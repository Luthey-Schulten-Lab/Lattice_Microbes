/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
 * Copyright 2012 Roberts Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.uiuc.edu/~schulten
 * 
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 * of the Software, and to permit persons to whom the Software is furnished to 
 * do so, subject to the following conditions:
 * 
 * - Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimers.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimers in the documentation 
 * and/or other materials provided with the distribution.
 * 
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Mike Hallock, Zane Thornburg
 */

#include <map>
#include <string>
#include <cstdlib>
#include "config.h"
#include <chrono>
#if defined(MACOSX)
#include <mach/mach_time.h>
#elif defined(LINUX)
#include <time.h>
#endif
#include "cuda/lm_cuda.h"
#include "cuda/ldg.h"
#include "core/Math.h"
#include "core/Print.h"
#include "cme/CMESolver.h"
#include "DiffusionModel.pb.h"
#include "Lattice.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "core/ResourceAllocator.h"
#include "rdme/ByteLattice.h"
#include "rdme/MGPUMpdRdmeSolver.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include "core/Timer.h"

#define MPD_WORDS_PER_SITE              (MPD_LATTICE_MAX_OCCUPANCY / 4)
#define MPD_APRON_SIZE                  1

#include "cuda/constant.cuh"
namespace lm {
namespace rdme {
namespace mgpumpdrdme_dev{
#include "rdme/dev/xor_random_dev.cu"
#include "rdme/dev/lattice_sim_1d_dev.cu"
#include "rdme/dev/byte_diffusion_1d_dev.cu"
#include "rdme/dev/byte_reaction_dev.cu"
}}}

#include "rdme/GPUMapper/MultiGPUMapper.h"
#include "rdme/GPUMapper/SegmentDescriptor.h"
#include "rdme/GPUMapper/ZDivMultiGPUMapper.h"

#define L3INDEX(_x,_y,_z, _g) ((_x)+((_y)*(_g.x))+((_z)*(_g.x)*(_g.y)))


extern bool globalAbort;


using std::map;
using lm::io::DiffusionModel;
using lm::rdme::Lattice;
using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {


volatile unsigned int mclkr_overflow_counter=0;
volatile unsigned int mclkr_total_overflows=0;
volatile unsigned int mclkr_overflow_reporters=0;
unsigned int mclkr_worker_count;
pthread_mutex_t mclkr_overflow_mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t mclkr_overflow_cond=PTHREAD_COND_INITIALIZER;

volatile int mclkr_generation = 0;
pthread_spinlock_t mclkr_spin;


void* gpu_worker_thread(void *arg)
{
	Print::printf(Print::VERBOSE_DEBUG,"thread %p launched with arg %p",pthread_self(), arg);
	gpu_worker_thread_params *params=(gpu_worker_thread_params*)arg;
	MGPUMpdRdmeSolver *solver=params->runner;

	return solver->run_thread(params->gpu);
}

inline void mclkr_barrier_spin(int overflows)
{
	PROF_BEGIN(MPD_MCLKR_BARRIER);

	int gen = mclkr_generation;
	pthread_spin_lock(&mclkr_spin);

    // Update shared counters
    mclkr_overflow_reporters++;
    mclkr_overflow_counter+=overflows;

	pthread_spin_unlock(&mclkr_spin);

    // simulate barrier
    if(mclkr_overflow_reporters == mclkr_worker_count)
    {
        mclkr_total_overflows=mclkr_overflow_counter;
        mclkr_overflow_reporters=0;
		mclkr_generation++;
    }
    else
    {
		while(mclkr_generation == gen) {} 
    }

	PROF_END(MPD_MCLKR_BARRIER);
}

inline void mclkr_barrier_cond(int overflows)
{
	PROF_BEGIN(MPD_MCLKR_BARRIER);

	pthread_mutex_lock(&mclkr_overflow_mutex);

    // Update shared counters
    mclkr_overflow_reporters++;
    mclkr_overflow_counter+=overflows;

    // simulate barrier
    if(mclkr_overflow_reporters == mclkr_worker_count)
    {
        mclkr_total_overflows=mclkr_overflow_counter;
        mclkr_overflow_reporters=0;
		pthread_cond_broadcast(&mclkr_overflow_cond);
    }
    else
    {
		pthread_cond_wait(&mclkr_overflow_cond, &mclkr_overflow_mutex);
    }

	pthread_mutex_unlock(&mclkr_overflow_mutex);

	PROF_END(MPD_MCLKR_BARRIER);
}

	

MGPUMpdRdmeSolver::MGPUMpdRdmeSolver()
:RDMESolver(lm::rng::RandomGenerator::NONE),
 seed(0), tau(0.0), threads(NULL),
 overflowTimesteps(0), overflowListUses(0),
 model_reactionOrders(NULL), model_reactionSites(NULL), model_reactionRates(NULL),
 model_D1(NULL), model_D2(NULL),
 model_S(NULL), model_T(NULL), model_RL(NULL),
 zeroOrder(NULL), firstOrder(NULL), secondOrder(NULL)
{
	// set defaults features
	aggcopy_x_unpack = true;
	aggcopy_r_pack = true;
	use_spin_barrier = false;
}

void MGPUMpdRdmeSolver::initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resourcesA)
{
	pthread_spin_init(&mclkr_spin, 0);

	RDMESolver::initialize(replicate, parameters, resourcesA);

    // Figure out the random seed.
    uint32_t seedTop=(unsigned int)atoi((*parameters)["seed"].c_str());
    if (seedTop == 0)
    {
        #if defined(MACOSX)
        seedTop = (uint32_t)mach_absolute_time();
        #elif defined(LINUX)
        struct timespec seed_timespec;
        if (clock_gettime(CLOCK_REALTIME, &seed_timespec) != 0) throw lm::Exception("Error getting time to use for random seed.");
        seedTop = seed_timespec.tv_nsec;
        #endif
    }
    seed = (seedTop<<16)|(replicate&0x0000FFFF);
    Print::printf(Print::INFO, "MPDRDME: Rng seed: top word %u, bottom word %u", (seed>>16)*0xffff, seed&0xFFFF);

	resources=resourcesA;

	//Verify that core count is at least as large as the gpu count.
	//Otherwise, threads compete and speedups are not realized.
	int cpu_count = resources->cpuCores.size();
	int gpu_count = resources->cudaDevices.size();

	if(cpu_count < gpu_count)
	{
		Print::printf(Print::FATAL, "GPU count (%d) exceeds CPU count (%d)", gpu_count, cpu_count);
		throw lm::Exception("Not enough cpu cores to run solver");
	}

	string &s = (*parameters)["rdme.mpd.overflowhandler"];
	if(s == "classic")
	{
		overflow_handling = OVERFLOW_MODE_CLASSIC;
		Print::printf(Print::DEBUG, "Overflow handler set to classic");
	}
	else if(s == "relaxed")
	{
		overflow_handling = OVERFLOW_MODE_RELAXED;
		Print::printf(Print::DEBUG, "Overflow handler set to relaxed");
	}
	else if(s == "")
	{
		overflow_handling = OVERFLOW_MODE_CLASSIC;
		Print::printf(Print::DEBUG, "Overflow handler set to default (classic)");
	}
	else
	{
		overflow_handling = OVERFLOW_MODE_CLASSIC;
		Print::printf(Print::WARNING, "Unknown overflow handler requested: '%s'", s.c_str());
	}

	if((*parameters)["mgpu.x_unpack"] != "")
		aggcopy_x_unpack = atoi((*parameters)["mgpu.x_unpack"].c_str());

	if((*parameters)["mgpu.r_pack"] != "")
		aggcopy_r_pack = atoi((*parameters)["mgpu.r_pack"].c_str());

	if((*parameters)["mgpu.spinlock"] != "")
		use_spin_barrier = atoi((*parameters)["mgpu.spinlock"].c_str());
}

void MGPUMpdRdmeSolver::initialize_decomposition()
{
	// Create mapper
	int latx, laty, latz;
	int apron, overlap;
	int gpus = resources->cudaDevices.size();

#if defined MPD_BOUNDARY_PERIODIC
	bool pz = true;
#else
	bool pz = false;
#endif

	apron = 2;
	overlap = 0;

	latx = lattice->getXSize();
	laty = lattice->getYSize();
	latz = lattice->getZSize();

	int site_size = sizeof(int);
	int pagecount = MPD_WORDS_PER_SITE;

	mapper = new ZDivMultiGPUMapper(latx, laty, latz,
	                                site_size, apron, overlap,
								    gpus, resources->cudaDevices.data(), pz, pagecount);

    pthread_barrier_init(&start_barrier,      NULL, gpus+1);
    pthread_barrier_init(&stop_barrier,       NULL, gpus+1);
    pthread_barrier_init(&simulation_barrier, NULL, gpus);

	mclkr_worker_count = gpus;
}

void MGPUMpdRdmeSolver::start_threads()
{
	// TODO FIXME: Think about this.  Perhaps each thread should be assigned a core
	// or perhaps we should be reusing the threading infrastructure already provided
	// to us via the core.  We also need NUMA awareness, which is another TODO.
	// For now, lets just open up our pinning to be over all the cores that we've
	// been graciously handed to us.
#if defined(LINUX) && !defined(ARM)
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);

	for (int c = 0; c < resources->cpuCores.size(); c++)
	{
		CPU_SET(resources->cpuCores[c], &cpuset);
	}
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif

	int gpus = resources->cudaDevices.size();
	threads = new gpu_worker_thread_params[gpus];

	for (int i = 0; i < gpus; i++)
	{
        threads[i].runner           = this;
        threads[i].mapper           = mapper;
        threads[i].gpu              = i;
        threads[i].ngpus            = gpus;
        threads[i].timesteps_to_run = 0;
		threads[i].segment          = mapper->getSegmentDescriptor(i);
		
        pthread_create(&(threads[i].thread), NULL,
		               &gpu_worker_thread, (void*)&(threads[i]));
	}
}

void MGPUMpdRdmeSolver::stop_threads()
{
	int gpus=resources->cudaDevices.size();
	for(int i=0; i<gpus; i++)
	{
		Print::printf(Print::DEBUG,"Waiting to join GPU thread %d",i);
		pthread_join(threads[i].thread, NULL);
	}
}

MGPUMpdRdmeSolver::~MGPUMpdRdmeSolver()
{
	// Free allocated memory in MGPUMpdRdmeSolver
	if (model_reactionOrders) {delete [] model_reactionOrders;}
	if (model_reactionSites)  {delete [] model_reactionSites;}
	if (model_reactionRates)  {delete [] model_reactionRates;}

	if (model_D1) {delete [] model_D1;}
	if (model_D2) {delete [] model_D2;}

	if (model_S)  {delete [] model_S;}
	if (model_T)  {delete [] model_T;}
	if (model_RL) {delete [] model_RL;}

	if (zeroOrder)   {delete [] zeroOrder;}
	if (firstOrder)  {delete [] firstOrder;}
	if (secondOrder) {delete [] secondOrder;}

	if (threads) {delete [] threads;}
}

void MGPUMpdRdmeSolver::allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing)
{
	assert(bytes_per_particle == 1);
    lattice = (Lattice *)new ByteLattice(latticeXSize, latticeYSize, latticeZSize, latticeSpacing, particlesPerSite);
}

void MGPUMpdRdmeSolver::buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionTypesA, const double * KA, const int * SA, const uint * DA, const uint kCols)
{
    CMESolver::buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypesA, KA, SA, DA, kCols);

    // Get the time step.
    tau=atof((*parameters)["timestep"].c_str());
    if (tau <= 0.0) throw InvalidArgException("timestep", "A positive timestep must be specified for the solver.");

    // Make sure we can support the reaction model.
    if (numberReactions > MPD_MAX_REACTION_TABLE_ENTRIES) throw Exception("The number of reaction table entries exceeds the maximum supported by the solver.");
#ifndef MPD_GLOBAL_S_MATRIX
    if (numberSpecies*numberReactions > MPD_MAX_S_MATRIX_ENTRIES) throw Exception("The number of S matrix entries exceeds the maximum supported by the solver.");
#endif

    // Setup the cuda reaction model.
    model_reactionOrders = new unsigned int[numberReactions];
    model_reactionSites = new unsigned int[numberReactions];
    model_D1 = new unsigned int[numberReactions];
    model_D2 = new unsigned int[numberReactions];
    for (uint i=0; i<numberReactions; i++)
    {
        if(reactionTypes[i] == ZerothOrderPropensityArgs::REACTION_TYPE) {
        	model_reactionOrders[i] = MPD_ZERO_ORDER_REACTION;
        	model_reactionSites[i] = 0;
        	model_D1[i] = 0; 
        	model_D2[i] = 0;
        }
    	else if (reactionTypes[i] == FirstOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionOrders[i] = MPD_FIRST_ORDER_REACTION;
    		model_reactionSites[i] = 0;
    		model_D1[i] = ((FirstOrderPropensityArgs *)propensityFunctionArgs[i])->si+1;
    		model_D2[i] = 0;
    	}
    	else if (reactionTypes[i] == SecondOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionOrders[i] = MPD_SECOND_ORDER_REACTION;
    		model_reactionSites[i] = 0;
    		model_D1[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s1i+1;
    		model_D2[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s2i+1;
    	}
    	else if (reactionTypes[i] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionOrders[i] = MPD_SECOND_ORDER_SELF_REACTION;
    		model_reactionSites[i] = 0;
    		model_D1[i] = ((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i])->si+1;
    		model_D2[i] = 0;
    	}
    	else
    	{
    		throw InvalidArgException("reactionTypeA", "the reaction type was not supported by the solver", reactionTypes[i]);
    	}
    }

    // Setup the cuda S matrix.
    model_S = new int8_t[numberSpecies*numberReactions];

	for(uint rx = 0; rx < numberReactions; rx++)
	{
		for (uint p=0; p<numberSpecies; p++)
		{
			//tmpS[numberSpecies*rx + p] = S[numberSpecies*rx + p];
			model_S[rx * numberSpecies + p] = S[numberReactions*p + rx];
		}
	}
}

void MGPUMpdRdmeSolver::buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData)
{
    RDMESolver::buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData);

    // Get the time step.
    tau=atof((*parameters)["timestep"].c_str());
    if (tau <= 0.0) throw InvalidArgException("timestep", "A positive timestep must be specified for the solver.");

    // Setup the cuda transition matrix.
    const size_t DFmatrixSize = numberSpecies*numberSiteTypes*numberSiteTypes;

#ifndef MPD_GLOBAL_T_MATRIX
    if (DFmatrixSize > MPD_MAX_TRANSITION_TABLE_ENTRIES) throw Exception("The number of transition table entries exceeds the maximum supported by the solver.");
#endif

#ifndef MPD_GLOBAL_S_MATRIX
    if (numberReactions*numberSiteTypes > MPD_MAX_RL_MATRIX_ENTRIES) throw Exception("The number of RL matrix entries exceeds the maximum supported by the solver.");
#endif

    // Calculate the probability from the diffusion coefficient and the lattice properties.
    /*
     * p0 = probability of staying at the site, q = probability of moving in plus or minus direction
     *
     * D=(1-p0)*lambda^2/2*tau
     * q=(1-p0)/2
     * D=2q*lambda^2/2*tau
     * q=D*tau/lambda^2
     */
    model_T = new float[DFmatrixSize];
    for (uint i=0; i<DFmatrixSize; i++)
    {
        float q=(float)(DF[i]*tau/pow(latticeSpacing,2));
        if (q > 0.50f) throw InvalidArgException("D", "The specified diffusion coefficient is too high for the diffusion model.");
        model_T[i] = q;
    }

    // Setup the cuda reaction location matrix.
    model_RL = new uint8_t[numberReactions*numberSiteTypes];
    for (uint i=0; i<numberReactions*numberSiteTypes; i++)
    {
    	model_RL[i] = RL[i];
    }

    // Set the cuda reaction model rates now that we have the subvolume size.
    model_reactionRates = new float[numberReactions];

	// Set up pre-configured propensity matrices
	zeroOrderSize=numberSiteTypes;
	firstOrderSize=numberSpecies*numberSiteTypes;
	secondOrderSize=numberSpecies*numberSpecies*numberSiteTypes;
	zeroOrder=new float[zeroOrderSize];
	firstOrder=new float[firstOrderSize];
	secondOrder=new float[secondOrderSize];

	computePropensities();
	reactionModelModified = false;
}

void MGPUMpdRdmeSolver::computePropensities()
{
	unsigned int latticeXSize = lattice->getXSize();
	unsigned int latticeYSize = lattice->getYSize();
	unsigned int latticeZSize = lattice->getZSize();
    for (uint i=0; i<numberReactions; i++)
    {
    	if (reactionTypes[i] == ZerothOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((ZerothOrderPropensityArgs *)propensityFunctionArgs[i])->k*tau/(latticeXSize*latticeYSize*latticeZSize);
    	}
    	else if (reactionTypes[i] == FirstOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((FirstOrderPropensityArgs *)propensityFunctionArgs[i])->k*tau;
    	}
    	else if (reactionTypes[i] == SecondOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->k*tau*latticeXSize*latticeYSize*latticeZSize;
    	}
    	else if (reactionTypes[i] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i])->k*tau*latticeXSize*latticeYSize*latticeZSize;
    	}
    	else
    	{
    		throw InvalidArgException("reactionTypeA", "the reaction type was not supported by the solver", reactionTypes[i]);
    	}
    }

    float scale=latticeXSize*latticeYSize*latticeZSize;
	for(uint site=0; site<numberSiteTypes; site++)
	{
		uint o1=site*numberSpecies;
		uint o2=site*numberSpecies*numberSpecies;

		zeroOrder[site]=0.0f;
		for (uint i=0; i<numberSpecies; i++)
		{
			firstOrder[o1 + i]=0.0f;
			for (uint j=0; j<numberSpecies; j++)
				secondOrder[o2 + i*numberSpecies + j]=0.0f;	
		}

		for (uint i=0; i<numberReactions; i++)
		{
			if(! RL[i*numberSiteTypes + site])
				continue;

			switch(reactionTypes[i])
			{
				case ZerothOrderPropensityArgs::REACTION_TYPE:
				{
					ZerothOrderPropensityArgs *rx=(ZerothOrderPropensityArgs *)propensityFunctionArgs[i];
					zeroOrder[site]+=rx->k*tau/scale;
				}   break;

				case FirstOrderPropensityArgs::REACTION_TYPE:
				{
					FirstOrderPropensityArgs *rx=(FirstOrderPropensityArgs *)propensityFunctionArgs[i];
					firstOrder[o1 + rx->si]+=rx->k*tau;
				}   break;
				
				case SecondOrderPropensityArgs::REACTION_TYPE:
				{
					SecondOrderPropensityArgs *rx=(SecondOrderPropensityArgs *)propensityFunctionArgs[i];
					secondOrder[o2 + (rx->s1i) * numberSpecies + (rx->s2i)]=rx->k*tau*scale;
					secondOrder[o2 + (rx->s2i) * numberSpecies + (rx->s1i)]=rx->k*tau*scale;
				}   break; 

				case SecondOrderSelfPropensityArgs::REACTION_TYPE:
				{
					SecondOrderSelfPropensityArgs *rx=(SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i];
					secondOrder[o2 + (rx->si) * numberSpecies + (rx->si)]+=rx->k*tau*scale*2;
				}
			}
		}
	}
}

void MGPUMpdRdmeSolver::copyModelsToDevice(int gpu)
{
    const size_t DFmatrixSize = numberSpecies*numberSiteTypes*numberSiteTypes;
	const unsigned int latticeXSize=threads[gpu].segment->local_dimensions.x;
	const unsigned int latticeYSize=threads[gpu].segment->local_dimensions.y;
	const unsigned int latticeZSize=threads[gpu].segment->local_dimensions.z;
	
	// the lattice{X,XY,XYZ}SizeC refer to sizes of the LOCAL lattice.
	// NOTE that since we have 1D diffusion only, the global X and XY sizes are the same.
    const unsigned int latticeXYSize = latticeXSize*latticeYSize;
    const unsigned int latticeXYZSize = latticeXSize*latticeYSize*latticeZSize;

	// gloal_latticeXYZSizeC is the full lattice XYZ size.  
    const unsigned int global_latticeZSize = lattice->getZSize();
    const unsigned int global_latticeXYZSize = latticeXSize*latticeYSize*lattice->getZSize();

    // Copy the reaction model and S matrix to constant memory on the GPU.
#ifdef MPD_GLOBAL_R_MATRIX
    // R matrix put in global memory
    //cudaMalloc(&numberReactionsG, sizeof(unsigned int));
    //cudaMemcpy(numberReactionsG, &numberReactions, sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(numberReactionsC, &numberReactions, sizeof(unsigned int));
    cudaMalloc(&(threads[gpu].reactionOrdersG), numberReactions*sizeof(unsigned int));
    cudaMemcpy(threads[gpu].reactionOrdersG, model_reactionOrders, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&(threads[gpu].reactionSitesG), numberReactions*sizeof(unsigned int));
    cudaMemcpy(threads[gpu].reactionSitesG, model_reactionSites, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&(threads[gpu].D1G), numberReactions*sizeof(unsigned int));
    cudaMemcpy(threads[gpu].D1G, model_D1, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&(threads[gpu].D2G), numberReactions*sizeof(unsigned int));
    cudaMemcpy(threads[gpu].D2G, model_D2, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
#else
    cudaMemcpyToSymbol(numberReactionsC, &numberReactions, sizeof(unsigned int));
    cudaMemcpyToSymbol(reactionOrdersC, model_reactionOrders, numberReactions*sizeof(unsigned int));
    cudaMemcpyToSymbol(reactionSitesC, model_reactionSites, numberReactions*sizeof(unsigned int));
    cudaMemcpyToSymbol(D1C, model_D1, numberReactions*sizeof(unsigned int));
    cudaMemcpyToSymbol(D2C, model_D2, numberReactions*sizeof(unsigned int));
#endif

#ifdef MPD_GLOBAL_S_MATRIX
	// If S matrix stored in global memory, allocate space and perform copy
	cudaMalloc(&(threads[gpu].SG), numberSpecies*numberReactions * sizeof(int8_t));
	cudaMemcpy(threads[gpu].SG, model_S, numberSpecies*numberReactions * sizeof(int8_t), cudaMemcpyHostToDevice);
#else
	// S matrix is in constant memory
    cudaMemcpyToSymbol(SC, model_S, numberSpecies*numberReactions*sizeof(int8_t));
#endif

#ifdef MPD_GLOBAL_T_MATRIX
	cudaMalloc(&(threads[gpu].TG), DFmatrixSize*sizeof(float));
	cudaMemcpy(threads[gpu].TG, model_T, DFmatrixSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(TC, &(threads[gpu].TG), sizeof(float*));
#else
    cudaMemcpyToSymbol(TC, model_T, DFmatrixSize*sizeof(float));
#endif

    cudaMemcpyToSymbol(numberSpeciesC, &numberSpecies, sizeof(numberSpeciesC));
    cudaMemcpyToSymbol(numberSiteTypesC, &numberSiteTypes, sizeof(numberSiteTypesC));
    cudaMemcpyToSymbol(latticeXSizeC, &latticeXSize, sizeof(latticeYSizeC));
    cudaMemcpyToSymbol(latticeYSizeC, &latticeYSize, sizeof(latticeYSizeC));
    cudaMemcpyToSymbol(latticeZSizeC, &latticeZSize, sizeof(latticeZSizeC));
    cudaMemcpyToSymbol(latticeXYSizeC, &latticeXYSize, sizeof(latticeXYSizeC));
    cudaMemcpyToSymbol(latticeXYZSizeC, &latticeXYZSize, sizeof(latticeXYZSizeC));
    cudaMemcpyToSymbol(global_latticeZSizeC, &global_latticeZSize, sizeof(global_latticeZSizeC));
    cudaMemcpyToSymbol(global_latticeXYZSizeC, &global_latticeXYZSize, sizeof(global_latticeXYZSizeC));
#ifdef MPD_GLOBAL_S_MATRIX
	// Store RL in global memory too, since I'm going to assume if S is too big, then RL is too.
	cudaMalloc(&(threads[gpu].RLG), numberReactions*numberSiteTypes * sizeof(uint8_t));
	cudaMemcpy(threads[gpu].RLG, model_RL, numberReactions*numberSiteTypes * sizeof(uint8_t), cudaMemcpyHostToDevice);
#else
	// RL is stored in constant memory
    cudaMemcpyToSymbol(RLC, model_RL, numberReactions*numberSiteTypes*sizeof(uint8_t));
#endif

#ifdef MPD_GLOBAL_R_MATRIX
    cudaMalloc(&(threads[gpu].reactionRatesG), numberReactions*sizeof(float));
    cudaMemcpy(threads[gpu].reactionRatesG, model_reactionRates, numberReactions*sizeof(float), cudaMemcpyHostToDevice);
#else
    cudaMemcpyToSymbol(reactionRatesC, model_reactionRates, numberReactions*sizeof(float));
#endif

	// Copy per-GPU global offsets
    //CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(lsizeXC, &(threads[gpu].local_dimensions.x), sizeof(unsigned int)));
    //CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(lsizeYC, &(threads[gpu].local_dimensions.y), sizeof(unsigned int)));
    //CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(lsizeZC, &(threads[gpu].local_dimensions.z), sizeof(unsigned int)));
    cudaMemcpyToSymbol(goffsetXC, &(threads[gpu].segment->global_offset.x), sizeof(int));
    cudaMemcpyToSymbol(goffsetYC, &(threads[gpu].segment->global_offset.y), sizeof(int));
    cudaMemcpyToSymbol(goffsetZC, &(threads[gpu].segment->global_offset.z), sizeof(int));

    cudaMemcpyToSymbol(gpuidC, &gpu, sizeof(unsigned int));

	cudaMalloc(&(threads[gpu].propZeroOrder), zeroOrderSize*sizeof(float));
	cudaMemcpy(threads[gpu].propZeroOrder, zeroOrder, zeroOrderSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc(&(threads[gpu].propFirstOrder), firstOrderSize*sizeof(float));
	cudaMemcpy(threads[gpu].propFirstOrder, firstOrder, firstOrderSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc(&(threads[gpu].propSecondOrder), secondOrderSize*sizeof(float));
	cudaMemcpy(threads[gpu].propSecondOrder, secondOrder, secondOrderSize*sizeof(float), cudaMemcpyHostToDevice);
	
	ZDivMultiGPUMapper *zmap = (ZDivMultiGPUMapper*)mapper;
	gpu_info *g = zmap->getinfo(gpu);
	unsigned int tmp = g->overlap_send[0].z;
	cudaMemcpyToSymbol(top_send_z, &tmp, sizeof(int));
	tmp = g->overlap_dim.z;
	cudaMemcpyToSymbol(top_dim_z, &tmp, sizeof(int));
	tmp = DIMSIZE(g->overlap_dim);
	cudaMemcpyToSymbol(top_dim_size, &tmp, sizeof(int));

	tmp = g->overlap_send[1].z;
	cudaMemcpyToSymbol(bot_send_z, &tmp, sizeof(int));
	tmp = g->overlap_dim.z;
	cudaMemcpyToSymbol(bot_dim_z, &tmp, sizeof(int));
	tmp = DIMSIZE(g->overlap_dim);
	cudaMemcpyToSymbol(bot_dim_size, &tmp, sizeof(int));

	tmp = g->overlap_recv[0].z;
	cudaMemcpyToSymbol(top_recv_z, &tmp, sizeof(int));

	tmp = g->overlap_recv[1].z;
	cudaMemcpyToSymbol(bot_recv_z, &tmp, sizeof(int));
}

void MGPUMpdRdmeSolver::setReactionRate(unsigned int rxid, float rate)
{
	if (reactionTypes[rxid] == ZerothOrderPropensityArgs::REACTION_TYPE)
	{
		((ZerothOrderPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
	else if (reactionTypes[rxid] == FirstOrderPropensityArgs::REACTION_TYPE)
	{
		((FirstOrderPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
	else if (reactionTypes[rxid] == SecondOrderPropensityArgs::REACTION_TYPE)
	{
		((SecondOrderPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
	else if (reactionTypes[rxid] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
	{
		((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
	else
    {
    	throw InvalidArgException("reactionTypeA", "the reaction type was not supported by the solver", reactionTypes[rxid]);
    }

	// Flag to trigger re-computation of rdme propensities
	reactionModelModified = true;
}

void MGPUMpdRdmeSolver::generateTrajectory()
{
	// Determine how to divide the lattice
	initialize_decomposition();

	// Start worker threads
	start_threads();

    // Shadow the lattice member as a cuda lattice.
    ByteLattice * lattice = (ByteLattice *)this->lattice;

    // Get the interval for writing species counts and lattices.
    uint32_t speciesCountsWriteInterval=atol((*parameters)["writeInterval"].c_str());
    uint32_t nextSpeciesCountsWriteTime = speciesCountsWriteInterval;
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpeciesToTrack);
    speciesCountsDataSet.set_number_entries(0);
    uint32_t latticeWriteInterval=atol((*parameters)["latticeWriteInterval"].c_str());
    uint32_t nextLatticeWriteTime = latticeWriteInterval;
    lm::io::Lattice latticeDataSet;

    // Get the simulation time limit. 
    double maxTime=atof((*parameters)["maxTime"].c_str());

    Print::printf(Print::INFO,
	              "Running mpd rdme simulation with %d species, %d reactions, %d site types for %e s with tau %e. Writing species at %d and lattice at %d intervals",
				  numberSpecies, numberReactions, numberSiteTypes,
				  maxTime, tau,
				  speciesCountsWriteInterval, latticeWriteInterval);
    
    // Set the initial time.
    double time = 0.0;
	current_timestep = 0;

	bool hookEnabled=false;
    uint32_t nextHookTime=0;
    uint32_t hookInterval=0;

    // Find out at what interval to hook simulations
	if((*parameters)["hookInterval"] != "")
    {
        hookInterval=atol((*parameters)["hookInterval"].c_str());
        hookEnabled=true;
        nextHookTime=hookInterval;
    }

	// Find out at what interval to write status messages
    printPerfInterval = 60;
    if((*parameters)["perfPrintInterval"] != "")
    {
        printPerfInterval=atof((*parameters)["perfPrintInterval"].c_str());
    }

    // Record the initial species counts.
    recordSpeciesCounts(time, lattice, &speciesCountsDataSet);

    // Write the initial lattice.
    writeLatticeData(time, lattice, &latticeDataSet);

	// Perform an initial hook check
    hookCheckSimulation(time, lattice);

    // Loop until we have finished the simulation.
    while (time < maxTime)
    {
		// compute number of steps to batch
		int nSteps = min(nextLatticeWriteTime-current_timestep,
			         min(nextSpeciesCountsWriteTime-current_timestep,
			         min(nextHookTime-current_timestep,
					 ((uint32_t) (ceil(maxTime/tau)-current_timestep)))));

		if(nSteps == 0) nSteps=1;
		assert(nSteps > 0);
		
		timesteps_to_run = nSteps;
	
		// wake up threads
		pthread_barrier_wait(&start_barrier);
		
		// ... and wait.
		pthread_barrier_wait(&stop_barrier);
		
		if(globalAbort)
		{
			printf("Global abort: terminating solver");
			break;
		}
	
        // Advance the timestep counter
		current_timestep += nSteps;

		// Update the time
        time = current_timestep*tau;

		// No need to sync from gpus first; each gpu thread does a stage_out
		// prior to entering the stop barrier
		// clear the flag for if reaction model has changed
		reactionModelModified = false;
		// Also no need to check the return value of the hook as a stage in is always done
		// at the beginning of the timestep batch

		// Check if we need to execute the hook
		if (hookEnabled && current_timestep >= nextHookTime)
        {
            // lattice->copyFromGPU();
            Print::printf(Print::INFO, "Hook time is %.14f, in steps is %d", time, nextHookTime);
            hookCheckSimulation(time, lattice);
            nextHookTime += hookInterval;
            Print::printf(Print::INFO, "Next hook time is %d", nextHookTime);
        }
			
        // See if we need to write the lattice.
        if (current_timestep >= nextLatticeWriteTime)
        {
            PROF_BEGIN(PROF_SERIALIZE_LATTICE);
            writeLatticeData(time, lattice, &latticeDataSet);
            nextLatticeWriteTime += latticeWriteInterval;
            PROF_END(PROF_SERIALIZE_LATTICE);
        }

        // See if we need to write the species counts.
        if (current_timestep >= nextSpeciesCountsWriteTime)
        {
            PROF_BEGIN(PROF_DETERMINE_COUNTS);
            recordSpeciesCounts(time, lattice, &speciesCountsDataSet);
            nextSpeciesCountsWriteTime += speciesCountsWriteInterval;
            PROF_END(PROF_DETERMINE_COUNTS);

            // See if we have accumulated enough species counts to send.
            if (speciesCountsDataSet.number_entries() >= TUNE_SPECIES_COUNTS_BUFFER_SIZE)
            {
                PROF_BEGIN(PROF_SERIALIZE_COUNTS);
                writeSpeciesCounts(&speciesCountsDataSet);
                PROF_END(PROF_SERIALIZE_COUNTS);
            }
        }
    }

    // Write any remaining species counts.
    writeSpeciesCounts(&speciesCountsDataSet);

	// take down threads
	timesteps_to_run=-1;
	Print::printf(Print::DEBUG, "all done, tear down");
	pthread_barrier_wait(&start_barrier);
	Print::printf(Print::DEBUG, "stop all threads");

	stop_threads();
}

void MGPUMpdRdmeSolver::writeLatticeData(double time, ByteLattice * lattice, lm::io::Lattice * latticeDataSet)
{
    Print::printf(Print::DEBUG, "Writing lattice at %e s", time);

    // Record the lattice data.
    latticeDataSet->Clear();
    latticeDataSet->set_lattice_x_size(lattice->getSize().x);
    latticeDataSet->set_lattice_y_size(lattice->getSize().y);
    latticeDataSet->set_lattice_z_size(lattice->getSize().z);
    latticeDataSet->set_particles_per_site(lattice->getMaxOccupancy());
    latticeDataSet->set_time(time);

    // Push it to the output queue.
    size_t payloadSize = size_t(lattice->getSize().x)*size_t(lattice->getSize().y)*size_t(lattice->getSize().z)*size_t(lattice->getMaxOccupancy())*sizeof(uint8_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::BYTE_LATTICE, replicate, latticeDataSet, lattice, payloadSize, &lm::rdme::ByteLattice::nativeSerialize);
}

void MGPUMpdRdmeSolver::writeLatticeSites(double time, ByteLattice * lattice)
{
    Print::printf(Print::DEBUG, "Writing lattice sites at %e s", time);

    lm::io::Lattice latticeDataSet;
    // Record the lattice data.
    latticeDataSet.Clear();
    latticeDataSet.set_lattice_x_size(lattice->getSize().x);
    latticeDataSet.set_lattice_y_size(lattice->getSize().y);
    latticeDataSet.set_lattice_z_size(lattice->getSize().z);
    latticeDataSet.set_particles_per_site(lattice->getMaxOccupancy());
    latticeDataSet.set_time(time);

    // Push it to the output queue.
    size_t payloadSize = lattice->getSize().x*lattice->getSize().y*lattice->getSize().z*sizeof(uint8_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SITE_LATTICE,
	                                                      replicate,
														  &latticeDataSet,
														  lattice,
														  payloadSize,
														  &lm::rdme::ByteLattice::nativeSerializeSites);
}

void MGPUMpdRdmeSolver::recordSpeciesCounts(double time, ByteLattice * lattice, lm::io::SpeciesCounts * speciesCountsDataSet)
{
    std::map<particle_t,uint> particleCounts = lattice->getParticleCounts();

    speciesCountsDataSet->set_number_entries(speciesCountsDataSet->number_entries()+1);
    speciesCountsDataSet->add_time(time);

    for (particle_t p=0; p<numberSpeciesToTrack; p++)
    {
        speciesCountsDataSet->add_species_count((particleCounts.count(p+1)>0)?particleCounts[p+1]:0);
    }
}

void MGPUMpdRdmeSolver::writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet)
{
    if (speciesCountsDataSet->number_entries() > 0)
    {
        // Push it to the output queue.
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, speciesCountsDataSet);

        // Reset the data set.
        speciesCountsDataSet->Clear();
        speciesCountsDataSet->set_number_species(numberSpeciesToTrack);
        speciesCountsDataSet->set_number_entries(0);
    }
}

void MGPUMpdRdmeSolver::hookCheckSimulation(double time, ByteLattice * lattice)
{
	switch(hookSimulation(time, lattice))
    {
        case 0:
            break;

        case 1:
            // lattice->copyToGPU();
            break;

        case 2:
            // lattice->copyToGPU();
            writeLatticeSites(time, lattice);
            break;

        case 3:
            printf("hook return value is 3, force to stop.\n");
			writeLatticeSites(time, lattice);
            return;

        default:
            throw("Unknown hook return value");
    }

	if(reactionModelModified)
		computePropensities();
}

uint64_t MGPUMpdRdmeSolver::getTimestepSeed(uint32_t timestep, uint32_t substep)
{
    uint64_t timestepHash = (((((uint64_t)seed)<<30)+timestep)<<2)+substep;
    timestepHash = timestepHash * 3202034522624059733ULL + 4354685564936845319ULL;
    timestepHash ^= timestepHash >> 20; timestepHash ^= timestepHash << 41; timestepHash ^= timestepHash >> 5;
    timestepHash *= 7664345821815920749ULL;
    return timestepHash;
}

int MGPUMpdRdmeSolver::run_next_timestep(int gpu, uint32_t timestep)
{
	int z_start = 0;

	unsigned int *dlat        = threads[gpu].dLattice;
	unsigned int *dtmp        = threads[gpu].dLatticeTmp;
	uint8_t      *sites       = threads[gpu].dSites;
	unsigned int *d_overflows = threads[gpu].d_overflows;

	cudaStream_t cudaStream = threads[gpu].stream1;
	cudaStream_t copyStream = threads[gpu].stream2;
	
    dim3 gridSize        = threads[gpu].grid_x;
	dim3 threadBlockSize = threads[gpu].threads_x;

	ZDivMultiGPUMapper* zmap = (ZDivMultiGPUMapper*) mapper;

	// Execute the diffusion kernel for the x direction.
	if (aggcopy_x_unpack)
	{
		unsigned int* rbuf_top = zmap->getrbuf(gpu, timestep, 0);
		unsigned int* rbuf_bot = zmap->getrbuf(gpu, timestep, 1);

		mgpumpdrdme_dev::MGPU_x_kernel_unpack<<<gridSize, threadBlockSize>>>
		(dlat, sites, dtmp, z_start,
		 getTimestepSeed(timestep, 0),
		 d_overflows,
		 rbuf_top, rbuf_bot);
	}
	else
	{
		mapper->schedule_recv(gpu, dlat, timestep, 0, cudaStream);
		mapper->schedule_recv(gpu, dlat, timestep, 1, cudaStream);

		CUDA_EXCEPTION_EXECUTE((
			mgpumpdrdme_dev::MGPU_x_kernel<<<gridSize, threadBlockSize, 0, cudaStream>>>
			(dlat, sites, dtmp, z_start,
			 getTimestepSeed(timestep, 0),
			 d_overflows)
		));
	}

    gridSize        = threads[gpu].grid_y;
	threadBlockSize = threads[gpu].threads_y;

    // Execute the diffusion kernel for the y direction.
    mgpumpdrdme_dev::MGPU_y_kernel<<<gridSize, threadBlockSize>>>
	(dtmp, sites, dlat,
	 getTimestepSeed(timestep, 1),
	 d_overflows);

    gridSize        = threads[gpu].grid_z;
	threadBlockSize = threads[gpu].threads_z;

    // Execute the diffusion kernel for the z direction.
    mgpumpdrdme_dev::MGPU_z_kernel<<<gridSize, threadBlockSize>>>
	(dlat, sites, dtmp,
	 getTimestepSeed(timestep, 2),
	 d_overflows,
	 threads[gpu].segment->active_offset.z);

    if (numberReactions > 0)
    {
		// Execute the kernel for the reaction.
		// This kernel updates the lattice in-place, so only the src pointer is passed.
		gridSize        = threads[gpu].grid_r;
		threadBlockSize = threads[gpu].threads_r;

		unsigned int* buf_top = zmap->gettbuf(gpu, timestep, 0);
		unsigned int* buf_bot = zmap->gettbuf(gpu, timestep, 1);

	#ifdef MPD_FREAKYFAST
		if (aggcopy_r_pack)
		{
			mgpumpdrdme_dev::MGPU_precomp_reaction_kernel_packing<<<gridSize, threadBlockSize>>>
			(dtmp, sites, dtmp,
			 getTimestepSeed(timestep, 3),
			 d_overflows,
			 threads[gpu].segment->active_offset.z,
	        #ifdef MPD_GLOBAL_S_MATRIX
			 threads[gpu].SG,
			 threads[gpu].RLG,
			#endif
			#ifdef MPD_GLOBAL_R_MATRIX
			 threads[gpu].reactionOrdersG,
			 threads[gpu].reactionSitesG,
			 threads[gpu].D1G,
			 threads[gpu].D2G,
			 threads[gpu].reactionRatesG,
			#endif
			 threads[gpu].propZeroOrder,
			 threads[gpu].propFirstOrder,
			 threads[gpu].propSecondOrder,
			 buf_top, buf_bot);
		}
		else
		{
			mgpumpdrdme_dev::MGPU_precomp_reaction_kernel<<<gridSize, threadBlockSize>>>
			(dtmp, sites, dtmp,
			 getTimestepSeed(timestep, 3),
			 d_overflows,
			 threads[gpu].segment->active_offset.z,
			#ifdef MPD_GLOBAL_S_MATRIX
			 threads[gpu].SG,
			 threads[gpu].RLG,
			#endif
			#ifdef MPD_GLOBAL_R_MATRIX
			 threads[gpu].reactionOrdersG,
			 threads[gpu].reactionSitesG,
			 threads[gpu].D1G,
			 threads[gpu].D2G,
			 threads[gpu].reactionRatesG,
			#endif
			 threads[gpu].propZeroOrder,
			 threads[gpu].propFirstOrder,
			 threads[gpu].propSecondOrder);
		}
	#else
    	mgpumpdrdme_dev::MGPU_reaction_kernel<<<gridSize, threadBlockSize>>>
		(dtmp, sites, dtmp,
		 getTimestepSeed(timestep, 3),
		 d_overflows,
		 threads[gpu].segment->active_offset.z,
		#ifdef MPD_GLOBAL_S_MATRIX
		 threads[gpu].SG,
		 threads[gpu].RLG,
		#endif
		#ifdef MPD_GLOBAL_R_MATRIX
         threads[gpu].reactionOrdersG,
		 threads[gpu].reactionSitesG,
		 threads[gpu].D1G,
		 threads[gpu].D2G,
		 threads[gpu].reactionRatesG
        #endif
		);
	#endif
    }

    if (overflow_handling == OVERFLOW_MODE_RELAXED)
		mgpumpdrdme_dev::correct_overflows_mgpu<<<dim3(1,1,1),
	                                              dim3(TUNE_MPD_MAX_PARTICLE_OVERFLOWS,1,1),
												  0, cudaStream>>>(dtmp, d_overflows);

	pthread_barrier_wait(&simulation_barrier);

	// send data
	if (aggcopy_r_pack && numberReactions > 0)
	{
		zmap->publish_state(gpu, timestep, cudaStream, copyStream);
	}
	else
	{
		zmap->publish_state(gpu, timestep, cudaStream, copyStream, dtmp);
	}

    if (overflow_handling == OVERFLOW_MODE_RELAXED)
	{
		uint unresolved = ((unsigned int*)d_overflows)[0];
		if (unresolved > 0)
		{
            Print::printf(Print::WARNING, "[gpu %d] %d unresolved overflows", gpu, unresolved);
		}
	}

	// Save correct pointers for current data
	threads[gpu].dLattice    = dtmp;
	threads[gpu].dLatticeTmp = dlat;

	if (overflow_handling == OVERFLOW_MODE_CLASSIC)
	{
	#ifndef MPD_MAPPED_OVERFLOWS
		cudaMemcpy(threads[gpu].h_overflows, d_overflows, MPD_OVERFLOW_LIST_SIZE, cudaMemcpyDeviceToHost);
	#endif

		return threads[gpu].h_overflows[0];
	}
	
	// relaxed overflow mode
	return 0;
}


int MGPUMpdRdmeSolver::handle_overflows(int gpu, void *hptr, void *dptr, int ts)
{
	// Copy data back from GPUs
	mapper->stage_out(gpu, hptr, dptr);

	pthread_barrier_wait(&simulation_barrier);

	// Resolve overflows
	if (gpu == 0)
	{
		handle_all_overflows();
		mclkr_total_overflows  = 0;
		mclkr_overflow_counter = 0;
	}

	pthread_barrier_wait(&simulation_barrier);

	cudaStream_t cudaStream = threads[gpu].stream1;
	cudaStream_t copyStream = threads[gpu].stream2;

	// Return to GPU
	mapper->stage_in(gpu, dptr, hptr);
	mapper->publish_state(gpu, ts, cudaStream, copyStream, dptr);
	cudaMemset(threads[gpu].d_overflows, 0, sizeof(unsigned int));

	pthread_barrier_wait(&simulation_barrier);

	return 0;
}

int MGPUMpdRdmeSolver::handle_all_overflows()
{
    // Handle any particle overflows.
    PROF_BEGIN(PROF_MPD_OVERFLOW);

    overflowTimesteps++;
    uint32_t overflowList[1+2*TUNE_MPD_MAX_PARTICLE_OVERFLOWS];
    uint numberExceptions = 0;

	// LOOP OVER EACH GPU
	int gpus=resources->cudaDevices.size();
	for(int i=0; i<gpus; i++)
	{
		uint* gpulist=threads[i].h_overflows;
		ssize_t min_index = mapper->get_authority_offset(i);
		ssize_t max_index = DIMSIZE(mapper->get_global_dim(i))+min_index;
		//printf("gpu %d: min %d max %d\n", i, min_index, max_index);

		for(int e=0; e<gpulist[0]; e++)
		{
			// DECODE INDEX
            size_t latticeIndex = gpulist[e*2 + 1];
            particle_t particle = gpulist[e*2 + 2];
			
			ssize_t globalIndex = latticeIndex + mapper->get_global_input_offset(i);
			//printf("Overflow of particle type %d at local index %d on gpu %d ==> index %d\n", particle, latticeIndex, i, globalIndex);
			if(latticeIndex < min_index || latticeIndex > max_index)
				continue;

			// FILL list
            overflowList[(numberExceptions*2)+1] = globalIndex;
			overflowList[(numberExceptions*2)+2] = particle;
			numberExceptions++;
		}
	}

    if (numberExceptions > 0)
    {
        Print::printf(Print::DEBUG, "%d overflows", numberExceptions);
        
        // Make sure we did not exceed the overflow buffer.
        if (numberExceptions > TUNE_MPD_MAX_PARTICLE_OVERFLOWS)
            throw Exception("Too many particle overflows for the available buffer", numberExceptions);
            
        // Go through each exception.
        for (uint i=0; i<numberExceptions; i++)
        {
            // Extract the index and particle type.
            lattice_size_t latticeIndex = overflowList[(i*2)+1];
            particle_t particle = overflowList[(i*2)+2];
            
            // Get the x, y, and z coordiantes.
            lattice_size_t x = latticeIndex%lattice->getXSize();
            lattice_size_t y = (latticeIndex/lattice->getXSize())%lattice->getYSize();
            lattice_size_t z = latticeIndex/(lattice->getXSize()*lattice->getYSize());
            
            // Put the particles back into a nearby lattice site.
            bool replacedParticle = false;
            for (uint searchRadius=0; !replacedParticle && searchRadius <= TUNE_MPD_MAX_OVERFLOW_REPLACEMENT_DIST; searchRadius++)
            {
                // Get the nearby sites.
                std::vector<lattice_coord_t> sites = lattice->getNearbySites(x,y,z,(searchRadius>0)?searchRadius-1:0,searchRadius);
                
                // TODO: Shuffle the sites.
                
                // Try to find one that in not fully occupied and of the same type.
                for (std::vector<lattice_coord_t>::iterator it=sites.begin(); it<sites.end(); it++)
                {
                    lattice_coord_t site = *it;
                    if (lattice->getOccupancy(site.x,site.y,site.z) < lattice->getMaxOccupancy() && lattice->getSiteType(site.x,site.y,site.z) == lattice->getSiteType(x,y,z))
                    {
                        lattice->addParticle(site.x, site.y, site.z, particle);
                        replacedParticle = true;
                        Print::printf(Print::VERBOSE_DEBUG, "Handled overflow of particle %d at site %d,%d,%d type=%d occ=%d by placing at site %d,%d,%d type=%d newocc=%d dist=%0.2f", particle, x, y, z, lattice->getSiteType(x,y,z), lattice->getOccupancy(x,y,z), site.x, site.y, site.z, lattice->getSiteType(site.x,site.y,site.z), lattice->getOccupancy(site.x,site.y,site.z), sqrt(pow((double)x-(double)site.x,2.0)+pow((double)y-(double)site.y,2.0)+pow((double)z-(double)site.z,2.0)));
                        break;
                    }
                }
            }
            
            // If we were not able to fix the exception, throw an error.
            if (!replacedParticle)
                throw Exception("Unable to find an available site to handle a particle overflow.");
        }
        
        // Track that we used the overflow list.
        overflowListUses++;
    }
    
    // If the overflow lsit is being used too often, print a warning.
    if (overflowTimesteps >= 1000)
    {
        if (overflowListUses > 10)
            Print::printf(Print::WARNING, "%d uses of the particle overflow list in the last 1000 timesteps, performance may be degraded.", overflowListUses);
        overflowTimesteps = 0;
        overflowListUses = 0;
    }
    PROF_END(PROF_MPD_OVERFLOW);

	return 0;
}

void* MGPUMpdRdmeSolver::run_thread(int gpu)
{
	cudaSetDevice(gpu);

	mapper->numa_bind_thread(gpu);
	mapper->initialize_gpu(gpu);
	pthread_barrier_wait(&simulation_barrier);
	mapper->set_affinity(gpu);
	
	gpu_worker_thread_params *p = &threads[gpu];

	// allocate device memory for lattice stuctures
	cudaMalloc(&p->dLattice,    mapper->get_local_size(gpu) * MPD_WORDS_PER_SITE);
	cudaMalloc(&p->dLatticeTmp, mapper->get_local_size(gpu) * MPD_WORDS_PER_SITE);
	cudaMalloc(&p->dSites,      DIMSIZE(mapper->get_local_dim(gpu)));

	// allocate device and host space for overflow lists
	if (overflow_handling == OVERFLOW_MODE_CLASSIC)
	{
	#ifdef MPD_MAPPED_OVERFLOWS
	    cudaHostAlloc(&p->h_overflows, MPD_OVERFLOW_LIST_SIZE,
					  cudaHostAllocPortable|cudaHostAllocMapped);
		p->d_overflows = p->h_overflows;
	#else
	    cudaHostAlloc(&p->h_overflows, MPD_OVERFLOW_LIST_SIZE, cudaHostAllocPortable);
	    cudaMalloc(&p->d_overflows, MPD_OVERFLOW_LIST_SIZE);
		cudaMemsetAsync(p->d_overflows, 0, MPD_OVERFLOW_LIST_SIZE);
	#endif

		memset(p->h_overflows, 0, MPD_OVERFLOW_LIST_SIZE);
	}
	else
	{
		p->h_overflows = NULL;
	    cudaMalloc(&p->d_overflows, MPD_OVERFLOW_LIST_SIZE);
		cudaMemsetAsync(p->d_overflows, 0, MPD_OVERFLOW_LIST_SIZE);
	}

    // Create data exchange streams and events
    cudaStreamCreate(&p->stream1);
	cudaStreamCreate(&p->stream2);

	dim3 ldim = mapper->get_local_dim(gpu);
	dim3 gdim = mapper->get_global_dim(gpu);

	// calculate launch params
    calculateXLaunchParameters(&(p->grid_x), &p->threads_x,
	                           TUNE_MPD_X_BLOCK_MAX_X_SIZE,
							   ldim.x, ldim.y, ldim.z);
    calculateYLaunchParameters(&p->grid_y, &p->threads_y,
	                           TUNE_MPD_Y_BLOCK_X_SIZE, TUNE_MPD_Y_BLOCK_Y_SIZE,
							   ldim.x, ldim.y, ldim.z);
    calculateZLaunchParameters(&p->grid_z, &p->threads_z,
	                           TUNE_MPD_Z_BLOCK_X_SIZE, TUNE_MPD_Z_BLOCK_Z_SIZE,
							   ldim.x, ldim.y, gdim.z);
	calculateReactionLaunchParameters(&p->grid_r, &p->threads_r,
	                                  TUNE_MPD_REACTION_BLOCK_X_SIZE, TUNE_MPD_REACTION_BLOCK_Y_SIZE,
									  ldim.x, ldim.y, gdim.z);

	// Now here, do not copy model every time
	copyModelsToDevice(gpu);

	// Copy the initialized data into GPU memory
	mapper->stage_in(gpu, p->dLattice, ((ByteLattice*) lattice)->getParticlesMemory());
	mapper->stage_in_sites(gpu, p->dSites, ((ByteLattice*) lattice)->getSitesMemory());
	mapper->publish_state(gpu, current_timestep + 1, p->stream1, p->stream2, p->dLattice);
	
	Timer timer;
    double lastT = 0;
    int lastSteps = 0;
	size_t maxTS = atof((*parameters)["maxTime"].c_str()) / tau;
	if (gpu == 0)
		timer.tick();

	do
	{
		// wait for master wakeup
		pthread_barrier_wait(&start_barrier);

		int steps = timesteps_to_run;
		if (steps < 0)
			break;

		// Don't copy models in on every step set
		// This was to accomodate changing lattice dims from lb changes
		// Unless the reaction model has changed
		if (reactionModelModified) 
		{
			copyModelsToDevice(gpu);
		}

		for (int ts = 0; ts < steps; ts++)
		{
			if (gpu == 0)
			{
				lastT     += timer.tock();
				lastSteps += 1;

				if (lastT >= printPerfInterval)
				{
					double stepTime = lastT / lastSteps;
					double completionTime = stepTime * (maxTS - current_timestep - ts);
					std::string units;
					if (completionTime > 60*60*24*365)
					{
						units = "weeks";
						completionTime /= 60*60*24*365;
					}
					else if (completionTime > 60*60*24*30)
					{
						units = "months";
						completionTime /= 60*60*24*30;
					}
					else if (completionTime > 60*60*24*7)
					{
						units = "weeks";
						completionTime /= 60*60*24*7;
					}
					else if (completionTime > 60*60*24)
					{
						units = "days";
						completionTime /= 60*60*24;
					}
					else if (completionTime > 60*60)
					{
						units = "hours";
						completionTime /= 60*60;
					}
					else if (completionTime > 60)
					{
						units = "minutes";
						completionTime /= 60;
					}
					else
					{
						units = "seconds";
					}

					Print::printf(Print::INFO, "Average walltime per timestep: %.2f ms. Progress: %.4fs/%.4fs (% .3g%% done / %.2g %s walltime remaining)",
								  1000.0*stepTime, (ts+current_timestep)*tau, maxTS*tau, 100.0*(ts+current_timestep)/maxTS, completionTime, units.c_str());

					lastT     = 0;
					lastSteps = 0;
				}

				timer.tick();
			}

			int overflows = run_next_timestep(gpu, ts + current_timestep);

			switch (overflow_handling)
			{
				case OVERFLOW_MODE_CLASSIC:
					if (use_spin_barrier)
						mclkr_barrier_spin(overflows);
					else
						mclkr_barrier_cond(overflows);

					if (mclkr_total_overflows > 0)
					{
						Print::printf(Print::DEBUG, "[GPU Thread %d] encountered %d overflows", gpu, overflows);
						handle_overflows(gpu, ((ByteLattice*) lattice)->getParticlesMemory(),
						                 p->dLattice, ts + current_timestep);
					}
					break;

				case OVERFLOW_MODE_RELAXED:
					// All overflows are handled by the gpu. No need to communicate.
					pthread_barrier_wait(&simulation_barrier);
					break;

				default:
					throw("Unimplemented overflow handler");
			}
		}

		mapper->stage_out(gpu, ((ByteLattice*) lattice)->getParticlesMemory(), p->dLattice);
		pthread_barrier_wait(&stop_barrier);

	} while(true);

	Print::printf(Print::DEBUG,"[GPU Thread %d] leaving.",gpu);

	return NULL;
}

int MGPUMpdRdmeSolver::hookSimulation(double time, ByteLattice *lattice)
{
	return 0;
}

void MGPUMpdRdmeSolver::calculateXLaunchParameters(dim3 * gridSize,
                                                   dim3 * threadBlockSize,
												   const unsigned int maxXBlockSize,
												   const unsigned int latticeXSize,
												   const unsigned int latticeYSize,
												   const unsigned int latticeZSize)
{
    unsigned int xBlockXSize = min(maxXBlockSize, latticeXSize);
    unsigned int gridXSize   = latticeXSize / xBlockXSize;
	
    if (gridXSize*xBlockXSize != latticeXSize)
	{
		// Find the largest number of warps that is divisible
		unsigned int tryx=32;
		while (tryx < maxXBlockSize)
		{
			if (latticeXSize % tryx == 0)
				xBlockXSize = tryx;
			tryx += 32;
		}
		gridXSize = latticeXSize / xBlockXSize;
	}
			
    (*gridSize).x = gridXSize;
    (*gridSize).y = latticeYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = xBlockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = 1;
}

void MGPUMpdRdmeSolver::calculateYLaunchParameters(dim3 * gridSize,
                                                   dim3 * threadBlockSize,
												   const unsigned int blockXSize,
												   const unsigned int blockYSize,
												   const unsigned int latticeXSize,
												   const unsigned int latticeYSize,
												   const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize/blockYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

void MGPUMpdRdmeSolver::calculateZLaunchParameters(dim3 * gridSize,
                                                   dim3 * threadBlockSize,
												   const unsigned int blockXSize,
												   const unsigned int blockZSize,
												   const unsigned int latticeXSize,
												   const unsigned int latticeYSize,
												   const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize;
    (*gridSize).z = latticeZSize/blockZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = blockZSize;
}

void MGPUMpdRdmeSolver::calculateReactionLaunchParameters(dim3 * gridSize,
                                                          dim3 * threadBlockSize,
														  const unsigned int blockXSize,
														  const unsigned int blockYSize,
														  const unsigned int latticeXSize,
														  const unsigned int latticeYSize,
														  const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize/blockYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

namespace mgpumpdrdme_dev {
/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void MGPU_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int z_start, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
	const unsigned int x = (blockDim.x * blockIdx.x) + threadIdx.x;
	const unsigned int y = (blockDim.y * blockIdx.y) + threadIdx.y;
	const unsigned int z = blockIdx.z+z_start;
	
    // Figure out the offset of this thread in the lattice and the lattice segment.
    const unsigned int latticeXIndex = x;
	unsigned int latticeIndex = local_index(x, y, z);
	unsigned int globalIndex = local_to_global(x, y, z); // Use the global lattice index for RNG
	const unsigned int windowIndex = threadIdx.x+MPD_APRON_SIZE;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment.
    __shared__ unsigned int window[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_X_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyXWindowFromLattice(blockIdx.x, inLattice, window, latticeIndex, latticeXIndex, windowIndex);
    copyXWindowFromSites(blockIdx.x, inSites, sitesWindow, latticeIndex, latticeXIndex, windowIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];

    // Make the choices.
    makeXDiffusionChoices(window, sitesWindow, choices, globalIndex, latticeXIndex, windowIndex, blockDim.x, timestepHash);
	__syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Propagate the choices to the new lattice segment.
    performPropagation(outLattice, window, choices, latticeIndex, windowIndex-1, windowIndex, windowIndex+1, MPD_X_WINDOW_SIZE, siteOverflowList);
}


/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void MGPU_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{

	const unsigned int x = (blockDim.x * blockIdx.x) + threadIdx.x;
	const unsigned int y = (blockDim.y * blockIdx.y) + threadIdx.y;
	const unsigned int z = blockIdx.z;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (blockIdx.y*blockDim.y) + threadIdx.y;
    unsigned int windowYIndex = threadIdx.y+MPD_APRON_SIZE;
    unsigned int windowIndex = (windowYIndex*blockDim.x) + threadIdx.x;

	unsigned int latticeIndex = local_index(x, y, z);
	unsigned int globalIndex = local_to_global(x, y, z); // Use the global lattice index for RNG

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[MPD_Y_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_Y_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyYWindowFromLattice(inLattice, window, latticeIndex, latticeYIndex, windowIndex, windowYIndex);
    copyYWindowFromSites(inSites, sitesWindow, latticeIndex, latticeYIndex, windowIndex, windowYIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_Y_WINDOW_SIZE*MPD_WORDS_PER_SITE];

    // Make the choices.
    makeYDiffusionChoices(window, sitesWindow, choices, globalIndex, latticeYIndex, windowIndex, windowYIndex, timestepHash);
	__syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagation(outLattice, window, choices, latticeIndex, windowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, siteOverflowList);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void MGPU_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start)
{
	const unsigned int x = (blockDim.x * blockIdx.x) + threadIdx.x;
	const unsigned int y = blockIdx.y;
	const unsigned int z = (blockDim.z * blockIdx.z) + threadIdx.z + z_start;

    // Figure out the offset of this thread in the lattice and the lattice segment.
	int global_z = z + goffsetZC;
    unsigned int windowZIndex = threadIdx.z+MPD_APRON_SIZE;
    unsigned int windowIndex = (windowZIndex*blockDim.x) + threadIdx.x;

	unsigned int latticeIndex = local_index(x, y, z);
	unsigned int globalIndex = local_to_global(x, y, z); // Use the global lattice index for RNG

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[MPD_Z_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_Z_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyZWindowFromLattice_MGPU(inLattice, window, latticeIndex, global_z, windowIndex, windowZIndex);
    copyZWindowFromSites_MGPU(inSites, sitesWindow, latticeIndex, global_z, windowIndex, windowZIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_Z_WINDOW_SIZE*MPD_WORDS_PER_SITE];

    // Make the choices.
    makeZDiffusionChoices(window, sitesWindow, choices, globalIndex, global_z, windowIndex, windowZIndex, timestepHash);
	__syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagation(outLattice, window, choices, latticeIndex, windowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, siteOverflowList);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void MGPU_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start,   const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
#else
__global__ void MGPU_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start,   const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG)
#endif
#else
__global__ void MGPU_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start)
#endif
{
	const unsigned int x = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int y = blockDim.y * blockIdx.y + threadIdx.y;
	const unsigned int z = blockIdx.z+z_start;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeIndex = local_index(x, y, z);
    unsigned int globalIndex = local_to_global(x, y, z); // Use the global lattice index for RNG

    ///////////////////////////////////////////
    // Load the particles and site.          //
    ///////////////////////////////////////////

    unsigned int particles[MPD_WORDS_PER_SITE];
    for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
        particles[w] = inLattice[latticeIndex+latticeOffset];
    uint8_t siteType = inSites[latticeIndex];

    ////////////////////////////////////////
    // Perform the reactions.             //
    ////////////////////////////////////////

    // Calculate the kinetic rate for each reaction at this site.
    float totalReactionPropensity = 0.0f;
    for (int i=0; i<numberReactionsC; i++)
    {
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        totalReactionPropensity += calculateReactionPropensity(siteType, (uint8_t*)particles, i, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        totalReactionPropensity += calculateReactionPropensity(siteType, (uint8_t*)particles, i, RLG);
#endif
#else
        totalReactionPropensity += calculateReactionPropensity(siteType, (uint8_t*)particles, i);
#endif
    }


    // See if a reaction occurred at the site.
    float reactionProbability = calculateReactionProbability(totalReactionPropensity);
    unsigned int reactionOccurred = checkForReaction(globalIndex, reactionProbability, timestepHash);

    // If there was a reaction, process it.
    if (reactionOccurred)
    {
        // Figure out which reaction occurred.
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash, RLG);
#endif
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash);
#endif

        // Construct the new site.
#ifdef MPD_GLOBAL_S_MATRIX
        evaluateReaction(latticeIndex, siteType, (uint8_t*)particles, reactionIndex, siteOverflowList, SG);
#else
        evaluateReaction(latticeIndex, siteType, (uint8_t*)particles, reactionIndex, siteOverflowList);
#endif

        // Copy the new particles back into the lattice.
        for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
             outLattice[latticeIndex+latticeOffset] = particles[w];
    }
}

#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) MGPU_precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) MGPU_precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2)
#endif
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) MGPU_precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2)
#endif
{
	const unsigned int x = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int y = blockDim.y * blockIdx.y + threadIdx.y;
	const unsigned int z = blockIdx.z+z_start;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeIndex = local_index(x, y, z);
    unsigned int globalIndex = local_to_global(x, y, z); // Use the global lattice index for RNG

    ///////////////////////////////////////////
    // Load the particles and site.          //
    ///////////////////////////////////////////

    unsigned int particles[MPD_WORDS_PER_SITE];
    for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
        particles[w] = inLattice[latticeIndex+latticeOffset];
    uint8_t siteType = inSites[latticeIndex];

    ////////////////////////////////////////
    // Perform the reactions.             //
    ////////////////////////////////////////

    // Calculate the kinetic rate for each reaction at this site.
    float totalReactionPropensity = read_element(qp0,siteType);
    //float totalReactionPropensity = qp0[siteType];
	for(uint i=0; i<MPD_PARTICLES_PER_SITE; i++)
	{
		uint8_t p1=((uint8_t*)particles)[i];
		if(p1 > 0)
		{
			totalReactionPropensity += read_element(qp1,(siteType * numberSpeciesC + (p1-1)));
			//totalReactionPropensity += qp1[(siteType * numberSpeciesC + (p1-1))];
			for(uint j=i+1; j<MPD_PARTICLES_PER_SITE; j++)
			{
				uint8_t p2=((uint8_t*)particles)[j];
				if(p2 > 0)
				{
					totalReactionPropensity += read_element(qp2,siteType*numberSpeciesC*numberSpeciesC + (p1-1)*numberSpeciesC + (p2-1));
					//totalReactionPropensity += qp2[siteType*numberSpeciesC*numberSpeciesC + (p1-1)*numberSpeciesC + (p2-1)];
				}
			}
		}
	}


	// If propensity is zero, no reaction can occur.
	if(totalReactionPropensity == 0.0f)
		return;

    // See if a reaction occurred at the site.
    float reactionProbability = calculateReactionProbability(totalReactionPropensity);
    unsigned int reactionOccurred = checkForReaction(globalIndex, reactionProbability, timestepHash);

    // If there was a reaction, process it.
    if (reactionOccurred)
    {
        // Figure out which reaction occurred.
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash, RLG);
#endif
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash);
#endif

        // Construct the new site.
#ifdef MPD_GLOBAL_S_MATRIX
        evaluateReaction(latticeIndex, siteType, (uint8_t*)particles, reactionIndex, siteOverflowList, SG);
#else
        evaluateReaction(latticeIndex, siteType, (uint8_t*)particles, reactionIndex, siteOverflowList);
#endif

        // Copy the new particles back into the lattice.
        for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
             outLattice[latticeIndex+latticeOffset] = particles[w];
    }
}

__global__ void correct_overflows_mgpu(unsigned int* lattice, unsigned int* siteOverflowList)
{
	// Remember: #define MPD_OVERFLOW_LIST_SIZE       1+2*TUNE_MPD_MAX_PARTICLE_OVERFLOWS*sizeof(uint32_t)
	// Requirement: one block with TUNE_MPD_MAX_PARTICLE_OVERFLOWS threads
	
	const unsigned int i = threadIdx.x;
	const unsigned int total = siteOverflowList[0];

	__shared__ unsigned int indexes[TUNE_MPD_MAX_PARTICLE_OVERFLOWS];
	__shared__ unsigned int maxround;

	// Do I have an overflow to look at?
	if(i >= total)
		return;

	// Abort if overflows have overflown
	if(threadIdx.x == 0) assert(total < TUNE_MPD_MAX_PARTICLE_OVERFLOWS);

	// load our index
	lattice_size_t latticeIndex = siteOverflowList[(i*2)+1];
	particle_t particle = siteOverflowList[(i*2)+2];

	indexes[i] = latticeIndex;

	// zero out list
	maxround=0;
	__syncthreads();
	siteOverflowList[0]=0;

	// Discover which round I should go in.  To prevent situations where two threads
	// will try and add a particle to the same site at the same time, and thus
	// negating the others, each thread determines what round it is allowed to make
	// edits in, by determining how many previous overflows are occuring at the 
	// same index.
	int round=0;
	for(int j=0; j < i; j++)
	{
		if(indexes[j] == latticeIndex)
			round++;
	}
	atomicMax(&maxround, round);
	__syncthreads();

	for(int r=0; r <= maxround; r++)
	{
		if(round == r)
		{
			unsigned int particles[MPD_WORDS_PER_SITE];
			for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
				particles[w] = lattice[latticeIndex+latticeOffset];
			
			uint8_t *p=(uint8_t*)particles;
			int ok=0;
			for(int pi=0; pi<MPD_PARTICLES_PER_SITE; pi++)
			{
				if(p[pi] == 0)
				{
					p[pi]=particle;
					ok=1;
					//printf("(round %d) Corrected overflow of particle %d at index %d\n", r, particle, latticeIndex);
					break;
				}
			}

			if(ok)
			{
				// Copy the new particles back into the lattice.
				for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
				     lattice[latticeIndex+latticeOffset] = particles[w];
			}
			else
			{
				int exceptionIndex = atomicAdd(siteOverflowList, 1);
				siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
				siteOverflowList[(exceptionIndex*2)+2]=particle;
				//printf("(round %d) Failed to correct overflow of particle %d at index %d\n", r, particle, latticeIndex);
			}
		}
		__syncthreads();
	}

	//if(i  == 0)
		//printf("%d) in: %d overflows, out: %d\n", gpuidC, total, siteOverflowList[0]);

}

#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) MGPU_precomp_reaction_kernel_packing(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, unsigned int* buf_top, unsigned int* buf_bot)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) MGPU_precomp_reaction_kernel_packing(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, unsigned int* buf_top, unsigned int* buf_bot)
#endif
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) MGPU_precomp_reaction_kernel_packing(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start,  const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, unsigned int* buf_top, unsigned int* buf_bot)
#endif
{
	const unsigned int x = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int y = blockDim.y * blockIdx.y + threadIdx.y;
	const unsigned int z = blockIdx.z+z_start;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeIndex = local_index(x, y, z);
    unsigned int globalIndex = local_to_global(x, y, z); // Use the global lattice index for RNG

    ///////////////////////////////////////////
    // Load the particles and site.          //
    ///////////////////////////////////////////

    unsigned int particles[MPD_WORDS_PER_SITE];
    for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
        particles[w] = inLattice[latticeIndex+latticeOffset];
    uint8_t siteType = inSites[latticeIndex];

    ////////////////////////////////////////
    // Perform the reactions.             //
    ////////////////////////////////////////

    // Calculate the kinetic rate for each reaction at this site.
    float totalReactionPropensity = read_element(qp0,siteType);
    //float totalReactionPropensity = qp0[siteType];
	for(uint i=0; i<MPD_PARTICLES_PER_SITE; i++)
	{
		uint8_t p1=((uint8_t*)particles)[i];
		if(p1 > 0)
		{
			totalReactionPropensity += read_element(qp1,(siteType * numberSpeciesC + (p1-1)));
			//totalReactionPropensity += qp1[(siteType * numberSpeciesC + (p1-1))];
			for(uint j=i+1; j<MPD_PARTICLES_PER_SITE; j++)
			{
				uint8_t p2=((uint8_t*)particles)[j];
				if(p2 > 0)
				{
					totalReactionPropensity += read_element(qp2,siteType*numberSpeciesC*numberSpeciesC + (p1-1)*numberSpeciesC + (p2-1));
					//totalReactionPropensity += qp2[siteType*numberSpeciesC*numberSpeciesC + (p1-1)*numberSpeciesC + (p2-1)];
				}
			}
		}
	}

	float reactionProbability;
	unsigned int reactionOccurred;

	// If propensity is zero, no reaction can occur.
	if(totalReactionPropensity == 0.0f)
		goto letspackitup;

    // See if a reaction occurred at the site.
    reactionProbability = calculateReactionProbability(totalReactionPropensity);
    reactionOccurred = checkForReaction(globalIndex, reactionProbability, timestepHash);

    // If there was a reaction, process it.
    if (reactionOccurred)
    {
        // Figure out which reaction occurred.
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash, RLG);
#endif
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash);
#endif

        // Construct the new site.
#ifdef MPD_GLOBAL_S_MATRIX
        evaluateReaction(latticeIndex, siteType, (uint8_t*)particles, reactionIndex, siteOverflowList, SG);
#else
        evaluateReaction(latticeIndex, siteType, (uint8_t*)particles, reactionIndex, siteOverflowList);
#endif

        // Copy the new particles back into the lattice.
        for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
             outLattice[latticeIndex+latticeOffset] = particles[w];
    }

letspackitup:

	bool in_top_send = (z < (top_send_z + top_dim_z));
	bool in_bot_send = (z >= bot_send_z);

/*
	if(x == 0 && y == 0 && in_top_send)
		printf("%d, %d, %d in top send send z =%d, dim z = %d, sz=%d\n", x, y, z, top_send_z, top_dim_z, top_dim_size);
	if(x == 0 && y == 0 && in_bot_send)
		printf("%d, %d, %d in bot send send z =%d, dim z = %d, sz=%d\n", x, y, z, bot_send_z, bot_dim_z, bot_dim_size);
*/

	if(in_top_send && buf_top != NULL)
	{
		unsigned int packidx = x + y*latticeXSizeC + (z-top_send_z)*latticeXYSizeC;
        for (uint w=0, packOffset=0; w<MPD_WORDS_PER_SITE; w++, packOffset+=top_dim_size)
             buf_top[packidx+packOffset] = particles[w];
	}

	if(in_bot_send && buf_bot != NULL)
	{
		unsigned int packidx = x + y*latticeXSizeC + (z-bot_send_z)*latticeXYSizeC;
        for (uint w=0, packOffset=0; w<MPD_WORDS_PER_SITE; w++, packOffset+=bot_dim_size)
             buf_bot[packidx+packOffset] = particles[w];
	}
}

__global__ void MGPU_x_kernel_unpack(const unsigned int* inLattice,
                                     const uint8_t * inSites,
									 unsigned int* outLattice,
									 const unsigned int z_start,
									 const unsigned long long timestepHash,
									 unsigned int* siteOverflowList,
									 unsigned int* buf_top,
									 unsigned int *buf_bot)
{
	const unsigned int x = (blockDim.x * blockIdx.x) + threadIdx.x;
	const unsigned int y = (blockDim.y * blockIdx.y) + threadIdx.y;
	const unsigned int z = blockIdx.z+z_start;
	
    // Figure out the offset of this thread in the lattice and the lattice segment.
    const unsigned int latticeXIndex = x;
	unsigned int latticeIndex = local_index(x, y, z);
	unsigned int globalIndex = local_to_global(x, y, z); // Use the global lattice index for RNG
	const unsigned int windowIndex = threadIdx.x+MPD_APRON_SIZE;

	// Check if I need to read from buffer and not lattice
	bool in_top_recv = (z < (top_recv_z + top_dim_z));
	bool in_bot_recv = (z >= bot_recv_z);

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment.
    __shared__ unsigned int window[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_X_WINDOW_SIZE];

	if(in_top_recv && buf_top != NULL)
	{
		unsigned int packidx = x + y*latticeXSizeC + (z-top_recv_z)*latticeXYSizeC;
		copyXWindowFromBuffer(blockIdx.x, buf_top, window, packidx, latticeXIndex, windowIndex, top_dim_size);
	}
	else if(in_bot_recv && buf_bot != NULL)
	{
		unsigned int packidx = x + y*latticeXSizeC + (z-bot_recv_z)*latticeXYSizeC;
		copyXWindowFromBuffer(blockIdx.x, buf_bot, window, packidx, latticeXIndex, windowIndex, bot_dim_size);
	}
	else
	{
		// Copy the x window from device memory into shared memory.
		copyXWindowFromLattice(blockIdx.x, inLattice, window, latticeIndex, latticeXIndex, windowIndex);
	}

    copyXWindowFromSites(blockIdx.x, inSites, sitesWindow, latticeIndex, latticeXIndex, windowIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];

    // Make the choices.
    makeXDiffusionChoices(window, sitesWindow, choices, globalIndex, latticeXIndex, windowIndex, blockDim.x, timestepHash);
	__syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Propagate the choices to the new lattice segment.
    performPropagation(outLattice, window, choices, latticeIndex, windowIndex-1, windowIndex, windowIndex+1, MPD_X_WINDOW_SIZE, siteOverflowList);
}

// For periodic lattices with a potentially negative offset
__device__ inline size_t local_to_global(unsigned int x, unsigned int y, unsigned int z)
{
    ssize_t ix=(int)x+goffsetXC;
    ix += ((int)y+goffsetYC)*(int)latticeXSizeC;
    ix += ((int)z+goffsetZC)*(int)latticeXYSizeC;
    size_t max_index = global_latticeXYZSizeC;
    if(ix < 0)
        return ix+max_index;
    else if(ix >= max_index)
        return ix-max_index;
    else
        return ix;
}
/*
// For non-periodic lattices with a strictly positive offset
__device__ inline size_t local_to_global(unsigned int x, unsigned int y, unsigned int z)
{
    size_t ix=x+goffsetXC;
    ix += (y+goffsetYC)*latticeXSizeC;
    ix += (z+goffsetZC)*latticeXYSizeC;
	return ix;
}
*/

__device__ inline size_t local_index(unsigned int x, unsigned int y, unsigned int z)
{
	return x + y*latticeXSizeC + z*latticeXYSizeC;
}

}
}
}
