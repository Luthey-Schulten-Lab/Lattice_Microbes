/*
 * University of Illinois Open Source License
 * Copyright 2014-2018 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.uiuc.edu/~schulten
 * 
 * Overflow algorithm in RDME solvers and CPU assignment (2012)
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
 * Author(s): Mike Hallock
 */

#include <map>
#include <string>
#include <cstdlib>
#include "config.h"
#if defined(MACOSX)
#include <mach/mach_time.h>
#elif defined(LINUX)
#include <time.h>
#endif
#include "cuda/lm_cuda.h"
#include "core/Math.h"
#include "core/Print.h"
#include "cme/CMESolver.h"
#include "DiffusionModel.pb.h"
#include "Lattice.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "core/ResourceAllocator.h"
#include "rdme/ByteLattice.h"
#include "rdme/MPIMpdRdmeSolver.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include "core/Timer.h"

#define MPD_WORDS_PER_SITE              (MPD_LATTICE_MAX_OCCUPANCY / 4)
#define MPD_APRON_SIZE                  1

#include "cuda/constant.cuh"

namespace lm {
namespace rdme {
namespace mpimpdrdme_dev {
#include "rdme/dev/xor_random_dev.cu"
#include "rdme/dev/lattice_sim_1d_dev.cu"
#include "rdme/dev/byte_diffusion_1d_dev.cu"
#include "rdme/dev/byte_reaction_dev.cu"
}}}

#include "rdme/GPUMapper/SegmentDescriptor.h"
#include "rdme/GPUMapper/ZDivMPIGPUMapper.h"
#include "rdme/GPUMapper/MultiGPUMapper.h"

#include "mpi.h"
#include "mpi/lm_mpi.h"

#define L3INDEX(_x,_y,_z, _g) ((_x)+((_y)*(_g.x))+((_z)*(_g.x)*(_g.y)))

extern int cudaDevicesPerNode;


using std::map;
using lm::io::DiffusionModel;
using lm::rdme::Lattice;
using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {


MPIMpdRdmeSolver::MPIMpdRdmeSolver()
:RDMESolver(lm::rng::RandomGenerator::NONE),seed(0),cudaOverflowList(NULL),tau(0.0),overflowTimesteps(0),overflowListUses(0)
{
}

MPIMpdRdmeSolver::~MPIMpdRdmeSolver()
{
}

void MPIMpdRdmeSolver::initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resourcesA)
{
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

	// Distribute seed to all ranks
	MPI_Bcast(&seedTop, 1, MPI_INT,  lm::MPI::MASTER, MPI_COMM_WORLD);
    }
    seed = (seedTop<<16)|(replicate&0x0000FFFF);
    Print::printf(Print::INFO, "MPDRDME: Rng seed: top word %u, bottom word %u", (seed>>16)*0xffff, seed&0xFFFF);

	resources=resourcesA;
}

// TODO make collective
void MPIMpdRdmeSolver::initialize_decomposition()
{
    if(MPI_Comm_size(MPI_COMM_WORLD, &world_size)!=MPI_SUCCESS)
    {
        printf("Failed to get comm size\n");
        exit(1);
    }

    if(MPI_Comm_rank(MPI_COMM_WORLD, &rank)!=MPI_SUCCESS)
    {
        printf("Failed to get rank\n");
        exit(1);
    }

	// TODO: This assumes we will use full nodes only and that ranks fill-up
	int gpu = rank % cudaDevicesPerNode;

	// Create mapper
	int latx,laty,latz;
	int apron,overlap;
	char *ids=NULL;
#if defined MPD_BOUNDARY_PERIODIC
	bool pz=true;
#else
	bool pz=false;
#endif
	apron=2;
	overlap=0;

	latx=lattice->getXSize();
	laty=lattice->getYSize();
	latz=lattice->getZSize();

	int site_size = sizeof(int);
	int pagecount = MPD_WORDS_PER_SITE;

	Print::printf(Print::DEBUG,"rank %d Create mapper: X=%d Y=%d Z=%d, sitesize=%d, apron=%d, overlap=%d, gpu=%d, gpu_ids=%s, periodic_z=%d, pagecount=%d\n",
			rank, latx, laty, latz,
			site_size, apron, overlap,
			gpu, ids, pz, pagecount);
	mapper=new ZDivMPIGPUMapper(latx, laty, latz, site_size, apron, overlap, gpu, pz, pagecount);
}

void MPIMpdRdmeSolver::allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing)
{
// TODO no allocate whole lattice!
    lattice = (Lattice *)new ByteLattice(latticeXSize, latticeYSize, latticeZSize, latticeSpacing, particlesPerSite);
}

void MPIMpdRdmeSolver::buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionTypesA, const double * KA, const int * SA, const uint * DA, const uint kCols)
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

void MPIMpdRdmeSolver::buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData)
{
    if(bytes_per_particle != 1) throw Exception("MPIMpdRdmeSolver only supports 1 byte per particle at this time.");

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
}

void MPIMpdRdmeSolver::copyModelsToDevice()
{
    const size_t DFmatrixSize = numberSpecies*numberSiteTypes*numberSiteTypes;
	const unsigned int latticeXSize=segment->local_dimensions.x;
	const unsigned int latticeYSize=segment->local_dimensions.y;
	const unsigned int latticeZSize=segment->local_dimensions.z;
	
	// the lattice{X,XY,XYZ}SizeC refer to sizes of the LOCAL lattice.
	// NOTE that since we have 1D diffusion only, the global X and XY sizes are the same.
    const unsigned int latticeXYSize = latticeXSize*latticeYSize;
    const unsigned int latticeXYZSize = latticeXSize*latticeYSize*latticeZSize;

	// gloal_latticeXYZSizeC is the full lattice XYZ size.  
    const unsigned int global_latticeZSize = lattice->getZSize();
    const unsigned int global_latticeXYZSize = latticeXSize*latticeYSize*lattice->getZSize();

    // Copy the reaction model and S matrix to constant memory on the GPU.
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberReactionsC, &numberReactions, sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionOrdersC, model_reactionOrders, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionSitesC, model_reactionSites, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(D1C, model_D1, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(D2C, model_D2, numberReactions*sizeof(unsigned int)));
#ifdef MPD_GLOBAL_S_MATRIX
	// If S matrix stored in global memory, allocate space and perform copy
	cudaMalloc(&SG, numberSpecies*numberReactions * sizeof(int8_t));
	cudaMemcpy(SG, model_S, numberSpecies*numberReactions * sizeof(int8_t), cudaMemcpyHostToDevice);
#else
	// S matrix is in constant memory
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(SC, model_S, numberSpecies*numberReactions*sizeof(int8_t)));
#endif

    // Copy the diffusion model to constant memory on the GPU.
#ifdef MPD_GLOBAL_T_MATRIX
	CUDA_EXCEPTION_CHECK(cudaMalloc(&TG, DFmatrixSize*sizeof(float)));
	CUDA_EXCEPTION_CHECK(cudaMemcpy(TG, model_T, DFmatrixSize*sizeof(float), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(TC, &TG, sizeof(float*)));
#else
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(TC, model_T, DFmatrixSize*sizeof(float)));
#endif

    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberSpeciesC, &numberSpecies, sizeof(numberSpeciesC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberSiteTypesC, &numberSiteTypes, sizeof(numberSiteTypesC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXSizeC, &latticeXSize, sizeof(latticeYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeYSizeC, &latticeYSize, sizeof(latticeYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeZSizeC, &latticeZSize, sizeof(latticeZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXYSizeC, &latticeXYSize, sizeof(latticeXYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXYZSizeC, &latticeXYZSize, sizeof(latticeXYZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(global_latticeZSizeC, &global_latticeZSize, sizeof(global_latticeZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(global_latticeXYZSizeC, &global_latticeXYZSize, sizeof(global_latticeXYZSizeC)));
#ifdef MPD_GLOBAL_S_MATRIX
	// Store RL in global memory too, since I'm going to assume if S is too big, then RL is too.
	cudaMalloc(&RLG, numberReactions*numberSiteTypes * sizeof(uint8_t));
	cudaMemcpy(RLG, model_RL, numberReactions*numberSiteTypes * sizeof(uint8_t), cudaMemcpyHostToDevice);
#else
	// RL is stored in constant memory
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(RLC, model_RL, numberReactions*numberSiteTypes*sizeof(uint8_t)));
#endif
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionRatesC, model_reactionRates, numberReactions*sizeof(float)));

	// Copy per-GPU global offsets
    //CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(lsizeXC, &(threads[gpu].local_dimensions.x), sizeof(unsigned int)));
    //CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(lsizeYC, &(threads[gpu].local_dimensions.y), sizeof(unsigned int)));
    //CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(lsizeZC, &(threads[gpu].local_dimensions.z), sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(goffsetXC, &(segment->global_offset.x), sizeof(int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(goffsetYC, &(segment->global_offset.y), sizeof(int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(goffsetZC, &(segment->global_offset.z), sizeof(int)));

//    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(gpuidC, &gpu, sizeof(unsigned int)));
	
}

void MPIMpdRdmeSolver::generateTrajectory()
{
	// Determine how to divide the lattice
	initialize_decomposition();

	//TODO mapper->numa_bind_thread(gpu);
	//TODO mapper->initialize_gpu(gpu);
	//TODO mapper->set_affinity(gpu);
	segment=mapper->initialize();
	mapper->initialize_gpu();

	prepare_gpu();

    // Shadow the lattice member as a cuda lattice.
    ByteLattice * lattice = (ByteLattice *)this->lattice;

    // Get the interval for writing species counts and lattices.
    double speciesCountsWriteInterval=atof((*parameters)["writeInterval"].c_str());
    double nextSpeciesCountsWriteTime = speciesCountsWriteInterval;
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpeciesToTrack);
    speciesCountsDataSet.set_number_entries(0);
    double latticeWriteInterval=atof((*parameters)["latticeWriteInterval"].c_str());
    double nextLatticeWriteTime = latticeWriteInterval;
    lm::io::Lattice latticeDataSet;

    // Get the simulation time limit.
    double maxTime=atof((*parameters)["maxTime"].c_str());

    Print::printf(Print::DEBUG, "Running mpd rdme simulation with %d species, %d reactions, %d site types for %e s with tau %e. Writing species at %e and lattice at %e intervals", numberSpecies, numberReactions, numberSiteTypes, maxTime, tau, speciesCountsWriteInterval, latticeWriteInterval);

    // Set the initial time.
    double time = 0.0;
	absolute_timestep=1;

    // Record the initial species counts.
    recordSpeciesCounts(time, lattice, &speciesCountsDataSet);

    // Write the initial lattice.
    writeLatticeData(time, lattice, &latticeDataSet);

	// Copy the parameter matricies to the device
	copyModelsToDevice();
	mapper->stage_in(dLattice, ((ByteLattice*) lattice)->getParticlesMemory());
	mapper->stage_in_sites(dSites, ((ByteLattice*) lattice)->getSitesMemory());

	// Make sure everyone is ready!
	MPI_Barrier(MPI_COMM_WORLD);
	mapper->communicate_edges(dLattice, absolute_timestep+1);

	// prime for correctness
    mapper->copy_edge_to_host(dLattice, absolute_timestep, Z_SHADOW_TOP, 0);
    mapper->copy_edge_to_host(dLattice, absolute_timestep, Z_SHADOW_BOTTOM, 0);

    Timer timer;
    timer.tick();
    double lastT=0;
    int lastSteps=0;

    // Loop until we have finished the simulation.
    while (time < maxTime)
    {
        lastT += timer.tock();
        lastSteps += 1;

        if ( lastT >= 5)
        {
            double stepTime = lastT/lastSteps;
            double completionTime = stepTime*(maxTime-time)/tau;
            std::string units;
            if (completionTime > 60*60*24*365) {
                units = "weeks";
                completionTime /= 60*60*24*365;
            } else if (completionTime > 60*60*24*30) {
                units = "months";
                completionTime /= 60*60*24*30;
            } else if (completionTime > 60*60*24*7) {
                units = "weeks";
                completionTime /= 60*60*24*7;
            } else if (completionTime > 60*60*24) {
                units = "days";
                completionTime /= 60*60*24;
            } else if (completionTime > 60*60) {
                units = "hours";
                completionTime /= 60*60;
            } else if (completionTime > 60) {
                units = "minutes";
                completionTime /= 60;
            } else {
                units = "seconds";
            }

            if(rank == 0) Print::printf(Print::INFO, "Average walltime per timestep: %.2f ms. Progress: %.4fs/%.4fs (% .3g%% done / %.2g %s walltime remaining)",
                                       1000.0*stepTime, time, maxTime, 100.0*time/maxTime, completionTime, units.c_str());

            lastT = 0;
            lastSteps = 0;
        }

        timer.tick();
		run_next_timestep(absolute_timestep);

		absolute_timestep ++;
        time += tau; 

        // See if we need to write out the any data.
        if (time >= nextLatticeWriteTime-EPS || time >= nextSpeciesCountsWriteTime-EPS)
        {
            // Synchronize the lattice.
			mapper->stage_out(((ByteLattice*) lattice)->getParticlesMemory(), dLattice);

            // See if we need to write the lattice.
            if (time >= nextLatticeWriteTime-EPS)
            {
                PROF_BEGIN(PROF_SERIALIZE_LATTICE);
                writeLatticeData(time, lattice, &latticeDataSet);
                nextLatticeWriteTime += latticeWriteInterval;
                PROF_END(PROF_SERIALIZE_LATTICE);
            }

            // See if we need to write the species counts.
            if (time >= nextSpeciesCountsWriteTime-EPS)
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
    }

    // Write any remaining species counts.
    writeSpeciesCounts(&speciesCountsDataSet);
}

void MPIMpdRdmeSolver::writeLatticeData(double time, ByteLattice * lattice, lm::io::Lattice * latticeDataSet)
{

	if (lm::MPI::worldRank == lm::MPI::MASTER)
	{
		// Collect lattice data from slave ranks
		for(int i=1; i < world_size; i++)
		{
	//		printf("Staging in from rank %d\n", i);
			mapper->stage_in_from_slave(i,((ByteLattice*) lattice)->getParticlesMemory());
		}
	}
	else
	{
	//	printf("rank %d staging out\n",rank);
		mapper->stage_out_to_master(((ByteLattice*) lattice)->getParticlesMemory());
		return;
	}

    Print::printf(Print::DEBUG, "Writing lattice at %e s", time);

    // Record the lattice data.
    latticeDataSet->Clear();
    latticeDataSet->set_lattice_x_size(lattice->getSize().x);
    latticeDataSet->set_lattice_y_size(lattice->getSize().y);
    latticeDataSet->set_lattice_z_size(lattice->getSize().z);
    latticeDataSet->set_particles_per_site(lattice->getMaxOccupancy());
    latticeDataSet->set_time(time);

    // Push it to the output queue.
    size_t payloadSize = lattice->getSize().x*lattice->getSize().y*lattice->getSize().z*lattice->getMaxOccupancy()*sizeof(uint8_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::BYTE_LATTICE, replicate, latticeDataSet, lattice, payloadSize, &lm::rdme::ByteLattice::nativeSerialize);
}

void MPIMpdRdmeSolver::recordSpeciesCounts(double time, ByteLattice * lattice, lm::io::SpeciesCounts * speciesCountsDataSet)
{
	if (lm::MPI::worldRank == lm::MPI::MASTER)
	{
		std::map<particle_t,uint> particleCounts = lattice->getParticleCounts();
		speciesCountsDataSet->set_number_entries(speciesCountsDataSet->number_entries()+1);
		speciesCountsDataSet->add_time(time);
		for (particle_t p=0; p<numberSpeciesToTrack; p++)
		{
			speciesCountsDataSet->add_species_count((particleCounts.count(p+1)>0)?particleCounts[p+1]:0);
		}
	}
}

void MPIMpdRdmeSolver::writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet)
{
    if ( (lm::MPI::worldRank == lm::MPI::MASTER) && speciesCountsDataSet->number_entries() > 0)
    {
        // Push it to the output queue.
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, speciesCountsDataSet);

        // Reset the data set.
        speciesCountsDataSet->Clear();
        speciesCountsDataSet->set_number_species(numberSpeciesToTrack);
        speciesCountsDataSet->set_number_entries(0);
    }
}

uint64_t MPIMpdRdmeSolver::getTimestepSeed(uint32_t timestep, uint32_t substep)
{
    uint64_t timestepHash = (((((uint64_t)seed)<<30)+timestep)<<2)+substep;
    timestepHash = timestepHash * 3202034522624059733ULL + 4354685564936845319ULL;
    timestepHash ^= timestepHash >> 20; timestepHash ^= timestepHash << 41; timestepHash ^= timestepHash >> 5;
    timestepHash *= 7664345821815920749ULL;
    return timestepHash;
}

int MPIMpdRdmeSolver::run_next_timestep(uint32_t timestep)
{
    PROF_BEGIN(PROF_MPD_TIMESTEP);

	unsigned int *dtmp=dLatticeTmp;
	cudaStream_t cudaStream=stream1;
	int z_start=0;
	
    dim3 gridSize=grid_x;
	dim3 threadBlockSize=threads_x;

    mapper->send_edge(timestep+1, Z_SHADOW_TOP);
    mapper->send_edge(timestep+1, Z_SHADOW_BOTTOM);
    mapper->recv_edge(timestep, Z_SHADOW_BOTTOM);
    mapper->recv_edge(timestep, Z_SHADOW_TOP);
    mapper->copy_edge_to_device(dLattice, timestep, Z_SHADOW_TOP, cudaStream);
    mapper->copy_edge_to_device(dLattice, timestep, Z_SHADOW_BOTTOM, cudaStream);

    // Execute the kernel for the x direction.
    PROF_CUDA_START(cudaStream);

	// TODO: edge splitting for overlap
    PROF_CUDA_BEGIN(PROF_MPD_X_DIFFUSION,cudaStream);
    CUDA_EXCEPTION_EXECUTE((
		mpimpdrdme_dev::MPI_x_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>(dLattice, dSites, dtmp, z_start,
			getTimestepSeed(timestep,0), d_overflows)
	));
    PROF_CUDA_END(PROF_MPD_X_DIFFUSION,cudaStream);

    gridSize=grid_y;
	threadBlockSize=threads_y;

    // Execute the kernel for the y direction.
    PROF_CUDA_BEGIN(PROF_MPD_Y_DIFFUSION,cudaStream);
    CUDA_EXCEPTION_EXECUTE((mpimpdrdme_dev::MPI_y_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>(
		dtmp, dSites, dLattice,
		getTimestepSeed(timestep,1),
		d_overflows)));
    PROF_CUDA_END(PROF_MPD_Y_DIFFUSION,cudaStream);

    gridSize=grid_z;
	threadBlockSize=threads_z;

    // Execute the kernel for the z direction.
    PROF_CUDA_BEGIN(PROF_MPD_Z_DIFFUSION,cudaStream);
    CUDA_EXCEPTION_EXECUTE((
	mpimpdrdme_dev::MPI_z_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>
	(
		dLattice, dSites, dtmp,
		getTimestepSeed(timestep,2),
		d_overflows,
		segment->active_offset.z
	)));
    PROF_CUDA_END(PROF_MPD_Z_DIFFUSION,cudaStream);
	
    if (numberReactions > 0)
    {
		gridSize=grid_r;
		threadBlockSize=threads_r;
        // Execute the kernel for the reaction, this kernel updates the lattice in-place, so only the src pointer is passed.
        PROF_CUDA_BEGIN(PROF_MPD_REACTION,cudaStream);
        CUDA_EXCEPTION_EXECUTE((
		mpimpdrdme_dev::MPI_reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>(
			dtmp, dSites, dtmp,
			getTimestepSeed(timestep,3),
			d_overflows,
			segment->active_offset.z
#ifdef MPD_GLOBAL_S_MATRIX
			, SG, RLG
#endif
		)));
        PROF_CUDA_END(PROF_MPD_REACTION,cudaStream);
    }

    mpimpdrdme_dev::mpi_correct_overflows<<<dim3(1,1,1), dim3(TUNE_MPD_MAX_PARTICLE_OVERFLOWS,1,1),0,cudaStream>>>(
		dtmp, d_overflows);

	// send data
    mapper->copy_edge_to_host(dtmp, timestep, Z_SHADOW_TOP, cudaStream);
    mapper->copy_edge_to_host(dtmp, timestep, Z_SHADOW_BOTTOM, cudaStream);
	
    // Wait for the kernels to complete.
    PROF_BEGIN(PROF_MPD_SYNCHRONIZE);
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(cudaStream));
    PROF_END(PROF_MPD_SYNCHRONIZE);
	
	uint unresolved = ((unsigned int*)d_overflows)[0];
	if (unresolved > 0)
	{
		Print::printf(Print::WARNING, "[rank %d] %d unresolved overflows", rank, unresolved);
	}
	

	// Save correct pointers for current data
	dLatticeTmp=dLattice;
	dLattice=dtmp;

	return 0;
}



void MPIMpdRdmeSolver::prepare_gpu()
{
	// allocate device memory for lattice stuctures
	cudaMalloc(&dLattice, mapper->get_local_size()*MPD_WORDS_PER_SITE);
	cudaMalloc(&dLatticeTmp, mapper->get_local_size()*MPD_WORDS_PER_SITE);
	cudaMalloc(&dSites, DIMSIZE(mapper->get_local_dim()));

	// Allocate zero copy memory for overflow buffer
	CUDA_EXCEPTION_CHECK(cudaHostAlloc(&d_overflows,
        MPD_OVERFLOW_LIST_SIZE,
        cudaHostAllocPortable|cudaHostAllocMapped));
	memset(d_overflows, 0, MPD_OVERFLOW_LIST_SIZE);

    // Create data exchange streams and events
    CUDA_EXCEPTION_CHECK(cudaStreamCreate(&stream1));
    CUDA_EXCEPTION_CHECK(cudaStreamCreate(&stream2));
#ifdef PROF_USE_CUEVENT
    CUDA_EXCEPTION_CHECK(cudaEventCreate(&x_finish));
    CUDA_EXCEPTION_CHECK(cudaEventCreate(&diffusion_finished));
    CUDA_EXCEPTION_CHECK(cudaEventCreate(&rx_finish));
#else
    CUDA_EXCEPTION_CHECK(cudaEventCreateWithFlags(&x_finish,
		cudaEventDisableTiming));
    CUDA_EXCEPTION_CHECK(cudaEventCreateWithFlags(&diffusion_finished,
		cudaEventDisableTiming));
    CUDA_EXCEPTION_CHECK(cudaEventCreateWithFlags(&rx_finish,
		cudaEventDisableTiming));
#endif

	dim3 ldim=mapper->get_local_dim();
	dim3 gdim=mapper->get_global_dim();
	
	// calculate launch params
    calculateXLaunchParameters(&grid_x, &threads_x,
		TUNE_MPD_X_BLOCK_MAX_X_SIZE,
		ldim.x, ldim.y, ldim.z);
    calculateYLaunchParameters(&grid_y, &threads_y,
		TUNE_MPD_Y_BLOCK_X_SIZE, TUNE_MPD_Y_BLOCK_Y_SIZE, 
		ldim.x, ldim.y, ldim.z);
    calculateZLaunchParameters(&grid_z, &threads_z,
		TUNE_MPD_Z_BLOCK_X_SIZE, TUNE_MPD_Z_BLOCK_Z_SIZE, 
		ldim.x, ldim.y, gdim.z);
	calculateReactionLaunchParameters(&grid_r, &threads_r,
		TUNE_MPD_REACTION_BLOCK_X_SIZE, TUNE_MPD_REACTION_BLOCK_Y_SIZE,
		ldim.x, ldim.y, gdim.z);
}


/**
 * Gets the launch parameters for launching an x diffusion kernel.
 */
void MPIMpdRdmeSolver::calculateXLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    unsigned int xBlockXSize = min(maxXBlockSize,latticeXSize);
    unsigned int gridXSize = latticeXSize/xBlockXSize;
    if (gridXSize*xBlockXSize != latticeXSize)
	{
		// Find the largest number of warps that is divisible
		unsigned int tryx=32;
		while(tryx < maxXBlockSize)
		{
			if (latticeXSize % tryx == 0)
				xBlockXSize = tryx;
			tryx +=32;
		}
		gridXSize = latticeXSize/xBlockXSize;
	}
    (*gridSize).x = gridXSize;
    (*gridSize).y = latticeYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = xBlockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a y diffusion kernel.
 */
void MPIMpdRdmeSolver::calculateYLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize/blockYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a z diffusion kernel.
 */
void MPIMpdRdmeSolver::calculateZLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize;
    (*gridSize).z = latticeZSize/blockZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = blockZSize;
}

/**
 * Gets the launch parameters for launching a y diffusion kernel.
 */
void MPIMpdRdmeSolver::calculateReactionLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize/blockYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}


namespace mpimpdrdme_dev {
/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void MPI_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int z_start, const unsigned long long timestepHash, unsigned int* siteOverflowList)
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
__global__ void MPI_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList)
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
__global__ void MPI_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start)
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
__global__ void MPI_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start,   const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG)
#else
__global__ void MPI_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start)
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
        totalReactionPropensity += calculateReactionPropensity(siteType, (uint8_t*)particles, i, RLG);
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
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, globalIndex, totalReactionPropensity, timestepHash, RLG);
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

__global__ void mpi_correct_overflows(unsigned int* lattice, unsigned int* siteOverflowList)
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
		//printf("in: %d overflows, out: %d\n", total, siteOverflowList[0]);

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
