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
 * Author(s): Elijah Roberts, Zane Thornburg
 */
#include <iostream>
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
#include "rdme/CudaByteLattice.h"
#include "rdme/MpdRdmeSolver.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include <sstream>
#include "core/Timer.h"
#include "io/ArbitraryH5.h"

#define MPD_WORDS_PER_SITE    (MPD_LATTICE_MAX_OCCUPANCY / 4) // 4 
#define MPD_APRON_SIZE        1

#include "cuda/constant.cuh"

namespace lm {
namespace rdme {
namespace mpdrdme_dev {
#include "rdme/dev/xor_random_dev.cu"
#include "rdme/dev/lattice_sim_1d_dev.cu"
#include "rdme/dev/byte_diffusion_1d_dev.cu"
#include "rdme/dev/byte_reaction_dev.cu"
}}}

extern bool globalAbort;

using std::map;
using lm::io::DiffusionModel;
using lm::rdme::Lattice;
using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {

MpdRdmeSolver::MpdRdmeSolver()
:RDMESolver(lm::rng::RandomGenerator::NONE),
 seed(0), tau(0.0),
 cudaOverflowList(NULL), cudaStream(0),
 overflowTimesteps(0), overflowListUses(0),
 model_reactionRates(NULL),
 zeroOrder(NULL), firstOrder(NULL), secondOrder(NULL)
{
}

void MpdRdmeSolver::initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources)
{
	RDMESolver::initialize(replicate, parameters, resources);

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

	
	#ifdef MPD_MAPPED_OVERFLOWS
	CUDA_EXCEPTION_CHECK(cudaHostAlloc(&cudaOverflowList,
		MPD_OVERFLOW_LIST_SIZE,
		cudaHostAllocPortable|cudaHostAllocMapped));
	memset(cudaOverflowList, 0, MPD_OVERFLOW_LIST_SIZE);
	#else
    // Allocate memory on the device for the exception list.
    CUDA_EXCEPTION_CHECK(cudaMalloc(&cudaOverflowList, MPD_OVERFLOW_LIST_SIZE)); //TODO: track memory usage.
    CUDA_EXCEPTION_CHECK(cudaMemset(cudaOverflowList, 0, MPD_OVERFLOW_LIST_SIZE));
	#endif

    // Create a stream for synchronizing the events.
    CUDA_EXCEPTION_CHECK(cudaStreamCreate(&cudaStream));
}

MpdRdmeSolver::~MpdRdmeSolver()
{
    // Free any device memory.
    if (cudaOverflowList != NULL)
    {
	#ifdef MPD_MAPPED_OVERFLOWS
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(cudaOverflowList)); //TODO: track memory usage.
	#else
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFree(cudaOverflowList)); //TODO: track memory usage.
	#endif
        cudaOverflowList = NULL;
    }

    // If we have created a stream, destroy it.
    if (cudaStream != 0)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaStreamDestroy(cudaStream));
        cudaStream = NULL;
    }

    // Free allocated memory in MpdRdmeSolver
    if (model_reactionRates)  {delete [] model_reactionRates;}

    if (zeroOrder)   {delete [] zeroOrder;}
	if (firstOrder)  {delete [] firstOrder;}
	if (secondOrder) {delete [] secondOrder;}
}

void MpdRdmeSolver::allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing)
{
	assert(bytes_per_particle == 1);
	if(particlesPerSite != MPD_LATTICE_MAX_OCCUPANCY)
	{
		Print::printf(Print::ERROR, "requested allocation for %d particles per site is not %d", particlesPerSite, MPD_LATTICE_MAX_OCCUPANCY);
		throw Exception("incorrect particle density");
	}

	if(latticeXSize % 32 != 0)
	{
		Print::printf(Print::ERROR, "The X dimension (%d) is not divisible by %d", latticeXSize, 32);
		throw Exception("incorrect lattice dimensions");
	}

	if(latticeYSize % TUNE_MPD_Y_BLOCK_Y_SIZE != 0)
	{
		Print::printf(Print::ERROR, "The Y dimension (%d) is not divisible by %d", latticeYSize, TUNE_MPD_Y_BLOCK_Y_SIZE);
		throw Exception("incorrect lattice dimensions");
	}

	if(latticeZSize % TUNE_MPD_Z_BLOCK_Z_SIZE != 0)
	{
		Print::printf(Print::ERROR, "The Z dimension (%d) is not divisible by %d", latticeZSize, TUNE_MPD_Z_BLOCK_Z_SIZE);
		throw Exception("incorrect lattice dimensions");
	}

	try
	{
		CudaByteLattice *cbl = new CudaByteLattice(latticeXSize, latticeYSize, latticeZSize, latticeSpacing, particlesPerSite);
		lattice = (Lattice *)cbl;
	}
	catch (CUDAException e)
	{
		Print::printf(Print::ERROR, e.what());
		throw Exception("Failed to create cuda lattice");
	}

}

void MpdRdmeSolver::buildModel(const uint numberSpeciesA,
                               const uint numberReactionsA,
                               const uint * initialSpeciesCountsA,
                               const uint * reactionTypesA,
                               const double * KA,
                               const int * SA,
                               const uint * DA,
                               const uint kCols)
{
    CMESolver::buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypesA, KA, SA, DA, kCols);

    // Get the time step.
    tau=atof((*parameters)["timestep"].c_str());
    if (tau <= 0.0) throw InvalidArgException("timestep", "A positive timestep must be specified for the solver.");

    // Make sure we can support the reaction model.
    if (numberReactions > MPD_MAX_REACTION_TABLE_ENTRIES) 
    {
        std::stringstream s;
        s << "The number of reaction table entries (" << numberReactions <<") exceeds the maximum supported by the solver (" << MPD_MAX_REACTION_TABLE_ENTRIES << ").";
        throw Exception(s.str().c_str());
    }
#ifndef MPD_GLOBAL_S_MATRIX
    if (numberSpecies*numberReactions > MPD_MAX_S_MATRIX_ENTRIES) 
    {
        std::stringstream s;
        s << "The number of S matrix entries (" << numberSpecies*numberReactions <<") exceeds the maximum supported by the solver (" << MPD_MAX_S_MATRIX_ENTRIES << ").";
        throw Exception(s.str().c_str());
    }
#endif


    // Setup the cuda reaction model.
    unsigned int * reactionOrders = new unsigned int[numberReactions];
    unsigned int * reactionSites = new unsigned int[numberReactions];
    unsigned int * D1 = new unsigned int[numberReactions];
    unsigned int * D2 = new unsigned int[numberReactions];

    for (uint i=0; i<numberReactions; i++)
    {
        if(reactionTypes[i] == ZerothOrderPropensityArgs::REACTION_TYPE)
        {
            reactionOrders[i] = MPD_ZERO_ORDER_REACTION;
        	reactionSites[i] = 0;
        	D1[i] = 0; 
        	D2[i] = 0;
        }
    	else if (reactionTypes[i] == FirstOrderPropensityArgs::REACTION_TYPE)
    	{
    		reactionOrders[i] = MPD_FIRST_ORDER_REACTION;
    		reactionSites[i] = 0;
    		D1[i] = ((FirstOrderPropensityArgs *)propensityFunctionArgs[i])->si+1;
    		D2[i] = 0;
    	}
    	else if (reactionTypes[i] == SecondOrderPropensityArgs::REACTION_TYPE)
    	{
    		reactionOrders[i] = MPD_SECOND_ORDER_REACTION;
    		reactionSites[i] = 0;
    		D1[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s1i+1;
    		D2[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s2i+1;
    	}
    	else if (reactionTypes[i] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
    	{
    		reactionOrders[i] = MPD_SECOND_ORDER_SELF_REACTION;
    		reactionSites[i] = 0;
    		D1[i] = ((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i])->si+1;
    		D2[i] = 0;
    	}
    	else
    	{
    		throw InvalidArgException("reactionTypeA", "the reaction type was not supported by the solver", reactionTypes[i]);
    	}
    }

    // Setup the cuda S matrix.
	// Transpose S on the device because we will want to read along the numSpecies axis
	// Before A1A2A3A4A5 B1B2B3B4B5 C1C2C3C4C5 ...
	// After A1B1C1 A2B2C2 ...
/*
    int8_t * tmpS = new int8_t[numberSpecies*numberReactions];
    for (uint i=0; i<numberSpecies*numberReactions; i++)
    {
    	tmpS[i] = S[i];
    }
*/
    int8_t * tmpS = new int8_t[numberSpecies*numberReactions];
	for(uint rx = 0; rx < numberReactions; rx++)
	{
		for (uint p=0; p<numberSpecies; p++)
		{
			//tmpS[numberSpecies*rx + p] = S[numberSpecies*rx + p];
			tmpS[rx * numberSpecies + p] = S[numberReactions*p + rx];
		}
	}

#ifdef MPD_GLOBAL_R_MATRIX
    // R matrix put in global memory
    //cudaMalloc(&numberReactionsG, sizeof(unsigned int));
    //cudaMemcpy(numberReactionsG, &numberReactions, sizeof(unsigned int), cudaMemcpyHostToDevice);
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberReactionsC, &numberReactions, sizeof(unsigned int)));
    cudaMalloc(&reactionOrdersG, numberReactions*sizeof(unsigned int));
    cudaMemcpy(reactionOrdersG, reactionOrders, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&reactionSitesG, numberReactions*sizeof(unsigned int));
    cudaMemcpy(reactionSitesG, reactionSites, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&D1G, numberReactions*sizeof(unsigned int));
    cudaMemcpy(D1G, D1, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&D2G, numberReactions*sizeof(unsigned int));
    cudaMemcpy(D2G, D2, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
#else
    // Copy the reaction model and S matrix to constant memory on the GPU.
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberReactionsC, &numberReactions, sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionOrdersC, reactionOrders, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionSitesC, reactionSites, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(D1C, D1, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(D2C, D2, numberReactions*sizeof(unsigned int)));
#endif

#ifdef MPD_GLOBAL_S_MATRIX
	// If S matrix stored in global memory, allocate space and perform copy
	cudaMalloc(&SG, numberSpecies*numberReactions * sizeof(int8_t));
	cudaMemcpy(SG, tmpS, numberSpecies*numberReactions * sizeof(int8_t), cudaMemcpyHostToDevice);
#else
	// S matrix is in constant memory
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(SC, tmpS, numberSpecies*numberReactions*sizeof(int8_t)));
#endif

    // Free any temporary resources.
    delete [] reactionSites;
    delete [] D1;
    delete [] D2;
    delete [] reactionOrders;
    delete [] tmpS;
}

void MpdRdmeSolver::buildDiffusionModel(const uint numberSiteTypesA,
                                        const double * DFA,
                                        const uint * RLA,
                                        lattice_size_t latticeXSize,
                                        lattice_size_t latticeYSize,
                                        lattice_size_t latticeZSize,
                                        site_size_t particlesPerSite,
                                        const unsigned int bytes_per_particle,
                                        si_dist_t latticeSpacing,
                                        const uint8_t * latticeData,
                                        const uint8_t * latticeSitesData,
                                        bool rowMajorData)
{
    RDMESolver::buildDiffusionModel(numberSiteTypesA, DFA, RLA,
                                    latticeXSize, latticeYSize, latticeZSize,
                                    particlesPerSite, bytes_per_particle,
                                    latticeSpacing, latticeData, latticeSitesData, rowMajorData);

    // Get the time step.
    tau=atof((*parameters)["timestep"].c_str());
    if (tau <= 0.0) throw InvalidArgException("timestep", "A positive timestep must be specified for the solver.");

    maxSiteCounts.resize(numberSiteTypes);
    currentMaxSiteCounts.resize(numberSiteTypes);
    std::fill(maxSiteCounts.begin(), maxSiteCounts.end(), 0);

    maxParticleCounts.resize(numberSpecies+1);
    currentMaxParticleCounts.resize(numberSpecies+1);
    std::fill(maxParticleCounts.begin(), maxParticleCounts.end(), 0);
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
    float * T = new float[DFmatrixSize];
    for (uint i=0; i<DFmatrixSize; i++)
    {
        float q=(float)(DF[i]*tau/pow(latticeSpacing,2));
        if (q > 0.50f) throw InvalidArgException("D", "The specified diffusion coefficient is too high for the diffusion model.");
        T[i] = q;
    }

    // Setup the cuda reaction location matrix.
    uint8_t * tmpRL = new uint8_t[numberReactions*numberSiteTypes];
    for (uint i=0; i<numberReactions*numberSiteTypes; i++)
    {
    	tmpRL[i] = RL[i];
    }

    // Copy the diffusion model to constant memory on the GPU.
#ifdef MPD_GLOBAL_T_MATRIX
	CUDA_EXCEPTION_CHECK(cudaMalloc(&TG, DFmatrixSize*sizeof(float)));
	CUDA_EXCEPTION_CHECK(cudaMemcpy(TG, T, DFmatrixSize*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(TC, &TG, sizeof(float*)));
#else
	CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(TC, T, DFmatrixSize*sizeof(float)));
#endif

    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberSpeciesC, &numberSpecies, sizeof(numberSpeciesC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberSiteTypesC, &numberSiteTypes, sizeof(numberSiteTypesC)));
    const unsigned int latticeXYSize = latticeXSize*latticeYSize;
    const unsigned int latticeXYZSize = latticeXSize*latticeYSize*latticeZSize;
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXSizeC, &latticeXSize, sizeof(latticeYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeYSizeC, &latticeYSize, sizeof(latticeYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeZSizeC, &latticeZSize, sizeof(latticeZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(global_latticeZSizeC, &latticeZSize, sizeof(latticeZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXYSizeC, &latticeXYSize, sizeof(latticeXYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXYZSizeC, &latticeXYZSize, sizeof(latticeXYZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(global_latticeXYZSizeC, &latticeXYZSize, sizeof(latticeXYZSizeC)));
#ifdef MPD_GLOBAL_S_MATRIX
	// Store RL in global memory too, since I'm going to assume if S is too big, then RL is too.
	cudaMalloc(&RLG, numberReactions*numberSiteTypes * sizeof(uint8_t));
	cudaMemcpy(RLG, tmpRL, numberReactions*numberSiteTypes * sizeof(uint8_t), cudaMemcpyHostToDevice);
#else
	// RL is stored in constant memory
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(RLC, tmpRL, numberReactions*numberSiteTypes*sizeof(uint8_t)));
#endif
    delete [] tmpRL;
    delete [] T;

    // Set the cuda reaction model rates now that we have the subvolume size.
    model_reactionRates = new float[numberReactions];

    // Set up pre-configured propensity matrices
	zeroOrderSize=numberSiteTypes;
	firstOrderSize=numberSpecies*numberSiteTypes;
	secondOrderSize=numberSpecies*numberSpecies*numberSiteTypes;
	zeroOrder=new float[zeroOrderSize];
	firstOrder=new float[firstOrderSize];
	secondOrder=new float[secondOrderSize];

    #ifdef MPD_GLOBAL_R_MATRIX
    cudaMalloc(&reactionRatesG, numberReactions*sizeof(float));
    #endif
		
	cudaMalloc(&propZeroOrder, zeroOrderSize*sizeof(float));
	cudaMalloc(&propFirstOrder, firstOrderSize*sizeof(float));
	cudaMalloc(&propSecondOrder, secondOrderSize*sizeof(float));

    computePropensities();
    copyModelsToDevice();
	reactionModelModified = false;
}

void MpdRdmeSolver::generateTrajectory()
{
    // Shadow the lattice member as a cuda lattice.
    CudaByteLattice * lattice = (CudaByteLattice *)this->lattice;

    // Get the interval for writing species counts and lattices.
    uint32_t speciesCountsWriteInterval=atof((*parameters)["writeInterval"].c_str());
    uint32_t nextSpeciesCountsWriteTime = speciesCountsWriteInterval;
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpeciesToTrack);
    speciesCountsDataSet.set_number_entries(0);
    uint32_t latticeWriteInterval=atof((*parameters)["latticeWriteInterval"].c_str());
    uint32_t nextLatticeWriteTime = latticeWriteInterval;
    lm::io::Lattice latticeDataSet;

    // Get the simulation time limit.
    double maxTime=atof((*parameters)["maxTime"].c_str());

    Print::printf(Print::INFO,
                  "Running mpd rdme simulation with %d species, %d reactions, %d site types for %e s with tau %e. Writing species at %e and lattice at %e intervals",
                  numberSpecies, numberReactions, numberSiteTypes,
                  maxTime, tau,
                  speciesCountsWriteInterval, latticeWriteInterval);

    // Set the initial time.
    double time = 0.0;
	uint32_t current_timestep=0;

	bool hookEnabled=false;
	uint32_t nextHookTime=0;
	uint32_t hookInterval=0;

    // Find out at what interval to hook simulations
	if((*parameters)["hookInterval"] != "")
	{
		hookInterval=atof((*parameters)["hookInterval"].c_str());
		hookEnabled=hookInterval >= 0.0;
		nextHookTime=hookInterval;
	}

	// Find out at what interval to write status messages
	double printPerfInterval = 60;	
	if((*parameters)["perfPrintInterval"] != "")
	{
		printPerfInterval=atof((*parameters)["perfPrintInterval"].c_str());
	}
	 

	// Call beginning hook.  Any modifications to the lattice will be
	// accounted for because we have not yet copied to device memory
	// thus we do not need to check return value
	onBeginTrajectory(lattice);

    // Synchronize the cuda memory.
    lattice->copyToGPU();

    // Record the initial species counts.
    recordSpeciesCounts(time, lattice, &speciesCountsDataSet);

    // Write the initial lattice.
    writeLatticeData(time, lattice, &latticeDataSet);

    // initialize max counts from initial particle/site lattice
    initMaxCounts(lattice);

    Timer timer;
    timer.tick();
    double lastT=0;
    int lastSteps=0;

    // Perform an initial hook check
    hookCheckSimulation(time, lattice);

    while (time < maxTime)
    {
        
        if(globalAbort)
	    {
		    printf("Global abort: terminating solver\n");
		    break;
	    }

	
        lastT += timer.tock();
        lastSteps += 1;

        if ( lastT >= printPerfInterval)
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

            Print::printf(Print::INFO, "Average walltime per timestep: %.2f ms. Progress: %.4fs/%.4fs (% .3g%% done / %.2g %s walltime remaining)",
                                       1000.0*stepTime, time, maxTime, 100.0*time/maxTime, completionTime, units.c_str());

            lastT = 0;
            lastSteps = 0;
        }

        timer.tick();

        // Run the next timestep.
        runTimestep(lattice, current_timestep++);

        // Update the time.
        time = current_timestep*tau;

	    // See if we need to write out the any data.
        if (current_timestep >= nextLatticeWriteTime
            || current_timestep >= nextSpeciesCountsWriteTime
            || (hookEnabled && current_timestep >= nextHookTime))
        {
            // Synchronize the lattice.
            lattice->copyFromGPU();

            Print::printf(Print::INFO, "Time is %.14f", time);

            // Check if we need to execute the hook
            if (hookEnabled && current_timestep >= nextHookTime)
            {
                Print::printf(Print::INFO, "Hook time is %.14f, in steps is %d", time, nextHookTime);
                hookCheckSimulation(time, lattice);
                nextHookTime += hookInterval;
                Print::printf(Print::INFO, "Next hook time is %d", nextHookTime);
            }

	        // See if we need to write the lattice.
            if (current_timestep >= nextLatticeWriteTime)
            {
                PROF_BEGIN(PROF_SERIALIZE_LATTICE);
                Print::printf(Print::INFO, "Lattice write time is %.14f, in steps is %.d", time, nextLatticeWriteTime);
                writeLatticeData(time, lattice, &latticeDataSet);
                nextLatticeWriteTime += latticeWriteInterval;
                Print::printf(Print::INFO, "Next lattice write time is %.d", nextLatticeWriteTime);
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
    }

	lattice->copyFromGPU();
	onEndTrajectory(lattice);

    // Write any remaining species counts.
    writeSpeciesCounts(&speciesCountsDataSet);
    writeMaxCounts();
}

void MpdRdmeSolver::setReactionRate(unsigned int rxid, float rate)
{
    if (reactionTypes[rxid] ==  ZerothOrderPropensityArgs::REACTION_TYPE)
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

void MpdRdmeSolver::computePropensities()
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

void MpdRdmeSolver::copyModelsToDevice()
{
#ifdef MPD_GLOBAL_R_MATRIX
    cudaMemcpy(reactionRatesG, model_reactionRates, numberReactions*sizeof(float), cudaMemcpyHostToDevice);
#else
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionRatesC, model_reactionRates, numberReactions*sizeof(float)));
#endif

	cudaMemcpy(propZeroOrder, zeroOrder, zeroOrderSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(propFirstOrder, firstOrder, firstOrderSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(propSecondOrder, secondOrder, secondOrderSize*sizeof(float), cudaMemcpyHostToDevice);
}

int MpdRdmeSolver::hookSimulation(double time, CudaByteLattice *lattice)
{
	// Overload this function in derivative classes
	// Return 0 if the lattice state is unchanged
	// Return 1 if the lattice state has been modified, 
    //          and it needs to be copied back to the GPU
	// Return 2 if lattice sites have changed, and should
	//			be copied back to the GPU *and* be recorded
	//			in the output file
	return 0;
}

int MpdRdmeSolver::onWriteLattice(double time, CudaByteLattice *lattice)
{
	// Overload this function in derivative classes
	// Return 0 if the lattice state is unchanged
	// Return 1 if the lattice state has been modified, 
    //          and it needs to be copied back to the GPU.
	// Return 2 if lattice sites have changed, and should
	//			be copied back to the GPU *and* be recorded
	//			in the output file
	return 0;
}

int MpdRdmeSolver::onBeginTrajectory(CudaByteLattice *lattice)
{
	return 0;
}

int MpdRdmeSolver::onEndTrajectory(CudaByteLattice *lattice)
{
	return 0;
}


void MpdRdmeSolver::initMaxCounts(CudaByteLattice *lattice)
{
    uint8_t *sites = reinterpret_cast<uint8_t*>(lattice->getSitesMemory());
    const size_t ns = lattice->getNumberSites();
    std::fill(currentMaxSiteCounts.begin(), currentMaxSiteCounts.end(), 0);
    for (size_t i=0; i < ns; i++) {
        if (sites[i] > 0) {
            assert(sites[i] < currentMaxSiteCounts.size());
            currentMaxSiteCounts[sites[i]]++;
        }
    }
    for (int i=0; i < maxSiteCounts.size(); i++) {
        maxSiteCounts[i] = std::max(maxSiteCounts[i], currentMaxSiteCounts[i]);
    }

    uint8_t *particles = reinterpret_cast<uint8_t*>(lattice->getParticlesMemory());
    const size_t np = MPD_LATTICE_MAX_OCCUPANCY * lattice->getNumberSites();
    std::fill(currentMaxParticleCounts.begin(), currentMaxParticleCounts.end(), 0);
    for (size_t i=0; i < np; i++) {
        if (particles[i] > 0) {
            currentMaxParticleCounts[particles[i]]++;
        }
    }
    for (int i=0; i < maxParticleCounts.size(); i++) {
        maxParticleCounts[i] = std::max(maxParticleCounts[i], currentMaxParticleCounts[i]);
    }
}



void MpdRdmeSolver::writeLatticeData(double time, CudaByteLattice * lattice, lm::io::Lattice * latticeDataSet)
{
    Print::printf(Print::DEBUG, "Writing lattice at %e s", time);

    // Record the lattice data.
    latticeDataSet->Clear();
    latticeDataSet->set_lattice_x_size(lattice->getSize().x);
    latticeDataSet->set_lattice_y_size(lattice->getSize().y);
    latticeDataSet->set_lattice_z_size(lattice->getSize().z);
    latticeDataSet->set_particles_per_site(lattice->getMaxOccupancy());
    latticeDataSet->set_time(time);

    // update max counts

    uint8_t *particles = reinterpret_cast<uint8_t*>(lattice->getParticlesMemory());
    const size_t np = MPD_LATTICE_MAX_OCCUPANCY * lattice->getNumberSites();
    std::fill(currentMaxParticleCounts.begin(), currentMaxParticleCounts.end(), 0);
    for (size_t i=0; i < np; i++) {
        if (particles[i] > 0) {
            currentMaxParticleCounts[particles[i]]++;
        }
    }
    for (int i=0; i < maxParticleCounts.size(); i++) {
        maxParticleCounts[i] = std::max(maxParticleCounts[i], currentMaxParticleCounts[i]);
    }

    // Push it to the output queue.
    size_t payloadSize = size_t(lattice->getSize().x)*size_t(lattice->getSize().y)*size_t(lattice->getSize().z)*size_t(lattice->getMaxOccupancy())*sizeof(uint8_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::BYTE_LATTICE, replicate, latticeDataSet, lattice, payloadSize, &lm::rdme::ByteLattice::nativeSerialize);
}


void MpdRdmeSolver::writeLatticeSites(double time, CudaByteLattice * lattice)
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

    // update max counts

    uint8_t *sites = reinterpret_cast<uint8_t*>(lattice->getSitesMemory());
    const size_t ns = lattice->getNumberSites();
    std::fill(currentMaxSiteCounts.begin(), currentMaxSiteCounts.end(), 0);
    for (size_t i=0; i < ns; i++) {
        if (sites[i] > 0) {
            assert(sites[i] < currentMaxSiteCounts.size());
            currentMaxSiteCounts[sites[i]]++;
        }
    }
    for (int i=0; i < maxSiteCounts.size(); i++) {
        maxSiteCounts[i] = std::max(maxSiteCounts[i], currentMaxSiteCounts[i]);
    }


    // Push it to the output queue.
    size_t payloadSize = lattice->getSize().x*lattice->getSize().y*lattice->getSize().z*sizeof(uint8_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SITE_LATTICE, replicate, &latticeDataSet, lattice, payloadSize, &lm::rdme::ByteLattice::nativeSerializeSites);
}

void 
MpdRdmeSolver::writeMaxCounts()
{
    std::vector<size_t> shape_p(1);
    shape_p[0] = maxParticleCounts.size();
    H5MetaData md_p = make_H5_meta<uint32_t>(H5MetaData::NEW_DATASET, shape_p, "MaxParticleCounts");
    lm::main::DataOutputQueue::getInstance()->pushDataSet(md_p, replicate, &(maxParticleCounts[0]));

    std::vector<size_t> shape_s(1);
    shape_s[0] = maxSiteCounts.size();
    H5MetaData md_s = make_H5_meta<uint32_t>(H5MetaData::NEW_DATASET, shape_s, "MaxSiteCounts");
    lm::main::DataOutputQueue::getInstance()->pushDataSet(md_s, replicate, &(maxSiteCounts[0]));
}
void MpdRdmeSolver::recordSpeciesCounts(double time, CudaByteLattice * lattice, lm::io::SpeciesCounts * speciesCountsDataSet)
{
    std::map<particle_t,uint> particleCounts = lattice->getParticleCounts();
    speciesCountsDataSet->set_number_entries(speciesCountsDataSet->number_entries()+1);
    speciesCountsDataSet->add_time(time);
    for (particle_t p=0; p<numberSpeciesToTrack; p++)
    {
        speciesCountsDataSet->add_species_count((particleCounts.count(p+1)>0)?particleCounts[p+1]:0);
    }
}

void MpdRdmeSolver::writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet)
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

void MpdRdmeSolver::hookCheckSimulation(double time, CudaByteLattice * lattice)
{
    switch(hookSimulation(time, lattice))
    {
        case 0:
            break; 

        case 1:
            lattice->copyToGPU();
            break;

        case 2:
            lattice->copyToGPU();
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

uint64_t MpdRdmeSolver::getTimestepSeed(uint32_t timestep, uint32_t substep)
{
    uint64_t timestepHash = (((((uint64_t)seed)<<30)+timestep)<<2)+substep;
    timestepHash = timestepHash * 3202034522624059733ULL + 4354685564936845319ULL;
    timestepHash ^= timestepHash >> 20; timestepHash ^= timestepHash << 41; timestepHash ^= timestepHash >> 5;
    timestepHash *= 7664345821815920749ULL;
    return timestepHash;
}

void MpdRdmeSolver::runTimestep(CudaByteLattice * lattice, uint32_t timestep)
{
    PROF_BEGIN(PROF_MPD_TIMESTEP);

    if(reactionModelModified)
        copyModelsToDevice();

    // Calculate some properties of the lattice.
    lattice_coord_t size = lattice->getSize();
    const unsigned int latticeXSize = size.x;
    const unsigned int latticeYSize = size.y;
    const unsigned int latticeZSize = size.z;

    dim3 gridSize, threadBlockSize;

    // Execute the kernel for the x direction.
    PROF_CUDA_START(cudaStream);

    PROF_CUDA_BEGIN(PROF_MPD_X_DIFFUSION,cudaStream);
    #ifdef MPD_CUDA_3D_GRID_LAUNCH
    calculateXLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::mpd_x_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,0), (unsigned int*)cudaOverflowList)));
    #else
    unsigned int gridXSize;
    calculateXLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::mpd_x_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), gridXSize, getTimestepSeed(timestep,0), (unsigned int*)cudaOverflowList)));
    #endif
    PROF_CUDA_END(PROF_MPD_X_DIFFUSION,cudaStream);
    lattice->swapSrcDest();

    // Execute the kernel for the y direction.
    PROF_CUDA_BEGIN(PROF_MPD_Y_DIFFUSION,cudaStream);
    #ifdef MPD_CUDA_3D_GRID_LAUNCH
    calculateYLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_Y_BLOCK_X_SIZE, TUNE_MPD_Y_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::mpd_y_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,1), (unsigned int*)cudaOverflowList)));
    #else
    calculateYLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_Y_BLOCK_X_SIZE, TUNE_MPD_Y_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::mpd_y_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), gridXSize, getTimestepSeed(timestep,1), (unsigned int*)cudaOverflowList)));
    #endif
    PROF_CUDA_END(PROF_MPD_Y_DIFFUSION,cudaStream);
    lattice->swapSrcDest();

    // Execute the kernel for the z direction.
    PROF_CUDA_BEGIN(PROF_MPD_Z_DIFFUSION,cudaStream);
    #ifdef MPD_CUDA_3D_GRID_LAUNCH
    calculateZLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_Z_BLOCK_X_SIZE, TUNE_MPD_Z_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::mpd_z_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,2), (unsigned int*)cudaOverflowList)));
    #else
    calculateZLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_Z_BLOCK_X_SIZE, TUNE_MPD_Z_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::mpd_z_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), gridXSize, getTimestepSeed(timestep,2), (unsigned int*)cudaOverflowList)));
    #endif
    PROF_CUDA_END(PROF_MPD_Z_DIFFUSION,cudaStream);
    lattice->swapSrcDest();

    if (numberReactions > 0)
    {
        // Execute the kernel for the reaction, this kernel updates the lattice in-place, so only the src pointer is passed.
        PROF_CUDA_BEGIN(PROF_MPD_REACTION,cudaStream);
        #ifdef MPD_CUDA_3D_GRID_LAUNCH
        calculateReactionLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_REACTION_BLOCK_X_SIZE, TUNE_MPD_REACTION_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
#ifdef MPD_GLOBAL_S_MATRIX

#ifdef MPD_GLOBAL_R_MATRIX
        #ifdef MPD_FREAKYFAST
                CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::precomp_reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG, propZeroOrder, propFirstOrder, propSecondOrder)));
        #else
                CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG)));
        #endif
#else
	#ifdef MPD_FREAKYFAST
		CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::precomp_reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG, propZeroOrder, propFirstOrder, propSecondOrder)));
	#else
		CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG)));
	#endif
#endif

#else
        CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList)));
#endif
        #else
        calculateReactionLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_REACTION_BLOCK_X_SIZE, TUNE_MPD_REACTION_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
#ifdef MPD_GLOBAL_S_MATRIX
        CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), gridXSize, getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG)));
#else
        CUDA_EXCEPTION_EXECUTE((mpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), gridXSize, getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList)));
#endif
        #endif
        PROF_CUDA_END(PROF_MPD_REACTION,cudaStream);
    }

    if(overflow_handling == OVERFLOW_MODE_RELAXED)
	mpdrdme_dev::correct_overflows<<<dim3(1,1,1), dim3(TUNE_MPD_MAX_PARTICLE_OVERFLOWS,1,1),0,cudaStream>>>((unsigned int*)lattice->getGPUMemorySrc(), (unsigned int*)cudaOverflowList);

    // Wait for the kernels to complete.
    PROF_BEGIN(PROF_MPD_SYNCHRONIZE);
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(cudaStream));
    PROF_END(PROF_MPD_SYNCHRONIZE);

    if(overflow_handling == OVERFLOW_MODE_RELAXED)
	{
		uint unresolved = ((unsigned int*)cudaOverflowList)[0];
		if (unresolved > 0)
		{
            Print::printf(Print::WARNING, "%d unresolved overflows", unresolved);
		}
	}

    if(overflow_handling == OVERFLOW_MODE_CLASSIC)
    {

    // Handle any particle overflows.
    PROF_BEGIN(PROF_MPD_OVERFLOW);
    
    overflowTimesteps++;
#ifndef MPD_MAPPED_OVERFLOWS
    uint32_t overflowList[1+2*TUNE_MPD_MAX_PARTICLE_OVERFLOWS];
    CUDA_EXCEPTION_CHECK(cudaMemcpy(overflowList, cudaOverflowList, MPD_OVERFLOW_LIST_SIZE, cudaMemcpyDeviceToHost));
#else
	uint32_t *overflowList = (uint32_t*)cudaOverflowList;
#endif

    uint numberExceptions = overflowList[0];
    if (numberExceptions > 0)
    {
        Print::printf(Print::DEBUG, "%d overflows", numberExceptions);

        // Make sure we did not exceed the overflow buffer.
        if (numberExceptions > TUNE_MPD_MAX_PARTICLE_OVERFLOWS)
            throw Exception("Too many particle overflows for the available buffer", numberExceptions);

        // Synchronize the lattice.
        lattice->copyFromGPU();
        
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
        
            
        // Copy the changes back to the GPU.
        lattice->copyToGPU();
            
        // Reset the overflow list.
        CUDA_EXCEPTION_CHECK(cudaMemset(cudaOverflowList, 0, MPD_OVERFLOW_LIST_SIZE));

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
    } // End classic overflow

    PROF_CUDA_FINISH(cudaStream);
    PROF_END(PROF_MPD_TIMESTEP);
}

#ifdef MPD_CUDA_3D_GRID_LAUNCH
/**
 * Gets the launch parameters for launching an x diffusion kernel.
 */
void MpdRdmeSolver::calculateXLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
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
void MpdRdmeSolver::calculateYLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
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
void MpdRdmeSolver::calculateZLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
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
void MpdRdmeSolver::calculateReactionLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize/blockYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

#else
/**
 * Gets the launch parameters for launching an x diffusion kernel.
 */
void MpdRdmeSolver::calculateXLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    unsigned int xBlockXSize = min(maxXBlockSize,latticeXSize);
    *gridXSize = latticeXSize/xBlockXSize;
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

    (*gridSize).x = (*gridXSize)*latticeYSize;
    (*gridSize).y = latticeZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = xBlockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a y diffusion kernel.
 */
void MpdRdmeSolver::calculateYLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    *gridXSize = latticeXSize/blockXSize;
    (*gridSize).x = (*gridXSize)*(latticeYSize/blockYSize);
    (*gridSize).y = latticeZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a z diffusion kernel.
 */
void MpdRdmeSolver::calculateZLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    *gridXSize = latticeXSize/blockXSize;
    (*gridSize).x = (*gridXSize)*(latticeYSize);
    (*gridSize).y = latticeZSize/blockZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = blockZSize;
}

/**
 * Gets the launch parameters for launching a reaction diffusion kernel.
 */
void MpdRdmeSolver::calculateReactionLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    *gridXSize = latticeXSize/blockXSize;
    (*gridSize).x = (*gridXSize)*(latticeYSize/blockYSize);
    (*gridSize).y = latticeZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}
#endif


namespace mpdrdme_dev {
/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
#ifdef MPD_CUDA_3D_GRID_LAUNCH
__global__ void __launch_bounds__(TUNE_MPD_X_BLOCK_MAX_X_SIZE,1) mpd_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
__global__ void __launch_bounds__(TUNE_MPD_X_BLOCK_MAX_X_SIZE,1) mpd_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeXIndex = (bx*blockDim.x) + threadIdx.x;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (by*latticeXSizeC) + latticeXIndex;
    unsigned int windowIndex = threadIdx.x+MPD_APRON_SIZE;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment.
    __shared__ unsigned int window[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_X_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyXWindowFromLattice(bx, inLattice, window, latticeIndex, latticeXIndex, windowIndex);
    copyXWindowFromSites(bx, inSites, sitesWindow, latticeIndex, latticeXIndex, windowIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];

    // Make the choices.
    makeXDiffusionChoices(window, sitesWindow, choices, latticeIndex, latticeXIndex, windowIndex, blockDim.x, timestepHash);
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
#ifdef MPD_CUDA_3D_GRID_LAUNCH
__global__ void __launch_bounds__(TUNE_MPD_Y_BLOCK_X_SIZE*TUNE_MPD_Y_BLOCK_Y_SIZE,1) mpd_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
__global__ void __launch_bounds__(TUNE_MPD_Y_BLOCK_X_SIZE*TUNE_MPD_Y_BLOCK_Y_SIZE,1) mpd_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (latticeYIndex*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;
    unsigned int windowYIndex = threadIdx.y+MPD_APRON_SIZE;
    unsigned int windowIndex = (windowYIndex*blockDim.x) + threadIdx.x;

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
    makeYDiffusionChoices(window, sitesWindow, choices, latticeIndex, latticeYIndex, windowIndex, windowYIndex, timestepHash);
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
#ifdef MPD_CUDA_3D_GRID_LAUNCH
__global__ void __launch_bounds__(TUNE_MPD_Z_BLOCK_X_SIZE*TUNE_MPD_Z_BLOCK_Z_SIZE,1) mpd_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
__global__ void __launch_bounds__(TUNE_MPD_Z_BLOCK_X_SIZE*TUNE_MPD_Z_BLOCK_Z_SIZE,1) mpd_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz*blockDim.z) + threadIdx.z;
    unsigned int latticeIndex = (latticeZIndex*latticeXYSizeC) + (by*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;
    unsigned int windowZIndex = threadIdx.z+MPD_APRON_SIZE;
    unsigned int windowIndex = (windowZIndex*blockDim.x) + threadIdx.x;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[MPD_Z_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_Z_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyZWindowFromLattice(inLattice, window, latticeIndex, latticeZIndex, windowIndex, windowZIndex);
    copyZWindowFromSites(inSites, sitesWindow, latticeIndex, latticeZIndex, windowIndex, windowZIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_Z_WINDOW_SIZE*MPD_WORDS_PER_SITE];

    // Make the choices.
    makeZDiffusionChoices(window, sitesWindow, choices, latticeIndex, latticeZIndex, windowIndex, windowZIndex, timestepHash);
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
#ifdef MPD_CUDA_3D_GRID_LAUNCH
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG,  const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG)
#endif
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList)
#endif
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG)
#endif
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
#endif
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (latticeYIndex*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;

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

	// If propensity is zero, no reaction can occur.
	if(totalReactionPropensity == 0.0f)
		return;

    // See if a reaction occurred at the site.
    float reactionProbability = calculateReactionProbability(totalReactionPropensity);
    unsigned int reactionOccurred = checkForReaction(latticeIndex, reactionProbability, timestepHash);

    // If there was a reaction, process it.
    if (reactionOccurred)
    {
        // Figure out which reaction occurred.
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, latticeIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, latticeIndex, totalReactionPropensity, timestepHash, RLG);
#endif
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, latticeIndex, totalReactionPropensity, timestepHash);
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
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2)
#endif
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2)
#endif
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (latticeYIndex*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;

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
    unsigned int reactionOccurred = checkForReaction(latticeIndex, reactionProbability, timestepHash);

    // If there was a reaction, process it.
    if (reactionOccurred)
    {
        // Figure out which reaction occurred.
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, latticeIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, latticeIndex, totalReactionPropensity, timestepHash, RLG);
#endif
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, (uint8_t*)particles, latticeIndex, totalReactionPropensity, timestepHash);
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

__global__ void correct_overflows(unsigned int* lattice, unsigned int* siteOverflowList)
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


	// Sanity test kernel to compare kernel outputs
	__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) sanity_check(const unsigned int* L1, const unsigned int* L2)
	{
		// Figure out the offset of this thread in the lattice and the lattice segment.
		unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
		unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
		unsigned int latticeIndex = (bz*latticeXYSizeC) + (latticeYIndex*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;

		///////////////////////////////////////////
		// Load the particles and site.          //
		///////////////////////////////////////////

		unsigned int particles1;
		unsigned int particles2;
		for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
		{
			particles1 = L1[latticeIndex+latticeOffset];
			particles2 = L2[latticeIndex+latticeOffset];

			if(particles1 != particles2)
			{
				printf("**\n ********* Sanity failure block(%d,%d,%d) thread(%d,%d,%d) w=%d L1=%8X L2=%8X\n**\n",
					bx,by,bz,
					threadIdx.x, threadIdx.y, threadIdx.z,
					w, particles1, particles2);
			}
		}
	}

}

}
}
