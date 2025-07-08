/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
 * Copyright 2012 Roberts Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 * 
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
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
 * Author(s): Elijah Roberts, Mike Hallock, Zane Thornburg
 */

#ifndef LM_RDME_MGPUMPDRDMESOLVER_H_
#define LM_RDME_MGPUMPDRDMESOLVER_H_

#include "cuda/lm_cuda.h"
#include "core/ResourceAllocator.h"
#include "rdme/RDMESolver.h"
#include "rdme/ByteLattice.h"
#include "rdme/GPUMapper/MultiGPUMapper.h"
#include "rdme/GPUMapper/SegmentDescriptor.h"

#if defined(MACOSX)
#include "GPUMapper/osx_barrier.h"
#endif

#define OVERFLOW_MODE_CLASSIC 0
#define OVERFLOW_MODE_RELAXED 1

using lm::main::ResourceAllocator;
using lm::rdme::RDMESolver;
using lm::rdme::Lattice;

namespace lm {

namespace io{
class Lattice;
class SpeciesCounts;
}
namespace rdme {

struct gpu_worker_thread_params;

class MGPUMpdRdmeSolver : public RDMESolver
{
    using CMESolver::hookSimulation;
    using RDMESolver::buildDiffusionModel;

public:
    MGPUMpdRdmeSolver();
    virtual ~MGPUMpdRdmeSolver();
    virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
    virtual bool needsReactionModel() {return true;}
    virtual bool needsDiffusionModel()  {return true;}
    virtual void buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionTypeA, const double * kA, const int * SA, const uint * DA, const uint kCols=1);
    virtual void buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData=true);
    virtual void generateTrajectory();
	virtual void setReactionRate(unsigned int rxid, float rate);

protected:
    virtual void allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing);
    virtual void writeLatticeData(double time, ByteLattice * lattice, lm::io::Lattice * latticeDataSet);
	virtual void writeLatticeSites(double time, ByteLattice * lattice);
    virtual void recordSpeciesCounts(double time, ByteLattice * lattice, lm::io::SpeciesCounts * speciesCountsDataSet);
    virtual void writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet);
	virtual void hookCheckSimulation(double time, ByteLattice * lattice);
    virtual int run_next_timestep(int gpu, uint32_t timestep);
    virtual uint64_t getTimestepSeed(uint32_t timestep, uint32_t substep);
	
	virtual void copyModelsToDevice(int gpu);
	virtual void initialize_decomposition();
	virtual void start_threads();
	virtual void stop_threads();
	virtual void computePropensities();

	// declare the child thread entry point as a friend function
	// so that it may call the run_thread method
	friend void* gpu_worker_thread(void *arg);
	virtual void* run_thread(int);
	virtual int handle_all_overflows();
	virtual int handle_overflows(int gpu, void *hptr, void *dptr, int ts);
	virtual int hookSimulation(double time, ByteLattice *lattice);

    virtual void calculateXLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateYLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateZLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateReactionLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);

protected:
    uint32_t seed;
    double tau;

    uint32_t overflowTimesteps;
    uint32_t overflowListUses;
    int      overflow_handling;

	// Stored model parameters for const memory
	unsigned int* model_reactionOrders;
	unsigned int* model_reactionSites;
	float*        model_reactionRates;

	unsigned int* model_D1;
	unsigned int* model_D2;

	int8_t*  model_S;
	float*   model_T;
	uint8_t* model_RL;

	size_t zeroOrderSize, firstOrderSize, secondOrderSize;
	float *zeroOrder, *firstOrder, *secondOrder;

	pthread_barrier_t start_barrier, stop_barrier, simulation_barrier;
	MultiGPUMapper *mapper;
	ResourceAllocator::ComputeResources *resources;
	gpu_worker_thread_params *threads;

	int timesteps_to_run;
	uint32_t current_timestep;
	bool reactionModelModified;
	double printPerfInterval;
	
	bool aggcopy_x_unpack;
	bool aggcopy_r_pack;
	bool use_spin_barrier;
};

struct gpu_worker_thread_params
{
    pthread_t thread;
    MGPUMpdRdmeSolver *runner;
    MultiGPUMapper *mapper;
    int gpu;
    int ngpus;
    int timesteps_to_run;
    bool lattice_synched;
    int load_balance_counter;

	// cuda objects
	unsigned int *dLattice, *dLatticeTmp;
	uint8_t *dSites;
	cudaStream_t stream1, stream2;
	unsigned int *h_overflows, *d_overflows;

	// kernel launch params
	dim3 grid_x, grid_y, grid_z, grid_r;
	dim3 threads_x, threads_y, threads_z, threads_r;

	// lattice segment geometry provided by the mapper
    SegmentDescriptor_s *segment;

#ifdef MPD_GLOBAL_S_MATRIX
	uint8_t *RLG;   // Device global memory pointer for RL matrix
	int8_t *SG;     // Device global memory pointer for S matrix
#endif
#ifdef MPD_GLOBAL_T_MATRIX
	float *TG;
#endif

#ifdef MPD_GLOBAL_R_MATRIX
    float *reactionRatesG;
    unsigned int *reactionOrdersG;
    unsigned int *reactionSitesG;
    unsigned int *D1G;
    unsigned int *D2G;
#endif

	float *propZeroOrder, *propFirstOrder, *propSecondOrder;
	

};



namespace mgpumpdrdme_dev {
__device__ inline size_t local_to_global(unsigned int x, unsigned int y, unsigned int z);
__device__ inline size_t local_index(unsigned int x, unsigned int y, unsigned int z);
__global__ void correct_overflows_mgpu(unsigned int* lattice, unsigned int* siteOverflowList);

#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void MGPU_precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrderG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2);
#else
__global__ void MGPU_precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2);
#endif
#else
__global__ void MGPU_precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2);
#endif


#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void MGPU_precomp_reaction_kernel_packing(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, unsigned int* buf_top, unsigned int* buf_bot);
#else
__global__ void MGPU_precomp_reaction_kernel_packing(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, unsigned int* buf_top, unsigned int* buf_bot);
#endif
#else
__global__ void MGPU_precomp_reaction_kernel_packing(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, int z_start, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, unsigned int* buf_top, unsigned int* buf_bot);
#endif

__global__ void MGPU_x_kernel_unpack(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int z_start, const unsigned long long timestepHash, unsigned int* siteOverflowList, unsigned int* buf_top, unsigned int *buf_bot);

__global__ void MGPU_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int z_start, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void MGPU_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void MGPU_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start);
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void MGPU_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start,   const int8_t* __restrict__ SG, const uint8_t* __restrict__ RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG);
#else
__global__ void MGPU_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start,   const int8_t* __restrict__ SG, const uint8_t* __restrict__ RLG);
#endif
#else
__global__ void MGPU_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start);
#endif


}

}
}

#endif
