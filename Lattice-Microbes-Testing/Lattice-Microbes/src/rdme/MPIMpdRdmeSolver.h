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
 * Overflow algorithm in RDME solvers and CPU assignment (2012)
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
 * Author(s): Mike Hallock
 */

#ifndef LM_RDME_MPIMPDRDMESOLVER_H_
#define LM_RDME_MPIMPDRDMESOLVER_H_

#include "cuda/lm_cuda.h"

#include "core/ResourceAllocator.h"
#include "rdme/RDMESolver.h"
#include "rdme/ByteLattice.h"
#include "rdme/GPUMapper/ZDivMPIGPUMapper.h"
#include "rdme/GPUMapper/SegmentDescriptor.h"

#if defined(MACOSX)
#include "rdme/GPUMapper/osx_barrier.h"
#endif

using lm::main::ResourceAllocator;
using lm::rdme::RDMESolver;
using lm::rdme::Lattice;

namespace lm {

namespace io{
class Lattice;
class SpeciesCounts;
}
namespace rdme {

class MPIMpdRdmeSolver : public RDMESolver
{
public:
    MPIMpdRdmeSolver();
    virtual ~MPIMpdRdmeSolver();
    virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
    virtual bool needsReactionModel() {return true;}
    virtual bool needsDiffusionModel()  {return true;}
    virtual void buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionTypeA, const double * kA, const int * SA, const uint * DA, const uint kCols=1);
    virtual void buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData=true);
    virtual void generateTrajectory();

protected:
    virtual void allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing);
    virtual void writeLatticeData(double time, ByteLattice * lattice, lm::io::Lattice * latticeDataSet);
    virtual void recordSpeciesCounts(double time, ByteLattice * lattice, lm::io::SpeciesCounts * speciesCountsDataSet);
    virtual void writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet);
    virtual int run_next_timestep(uint32_t timestep);
    virtual uint64_t getTimestepSeed(uint32_t timestep, uint32_t substep);
	virtual void copyModelsToDevice();
	virtual void initialize_decomposition();

	virtual void prepare_gpu();


    virtual void calculateXLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateYLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateZLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateReactionLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);

protected:
    uint32_t seed;
    void * cudaOverflowList;
    double tau;
    uint32_t overflowTimesteps;
    uint32_t overflowListUses;

	// Stored model parameters for const memory
	unsigned int* model_reactionOrders;
	unsigned int* model_reactionSites;
	unsigned int* model_D1;
	unsigned int* model_D2;
	int8_t* model_S;
	float* model_T;
	uint8_t* model_RL;
	float* model_reactionRates;

	int world_size, rank;

#ifdef MPD_GLOBAL_S_MATRIX
	uint8_t *RLG;   // Device global memory pointer for RL matrix
	int8_t *SG;     // Device global memory pointer for S matrix
#endif

#ifdef MPD_GLOBAL_T_MATRIX
	float *TG;
#endif
	
	ZDivMPIGPUMapper *mapper;
	ResourceAllocator::ComputeResources *resources;

	int timesteps_to_run;
	uint32_t absolute_timestep;

    int gpu;
    int ngpus;
    bool lattice_synched;

	// cuda objects
	unsigned int *dLattice, *dLatticeTmp;
	uint8_t *dSites;
	cudaStream_t stream1, stream2;
	cudaEvent_t x_finish, diffusion_finished, rx_finish;
	unsigned int *h_overflows, *d_overflows;

	// kernel launch params
	dim3 grid_x, grid_y, grid_z, grid_r;
	dim3 threads_x, threads_y, threads_z, threads_r;

	// lattice segment geometry provided by the mapper
    SegmentDescriptor_s *segment;
};

__device__ inline size_t local_to_global(unsigned int x, unsigned int y, unsigned int z);
__device__ inline size_t local_index(unsigned int x, unsigned int y, unsigned int z);
__global__ void MPI_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int z_start, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void MPI_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void MPI_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start);
#ifdef MPD_GLOBAL_S_MATRIX
__global__ void MPI_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start,   const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG);
#else
__global__ void MPI_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const unsigned int z_start);
#endif
__global__ void mpi_correct_overflows(unsigned int* inLattice, unsigned int* siteOverflowList);

}
}

#endif
