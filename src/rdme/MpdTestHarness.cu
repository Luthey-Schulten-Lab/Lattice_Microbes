/*
 * University of Illinois Open Source License
 * Copyright 2015-2018 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.uiuc.edu/~schulten
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
#include <sstream>
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
#include "rdme/CudaByteLattice.h"
#include "rdme/MpdRdmeSolver.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include "rdme/MpdTestHarness.h"

#include <nvrtc.h>

#define MPD_WORDS_PER_SITE              2
#define MPD_APRON_SIZE                  1

using std::map;
using lm::io::DiffusionModel;
using lm::rdme::Lattice;
using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {

__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) sanity_check(const unsigned int* L1, const unsigned int* L2);

extern __global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList
#ifdef MPD_GLOBAL_S_MATRIX
 , const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const float* qp0, const float* qp1, const float* qp2
#endif
);

MpdTestHarness::MpdTestHarness()
:MpdRdmeSolver()
{
}

MpdTestHarness::~MpdTestHarness()
{
}


void MpdTestHarness::generateTrajectory()
{
    // Shadow the lattice member as a cuda lattice.
    CudaByteLattice * lattice = (CudaByteLattice *)this->lattice;

    // Synchronize the cuda memory.
    lattice->copyToGPU();

    // Get the interval for writing species counts and lattices.
    // Get the simulation time limit.
    double maxTime=atof((*parameters)["maxTime"].c_str());

    Print::printf(Print::INFO, "Running harness with %d species, %d reactions, %d site types for %e s with tau %e.", numberSpecies, numberReactions, numberSiteTypes, maxTime, tau);

    // Set the initial time.
    double time = 0.0;
    uint32_t timestep=1;

	total_orig = total_jit = 0.0f;

	cudaEventCreate(&original_start);
	cudaEventCreate(&original_end);
	cudaEventCreate(&jit_start);
	cudaEventCreate(&jit_end);
    // Loop until we have finished the simulation.
    while (time < maxTime)
    {
        // Run the next timestep.
        runTimestep(lattice, timestep++);

        // Update the time.
        time += tau;

    }

	timestep--;
	float avg_orig = total_orig / timestep;
	float avg_jit = total_jit / timestep;
	printf("FINAL Steps %d avg_orig %f avg_alt %f\n", timestep, avg_orig, avg_jit);
}

void MpdTestHarness::runTimestep(CudaByteLattice * lattice, uint32_t timestep)
{
//	printf("*@ Timestep %d @\n", timestep);

    // Calculate some properties of the lattice.
    lattice_coord_t size = lattice->getSize();
    const unsigned int latticeXSize = size.x;
    const unsigned int latticeYSize = size.y;
    const unsigned int latticeZSize = size.z;

    dim3 gridSize, threadBlockSize;

    // Execute the kernel for the x direction.


    calculateXLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((MpdRdmeSolver_x_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,0), (unsigned int*)cudaOverflowList)));
    lattice->swapSrcDest();

    // Execute the kernel for the y direction.
    calculateYLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_Y_BLOCK_X_SIZE, TUNE_MPD_Y_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((MpdRdmeSolver_y_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,1), (unsigned int*)cudaOverflowList)));
    lattice->swapSrcDest();

    // Execute the kernel for the z direction.
    calculateZLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_Z_BLOCK_X_SIZE, TUNE_MPD_Z_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((MpdRdmeSolver_z_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,2), (unsigned int*)cudaOverflowList)));
    lattice->swapSrcDest();

	// Execute the kernel for the reaction, this kernel updates the lattice in-place, so only the src pointer is passed.
	calculateReactionLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_REACTION_BLOCK_X_SIZE, TUNE_MPD_REACTION_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);

	// Copy SRC to DST so we can co-run
	CUDA_EXCEPTION_EXECUTE(cudaMemcpyAsync(lattice->getGPUMemoryDest(), lattice->getGPUMemorySrc(), lattice->getParticleMemorySize(), cudaMemcpyDeviceToDevice, cudaStream));

	// Run original

	// Run static alternate test kernel

	cudaEventRecord(original_start, cudaStream);
	CUDA_EXCEPTION_EXECUTE((MpdRdmeSolver_reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList
#ifdef MPD_GLOBAL_S_MATRIX
, SG, RLG
#endif
)));
	cudaEventRecord(original_end, cudaStream);

	cudaEventRecord(jit_start, cudaStream);
	// Defined in MpdRdmeSolver.cu to make sure we get same constants
	CUDA_EXCEPTION_EXECUTE((precomp_reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemoryDest(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList
#ifdef MPD_GLOBAL_S_MATRIX
, SG, RLG, propZeroOrder, propFirstOrder, propSecondOrder
#endif
)));

	cudaEventRecord(jit_end, cudaStream);

	// Check sanity
	sanity_check<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (unsigned int *)lattice->getGPUMemoryDest());

    // Wait for the kernels to complete.
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(cudaStream));

	float otime, jtime;
	cudaEventElapsedTime(&otime, original_start, original_end);
	cudaEventElapsedTime(&jtime, jit_start, jit_end);
	if (timestep % 5000 == 0) printf("TS %d Orignial time: %f ms, Alt time: %f ms\n", timestep, otime, jtime);
	total_orig += otime;
	total_jit  += jtime;

    uint32_t overflowList[1+2*TUNE_MPD_MAX_PARTICLE_OVERFLOWS];
    CUDA_EXCEPTION_CHECK(cudaMemcpy(overflowList, cudaOverflowList, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    uint numberExceptions = overflowList[0];
    if (numberExceptions > 0)
    {
        Print::printf(Print::DEBUG, "%d overflows (not resolving)", numberExceptions);

        // Reset the overflow list.
        CUDA_EXCEPTION_CHECK(cudaMemset(cudaOverflowList, 0, sizeof(uint32_t)));
    }
}


}
}
