/*
 * University of Illinois Open Source License
 * Copyright 2010 Luthey-Schulten Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
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
 * Urbana-Champaign, nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdint.h>
#include <cuda.h>
#include "lptf/Profile.h"
#include "lm/Cuda.h"
#include "lm/Math.h"
#include "TimingConstants.h"

#define LS_WORDS_PER_SITE               2
#define LS_APRON_SIZE                   1


#if !defined LS_X_BLOCK_MAX_X_SIZE
#define LS_X_BLOCK_MAX_X_SIZE           256
#endif
#if !defined LS_Y_BLOCK_X_SIZE
#define LS_Y_BLOCK_X_SIZE               32
#endif
#if !defined LS_Y_BLOCK_Y_SIZE
#define LS_Y_BLOCK_Y_SIZE               4
#endif
#if !defined LS_Z_BLOCK_X_SIZE
#define LS_Z_BLOCK_X_SIZE               32
#endif
#if !defined LS_Z_BLOCK_Z_SIZE
#define LS_Z_BLOCK_Z_SIZE               4
#endif

#define LS_PACKED_SITES
#define LS_PACKED_LAST_OBJECT_MASK      0xFF000000


#define MPD_MAX_PARTICLE_OVERFLOWS      512
#define MPD_OVERFLOW_LIST_ENTRIES       1+2*MPD_MAX_PARTICLE_OVERFLOWS


#include "lm/rdme/dev/xor_random_dev.cu"
#include "lm/rdme/dev/lattice_sim_1d_dev.cu"
#include "lm/rdme/dev/byte_diffusion_1d_dev.cu"

// Allocate the profile space.
PROF_ALLOC;

#define X_SIZE          128
#define Y_SIZE          128
#define Z_SIZE          64
#define PARTICLE_COUNT  216720      //   1 mM

__global__ void x_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void y_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void z_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList);
void runTimestep(cudaStream_t stream, void* inLattice, void* outLattice, void* siteOverflowList, uint64_t xseed, uint64_t yseed, uint64_t zseed) throw(lm::CUDAException);

int main(int argc, char **argv)
{
    try
    {
        PROF_INIT;
        PROF_BEGIN(PROF_MAIN_RUN);

        // Allocate the cuda resources.
        cudaStream_t stream;
        unsigned int* startLattice;
        unsigned int* startLatticeCounts;
        void* inLattice;
        void* outLattice;
        void* overflowList;
        startLattice = new unsigned int[X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE];
        startLatticeCounts = new unsigned int[X_SIZE*Y_SIZE*Z_SIZE];
        memset(startLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int));
        memset(startLatticeCounts, 0, X_SIZE*Y_SIZE*Z_SIZE*sizeof(unsigned int));
        CUDA_EXCEPTION_CHECK(cudaStreamCreate(&stream));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&overflowList, MPD_OVERFLOW_LIST_ENTRIES*sizeof(unsigned int)));

        // Fill in some random particles.
        srand(2010);
        for (unsigned int i=0; i<PARTICLE_COUNT; i++)
        {
            unsigned int r = (unsigned int)((((double)rand())/((double)RAND_MAX))*((double)X_SIZE)*((double)Y_SIZE)*((double)Z_SIZE));
            if (startLatticeCounts[r] < 4)
            {
                ((unsigned char*)&startLattice[r])[startLatticeCounts[r]] = (rand()%255)+1;
                startLatticeCounts[r]++;
            }
            else if (startLatticeCounts[r] < 8)
            {
                ((unsigned char*)&startLattice[r+(X_SIZE*Y_SIZE*Z_SIZE)])[startLatticeCounts[r]] = (rand()%255)+1;
                startLatticeCounts[r]++;
            }
            else
            {
                printf("Warning: skipped adding particle to fully occupied site.\n");
            }
        }

        // Start timings the kernels.
        PROF_BEGIN(PROF_SUBMIT_KERNELS);
        PROF_CUDA_START(stream);

        // Launch the kernels.
        int NUM_LAUNCHES=100;
        for (int i=0; i<NUM_LAUNCHES; i++)
        {
            // Reset the memory.
            CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, startLattice, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyHostToDevice));
            CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
            CUDA_EXCEPTION_CHECK(cudaMemset(overflowList, 0, MPD_OVERFLOW_LIST_ENTRIES*sizeof(unsigned int)));

            // Run the timestep.
            PROF_CUDA_BEGIN(PROF_TIMESTEP_RUNNING,stream);
            runTimestep(stream, inLattice, outLattice, overflowList, 1, 2, 3);
            PROF_CUDA_END(PROF_TIMESTEP_RUNNING,stream);
        }

        // Wait for all of the kernels to finish.
        CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));

        // Record the timings.
        PROF_CUDA_FINISH(stream);
        CUDA_EXCEPTION_CHECK(cudaFree(overflowList));
        CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
        CUDA_EXCEPTION_CHECK(cudaFree(inLattice));
        delete[] startLatticeCounts;
        delete[] startLattice;
        CUDA_EXCEPTION_CHECK(cudaStreamDestroy(stream));
        PROF_END(PROF_SUBMIT_KERNELS);

        printf("Profile file saved as: %s\n",PROF_MAKE_STR(PROF_OUT_FILE));
        PROF_END(PROF_MAIN_RUN);
        PROF_WRITE;
        return 0;
    }
    catch (lm::CUDAException& e)
    {
        std::cerr << "CUDA Exception during execution: " << e.what() << std::endl;
    }
    catch (std::exception& e)
    {
        std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Unknown Exception during execution." << std::endl;
    }
    PROF_END(PROF_MAIN_RUN);
    PROF_WRITE;
    return -1;
}

void runTimestep(cudaStream_t stream, void* inLattice, void* outLattice, void* siteOverflowList, uint64_t xseed, uint64_t yseed, uint64_t zseed)
throw(lm::CUDAException)
{
    // Calculate some properties of the lattice.
    const unsigned int latticeXSize = X_SIZE;
    const unsigned int latticeYSize = Y_SIZE;
    const unsigned int latticeZSize = Z_SIZE;
    const unsigned int latticeXYSize = X_SIZE*Y_SIZE;
    const unsigned int latticeXYZSize = X_SIZE*Y_SIZE*Z_SIZE;

    // Execute the kernel for the x direction.
    PROF_CUDA_BEGIN(PROF_X_DIFFUSION,stream);
    unsigned int gridXSize;
    dim3 gridSize, threadBlockSize;
    if (!calculateXLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeXSize, latticeZSize))
        throw lm::InvalidArgException("Unable to calculate correct x launch parameters, the lattice size is incompatible.");
    CUDA_EXCEPTION_EXECUTE((x_kernel<<<gridSize,threadBlockSize,0,stream>>>((unsigned int*)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeXYSize, latticeXYZSize, xseed, (unsigned int*)siteOverflowList)));
    PROF_CUDA_END(PROF_X_DIFFUSION,stream);

    // Execute the kernel for the y direction.
    PROF_CUDA_BEGIN(PROF_Y_DIFFUSION,stream);
    if (!calculateYLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_Y_BLOCK_X_SIZE, LS_Y_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize))
        throw lm::InvalidArgException("Unable to calculate correct y launch parameters, the lattice size is incompatible.");
    CUDA_EXCEPTION_EXECUTE((y_kernel<<<gridSize,threadBlockSize,0,stream>>>((unsigned int*)outLattice, (unsigned int*)inLattice, gridXSize, latticeXSize, latticeYSize, latticeXYSize, latticeXYZSize, yseed, (unsigned int*)siteOverflowList)));
    PROF_CUDA_END(PROF_Y_DIFFUSION,stream);

    // Execute the kernel for the z direction.
    PROF_CUDA_BEGIN(PROF_Z_DIFFUSION,stream);
    if (!calculateZLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_Z_BLOCK_X_SIZE, LS_Z_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize))
        throw lm::InvalidArgException("Unable to calculate correct z launch parameters, the lattice size is incompatible.");
    CUDA_EXCEPTION_EXECUTE((z_kernel<<<gridSize,threadBlockSize,0,stream>>>((unsigned int*)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeZSize, latticeXYSize, latticeXYZSize, zseed, (unsigned int*)siteOverflowList)));
    PROF_CUDA_END(PROF_Z_DIFFUSION,stream);
}

__global__ void __launch_bounds__(LS_X_BLOCK_MAX_X_SIZE,1) x_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeXIndex = (bx*blockDim.x) + threadIdx.x;
    unsigned int latticeIndex = (bz*latticeXYSize) + (by*latticeXSize) + latticeXIndex;
    unsigned int windowIndex = threadIdx.x+LS_APRON_SIZE;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_X_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    copyXWindowFromLattice(bx, inLattice, window, latticeIndex, latticeXIndex, latticeXSize, latticeXYZSize, windowIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[LS_X_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Make the choices.
    makeXDiffusionChoices(window, choices, latticeIndex, latticeXIndex, latticeXSize, latticeXYZSize, windowIndex, blockDim.x, timestepHash);
    __syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Propagate the choices to the new lattice segment.
    performPropagation(outLattice, window, choices, latticeIndex, latticeXYZSize, windowIndex-1, windowIndex, windowIndex+1, LS_X_WINDOW_SIZE, siteOverflowList);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void __launch_bounds__(LS_Y_BLOCK_X_SIZE*LS_Y_BLOCK_Y_SIZE,1) y_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
    unsigned int latticeIndex = (bz*latticeXYSize) + (latticeYIndex*latticeXSize) + (bx*blockDim.x) + threadIdx.x;
    unsigned int windowYIndex = threadIdx.y+LS_APRON_SIZE;
    unsigned int windowIndex = (windowYIndex*blockDim.x) + threadIdx.x;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_Y_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    copyYWindowFromLattice(inLattice, window, latticeIndex, latticeYIndex, latticeXSize, latticeYSize, latticeXYZSize, windowIndex, windowYIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[LS_Y_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Make the choices.
    makeYDiffusionChoices(window, choices, latticeIndex, latticeYIndex, latticeXSize, latticeYSize, latticeXYSize, windowIndex, windowYIndex, timestepHash);
    __syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagation(outLattice, window, choices, latticeIndex, latticeXYZSize, windowIndex-LS_Y_BLOCK_X_SIZE, windowIndex, windowIndex+LS_Y_BLOCK_X_SIZE, LS_Y_WINDOW_SIZE, siteOverflowList);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void __launch_bounds__(LS_Z_BLOCK_X_SIZE*LS_Z_BLOCK_Z_SIZE,1) z_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz*blockDim.z) + threadIdx.z;
    unsigned int latticeIndex = (latticeZIndex*latticeXYSize) + (by*latticeXSize) + (bx*blockDim.x) + threadIdx.x;
    unsigned int windowZIndex = threadIdx.z+LS_APRON_SIZE;
    unsigned int windowIndex = (windowZIndex*blockDim.x) + threadIdx.x;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_Z_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    copyZWindowFromLattice(inLattice, window, latticeIndex, latticeZIndex, latticeZSize, latticeXYSize, latticeXYZSize, windowIndex, windowZIndex);
    __syncthreads();

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[LS_Z_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Make the choices.
    makeZDiffusionChoices(window, choices, latticeIndex, latticeZIndex, latticeZSize, latticeXYSize, latticeXYZSize, windowIndex, windowZIndex, timestepHash);
    __syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagation(outLattice, window, choices, latticeIndex, latticeXYZSize, windowIndex-LS_Z_BLOCK_X_SIZE, windowIndex, windowIndex+LS_Z_BLOCK_X_SIZE, LS_Z_WINDOW_SIZE, siteOverflowList);
}
