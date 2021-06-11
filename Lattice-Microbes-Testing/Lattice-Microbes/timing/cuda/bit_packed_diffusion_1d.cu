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

#include "lm/rdme/dev/xor_random_dev.cu"
#include "lm/rdme/dev/bit_packed_diffusion_1d_dev.cu"

// Allocate the profile space.
PROF_ALLOC;

#define X_SIZE          128
#define Y_SIZE          128
#define Z_SIZE          64
#define PARTICLE_COUNT  216720      //   1 mM

__global__ void x_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int xBlockMaskMod, const unsigned int yBlockShiftMult, const unsigned int latticeXSize, const unsigned int latticeXSizeShiftMult, const unsigned int latticeXYSizeShiftMult, const unsigned int blockXSize, const unsigned int blockXSizeShiftMult, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void y_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int xBlockMaskMod, const unsigned int yBlockShiftMult, const unsigned int latticeXSize, const unsigned int latticeXSizeShiftMult, const unsigned int latticeYSize, const unsigned int latticeXYSizeShiftMult, const unsigned long long timestepHash, unsigned int* siteOverflowList);
__global__ void z_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int xBlockMaskMod, const unsigned int zBlockShiftDiv, const unsigned int latticeXSizeShiftMult, const unsigned int latticeXYSize, const unsigned int latticeXYSizeShiftMult, const unsigned int latticeZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList);
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
        startLattice = new unsigned int[X_SIZE*Y_SIZE*Z_SIZE];
        startLatticeCounts = new unsigned int[X_SIZE*Y_SIZE*Z_SIZE];
        memset(startLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*sizeof(unsigned int));
        memset(startLatticeCounts, 0, X_SIZE*Y_SIZE*Z_SIZE*sizeof(unsigned int));
        CUDA_EXCEPTION_CHECK(cudaStreamCreate(&stream));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, X_SIZE*Y_SIZE*Z_SIZE*sizeof(unsigned int)));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, X_SIZE*Y_SIZE*Z_SIZE*sizeof(unsigned int)));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&overflowList, LKCUDA_OVERFLOW_LIST_ENTRIES*sizeof(unsigned int)));

        // Fill in some random particles.
        srand(2010);
        for (unsigned int i=0; i<PARTICLE_COUNT; i++)
        {
            unsigned int r = (unsigned int)((((double)rand())/((double)RAND_MAX))*((double)X_SIZE)*((double)Y_SIZE)*((double)Z_SIZE));
            if (startLatticeCounts[r] < MPD_PARTICLE_COUNT)
            {
                startLattice[r] |= ((rand()%15)+1)<<(MPD_PARTICLE_SHIFT*startLatticeCounts[r]++);
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
            CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, startLattice, X_SIZE*Y_SIZE*Z_SIZE*sizeof(unsigned int), cudaMemcpyHostToDevice));
            CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*sizeof(unsigned int)));
            CUDA_EXCEPTION_CHECK(cudaMemset(overflowList, 0, LKCUDA_OVERFLOW_LIST_ENTRIES*sizeof(unsigned int)));

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
    const unsigned int latticeXYSize = latticeXSize*latticeYSize;
    const unsigned int latticeXSizeShiftMult = log2(latticeXSize);
    const unsigned int latticeYSizeShiftMult = log2(latticeYSize);
    const unsigned int latticeZSizeShiftMult = log2(latticeZSize);
    const unsigned int latticeXYSizeShiftMult = latticeXSizeShiftMult+latticeYSizeShiftMult;

    // Execute the kernel for the x direction.
    PROF_CUDA_BEGIN(PROF_X_DIFFUSION,stream);
    unsigned int xBlockXSize = min(MPD_X_BLOCK_MAX_X_SIZE,latticeXSize);
    unsigned int xGrid = latticeXSize/xBlockXSize;
    unsigned int yGrid = latticeYSize;
    unsigned int zGrid = latticeZSize;
    dim3 grid(xGrid*yGrid, zGrid);
    dim3 threads(xBlockXSize, 1, 1);
    CUDA_EXCEPTION_EXECUTE((x_kernel<<<grid,threads,0,stream>>>((unsigned int*)inLattice, (unsigned int*)outLattice, xGrid-1, log2(xGrid), latticeXSize, latticeXSizeShiftMult, latticeXYSizeShiftMult, xBlockXSize, log2(xBlockXSize), xseed, (unsigned int*)siteOverflowList)));
    PROF_CUDA_END(PROF_X_DIFFUSION,stream);

    // Execute the kernel for the y direction.
    PROF_CUDA_BEGIN(PROF_Y_DIFFUSION,stream);
    xGrid = latticeXSize/MPD_Y_BLOCK_X_SIZE;
    yGrid = latticeYSize/MPD_Y_BLOCK_Y_SIZE;
    zGrid = latticeZSize;
    grid.x = xGrid*yGrid;
    grid.y = zGrid;
    threads.x = MPD_Y_BLOCK_X_SIZE;
    threads.y = MPD_Y_BLOCK_Y_SIZE;
    threads.z = 1;
    CUDA_EXCEPTION_EXECUTE((y_kernel<<<grid,threads,0,stream>>>((unsigned int*)outLattice, (unsigned int*)inLattice, xGrid-1, log2(xGrid), latticeXSize, latticeXSizeShiftMult, latticeYSize, latticeXYSizeShiftMult, yseed, (unsigned int*)siteOverflowList)));
    PROF_CUDA_END(PROF_Y_DIFFUSION,stream);

    // Execute the kernel for the z direction.
    PROF_CUDA_BEGIN(PROF_Z_DIFFUSION,stream);
    xGrid = latticeXSize/MPD_Z_BLOCK_X_SIZE;
    yGrid = latticeYSize;
    zGrid = latticeZSize/MPD_Z_BLOCK_Z_SIZE;
    grid.x = xGrid*zGrid;
    grid.y = yGrid;
    threads.x = MPD_Z_BLOCK_X_SIZE;
    threads.y = 1;
    threads.z = MPD_Z_BLOCK_Z_SIZE;
    CUDA_EXCEPTION_EXECUTE((z_kernel<<<grid,threads,0,stream>>>((unsigned int*)inLattice, (unsigned int*)outLattice, xGrid-1, log2(xGrid), latticeXSizeShiftMult, latticeXYSize, latticeXYSizeShiftMult, latticeZSize, zseed, (unsigned int*)siteOverflowList)));
    PROF_CUDA_END(PROF_Z_DIFFUSION,stream);
}

__global__ void x_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int xBlockMaskMod, const unsigned int yBlockShiftMult, const unsigned int latticeXSize, const unsigned int latticeXSizeShiftMult, const unsigned int latticeXYSizeShiftMult, const unsigned int blockXSize, const unsigned int blockXSizeShiftMult, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    const unsigned int bx = blockIdx.x&xBlockMaskMod;
    const unsigned int by = blockIdx.x>>yBlockShiftMult;
    const unsigned int bz = blockIdx.y;
    const unsigned int x = threadIdx.x;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeXIndex = (bx<<blockXSizeShiftMult) + x;
    unsigned int latticeIndex = (bz<<latticeXYSizeShiftMult) + (by<<latticeXSizeShiftMult) + latticeXIndex;
    unsigned int latticeSegmentIndex = x+MPD_WINDOW_APRON_SIZE;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int latticeSegment[MPD_X_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyXWindowFromLattice(inLattice, latticeIndex, latticeXIndex, latticeXSize, latticeSegment, latticeSegmentIndex, blockXSize, timestepHash);


    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_X_WINDOW_SIZE];

    // Make the choices.
    makeXDiffusionChoices(latticeIndex, latticeXIndex, latticeXSize, latticeSegment, latticeSegmentIndex, choices, blockXSize, timestepHash);

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagation(outLattice, latticeIndex, latticeSegment, latticeSegmentIndex-1, latticeSegmentIndex, latticeSegmentIndex+1, choices, siteOverflowList);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void y_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int xBlockMaskMod, const unsigned int yBlockShiftMult, const unsigned int latticeXSize, const unsigned int latticeXSizeShiftMult, const unsigned int latticeYSize, const unsigned int latticeXYSizeShiftMult, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    const unsigned int bx = blockIdx.x&xBlockMaskMod;
    const unsigned int by = blockIdx.x>>yBlockShiftMult;
    const unsigned int bz = blockIdx.y;
    const unsigned int x = threadIdx.x;
    const unsigned int y = threadIdx.y;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by<<MPD_Y_BLOCK_Y_SIZE_SHIFT_MULT) + y;
    unsigned int latticeIndex = (bz<<latticeXYSizeShiftMult) + (latticeYIndex<<latticeXSizeShiftMult) + (bx<<MPD_Y_BLOCK_X_SIZE_SHIFT_MULT) + x;
    unsigned int latticeYSegmentIndex = y+MPD_WINDOW_APRON_SIZE;
    unsigned int latticeSegmentIndex = (latticeYSegmentIndex<<MPD_Y_BLOCK_X_SIZE_SHIFT_MULT) + x;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int latticeSegment[MPD_Y_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyYWindowFromLattice(inLattice, latticeIndex, latticeYIndex, latticeXSize, latticeYSize, latticeSegment, latticeSegmentIndex, latticeYSegmentIndex, timestepHash);


    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_Y_WINDOW_SIZE];

    // Make the choices.
    makeYDiffusionChoices(latticeIndex, latticeYIndex, latticeXSize, latticeYSize, latticeSegment, latticeSegmentIndex, latticeYSegmentIndex, choices, timestepHash);

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagation(outLattice, latticeIndex, latticeSegment, latticeSegmentIndex-MPD_Y_BLOCK_X_SIZE, latticeSegmentIndex, latticeSegmentIndex+MPD_Y_BLOCK_X_SIZE, choices, siteOverflowList);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
__global__ void z_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int xBlockMaskMod, const unsigned int zBlockShiftDiv, const unsigned int latticeXSizeShiftMult, const unsigned int latticeXYSize, const unsigned int latticeXYSizeShiftMult, const unsigned int latticeZSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
{
    const unsigned int bx = blockIdx.x&xBlockMaskMod;
    const unsigned int by = blockIdx.y;
    const unsigned int bz = blockIdx.x>>zBlockShiftDiv;
    const unsigned int x = threadIdx.x;
    const unsigned int z = threadIdx.z;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz<<MPD_Z_BLOCK_Z_SIZE_SHIFT_MULT) + z;
    unsigned int latticeIndex = (latticeZIndex<<latticeXYSizeShiftMult) + (by<<latticeXSizeShiftMult) + (bx<<MPD_Z_BLOCK_X_SIZE_SHIFT_MULT) + x;
    unsigned int latticeZSegmentIndex = z+MPD_WINDOW_APRON_SIZE;
    unsigned int latticeSegmentIndex = (latticeZSegmentIndex<<MPD_Z_BLOCK_X_SIZE_SHIFT_MULT) + x;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int latticeSegment[MPD_Z_WINDOW_SIZE];

    // Copy the x window from device memory into shared memory.
    copyZWindowFromLattice(inLattice, latticeIndex, latticeZIndex, latticeXYSize, latticeZSize, latticeSegment, latticeSegmentIndex, latticeZSegmentIndex, timestepHash);


    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ unsigned int choices[MPD_Z_WINDOW_SIZE];

    // Make the choices.
    makeZDiffusionChoices(latticeIndex, latticeZIndex, latticeXYSize, latticeZSize, latticeSegment, latticeSegmentIndex, latticeZSegmentIndex, choices, timestepHash);

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagation(outLattice, latticeIndex, latticeSegment, latticeSegmentIndex-MPD_Z_BLOCK_X_SIZE, latticeSegmentIndex, latticeSegmentIndex+MPD_Z_BLOCK_X_SIZE, choices, siteOverflowList);
}
