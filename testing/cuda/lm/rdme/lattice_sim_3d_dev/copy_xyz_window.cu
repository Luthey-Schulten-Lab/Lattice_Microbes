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

#include "lm/Cuda.h"

#define LS_WORDS_PER_SITE               2
#define LS_APRON_SIZE                   2
#define LS_XYZ_BLOCK_X_SIZE             16
#define LS_XYZ_BLOCK_Y_SIZE             8
#define LS_XYZ_BLOCK_Z_SIZE             4
#define LS_XYZ_Z_THREADS                2
#define LS_BOUNDARY_VALUE               0xFFEEDDCC

#include "lm/rdme/dev/lattice_sim_3d_dev.cu"

__global__ void cu_CopyXYZWindowSites_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize);
__global__ void cu_CopyXYZWindowAprons_kernel(int mode, const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize);

void cu_CopyXYZWindowSites(unsigned int * host_inLattice, unsigned int * host_outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize)
{
    void* inLattice;
    void* outLattice;
    CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, host_inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0xFF, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));

    unsigned int gridXSize;
    dim3 gridSize, threadBlockSize;
    if (!calculateLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_XYZ_BLOCK_X_SIZE, LS_XYZ_BLOCK_Y_SIZE, LS_XYZ_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize, LS_XYZ_X_THREADS, LS_XYZ_Y_THREADS, LS_XYZ_Z_THREADS))
        throw lm::InvalidArgException("Unable to calculate correct launch parameters, the lattice, block, and thread sizes are incompatible.");
    CUDA_EXCEPTION_EXECUTE((cu_CopyXYZWindowSites_kernel<<<gridSize,threadBlockSize>>>((unsigned int*)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeYSize, latticeZSize, latticeXSize*latticeYSize, latticeXSize*latticeYSize*latticeZSize)));

    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(0));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(host_outLattice, outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyDeviceToHost));

    CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
    CUDA_EXCEPTION_CHECK(cudaFree(inLattice));
}

__global__ void cu_CopyXYZWindowSites_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize)
{
    __shared__ unsigned int bx, by, bz, gx, gy;
    calculateBlockIndices(&bx, &by, &bz, &gx, &gy, gridXSize);

    // Figure out the indices of this thread.
    unsigned int latticeXIndex, latticeYIndex, latticeZIndex, latticeIndex;
    unsigned int windowXIndex, windowYIndex, windowZIndex, windowIndex;
    calculateThreadIndices(bx, by, bz, latticeXSize, latticeXYSize, &latticeXIndex, &latticeYIndex, &latticeZIndex, &latticeIndex, &windowXIndex, &windowYIndex, &windowZIndex, &windowIndex);

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_XYZ_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    //copyXYZWindowFromLatticeStriped(inLattice, window, latticeXSize, latticeXYSize, latticeXYZSize, bx, by, bz, gx, gy);
    //copyXYZWindowFromLatticeQuadrant(inLattice, window, latticeIndex, latticeXSize, latticeXYSize, latticeXYZSize, windowIndex, bx, by, gx, gy);
    __syncthreads();

    // Copy the x window from shared memory to device memory.
    copyXYZWindowToLattice(outLattice, window, latticeIndex, latticeXYSize, latticeXYZSize, windowIndex);
}

void cu_CopyXYZWindowAprons(int mode, unsigned int * host_inLattice, unsigned int * host_outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize)
{
    void* inLattice;
    void* outLattice;
    CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, host_inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0xFF, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));

    unsigned int gridXSize;
    dim3 gridSize, threadBlockSize;
    if (!calculateLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_XYZ_BLOCK_X_SIZE, LS_XYZ_BLOCK_Y_SIZE, LS_XYZ_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize, LS_XYZ_X_THREADS, LS_XYZ_Y_THREADS, LS_XYZ_Z_THREADS))
        throw lm::InvalidArgException("Unable to calculate correct launch parameters, the lattice, block, and thread sizes are incompatible.");
    CUDA_EXCEPTION_EXECUTE((cu_CopyXYZWindowAprons_kernel<<<gridSize,threadBlockSize>>>(mode, (unsigned int*)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeYSize, latticeZSize, latticeXSize*latticeYSize, latticeXSize*latticeYSize*latticeZSize)));

    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(0));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(host_outLattice, outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyDeviceToHost));

    CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
    CUDA_EXCEPTION_CHECK(cudaFree(inLattice));
}

__global__ void cu_CopyXYZWindowAprons_kernel(int mode, const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize)
{
    __shared__ unsigned int bx, by, bz, gx, gy;
    calculateBlockIndices(&bx, &by, &bz, &gx, &gy, gridXSize);

    // Figure out the indices of this thread.
    unsigned int latticeXIndex, latticeYIndex, latticeZIndex, latticeIndex;
    unsigned int windowXIndex, windowYIndex, windowZIndex, windowIndex;
    calculateThreadIndices(bx, by, bz, latticeXSize, latticeXYSize, &latticeXIndex, &latticeYIndex, &latticeZIndex, &latticeIndex, &windowXIndex, &windowYIndex, &windowZIndex, &windowIndex);

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_XYZ_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    //copyXYZWindowFromLatticeStriped(inLattice, window, latticeXSize, latticeXYSize, latticeXYZSize, bx, by, bz, gx, gy);
    copyXYZWindowFromLatticeQuadrant(inLattice, window, latticeIndex, latticeXSize, latticeXYSize, latticeXYZSize, windowIndex, bx, by, gx, gy);
    __syncthreads();

    #if LS_WORDS_PER_SITE >= 2
    unsigned int* window2 = window+LS_XYZ_WINDOW_SIZE;
    #endif

    // Load the block.
    for (int i=0; i<LS_XYZ_Z_LOOPS; i++)
    {
        unsigned int loopWindowIndex = LS_Z_LOOP_WINDOW_INDEX(windowIndex,i);

        // -x,-y
        if (mode == 1 && threadIdx.x < LS_APRON_SIZE && threadIdx.y < LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex-LS_APRON_SIZE-(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex-LS_APRON_SIZE-(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // -y
        else if (mode == 2 && threadIdx.y < LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex-(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex-(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // +x,-y
        else if (mode == 3 && threadIdx.x >= LS_XYZ_BLOCK_X_SIZE-LS_APRON_SIZE && threadIdx.y < LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex+LS_APRON_SIZE-(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex+LS_APRON_SIZE-(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // -x
        else if (mode == 4 && threadIdx.x < LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex-LS_APRON_SIZE];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex-LS_APRON_SIZE];
            #endif
        }

        // +x
        else if (mode == 5 && threadIdx.x >= LS_XYZ_BLOCK_X_SIZE-LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex+LS_APRON_SIZE];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex+LS_APRON_SIZE];
            #endif
        }

        // -x,+y
        else if (mode == 6 && threadIdx.x < LS_APRON_SIZE && threadIdx.y >= LS_XYZ_BLOCK_Y_SIZE-LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex-LS_APRON_SIZE+(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex-LS_APRON_SIZE+(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // +y
        else if (mode == 7 && threadIdx.y >= LS_XYZ_BLOCK_Y_SIZE-LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex+(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex+(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // +x,+y
        else if (mode == 8 && threadIdx.x >= LS_XYZ_BLOCK_X_SIZE-LS_APRON_SIZE && threadIdx.y >= LS_XYZ_BLOCK_Y_SIZE-LS_APRON_SIZE)
        {
            window[loopWindowIndex]=window[loopWindowIndex+LS_APRON_SIZE+(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=window2[loopWindowIndex+LS_APRON_SIZE+(LS_XYZ_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        else
        {
            window[loopWindowIndex]=0;
            #if LS_WORDS_PER_SITE >= 2
            window2[loopWindowIndex]=0;
            #endif
        }
    }
    __syncthreads();


    // Copy the x window from shared memory to device memory.
    copyXYZWindowToLattice(outLattice, window, latticeIndex, latticeXYSize, latticeXYZSize, windowIndex);
}
