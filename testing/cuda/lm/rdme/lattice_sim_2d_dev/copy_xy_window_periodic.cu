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
#define LS_XY_BLOCK_X_SIZE              16
#define LS_XY_BLOCK_Y_SIZE              8
#define LS_Z_BLOCK_X_SIZE               8
#define LS_Z_BLOCK_Z_SIZE               8
#define LS_BOUNDARY_PERIODIC            1

#include "lm/rdme/dev/lattice_sim_2d_dev.cu"

__global__ void cu_CopyXYWindowPeriodicSites_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize);
__global__ void cu_CopyXYWindowPeriodicAprons_kernel(int mode, const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize);


void cu_CopyXYWindowPeriodicSites(unsigned int * host_inLattice, unsigned int * host_outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize)
{
    void* inLattice;
    void* outLattice;
    CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, host_inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0xFF, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));

    unsigned int gridXSize;
    dim3 gridSize, threadBlockSize;
    if (!calculateXYLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_XY_BLOCK_X_SIZE, LS_XY_BLOCK_Y_SIZE, LS_APRON_SIZE, latticeXSize, latticeYSize, latticeZSize))
        throw lm::InvalidArgException("Unable to calculate correct launch parameters, the lattice, block, and thread sizes are incompatible.");
    CUDA_EXCEPTION_EXECUTE((cu_CopyXYWindowPeriodicSites_kernel<<<gridSize,threadBlockSize>>>((unsigned int *)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeYSize, latticeXSize*latticeYSize, latticeXSize*latticeYSize*latticeZSize)));

    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(0));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(host_outLattice, outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyDeviceToHost));

    CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
    CUDA_EXCEPTION_CHECK(cudaFree(inLattice));
}


__global__ void cu_CopyXYWindowPeriodicSites_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize)
{
    __shared__ unsigned int bx, by, bz, gx, gy;
    calculateXYBlockIndices(&bx, &by, &bz, &gx, &gy, gridXSize);

    // Figure out the indices for this thread in the block, lattice, and window.
    int blockXIndex, blockYIndex, latticeXIndex, latticeYIndex;
    unsigned int latticeIndex,  windowXIndex,  windowYIndex,  windowIndex;
    calculateXYThreadIndices(bx, by, bz, latticeXSize, latticeYSize, latticeXYSize, &blockXIndex, &blockYIndex, &latticeXIndex, &latticeYIndex, &latticeIndex,  &windowXIndex,  &windowYIndex,  &windowIndex);

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_XY_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    copyXYWindowFromLattice(inLattice, window, latticeIndex, latticeYIndex, latticeXSize, latticeYSize, latticeXYZSize, windowIndex, windowXIndex, bx, gx);


    // Copy the x window from shared memory to device memory.
    copyXYWindowToLattice(outLattice, window, latticeIndex, latticeXYZSize, windowIndex, blockYIndex);
}


void cu_CopyXYWindowPeriodicAprons(int mode, unsigned int * host_inLattice, unsigned int * host_outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize)
{
    void* inLattice;
    void* outLattice;
    CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, host_inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0xFF, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));

    unsigned int gridXSize;
    dim3 gridSize, threadBlockSize;
    if (!calculateXYLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_XY_BLOCK_X_SIZE, LS_XY_BLOCK_Y_SIZE, LS_APRON_SIZE, latticeXSize, latticeYSize, latticeZSize))
        throw lm::InvalidArgException("Unable to calculate correct launch parameters, the lattice, block, and thread sizes are incompatible.");
    CUDA_EXCEPTION_EXECUTE((cu_CopyXYWindowPeriodicAprons_kernel<<<gridSize,threadBlockSize>>>(mode, (unsigned int*)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeYSize, latticeXSize*latticeYSize, latticeXSize*latticeYSize*latticeZSize)));

    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(0));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(host_outLattice, outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyDeviceToHost));

    CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
    CUDA_EXCEPTION_CHECK(cudaFree(inLattice));
}

__global__ void cu_CopyXYWindowPeriodicAprons_kernel(int mode, const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize)
{
    __shared__ unsigned int bx, by, bz, gx, gy;
    calculateXYBlockIndices(&bx, &by, &bz, &gx, &gy, gridXSize);

    // Figure out the indices for this thread in the block, lattice, and window.
    int blockXIndex, blockYIndex, latticeXIndex, latticeYIndex;
    unsigned int latticeIndex,  windowXIndex,  windowYIndex,  windowIndex;
    calculateXYThreadIndices(bx, by, bz, latticeXSize, latticeYSize, latticeXYSize, &blockXIndex, &blockYIndex, &latticeXIndex, &latticeYIndex, &latticeIndex,  &windowXIndex,  &windowYIndex,  &windowIndex);

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_XY_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    copyXYWindowFromLattice(inLattice, window, latticeIndex, latticeYIndex, latticeXSize, latticeYSize, latticeXYZSize, windowIndex, windowXIndex, bx, gx);
    __syncthreads();

    // Put the apron into the lattice block
    if (blockXIndex >= 0 && blockXIndex < LS_XY_BLOCK_X_SIZE && blockYIndex >= 0 && blockYIndex < LS_XY_BLOCK_Y_SIZE)
    {
        #if LS_WORDS_PER_SITE >= 2
        unsigned int* window2 = window+LS_XY_WINDOW_SIZE;
        #endif

        // -x,-y
        if (mode == 1 && blockXIndex < LS_APRON_SIZE && blockYIndex < LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex-LS_APRON_SIZE-(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex-LS_APRON_SIZE-(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // -y
        else if (mode == 2 && blockYIndex < LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex-(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex-(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // +x,-y
        else if (mode == 3 && blockXIndex >= LS_XY_BLOCK_X_SIZE-LS_APRON_SIZE && blockYIndex < LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex+LS_APRON_SIZE-(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex+LS_APRON_SIZE-(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // -x
        else if (mode == 4 && blockXIndex < LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex-LS_APRON_SIZE];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex-LS_APRON_SIZE];
            #endif
        }

        // +x
        else if (mode == 5 && blockXIndex >= LS_XY_BLOCK_X_SIZE-LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex+LS_APRON_SIZE];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex+LS_APRON_SIZE];
            #endif
        }

        // -x,+y
        else if (mode == 6 && blockXIndex < LS_APRON_SIZE && blockYIndex >= LS_XY_BLOCK_Y_SIZE-LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex-LS_APRON_SIZE+(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex-LS_APRON_SIZE+(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // +y
        else if (mode == 7 && blockYIndex >= LS_XY_BLOCK_Y_SIZE-LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex+(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex+(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        // +x,+y
        else if (mode == 8 && blockXIndex >= LS_XY_BLOCK_X_SIZE-LS_APRON_SIZE && blockYIndex >= LS_XY_BLOCK_Y_SIZE-LS_APRON_SIZE)
        {
            window[windowIndex]=window[windowIndex+LS_APRON_SIZE+(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=window2[windowIndex+LS_APRON_SIZE+(LS_XY_WINDOW_X_SIZE*LS_APRON_SIZE)];
            #endif
        }

        else
        {
            window[windowIndex]=0;
            #if LS_WORDS_PER_SITE >= 2
            window2[windowIndex]=0;
            #endif
        }
    }
    __syncthreads();

    // Copy the x window from shared memory to device memory.
    copyXYWindowToLattice(outLattice, window, latticeIndex, latticeXYZSize, windowIndex, blockYIndex);
}
