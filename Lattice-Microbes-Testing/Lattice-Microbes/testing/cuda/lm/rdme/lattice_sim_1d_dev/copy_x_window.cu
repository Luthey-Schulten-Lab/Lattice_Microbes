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
#define LS_X_BLOCK_MAX_X_SIZE           128
#define LS_Y_BLOCK_X_SIZE               16
#define LS_Y_BLOCK_Y_SIZE               8
#define LS_Z_BLOCK_X_SIZE               16
#define LS_Z_BLOCK_Z_SIZE               8
#define LS_BOUNDARY_VALUE               0xFFEEDDCC

#include "lm/rdme/dev/lattice_sim_1d_dev.cu"

__global__ void cu_CopyXWindowSites_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYZSize);
__global__ void cu_CopyXWindowAprons_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYZSize);

void cu_CopyXWindowSites(unsigned int * host_inLattice, unsigned int * host_outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize)
{
    void* inLattice;
    void* outLattice;
    CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, host_inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0xFF, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));

    unsigned int gridXSize;
    dim3 gridSize, threadBlockSize;
    calculateXLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((cu_CopyXWindowSites_kernel<<<gridSize,threadBlockSize>>>((unsigned int*)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeYSize, latticeXSize*latticeYSize*latticeZSize)));

    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(0));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(host_outLattice, outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyDeviceToHost));

    CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
    CUDA_EXCEPTION_CHECK(cudaFree(inLattice));
}

__global__ void cu_CopyXWindowSites_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYZSize)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeXIndex = (bx*blockDim.x) + threadIdx.x;
    unsigned int latticeIndex = (bz*latticeXSize*latticeYSize) + (by*latticeXSize) + latticeXIndex;
    unsigned int windowIndex = threadIdx.x+LS_APRON_SIZE;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_X_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    copyXWindowFromLattice(bx, inLattice, window, latticeIndex, latticeXIndex, latticeXSize, latticeXYZSize, windowIndex);

    // Copy the x window from shared memory to device memory.
    copyXWindowToLattice(outLattice, window, latticeIndex, latticeXIndex, latticeXSize, latticeXYZSize, windowIndex);
}

void cu_CopyXWindowAprons(unsigned int * host_inLattice, unsigned int * host_outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize)
{
    void* inLattice;
    void* outLattice;
    CUDA_EXCEPTION_CHECK(cudaMalloc(&inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(inLattice, host_inLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0xFF, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int)));

    unsigned int gridXSize;
    dim3 gridSize, threadBlockSize;
    calculateXLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((cu_CopyXWindowAprons_kernel<<<gridSize,threadBlockSize>>>((unsigned int*)inLattice, (unsigned int*)outLattice, gridXSize, latticeXSize, latticeYSize, latticeXSize*latticeYSize*latticeZSize)));

    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(0));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(host_outLattice, outLattice, latticeXSize*latticeYSize*latticeZSize*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyDeviceToHost));

    CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
    CUDA_EXCEPTION_CHECK(cudaFree(inLattice));
}

__global__ void cu_CopyXWindowAprons_kernel(const unsigned int* inLattice, unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYZSize)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeXIndex = (bx*blockDim.x) + threadIdx.x;
    unsigned int latticeIndex = (bz*latticeXSize*latticeYSize) + (by*latticeXSize) + latticeXIndex;
    unsigned int windowIndex = threadIdx.x+LS_APRON_SIZE;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[LS_X_WINDOW_SIZE*LS_WORDS_PER_SITE];

    // Copy the x window from device memory into shared memory.
    copyXWindowFromLattice(bx, inLattice, window, latticeIndex, latticeXIndex, latticeXSize, latticeXYZSize, windowIndex);

    __syncthreads();

    if (latticeXIndex < latticeXSize)
    {
        outLattice[latticeIndex] = 0;
        outLattice[latticeIndex+latticeXYZSize] = 0;

        // If this is the first part of the block, load the leading apron.
        if (windowIndex < 2*LS_APRON_SIZE)
        {
            outLattice[latticeIndex] = window[windowIndex-LS_APRON_SIZE];
            outLattice[latticeIndex+latticeXYZSize] = window[windowIndex-LS_APRON_SIZE+LS_X_WINDOW_SIZE];
        }

        // If this is the last part of the block, load the trailing apron.
        if (windowIndex >= blockDim.x)
        {
            outLattice[latticeIndex] = window[windowIndex+LS_APRON_SIZE];
            outLattice[latticeIndex+latticeXYZSize] = window[windowIndex+LS_APRON_SIZE+LS_X_WINDOW_SIZE];
        }
    }
}
