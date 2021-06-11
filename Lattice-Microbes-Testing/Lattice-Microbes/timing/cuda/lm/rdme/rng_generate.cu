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
#include <cstdio>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#define PROF_ENABLE
#define PROF_MAX_THREADS                1
#define PROF_MAX_EVENTS                 1000
#define PROF_MAX_CUDA_EVENT_BUFFER      1000
#include "lptf/Profile.h"
#undef PROF_ENABLE
#include "lm/Cuda.h"
#include "TimingConstants.h"

#define LS_WORDS_PER_SITE               2
#define LS_APRON_SIZE                   1
#define LS_X_BLOCK_MAX_X_SIZE           256
#define LS_Y_BLOCK_X_SIZE               32
#define LS_Y_BLOCK_Y_SIZE               4
#define LS_Z_BLOCK_X_SIZE               32
#define LS_Z_BLOCK_Z_SIZE               4
#include "lm/rdme/dev/xor_random_dev.cu"
#include "lm/rdme/dev/lattice_sim_1d_dev.cu"

#include <curand_kernel.h>


// Allocate the profile space.
PROF_ALLOC;

#define X_SIZE          128
#define Y_SIZE          128
#define Z_SIZE          64
#define NUM_LAUNCHES    100

__global__ void xorshift_int_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash);
__global__ void xorshift_float_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash);
__global__ void xorwow_init_kernel(unsigned int* outLattice, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, curandState *rngState);
__global__ void xorwow_int_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, curandState *state);
__global__ void xorwow_float_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, curandState *state);

int main(int argc, char **argv)
{
    try
    {
        PROF_INIT;
        PROF_BEGIN(PROF_MAIN_RUN);

        // Allocate the cuda resources.
        cudaStream_t stream;
        unsigned int * hostOutLattice = new unsigned int[X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE];
        void* outLattice;
        void* rngState;
        CUDA_EXCEPTION_CHECK(cudaStreamCreate(&stream));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&outLattice, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&rngState, X_SIZE*Y_SIZE*Z_SIZE*sizeof(curandState)));

        // Start timings the kernels.
        PROF_BEGIN(PROF_SUBMIT_KERNELS);

        // Calculate some properties of the lattice.
        const unsigned int latticeXSize = X_SIZE;
        const unsigned int latticeYSize = Y_SIZE;
        const unsigned int latticeZSize = Z_SIZE;
        const unsigned int latticeXYSize = X_SIZE*Y_SIZE;
        const unsigned int latticeXYZSize = X_SIZE*Y_SIZE*Z_SIZE;
        unsigned int gridXSize;
        dim3 gridSize, threadBlockSize;
        if (!calculateXLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, LS_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize))
            throw lm::InvalidArgException("Unable to calculate correct launch parameters, the lattice size is incompatible.");

        // Launch the xorshift kernels.
        PROF_CUDA_START(stream);
        for (int i=0; i<NUM_LAUNCHES; i++)
        {
            // Execute the kernel.
            unsigned long long hash = (i+1);
            CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
            PROF_CUDA_BEGIN(PROF_XORSHIFT_INT,stream);
            CUDA_EXCEPTION_EXECUTE((xorshift_int_kernel<<<gridSize,threadBlockSize,0,stream>>>((unsigned int *)outLattice, gridXSize, latticeXSize, latticeXYSize, latticeXYZSize, hash)));
            PROF_CUDA_END(PROF_XORSHIFT_INT,stream);
        }
        for (int i=0; i<NUM_LAUNCHES; i++)
        {
            // Execute the kernel.
            unsigned long long hash = (i+1);
            CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
            PROF_CUDA_BEGIN(PROF_XORSHIFT_FLOAT,stream);
            CUDA_EXCEPTION_EXECUTE((xorshift_float_kernel<<<gridSize,threadBlockSize,0,stream>>>((unsigned int *)outLattice, gridXSize, latticeXSize, latticeXYSize, latticeXYZSize, hash)));
            PROF_CUDA_END(PROF_XORSHIFT_FLOAT,stream);
        }
        CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));
        PROF_CUDA_FINISH(stream);

        // Initialize the xorwow state.
        //printf("State init\n");
        CUDA_EXCEPTION_CHECK(cudaThreadSetLimit(cudaLimitStackSize, 16384));
        PROF_CUDA_START(stream);
        CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
        CUDA_EXCEPTION_CHECK(cudaMemset(rngState, 0, X_SIZE*Y_SIZE*Z_SIZE*sizeof(curandState)));
        PROF_CUDA_BEGIN(PROF_XORWOW_INIT,stream);
        dim3 gridSizeInit(X_SIZE/32,1,1), threadBlockSizeInit(32,1,1);
        CUDA_EXCEPTION_EXECUTE((xorwow_init_kernel<<<gridSizeInit,threadBlockSizeInit,0,stream>>>((unsigned int *)outLattice, latticeXSize, latticeYSize, latticeZSize, (curandState *)rngState)));
        PROF_CUDA_END(PROF_XORWOW_INIT,stream);
        CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));
        PROF_CUDA_FINISH(stream);
        CUDA_EXCEPTION_CHECK(cudaThreadSetLimit(cudaLimitStackSize, 1024));
        //printf("Done state init\n");

        // Make sure the init was done correctly.
        CUDA_EXCEPTION_CHECK(cudaMemcpy(hostOutLattice, outLattice, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int), cudaMemcpyDeviceToHost));
        /*
        for (int z=0, i=0; z<3; z++)
        {
            printf("%2d--------------------------------------------------------------\n",z);
            for (int y=0; y<Y_SIZE; y++)
            {
                for (int x=0; x<X_SIZE; x++, i++)
                {
                    printf("%d ",hostOutLattice[i]);
                }
                printf("\n");
            }
        }
        int totalL=0;
        int totalU=0;
        for (int z=0, i=0; z<Z_SIZE; z++)
        {
            for (int y=0; y<Y_SIZE; y++)
            {
                for (int x=0; x<X_SIZE; x++, i++)
                {
                    totalL += hostOutLattice[i];
                    totalU += hostOutLattice[i+latticeXYZSize];
                }
            }
        }
        printf("Total initialized sites: %d (lower), %d (upper)\n",totalL, totalU);
        */

        // Launch the xorwow kernels.
        PROF_CUDA_START(stream);
        for (int i=0; i<NUM_LAUNCHES; i++)
        {
            // Execute the kernel.
            CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
            PROF_CUDA_BEGIN(PROF_XORWOW_INT,stream);
            CUDA_EXCEPTION_EXECUTE((xorwow_int_kernel<<<gridSize,threadBlockSize,0,stream>>>((unsigned int *)outLattice, gridXSize, latticeXSize, latticeXYSize, latticeXYZSize, (curandState *)rngState)));
            PROF_CUDA_END(PROF_XORWOW_INT,stream);
        }
        for (int i=0; i<NUM_LAUNCHES; i++)
        {
            // Execute the kernel.
            CUDA_EXCEPTION_CHECK(cudaMemset(outLattice, 0, X_SIZE*Y_SIZE*Z_SIZE*LS_WORDS_PER_SITE*sizeof(unsigned int)));
            PROF_CUDA_BEGIN(PROF_XORWOW_FLOAT,stream);
            CUDA_EXCEPTION_EXECUTE((xorwow_float_kernel<<<gridSize,threadBlockSize,0,stream>>>((unsigned int *)outLattice, gridXSize, latticeXSize, latticeXYSize, latticeXYZSize, (curandState *)rngState)));
            PROF_CUDA_END(PROF_XORWOW_FLOAT,stream);
        }
        CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));
        PROF_CUDA_FINISH(stream);

        // Free any resources.
        CUDA_EXCEPTION_CHECK(cudaFree(rngState));
        CUDA_EXCEPTION_CHECK(cudaFree(outLattice));
        CUDA_EXCEPTION_CHECK(cudaStreamDestroy(stream));
        delete[] hostOutLattice;
        PROF_END(PROF_SUBMIT_KERNELS);

        printf("Profile file saved as: %s\n",PROF_MAKE_STR(PROF_OUT_FILE));
        PROF_END(PROF_MAIN_RUN);
        PROF_WRITE;
        return 0;
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

__global__ void xorshift_int_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz*blockDim.z) + threadIdx.z;
    unsigned int latticeIndex = (latticeZIndex*latticeXYSize) + (by*latticeXSize) + (bx*blockDim.x) + threadIdx.x;

    int sum=0;
    for (int i=0; i<16; i++)
    {
        unsigned int randomValue = getRandomHash(latticeIndex, 4, i, timestepHash);
        sum += (randomValue>2147483648)?1:0;
    }

    outLattice[latticeIndex] = sum;
}


__global__ void xorshift_float_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned long long timestepHash)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz*blockDim.z) + threadIdx.z;
    unsigned int latticeIndex = (latticeZIndex*latticeXYSize) + (by*latticeXSize) + (bx*blockDim.x) + threadIdx.x;

    int sum=0;
    for (int i=0; i<16; i++)
    {
        float randomValue = getRandomHashFloat(latticeIndex, 4, i, timestepHash);
        sum += (randomValue>0.5)?1:0;
    }

    outLattice[latticeIndex] = sum;
}


__global__ void xorwow_init_kernel(unsigned int* outLattice, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, curandState *rngState)
{
    int latticeXIndex = (blockIdx.x*blockDim.x)+threadIdx.x;

    for (int latticeZIndex=0; latticeZIndex<latticeZSize; latticeZIndex++)
    {
        for (int latticeYIndex=0; latticeYIndex<latticeYSize; latticeYIndex++)
        {
            unsigned int latticeIndex = (latticeZIndex*latticeXSize*latticeYSize) + (latticeYIndex*latticeXSize) + latticeXIndex;
            curand_init(1234, latticeIndex, 0, &rngState[latticeIndex]);
            outLattice[latticeIndex] += 1;
        }
    }
}

__global__ void xorwow_int_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, curandState *rngState)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz*blockDim.z) + threadIdx.z;
    unsigned int latticeIndex = (latticeZIndex*latticeXYSize) + (by*latticeXSize) + (bx*blockDim.x) + threadIdx.x;

    curandState localRngState = rngState[latticeIndex];
    int sum=0;
    for (int i=0; i<16; i++)
    {
        unsigned int randomValue = curand(&localRngState);
        sum += (randomValue>2147483648)?1:0;
    }
    rngState[latticeIndex] = localRngState;

    outLattice[latticeIndex] = sum;
}

__global__ void xorwow_float_kernel(unsigned int* outLattice, const unsigned int gridXSize, const unsigned int latticeXSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, curandState *rngState)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz*blockDim.z) + threadIdx.z;
    unsigned int latticeIndex = (latticeZIndex*latticeXYSize) + (by*latticeXSize) + (bx*blockDim.x) + threadIdx.x;

    curandState localRngState = rngState[latticeIndex];
    int sum=0;
    for (int i=0; i<16; i++)
    {
        float randomValue = curand_uniform(&localRngState);
        sum += (randomValue>0.5)?1:0;
    }
    rngState[latticeIndex] = localRngState;

    outLattice[latticeIndex] = sum;
}
