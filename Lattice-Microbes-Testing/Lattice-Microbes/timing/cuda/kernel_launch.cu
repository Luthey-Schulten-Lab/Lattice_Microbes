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
#include "lptf/Profiler.h"
#include "lm/Cuda.h"
#include "TimingConstants.h"

// Allocate the profile space.
PROF_ALLOC;

__global__ void kernel(int* mem);

int main(int argc, char **argv)
{
    try
    {
        PROF_INIT;
        PROF_BEGIN(PROF_MAIN_RUN);

        dim3 grid(1);
        dim3 threads(32);
        cudaStream_t stream;
        int* deviceMem;
        CUDA_EXCEPTION_CHECK(cudaStreamCreate(&stream));
        CUDA_EXCEPTION_CHECK(cudaMalloc(&deviceMem, threads.x*sizeof(int)));

        // Launch a few kernels to warm up the system.
        for (int i=0; i<10; i++)
        {
            CUDA_EXCEPTION_EXECUTE((kernel<<<grid,threads,0,stream>>>(deviceMem)));
        }
        CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));

        // Start timings the kernels.
        PROF_BEGIN(PROF_SUBMIT_KERNELS);
        PROF_CUDA_START(stream);

        // Launch the kernels.
        int NUM_LAUNCHES=100;
        for (int i=0; i<NUM_LAUNCHES; i++)
        {
            PROF_CUDA_BEGIN(PROF_KERNEL_RUNNING,stream);
            CUDA_EXCEPTION_EXECUTE((kernel<<<grid,threads,0,stream>>>(deviceMem)));
            PROF_CUDA_END(PROF_KERNEL_RUNNING,stream);
        }

        // Wait for all of the kernels to finish.
        CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));

        // record the timings.
        PROF_CUDA_FINISH(stream);
        CUDA_EXCEPTION_CHECK(cudaFree(deviceMem));
        CUDA_EXCEPTION_CHECK(cudaStreamDestroy(stream));
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

__global__ void kernel(int* mem)
{
    mem[threadIdx.x] = threadIdx.x;
}
