    /*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
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

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "config.h"
#include "cuda/lm_cuda.h"
#include "core/Types.h"
#include "rng/RandomGenerator.h"
#include "rng/XORWow.h"
#include "lptf/Profile.h"

namespace lm {
namespace rng {

XORWow::XORWow(int cudaDevice, uint32_t seedTop, uint32_t seedBottom, Distributions availableDists)
:RandomGenerator(seedTop,seedBottom,availableDists),cudaDevice(cudaDevice),state(NULL),stream(NULL),randomValues(NULL),nextRandomValues(NULL),randomValuesDev(NULL),expRandomValues(NULL),nextExpRandomValues(NULL),expRandomValuesDev(NULL),normRandomValues(NULL),nextNormRandomValues(NULL),normRandomValuesDev(NULL),numberValues(RNG_TUNE_XORWOW_BLOCK_SIZE*RNG_TUNE_XORWOW_GRID_SIZE*RNG_TUNE_XORWOW_THREAD_ITERATIONS),nextValue(0)
{
    // Create the cuda stream.
    CUDA_EXCEPTION_CHECK(cudaStreamCreate(&stream));

    // Allocate the memory for the rng state.
    CUDA_EXCEPTION_CHECK(cudaMalloc((void **)&state, RNG_TUNE_XORWOW_BLOCK_SIZE*RNG_TUNE_XORWOW_GRID_SIZE*sizeof(curandState)));

    // Allocate memory for the random values.
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&randomValues, numberValues*sizeof(double), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&nextRandomValues, numberValues*sizeof(double), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaMalloc((void **)&randomValuesDev, numberValues*sizeof(double)));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&expRandomValues, numberValues*sizeof(double), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&nextExpRandomValues, numberValues*sizeof(double), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaMalloc((void **)&expRandomValuesDev, numberValues*sizeof(double)));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&normRandomValues, numberValues*sizeof(double), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&nextNormRandomValues, numberValues*sizeof(double), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaMalloc((void **)&normRandomValuesDev, numberValues*sizeof(double)));
    #else
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&randomValues, numberValues*sizeof(float), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&nextRandomValues, numberValues*sizeof(float), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaMalloc((void **)&randomValuesDev, numberValues*sizeof(float)));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&expRandomValues, numberValues*sizeof(float), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&nextExpRandomValues, numberValues*sizeof(float), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaMalloc((void **)&expRandomValuesDev, numberValues*sizeof(float)));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&normRandomValues, numberValues*sizeof(float), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaHostAlloc((void **)&nextNormRandomValues, numberValues*sizeof(float), cudaHostAllocDefault));
    CUDA_EXCEPTION_CHECK(cudaMalloc((void **)&normRandomValuesDev, numberValues*sizeof(float)));
    #endif

    PROF_CUDA_START(stream);

    // Generate the initial state for the rng.
    PROF_CUDA_BEGIN(PROF_INIT_XORWOW_RNG,stream);
    CUDA_EXCEPTION_EXECUTE((xorwow_init_kernel<<<RNG_TUNE_XORWOW_BLOCK_SIZE,RNG_TUNE_XORWOW_GRID_SIZE,0,stream>>>(seed, state)));
    PROF_CUDA_END(PROF_INIT_XORWOW_RNG,stream);

    // Start the first round of random values.
    launchGenerateKernel();
    nextValue = numberValues;

    // Wait for the first rng generation to finish.
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));
}

XORWow::~XORWow()
{
    // Close the streams.
    if (stream != NULL)
    {
        // Wait for any kernels or mem copies to finish so they aren't using invalid memory.
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaStreamSynchronize(stream));
        PROF_CUDA_FINISH(stream);

        CUDA_EXCEPTION_CHECK_NOTHROW(cudaStreamDestroy(stream));
        stream = NULL;
    }

    // Free any memory.
    if (randomValues != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(randomValues));
        randomValues = NULL;
    }
    if (nextRandomValues != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(nextRandomValues));
        nextRandomValues = NULL;
    }
    if (randomValuesDev != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFree(randomValuesDev));
        randomValuesDev = NULL;
    }
    if (expRandomValues != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(expRandomValues));
        expRandomValues = NULL;
    }
    if (nextExpRandomValues != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(nextExpRandomValues));
        nextExpRandomValues = NULL;
    }
    if (expRandomValuesDev != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFree(expRandomValuesDev));
        expRandomValuesDev = NULL;
    }
    if (normRandomValues != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(normRandomValues));
        normRandomValues = NULL;
    }
    if (nextNormRandomValues != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(nextNormRandomValues));
        nextNormRandomValues = NULL;
    }
    if (normRandomValuesDev != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFree(normRandomValuesDev));
        normRandomValuesDev = NULL;
    }
    if (state != NULL)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFree(state));
        state = NULL;
    }
}

void XORWow::launchGenerateKernel()
{
    PROF_CUDA_START(stream);

    // Generate the rngs.
    PROF_CUDA_BEGIN(PROF_GENERATE_XORWOW_RNG,stream);
    CUDA_EXCEPTION_EXECUTE((xorwow_generate_kernel<<<RNG_TUNE_XORWOW_BLOCK_SIZE,RNG_TUNE_XORWOW_GRID_SIZE,0,stream>>>(state, randomValuesDev, expRandomValuesDev, normRandomValuesDev, RNG_TUNE_XORWOW_THREAD_ITERATIONS)));
    PROF_CUDA_END(PROF_GENERATE_XORWOW_RNG,stream);

    // Copy the rng values back to the host.
    PROF_CUDA_BEGIN(PROF_COPY_XORWOW_RNG,stream);
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(nextRandomValues, randomValuesDev, numberValues*sizeof(double), cudaMemcpyDeviceToHost, stream));
    #else
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(nextRandomValues, randomValuesDev, numberValues*sizeof(float), cudaMemcpyDeviceToHost, stream));
    #endif
    PROF_CUDA_END(PROF_COPY_XORWOW_RNG,stream);

    // Copy the exp rng values back to the host.
    PROF_CUDA_BEGIN(PROF_COPY_XORWOW_EXP_RNG,stream);
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(nextExpRandomValues, expRandomValuesDev, numberValues*sizeof(double), cudaMemcpyDeviceToHost, stream));
    #else
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(nextExpRandomValues, expRandomValuesDev, numberValues*sizeof(float), cudaMemcpyDeviceToHost, stream));
    #endif
    PROF_CUDA_END(PROF_COPY_XORWOW_EXP_RNG,stream);

    // Copy the norm rng values back to the host.
    PROF_CUDA_BEGIN(PROF_COPY_XORWOW_NORM_RNG,stream);
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(nextNormRandomValues, normRandomValuesDev, numberValues*sizeof(double), cudaMemcpyDeviceToHost, stream));
    #else
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(nextNormRandomValues, normRandomValuesDev, numberValues*sizeof(float), cudaMemcpyDeviceToHost, stream));
    #endif
    PROF_CUDA_END(PROF_COPY_XORWOW_NORM_RNG,stream);
}

void XORWow::generateRandomValues()
{
    PROF_BEGIN(PROF_LAUNCH_XORWOW_RNG);

    // Wait for the previous rng generation launch to finish, if it hasn't yet.
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream));
    PROF_CUDA_FINISH(stream);

    // Swap the random value pointers.
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    double * tmp;
    #else
    float * tmp;
    #endif
    tmp = randomValues;
    randomValues = nextRandomValues;
    nextRandomValues = tmp;

    // Swap the exp random value pointers.
    tmp = expRandomValues;
    expRandomValues = nextExpRandomValues;
    nextExpRandomValues = tmp;

    // Swap the norm random value pointers.
    tmp = normRandomValues;
    normRandomValues = nextNormRandomValues;
    nextNormRandomValues = tmp;

    // Mark that we are at the beginning again.
    nextValue = 0;

    // Start the next round of random value generation, but don't wait for it.
    launchGenerateKernel();
    PROF_END(PROF_LAUNCH_XORWOW_RNG);
}

uint32_t XORWow::getRandom()
{
    return 1000000;
}

double XORWow::getRandomDouble()
{
    if (nextValue >= numberValues)
    {
        generateRandomValues();
    }
    return randomValues[nextValue++];
}

void XORWow::getRandomDoubles(double * rngs, int numberRNGs)
{
    // If we don't have enough values, copy what is left and then generate more.
    if (numberRNGs > numberValues-nextValue)
    {
        int availableRNGs = numberValues-nextValue;

        #ifdef RNG_CUDA_DOUBLE_PRECISION
        memcpy(rngs, &randomValues[nextValue], sizeof(double)*availableRNGs);
        #else
        for (int i=0; i<availableRNGs; i++) rngs[i] = randomValues[nextValue+i];
        #endif

        rngs = &rngs[availableRNGs];
        numberRNGs -= availableRNGs;
        generateRandomValues();
    }

    // Copy the values into the buffer.
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    memcpy(rngs, &randomValues[nextValue], sizeof(double)*numberRNGs);
    #else
    for (int i=0; i<numberRNGs; i++) rngs[i] = randomValues[nextValue+i];
    #endif
    nextValue += numberRNGs;
}

double XORWow::getExpRandomDouble()
{
    if (nextValue >= numberValues)
    {
        generateRandomValues();
    }
    return expRandomValues[nextValue++];
}

void XORWow::getExpRandomDoubles(double * rngs, int numberRNGs)
{
    // If we don't have enough values, copy what is left and then generate more.
    if (numberRNGs > numberValues-nextValue)
    {
        int availableRNGs = numberValues-nextValue;

        #ifdef RNG_CUDA_DOUBLE_PRECISION
        memcpy(rngs, &expRandomValues[nextValue], sizeof(double)*availableRNGs);
        #else
        for (int i=0; i<availableRNGs; i++) rngs[i] = expRandomValues[nextValue+i];
        #endif

        rngs = &rngs[availableRNGs];
        numberRNGs -= availableRNGs;
        generateRandomValues();
    }

    // Copy the values into the buffer.
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    memcpy(rngs, &expRandomValues[nextValue], sizeof(double)*numberRNGs);
    #else
    for (int i=0; i<numberRNGs; i++) rngs[i] = expRandomValues[nextValue+i];
    #endif
    nextValue += numberRNGs;
}

double XORWow::getNormRandomDouble()
{
    if (nextValue >= numberValues)
    {
        generateRandomValues();
    }
    return normRandomValues[nextValue++];
}

void XORWow::getNormRandomDoubles(double * rngs, int numberRNGs)
{
    // If we don't have enough values, copy what is left and then generate more.
    if (numberRNGs > numberValues-nextValue)
    {
        int availableRNGs = numberValues-nextValue;

        #ifdef RNG_CUDA_DOUBLE_PRECISION
        memcpy(rngs, &normRandomValues[nextValue], sizeof(double)*availableRNGs);
        #else
        for (int i=0; i<availableRNGs; i++) rngs[i] = normRandomValues[nextValue+i];
        #endif

        rngs = &rngs[availableRNGs];
        numberRNGs -= availableRNGs;
        generateRandomValues();
    }

    // Copy the values into the buffer.
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    memcpy(rngs, &normRandomValues[nextValue], sizeof(double)*numberRNGs);
    #else
    for (int i=0; i<numberRNGs; i++) rngs[i] = normRandomValues[nextValue+i];
    #endif
    nextValue += numberRNGs;
}

__global__ void xorwow_init_kernel(const unsigned long long seed, curandState *state)
{
    int id = (blockIdx.x*blockDim.x)+threadIdx.x;
    curand_init(seed, id, 0, &state[id]);
}

#ifdef RNG_CUDA_DOUBLE_PRECISION
__global__ void xorwow_generate_kernel(curandState *state, double * randomValues, double * expRandomValues, double * normRandomValues, uint iterations)
#else
__global__ void xorwow_generate_kernel(curandState *state, float * randomValues, float * expRandomValues, float * normRandomValues, uint iterations)
#endif
{
    int id = (blockIdx.x*blockDim.x)+threadIdx.x;
    int valueIndex = id;
    int valueIndex2 = id+blockDim.x*gridDim.x;
    int offset = 2*blockDim.x*gridDim.x;

    // Get a local copy of the state.
    curandState localState = state[id];

    for (uint i=0; i<iterations; i+=2, valueIndex+=offset, valueIndex2+=offset)
    {
        #ifdef RNG_CUDA_DOUBLE_PRECISION
        unsigned int x1, y1, x2, y2;
        x1 = curand(&localState);
        y1 = curand(&localState);
        x2 = curand(&localState);
        y2 = curand(&localState);
        double u1 = _curand_uniform_double_hq(x1, y1);
        double u2 = _curand_uniform_double_hq(x2, y2);
        randomValues[valueIndex] = u1;
        randomValues[valueIndex2] = u2;
        expRandomValues[valueIndex] = (u1>=1.0)?(CURAND_2POW32_INV_DOUBLE):(-log(u1));
        expRandomValues[valueIndex2] = (u2>=1.0)?(CURAND_2POW32_INV_DOUBLE):(-log(u2));
        double2 n = _curand_box_muller_double(x1, y1, x2, y2);
        normRandomValues[valueIndex] = n.x;
        normRandomValues[valueIndex2] = n.y;
        #else
        unsigned int x1, x2;
        x1 = curand(&localState);
        x2 = curand(&localState);
        float u1 = _curand_uniform(x1);
        float u2 = _curand_uniform(x2);
        randomValues[valueIndex] = u1;
        randomValues[valueIndex2] = u2;
        expRandomValues[valueIndex] = (u1>=1.0f)?(CURAND_2POW32_INV):(-logf(u1));
        expRandomValues[valueIndex2] = (u2>=1.0f)?(CURAND_2POW32_INV):(-logf(u2));
        float2 n = _curand_box_muller(x1, x2);
        normRandomValues[valueIndex] = n.x;
        normRandomValues[valueIndex2] = n.y;
        #endif
    }

    // Copy the state back to global memory.
    state[id] = localState;
}

}
}
