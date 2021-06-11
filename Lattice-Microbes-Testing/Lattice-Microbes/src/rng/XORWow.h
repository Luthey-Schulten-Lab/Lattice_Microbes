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

#ifndef LM_RNG_XORWOW_H_
#define LM_RNG_XORWOW_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "core/Types.h"
#include "rng/RandomGenerator.h"

namespace lm {
namespace rng {

class XORWow: public RandomGenerator
{
public:
    XORWow(int cudaDevice, uint32_t seedTop, uint32_t seedBottom, Distributions availableDists);
    virtual ~XORWow();

    virtual uint32_t getRandom();
    virtual double getRandomDouble();
    virtual double getExpRandomDouble();
    virtual double getNormRandomDouble();
    virtual void getRandomDoubles(double * rngs, int numberRNGs);
    virtual void getExpRandomDoubles(double * rngs, int numberRNGs);
    virtual void getNormRandomDoubles(double * rngs, int numberRNGs);

protected:
    virtual void generateRandomValues();
    virtual void launchGenerateKernel();
    int cudaDevice;
    curandState * state;
    cudaStream_t stream;
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    double * randomValues;
    double * nextRandomValues;
    double * randomValuesDev;
    double * expRandomValues;
    double * nextExpRandomValues;
    double * expRandomValuesDev;
    double * normRandomValues;
    double * nextNormRandomValues;
    double * normRandomValuesDev;
    #else
    float * randomValues;
    float * nextRandomValues;
    float * randomValuesDev;
    float * expRandomValues;
    float * nextExpRandomValues;
    float * expRandomValuesDev;
    float * normRandomValues;
    float * nextNormRandomValues;
    float * normRandomValuesDev;
    #endif
    const size_t numberValues;
    size_t nextValue;
};

__global__ void xorwow_init_kernel(const unsigned long long seed, curandState *rngState);
#ifdef RNG_CUDA_DOUBLE_PRECISION
__global__ void xorwow_generate_kernel(curandState *state, double * randomValues, double * expRandomValues, double * normRandomValues, uint iterations);
#else
__global__ void xorwow_generate_kernel(curandState *state, float * randomValues, float * expRandomValues, float * normRandomValues, uint iterations);
#endif

}
}

#endif
