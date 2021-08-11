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

#include <cmath>
#include "config.h"
#include "core/Types.h"
#include "rng/RandomGenerator.h"
#include "rng/XORShift.h"

namespace lm {
namespace rng {

XORShift::XORShift(uint32_t seedTop, uint32_t seedBottom)
:RandomGenerator(seedTop,seedBottom),state(1),isNextGaussianValid(false)
{
}

uint32_t XORShift::getRandom()
{
    unsigned long long v = state++;
    v ^= seed;
    v = v * 3935559000370003845ULL + 2691343689449507681ULL;
    v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
    v *= 2685821657736338717ULL;
    return v>>32;
}

/**
 * Returns a float value in the range [0.0 1.0).
 */
double XORShift::getRandomDouble()
{
    uint32_t r = getRandom();
    return ((double)r)*(2.328306436539e-10); //1/(2^32)
}

/**
 * Returns an exponentially distributed value.
 */
double XORShift::getExpRandomDouble()
{
    uint64_t r;
    while ((r=(uint64_t)getRandom()) == 0);
    double d = ((double)r)*(2.328306435997e-10); //1/((2^32)+1) range (0.0 1.0)
    return -log(d);
}

/**
 * Returns an normally distributed value.
 */
double XORShift::getNormRandomDouble()
{
    if (isNextGaussianValid)
    {
        isNextGaussianValid = false;
        return nextGaussian;
    }

    // Generate two uniform random numbers.
    uint64_t r;
    while ((r=(uint64_t)getRandom()) == 0);
    double d1 = ((double)r)*(2.328306435997e-10); //1/((2^32)+1) range (0.0 1.0)
    while ((r=(uint64_t)getRandom()) == 0);
    double d2 = ((double)r)*(2.328306435997e-10 * 6.2831853071795860); //1/((2^32)+1)*PI range (0.0 2PI)

    // Transform using Box-Muller.
    isNextGaussianValid = true;
    double s = sqrt(-2.0 * log(d1));
    #if defined(LINUX)
    double ret;
    sincos(d2, &nextGaussian, &ret);
    nextGaussian *= s;
    return s*ret;
    #else
    nextGaussian = s * sin(d2);
    return s * cos(d2);
    #endif
}

void XORShift::getRandomDoubles(double * rngs, int numberRNGs)
{
    for (int i=0; i<numberRNGs; i++)
    {
        uint32_t r = getRandom();
        rngs[i] = ((double)r)*(2.328306436539e-10); //1/(2^32)
    }
}

void XORShift::getExpRandomDoubles(double * rngs, int numberRNGs)
{
    for (int i=0; i<numberRNGs; i++)
    {
        uint64_t r;
        while ((r=(uint64_t)getRandom()) == 0);
        double d = ((double)r)*(2.328306435997e-10); //1/((2^32)+1) range (0.0 1.0)
        rngs[i] = -log(d);
    }
}

void XORShift::getNormRandomDoubles(double * rngs, int numberRNGs)
{
    // Generate an even number of rngs that does not exceed the buffer size.
    int i;
    for (i=0; i<(numberRNGs>>1)<<1; i+=2)
    {
        // Generate two uniform random numbers.
        uint64_t r;
        while ((r=(uint64_t)getRandom()) == 0);
        double d1 = ((double)r)*(2.328306435997e-10); //1/((2^32)+1) range (0.0 1.0)
        while ((r=(uint64_t)getRandom()) == 0);
        double d2 = ((double)r)*(2.328306435997e-10 * 6.2831853071795860); //1/((2^32)+1)*PI range (0.0 2PI)

        // Transform using Box-Muller.
        double s = sqrt(-2.0 * log(d1));
        #if defined(LINUX)
        sincos(d2, &rngs[i], &rngs[i+1]);
        rngs[i] *= s;
        rngs[i+1] *= s;
        #else
        rngs[i] = s * sin(d2);
        rngs[i+1] = s * cos(d2);
        #endif
    }

    // If there is one more number to generate, do so and discard the second.
    if (i<numberRNGs)
    {
        // Generate two uniform random numbers.
        uint64_t r;
        while ((r=(uint64_t)getRandom()) == 0);
        double d1 = ((double)r)*(2.328306435997e-10); //1/((2^32)+1) range (0.0 1.0)
        while ((r=(uint64_t)getRandom()) == 0);
        double d2 = ((double)r)*(2.328306435997e-10 * 6.2831853071795860); //1/((2^32)+1)*PI range (0.0 2PI)

        // Transform using Box-Muller.
        isNextGaussianValid = true;
        double s = sqrt(-2.0 * log(d1));
        #if defined(LINUX)
        sincos(d2, &rngs[i], &rngs[i+1]);
        rngs[i] *= s;
        rngs[i+1] *= s;
        #else
        rngs[i] = s * sin(d2);
        rngs[i+1] = s * cos(d2);
        #endif
    }
}

}
}
