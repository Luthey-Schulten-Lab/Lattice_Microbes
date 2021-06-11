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
#if defined(MACOSX)
#include <mach/mach_time.h>
#elif defined(LINUX)
#include <time.h>
#endif
#include "core/Exceptions.h"
#include "core/Types.h"
#include "rng/RandomGenerator.h"

namespace lm {
namespace rng {

RandomGenerator::RandomGenerator(uint32_t seedTop, uint32_t seedBottom, Distributions availableDists)
:seed((((uint64_t)(seedTop))<<32)|(uint64_t)(seedBottom)),availableDists(availableDists)
{
    if (seedBottom == 0)
    {
        #if defined(MACOSX)
        seed |= (uint32_t)mach_absolute_time();
        #elif defined(LINUX)
        struct timespec seed_timespec;
        if (clock_gettime(CLOCK_REALTIME, &seed_timespec) != 0) throw lm::Exception("Error getting time to use for random seed.");
        seed |= seed_timespec.tv_nsec;
        #endif
    }
}

void RandomGenerator::getRandomDoubles(double * rngs, int numberRNGs)
{
    for (int i=0; i<numberRNGs; i++)
    {
        rngs[i] = getRandomDouble();
    }
}

void RandomGenerator::getExpRandomDoubles(double * rngs, int numberRNGs)
{
    for (int i=0; i<numberRNGs; i++)
    {
        rngs[i] = getExpRandomDouble();
    }
}

void RandomGenerator::getNormRandomDoubles(double * rngs, int numberRNGs)
{
    for (int i=0; i<numberRNGs; i++)
    {
        rngs[i] = getNormRandomDouble();
    }
}

}
}
