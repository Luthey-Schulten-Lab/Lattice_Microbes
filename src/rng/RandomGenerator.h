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

#ifndef LM_RNG_RANDOMGENERATOR_H_
#define LM_RNG_RANDOMGENERATOR_H_

#include "core/Types.h"

namespace lm {
namespace rng {

/// @class RandomGenerator
/// @brief Base class for random number generators in Lattice Microbes
class RandomGenerator
{
public:
    /// @enum Distributions
    /// @brief Types of random number generators that are allowed
    enum Distributions
    {
       NONE        = 0x00,
       UNIFORM     = 0x01,
       EXPONENTIAL = 0x02,
       NORMAL      = 0x04,
       ALL         = 0xFF
    };

public:
    RandomGenerator(uint32_t seedTop, uint32_t seedBottom, Distributions availableDists=(Distributions)(ALL));
    virtual ~RandomGenerator() {}
    
    /// @brief Get the current seed
    /// @return seed
    virtual uint64_t getSeed() {return seed;}

    /// @brief Get a random integer
    /// @return random A random integer
    virtual uint32_t getRandom()=0;
    /// @brief Get a random double
    /// @return random A random double
    virtual double getRandomDouble()=0;
    /// @brief Get a random exponentially distributed double
    /// @return random A random double
    virtual double getExpRandomDouble()=0;
    /// @brief Get a random normally distributed double
    /// @return random A random double
    virtual double getNormRandomDouble()=0;
    
    /// @brief Get a number of random doubles
    /// @param rngs A preallocated array for random numbers
    /// @param numberRNGs Number or randoms to put in array
    virtual void getRandomDoubles(double * rngs, int numberRNGs);
    /// @brief Get a number of random exponentially distiributed doubles
    /// @param rngs A preallocated array for random numbers
    /// @param numberRNGs Number or randoms to put in array
    virtual void getExpRandomDoubles(double * rngs, int numberRNGs);
    /// @brief Get a number of random normally distributed doubles
    /// @param rngs A preallocated array for random numbers
    /// @param numberRNGs Number or randoms to put in array
    virtual void getNormRandomDoubles(double * rngs, int numberRNGs);

protected:
    uint64_t seed;                  // seed for the simulation
    Distributions availableDists;   // Which distributions are available to the simulation
};

}
}

#endif
