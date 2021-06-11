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

#ifndef LM_RNG_XORSHIFT_H_
#define LM_RNG_XORSHIFT_H_

#include "core/Types.h"
#include "rng/RandomGenerator.h"

namespace lm {
namespace rng {

/// @class XORShift
/// @brief A random number generator that works for a grid based on the index
class XORShift: public RandomGenerator
{
public:
    /// @brief Create the XORShift random number generator
    /// @param seedTop Top 32 bits of the random seed
    /// @param seedBottom Bottom 32 bits of the random seed
    XORShift(uint32_t seedTop, uint32_t seedBottom);
    /// @brief Destory the XORShift object
    virtual ~XORShift() {}

    /// @brief Get a random integer
    /// @return random A random integer
    virtual uint32_t getRandom();
    /// @brief Get a random double
    /// @return random A random double
    virtual double getRandomDouble();
    /// @brief Get a random exponentially distributed double
    /// @return random A random double
    virtual double getExpRandomDouble();
    /// @brief Get a random normally distributed double
    /// @return random A random double
    virtual double getNormRandomDouble();
    
    /// @brief Get a set of random doubles
    /// @param rngs Pointer to memory for random numbers
    /// @param numberRNGs The size of the rngs variable
    virtual void getRandomDoubles(double * rngs, int numberRNGs);
    /// @brief Get a set of random exponentially distributed doubles
    /// @param rngs Pointer to memory for random numbers
    /// @param numberRNGs The size of the rngs variable
    virtual void getExpRandomDoubles(double * rngs, int numberRNGs);
    /// @brief Get a set of random normally distributed doubles
    /// @param rngs Pointer to memory for random numbers
    /// @param numberRNGs The size of the rngs variable
    virtual void getNormRandomDoubles(double * rngs, int numberRNGs);


protected:
    unsigned long long state;
    bool isNextGaussianValid;
    double nextGaussian;
};

}
}

#endif
