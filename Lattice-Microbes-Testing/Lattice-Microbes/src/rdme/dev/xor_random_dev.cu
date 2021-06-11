/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group, 
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

#include <cuda.h>
#include <cuda_runtime.h>
#include "config.h"

#define UNIFORM_RAND_STEP_SIZE	2.32830643654e-10f		// 1/(2^32)

inline __device__ unsigned int getRandomHash(unsigned long long v, const unsigned long long seed)
{
	v ^= seed;
	v = v * 3935559000370003845ULL + 2691343689449507681ULL;
	v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
	v *= 2685821657736338717ULL;
	return v>>32;
}

inline __device__ unsigned int getRandomHash(const unsigned int v1, const unsigned int v1Shift, const unsigned int v2, const unsigned long long seed)
{
	unsigned long long v = v1;
	v <<= v1Shift;
	v |= v2;
	return getRandomHash(v, seed);
}

/**
 * Returns a float value in the range [0.0 1.0).
 */
inline __device__ float getRandomHashFloat(unsigned long long v, const unsigned long long seed)
{
	unsigned int r = getRandomHash(v, seed);
	return ((float)r)*UNIFORM_RAND_STEP_SIZE;
}

/**
 * Returns a float value in the range [0.0 1.0).
 */
inline __device__ float getRandomHashFloat(const unsigned int v1, const unsigned int v1Shift, const unsigned int v2, const unsigned long long seed)
{
	unsigned long long v = v1;
	v <<= v1Shift;
	v |= v2;
	return getRandomHashFloat(v, seed);
}

