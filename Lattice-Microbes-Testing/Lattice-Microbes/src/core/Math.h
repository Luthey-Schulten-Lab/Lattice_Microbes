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

#ifndef LM_MATH_H_
#define LM_MATH_H_

#include <algorithm>
#include <cmath>

#define TWOPI       6.283185307179586476925287
#define PI          3.141592653589793238462643
#define PID2        1.570796326794896619231322
#define PID4        0.7853981634
#define PIOVER180   0.0174532925
#define NA          6.02214179e23
#define EPS         1e-12



inline bool isPower2(unsigned int x)
{
    return !(x&(x - 1)) && x;
}

inline bool isPower2(unsigned long x)
{
    return !(x&(x - 1)) && x;
}

inline bool isPower2(unsigned long long x)
{
    return !(x&(x - 1)) && x;
}

inline unsigned int log2(unsigned int x)
{
    unsigned int r = 0;
    while (x>>=1)
        r++;
    return r;
}

inline unsigned int log2(unsigned long x)
{
    unsigned int r = 0;
    while (x>>=1)
        r++;
    return r;
}

inline unsigned int log2(unsigned long long x)
{
    unsigned int r = 0;
    while (x>>=1)
        r++;
    return r;
}


#ifndef __cuda_cuda_h__
using std::min;
using std::max;
/*
template <class T> inline T min(T x, T y) {return x<y?x:y;}
template <class T> inline T min(T x, T y, T z) {return x<y?(x<z?x:z):(y<z?y:z);}
template <class T> inline T max(T x, T y) {return x>y?x:y;}
template <class T> inline T max(Tt x, T y, T z) {return x>y?(x>z?x:z):(y>z?y:z);}
*/
#endif

#endif
