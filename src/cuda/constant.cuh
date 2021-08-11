/*
 * University of Illinois Open Source License
 * Copyright 2018-2018 Luthey-Schulten Group,
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
 * Author(s): Tyler M. Earnest
 */

#ifndef _CONSTANT_CUH
#define _CONSTANT_CUH
#include "config.h"

#ifndef _CONSTANT_CU_
extern __constant__ unsigned int gpuidC;
extern __constant__ unsigned int numberReactionsC;
extern __constant__ unsigned int numberSpeciesC;
extern __constant__ unsigned int numberSiteTypesC;
extern __constant__ unsigned int latticeXSizeC;
extern __constant__ unsigned int latticeYSizeC;
extern __constant__ unsigned int latticeZSizeC;
extern __constant__ unsigned int latticeXYSizeC;
extern __constant__ unsigned int latticeXYZSizeC;
extern __constant__ unsigned int global_latticeZSizeC;
extern __constant__ unsigned int global_latticeXYZSizeC;
extern __constant__ int goffsetXC;
extern __constant__ int goffsetYC;
extern __constant__ int goffsetZC;
extern __constant__ unsigned int top_send_z;
extern __constant__ unsigned int top_recv_z;
extern __constant__ unsigned int top_dim_z;
extern __constant__ unsigned int top_dim_size;
extern __constant__ unsigned int bot_send_z;
extern __constant__ unsigned int bot_recv_z;
extern __constant__ unsigned int bot_dim_z;
extern __constant__ unsigned int bot_dim_size;

#ifndef MPD_GLOBAL_R_MATRIX
extern __constant__ unsigned int reactionOrdersC[];
extern __constant__ unsigned int reactionSitesC[];
extern __constant__ unsigned int D1C[];
extern __constant__ unsigned int D2C[];
extern __constant__ float reactionRatesC[];
#endif

#ifndef MPD_GLOBAL_S_MATRIX
extern __constant__ int8_t SC[];
extern __constant__ uint8_t RLC[];
#endif

#ifndef MPD_GLOBAL_T_MATRIX
extern __constant__ float TC[];
#else
extern __constant__ float* TC;
#endif

#endif /* _CONSTANT_CU */

#endif /* _CONSTANT_CUH */
