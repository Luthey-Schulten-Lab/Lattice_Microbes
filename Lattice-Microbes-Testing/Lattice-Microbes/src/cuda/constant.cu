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

#define _CONSTANT_CU_
#include "config.h"
#include <cstdint>
#include "cuda/constant.cuh"
__constant__ unsigned int gpuidC;
__constant__ unsigned int numberReactionsC;
__constant__ unsigned int numberSpeciesC;
__constant__ unsigned int numberSiteTypesC;
__constant__ unsigned int latticeXSizeC;
__constant__ unsigned int latticeYSizeC;
__constant__ unsigned int latticeZSizeC;
__constant__ unsigned int latticeXYSizeC;
__constant__ unsigned int latticeXYZSizeC;
__constant__ unsigned int global_latticeZSizeC;
__constant__ unsigned int global_latticeXYZSizeC;


// goffset - the offset of the local lattice as it appears in the global space
// This can be negative for periofic boundaries.  It will be interpreted as 
// offsets from the other end
__constant__ int goffsetXC;
__constant__ int goffsetYC;
__constant__ int goffsetZC;

__constant__ unsigned int top_send_z;
__constant__ unsigned int top_recv_z;
__constant__ unsigned int top_dim_z;
__constant__ unsigned int top_dim_size;
__constant__ unsigned int bot_send_z;
__constant__ unsigned int bot_recv_z;
__constant__ unsigned int bot_dim_z;
__constant__ unsigned int bot_dim_size;

#ifndef MPD_GLOBAL_R_MATRIX
__constant__ unsigned int reactionOrdersC[MPD_MAX_REACTION_TABLE_ENTRIES];
__constant__ unsigned int reactionSitesC[MPD_MAX_REACTION_TABLE_ENTRIES];
__constant__ unsigned int D1C[MPD_MAX_REACTION_TABLE_ENTRIES];
__constant__ unsigned int D2C[MPD_MAX_REACTION_TABLE_ENTRIES];
__constant__ float reactionRatesC[MPD_MAX_REACTION_TABLE_ENTRIES];
#endif

#ifndef MPD_GLOBAL_S_MATRIX
__constant__ int8_t SC[MPD_MAX_S_MATRIX_ENTRIES]; // The stoichiometric matrix: numberSpecies x numberReactions
__constant__ uint8_t RLC[MPD_MAX_RL_MATRIX_ENTRIES]; // The reaction location matrix: numberReaction x numberSiteTypes
#endif

#ifndef MPD_GLOBAL_T_MATRIX
__constant__ float TC[MPD_MAX_TRANSITION_TABLE_ENTRIES];
#else
__constant__ float* TC;
#endif
