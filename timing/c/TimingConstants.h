/*
 * University of Illinois Open Source License
 * Copyright 2010 Luthey-Schulten Group,
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
#ifndef TIMINGCONSTANTS_H_
#define TIMINGCONSTANTS_H_

#define PROF_MAIN_RUN                               1
#define PROF_SUBMIT_KERNELS                         2

// kernel_launch.cu
#define PROF_KERNEL_RUNNING                         3

// bit_packed_diffusion.cu
#define PROF_TIMESTEP_RUNNING                       4
#define PROF_X_DIFFUSION                            5
#define PROF_Y_DIFFUSION                            6
#define PROF_Z_DIFFUSION                            7
#define PROF_XYZ_DIFFUSION                          8
#define PROF_XY_DIFFUSION                           9


// rng_generate
#define PROF_XORSHIFT_INT                           300
#define PROF_XORSHIFT_FLOAT                         301
#define PROF_XORWOW_INIT                            302
#define PROF_XORWOW_INT                             303
#define PROF_XORWOW_FLOAT                           304




#endif /* TIMINGCONSTANTS_H_ */
