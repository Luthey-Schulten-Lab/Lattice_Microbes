/*
 * University of Illinois Open Source License
 * Copyright 2010-2018 Luthey-Schulten Group,
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

#define PROF_THREAD_VARIABLE_START                  3

#define PROF_MAIN_RUN                               1
#define PROF_SIM_RUN                                2

#define PROF_MASTER_READ_STATIC_MSG                 3
#define PROF_MASTER_READ_FINISHED_MSG               4
#define PROF_MASTER_FINISHED_THREAD                 5
#define PROF_MASTER_SLEEP                           6

#define PROF_REPLICATE_EXECUTE                      10
#define PROF_REPLICATE_WRITE_DATASET                11

#define PROF_DATAOUTPUT_RUN                         100
#define PROF_DATAOUTPUT_WRITE_DATASET               101
#define PROF_DATAOUTPUT_HDF_WRITE_COUNTS            110
#define PROF_DATAOUTPUT_HDF_WRITE_FPT               111
#define PROF_DATAOUTPUT_HDF_WRITE_PV                112
#define PROF_DATAOUTPUT_HDF_WRITE_LATTICE           113

#define PROF_SLAVE_SLEEP                            200

#define PROF_SIM_EXECUTE                            299
#define PROF_SERIALIZE_COUNTS                       300
#define PROF_SERIALIZE_FPT                          301
#define PROF_DESERIALIZE_COUNTS                     302
#define PROF_DESERIALIZE_FPT                        303
#define PROF_SERIALIZE_PV                           304
#define PROF_DESERIALIZE_PV                         305
#define PROF_SERIALIZE_LATTICE                      306
#define PROF_DESERIALIZE_LATTICE                    307
#define PROF_DETERMINE_COUNTS                       308

#define PROF_INIT_XORWOW_RNG                        324
#define PROF_GENERATE_XORWOW_RNG                    325
#define PROF_CACHE_RNG                              326
#define PROF_LAUNCH_XORWOW_RNG                      327
#define PROF_COPY_XORWOW_RNG                        328
#define PROF_COPY_XORWOW_EXP_RNG                    329
#define PROF_COPY_XORWOW_NORM_RNG                   330
#define PROF_CACHE_EXP_RNG                          331

#define PROF_MPD_TIMESTEP                           500
#define PROF_MPD_X_DIFFUSION                        501
#define PROF_MPD_Y_DIFFUSION                        502
#define PROF_MPD_Z_DIFFUSION                        503
#define PROF_MPD_REACTION                           504
#define PROF_MPD_SYNCHRONIZE                        505
#define PROF_MPD_OVERFLOW                           506

#define PROF_BARRIER								550

#define PROF_MCLKR_TS_LOOP							551
#define PROF_MCLKR_PUBLISH							552
#define PROF_MCLKR_OVERFLOWS						553
#define PROF_MCLKR_DDX								554

#define PROF_D2H									555
#define PROF_H2D									556
#define PROF_D2D									557
#define PROF_CONST									558
#define PROF_FETCH_EX								569

#define PROF_XD1									570
#define PROF_XD2									571
#define PROF_XD3									572
#define PROF_RX1									573
#define PROF_RX2									574
#define PROF_RX3									575

#define PROF_MRBASE									576
#define PROF_MR1									576
#define PROF_MR2									577
#define PROF_MSBASE									578
#define PROF_MS1									578
#define PROF_MS2									579


#define PROF_NSM_INIT_QUEUE                         600
#define PROF_NSM_BUILD_QUEUE                        601
#define PROF_NSM_LOOP		                        602


#endif /* TIMINGCONSTANTS_H_ */
