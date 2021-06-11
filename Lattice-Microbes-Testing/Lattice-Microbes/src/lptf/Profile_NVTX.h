/*
 * University of Illinois Open Source License
 * Copyright 2014-2018 Luthey-Schulten Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.illinois.edu/schulten/lm
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
 * Author(s): Mike Hallock
 */
#ifndef PROF_NVTX_H
#define PROF_NVTX_H

enum PROFILE_CODE {
PROF_THREAD_VARIABLE_START=0 ,

PROF_MAIN_RUN                               ,
PROF_SIM_RUN                                ,

PROF_MASTER_READ_STATIC_MSG                 ,
PROF_MASTER_READ_FINISHED_MSG               ,
PROF_MASTER_FINISHED_THREAD                 ,
PROF_MASTER_SLEEP                           ,

PROF_REPLICATE_EXECUTE                      ,
PROF_REPLICATE_WRITE_DATASET                ,

PROF_DATAOUTPUT_RUN                         ,
PROF_DATAOUTPUT_WRITE_DATASET               ,
PROF_DATAOUTPUT_HDF_WRITE_COUNTS            ,
PROF_DATAOUTPUT_HDF_WRITE_FPT               ,
PROF_DATAOUTPUT_HDF_WRITE_PV                ,
PROF_DATAOUTPUT_HDF_WRITE_LATTICE           ,

PROF_SLAVE_SLEEP                            ,

PROF_SIM_EXECUTE                            ,
PROF_SERIALIZE_COUNTS                       ,
PROF_SERIALIZE_FPT                          ,
PROF_DESERIALIZE_COUNTS                     ,
PROF_DESERIALIZE_FPT                        ,
PROF_SERIALIZE_PV                           ,
PROF_DESERIALIZE_PV                         ,
PROF_SERIALIZE_LATTICE                      ,
PROF_DESERIALIZE_LATTICE                    ,
PROF_DETERMINE_COUNTS                       ,

PROF_INIT_XORWOW_RNG                        ,
PROF_GENERATE_XORWOW_RNG                    ,
PROF_CACHE_RNG                              ,
PROF_LAUNCH_XORWOW_RNG                      ,
PROF_COPY_XORWOW_RNG                        ,
PROF_COPY_XORWOW_EXP_RNG                    ,
PROF_COPY_XORWOW_NORM_RNG                   ,
PROF_CACHE_EXP_RNG                          ,

PROF_MPD_TIMESTEP                           ,
PROF_MPD_X_DIFFUSION                        ,
PROF_MPD_Y_DIFFUSION                        ,
PROF_MPD_Z_DIFFUSION                        ,
PROF_MPD_REACTION                           ,
PROF_MPD_SYNCHRONIZE                        ,
PROF_MPD_OVERFLOW                           ,

PROF_NSM_INIT_QUEUE                         ,
PROF_NSM_BUILD_QUEUE                        ,
PROF_NSM_LOOP		                        ,

PROF_H2D,
PROF_D2H,
PROF_D2D,
PROF_P2P,

MPD_MCLKR_BARRIER
};

static const char *profile_description[] = {
"THREAD_VARIABLE_START",

"MAIN_RUN",
"SIM_RUN",

"MASTER_READ_STATIC_MSG",
"MASTER_READ_FINISHED_MSG",
"MASTER_FINISHED_THREAD",
"MASTER_SLEEP",

"REPLICATE_EXECUTE",
"REPLICATE_WRITE_DATASET",

"DATAOUTPUT_RUN",
"DATAOUTPUT_WRITE_DATASET",
"DATAOUTPUT_HDF_WRITE_COUNT",
"DATAOUTPUT_HDF_WRITE_FPT",
"DATAOUTPUT_HDF_WRITE_PV",
"DATAOUTPUT_HDF_WRITE_LATTICE",

"SLAVE_SLEEP",

"SIM_EXECUTE",
"SERIALIZE_COUNTS",
"SERIALIZE_FPT",
"DESERIALIZE_COUNTS",
"DESERIALIZE_FPT",
"SERIALIZE_PV",
"DESERIALIZE_PV",
"SERIALIZE_LATTICE",
"DESERIALIZE_LATTICE",
"DETERMINE_COUNTS",

"INIT_XORWOW_RN",
"GENERATE_XORWOW_RNG",
"CACHE_RNG",
"LAUNCH_XORWOW_RNG",
"COPY_XORWOW_RNG",
"COPY_XORWOW_EXP_RNG",
"COPY_XORWOW_NORM_RNG",
"CACHE_EXP_RNG",

"MPD_TIMESTEP",
"MPD_X_DIFFUSION",
"MPD_Y_DIFFUSION ",
"MPD_Z_DIFFUSION",
"MPD_REACTION",
"MPD_SYNCHRONIZE",
"MPD_OVERFLOW",

"NSM_INIT_QUEUE",
"NSM_BUILD_QUEUE",
"NSM_LOOP",

"Memcopy H->D",
"Memcopy D->H",
"Memcopy D->D",
"Memcopy P->P",

"MCLKR Barrier"

};



#endif /* PROF_NVTX_H */
