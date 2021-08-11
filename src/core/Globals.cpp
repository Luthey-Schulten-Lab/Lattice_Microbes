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

#include <string>
#include <list>
#include <vector>
#include <ctime>
#include "config.h"
#include "core/Types.h"
#include "core/Exceptions.h"
#include "me/MESolverFactory.h"


using std::string;
using std::vector;

/**
 * The function being performed.
 */
string functionOption = "interpreter";

/**
 * The name of the file containing the simulation.
 */
string simulationFilename;

/**
 * The number of replicates of the simulation that should be performed.
 */
vector<int> replicates;

/**
 * The interval at which the results file should be checkpointed.
 */
time_t checkpointInterval = 0;

/**
 * If a global abort signal has been received.
 */
volatile bool globalAbort = false;

/**
 * The solver to use for the simulations.
 */
lm::me::MESolverFactory solverFactory;

/**
 * The number of cpu cores assigned to each process.
 */
int numberCpuCores;

/**
 * The number of cpu cores to assign per replicate (can be a fraction, e.g., 1/2, 1/4, etc).
 */
float cpuCoresPerReplicate;

#ifdef OPT_CUDA

/**
 * The cuda devices assigned to each process.
 */
vector<int> cudaDevices;

/**
 * The number of cuda devices to assign per replicate (can be a fraction, e.g., 1/2, 1/4, etc).
 */
float cudaDevicesPerReplicate;

/**
 * Whether we should print the cuda device capabilities on startup.
 */
bool shouldPrintCudaCapabilities;

/** 
 * Boolean to signify if peering should be enabled for multi-GPU communication
 */
bool mgpu_disablePeering=false;

#endif

/**
 * Whether we should reserve a core for the output thread.
 */
bool shouldReserveOutputCore;


int cudaDevicesPerNode;

// Parameters given on command line
vector<string> cmdline_parameters;

