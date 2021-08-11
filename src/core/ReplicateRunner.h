/*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
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

#ifndef LM_MAIN_REPLICATERUNNER_H_
#define LM_MAIN_REPLICATERUNNER_H_

#include <map>
#include <string>
#include "core/ResourceAllocator.h"
#include "me/MESolverFactory.h"
#include "ReactionModel.pb.h"
#include "DiffusionModel.pb.h"
#include "thread/Thread.h"
#include "thread/Worker.h"

using std::map;
using std::string;
using lm::thread::PthreadException;
using lm::thread::Worker;
using lm::me::MESolverFactory;

namespace lm {
namespace main {

/// @class ReplicateRunner
/// @brief A thread that launches all the various replicates requested
class ReplicateRunner : public Worker
{
public:
    /// @brief Create a new replicate runner
    /// @param replicate The number of the first replicate
    /// @param solverFactory The factory used to create a solver object for each replicate
    /// @param parameters A map of parameters for the simulation
    /// @param reactionModel An object to take care of reactions
    /// @param diffusionModel An object to take care of diffusion
    /// @param lattice The actual data representing the lattice
    /// @param latticeSize The size of the lattice in bytes
    /// @param latticeSitesSize The size of the lattice sites in bytes
    /// @param resources A manager for the resources that handles GPUs and CPUs
    ReplicateRunner(int replicate, MESolverFactory solverFactory, map<string,string> * parameters, lm::io::ReactionModel * reactionModel, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize, ResourceAllocator::ComputeResources resources);
    virtual ~ReplicateRunner();
    
    /// @brief Wake the thread from sleep
    virtual void wake();
    /// @brief Run the thread
    /// @return 0 on success, -1 on failure
    virtual int run();

    /// @brief Get the current replicate number
    virtual int getReplicate() {return replicate;}
    /// @brief Tell whether the previously running replicate has finished running
    virtual bool hasReplicateFinished();
    /// @brief Get the exit code for the last replicate to finish
    virtual int getReplicateExitCode();
	
protected:
    int replicate;
    MESolverFactory solverFactory;
    map<string,string> * parameters;
    lm::io::ReactionModel * reactionModel;
    lm::io::DiffusionModel * diffusionModel;
    uint8_t * lattice;
    size_t latticeSize;
    uint8_t * latticeSites;
    size_t latticeSitesSize;
    ResourceAllocator::ComputeResources resources;
    volatile bool replicateFinished;
    volatile int replicateExitCode;
};

}
}

#endif
