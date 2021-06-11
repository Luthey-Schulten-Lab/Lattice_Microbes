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

#ifndef LM_ME_MESOLVER_H
#define LM_ME_MESOLVER_H

#include <map>
#include "core/ResourceAllocator.h"

using std::map;
using lm::main::ResourceAllocator;

namespace lm {

namespace io {
class ReactionModel;
}
namespace me {

/// @class MESolver
/// @brief An abstract base class for all Master Equation solvers, this is essentially a representation of "the simulation instance"
class MESolver
{

public:
    /// @brief Create the MESolver
    MESolver();
    virtual ~MESolver();
    
    /// @brief Initialize the simulation
    /// @param replicate Replicate number out of total replicates
    /// @param parameters A map of all the parameters for the simulation
    /// @param A list of resources assigned to the simulation
    virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources)=0;
    
    /// @brief Tells whether the solver needs a reaction model
    virtual bool needsReactionModel()=0;
    /// @brief Tells whether the solver needs a reaction model
    virtual bool needsDiffusionModel()=0;
    /// @brief Actually run the simulation
    virtual void generateTrajectory()=0;
};

}
}

#endif
