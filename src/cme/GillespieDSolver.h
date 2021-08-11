/*
 * University of Illinois Open Source License
 * Copyright 2010-2018 Luthey-Schulten Group,
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

#ifndef LM_CME_GILLESPIEDSOLVER_H_
#define LM_CME_GILLESPIEDSOLVER_H_

#include <map>
#include <list>
#include <string>
#include "cme/CMESolver.h"
#include "FirstPassageTimes.pb.h"
#include "core/ResourceAllocator.h"
#include "rng/RandomGenerator.h"

using std::map;
using std::list;
using std::string;
using lm::main::ResourceAllocator;
using lm::rng::RandomGenerator;

namespace lm {
namespace cme {

class GillespieDSolver : public CMESolver
{
public:
    GillespieDSolver();
    virtual ~GillespieDSolver();
    virtual bool needsReactionModel() {return true;}
    virtual bool needsDiffusionModel()  {return false;}
    virtual void buildModel(const uint numberSpecies, const uint numberReactions, const uint * initialSpeciesCounts, const uint * reactionType, const double * k, const int * S, const uint * D, const uint kCols=1);
    virtual void generateTrajectory();

protected:
    virtual void destroyModel();
    inline void updateAllPropensities(double time);
    inline void updatePropensities(double time, uint r);

protected:
    double * propensities;

protected:
    virtual int hookSimulation(double time);
    virtual int onBeginTrajectory();
    virtual int onEndTrajectory();
};

}
}

#endif
