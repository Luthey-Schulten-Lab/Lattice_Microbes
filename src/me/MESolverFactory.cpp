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

#include <string>
#include "config.h"
#include "core/Exceptions.h"
#include "cme/FluctuatingNRSolver.h"
#include "cme/GillespieDSolver.h"
#include "cme/HillSwitch.h"
#include "cme/LacHillSwitch.h"
#include "cme/SelfRegulatingGeneSwitch.h"
#include "cme/TwoStateExpression.h"
#include "cme/TwoStateHillSwitch.h"
#include "cme/TwoStateHillLoopSwitch.h"
#include "cme/GillespieDSolver.h"
#include "me/MESolver.h"
#include "me/MESolverFactory.h"
#ifdef OPT_CUDA
#include "rdme/MpdRdmeSolver.h"
#include "rdme/MGPUMpdRdmeSolver.h"
#include "rdme/IntMpdRdmeSolver.h"
#ifdef OPT_MPI
#include "rdme/MPIMpdRdmeSolver.h"
#endif
#endif
#include "rdme/NextSubvolumeSolver.h"

using std::string;

namespace lm {
namespace me {

MESolverFactory::MESolverFactory()
:solver("")
{
}

void MESolverFactory::setSolver(string solver)
{
    if (solver == "lm::cme::FluctuatingNRSolver" || \
        solver == "lm::cme::GillespieDSolver" || \
        solver == "lm::cme::HillSwitch" || \
        solver == "lm::cme::LacHillSwitch" || \
        solver == "lm::cme::NextReactionSolver" || \
        solver == "lm::cme::SelfRegulatingGeneSwitch" || \
        solver == "lm::cme::TwoStateExpression" || \
        solver == "lm::cme::TwoStateHillSwitch" || \
        solver == "lm::cme::TwoStateHillLoopSwitch" || \
        solver == "lm::rdme::MpdRdmeSolver" || \
        solver == "lm::rdme::MGPUMpdRdmeSolver" || \
        solver == "lm::rdme::MPIMpdRdmeSolver" || \
        solver == "lm::rdme::IntMpdRdmeSolver" || \
        solver == "lm::rdme::NextSubvolumeSolver")
    {
        this->solver = solver;
        return;
    }
    throw lm::InvalidArgException("solver", "The specified solver is not known", solver.c_str());
}

bool MESolverFactory::needsReactionModel()
{
    MESolver * solver = instantiate();
    bool ret = solver->needsReactionModel();
    delete solver;
    return ret;
}

bool MESolverFactory::needsDiffusionModel()
{
    MESolver * solver = instantiate();
    bool ret = solver->needsDiffusionModel();
    delete solver;
    return ret;
}

MESolver * MESolverFactory::instantiate()
{
    // Run the simulation using the specified model.
    if (solver == "lm::cme::FluctuatingNRSolver")
    {
        return new lm::cme::FluctuatingNRSolver;
    }
    else if (solver == "lm::cme::GillespieDSolver")
    {
        return new lm::cme::GillespieDSolver;
    }
    else if (solver == "lm::cme::HillSwitch")
    {
        return new lm::cme::HillSwitch;
    }
    else if (solver == "lm::cme::LacHillSwitch")
    {
        return new lm::cme::LacHillSwitch;
    }
    else if (solver == "lm::cme::NextReactionSolver")
    {
        return new lm::cme::NextReactionSolver;
    }
    else if (solver == "lm::cme::SelfRegulatingGeneSwitch")
    {
        return new lm::cme::SelfRegulatingGeneSwitch;
    }
    else if (solver == "lm::cme::TwoStateExpression")
    {
        return new lm::cme::TwoStateExpression;
    }
    else if (solver == "lm::cme::TwoStateHillLoopSwitch")
    {
        return new lm::cme::TwoStateHillLoopSwitch;
    }
    else if (solver == "lm::cme::TwoStateHillSwitch")
    {
        return new lm::cme::TwoStateHillSwitch;
    }
    #ifdef OPT_CUDA
    else if (solver == "lm::rdme::MpdRdmeSolver")
    {
        return new lm::rdme::MpdRdmeSolver;
    }
	else if (solver == "lm::rdme::MGPUMpdRdmeSolver")
	{
		return new lm::rdme::MGPUMpdRdmeSolver;
	}
	else if (solver == "lm::rdme::IntMpdRdmeSolver")
	{
		return new lm::rdme::IntMpdRdmeSolver;
	}
	#ifdef OPT_MPI
	else if (solver == "lm::rdme::MPIMpdRdmeSolver")
	{
		return new lm::rdme::MPIMpdRdmeSolver;
	}
	#endif
    #else
    else if (solver == "lm::rdme::MpdRdmeSolver")
    {
        throw lm::Exception("The specified solver is only available using CUDA:", solver.c_str());
    }
    #endif
    else if (solver == "lm::rdme::NextSubvolumeSolver")
    {
        return new lm::rdme::NextSubvolumeSolver;
    }
    else
    {
        throw lm::Exception("The specified solver is unknown:", solver.c_str());
    }
}

}
}
