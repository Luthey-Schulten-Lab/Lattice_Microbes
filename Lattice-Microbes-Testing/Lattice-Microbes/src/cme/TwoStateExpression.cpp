/*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
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

#include <string>
#include <map>
#include <cmath>
#include "config.h"
#if defined(MACOSX)
#elif defined(LINUX)
#include <time.h>
#endif
#include "core/Math.h"
#include "core/Print.h"
#include "cme/NextReactionSolver.h"
#include "cme/TwoStateExpression.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "lptf/Profile.h"

using std::string;
using std::map;

namespace lm {
namespace cme {

TwoStateExpression::TwoStateExpression()
{

}

TwoStateExpression::~TwoStateExpression()
{
}

void TwoStateExpression::generateTrajectory()
{
    uint ns         = 4;
    uint nr         = 6;
    uint c[]        = {  1,   0,   0,   0};
    uint o[]        = {  1,   1,   1,   1,   1,   1};
    double k[]      = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,
                         1,  -1,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,
                         0,   0,   0,   1,   0,  -1};
    uint C[]        = {  1,   0,   0,   0,   0,   0,
                         0,   1,   1,   0,   0,   0,
                         0,   0,   0,   1,   1,   0,
                         0,   0,   0,   0,   0,   1};

    // Get the simulation (*parameters).
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double kap0=atof((*parameters)["kap0"].c_str());
    double kap1=atof((*parameters)["kap1"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    // Calculate the rate constants,
    double k0=kap0*d1, k1=kap1*d1, v0=a*d1, v1=b*gamma*d1, d0=gamma*d1;

    // Get the initial conditions.
    c[0] = 1;
    c[3] = (unsigned int)atoi((*parameters)["initialProtein"].c_str());

    // Set the rates.
    k[0] = k0;
    k[1] = k1;
    k[2] = v0;
    k[3] = v1;
    k[4] = d0;
    k[5] = d1;

    Print::printf(Print::DEBUG, "Running two-state gene expression with a=%e b=%e kap0=%e kap1=%e gamma=%e d1=%e", a, b, kap0, kap1, gamma, d1);

    // Build the model.
    buildModel(ns, nr, c, o, k, S, C);

    // Set the protein limits.
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());
    if (maxProtein > 0) setSpeciesUpperLimit(3, maxProtein);
    if (minProtein > 0) setSpeciesLowerLimit(3, minProtein);

    // Set the fpt tracking list.
    setFptTrackingList(list<uint>(1,3));

    // Generate the trajectory.
    GillespieDSolver::generateTrajectory();
}

}
}
