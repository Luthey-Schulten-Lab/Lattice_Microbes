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

#include <string>
#include <list>
#include <map>
#include <utility>
#include <cmath>
#include "config.h"
#if defined(MACOSX)
#elif defined(LINUX)
#include <time.h>
#endif
#include "core/Math.h"
#include "core/Print.h"
#include "cme/GillespieDSolver.h"
#include "cme/SelfRegulatingGeneSwitch.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "lptf/Profile.h"

using std::pair;
using std::string;
using std::list;
using std::map;

namespace lm {
namespace cme {

SelfRegulatingGeneSwitch::SelfRegulatingGeneSwitch()
{

}

SelfRegulatingGeneSwitch::~SelfRegulatingGeneSwitch()
{
}

void SelfRegulatingGeneSwitch::generateTrajectory()
{
    uint ns         = 3;
    uint nr         = 3;
    uint c[]        = {  0,   0,   1};
    uint o[]        = { 99,   1,   1};
    double k[]      = {0.0, 1.0, 0.0};
    int S[]         = {  1,  -1,   0,       //
                         0,   0,   1,
                         0,   0,   0};
    uint C[]        = {  1,   1,   0,
                         1,   0,   0,
                         0,   0,   1};

    // Get the simulation (*parameters).
    double a0=atof((*parameters)["a0"].c_str());
    double a1=atof((*parameters)["a1"].c_str());
    double p0=atof((*parameters)["p0"].c_str());
    double nvar=atof((*parameters)["nvar"].c_str());
    double tau=atof((*parameters)["tau"].c_str());
    double tauRecalc=atof((*parameters)["tau_recalc"].c_str());

    // Get the initial conditions.
    c[0] = (unsigned int)atoi((*parameters)["initialProtein"].c_str());

    Print::printf(Print::INFO, "Running self regulating gene switch with a0=%e a1=%e p0=%e nvar=%e tau=%e tau_recalc=%e", a0, a1, p0, nvar, tau, tauRecalc);

    // Calculate the rate constant for the ou update pseudo-reaction.
    if (tauRecalc > 0.0) k[2] = 1/tauRecalc;

    // Build the model.
    buildModel(ns, nr, c, o, k, S, C);

    // Set the ou propensity function.
    OUKHillPropensityArgs ka1(0,1,a0,a1,p0, nvar, tau, rng);
    setModelPropensityFunction(0, &ouKHillPropensity, &ka1);

    // Set the protein limits.
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());
    if (maxProtein > 0) setSpeciesUpperLimit(0, maxProtein);
    if (minProtein > 0) setSpeciesLowerLimit(0, minProtein);

    // Set the fpt tracking list.
    setFptTrackingList(list<uint>(1,0));

    // Turn on tracking of the noise.
    addToParameterTrackingList(pair<string,double*>("noise",&(ka1.noise)));

    // Generate the trajectory.
    NextReactionSolver::generateTrajectory();

    // Free the model since our propensity arg pointers are no longer valid.
    destroyModel();
}


double SelfRegulatingGeneSwitch::ouKHillPropensity(double time, uint * speciesCounts, void * pargs)
{
    OUKHillPropensityArgs * args = (OUKHillPropensityArgs *)pargs;

    // If we are using a noisy rate, calculate the new noise term.
    if (args->noiseVariance > 0.0)
    {
        // Get the current jump number.
        uint ouJumpNumber = speciesCounts[args->oui];

        // If this is the start, choose a properly distributed noise term to begin with.
        if (time == 0.0)
        {
            if (args->rngNext >= TUNE_LOCAL_RNG_CACHE_SIZE)
            {
                args->rng->getNormRandomDoubles(args->normRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                args->rngNext=0;
            }
            args->noise = sqrt(args->noiseVariance)*args->normRngValues[args->rngNext++];
            args->previousTime = time;
            args->lastOUJumpNumber = ouJumpNumber;
            Print::printf(Print::VERBOSE_DEBUG, "Initialized OU noise at %0.4e (JN %d): %0.4e", time, ouJumpNumber, args->noise);
        }

        // Otherwise, see if it is time to update the noise term.
        else if (ouJumpNumber != args->lastOUJumpNumber)
        {
            if (args->rngNext >= TUNE_LOCAL_RNG_CACHE_SIZE)
            {
                args->rng->getNormRandomDoubles(args->normRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                args->rngNext=0;
            }

            // Update the noise as an Ornstein-Uhlenbeck process to get the right autocorrelation.
            double dt = time-args->previousTime;
            double mu = exp(-dt/args->noiseTau);
            double sigma = sqrt(args->noiseVariance*(1-(mu*mu)));
            args->noise = args->noise*mu + sigma*args->normRngValues[args->rngNext++];
            args->previousTime = time;
            args->lastOUJumpNumber = ouJumpNumber;
            Print::printf(Print::VERBOSE_DEBUG, "Recalculated OU noise at %0.4e (JN %d): %0.4e", time, ouJumpNumber, args->noise);
        }
    }

    double x = (double)speciesCounts[args->xi];
    double xsq = x*x;
    double p = args->kmin+((args->dk*xsq/(args->x50sq+xsq)))+args->noise;
    Print::printf(Print::VERBOSE_DEBUG, "Recalculated OU KHill Propensity at %0.4e: %0.4e (noise %0.4e)", time, p, args->noise);
    return p>0.0?p:0.0;
}

}
}
