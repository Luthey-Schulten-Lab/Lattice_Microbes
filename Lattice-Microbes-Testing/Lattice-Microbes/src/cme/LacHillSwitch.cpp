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
#include <map>
#include <cmath>
#include "config.h"
#if defined(MACOSX)
#elif defined(LINUX)
#include <time.h>
#endif
#include "core/Math.h"
#include "core/Print.h"
#include "cme/FluctuatingNRSolver.h"
#include "cme/LacHillSwitch.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "lptf/Profile.h"

using std::string;
using std::map;

namespace lm {
namespace cme {

LacHillSwitch::LacHillSwitch()
{

}

LacHillSwitch::~LacHillSwitch()
{
}

void LacHillSwitch::setReactionModel(lm::io::ReactionModel * rm)
{
    FluctuatingNRSolver::setReactionModel(rm);

    // See if we should be using a fluctuating indcuer concentration.
    if (atoi((*parameters)["fluctuatingInducer"].c_str()))
    {
    	// Make the necessary adjustments to the fluctuating reactions.

    }

    if (atoi((*parameters)["looping"].c_str()))
    {
        Print::printf(Print::DEBUG, "Running looped lac switch.");
        // Set the protein limits.
        unsigned int Ymax=(unsigned int)atoi((*parameters)["Ymax"].c_str());
        unsigned int Ymin=(unsigned int)atoi((*parameters)["Ymin"].c_str());
        if (Ymax > 0) FluctuatingNRSolver::setSpeciesUpperLimit(4, Ymax);
        if (Ymin > 0) FluctuatingNRSolver::setSpeciesLowerLimit(4, Ymin);

        // Set the fpt tracking list.
        FluctuatingNRSolver::setFptTrackingList(list<uint>(1,4));
    }
    else
    {
        Print::printf(Print::DEBUG, "Running lac switch.");
        // Set the protein limits.
        unsigned int Ymax=(unsigned int)atoi((*parameters)["Ymax"].c_str());
        unsigned int Ymin=(unsigned int)atoi((*parameters)["Ymin"].c_str());
        if (Ymax > 0) FluctuatingNRSolver::setSpeciesUpperLimit(3, Ymax);
        if (Ymin > 0) FluctuatingNRSolver::setSpeciesLowerLimit(3, Ymin);

        // Set the fpt tracking list.
        FluctuatingNRSolver::setFptTrackingList(list<uint>(1,3));
    }

}

/*void LacHillSwitch::generateTrajectory()
{
    if (atoi((*parameters)["looping"].c_str()))
        generateTrajectoryLooping();
    else
        generateTrajectoryNoLooping();
}

void LacHillSwitch::generateTrajectoryNoLooping()
{
    uint ns         = 4;                                //I A M Y
    uint nr         = 6;
    uint c[]        = {  1,   0,   0,   0};
    uint o[]        = { 99,  99,   1,   1,   1,   1};
    double k[]      = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,
                         1,  -1,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,
                         0,   0,   0,   1,   0,  -1};
    uint C[]        = {  1,   0,   0,   0,   0,   0,
                         0,   1,   1,   0,   0,   0,
                         0,   0,   0,   1,   1,   0,
                         1,   1,   0,   0,   0,   1};

    // Get the transcription and translation parameters.
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    // Get the hill function parameters.
    double kap0min=atof((*parameters)["kap0min"].c_str());
    double kap0max=atof((*parameters)["kap0max"].c_str());
    double I050=atof((*parameters)["I050"].c_str());
    double h0=atof((*parameters)["h0"].c_str());
    double kap1min=atof((*parameters)["kap1min"].c_str());
    double kap1max=atof((*parameters)["kap1max"].c_str());
    double I150=atof((*parameters)["I150"].c_str());
    double h1=atof((*parameters)["h1"].c_str());

    // Get some other system parameters.
    double Iex=atof((*parameters)["Iex"].c_str());
    double kit=atof((*parameters)["kit"].c_str());
    double kid=atof((*parameters)["kid"].c_str());
    double KM=atof((*parameters)["KM"].c_str());
    double V=atof((*parameters)["V"].c_str());

    Print::printf(Print::DEBUG, "Running lac switch with a=%e b=%e kap0min=%e kap0max=%e I050=%e h0=%e kap1min=%e kap1max=%e I150=%e h1=%e gamma=%e d1=%e Iex=%e kit=%e kid=%e KM=%e V=%e", a, b, kap0min, kap0max, I050, h0, kap1min, kap1max, I150, h1, gamma, d1, Iex, kit, kid, KM, V);

    // Calculate the rate constants.
    double k0min=kap0min*d1, k0max=kap0max*d1, k1min=kap1min*d1, k1max=kap1max*d1, v0=a*d1, v1=b*gamma*d1, d0=gamma*d1;

    // Get the initial conditions.
    c[0] = 1;
    c[3] = (unsigned int)atoi((*parameters)["Y0"].c_str());

    // Set the rates.
    k[2] = v0;
    k[3] = v1;
    k[4] = d0;
    k[5] = d1;

    // Build the model.
    GillespieDSolver::buildModel(ns, nr, c, o, k, S, C);

    // Set the hill propensity functions.
    IHillPropensityArgs ka1(0,3,k0min,k0max,I050,Iex,kit,kid,KM,h0,V);
    GillespieDSolver::setModelPropensityFunction(0, &activationPropensity, &ka1);
    IHillPropensityArgs ka2(1,3,k1min,k1max,I150,Iex,kit,kid,KM,h1,V);
    GillespieDSolver::setModelPropensityFunction(1, &inactivationPropensity, &ka2);

    // Set the protein limits.
    unsigned int Ymax=(unsigned int)atoi((*parameters)["Ymax"].c_str());
    unsigned int Ymin=(unsigned int)atoi((*parameters)["Ymin"].c_str());
    if (Ymax > 0) GillespieDSolver::setSpeciesUpperLimit(3, Ymax);
    if (Ymin > 0) GillespieDSolver::setSpeciesLowerLimit(3, Ymin);

    // Set the fpt tracking list.
    GillespieDSolver::setFptTrackingList(list<uint>(1,3));

    // Generate the trajectory.
    GillespieDSolver::generateTrajectory();

    // Free the model since our propensity arg pointers are no longer valid.
    GillespieDSolver::destroyModel();
}

void LacHillSwitch::generateTrajectoryLooping()
{
    uint ns         = 5;                                //I A M Y L
    uint nr         = 9;
    uint c[]        = {  1,   0,   0,   0,   0};
    uint o[]        = { 99,  99,   1,   1,   1,   1,   1,  99,   1};
    double k[]      = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,  -1,   1,   0,
                         1,  -1,   0,   0,   0,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,   0,   0,   1,
                         0,   0,   0,   1,   0,  -1,   0,   0,   0,
                         0,   0,   0,   0,   0,   0,   1,  -1,   0};
    uint C[]        = {  1,   0,   0,   0,   0,   0,   1,   0,   0,
                         0,   1,   1,   0,   0,   0,   0,   0,   0,
                         0,   0,   0,   1,   1,   0,   0,   0,   0,
                         1,   1,   0,   0,   0,   1,   0,   1,   0,
                         0,   0,   0,   0,   0,   0,   0,   1,   1,};

    // Get the transcription and translation parameters.
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    // Get the hill function parameters.
    double kap0min=atof((*parameters)["kap0min"].c_str());
    double kap0max=atof((*parameters)["kap0max"].c_str());
    double I050=atof((*parameters)["I050"].c_str());
    double h0=atof((*parameters)["h0"].c_str());
    double kap1min=atof((*parameters)["kap1min"].c_str());
    double kap1max=atof((*parameters)["kap1max"].c_str());
    double I150=atof((*parameters)["I150"].c_str());
    double h1=atof((*parameters)["h1"].c_str());

    // Get looping function parameters.
    double kapfl=atof((*parameters)["kapfl"].c_str());
    double kaplfmin=atof((*parameters)["kaplfmin"].c_str());
    double kaplfmax=atof((*parameters)["kaplfmax"].c_str());
    double Ilf50=atof((*parameters)["Ilf50"].c_str());
    double hlf=atof((*parameters)["hlf"].c_str());
    double epsilon=atof((*parameters)["epsilon"].c_str());

    // Get some other system parameters.
    double Iex=atof((*parameters)["Iex"].c_str());
    double kit=atof((*parameters)["kit"].c_str());
    double kid=atof((*parameters)["kid"].c_str());
    double KM=atof((*parameters)["KM"].c_str());
    double V=atof((*parameters)["V"].c_str());

    Print::printf(Print::DEBUG, "Running lac switch with a=%e b=%e kap0min=%e kap0max=%e I050=%e h0=%e kap1min=%e kap1max=%e I150=%e h1=%e gamma=%e d1=%e kapfl=%e kaplfmin=%e kaplfmax=%e Ilf50=%e hlf=%e epsilon=%e Iex=%e kit=%e kid=%e KM=%e V=%e", a, b, kap0min, kap0max, I050, h0, kap1min, kap1max, I150, h1, gamma, d1, kapfl, kaplfmin, kaplfmax, Ilf50, hlf, epsilon, Iex, kit, kid, KM, V);

    // Calculate the rate constants.
    double k0min=kap0min*d1, k0max=kap0max*d1, k1min=kap1min*d1, k1max=kap1max*d1, v0=a*d1, v1=b*gamma*d1, d0=gamma*d1;
    double kfl=kapfl*d1, klfmin=kaplfmin*d1, klfmax=kaplfmax*d1;

    // Get the initial conditions.
    c[0] = (unsigned int)atoi((*parameters)["I0"].c_str());
    c[1] = (unsigned int)atoi((*parameters)["A0"].c_str());
    c[2] = (unsigned int)atoi((*parameters)["M0"].c_str());
    c[3] = (unsigned int)atoi((*parameters)["Y0"].c_str());
    c[4] = (unsigned int)atoi((*parameters)["L0"].c_str());

    // Set the rates.
    k[2] = v0;
    k[3] = v1;
    k[4] = d0;
    k[5] = d1;
    k[6] = kfl;
    k[8] = epsilon*v0;

    // Build the model.
    GillespieDSolver::buildModel(ns, nr, c, o, k, S, C);

    // Set the hill propensity functions.
    IHillPropensityArgs ka1(0,3,k0min,k0max,I050,Iex,kit,kid,KM,h0,V);
    GillespieDSolver::setModelPropensityFunction(0, &activationPropensity, &ka1);
    IHillPropensityArgs ka2(1,3,k1min,k1max,I150,Iex,kit,kid,KM,h1,V);
    GillespieDSolver::setModelPropensityFunction(1, &inactivationPropensity, &ka2);
    IHillPropensityArgs ka3(4,3,klfmin,klfmax,Ilf50,Iex,kit,kid,KM,hlf,V);
    GillespieDSolver::setModelPropensityFunction(7, &unloopingPropensity, &ka3);

    // Set the protein limits.
    unsigned int Ymax=(unsigned int)atoi((*parameters)["Ymax"].c_str());
    unsigned int Ymin=(unsigned int)atoi((*parameters)["Ymin"].c_str());
    if (Ymax > 0) GillespieDSolver::setSpeciesUpperLimit(3, Ymax);
    if (Ymin > 0) GillespieDSolver::setSpeciesLowerLimit(3, Ymin);

    // Set the fpt tracking list.
    GillespieDSolver::setFptTrackingList(list<uint>(1,3));

    // Generate the trajectory.
    GillespieDSolver::generateTrajectory();

    // Free the model since our propensity arg pointers are no longer valid.
    GillespieDSolver::destroyModel();
}

double LacHillSwitch::activationPropensity(double time, uint * speciesCounts, void * pargs)
{
    IHillPropensityArgs * args = (IHillPropensityArgs *)pargs;

    uint s = speciesCounts[args->si];
    if (s == 0) return 0.0;

    double x = (double)speciesCounts[args->xi];
    double xh = pow(1+(args->ITp*x),args->h);
    return args->kmin+((args->dk*xh)/(args->IRh+xh));
}

double LacHillSwitch::inactivationPropensity(double time, uint * speciesCounts, void * pargs)
{
    IHillPropensityArgs * args = (IHillPropensityArgs *)pargs;

    uint s = speciesCounts[args->si];
    if (s == 0) return 0.0;

    double x = (double)speciesCounts[args->xi];
    double xh = pow(1+(args->ITp*x),args->h);
    return args->kmax-((args->dk*xh)/(args->IRh+xh));
}

double LacHillSwitch::unloopingPropensity(double time, uint * speciesCounts, void * pargs)
{
    IHillPropensityArgs * args = (IHillPropensityArgs *)pargs;

    uint s = speciesCounts[args->si];
    if (s == 0) return 0.0;

    double x = (double)speciesCounts[args->xi];
    double xh = pow(1+(args->ITp*x),args->h);
    return args->kmin+((args->dk*xh)/(args->IRh+xh));
}*/

}
}
