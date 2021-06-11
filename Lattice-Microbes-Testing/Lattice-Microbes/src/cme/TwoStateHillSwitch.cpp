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
#include "cme/GillespieDSolver.h"
#include "cme/TwoStateHillSwitch.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "lptf/Profile.h"

using std::string;
using std::map;

namespace lm {
namespace cme {

TwoStateHillSwitch::TwoStateHillSwitch()
{

}

TwoStateHillSwitch::~TwoStateHillSwitch()
{
}

void TwoStateHillSwitch::generateTrajectory()
{
    uint ns         = 4;
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

    // Get the simulation (*parameters).
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double kap0=atof((*parameters)["kap0"].c_str());
    double kap00=atof((*parameters)["kap00"].c_str());
    double n050=atof((*parameters)["n050"].c_str());
    double h0=atof((*parameters)["h0"].c_str());
    double kap1=atof((*parameters)["kap1"].c_str());
    double kap11=atof((*parameters)["kap11"].c_str());
    double n150=atof((*parameters)["n150"].c_str());
    double h1=atof((*parameters)["h1"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    // Calculate the rate constants,
    double k0=kap0*d1, k00=kap00*d1, k1=kap1*d1, k11=kap11*d1, v0=a*d1, v1=b*gamma*d1, d0=gamma*d1;

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

    Print::printf(Print::DEBUG, "Running two-state Hill switch with a=%e b=%e kap0=%e kap00=%e h0=%e n050=%e kap11=%e kap1=%e h1=%e n150=%e gamma=%e d1=%e", a, b, kap0, kap00, h0, n050, kap11, kap1, h1, n150, gamma, d1);

    GillespieDSolver::buildModel(ns, nr, c, o, k, S, C);

    list<FirstOrderPropensityArgs *> kHillArgs;
    KHillPropensityArgs ka1(0,3,k0,k00,n050,h0);
    GillespieDSolver::setModelPropensityFunction(0, &kHillPropensity, &ka1);
    KHillPropensityArgs ka2(1,3,k11,k1,n150,h1);
    GillespieDSolver::setModelPropensityFunction(1, &negativeKHillPropensity, &ka2);

    // Set the protein limits.
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());
    if (maxProtein > 0) GillespieDSolver::setSpeciesUpperLimit(3, maxProtein);
    if (minProtein > 0) GillespieDSolver::setSpeciesLowerLimit(3, minProtein);

    // Set the fpt tracking list.
    GillespieDSolver::setFptTrackingList(list<uint>(1,3));

    // Generate the trajectory.
    GillespieDSolver::generateTrajectory();

    // Free the model since our propensity arg pointers are no longer valid.
    GillespieDSolver::destroyModel();
}


double TwoStateHillSwitch::kHillPropensity(double time, uint * speciesCounts, void * pargs)
{
    KHillPropensityArgs * args = (KHillPropensityArgs *)pargs;

    uint s = speciesCounts[args->si];
    if (s == 0) return 0.0;

    double x = (double)speciesCounts[args->xi];
    double xph = pow(x,args->h);
    return args->kmin+((args->dk*xph/(args->x50ph+xph)));
}

double TwoStateHillSwitch::negativeKHillPropensity(double time, uint * speciesCounts, void * pargs)
{
    KHillPropensityArgs * args = (KHillPropensityArgs *)pargs;

    uint s = speciesCounts[args->si];
    if (s == 0) return 0.0;

    double x = (double)speciesCounts[args->xi];
    double xph = pow(x,args->h);
    return args->kmax-((args->dk*xph/(args->x50ph+xph)));
}

/*void TwoStateHillSwitch::generateTrajectory2()
{
    if ((*parameters)["proteinModel"] == "constant")
    {
        runConstantProtein();
        return;
    }
    else if ((*parameters)["proteinModel"] == "geometric")
    {
        runGeometricProtein();
        return;
    }

    // Get the simulation (*parameters).
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double kap0=atof((*parameters)["kap0"].c_str());
    double kap00=atof((*parameters)["kap00"].c_str());
    double n050=atof((*parameters)["n050"].c_str());
    double h0=atof((*parameters)["h0"].c_str());
    double kap1=atof((*parameters)["kap1"].c_str());
    double kap11=atof((*parameters)["kap11"].c_str());
    double n150=atof((*parameters)["n150"].c_str());
    double h1=atof((*parameters)["h1"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    double writeInterval=atof((*parameters)["writeInterval"].c_str());
    double maxTime=atof((*parameters)["maxTime"].c_str());
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());

    // Calculate the rate constants.
    double k0=kap0*d1, k00=kap00*d1, k1=kap1*d1, k11=kap11*d1, v0=a*d1, v1=b*gamma*d1, d0=gamma*d1;

    // Get the initial conditions.
    const uint numberSpecies = 4;
    unsigned int speciesCounts[numberSpecies];  // I A M Y
    speciesCounts[0] = 1;
    speciesCounts[1] = 0;
    speciesCounts[2] = 0;
    speciesCounts[3] = (unsigned int)atoi((*parameters)["initialProtein"].c_str());

    // Create the species counts data set to track during the simulation.
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpecies);
    speciesCountsDataSet.set_number_entries(1);
    speciesCountsDataSet.add_time(0.0);
    for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    double nextSpeciesCountsWriteTime = writeInterval;

    // Track the first passage times.
    unsigned int minValueAchieved = speciesCounts[3];
    unsigned int maxValueAchieved = speciesCounts[3];
    lm::io::FirstPassageTimes fptDataSet;
    fptDataSet.set_species(3);
    fptDataSet.add_species_count(speciesCounts[3]);
    fptDataSet.add_first_passage_time(0.0);

    Print::printf(Print::DEBUG, "Running simulation with a=%e b=%e kap0=%e kap00=%e h0=%e n050=%e kap11=%e kap1=%e h1=%e n150=%e gamma=%e d1=%e", a, b, kap0, kap00, h0, n050, kap11, kap1, h1, n150, gamma, d1);

    // Run the direct method.
    double time = 0.0;
    double propensities[6];
    double totalPropensity;
    unsigned long long steps=0;

    // Local cache of random numbers.
    int rngNext=TUNE_LOCAL_RNG_CACHE_SIZE;
    double rngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
    int expRngNext=TUNE_LOCAL_RNG_CACHE_SIZE;
    double expRngValues[TUNE_LOCAL_RNG_CACHE_SIZE];

    // Calculate the initial propensities.
    propensities[0] = ((double)speciesCounts[0])*KHillPropensity(((double)speciesCounts[3]), k0, k00, n050, h0);
    propensities[1] = ((double)speciesCounts[1])*negativeKHillPropensity(((double)speciesCounts[3]), k11, k1, n150, h1);
    propensities[2] = ((double)speciesCounts[1])*v0;
    propensities[3] = ((double)speciesCounts[2])*v1;
    propensities[4] = ((double)speciesCounts[2])*d0;
    propensities[5] = ((double)speciesCounts[3])*d1;
    totalPropensity = propensities[0]+propensities[1]+propensities[2]+propensities[3]+propensities[4]+propensities[5];

    PROF_BEGIN(PROF_SIM_EXECUTE);
    bool addedSpeciesCounts;
    bool addedFpt;
    while (totalPropensity > 0)
    {
        addedSpeciesCounts = false;
        addedFpt = false;
        steps++;

        // See if we need to update our rng caches.
        if (rngNext>=TUNE_LOCAL_RNG_CACHE_SIZE)
        {
            if (steps < 20000)
            {
                PROF_BEGIN(PROF_CACHE_RNG);
            }
            rng->getRandomDoubles(rngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            rngNext=0;
            if (steps < 20000)
            {
                PROF_END(PROF_CACHE_RNG);
            }
        }
        if (expRngNext>=TUNE_LOCAL_RNG_CACHE_SIZE)
        {
            if (steps < 20000)
            {
                PROF_BEGIN(PROF_CACHE_EXP_RNG);
            }
            rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            expRngNext=0;
            if (steps < 20000)
            {
                PROF_END(PROF_CACHE_EXP_RNG);
            }
        }

        // Calculate the time to the next reaction.
        double expR = expRngValues[expRngNext++];
        time += expR/totalPropensity;

        // If the new time is past the end time, we are done.
        if (time >= maxTime)
        {
            time = maxTime;
            break;
        }

        // See if the new time is past the write time.
        while (nextSpeciesCountsWriteTime <= (time+1e-9))
        {
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
            nextSpeciesCountsWriteTime += writeInterval;
            addedSpeciesCounts = true;
        }

        // Calculate which reaction it was.
        double r = rngValues[rngNext++]*totalPropensity;
        int i=0;
        for (; i<(6-1); i++)
        {
            if (r < propensities[i])
                break;
            else
                r -= propensities[i];
        }

        // Update the species counts and propensities.
        switch (i)
        {

        // Switch to active.
        case 0:
            speciesCounts[0] = 0;
            speciesCounts[1] = 1;
            propensities[0] = 0.0;
            propensities[1] = negativeKHillPropensity(((double)speciesCounts[3]), k11, k1, n150, h1);
            propensities[2] = v0;
            break;

        // Switch to inactive.
        case 1:
            speciesCounts[0] = 1;
            speciesCounts[1] = 0;
            propensities[0] = KHillPropensity(((double)speciesCounts[3]), k0, k00, n050, h0);
            propensities[1] = 0.0;
            propensities[2] = 0.0;
            break;

        // mRNA produced.
        case 2:
            speciesCounts[2]++;
            propensities[3] = ((double)speciesCounts[2])*v1;
            propensities[4] = ((double)speciesCounts[2])*d0;
            break;

        // Protein produced.
        case 3:
            speciesCounts[3]++;
            propensities[0] = (speciesCounts[0]==1)?(KHillPropensity(((double)speciesCounts[3]), k0, k00, n050, h0)):(0.0);
            propensities[1] = (speciesCounts[1]==1)?(negativeKHillPropensity(((double)speciesCounts[3]), k11, k1, n150, h1)):(0.0);
            propensities[5] = ((double)speciesCounts[3])*d1;
            break;

        // mRNA decayed.
        case 4:
            speciesCounts[2]--;
            propensities[3] = ((double)speciesCounts[2])*v1;
            propensities[4] = ((double)speciesCounts[2])*d0;
            break;

        // Protein decayed.
        case 5:
            speciesCounts[3]--;
            propensities[0] = (speciesCounts[0]==1)?(KHillPropensity(((double)speciesCounts[3]), k0, k00, n050, h0)):(0.0);
            propensities[1] = (speciesCounts[1]==1)?(negativeKHillPropensity(((double)speciesCounts[3]), k11, k1, n150, h1)):(0.0);
            propensities[5] = ((double)speciesCounts[3])*d1;
            break;
        }

        // Recalculate the total propensity.
        totalPropensity = propensities[0]+propensities[1]+propensities[2]+propensities[3]+propensities[4]+propensities[5];

        // Update the first passage time table.
        while (speciesCounts[3] < minValueAchieved)
        {
            fptDataSet.add_species_count(--minValueAchieved);
            fptDataSet.add_first_passage_time(time);
            addedFpt = true;
        }
        while (speciesCounts[3] > maxValueAchieved)
        {
            fptDataSet.add_species_count(++maxValueAchieved);
            fptDataSet.add_first_passage_time(time);
            addedFpt = true;
        }

        // If a species is past the limit, we are done.
        if (maxProtein > 0 && speciesCounts[3] >= maxProtein)
        {
            break;
        }
        else if (minProtein > 0 && speciesCounts[3] <= minProtein)
        {
            break;
        }

        // See if we have accumulated enough data to send.
        if (addedSpeciesCounts && speciesCountsDataSet.number_entries() >= TUNE_SPECIES_COUNTS_BUFFER_SIZE)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_COUNTS);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
            PROF_END(PROF_SERIALIZE_COUNTS);

            // Reset the data set.
            speciesCountsDataSet.Clear();
            speciesCountsDataSet.set_number_species(numberSpecies);
            speciesCountsDataSet.set_number_entries(0);
        }

        if (addedFpt && fptDataSet.first_passage_time_size() >= TUNE_FIRST_PASSAGE_TIME_BUFFER_SIZE)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_FPT);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptDataSet);
            PROF_END(PROF_SERIALIZE_FPT);

            // Reset the data set.
            fptDataSet.Clear();
            fptDataSet.set_species(3);
        }
    }
    PROF_END(PROF_SIM_EXECUTE);

    Print::printf(Print::DEBUG, "Generated trajectory for replicate %d in %llu steps.", replicate, steps);

    // If we finished the total time or ran out of reactions, write out the remaining time steps.
    if (time >= maxTime || totalPropensity <= 0)
    {
        Print::printf(Print::DEBUG, "Finished with time %e (%e)", time, maxTime);
        while (nextSpeciesCountsWriteTime <= (maxTime+1e-9))
        {
            Print::printf(Print::DEBUG, "Recording event at time %e (%e)", nextSpeciesCountsWriteTime, maxTime);
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
            nextSpeciesCountsWriteTime += writeInterval;
        }
        Print::printf(Print::DEBUG, "Done recording events at time %e (%e)", nextSpeciesCountsWriteTime, maxTime);
    }

    // Otherwise we must have finished because of a species limit, so just write out the last time.
    else
    {
        // Record the species counts.
        speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
        speciesCountsDataSet.add_time(time);
        for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    }

    // Send any remaining species counts to the queue.
    if (speciesCountsDataSet.number_entries() > 0)
    {
        PROF_BEGIN(PROF_SERIALIZE_COUNTS);
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
        PROF_END(PROF_SERIALIZE_COUNTS);
    }

    // Send any remaining first passage times to the queue.
    if (fptDataSet.first_passage_time_size() > 0)
    {
        // Push it to the output queue.
        PROF_BEGIN(PROF_SERIALIZE_FPT);
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptDataSet);
        PROF_END(PROF_SERIALIZE_FPT);
    }
}*/

void TwoStateHillSwitch::runConstantProtein()
{
    /*// Get the simulation (*parameters).
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double kap0=atof((*parameters)["kap0"].c_str());
    double kap00=atof((*parameters)["kap00"].c_str());
    double n050=atof((*parameters)["n050"].c_str());
    double h0=atof((*parameters)["h0"].c_str());
    double kap1=atof((*parameters)["kap1"].c_str());
    double kap11=atof((*parameters)["kap11"].c_str());
    double n150=atof((*parameters)["n150"].c_str());
    double h1=atof((*parameters)["h1"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    double writeInterval=atof((*parameters)["writeInterval"].c_str());
    double maxTime=atof((*parameters)["maxTime"].c_str());
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());

    // Calculate the rate constants,
    double k0=kap0*d1, k00=kap00*d1, k1=kap1*d1, k11=kap11*d1, v0=a*d1;

    // Get the initial conditions.
    const uint numberSpecies = 4;
    unsigned int speciesCounts[numberSpecies];  // I A M Y
    speciesCounts[0] = 1;
    speciesCounts[1] = 0;
    speciesCounts[2] = 0;
    speciesCounts[3] = (unsigned int)atoi((*parameters)["initialProtein"].c_str());

    // Create the species counts data set to track during the simulation.
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpecies);
    speciesCountsDataSet.set_number_entries(1);
    speciesCountsDataSet.add_time(0.0);
    for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    double nextSpeciesCountsWriteTime = writeInterval;

    // Track the first passage times.
    unsigned int minValueAchieved = speciesCounts[3];
    unsigned int maxValueAchieved = speciesCounts[3];
    lm::io::FirstPassageTimes fptDataSet;
    fptDataSet.set_species(3);
    fptDataSet.add_species_count(speciesCounts[3]);
    fptDataSet.add_first_passage_time(0.0);

    Print::printf(Print::DEBUG, "Running constant protein simulation with a=%e b=%e kap0=%e kap00=%e h0=%e n050=%e kap11=%e kap1=%e h1=%e n150=%e gamma=%e d1=%e", a, b, kap0, kap00, h0, n050, kap11, kap1, h1, n150, gamma, d1);

    // Run the direct method.
    double time = 0.0;
    double propensities[6];
    double totalPropensity;
    while (true)
    {
        // Calculate the propensities.
        propensities[0] = ((double)speciesCounts[0])*KHillPropensity(((double)speciesCounts[3]), k0, k00, n050, h0);
        propensities[1] = ((double)speciesCounts[1])*negativeKHillPropensity(((double)speciesCounts[3]), k11, k1, n150, h1);
        propensities[2] = ((double)speciesCounts[1])*v0;
        propensities[3] = 0.0;
        propensities[4] = 0.0;
        propensities[5] = ((double)speciesCounts[3])*d1;
        totalPropensity = propensities[0]+propensities[1]+propensities[2]+propensities[3]+propensities[4]+propensities[5];

        // If the total propensity is zero, we are done.
        if (totalPropensity <= 0)
        {
            time = maxTime;
            break;
        }

        // Calculate the time to the next reaction.
        double r1 = rng->getRandomDoubleEE();
        time = time + (-log(r1))/totalPropensity;

        // If the new time is past the end time, we are done.
        if (time >= maxTime)
        {
            time = maxTime;
            break;
        }

        // See if the new time is past the write time.
        while (time >= (nextSpeciesCountsWriteTime-1e-9))
        {
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
            nextSpeciesCountsWriteTime += writeInterval;
        }

        // Calculate which reaction it was.
        double r2 = rng->getRandomDoubleIE()*totalPropensity;
        int i=0;
        for (; i<(6-1); i++)
        {
            if (r2 < propensities[i])
                break;
            else
                r2 -= propensities[i];
        }

        // Update the species counts;
        switch (i)
        {
        case 0:
            speciesCounts[0]--;
            speciesCounts[1]++;
            break;
        case 1:
            speciesCounts[0]++;
            speciesCounts[1]--;
            break;
        case 2:
            speciesCounts[3]+= (uint)floor(b);
            break;
        case 3:
            break;
        case 4:
            break;
        case 5:
            speciesCounts[3]--;
            break;
        }

        // Update the first passage time table.
        while (speciesCounts[3] < minValueAchieved)
        {
            fptDataSet.add_species_count(--minValueAchieved);
            fptDataSet.add_first_passage_time(time);
        }
        while (speciesCounts[3] > maxValueAchieved)
        {
            fptDataSet.add_species_count(++maxValueAchieved);
            fptDataSet.add_first_passage_time(time);
        }

        // If a species is past the limit, we are done.
        if (maxProtein > 0 && speciesCounts[3] >= maxProtein)
        {
            break;
        }
        else if (minProtein > 0 && speciesCounts[3] <= minProtein)
        {
            break;
        }

        // See if we have accumulated enough data to send.
        if (speciesCountsDataSet.number_entries() >= TUNE_SPECIES_COUNTS_BUFFER_SIZE)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_COUNTS);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
            PROF_END(PROF_SERIALIZE_COUNTS);

            // Reset the data set.
            speciesCountsDataSet.Clear();
            speciesCountsDataSet.set_number_species(numberSpecies);
            speciesCountsDataSet.set_number_entries(0);
        }

        if (fptDataSet.first_passage_time_size() >= TUNE_FIRST_PASSAGE_TIME_BUFFER_SIZE)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_FPT);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptDataSet);
            PROF_END(PROF_SERIALIZE_FPT);

            // Reset the data set.
            fptDataSet.Clear();
            fptDataSet.set_species(3);
        }
    }

    // If we finish the total time, write out the remaining time steps.
    if (time == maxTime)
    {
        while (time >= (nextSpeciesCountsWriteTime-1e-9))
        {
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
            nextSpeciesCountsWriteTime += writeInterval;
        }
    }

    // Otherwise we must have finished because of a species limit, so just write out the last time.
    else
    {
        // Record the species counts.
        speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
        speciesCountsDataSet.add_time(time);
        for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    }

    // Send any remaining species counts to the queue.
    if (speciesCountsDataSet.number_entries() > 0)
    {
        PROF_BEGIN(PROF_SERIALIZE_COUNTS);
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
        PROF_END(PROF_SERIALIZE_COUNTS);
    }

    // Send any remaining first passage times to the queue.
    if (fptDataSet.first_passage_time_size() > 0)
    {
        // Push it to the output queue.
        PROF_BEGIN(PROF_SERIALIZE_FPT);
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptDataSet);
        PROF_END(PROF_SERIALIZE_FPT);
    }*/
}

void TwoStateHillSwitch::runGeometricProtein()
{
    /*
    // Get the simulation (*parameters).
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double kap0=atof((*parameters)["kap0"].c_str());
    double kap00=atof((*parameters)["kap00"].c_str());
    double n050=atof((*parameters)["n050"].c_str());
    double h0=atof((*parameters)["h0"].c_str());
    double kap1=atof((*parameters)["kap1"].c_str());
    double kap11=atof((*parameters)["kap11"].c_str());
    double n150=atof((*parameters)["n150"].c_str());
    double h1=atof((*parameters)["h1"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    double writeInterval=atof((*parameters)["writeInterval"].c_str());
    double maxTime=atof((*parameters)["maxTime"].c_str());
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());

    // Calculate the rate constants,
    double k0=kap0*d1, k00=kap00*d1, k1=kap1*d1, k11=kap11*d1, v0=a*d1;

    // Get the initial conditions.
    const uint numberSpecies = 4;
    unsigned int speciesCounts[numberSpecies];  // I A M Y
    speciesCounts[0] = 1;
    speciesCounts[1] = 0;
    speciesCounts[2] = 0;
    speciesCounts[3] = (unsigned int)atoi((*parameters)["initialProtein"].c_str());

    // Create the species counts data set to track during the simulation.
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpecies);
    speciesCountsDataSet.set_number_entries(1);
    speciesCountsDataSet.add_time(0.0);
    for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    double nextSpeciesCountsWriteTime = writeInterval;

    // Track the first passage times.
    unsigned int minValueAchieved = speciesCounts[3];
    unsigned int maxValueAchieved = speciesCounts[3];
    lm::io::FirstPassageTimes fptDataSet;
    fptDataSet.set_species(3);
    fptDataSet.add_species_count(speciesCounts[3]);
    fptDataSet.add_first_passage_time(0.0);

    Print::printf(Print::DEBUG, "Running geometric protein simulation with a=%e b=%e kap0=%e kap00=%e h0=%e n050=%e kap11=%e kap1=%e h1=%e n150=%e gamma=%e d1=%e", a, b, kap0, kap00, h0, n050, kap11, kap1, h1, n150, gamma, d1);

    // Run the direct method.
    double time = 0.0;
    double propensities[6];
    double totalPropensity;
    while (true)
    {
        // Calculate the propensities.
        propensities[0] = ((double)speciesCounts[0])*KHillPropensity(((double)speciesCounts[3]), k0, k00, n050, h0);
        propensities[1] = ((double)speciesCounts[1])*negativeKHillPropensity(((double)speciesCounts[3]), k11, k1, n150, h1);
        propensities[2] = ((double)speciesCounts[1])*v0;
        propensities[3] = 0.0;
        propensities[4] = 0.0;
        propensities[5] = ((double)speciesCounts[3])*d1;
        totalPropensity = propensities[0]+propensities[1]+propensities[2]+propensities[3]+propensities[4]+propensities[5];

        // If the total propensity is zero, we are done.
        if (totalPropensity <= 0)
        {
            time = maxTime;
            break;
        }

        // Calculate the time to the next reaction.
        double r1 = rng->getRandomDoubleEE();
        time = time + (-log(r1))/totalPropensity;

        // If the new time is past the end time, we are done.
        if (time >= maxTime)
        {
            time = maxTime;
            break;
        }

        // See if the new time is past the write time.
        while (time >= (nextSpeciesCountsWriteTime-1e-9))
        {
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
            nextSpeciesCountsWriteTime += writeInterval;
        }

        // Calculate which reaction it was.
        double r2 = rng->getRandomDoubleIE()*totalPropensity;
        int i=0;
        for (; i<(6-1); i++)
        {
            if (r2 < propensities[i])
                break;
            else
                r2 -= propensities[i];
        }

        // Update the species counts;
        double r3;
        uint geoR;
        switch (i)
        {
        case 0:
            speciesCounts[0]--;
            speciesCounts[1]++;
            break;
        case 1:
            speciesCounts[0]++;
            speciesCounts[1]--;
            break;
        case 2:
            r3 = rng->getRandomDoubleEE();
            geoR = (uint)floor(log(r3)/log(1-(1/(b+1))));
            speciesCounts[3]+= geoR;
            break;
        case 3:
            break;
        case 4:
            break;
        case 5:
            speciesCounts[3]--;
            break;
        }

        // Update the first passage time table.
        while (speciesCounts[3] < minValueAchieved)
        {
            fptDataSet.add_species_count(--minValueAchieved);
            fptDataSet.add_first_passage_time(time);
        }
        while (speciesCounts[3] > maxValueAchieved)
        {
            fptDataSet.add_species_count(++maxValueAchieved);
            fptDataSet.add_first_passage_time(time);
        }

        // If a species is past the limit, we are done.
        if (maxProtein > 0 && speciesCounts[3] >= maxProtein)
        {
            break;
        }
        else if (minProtein > 0 && speciesCounts[3] <= minProtein)
        {
            break;
        }

        // See if we have accumulated enough data to send.
        if (speciesCountsDataSet.number_entries() >= TUNE_SPECIES_COUNTS_BUFFER_SIZE)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_COUNTS);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
            PROF_END(PROF_SERIALIZE_COUNTS);

            // Reset the data set.
            speciesCountsDataSet.Clear();
            speciesCountsDataSet.set_number_species(numberSpecies);
            speciesCountsDataSet.set_number_entries(0);
        }

        if (fptDataSet.first_passage_time_size() >= TUNE_FIRST_PASSAGE_TIME_BUFFER_SIZE)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_FPT);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptDataSet);
            PROF_END(PROF_SERIALIZE_FPT);

            // Reset the data set.
            fptDataSet.Clear();
            fptDataSet.set_species(3);
        }
    }

    // If we finish the total time, write out the remaining time steps.
    if (time == maxTime)
    {
        while (time >= (nextSpeciesCountsWriteTime-1e-9))
        {
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
            nextSpeciesCountsWriteTime += writeInterval;
        }
    }

    // Otherwise we must have finished because of a species limit, so just write out the last time.
    else
    {
        // Record the species counts.
        speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
        speciesCountsDataSet.add_time(time);
        for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    }

    // Send any remaining species counts to the queue.
    if (speciesCountsDataSet.number_entries() > 0)
    {
        PROF_BEGIN(PROF_SERIALIZE_COUNTS);
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
        PROF_END(PROF_SERIALIZE_COUNTS);
    }

    // Send any remaining first passage times to the queue.
    if (fptDataSet.first_passage_time_size() > 0)
    {
        // Push it to the output queue.
        PROF_BEGIN(PROF_SERIALIZE_FPT);
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptDataSet);
        PROF_END(PROF_SERIALIZE_FPT);
    }*/
}

}
}
