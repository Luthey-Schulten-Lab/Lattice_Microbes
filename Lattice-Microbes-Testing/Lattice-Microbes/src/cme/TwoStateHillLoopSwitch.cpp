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
#include "cme/TwoStateHillLoopSwitch.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "lptf/Profile.h"

using std::string;
using std::map;

namespace lm {
namespace cme {

void TwoStateHillLoopSwitch::generateTrajectory()
{
    // Two-state (*parameters).
    double a=atof((*parameters)["a"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double kl=atof((*parameters)["kl"].c_str());
    double kh=atof((*parameters)["kh"].c_str());
    double h=atof((*parameters)["h"].c_str());
    double n50=atof((*parameters)["n50"].c_str());
    double klL=atof((*parameters)["klL"].c_str());
    double khL=atof((*parameters)["khL"].c_str());
    double kfl0=atof((*parameters)["kfl0"].c_str());
    double hL=atof((*parameters)["hL"].c_str());
    double n50L=atof((*parameters)["n50L"].c_str());
    double epsilon=atof((*parameters)["epsilon"].c_str());
    double d1=atof((*parameters)["d1"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());

    double writeInterval=atof((*parameters)["writeInterval"].c_str());
    double maxTime=atof((*parameters)["maxTime"].c_str());
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());

    // Calculate the rate constants.
    double v0=a*d1, v0L=epsilon*a*d1, v1=b*gamma*d1, d0=gamma*d1;

    // Get the initial conditions.
    const uint numberSpecies = 5;
    unsigned int speciesCounts[numberSpecies];  // L I A M N
    speciesCounts[0] = 1;
    speciesCounts[1] = 0;
    speciesCounts[2] = 0;
    speciesCounts[3] = 0;
    speciesCounts[4] = (unsigned int)atoi((*parameters)["initialProtein"].c_str());

    // Create the species counts data set to track during the simulation.
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpecies);
    speciesCountsDataSet.set_number_entries(1);
    speciesCountsDataSet.add_time(0.0);
    for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    double nextSpeciesCountsWriteTime = writeInterval;

    // Track the first passage times.
    unsigned int minValueAchieved = speciesCounts[4];
    unsigned int maxValueAchieved = speciesCounts[4];
    lm::io::FirstPassageTimes fptDataSet;
    fptDataSet.set_species(4);
    fptDataSet.add_species_count(speciesCounts[4]);
    fptDataSet.add_first_passage_time(0.0);

    Print::printf(Print::DEBUG, "Running simulation with v0=%e v0L=%e v1=%e d0=%e d1=%e", v0, v0L, v1, d0, d1);

    // Run the direct method.
    double time = 0.0;
    double propensities[9];
    double totalPropensity;
    unsigned long long steps=0;

    // Local cache of random numbers.
    int rngNext=TUNE_LOCAL_RNG_CACHE_SIZE;
    double rngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
    int expRngNext=TUNE_LOCAL_RNG_CACHE_SIZE;
    double expRngValues[TUNE_LOCAL_RNG_CACHE_SIZE];

    // Calculate the initial propensities.
    //0: I -> A
    propensities[0] = ((double)speciesCounts[1])*KHillPropensity(((double)speciesCounts[4]), kl, kh, n50, h);

    //1: A -> I
    propensities[1] = ((double)speciesCounts[2])*negativeKHillPropensity(((double)speciesCounts[4]), kl, kh, n50, h);

    //2: A -> A + M
    propensities[2] = ((double)speciesCounts[2])*v0;

    //3: M -> M + N
    propensities[3] = ((double)speciesCounts[3])*v1;

    //4: M -> 0
    propensities[4] = ((double)speciesCounts[3])*d0;

    //5: N -> 0
    propensities[5] = ((double)speciesCounts[4])*d1;

    //6: I -> L
    propensities[6] = ((double)speciesCounts[1])*kfl0;

    //7: L -> I
    propensities[7] = ((double)speciesCounts[0])*KHillPropensity(((double)speciesCounts[4]), klL, khL, n50L, hL);

    //8: L -> L + M
    propensities[8] = ((double)speciesCounts[0])*v0L;

    // Sum the total propensity.
    totalPropensity = propensities[0]+propensities[1]+propensities[2]+propensities[3]+propensities[4]+propensities[5]+propensities[6]+propensities[7]+propensities[8];

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
        for (; i<(9-1); i++)
        {
            if (r < propensities[i])
                break;
            else
                r -= propensities[i];
        }

        // 0:L  1:I  2:A  3:M  4:N
        // 0: I -> A
        // 1: A -> I
        // 2: A -> A + M
        // 3: M -> M + N
        // 4: M -> 0
        // 5: N -> 0
        // 6: I -> L
        // 7: L -> I
        // 8: L -> L + M
        // Update the species counts;
        switch (i)
        {
        case 0:
            speciesCounts[1]--;
            speciesCounts[2]++;
            break;
        case 1:
            speciesCounts[1]++;
            speciesCounts[2]--;
            break;
        case 2:
            speciesCounts[3]++;
            break;
        case 3:
            speciesCounts[4]++;
            break;
        case 4:
            speciesCounts[3]--;
            break;
        case 5:
            speciesCounts[4]--;
            break;
        case 6:
            speciesCounts[0]++;
            speciesCounts[1]--;
            break;
        case 7:
            speciesCounts[0]--;
            speciesCounts[1]++;
            break;
        case 8:
            speciesCounts[3]++;
            break;
        }

        // Update the propensities.
        //0: I -> A
        propensities[0] = ((double)speciesCounts[1])*KHillPropensity(((double)speciesCounts[4]), kl, kh, n50, h);

        //1: A -> I
        propensities[1] = ((double)speciesCounts[2])*negativeKHillPropensity(((double)speciesCounts[4]), kl, kh, n50, h);

        //2: A -> A + M
        propensities[2] = ((double)speciesCounts[2])*v0;

        //3: M -> M + N
        propensities[3] = ((double)speciesCounts[3])*v1;

        //4: M -> 0
        propensities[4] = ((double)speciesCounts[3])*d0;

        //5: N -> 0
        propensities[5] = ((double)speciesCounts[4])*d1;

        //6: I -> L
        propensities[6] = ((double)speciesCounts[1])*kfl0;

        //7: L -> I
        propensities[7] = ((double)speciesCounts[0])*KHillPropensity(((double)speciesCounts[4]), klL, khL, n50L, hL);

        //8: L -> L + M
        propensities[8] = ((double)speciesCounts[0])*v0L;

        // Sum the total propensity.
        totalPropensity = propensities[0]+propensities[1]+propensities[2]+propensities[3]+propensities[4]+propensities[5]+propensities[6]+propensities[7]+propensities[8];


        // Update the first passage time table.
        while (speciesCounts[4] < minValueAchieved)
        {
            fptDataSet.add_species_count(--minValueAchieved);
            fptDataSet.add_first_passage_time(time);
            addedFpt = true;
        }
        while (speciesCounts[4] > maxValueAchieved)
        {
            fptDataSet.add_species_count(++maxValueAchieved);
            fptDataSet.add_first_passage_time(time);
            addedFpt = true;
        }

        // If a species is past the limit, we are done.
        if (maxProtein > 0 && speciesCounts[4] >= maxProtein)
        {
            break;
        }
        else if (minProtein > 0 && speciesCounts[4] <= minProtein)
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
            fptDataSet.set_species(4);
        }
    }

    Print::printf(Print::DEBUG, "Generated trajectory for replicate %d in %llu steps.", replicate, steps);

    // If we finished the total time or ran out of reactions, write out the remaining time steps.
    if (time >= maxTime || totalPropensity <= 0)
    {
        while (nextSpeciesCountsWriteTime <= (maxTime+1e-9))
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
    }
}

}
}
