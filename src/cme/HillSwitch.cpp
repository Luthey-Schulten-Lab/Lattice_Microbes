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
#include "cme/GillespieDSolver.h"
#include "cme/HillSwitch.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "lptf/Profile.h"

using std::string;
using std::map;

namespace lm {
namespace cme {

void HillSwitch::generateTrajectory()
{
    // Get the simulation parameters.
    double a0=atof((*parameters)["a0"].c_str());
    double a00=atof((*parameters)["a00"].c_str());
    double n50=atof((*parameters)["n50"].c_str());
    double h=atof((*parameters)["h"].c_str());
    double b=atof((*parameters)["b"].c_str());
    double gamma=atof((*parameters)["gamma"].c_str());
    double d1=atof((*parameters)["d1"].c_str());

    double maxTime=atof((*parameters)["maxTime"].c_str());
    unsigned int maxProtein=(unsigned int)atoi((*parameters)["maxProtein"].c_str());
    unsigned int minProtein=(unsigned int)atoi((*parameters)["minProtein"].c_str());

    // Calculate the rate constants,
    double v0=a0*d1, v00=a00*d1, v1=b*gamma*d1, d0=gamma*d1;

    // Get the initial conditions.
    const uint numberSpecies = 2;
    const uint MRNA=0;
    const uint PROTEIN=1;
    unsigned int speciesCounts[numberSpecies];  // M Y
    speciesCounts[MRNA] = (unsigned int)atoi((*parameters)["initialMRNA"].c_str());
    speciesCounts[PROTEIN] = (unsigned int)atoi((*parameters)["initialProtein"].c_str());

    // Create the species counts data set to track during the simulation.
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpecies);
    speciesCountsDataSet.set_number_entries(1);
    speciesCountsDataSet.add_time(0.0);
    for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);

    // Get the interval for writing species counts.
    double writeInterval=atof((*parameters)["writeInterval"].c_str());
    bool writeTimeSteps = (writeInterval > 0.0);
    double nextSpeciesCountsWriteTime = writeInterval;

    // Track the first passage times.
    unsigned int minValueAchieved = speciesCounts[PROTEIN];
    unsigned int maxValueAchieved = speciesCounts[PROTEIN];
    lm::io::FirstPassageTimes fptDataSet;
    fptDataSet.set_species(PROTEIN);
    fptDataSet.add_species_count(speciesCounts[PROTEIN]);
    fptDataSet.add_first_passage_time(0.0);

    Print::printf(Print::DEBUG, "Running simulation with a0=%e a00=%e h=%e n50=%e b=%e gamma=%e d1=%e and write mode %d", a0, a00, h, n50, b, gamma, d1, writeTimeSteps);

    // Run the direct method.
    double time = 0.0;
    const uint numberReactions=4;
    double propensities[numberReactions];
    double totalPropensity;
    unsigned long long steps=0;
    while (true)
    {
        steps++;

        // Calculate the propensities.
        propensities[0] = KHillPropensity(((double)speciesCounts[PROTEIN]), v0, v00, n50, h);
        propensities[1] = ((double)speciesCounts[MRNA])*v1;
        propensities[2] = ((double)speciesCounts[MRNA])*d0;
        propensities[3] = ((double)speciesCounts[PROTEIN])*d1;
        totalPropensity = propensities[0]+propensities[1]+propensities[2]+propensities[3];

        // If the total propensity is zero, we are done.
        if (totalPropensity <= 0)
        {
            time = maxTime;
            break;
        }

        // Calculate the time to the next reaction.
        double expr = rng->getExpRandomDouble();
        time += expr/totalPropensity;

        // If the new time is past the end time, we are done.
        if (time >= maxTime)
        {
            time = maxTime;
            break;
        }

        // If we are writing time steps, write out any time steps before this event occurred.
        if (writeTimeSteps)
        {
            // Write time steps until the next write time is past the current time.
            while (time >= (nextSpeciesCountsWriteTime-1e-9))
            {
                // Record the species counts.
                speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
                speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
                for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
                nextSpeciesCountsWriteTime += writeInterval;
            }
        }

        // Calculate which reaction it was.
        double r = rng->getRandomDouble()*totalPropensity;
        uint i=0;
        for (; i<(numberReactions-1); i++)
        {
            if (r < propensities[i])
                break;
            else
                r -= propensities[i];
        }

        // Update the species counts.
        switch (i)
        {
        case 0:
            speciesCounts[MRNA]++;
            break;
        case 1:
            speciesCounts[PROTEIN]++;
            break;
        case 2:
            speciesCounts[MRNA]--;
            break;
        case 3:
            speciesCounts[PROTEIN]--;
            break;
        }

        // We are recording every event, so add it.
        if (!writeTimeSteps)
        {
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(time);
            for (uint j=0; j<numberSpecies; j++) speciesCountsDataSet.add_species_count(speciesCounts[j]);
        }

        // Update the first passage time table.
        while (speciesCounts[PROTEIN] < minValueAchieved)
        {
            fptDataSet.add_species_count(--minValueAchieved);
            fptDataSet.add_first_passage_time(time);
        }
        while (speciesCounts[PROTEIN] > maxValueAchieved)
        {
            fptDataSet.add_species_count(++maxValueAchieved);
            fptDataSet.add_first_passage_time(time);
        }

        // If a species is past the limit, we are done.
        if (maxProtein > 0 && speciesCounts[PROTEIN] >= maxProtein)
        {
            break;
        }
        else if (minProtein > 0 && speciesCounts[PROTEIN] <= minProtein)
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
            fptDataSet.set_species(PROTEIN);
        }
    }

    Print::printf(Print::DEBUG, "Generated trajectory in %llu steps.", steps);

    // See if we finished the total simulation time.
    if (time == maxTime)
    {
        // We are recording every event, so add it.
        if (!writeTimeSteps)
        {
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(time);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
        }

        // Otherwise, write out the remaining time steps.
        else
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
    }

    // Otherwise we must have finished because of a species limit.
    else
    {
        // If we are writing time steps, just write out the last time.
        if (writeTimeSteps)
        {
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(time);
            for (uint i=0; i<numberSpecies; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
        }
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

double HillSwitch::KHillPropensity(double x, double kmin, double kmax, double x50, double h)
{
    return kmin+(((kmax-kmin)*(pow(x,h)))/((pow(x50,h))+(pow(x,h))));
}

}
}
