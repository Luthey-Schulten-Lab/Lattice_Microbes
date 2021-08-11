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

#include <string>
#include <list>
#include <map>
#include <cmath>
#include <iostream>
#include <climits>
#include "config.h"
#if defined(MACOSX)
#elif defined(LINUX)
#include <time.h>
#endif
#include "core/Math.h"
#include "core/Print.h"
#include "cme/CMESolver.h"
#include "cme/GillespieDSolver.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "core/ResourceAllocator.h"
#include "rng/RandomGenerator.h"
#include "rng/XORShift.h"
#ifdef OPT_CUDA
#include "rng/XORWow.h"
#endif
#include "thread/Thread.h"
#include "thread/Worker.h"
#include "lptf/Profile.h"




using std::string;
using std::list;
using std::map;
using lm::rng::RandomGenerator;

namespace lm {
namespace cme {

GillespieDSolver::GillespieDSolver():CMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|RandomGenerator::UNIFORM)),propensities(NULL)
{
}

GillespieDSolver::~GillespieDSolver()
{
}

void GillespieDSolver::buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionType, const double * K, const int * SA, const uint * DA, const uint kCols)
{
    CMESolver::buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionType, K, SA, DA, kCols);

    // Allocate reaction propensities table.
    propensities = new double[numberReactions];
    for (uint i=0; i<numberReactions; i++)
    {
        propensities[i] = 0.0;
    }
}

void GillespieDSolver::destroyModel()
{
    CMESolver::destroyModel();

    // Free the propensities.
    if (propensities != NULL) {delete[] propensities; propensities = NULL;}
}

void GillespieDSolver::generateTrajectory()
{
    // Make sure we have propensity functions for every reaction.
    for (uint i=0; i<numberReactions; i++)
        if (propensityFunctions[i] == NULL || propensityFunctionArgs[i] == NULL)
            throw Exception("A reaction did not have a valid propensity function",i);

    // Create local copies of the data for efficiency.
    uint numberSpecies = this->numberSpecies;
    uint numberReactions = this->numberReactions;
    uint * speciesCounts = this->speciesCounts;
    double * propensities = this->propensities;

    // Initialize the species counts.
    for (uint i=0; i<numberSpecies; i++) speciesCounts[i] = initialSpeciesCounts[i];

    // Initialize the propensities.
    updateAllPropensities(0.0);
    double totalPropensity = 0.0;
    for (uint i=0; i<numberReactions; i++) totalPropensity += propensities[i];

    // Create the species counts data set to track during the simulation.
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpeciesToTrack);
    speciesCountsDataSet.set_number_entries(1);
    speciesCountsDataSet.add_time(0.0);
    for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);

    // Initialize tracking of the first passage times.
    for (uint i=0; i<numberFptTrackedSpecies; i++)
    {
        fptTrackedSpecies[i].minValueAchieved = speciesCounts[fptTrackedSpecies[i].species];
        fptTrackedSpecies[i].maxValueAchieved = speciesCounts[fptTrackedSpecies[i].species];
        fptTrackedSpecies[i].dataSet.Clear();
        fptTrackedSpecies[i].dataSet.set_species(fptTrackedSpecies[i].species);
        fptTrackedSpecies[i].dataSet.add_species_count(speciesCounts[fptTrackedSpecies[i].species]);
        fptTrackedSpecies[i].dataSet.add_first_passage_time(0.0);
    }

    // Get the interval for writing species counts.
    double writeInterval=atof((*parameters)["writeInterval"].c_str());
    bool writeTimeSteps = (writeInterval > 0.0);
    double nextSpeciesCountsWriteTime = writeInterval;

    // Get the interval for writing parameters.
    double nextParameterWriteTime = INFINITY;
    double parameterWriteInterval = atof((*parameters)["parameterWriteInterval"].c_str());
    if (trackedParameters.size() > 0 && parameterWriteInterval > 0.0)
        nextParameterWriteTime = recordParameters(0.0, parameterWriteInterval, 0.0);

    // Get the simulation time limit.
    double maxTime=atof((*parameters)["maxTime"].c_str());

    // Local cache of random numbers.
    double rngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
    double expRngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
    rng->getRandomDoubles(rngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
    rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
    int rngNext=0;

    // Run the direct method.
    Print::printf(Print::DEBUG, "Running Gillespie direct simulation with %d species, %d reactions, %d species limits, and write mode %d", numberSpecies, numberReactions, numberSpeciesLimits, writeTimeSteps);
    PROF_BEGIN(PROF_SIM_EXECUTE);
    bool addedSpeciesCounts;
    bool addedFpt;
    double time = 0.0;
    unsigned long long steps=0;
    unsigned long long maxSteps = atoll((*parameters)["maxSteps"].c_str());
    if (maxSteps == 0) maxSteps = ULLONG_MAX;

	bool hookEnabled=false;
	double nextHookTime=0.0;
	double hookInterval=0.0;

	if((*parameters)["hookInterval"] != "")
	{
		hookEnabled=true;
		hookInterval=atof((*parameters)["hookInterval"].c_str());
		nextHookTime=hookInterval;
	}

	onBeginTrajectory();

    while (totalPropensity > 0 && steps < maxSteps && !reachedSpeciesLimit())
    {
        addedSpeciesCounts = false;
        addedFpt = false;
        steps++;

        // See if we need to update our rng caches.
        if (rngNext >= TUNE_LOCAL_RNG_CACHE_SIZE)
        {
            rng->getRandomDoubles(rngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            rngNext=0;
        }

        // Calculate the time to the next reaction.
        double expR = expRngValues[rngNext];
        time += expR/totalPropensity;

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
            while (nextSpeciesCountsWriteTime <= (time+1e-9))
            {
                // Record the species counts.
                speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
                speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
                for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
                nextSpeciesCountsWriteTime += writeInterval;
                addedSpeciesCounts = true;
            }
        }

        // If we are recording parameter values, write out the values before this event occurred.
        if (nextParameterWriteTime <= (time+1e-9))
        {
            nextParameterWriteTime = recordParameters(nextParameterWriteTime, parameterWriteInterval, time);
        }

        // Calculate which reaction it was.
        double rngValue = rngValues[rngNext]*totalPropensity;
        uint r=0;
        for (; r<(numberReactions-1); r++)
        {
            if (rngValue < propensities[r])
                break;
            else
                rngValue -= propensities[r];
        }

        // Update species counts and propensities given the reaction that occurred.
        updateSpeciesCounts(r);
        updatePropensities(time, r);

		if(hookEnabled && nextHookTime <= (time+1e-9))
		{
				// Hook the simulation so that the user can handle this
				if(hookSimulation(time)) {
        		  updateAllPropensities(time);
        		}
                nextHookTime += hookInterval;
		}

        // Recalculate the total propensity.
        totalPropensity = 0.0;
        for (uint i=0; i<numberReactions; i++) {
            totalPropensity += propensities[i];
        }

        // If we are recording every event, add it.
        if (!writeTimeSteps)
        {
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(time);
            for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
        }

        // Update the first passage time tables.
        for (uint i=0; i<numberFptTrackedSpecies; i++)
        {
            uint speciesCount = speciesCounts[fptTrackedSpecies[i].species];
            while (fptTrackedSpecies[i].minValueAchieved > speciesCount)
            {
                fptTrackedSpecies[i].dataSet.add_species_count(--fptTrackedSpecies[i].minValueAchieved);
                fptTrackedSpecies[i].dataSet.add_first_passage_time(time);
                addedFpt = true;
            }
            while (fptTrackedSpecies[i].maxValueAchieved < speciesCount)
            {
                fptTrackedSpecies[i].dataSet.add_species_count(++fptTrackedSpecies[i].maxValueAchieved);
                fptTrackedSpecies[i].dataSet.add_first_passage_time(time);
                addedFpt = true;
            }

            // See if we have accumulated enough fpt data to send.
            if (addedFpt && fptTrackedSpecies[i].dataSet.first_passage_time_size() >= TUNE_FIRST_PASSAGE_TIME_BUFFER_SIZE)
            {
                // Push it to the output queue.
                PROF_BEGIN(PROF_SERIALIZE_FPT);
                lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptTrackedSpecies[i].dataSet);
                PROF_END(PROF_SERIALIZE_FPT);

                // Reset the data set.
                fptTrackedSpecies[i].dataSet.Clear();
                fptTrackedSpecies[i].dataSet.set_species(fptTrackedSpecies[i].species);
            }
        }

        // See if we have accumulated enough species counts to send.
        if (addedSpeciesCounts && speciesCountsDataSet.number_entries() >= TUNE_SPECIES_COUNTS_BUFFER_SIZE)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_COUNTS);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
            PROF_END(PROF_SERIALIZE_COUNTS);

            // Reset the data set.
            speciesCountsDataSet.Clear();
            speciesCountsDataSet.set_number_species(numberSpeciesToTrack);
            speciesCountsDataSet.set_number_entries(0);
        }

         // Go to the next rng pair.
        rngNext++;

    }
    PROF_END(PROF_SIM_EXECUTE);

	onEndTrajectory();

    Print::printf(Print::DEBUG, "Generated trajectory for replicate %d in %llu steps.", replicate, steps);

    // If we finished the total time or ran out of reactions, write out the remaining time steps.
    if (time >= maxTime || totalPropensity <= 0)
    {
        Print::printf(Print::DEBUG, "Finished with time %e (%e)", time, maxTime);
        while (nextSpeciesCountsWriteTime <= (maxTime+1e-9))
        {
            Print::printf(Print::VERBOSE_DEBUG, "Recording event at time %e (%e)", nextSpeciesCountsWriteTime, maxTime);
            // Record the species counts.
            speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
            speciesCountsDataSet.add_time(nextSpeciesCountsWriteTime);
            for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
            nextSpeciesCountsWriteTime += writeInterval;
        }
        Print::printf(Print::VERBOSE_DEBUG, "Done recording events at time %e (%e)", nextSpeciesCountsWriteTime, maxTime);

        // If we are recording parameter values, write out the remaining value intervals.
        if (nextParameterWriteTime <= (maxTime+1e-9))
        {
            recordParameters(nextParameterWriteTime, parameterWriteInterval, maxTime);
        }
    }

    // Otherwise we must have finished because of a species limit or step, so just write out the last time.
    else
    {
        // Record the species counts.
        speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
        speciesCountsDataSet.add_time(time);
        for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);
    }

    // Send any remaining species counts to the queue.
    if (speciesCountsDataSet.number_entries() > 0)
    {
        PROF_BEGIN(PROF_SERIALIZE_COUNTS);
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, &speciesCountsDataSet);
        PROF_END(PROF_SERIALIZE_COUNTS);
    }

    // Send any remaining first passage times to the queue.
    for (uint i=0; i<numberFptTrackedSpecies; i++)
    {
        if (fptTrackedSpecies[i].dataSet.first_passage_time_size() > 0)
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_FPT);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::FIRST_PASSAGE_TIMES, replicate, &fptTrackedSpecies[i].dataSet);
            PROF_END(PROF_SERIALIZE_FPT);
        }
    }

    // Send any remaining parameter values to the queue.
    queueRecordedParameters(true);
}

void GillespieDSolver::updateAllPropensities(double time)
{
    // Update the propensities.
    for (uint i=0; i<numberReactions; i++)
    {
        double (*propensityFunction)(double, uint * speciesCounts, void * args) = (double (*)(double, uint*, void*))propensityFunctions[i];
        propensities[i] = (*propensityFunction)(time, speciesCounts, propensityFunctionArgs[i]);
    }
}

void GillespieDSolver::updatePropensities(double time, uint sourceReaction)
{
    // Update the propensities of the dependent reactions.
    for (uint i=0; i<numberDependentReactions[sourceReaction]; i++)
    {
        uint r = dependentReactions[sourceReaction][i];
        double (*propensityFunction)(double, uint * speciesCounts, void * args) = (double (*)(double, uint*, void*))propensityFunctions[r];
        propensities[r] = (*propensityFunction)(time, speciesCounts, propensityFunctionArgs[r]);
    }
}


int GillespieDSolver::hookSimulation(double time)
{
    // Overload this function in derivative classes
    // Return 0 if the lattice state is unchanged
    // Return 1 if the lattice state has been modified, 
    //          and it needs to be copied back to the GPU.
    return 0;
}

int GillespieDSolver::onBeginTrajectory()
{
	return 0;
}

int GillespieDSolver::onEndTrajectory()
{
	return 0;
}
	

}
}
