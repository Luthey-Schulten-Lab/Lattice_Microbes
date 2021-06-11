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
#include "config.h"
#if defined(MACOSX)
#elif defined(LINUX)
#include <time.h>
#endif
#include "core/Math.h"
#include "core/Print.h"
#include "cme/CMESolver.h"
#include "cme/NextReactionSolver.h"
#include "FirstPassageTimes.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "core/ResourceAllocator.h"
#include "reaction/ReactionQueue.h"
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
using lm::reaction::ReactionQueue;
using lm::rng::RandomGenerator;

namespace lm {
namespace cme {

NextReactionSolver::NextReactionSolver():CMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL)),reactionQueue(NULL)
{
}

NextReactionSolver::NextReactionSolver(RandomGenerator::Distributions neededDists):CMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|neededDists)),reactionQueue(NULL)
{
}

NextReactionSolver::~NextReactionSolver()
{
}

void NextReactionSolver::buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionType, const double * K, const int * SA, const uint * DA, const uint kCols)
{
    CMESolver::buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionType, K, SA, DA, kCols);

    // Allocate reaction queue.
    reactionQueue = new ReactionQueue(numberReactions);
}

void NextReactionSolver::destroyModel()
{
    CMESolver::destroyModel();

    // Free the reaction queue.
    if (reactionQueue != NULL) {delete reactionQueue; reactionQueue = NULL;}
}

void NextReactionSolver::generateTrajectory()
{
    // Make sure we have propensity functions for every reaction.
    for (uint i=0; i<numberReactions; i++)
        if (propensityFunctions[i] == NULL || propensityFunctionArgs[i] == NULL)
            throw Exception("A reaction did not have a valid propensity function",i);

    // Create local copies of the data for efficiency.
    uint numberSpecies = this->numberSpecies;
    uint numberReactions = this->numberReactions;
    uint * speciesCounts = this->speciesCounts;

    // Initialize the species counts.
    for (uint i=0; i<numberSpecies; i++) speciesCounts[i] = initialSpeciesCounts[i];

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
    double nextParameterWriteTime = INFINITY; //std::numeric_limits<double>::infinity();
    double parameterWriteInterval = atof((*parameters)["parameterWriteInterval"].c_str());
    if (trackedParameters.size() > 0 && parameterWriteInterval > 0.0)
        nextParameterWriteTime = recordParameters(0.0, parameterWriteInterval, 0.0);

    // Get the simulation time limit.
    double maxTime=atof((*parameters)["maxTime"].c_str());

    // Local cache of random numbers.
    double expRngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
    rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
    int rngNext=0;

    // Initialize the reaction queue.
    double time = 0.0;
    rngNext=updateAllReactionEvents(time, rngNext, expRngValues);

    // Run the next reaction method.
    Print::printf(Print::DEBUG, "Running next reaction simulation with %d species, %d reactions, %d species limits, and write mode %d", numberSpecies, numberReactions, numberSpeciesLimits, writeTimeSteps);
    PROF_BEGIN(PROF_SIM_EXECUTE);
    bool addedSpeciesCounts;
    bool addedFpt;
    unsigned long long steps=0;
    unsigned long long maxSteps = atoll((*parameters)["maxSteps"].c_str());
    if (maxSteps == 0) maxSteps = ULLONG_MAX;

	// Call Hook Simulation with initla time
    hookSimulation(time);

    while (steps < maxSteps && !reachedSpeciesLimit())
    {
        addedSpeciesCounts = false;
        addedFpt = false;
        steps++;

        //if (steps%100000 == 0)
        //{
			//Print::printf(Print::INFO, "Before step %3d reaction times: %12.6e, %12.6e, %12.6e, %12.6e, %12.6e, %12.6e",steps,reactionQueue->getReactionEvent(0).time,reactionQueue->getReactionEvent(1).time,reactionQueue->getReactionEvent(2).time,reactionQueue->getReactionEvent(3).time,reactionQueue->getReactionEvent(4).time,reactionQueue->getReactionEvent(5).time);
			//Print::printf(Print::DEBUG, "                         props: %12.6e, %12.6e, %12.6e, %12.6e, %12.6e, %12.6e",reactionQueue->getReactionEvent(0).propensity,reactionQueue->getReactionEvent(1).propensity,reactionQueue->getReactionEvent(2).propensity,reactionQueue->getReactionEvent(3).propensity,reactionQueue->getReactionEvent(4).propensity,reactionQueue->getReactionEvent(5).propensity);
			//Print::printf(Print::DEBUG, "                        counts: %12d, %12d, %12d, %12d",speciesCounts[0],speciesCounts[1],speciesCounts[2],speciesCounts[3]);
        //}

        // Get the next reaction.
        uint r = reactionQueue->getNextReaction();
        time = reactionQueue->getReactionEvent(r).time;

        //Print::printf(Print::DEBUG, "Step %d was reaction %d at time %0.6e",steps,r,time);

         // If the new time is past the end time, we are done.
        if (time >= maxTime)
        {
            time = maxTime;
            break;
        }

        // If we are writing time steps, write out any time steps before this event occurred.
        bool hookChangedState = false;
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

				// Hook the simulation so that the user can handle this
				if(hookSimulation(time)) {
                    hookChangedState = true;
                }
            }
        }

        // If we are recording parameter values, write out the values before this event occurred.
        if (nextParameterWriteTime <= (time+1e-9))
        {
            nextParameterWriteTime = recordParameters(nextParameterWriteTime, parameterWriteInterval, time);
        }

        // Update species counts and propensities given the reaction that occurred.
        updateSpeciesCounts(r);
        if(hookChangedState)
            rngNext=updateAllReactionEvents(time, rngNext, expRngValues);
        else
            rngNext=updateReactionEvents(r, time, rngNext, expRngValues);

        // See if we need to update our rng caches.
        if (rngNext >= TUNE_LOCAL_RNG_CACHE_SIZE)
        {
            rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            rngNext=0;
        }

        // Generate a new time for the reaction that occurred.
        double propensity = reactionQueue->getReactionEvent(r).propensity;
        double newTime = time+(expRngValues[rngNext++]/propensity);
        reactionQueue->updateReactionEvent(r, newTime, propensity);

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
    }
    PROF_END(PROF_SIM_EXECUTE);
    Print::printf(Print::DEBUG, "Generated trajectory for replicate %d in %llu steps.", replicate, steps);

    // If we finished the total time, write out the remaining time steps.
    if (time >= maxTime)
    {
    	if (writeTimeSteps)
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
    	}
    	else
    	{
			// We are tracking events, so just write out the last time.
			speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
			speciesCountsDataSet.add_time(maxTime);
			for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);

    	}

		// If we are recording parameter values, write out the remaining value intervals.
		if (nextParameterWriteTime <= (maxTime+1e-9))
		{
			recordParameters(nextParameterWriteTime, parameterWriteInterval, maxTime);
		}
    }

    // Otherwise we must have finished because of a species limit or step, so just write out the last time and last parameter value.
    else
    {
        // Record the species counts.
        speciesCountsDataSet.set_number_entries(speciesCountsDataSet.number_entries()+1);
        speciesCountsDataSet.add_time(time);
        for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet.add_species_count(speciesCounts[i]);

        // Record the final parameter value.
		recordParameters(time, parameterWriteInterval, time);
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

int NextReactionSolver::updateAllReactionEvents(double time, int rngNext, double * expRngValues)
{
    // Update all of the reactions.
    for (uint r=0; r<numberReactions; r++)
    {
        double (*propensityFunction)(double, uint * speciesCounts, void * args) = (double (*)(double, uint*, void*))propensityFunctions[r];
        double propensity = (*propensityFunction)(time, speciesCounts, propensityFunctionArgs[r]);
        double newTime = INFINITY;
        if (propensity > 0.0)
        {
            if (rngNext >= TUNE_LOCAL_RNG_CACHE_SIZE)
            {
                rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                rngNext=0;
            }
            newTime = time+expRngValues[rngNext++]/propensity;
        }
        reactionQueue->updateReactionEvent(r, newTime, propensity);
    }
    return rngNext;
}

int NextReactionSolver::updateReactionEvents(uint sourceReaction, double time, int rngNext, double * expRngValues)
{
    // Update the propensities of the dependent reactions.
    for (uint i=0; i<numberDependentReactions[sourceReaction]; i++)
    {
        uint r = dependentReactions[sourceReaction][i];
        double (*propensityFunction)(double, uint * speciesCounts, void * args) = (double (*)(double, uint*, void*))propensityFunctions[r];
        double newPropensity = (*propensityFunction)(time, speciesCounts, propensityFunctionArgs[r]);
        double newTime = INFINITY;

        // If there is some propensity and this is not the reaction that occurred, figure out the new time.
        if (newPropensity > 0.0 && r != sourceReaction)
        {
            #if TUNE_NRM_REUSE_PROPENSITIES == 1
            double oldTime = reactionQueue->getReactionEvent(r).time;
            double oldPropensity = reactionQueue->getReactionEvent(r).propensity;

            // If this reaction had a previous time, reuse it.
            if (oldPropensity > 0.0)
            {
                newTime = time+((oldPropensity/newPropensity)*(oldTime-time));
            }

            // Otherwise choose a new time.
            else
            {
            #endif

                if (rngNext >= TUNE_LOCAL_RNG_CACHE_SIZE)
                {
                    rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                    rngNext=0;
                }
                newTime = time+expRngValues[rngNext++]/newPropensity;
            #if TUNE_NRM_REUSE_PROPENSITIES == 1
            }
            #endif
        }
        reactionQueue->updateReactionEvent(r, newTime, newPropensity);
    }

    return rngNext;
}

int NextReactionSolver::hookSimulation(double time)
{
    // Overload this function in derivative classes
    // Return 0 if the lattice state is unchanged
    // Return 1 if the lattice state has been modified, 
    //          and it needs to be copied back to the GPU.
    return 0;
}




}
}
