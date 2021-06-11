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

#include "config.h"
#include "core/Exceptions.h"
#include "core/Print.h"
#include "cme/CMESolver.h"
#include "FirstPassageTimes.pb.h"
#include "Lattice.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "rdme/Lattice.h"
#include "rdme/ByteLattice.h"
#include "rdme/NextSubvolumeSolver.h"
#include "reaction/ReactionQueue.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"

using lm::rdme::Lattice;
using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {

NextSubvolumeSolver::NextSubvolumeSolver():RDMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|RandomGenerator::UNIFORM)),numberSubvolumes(0),reactionPropensities(NULL),latticeSpacingSquared(0.0),reactionQueue(NULL),currentSubvolumeSpeciesCounts(NULL)
{
}

NextSubvolumeSolver::NextSubvolumeSolver(RandomGenerator::Distributions neededDists):RDMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|RandomGenerator::UNIFORM|neededDists)),numberSubvolumes(0),reactionPropensities(NULL),latticeSpacingSquared(0.0),reactionQueue(NULL)
{
}


NextSubvolumeSolver::~NextSubvolumeSolver()
{
}

void NextSubvolumeSolver::allocateModel(uint numberSpeciesA, uint numberReactionsA)
{
	RDMESolver::allocateModel(numberSpeciesA, numberReactionsA);

    // Allocate subvolume species counts.
	currentSubvolumeSpeciesCounts = new uint[numberSpeciesA];
}

void NextSubvolumeSolver::destroyModel()
{
	RDMESolver::destroyModel();
    if (currentSubvolumeSpeciesCounts != NULL) {delete[] currentSubvolumeSpeciesCounts; currentSubvolumeSpeciesCounts = NULL;}
}

void NextSubvolumeSolver::buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData)
{
    RDMESolver::buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData);

    // Allocate reaction queue.
    PROF_BEGIN(PROF_NSM_INIT_QUEUE);
    numberSubvolumes = lattice->getNumberSites();
    reactionQueue = new ReactionQueue(numberSubvolumes);
    PROF_END(PROF_NSM_INIT_QUEUE);

    // Update the propensity functions with the subvolume size.
    for (uint i=0; i<numberReactions; i++)
    {
        if (reactionTypes[i] == SecondOrderPropensityArgs::REACTION_TYPE)
        {
            ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->k *= numberSubvolumes;
        }
        else if (reactionTypes[i] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
        {
            ((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i])->k *= numberSubvolumes;
        }
    }
    latticeSpacingSquared = latticeSpacing*latticeSpacing;
}

void NextSubvolumeSolver::destroyDiffusionModel()
{
    RDMESolver::destroyDiffusionModel();

    // Free the reaction queue.
    if (reactionQueue != NULL) {delete reactionQueue; reactionQueue = NULL;}
    numberSubvolumes = 0;
    latticeSpacingSquared = 0.0;
}

void NextSubvolumeSolver::generateTrajectory()
{
    // Make sure we have propensity functions for every reaction.
    for (uint i=0; i<numberReactions; i++)
        if (propensityFunctions[i] == NULL || propensityFunctionArgs[i] == NULL)
            throw Exception("A reaction did not have a valid propensity function",i);

    // Initialize the species counts.
    for (uint i=0; i<numberSpecies; i++) speciesCounts[i] = initialSpeciesCounts[i];

    // Shadow the lattice member as a byte lattice.
    ByteLattice * lattice = (ByteLattice *)this->lattice;

    // Make sure that the initial species counts agree with the actual number in the lattice.
    checkSpeciesCountsAgainstLattice();

    // Get the interval for writing species counts and lattices.
    double speciesCountsWriteInterval=atof((*parameters)["writeInterval"].c_str());
    double nextSpeciesCountsWriteTime = speciesCountsWriteInterval;
    double latticeWriteInterval=atof((*parameters)["latticeWriteInterval"].c_str());
    double nextLatticeWriteTime = latticeWriteInterval;

    // Create the species counts data set to use during the simulation.
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpeciesToTrack);
    speciesCountsDataSet.set_number_entries(0);

    // Create the lattice data set to use during the simulation.
    lm::io::Lattice latticeDataSet;

    // Create the first passage time data set to use during the simulation.
    for (uint i=0; i<numberFptTrackedSpecies; i++)
    {
        fptTrackedSpecies[i].minValueAchieved = speciesCounts[fptTrackedSpecies[i].species];
        fptTrackedSpecies[i].maxValueAchieved = speciesCounts[fptTrackedSpecies[i].species];
        fptTrackedSpecies[i].dataSet.Clear();
        fptTrackedSpecies[i].dataSet.set_species(fptTrackedSpecies[i].species);
        fptTrackedSpecies[i].dataSet.add_species_count(speciesCounts[fptTrackedSpecies[i].species]);
        fptTrackedSpecies[i].dataSet.add_first_passage_time(0.0);
    }

    // Get the simulation time limit.
    si_time_t maxTime=atof((*parameters)["maxTime"].c_str());

    // Local cache of random numbers.
    double expRngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
    double uniRngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
    rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
    rng->getRandomDoubles(uniRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
    int expRngNext=0;
    int uniRngNext=0;

    // Initialize the reaction queue.
    double time = 0.0;
    expRngNext=updateAllSubvolumePropensities(time, expRngNext, expRngValues);

    // Record the initial species counts.
    recordSpeciesCounts(time, &speciesCountsDataSet);

    // Write the initial lattice.
    if (nextLatticeWriteTime > 0.0)
    	writeLatticeData(time, lattice, &latticeDataSet);

    // Run the next subvolume method.
    Print::printf(Print::DEBUG, "Running next subvolume simulation with %d species, %d reactions, %d subvolumes, %d site types for %e s. Writing species at %e and lattice at %e intervals", numberSpecies, numberReactions, numberSubvolumes, numberSiteTypes, maxTime, speciesCountsWriteInterval, latticeWriteInterval);
    PROF_BEGIN(PROF_SIM_EXECUTE);
    bool addedSpeciesCounts;
    bool addedFpt;
    unsigned long long steps=0;
    unsigned long long reactionSteps=0;
    unsigned long long maxSteps = atoll((*parameters)["maxSteps"].c_str());
    if (maxSteps == 0) maxSteps = ULLONG_MAX;
    bool affectedNeighbor;
    lattice_size_t subvolume;
    lattice_size_t neighborSubvolume;
    while (reactionSteps < maxSteps)
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Running step %d at time %f.", reactionSteps, time);
        addedSpeciesCounts = false;
        addedFpt = false;
        steps++;

        // Get the next subvolume with a reaction and the reaction time.
        subvolume = reactionQueue->getNextReaction();
        time = reactionQueue->getReactionEvent(subvolume).time;

        // If the new time is past the end time, we are done.
       if (time >= maxTime)
       {
           time = maxTime;
           break;
       }

       // Write species counts until the next write time is past the current time.
       while (nextSpeciesCountsWriteTime <= (time+EPS))
       {
           recordSpeciesCounts(nextSpeciesCountsWriteTime, &speciesCountsDataSet);
           nextSpeciesCountsWriteTime += speciesCountsWriteInterval;
           addedSpeciesCounts = true;
       }

       // Write lattice frames until the next write time is past the current time.
       while (nextLatticeWriteTime > 0.0 && nextLatticeWriteTime <= (time+EPS))
       {
           writeLatticeData(nextLatticeWriteTime, lattice, &latticeDataSet);
           nextLatticeWriteTime += latticeWriteInterval;
       }

       // Update the system with the reaction.
       affectedNeighbor = false;
       uniRngNext=performSubvolumeReaction(time, subvolume, uniRngNext, uniRngValues, &affectedNeighbor, &neighborSubvolume);

       // If the event didn't affect a neighbor, it must have been a reaction.
       if (!affectedNeighbor) reactionSteps++;

       // Update the propensity in the affected subvolumes.
       expRngNext=updateSubvolumePropensity(time, subvolume, expRngNext, expRngValues);
       if (affectedNeighbor) expRngNext=updateSubvolumePropensity(time, neighborSubvolume, expRngNext, expRngValues);

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
           writeSpeciesCounts(&speciesCountsDataSet);
       }
    }

    PROF_END(PROF_SIM_EXECUTE);
    Print::printf(Print::DEBUG, "Generated trajectory for replicate %d in %llu steps (%llu reaction events).", replicate, steps, reactionSteps);

    // If we finished the total time, write out the remaining time steps.
    if (time >= maxTime)
    {
        while (nextSpeciesCountsWriteTime <= (maxTime+EPS))
        {
            // Record the species counts.
            recordSpeciesCounts(nextSpeciesCountsWriteTime, &speciesCountsDataSet);
            nextSpeciesCountsWriteTime += speciesCountsWriteInterval;
        }
        while (nextLatticeWriteTime > 0.0 && nextLatticeWriteTime <= (maxTime+EPS))
        {
            writeLatticeData(nextLatticeWriteTime, lattice, &latticeDataSet);
            nextLatticeWriteTime += latticeWriteInterval;
        }
    }

    // Otherwise we must have finished because of a species limit or step, so just write out the last time.
    else
    {
        // Record the species counts.
        recordSpeciesCounts(time, &speciesCountsDataSet);
        if (nextLatticeWriteTime > 0.0)
        	writeLatticeData(time, lattice, &latticeDataSet);
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

    // Send any remaining species counts to the queue.
    writeSpeciesCounts(&speciesCountsDataSet);

    // Make sure that the final species counts agree with the actual number in the lattice.
    checkSpeciesCountsAgainstLattice();
}

void NextSubvolumeSolver::checkSpeciesCountsAgainstLattice()
{
	std::map<particle_t,uint> particleCounts = lattice->getParticleCounts();
	for (uint i=0; i<numberSpecies; i++)
	{
		if (speciesCounts[i] != ((particleCounts.count(i+1)>0)?particleCounts[i+1]:0))
			throw lm::Exception("Consistency error between species counts and lattice data", i, speciesCounts[i], ((particleCounts.count(i+1)>0)?particleCounts[i+1]:0));
	}
}

void NextSubvolumeSolver::writeLatticeData(double time, ByteLattice * lattice, lm::io::Lattice * latticeDataSet)
{
    Print::printf(Print::DEBUG, "Writing lattice at %e s", time);

    // Record the lattice data.
    latticeDataSet->Clear();
    latticeDataSet->set_lattice_x_size(lattice->getSize().x);
    latticeDataSet->set_lattice_y_size(lattice->getSize().y);
    latticeDataSet->set_lattice_z_size(lattice->getSize().z);
    latticeDataSet->set_particles_per_site(lattice->getMaxOccupancy());
    latticeDataSet->set_time(time);

    // Push it to the output queue.
    size_t payloadSize = lattice->getSize().x*lattice->getSize().y*lattice->getSize().z*lattice->getMaxOccupancy()*sizeof(uint8_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::BYTE_LATTICE, replicate, latticeDataSet, lattice, payloadSize, &lm::rdme::ByteLattice::nativeSerialize);
}

void NextSubvolumeSolver::recordSpeciesCounts(double time, lm::io::SpeciesCounts * speciesCountsDataSet)
{
    speciesCountsDataSet->set_number_entries(speciesCountsDataSet->number_entries()+1);
    speciesCountsDataSet->add_time(time);
    for (uint i=0; i<numberSpeciesToTrack; i++) speciesCountsDataSet->add_species_count(speciesCounts[i]);
}

void NextSubvolumeSolver::writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet)
{
    if (speciesCountsDataSet->number_entries() > 0)
    {
        PROF_BEGIN(PROF_SERIALIZE_COUNTS);
        // Push it to the output queue.
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, speciesCountsDataSet);

        // Reset the data set.
        speciesCountsDataSet->Clear();
        speciesCountsDataSet->set_number_species(numberSpeciesToTrack);
        speciesCountsDataSet->set_number_entries(0);
        PROF_END(PROF_SERIALIZE_COUNTS);
    }
}

int NextSubvolumeSolver::updateAllSubvolumePropensities(si_time_t time, int rngNext, double * expRngValues)
{
    // Update all of the subvolumes.
    PROF_BEGIN(PROF_NSM_BUILD_QUEUE);
    for (lattice_size_t s=0; s<numberSubvolumes; s++)
    {
        double propensity = calculateSubvolumePropensity(time, s);
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
        reactionQueue->updateReactionEvent(s, newTime, propensity);
    }
    PROF_END(PROF_NSM_BUILD_QUEUE);

    return rngNext;
}

int NextSubvolumeSolver::updateSubvolumePropensity(si_time_t time, lattice_size_t subvolume, int rngNext, double * expRngValues)
{
    double propensity = calculateSubvolumePropensity(time, subvolume);

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
    reactionQueue->updateReactionEvent(subvolume, newTime, propensity);

    return rngNext;
}

double NextSubvolumeSolver::calculateSubvolumePropensity(si_time_t time, lattice_size_t subvolume)
{
	// Update the species counts member for this subvolume.
	updateSpeciesCountsForSubvolume(subvolume);

	// Get the site type for this subvolume.
    site_t sourceSite = lattice->getSiteType(subvolume);

    // Calculate all of the reaction propensities.
	double subvolumePropensity = 0.0;
    for (uint i=0; i<numberReactions; i++)
    {
    	// Make sure the reaction can occur in this subvolume.
    	if (RL[i*numberSiteTypes+sourceSite])
    	{
			double (*propensityFunction)(double, uint *, void *) = (double (*)(double, uint*, void*))propensityFunctions[i];
			subvolumePropensity += (*propensityFunction)(time, currentSubvolumeSpeciesCounts, propensityFunctionArgs[i]);
    	}
    }

    // Calculate all of the diffusion propensities.
    const uint NUM_DEST_SITES=6;
    lattice_size_t neighboringSubvolumes[NUM_DEST_SITES];
    lattice->getNeighboringSites(subvolume, neighboringSubvolumes);
    for (uint i=0; i<numberSpecies; i++)
    {
    	if (currentSubvolumeSpeciesCounts[i] > 0)
    	{
    		for (uint j=0; j<NUM_DEST_SITES; j++)
    		{
    			subvolumePropensity += currentSubvolumeSpeciesCounts[i] * (DF[sourceSite*numberSiteTypes*numberSpecies + lattice->getSiteType(neighboringSubvolumes[j])*numberSpecies + i]/latticeSpacingSquared);
    		}
    	}
    }

    return subvolumePropensity;
}

void  NextSubvolumeSolver::updateSpeciesCountsForSubvolume(lattice_size_t subvolume)
{
	// Reset the species counts.
	for (uint i=0; i<numberSpecies; i++)
		currentSubvolumeSpeciesCounts[i] = 0;

	// Count the species that are in this subvolume.
	site_size_t numberParticles = lattice->getOccupancy(subvolume);
	for (site_size_t i=0; i<numberParticles; i++)
		currentSubvolumeSpeciesCounts[lattice->getParticle(subvolume, i)-1]++;
}

void NextSubvolumeSolver::updateSubvolumeWithSpeciesCounts(lattice_size_t subvolume)
{
	lattice->removeParticles(subvolume);
	for (uint i=0; i<numberSpecies; i++)
			addParticles(subvolume, i+1, currentSubvolumeSpeciesCounts[i]);
}

int NextSubvolumeSolver::performSubvolumeReaction(si_time_t time, lattice_size_t subvolume, int rngNext, double * uniRngValues, bool * affectedNeighbor, lattice_size_t * neighborSubvolume)
{
	// Only set the affected neighbor if this was a diffusion event.
    *affectedNeighbor = false;

	// Update the species counts member for this subvolume.
	updateSpeciesCountsForSubvolume(subvolume);

    // Stretch the random value across the total propensity range.
    if (rngNext >= TUNE_LOCAL_RNG_CACHE_SIZE)
    {
        rng->getRandomDoubles(uniRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
        rngNext=0;
    }
    double rngValue = uniRngValues[rngNext++]*reactionQueue->getReactionEvent(subvolume).propensity;

    // Get the site type for this subvolume.
    site_t sourceSite = lattice->getSiteType(subvolume);

    // See if it was a reaction that occurred.
    for (uint r=0; r<numberReactions; r++)
    {
    	// Make sure the reaction can occur in this subvolume.
    	if (RL[r*numberSiteTypes+sourceSite])
    	{
			double (*propensityFunction)(double, uint *, void *) = (double (*)(double, uint*, void*))propensityFunctions[r];
			double reactionPropensity = (*propensityFunction)(time, currentSubvolumeSpeciesCounts, propensityFunctionArgs[r]);
			if (reactionPropensity > 0.0)
			{
				if (rngValue <= reactionPropensity)
				{
					updateSpeciesCounts(r);
					updateCurrentSubvolumeSpeciesCounts(r);
					updateSubvolumeWithSpeciesCounts(subvolume);
					return rngNext;
				}
				else
				{
					rngValue -= reactionPropensity;
				}
			}
    	}
    }

    // See if it was a diffusion event that occurred.
    const uint NUM_DEST_SITES=6;
    lattice_size_t neighboringSubvolumes[NUM_DEST_SITES];
    lattice->getNeighboringSites(subvolume, neighboringSubvolumes);
    for (uint i=0; i<numberSpecies; i++)
    {
    	if (currentSubvolumeSpeciesCounts[i] > 0)
    	{
    		for (uint j=0; j<NUM_DEST_SITES; j++)
    		{
    			double diffusionPropensity = ((double)currentSubvolumeSpeciesCounts[i]) * (DF[sourceSite*numberSiteTypes*numberSpecies + lattice->getSiteType(neighboringSubvolumes[j])*numberSpecies + i]/latticeSpacingSquared);
    			if (rngValue <= diffusionPropensity)
    			{
    				currentSubvolumeSpeciesCounts[i]--;
    				updateSubvolumeWithSpeciesCounts(subvolume);
    				addParticles(neighboringSubvolumes[j], i+1, 1);
    				*affectedNeighbor = true;
    				*neighborSubvolume = neighboringSubvolumes[j];
    				return rngNext;
    			}
    			else
    			{
    	            rngValue -= diffusionPropensity;
    			}
    		}
    	}
    }

    throw Exception("Unable to determine correct reaction or diffusion event in subvolume.");
}

void NextSubvolumeSolver::updateCurrentSubvolumeSpeciesCounts(uint r)
{
    for (uint i=0; i<numberDependentSpecies[r]; i++)
    {
    	currentSubvolumeSpeciesCounts[dependentSpecies[r][i]] += dependentSpeciesChange[r][i];
    }
}

void NextSubvolumeSolver::addParticles(lattice_size_t subvolume, particle_t particle, uint count)
{
	for (uint i=0; i<count; i++)
	{
		try
		{
				lattice->addParticle(subvolume, particle);
		}
		catch (InvalidParticleException & e)
		{
			// We need to perform some overflow processing.
			const uint NUM_NEIGHBORS=6;
			lattice_size_t neighboringSubvolumes[NUM_NEIGHBORS];
			lattice->getNeighboringSites(subvolume, neighboringSubvolumes);
			bool handled = false;
			for (uint i=0; i<NUM_NEIGHBORS && !handled; i++)
			{
				if (lattice->getOccupancy(neighboringSubvolumes[i]) < lattice->getMaxOccupancy() && lattice->getSiteType(neighboringSubvolumes[i]) == lattice->getSiteType(subvolume))
				{
					lattice->addParticle(neighboringSubvolumes[i], particle);
					handled = true;
					Print::printf(Print::WARNING, "Handled overflow of particle type %d (%d total) from subvolume %d (type %d,occupancy %d) by moving to subvolume %d (type %d,occupancy %d).", particle, count, subvolume, lattice->getSiteType(subvolume), lattice->getOccupancy(subvolume), neighboringSubvolumes[i], lattice->getSiteType(neighboringSubvolumes[i]), lattice->getOccupancy(neighboringSubvolumes[i]));
				}
			}
			if (!handled) throw Exception("Unable to handle overflow at site", subvolume);
		}
	}
}

}
}
