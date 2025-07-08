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
#include "FirstPassageTimes.pb.h"
#include "ReactionModel.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "core/ResourceAllocator.h"
#include "rng/RandomGenerator.h"
#include "rng/XORShift.h"
#ifdef OPT_CUDA
#include "rng/XORWow.h"
#include "cuda/lm_cuda.h"
#endif
#include "thread/Thread.h"
#include "thread/Worker.h"
#include "lptf/Profile.h"




using std::string;
using std::list;
using std::map;

namespace lm {
namespace cme {

CMESolver::CMESolver(RandomGenerator::Distributions neededDists)
:neededDists(neededDists),replicate(-1),parameters(NULL),resources(NULL),rng(NULL),numberSpecies(0),numberSpeciesToTrack(0),numberReactions(0),initialSpeciesCounts(NULL),speciesCounts(NULL),reactionTypes(NULL),S(NULL),D(NULL),propensityFunctions(NULL),propensityFunctionArgs(NULL),numberSpeciesLimits(0),speciesLimits(NULL),numberFptTrackedSpecies(0),fptTrackedSpecies(NULL),numberDependentSpecies(NULL),dependentSpecies(NULL),dependentSpeciesChange(NULL),numberDependentReactions(NULL),dependentReactions(NULL)
{
}

CMESolver::~CMESolver()
{
    // Free any model memory.
    destroyModel();

    // Free any other memory.
    if (rng != NULL) {delete rng; rng = NULL;}
    if (speciesLimits != NULL) {delete[] speciesLimits; speciesLimits = NULL;}
    if (fptTrackedSpecies != NULL) {delete[] fptTrackedSpecies; fptTrackedSpecies = NULL;}
}



void CMESolver::initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources)
{
    this->replicate = replicate;
    this->parameters = parameters;
    this->resources = resources;

    if (neededDists != RandomGenerator::NONE)
    {
        #ifdef OPT_CUDA
        // Create the cuda based rng, if we are using cuda and we have a cuda device assigned.
        if (resources->cudaDevices.size() > 0)
        {
			try
			{
				rng = new lm::rng::XORWow(resources->cudaDevices[0], replicate, atoi((*parameters)["seed"].c_str()), neededDists);
				Print::printf(Print::INFO, "Seeding xorwow rng with top word %u and bottom word %u", (unsigned int)(rng->getSeed()>>32), (unsigned int)(rng->getSeed()&0xFFFFFFFFLL));
			}
			catch(lm::CUDAException e)
			{
				Print::printf(Print::ERROR, "Initializing GPU random number generator failed: %s.  Falling back to CPU random number generator.", e.what());
				rng=NULL;
			}
		
        }
        #endif

        if (rng == NULL)
        {
            rng = new lm::rng::XORShift(replicate, atoi((*parameters)["seed"].c_str()));
            Print::printf(Print::INFO, "Seeding xorshift rng with top word %u and bottom word %u", (unsigned int)(rng->getSeed()>>32), (unsigned int)(rng->getSeed()&0xFFFFFFFFLL));
        }
    }

    // Set the fpt tracked species from the parameters.
    string listString = (*parameters)["fptTrackingList"];
    list<uint> fptList;
    size_t start=0, end=0;
    while (end != string::npos)
    {
        end = listString.find(',', start);
        string trackedSpecies = listString.substr(start, (end == string::npos) ? string::npos : end - start);
        fptList.push_back(atoi(trackedSpecies.c_str()));
        start = end+1;
    }
    setFptTrackingList(fptList);

    // Set the species lower limits from the parameters.
    listString = (*parameters)["speciesLowerLimitList"];
    start=0, end=0;
    while (end != string::npos)
    {
        end = listString.find(',', start);
        string speciesLowerLimit = listString.substr(start, (end == string::npos) ? string::npos : end - start);

        size_t equalsPos=0;
        equalsPos = speciesLowerLimit.find(':', 0);
        if (equalsPos > 0 && equalsPos < speciesLowerLimit.length()-1)
        {
        	uint parsedSpecies = atoi(speciesLowerLimit.substr(0, equalsPos).c_str());
        	uint parsedLimit = atoi(speciesLowerLimit.substr(equalsPos+1, string::npos).c_str());
        	setSpeciesLowerLimit(parsedSpecies, parsedLimit);
        	Print::printf(Print::DEBUG, "Parsed lower limit %s to: %d => %d", speciesLowerLimit.c_str(), parsedSpecies, parsedLimit);
        }
        start = end+1;
    }

    // Set the species upper limits from the parameters.
    listString = (*parameters)["speciesUpperLimitList"];
    start=0, end=0;
    while (end != string::npos)
    {
        end = listString.find(',', start);
        string speciesUpperLimit = listString.substr(start, (end == string::npos) ? string::npos : end - start);

        size_t equalsPos=0;
        equalsPos = speciesUpperLimit.find(':', 0);
        if (equalsPos > 0 && equalsPos < speciesUpperLimit.length()-1)
        {
        	uint parsedSpecies = atoi(speciesUpperLimit.substr(0, equalsPos).c_str());
        	uint parsedLimit = atoi(speciesUpperLimit.substr(equalsPos+1, string::npos).c_str());
        	setSpeciesUpperLimit(parsedSpecies, parsedLimit);
        	Print::printf(Print::DEBUG, "Parsed upper limit %s to: %d <= %d", speciesUpperLimit.c_str(), parsedSpecies, parsedLimit);
        }
        start = end+1;
    }
}

void CMESolver::allocateModel(uint numberSpeciesA, uint numberReactionsA)
{
    // Set the number of species and reactions.
    numberSpecies = numberSpeciesA;
    numberSpeciesToTrack = numberSpeciesA;
    numberReactions = numberReactionsA;

    // Allocate species counts.
    initialSpeciesCounts = new uint[numberSpecies];
    speciesCounts = new uint[numberSpecies];
    for (uint i=0; i<numberSpecies; i++)
    {
        initialSpeciesCounts[i] = 0;
        speciesCounts[i] = 0;
    }

    if (numberReactions > 0)
    {
		// Allocate reaction/species matrices.
		reactionTypes = new uint [numberReactions];
		for (uint i=0; i<numberReactions; i++) reactionTypes[i] = 0;
		S = new int [numberSpecies*numberReactions];
		D = new uint [numberSpecies*numberReactions];
		for (uint i=0; i<numberSpecies*numberReactions; i++)
		{
			S[i] = 0;
			D[i] = 0;
		}

		// Allocate propensity function tables.
		propensityFunctions = new void *[numberReactions];
		propensityFunctionArgs = new void *[numberReactions];
		for (uint i=0; i<numberReactions; i++)
		{
			propensityFunctions[i] = NULL;
			propensityFunctionArgs[i] = NULL;
		}

		// Allocate the species dependency tables.
		numberDependentSpecies = new uint[numberReactions];
		dependentSpecies = new uint*[numberReactions];
		dependentSpeciesChange = new int*[numberReactions];

		// Allocate the reaction dependency tables.
		numberDependentReactions = new uint[numberReactions];
		dependentReactions = new uint*[numberReactions];
    }
}

void CMESolver::destroyModel()
{
    if (initialSpeciesCounts != NULL) {delete[] initialSpeciesCounts; initialSpeciesCounts = NULL;}
    if (speciesCounts != NULL) {delete[] speciesCounts; speciesCounts = NULL;}
    if (reactionTypes != NULL) {delete[] reactionTypes; reactionTypes = NULL;}
    if (S != NULL) {delete[] S; S = NULL;}
    if (D != NULL) {delete[] D; D = NULL;}
    if (propensityFunctions != NULL) {delete[] propensityFunctions; propensityFunctions = NULL;}

    // Free the propensity function arguments.
    if (propensityFunctionArgs != NULL) {delete[] propensityFunctionArgs; propensityFunctionArgs = NULL;}
    for (list<PropensityArgs *>::iterator it = propensityArgs.begin(); it != propensityArgs.end(); it++) delete *it;
    propensityArgs.clear();

    // Free the species dependency tables.
    if (numberDependentSpecies != NULL) {delete[] numberDependentSpecies; numberDependentSpecies = NULL;}
    if (dependentSpecies != NULL)
    {
        for (uint i=0; i<numberReactions; i++)
            delete[] dependentSpecies[i];
        delete[] dependentSpecies;
        dependentSpecies = NULL;
    }
    if (dependentSpeciesChange != NULL)
    {
        for (uint i=0; i<numberReactions; i++)
            delete[] dependentSpeciesChange[i];
        delete[] dependentSpeciesChange;
        dependentSpeciesChange = NULL;
    }

    // Free the reaction dependency tables.
    if (numberDependentReactions != NULL) {delete[] numberDependentReactions; numberDependentReactions = NULL;}
    if (dependentReactions != NULL)
    {
        for (uint i=0; i<numberReactions; i++)
            delete[] dependentReactions[i];
        delete[] dependentReactions;
        dependentReactions = NULL;
    }

    // Reset the species and reaction counts.
    numberSpecies = 0;
    numberSpeciesToTrack = 0;
    numberReactions = 0;
}

void CMESolver::setReactionModel(lm::io::ReactionModel * rm)
{
    if (rm->number_reactions() != (uint)rm->reaction_size()) throw InvalidArgException("rm", "number of reaction does not agree with reaction list size");

    // Figure out the max number of columns we need in the k matrix.
    uint kCols = 0;
    for (uint i=0; i<rm->number_reactions(); i++)
        kCols = max(kCols,(uint)rm->reaction(i).rate_constant_size());

    // Set the K and reaction type tables.
    uint * reactionType = new uint[rm->number_reactions()];
    double * K = new double[rm->number_reactions()*kCols];
    for (uint i=0; i<rm->number_reactions(); i++)
    {
        reactionType[i] = rm->reaction(i).type();
        for (uint j=0; j<(uint)rm->reaction(i).rate_constant_size(); j++)
        {
            K[i*kCols+j] = rm->reaction(i).rate_constant(j);
        }
    }

    // Build the model.
    buildModel(rm->number_species(), rm->number_reactions(), rm->initial_species_count().data(), reactionType, K, rm->stoichiometric_matrix().data(), rm->dependency_matrix().data(), kCols);

    // Free any resources.
    if (reactionType !=  NULL) {delete [] reactionType; reactionType = NULL;}
    if (K !=  NULL) {delete [] K; K = NULL;}
}

void CMESolver::buildModel(const uint numberSpeciesA,
                           const uint numberReactionsA,
                           const uint * initialSpeciesCountsA,
                           const uint * reactionTypesA,
                           const double * K,
                           const int * SA,
                           const uint * DA,
                           const uint kCols)
{
    if (numberReactionsA > 0 && kCols == 0) throw InvalidArgException("K", "must have at least 1 column");

    // Destroy the previous model, if we have one.
    destroyModel();

    // Allocate space for the new model.
    allocateModel(numberSpeciesA, numberReactionsA);

    // Set the initial species counts.
    for (uint i=0; i<numberSpecies; i++)
    {
        initialSpeciesCounts[i] = initialSpeciesCountsA[i];
    }

    // Set the reaction types.
    for (uint i=0; i<numberReactions; i++)
    {
    	reactionTypes[i] = reactionTypesA[i];
    }

    // Set the stoichiometric and dependency matrices.
    for (uint i=0; i<numberSpecies*numberReactions; i++)
    {
        S[i] = SA[i];
        D[i] = DA[i];
    }

    // Create the propensity functions table.
    for (uint i=0; i<numberReactions; i++)
    {
        switch(reactionTypes[i]) {
			// Zeroth order reaction
            case ZerothOrderPropensityArgs::REACTION_TYPE:
            {
                // Find the dependencies.
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) 
                    if (D[j*numberReactions+i] == 1) 
                        numberDependencies++;

                if (numberDependencies > 0) {
                    throw InvalidArgException("D", "zeroth order reaction cannot have any dependencies",numberDependencies);
                } else {
                    propensityFunctions[i] = (void *)&zerothOrderPropensity;
                    propensityFunctionArgs[i] =  (void *)new ZerothOrderPropensityArgs(K[i*kCols]);
                    propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                }
            }
            break;

            // First order reaction
            case FirstOrderPropensityArgs::REACTION_TYPE:
            {
                // Find the dependencies.
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1) {
                        numberDependencies++;

                        // Set the table entry to the first non-zero dependency.
                        if (numberDependencies == 1) {
                            propensityFunctions[i] = (void *)&firstOrderPropensity;
                            propensityFunctionArgs[i] =  (void *)new FirstOrderPropensityArgs(j, K[i*kCols]);
                            propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                        } else {
                            throw InvalidArgException("D", "first order reaction had invalid number of dependencies",numberDependencies);
                        }
                    }
                }
            }
            break;

            // Second order reaction
            case SecondOrderPropensityArgs::REACTION_TYPE:
            {
                // Find the dependencies.
                uint firstDependency = 0;
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1) {
                        numberDependencies++;

                        // Set the table entry to the first two non-zero dependencies.
                        if (numberDependencies == 1) {
                            firstDependency = j;
                        }
                        else if (numberDependencies == 2) {
                            propensityFunctions[i] = (void *)&secondOrderPropensity;
                            propensityFunctionArgs[i] =  (void *)new SecondOrderPropensityArgs(firstDependency, j, K[i*kCols]);
                            propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                        } else {
                            throw InvalidArgException("D", "second order reaction had invalid number of dependencies",numberDependencies);
                        }
                    }
                }
            } 
            break;

            // Second order self-reaction
            case SecondOrderSelfPropensityArgs::REACTION_TYPE:
            {
                // Find the dependencies.
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1) {
                        numberDependencies++;

                        // Set the table entry to the first non-zero dependency.
                        if (numberDependencies == 1) {
                            propensityFunctions[i] = (void *)&secondOrderSelfPropensity;
                            propensityFunctionArgs[i] =  (void *)new SecondOrderSelfPropensityArgs(j, K[i*kCols]);
                            propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                        } else {
                            throw InvalidArgException("D", "second order self reaction had invalid number of dependencies",numberDependencies);
                        }
                    }
                }
            }
            break;

            // Hill Function
            case KHillPropensityArgs::REACTION_TYPE:
            {
                // Find the dependencies.
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1) {
                        numberDependencies++;

                        // Set the table entry to the first non-zero dependency.
                        if (numberDependencies == 1) {
                            propensityFunctions[i] = (void *)&kHillPropensity;
                            propensityFunctionArgs[i] =  (void *)new KHillPropensityArgs(j, K[i*kCols], K[i*kCols+1], K[i*kCols+2], K[i*kCols+3], K[i*kCols+4]);
                            propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                        } else {
                            throw InvalidArgException("D", "khill reaction had invalid number of dependencies",numberDependencies);
                        }
                    }
                }
            }
            break;

            // Hill-Transport Function
            case KHillTransportPropensityArgs::REACTION_TYPE:
            {
                // Find the dependencies.
                uint firstDependency = 0, secondDependency = 0;
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1) {
                        numberDependencies++;

                        // Find the two dependencies.
                        if (numberDependencies == 1)
                            firstDependency = j;
                        else if (numberDependencies == 2)
                            secondDependency = j;
                        else
                            throw InvalidArgException("D", "kinetic hill transport reaction had invalid number of dependencies",numberDependencies);
                    }
                }

                // Figure out which dependency is the controlled species and which is the controlling.
                uint si, xi;
                if (S[firstDependency*numberReactions+i] == -1 && S[secondDependency*numberReactions+i] == 0) {
                    si=firstDependency;
                    xi=secondDependency;
                } else if (S[firstDependency*numberReactions+i] == 0 && S[secondDependency*numberReactions+i] == -1) {
                    si=secondDependency;
                    xi=firstDependency;
                } else {
                    throw InvalidArgException("D", "kinetic hill transport reaction cannot be parsed",numberDependencies);
                }

                // Set up the propensity function.
                if (kCols < 9) throw InvalidArgException("kCols", "kinetic hill transport reaction requires nine K values",kCols);
                propensityFunctions[i] = (void *)&kHillTransportPropensity;
                propensityFunctionArgs[i] =  (void *)new KHillTransportPropensityArgs(si, xi, K[i*kCols], K[i*kCols+1], K[i*kCols+2], K[i*kCols+3], K[i*kCols+4], K[i*kCols+5], K[i*kCols+6], K[i*kCols+7], K[i*kCols+8]);
                propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
            }
            break;

            //  Heaviside Function
            case ZerothOrderHeavisidePropensityArgs::REACTION_TYPE:
            {
                // Find the dependency.
                int xi=-1;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1) {
                    	if (xi != -1) throw InvalidArgException("D", "zeroth order Heaviside reaction can only have one dependency");
                    	xi = j;
                    }
                }

                // Make sure we found the right dependencies.
		    	if (xi == -1) throw InvalidArgException("D", "zeroth order Heaviside reaction must have one dependency");

		    	// Set the table entry.
		    	propensityFunctions[i] = (void *)&zerothOrderHeavisidePropensity;
		    	propensityFunctionArgs[i] =  (void *)new ZerothOrderHeavisidePropensityArgs(xi, (uint)round(K[i*kCols]), K[i*kCols+1], K[i*kCols+2]);
		    	propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
            }
            break;

            // Michaelis-Menten reaction
            case MichaelisMentenPropensityArgs::REACTION_TYPE:
            {
                // Find the dependency
                uint firstDependency = 0, secondDependency = 0;
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1) {
                        numberDependencies++;

                        // Find the two dependencies.
                        if (numberDependencies == 1)
                            firstDependency = j;
                        else if (numberDependencies == 2)
                            secondDependency = j;
                        else
                            throw InvalidArgException("D", "michaelis menten reaction had invalid number of dependencies",numberDependencies);
                    }
                }

                // Figure out which dependency is the controlled species and which is the controlling.
                uint si, ei;
                if (S[firstDependency*numberReactions+i] == -1 && S[secondDependency*numberReactions+i] == 0) {
                    si=firstDependency;
                    ei=secondDependency;
                } else if (S[firstDependency*numberReactions+i] == 0 && S[secondDependency*numberReactions+i] == -1) {
                    si=secondDependency;
                    ei=firstDependency;
                } else {
                    throw InvalidArgException("D", "michaelis menten reaction cannot be parsed",numberDependencies);
                }

                // Set up the propensity function.
                if (kCols < 2) throw InvalidArgException("kCols", "kinetic hill transport reaction requires two K values",kCols);
                propensityFunctions[i] = (void *)&michaelisMentenPropensity;
                propensityFunctionArgs[i] =  (void *)new MichaelisMentenPropensityArgs(si, ei, K[i*kCols], K[i*kCols+1]);
                propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
            }
            break;

            // Competitive Michaelis-Menten reaction
            case CompetitiveMMPropensityArgs::REACTION_TYPE:
            // Uncompetative Michaelis-Menten reaction
            case UncompetitiveMMPropensityArgs::REACTION_TYPE:
            // Noncompetative Michaelis-Menten reaction
            case NoncompetitiveMMPropensityArgs::REACTION_TYPE:
            {
                // Find the dependency
                uint firstDependency = 0, secondDependency = 0, thirdDependency = 0;
                uint numberDependencies = 0;
                for (uint j=0; j<numberSpecies; j++) {
                    if (D[j*numberReactions+i] == 1 or D[j*numberReactions+i] == 2) {
                        numberDependencies++;

                        // Find the two dependencies.
                        if (numberDependencies == 1)
                            firstDependency = j;
                        else if (numberDependencies == 2)
                            secondDependency = j;
                        else if (numberDependencies == 3)
                            thirdDependency = j;
                        else
                            throw InvalidArgException("D", "* michaelis menten reaction had invalid number of dependencies",numberDependencies);
                    }
                }

                // Figure out which dependency is the controlled species and which is the controlling.
                uint si, ei, ii;
                if (S[firstDependency*numberReactions+i] == -1) {
                    si=firstDependency;
                    if(D[secondDependency*numberReactions+i] == 2) {
                        ei=thirdDependency;
                        ii=secondDependency;
                    } else {
                        ei=thirdDependency;
                        ii=secondDependency;
                    }
                } else if (S[secondDependency*numberReactions+i] == -1) {
                    si=secondDependency;
                    if(D[firstDependency*numberReactions+i] == 2) {
                        ei=thirdDependency;
                        ii=firstDependency;
                    } else {
                        ei=firstDependency;
                        ii=thirdDependency;
                    }
                } else if (S[thirdDependency*numberReactions+i] == -1) {
                    si=thirdDependency;
                    if(D[firstDependency*numberReactions+i] == 2) {
                        ei=secondDependency;
                        ii=firstDependency;
                    } else {
                        ei=firstDependency;
                        ii=secondDependency;
                    }
                } else {
                    throw InvalidArgException("D", "* michaelis menten reaction cannot be parsed",numberDependencies);
                }

                // Set up the propensity function.
                if (kCols < 3) throw InvalidArgException("kCols", "competitive michaelis menten reaction requires three K values",kCols);
                if (reactionTypes[i] == CompetitiveMMPropensityArgs::REACTION_TYPE) {
                    propensityFunctions[i] = (void *)&competitiveMMPropensity;
                    propensityFunctionArgs[i] =  (void *)new CompetitiveMMPropensityArgs(si, ei, ii, K[i*kCols], K[i*kCols+1], K[i*kCols+2]);
                    propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                } else if(reactionTypes[i] == NoncompetitiveMMPropensityArgs::REACTION_TYPE) {
                    propensityFunctions[i] = (void *)&uncompetitiveMMPropensity;
                    propensityFunctionArgs[i] =  (void *)new UncompetitiveMMPropensityArgs(si, ei, ii, K[i*kCols], K[i*kCols+1], K[i*kCols+2]);
                    propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                } else if(reactionTypes[i] == NoncompetitiveMMPropensityArgs::REACTION_TYPE) {
                    propensityFunctions[i] = (void *)&noncompetitiveMMPropensity;
                    propensityFunctionArgs[i] =  (void *)new NoncompetitiveMMPropensityArgs(si, ei, ii, K[i*kCols], K[i*kCols+1], K[i*kCols+2]);
                    propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[i]);
                } else {
                    throw InvalidArgException("D", "* michaelis menten reaction");
                }
            }
            break;

            default:
                throw InvalidArgException("reactionType", "Unknown reaction type identifier.");

        } // End switch(reactionTypes[i])
    }

    // Create the species dependency tables from the S matrix.
    for (uint i=0; i<numberReactions; i++)
    {
        numberDependentSpecies[i]=0;
        for (uint j=0, index=i; j<numberSpecies; j++, index+=numberReactions)
            if (S[index] != 0)
                numberDependentSpecies[i]++;
        dependentSpecies[i] = new uint[numberDependentSpecies[i]];
        dependentSpeciesChange[i] = new int[numberDependentSpecies[i]];
        for (uint j=0, index=i, k=0; j<numberSpecies; j++, index+=numberReactions)
        {
            if (S[index] != 0 && k < numberDependentSpecies[i])
            {
                dependentSpecies[i][k] = j;
                dependentSpeciesChange[i][k] = S[index];
                k++;
            }
        }
    }

    // Create the reaction dependency tables from the other tables.
    for (uint r=0; r<numberReactions; r++)
    {
        list<uint> dependentReactionList;

        // Go through all of the species changed by this reaction.
        for (uint d=0; d<numberDependentSpecies[r]; d++)
        {
            uint s = dependentSpecies[r][d];

            // Find all of the reactions that depend on this species.
            for (uint i=0, index=s*numberReactions; i<numberReactions; i++, index++)
            {
                if (D[index] > 0) dependentReactionList.push_back(i);
            }
        }

        // Eliminate any duplicates from the list.
        dependentReactionList.sort();
        dependentReactionList.unique();

        // Create the table.
        numberDependentReactions[r] = dependentReactionList.size();
        dependentReactions[r] = new uint[numberDependentReactions[r]];
        uint i=0;
        for (list<uint>::iterator it=dependentReactionList.begin(); it != dependentReactionList.end() && i<numberDependentReactions[r]; it++, i++)
        {
            dependentReactions[r][i] = *it;
        }
    }
}

void CMESolver::getSpeciesCountView(uint **counts, int *number)
{
    *counts = speciesCounts;
    *number = (int)numberSpecies;
}

void CMESolver::getReactionRateConstantsView(int reactionNumber, double **rates, int *rateConstantCount)
{
	if(reactionNumber < 0 || reactionNumber >= (int)numberReactions)
        throw InvalidArgException("reactionNumber","reaction index is out of bounds.", reactionNumber);
    PropensityArgs *prop = (PropensityArgs*)propensityFunctionArgs[reactionNumber];

    switch(reactionTypes[reactionNumber]) {
		case ZerothOrderPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((ZerothOrderPropensityArgs *)prop)->k);
		    *rateConstantCount = 1;
		}
		break;

		case FirstOrderPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((FirstOrderPropensityArgs *)prop)->k);
		    *rateConstantCount = 1;

		}
		break;

		case SecondOrderPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((SecondOrderPropensityArgs *)prop)->k);
		    *rateConstantCount = 1;
		}
		break;

		case SecondOrderSelfPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((SecondOrderSelfPropensityArgs *)prop)->k);
		    *rateConstantCount = 1;
		}
		break;

		case KHillPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((KHillPropensityArgs *)prop)->k);
		    *rateConstantCount = 1;
		}
		break;

		case KHillTransportPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((KHillTransportPropensityArgs *)prop)->k0);
		    *rateConstantCount = 5;
		}
		break;

		case ZerothOrderHeavisidePropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((ZerothOrderHeavisidePropensityArgs *)prop)->k0);
		    *rateConstantCount = 2;
		}
		break;

		case MichaelisMentenPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((MichaelisMentenPropensityArgs *)prop)->kcat);
		    *rateConstantCount = 2;
		}
		break;

		case CompetitiveMMPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((CompetitiveMMPropensityArgs *)prop)->kcat);
		    *rateConstantCount = 3;
		}
		break;

		case UncompetitiveMMPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((UncompetitiveMMPropensityArgs *)prop)->kcat);
		    *rateConstantCount = 3;
		}
		break;

		case NoncompetitiveMMPropensityArgs::REACTION_TYPE:
		{
		    *rates = &(((NoncompetitiveMMPropensityArgs *)prop)->kcat);
		    *rateConstantCount = 3;
		}
		break;

		default:
			throw InvalidArgException("reactionType","Unknown reaction type identifier.");
	}
}


////////////////////////////
// Propensity Calculators //
////////////////////////////
double CMESolver::zerothOrderPropensity(double time, uint * speciesCounts, void * pargs)
{
    ZerothOrderPropensityArgs * args = (ZerothOrderPropensityArgs *)pargs;
    return args->k;
}

double CMESolver::firstOrderPropensity(double time, uint * speciesCounts, void * pargs)
{
    FirstOrderPropensityArgs * args = (FirstOrderPropensityArgs *)pargs;
    return args->k * (double)speciesCounts[args->si];
}

double CMESolver::secondOrderPropensity(double time, uint * speciesCounts, void * pargs)
{
    SecondOrderPropensityArgs * args = (SecondOrderPropensityArgs *)pargs;
    return args->k * ((double)speciesCounts[args->s1i]) * ((double)speciesCounts[args->s2i]);
}

double CMESolver::secondOrderSelfPropensity(double time, uint * speciesCounts, void * pargs)
{
    SecondOrderSelfPropensityArgs * args = (SecondOrderSelfPropensityArgs *)pargs;
    uint count = speciesCounts[args->si];
    if (count >= 2)
		// We divide this number by 2 as the two particles cannot be distinguished
        return 0.5 * args->k * ((double)count) * ((double)(count-1));
    return 0.0;
}

double CMESolver::kHillPropensity(double time, uint * speciesCounts, void * pargs)
{
	KHillPropensityArgs * args = (KHillPropensityArgs *)pargs;
    return args->k * (double)speciesCounts[args->si];
}

double CMESolver::kHillTransportPropensity(double time, uint * speciesCounts, void * pargs)
{
    KHillTransportPropensityArgs * args = (KHillTransportPropensityArgs *)pargs;

    uint s = speciesCounts[args->si];
    if (s == 0) return 0.0;

    double x = (double)speciesCounts[args->xi];
    double xh = pow(1+(args->ITp*x),args->h);
    double p = args->k0+((args->dk*xh)/(args->IRh+xh));

    if (time == 0.0)
        Print::printf(Print::DEBUG, "Recalculating hill transport propensity for %d,%d with %e,%e,%e,%e,%e = %e", args->si, args->xi, args->k0, args->dk, args->IRh, args->ITp, args->h, p);

    return p;
}

double CMESolver::zerothOrderHeavisidePropensity(double time, uint * speciesCounts, void * pargs)
{
	ZerothOrderHeavisidePropensityArgs * args = (ZerothOrderHeavisidePropensityArgs *)pargs;

	//Print::printf(Print::DEBUG, "Recalculating zerothOrderHeavisidePropensity for %d (count=%d) with %d,%e,%e = %e", args->xi, speciesCounts[args->xi], args->x0, args->k0, args->k1, ((speciesCounts[args->xi]<args->x0)?(args->k0):(args->k1)));

	return ((speciesCounts[args->xi]<args->x0)?(args->k0):(args->k1));
}

double CMESolver::michaelisMentenPropensity(double time, uint * speciesCounts, void *pargs) 
{
    MichaelisMentenPropensityArgs *args = (MichaelisMentenPropensityArgs *)pargs;

    double E = (double)speciesCounts[args->si];
    double S = (double)speciesCounts[args->ei];

    return args->kcat*E*S/ (args->Km + S);
}

double CMESolver::competitiveMMPropensity(double time, uint * speciesCounts, void *pargs)
{
    CompetitiveMMPropensityArgs *args = (CompetitiveMMPropensityArgs *)pargs;

    double S = (double)speciesCounts[args->si];
    double E = (double)speciesCounts[args->ei];
    double I = (double)speciesCounts[args->ii];

    return args->kcat*E*S/(args->Km * (1.0 + I/args->Ki) + S);
}

double CMESolver::uncompetitiveMMPropensity(double time, uint * speciesCounts, void *pargs)
{
    UncompetitiveMMPropensityArgs *args = (UncompetitiveMMPropensityArgs * )pargs;

    double S = (double)speciesCounts[args->si];
    double E = (double)speciesCounts[args->ei];
    double I = (double)speciesCounts[args->ii];

    return args->kcat*E*S/(args->Km + S + S*I/args->Ki);
}

double CMESolver::noncompetitiveMMPropensity(double time, uint * speciesCounts, void *pargs)
{
    NoncompetitiveMMPropensityArgs *args = (NoncompetitiveMMPropensityArgs *)pargs;

    double S = (double)speciesCounts[args->si];
    double E = (double)speciesCounts[args->ei];
    double I = (double)speciesCounts[args->ii];

    return args->kcat*(E/(1.0+I/args->Ki))*S/(args->Km + S);
}


void CMESolver::setModelPropensityFunction(uint reaction, double (*propensityFunction)(double time, uint * speciesCounts, void * args), void * propensityFunctionArg)
{
    if (reaction >= numberReactions) throw InvalidArgException("reaction", "reaction index exceeded the number of reactions");
    propensityFunctions[reaction] = (void *)propensityFunction;
    propensityFunctionArgs[reaction] = propensityFunctionArg;
}

void CMESolver::setSpeciesUpperLimit(uint species, uint limit)
{
    // Allocate a larger list for the limits.
    SpeciesLimit * newSpeciesLimits = new SpeciesLimit[numberSpeciesLimits++];
    if (numberSpeciesLimits > 1)
    {
        for (uint i=0; i<numberSpeciesLimits-1; i++)
            newSpeciesLimits[i] = speciesLimits[i];
        delete[] speciesLimits;
    }
    speciesLimits = newSpeciesLimits;
    speciesLimits[numberSpeciesLimits-1].type = 1;
    speciesLimits[numberSpeciesLimits-1].species = species;
    speciesLimits[numberSpeciesLimits-1].limit = limit;
}

void CMESolver::setSpeciesLowerLimit(uint species, uint limit)
{
    // Allocate a larger list for the limits.
    SpeciesLimit * newSpeciesLimits = new SpeciesLimit[numberSpeciesLimits++];
    if (numberSpeciesLimits > 1)
    {
        for (uint i=0; i<numberSpeciesLimits-1; i++)
            newSpeciesLimits[i] = speciesLimits[i];
        delete[] speciesLimits;
    }
    speciesLimits = newSpeciesLimits;
    speciesLimits[numberSpeciesLimits-1].type = -1;
    speciesLimits[numberSpeciesLimits-1].species = species;
    speciesLimits[numberSpeciesLimits-1].limit = limit;
}



/////////////////////////////
// First Passage Time Code //
/////////////////////////////
void CMESolver::setFptTrackingList(list<uint> speciesList)
{
    // If we already had a list, free it.
    if (fptTrackedSpecies != NULL) {delete[] fptTrackedSpecies; fptTrackedSpecies = NULL;}

    // Allocate a new list.
    numberFptTrackedSpecies = speciesList.size();
    fptTrackedSpecies = new FPTTracking[numberFptTrackedSpecies];

    // Initialize the list.
    int i=0;
    for (list<uint>::iterator it = speciesList.begin(); it != speciesList.end(); it++, i++)
    {
        fptTrackedSpecies[i].species = *it;
        fptTrackedSpecies[i].minValueAchieved = 0;
        fptTrackedSpecies[i].maxValueAchieved = 0;
        fptTrackedSpecies[i].dataSet.set_species(fptTrackedSpecies[i].species);
    }
}

void CMESolver::addToParameterTrackingList(pair<string,double*> parameter)
{
    trackedParameters.push_back(TrackedParameter(parameter.first, parameter.second));
}

double CMESolver::recordParameters(double nextRecordTime, double recordInterval, double simulationTime)
{
	if (recordInterval > 0.0)
	{
		// Write parameter values until the next write time is past the current time.
		do
		{
			// Update the parameter values.
			for (list<TrackedParameter>::iterator it = trackedParameters.begin(); it != trackedParameters.end(); it++)
			{
				// Record the species counts.
				it->dataSet.add_time(nextRecordTime);
				it->dataSet.add_value(*(it->valuePointer));
			}
			nextRecordTime += recordInterval;
		}
		while (nextRecordTime <= (simulationTime+1e-9));

		// Output the data, if necessary.
		queueRecordedParameters(false);
	}

    return nextRecordTime;
}

void CMESolver::queueRecordedParameters(bool flush)
{
    for (list<TrackedParameter>::iterator it = trackedParameters.begin(); it != trackedParameters.end(); it++)
    {
        if (it->dataSet.value_size() >= TUNE_PARAMETER_VALUES_BUFFER_SIZE || (flush && it->dataSet.value_size() > 0))
        {
            // Push it to the output queue.
            PROF_BEGIN(PROF_SERIALIZE_PV);
            lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::PARAMETER_VALUES, replicate, &it->dataSet);
            PROF_END(PROF_SERIALIZE_PV);

            // Reset the data set.
            it->dataSet.Clear();
            it->dataSet.set_parameter(it->name);
        }
    }
}

int CMESolver::hookSimulation(double time)
{
    // Overload this function in derivative classes
    // Return 0 if the lattice state is unchanged
    // Return 1 if the lattice state has been modified, 
    //          and it needs to be copied back to the GPU.
    return 0;
}

int CMESolver::onBeginTrajectory()
{
	return 0;
}

int CMESolver::onEndTrajectory()
{
	return 0;
}

}
}
