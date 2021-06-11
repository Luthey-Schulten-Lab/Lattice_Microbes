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
#include "cme/FluctuatingNRSolver.h"
#include "ReactionModel.pb.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"




using std::string;
using std::list;
using std::map;
using lm::reaction::ReactionQueue;
using lm::rng::RandomGenerator;

namespace lm {
namespace cme {

FluctuatingNRSolver::FluctuatingNRSolver():NextReactionSolver((RandomGenerator::Distributions)(RandomGenerator::NORMAL))
{
}

FluctuatingNRSolver::~FluctuatingNRSolver()
{
}

void FluctuatingNRSolver::setReactionModel(lm::io::ReactionModel * rm)
{
    if (rm->number_reactions() != (uint)rm->reaction_size()) throw InvalidArgException("rm", "number of reaction does not agree with reaction list size");

    // Figure out the max number of columns we need in the k matrix.
    uint kCols = 0;
    for (uint i=0; i<rm->number_reactions(); i++)
        kCols = max(kCols,(uint)rm->reaction(i).rate_constant_size());

    // Set the K and reaction type tables.
    uint * reactionType = new uint[rm->number_reactions()];
    double * K = new double[rm->number_reactions()*kCols];
    double * nvar = new double[rm->number_reactions()];
    double * ntau = new double[rm->number_reactions()];
    bool hasNoisyRates = false;
    for (uint i=0; i<rm->number_reactions(); i++)
    {
        reactionType[i] = rm->reaction(i).type();
        for (uint j=0; j<(uint)rm->reaction(i).rate_constant_size(); j++)
        {
            K[i*kCols+j] = rm->reaction(i).rate_constant(j);
        }
        if (rm->reaction(i).rate_has_noise())
        {
            hasNoisyRates = true;
            nvar[i] = rm->reaction(i).rate_noise_variance();
            ntau[i] = rm->reaction(i).rate_noise_tau();
        }
        else
        {
            nvar[i] = 0.0;
            ntau[i] = 0.0;
        }
    }

    // If there are any noisy rates, call the alternate builder.
    if (hasNoisyRates)
    {
        // Make sure we have a non-zero recalc rate.
        double noiseRecalcFraction=atof((*parameters)["noiseRecalcFraction"].c_str());;
        if (noiseRecalcFraction <= 0.0) throw Exception("The noise recalculation rate must be greater than zero.");

        // Build the model.
        buildModel(rm->number_species(), rm->number_reactions(), rm->initial_species_count().data(), reactionType, K, rm->stoichiometric_matrix().data(), rm->dependency_matrix().data(), nvar, ntau, noiseRecalcFraction, kCols);
    }
    else
    {
        // Build the model.
        buildModel(rm->number_species(), rm->number_reactions(), rm->initial_species_count().data(), reactionType, K, rm->stoichiometric_matrix().data(), rm->dependency_matrix().data(), kCols);
    }

    // Free any resources.
    if (ntau !=  NULL) {delete [] ntau; ntau = NULL;}
    if (nvar !=  NULL) {delete [] nvar; nvar = NULL;}
    if (K !=  NULL) {delete [] K; K = NULL;}
    if (reactionType !=  NULL) {delete [] reactionType; reactionType = NULL;}
}

void FluctuatingNRSolver::buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionTypeA, const double * KA, const int * SA, const uint * DA, const uint kCols)
{
    NextReactionSolver::buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypeA, KA, SA, DA, kCols);
}

void FluctuatingNRSolver::buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionTypeA, const double * KA, const int * SA, const uint * DA, const double * nvar, const double * ntau, const double noiseRecalcFraction, const uint kCols)
{
    // Figure out how many fluctuating reactions we have.
    uint numberFluctuatingReactions=0;
    for (uint i=0; i<numberReactionsA; i++)
    {
        if (nvar[i] > 0.0 && ntau[i] > 0.0)
            numberFluctuatingReactions++;
    }

    // Figure out the reaction index for each fluctuating reaction.
    uint * fluctuatingReactionIndex = new uint[numberFluctuatingReactions];
    numberFluctuatingReactions=0;
    for (uint i=0; i<numberReactionsA; i++)
    {
        if (nvar[i] > 0.0 && ntau[i] > 0.0)
            fluctuatingReactionIndex[numberFluctuatingReactions++] = i;
    }

    // Copy the parameters while expanding them to include space for the pseudo-reactions.
    uint numberSpecies = numberSpeciesA + numberFluctuatingReactions;
    uint * initialSpeciesCounts = new uint[numberSpecies];
    for (uint j=0; j<numberSpecies; j++)
    {
        if (j <numberSpeciesA)
            initialSpeciesCounts[j] = initialSpeciesCountsA[j];
        else
            initialSpeciesCounts[j] = 0;
    }
    uint numberReactions = numberReactionsA + numberFluctuatingReactions;
    uint * reactionType = new uint[numberReactions];
    double * K = new double[numberReactions*kCols];
    for (uint i=0; i<numberReactions; i++)
    {
        if (i<numberReactionsA)
        {
        	// Copy the reaction parameters.
			reactionType[i] = reactionTypeA[i];
			for (uint j=0; j<kCols; j++)
				K[i*kCols+j] = KA[i*kCols+j];
        }
        else
        {
            reactionType[i] = ZerothOrderPropensityArgs::REACTION_TYPE;
            K[i*kCols] = 1/(ntau[fluctuatingReactionIndex[i-numberReactionsA]]*noiseRecalcFraction);
        }
    }
    int * S = new int [numberReactions*numberSpecies];
    uint * D = new uint [numberReactions*numberSpecies];
    for (uint i=0; i<numberReactions; i++)
    {
        for (uint j=0; j<numberSpecies; j++)
        {
            if (i < numberReactionsA && j < numberSpeciesA)
            {
                S[j*numberReactions+i] = SA[j*numberReactionsA+i];
                D[j*numberReactions+i] = DA[j*numberReactionsA+i];
            }
            else
            {
                S[j*numberReactions+i] = 0;
                D[j*numberReactions+i] = 0;
            }
        }
    }

    // Fill in the diagonal of the new portion of the S matrix to associate each OU reaction with increasing its OU counter species.
    for (uint i=numberReactionsA, j=numberSpeciesA; i<numberReactions && j<numberSpecies; i++, j++)
    {
        S[j*numberReactions+i] = 1;
    }

    // Fill in the dependency matrix such that each fluctuating reaction will have its propensity recomputed when its OU counter species is increased.
    for (uint i=0; i<numberFluctuatingReactions; i++)
    {
        uint reactionIndex = fluctuatingReactionIndex[i];

        // Make these be secondary dependencies, such that they will not be used in propensity calculations.
        D[(numberSpeciesA+i)*numberReactions+reactionIndex] = 2;
    }

    // Build the model with the expanded parameters.
    NextReactionSolver::buildModel(numberSpecies, numberReactions, initialSpeciesCounts, reactionType, K, S, D, kCols);

    // Add the custom reaction propensity arguments and functions.
    for (uint i=0; i<numberFluctuatingReactions; i++)
    {
        uint reactionIndex = fluctuatingReactionIndex[i];

        // Figure out the variable to which the noise should be applied.
        double * noisyK = NULL;
        if (reactionType[reactionIndex] == ZerothOrderPropensityArgs::REACTION_TYPE)
        {
        	noisyK = &(((ZerothOrderPropensityArgs *)propensityFunctionArgs[reactionIndex])->k);
        }
        else if (reactionType[reactionIndex] == FirstOrderPropensityArgs::REACTION_TYPE)
        {
        	noisyK = &(((FirstOrderPropensityArgs *)propensityFunctionArgs[reactionIndex])->k);
        }
        else if (reactionType[reactionIndex] == ZerothOrderHeavisidePropensityArgs::REACTION_TYPE)
        {
        	noisyK = &(((ZerothOrderHeavisidePropensityArgs *)propensityFunctionArgs[reactionIndex])->k1);
        }
        else
        {
            throw InvalidArgException("nvar","noise on reactions of the specified type is not currently supported", reactionType[reactionIndex]);
        }

        // Override the existing propensity function and arg pointer with an ou propensity wrapper.
        propensityFunctionArgs[reactionIndex] =  (void *)new OUPropensityArgs(propensityFunctions[reactionIndex], propensityFunctionArgs[reactionIndex], noisyK, *noisyK, numberSpeciesA+i, nvar[reactionIndex], ntau[reactionIndex], rng);
        propensityArgs.push_back((PropensityArgs *)propensityFunctionArgs[reactionIndex]);
        propensityFunctions[reactionIndex] = (void *)&ouPropensity;
        Print::printf(Print::INFO, "Set reaction %d to be a noisy reaction with var=%e, tau=%e, and recalc fraction %e. Recalcs tracked in species %d", reactionIndex, nvar[reactionIndex], ntau[reactionIndex], noiseRecalcFraction, numberSpeciesA+i);

        // Add the noise variable to the parameter tracking list.
        char nameBuffer[8];
        snprintf(nameBuffer, sizeof(nameBuffer), "%07d", reactionIndex);
        addToParameterTrackingList(pair<string,double*>(nameBuffer,&(((OUPropensityArgs *)propensityFunctionArgs[reactionIndex])->noise)));
    }

    // Change the number of species to track so we don't track the pseduo-species.
    numberSpeciesToTrack = numberSpecies-numberFluctuatingReactions;

    // Free any resources.
    delete [] initialSpeciesCounts;
    delete [] reactionType;
    delete [] K;
    delete [] S;
    delete [] D;
    delete [] fluctuatingReactionIndex;
}

void FluctuatingNRSolver::destroyModel()
{
    NextReactionSolver::destroyModel();

    // Free propensity args.
}

double FluctuatingNRSolver::ouPropensity(double time, uint * speciesCounts, void * pargs)
{
    OUPropensityArgs * args = (OUPropensityArgs *)pargs;

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
            Print::printf(Print::DEBUG, "Initialized OU 1st order noise at %0.4e (JN %d): %0.4e", time, ouJumpNumber, args->noise);
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

            // Update the noisy k with its new value.
            *(args->noisyK) = args->noisyKInitialValue + args->noise;
            if (*(args->noisyK) < 0.0) *(args->noisyK) = 0.0;
            Print::printf(Print::DEBUG, "Recalculated OU 1st order noise at %0.4e (JN %d): k=%0.4e (%0.4e+%0.4e)", time, ouJumpNumber, *(args->noisyK), args->noisyKInitialValue, args->noise);
        }
    }

    // Perform the base propensity calculation.
    double (*propensityFunction)(double, uint * speciesCounts, void * args) = (double (*)(double, uint*, void*))args->basePropensityFunction;
    return (*propensityFunction)(time, speciesCounts, args->basePropensityFunctionArgs);
}

}
}
