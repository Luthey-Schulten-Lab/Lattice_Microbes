/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
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

#include <cuda.h>
#include <cuda_runtime.h>
#include "config.h"
#include "cuda/ldg.h"

#if !defined MPD_WORDS_PER_SITE
#error "Must define the number of words per site."
#endif

#define MPD_PARTICLES_PER_SITE MPD_WORDS_PER_SITE

#if !defined TUNE_MPD_MAX_PARTICLE_OVERFLOWS
#error "Must define the maximum size of overflow particle list."
#endif

#define MPD_ZERO_ORDER_REACTION         1
#define MPD_FIRST_ORDER_REACTION        2
#define MPD_SECOND_ORDER_REACTION       3
#define MPD_SECOND_ORDER_SELF_REACTION  4

#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
inline __device__ float calculateReactionPropensity(const uint8_t siteType, const unsigned int * __restrict__ particles, const unsigned int reactionIndex, const uint8_t __restrict__ *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
#else
inline __device__ float calculateReactionPropensity(const uint8_t siteType, const unsigned int * __restrict__ particles, const unsigned int reactionIndex, const uint8_t __restrict__ *RLG)
#endif
#else
inline __device__ float calculateReactionPropensity(const uint8_t siteType, const unsigned int * __restrict__ particles, const unsigned int reactionIndex)
#endif
{
    // Make sure that the reaction is valid for this site type.
#ifdef MPD_GLOBAL_S_MATRIX
    if (!read_element(RLG,reactionIndex*numberSiteTypesC+siteType)) return 0.0f;
#else
    if (!RLC[reactionIndex*numberSiteTypesC+siteType]) return 0.0f;
#endif

#ifdef MPD_GLOBAL_R_MATRIX
    // Get the number of each reactant.
    float numberParticles1 = 0.0f;
    float numberParticles2 = 0.0f;
    for (int i=0; i<MPD_PARTICLES_PER_SITE; i++)
    {
        if (particles[i] > 0)
        {
            numberParticles1 += (particles[i] == D1G[reactionIndex])?(1.0f):(0.0f);
            numberParticles2 += (particles[i] == D2G[reactionIndex])?(1.0f):(0.0f);
        }
    }

    // Calculate the propensity according to the reaction order.
    if (reactionOrdersG[reactionIndex] == MPD_ZERO_ORDER_REACTION)
        return reactionRatesG[reactionIndex]; // NOTE: The reactionRatesC is already converted to molecules/sec.  This is performed in "MpdRdmeSolver::buildDiffusionModel"
    else if (reactionOrdersG[reactionIndex] == MPD_FIRST_ORDER_REACTION)
        return reactionRatesG[reactionIndex]*numberParticles1;
    else if (reactionOrdersG[reactionIndex] == MPD_SECOND_ORDER_REACTION)
        return reactionRatesG[reactionIndex]*numberParticles1*numberParticles2;
    else if (reactionOrdersG[reactionIndex] == MPD_SECOND_ORDER_SELF_REACTION)
        return reactionRatesG[reactionIndex]*numberParticles1*(numberParticles1-1.0f);
    return 0.0f;

#else

// Get the number of each reactant.
    float numberParticles1 = 0.0f;
    float numberParticles2 = 0.0f;
    for (int i=0; i<MPD_PARTICLES_PER_SITE; i++)
    {
        if (particles[i] > 0)
        {
            numberParticles1 += (particles[i] == D1C[reactionIndex])?(1.0f):(0.0f);
            numberParticles2 += (particles[i] == D2C[reactionIndex])?(1.0f):(0.0f);
        }
    }

    // Calculate the propensity according to the reaction order.
    if (reactionOrdersC[reactionIndex] == MPD_ZERO_ORDER_REACTION)
        return reactionRatesC[reactionIndex]; // NOTE: The reactionRatesC is already converted to molecules/sec.  This is performed in "MpdRdmeSolver::buildDiffusionModel"
    else if (reactionOrdersC[reactionIndex] == MPD_FIRST_ORDER_REACTION)
        return reactionRatesC[reactionIndex]*numberParticles1;
    else if (reactionOrdersC[reactionIndex] == MPD_SECOND_ORDER_REACTION)
        return reactionRatesC[reactionIndex]*numberParticles1*numberParticles2;
    else if (reactionOrdersC[reactionIndex] == MPD_SECOND_ORDER_SELF_REACTION)
        return reactionRatesC[reactionIndex]*numberParticles1*(numberParticles1-1.0f);
    return 0.0f;

#endif

}

inline __device__ float calculateReactionProbability(const float rate)
{
    #ifdef RNG_CUDA_DOUBLE_PRECISION
    return (float)(1.0-exp(-(double)rate));
    #else
    return (rate > 2e-4f)?(1.0f-__expf(-rate)):(rate);
    #endif
}

inline __device__ unsigned int checkForReaction(const unsigned int latticeIndex, const float reactionProbability, const unsigned long long timestepHash)
{
    return reactionProbability > 0.0f && getRandomHashFloat(latticeIndex, 1, 0, timestepHash) <= reactionProbability;
}

#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
inline __device__ unsigned int determineReactionIndex(const uint8_t siteType, const unsigned int * __restrict__ particles, const unsigned int latticeIndex, const float totalReactionPropensity, const unsigned long long timestepHash, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
#else
inline __device__ unsigned int determineReactionIndex(const uint8_t siteType, const unsigned int * __restrict__ particles, const unsigned int latticeIndex, const float totalReactionPropensity, const unsigned long long timestepHash, const __restrict__ uint8_t *RLG)
#endif	
#else
inline __device__ unsigned int determineReactionIndex(const uint8_t siteType, const unsigned int * __restrict__ particles, const unsigned int latticeIndex, const float totalReactionPropensity, const unsigned long long timestepHash)
#endif
	{
		float randomPropensity = getRandomHashFloat(latticeIndex, 1, 1, timestepHash)*totalReactionPropensity;
		unsigned int reactionIndex = 0;
		for (int i=0; i<numberReactionsC; i++)
		{
	#ifdef MPD_GLOBAL_S_MATRIX
        #ifdef MPD_GLOBAL_R_MATRIX
                        float propensity = calculateReactionPropensity(siteType, particles, i, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
        #else
			float propensity = calculateReactionPropensity(siteType, particles, i, RLG);
	#endif
        #else
			float propensity = calculateReactionPropensity(siteType, particles, i);
	#endif
			if (propensity > 0.0f)
			{
				if (randomPropensity > 0.0f)
					reactionIndex = i;
				randomPropensity -= propensity;
			}

		}
		return reactionIndex;
	}

	#ifdef MPD_GLOBAL_S_MATRIX
	inline __device__ void evaluateReaction(const unsigned int latticeIndex, const uint8_t siteType, unsigned int * __restrict__ particles, const unsigned int reactionIndex, unsigned int * siteOverflowList, const int8_t* __restrict__ SG)
	#else
	inline __device__ void evaluateReaction(const unsigned int latticeIndex, const uint8_t siteType, unsigned int * __restrict__ particles, const unsigned int reactionIndex, unsigned int * siteOverflowList)
	#endif
	{

		int8_t *S;
		// Copy the S matrix entries for this reaction.
		#if __CUDA_ARCH__ >= 200
		S = (int8_t *)malloc(numberSpeciesC);
		#else
			#pragma error CUDA_ARCH >= 200 required
    #endif
	//for (uint i=0, index=reactionIndex; i<numberSpeciesC; i++, index+=numberReactionsC)
	for (uint i=0, index=reactionIndex*numberSpeciesC; i<numberSpeciesC; i++, index++)
	{
#ifdef MPD_GLOBAL_S_MATRIX
		S[i] = read_element((char*)SG,index);
#else
		S[i] = SC[index];
#endif
	}

    // Build the new site, copying particles that didn't react and removing those that did.
    int nextParticle=0;
    for (uint i=0; i<MPD_PARTICLES_PER_SITE; i++)
    {
    	unsigned int particle = particles[i];
        if (particle > 0)
        {
        	// If this particle was unaffected, copy it.
        	if (S[particle-1] >= 0)
        	{
        		particles[nextParticle++] = particle;
        	}

            // Otherwise, don't copy the particle and mark that we destroyed it.
        	else
        	{
        		S[particle-1]++;
        	}
        }
    }

    // Go through the S matrix and add in any new particles that were created.
    for (uint i=0; i<numberSpeciesC; i++)
    {
		for (uint j=0; j<S[i]; j++)
		{
			// If the particle will fit into the site, add it.
			if (nextParticle < MPD_PARTICLES_PER_SITE)
			{
				particles[nextParticle++] = i+1;
			}

			// Otherwise add it to the exception list.
			else
			{
				int exceptionIndex = atomicAdd(siteOverflowList, 1);
				if (exceptionIndex < TUNE_MPD_MAX_PARTICLE_OVERFLOWS)
				{
					siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
					siteOverflowList[(exceptionIndex*2)+2]=i+1;
				}
			}

		}
    }

    // Clear any remaining particles in the site.
    while (nextParticle < MPD_PARTICLES_PER_SITE)
    	particles[nextParticle++] = 0;

    // Free any allocated memory.
    #if __CUDA_ARCH__ >= 200
    free(S);
    #endif
}
