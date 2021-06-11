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

#include <cuda.h>
#include <cuda_runtime.h>
#include "config.h"
#include "cuda/ldg.h"

#if !defined TUNE_MPD_MAX_PARTICLE_OVERFLOWS
#error "Must define the maximum size of overflow particle list."
#endif

#define MPD_OVERFLOW_LIST_SIZE       1+2*TUNE_MPD_MAX_PARTICLE_OVERFLOWS*sizeof(uint32_t)

#define MPD_PARTICLES_PER_SITE  MPD_WORDS_PER_SITE*4
#if MPD_PARTICLES_PER_SITE <= 4
#define MPD_PARTICLE_COUNT_BITS 2
#elif MPD_PARTICLES_PER_SITE <= 8
#define MPD_PARTICLE_COUNT_BITS 3
#elif 1
#define MPD_PARTICLE_COUNT_BITS 4
#else
#error "Unsupported number of particles per site."
#endif

#define MPD_MOVE_NONE           0
#define MPD_MOVE_STAY           1
#define MPD_MOVE_PLUS           2
#define MPD_MOVE_MINUS          3


#include "cuda/constant.cuh"

inline __device__ float lookupTransitionProbability(const unsigned int particleType, const unsigned char sourceSite, const unsigned char destSite)
{
    if ((particleType-1) >= numberSpeciesC || destSite >= numberSiteTypesC || sourceSite >= numberSiteTypesC) return 0.0;
#ifdef MPD_GLOBAL_T_MATRIX
    return read_element(TC, sourceSite*numberSiteTypesC*numberSpeciesC + destSite*numberSpeciesC + (particleType-1));
#else
    return TC[sourceSite*numberSiteTypesC*numberSpeciesC + destSite*numberSpeciesC + (particleType-1)];
#endif
}

/**
 * Gets the random choices for a one-dimension diffusion move.
 *
 * @return  For each bit-packed each particle: 0 = no particle, 1 = stay in place, 2 = move in the positive direction,
 *          4 = move in the minus direction.
 */
inline __device__ void makeDiffusionChoices(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int windowIndexMinus, const unsigned int windowIndex, const unsigned int windowIndexPlus, const unsigned int windowSize, const unsigned long long timestepHash)
{
	const unsigned char site = sitesWindow[windowIndex];
	const unsigned char siteMinus = sitesWindow[windowIndexMinus];
	const unsigned char sitePlus = sitesWindow[windowIndexPlus];

	for(int w=0; w < MPD_WORDS_PER_SITE; w++, window = window+windowSize, choices = choices+windowSize)
	{
		// Set the default choice to none.
		choices[windowIndex]=0;

		// If there are no particles, we are done.
		if (window[windowIndex] > 0)
		{
			// Get the particles as a byte pointer.
			const unsigned char* particles = (const unsigned char*)&window[windowIndex];
			unsigned char* particleChoices = (unsigned char*)&choices[windowIndex];

			// Make choices for all of the particles in the site.
			for (int i=0; i<4; i++)
			{
				if (particles[i] > 0)
				{
					// Get the probability of moving plus and minus.
					float probMinus=lookupTransitionProbability(particles[i], site, siteMinus);
					float probPlus=lookupTransitionProbability(particles[i], site, sitePlus);

					// Get the random value and see which direction we should move in.
					float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, i + w*4, timestepHash);
					particleChoices[i] = (randomValue < probMinus)?(MPD_MOVE_MINUS):(MPD_MOVE_STAY);
					particleChoices[i] = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_PLUS):(particleChoices[i]);

				}
			}
		}
	}
}

/**
 * Makes the X diffusion choices for all particles in loaded sites further than one from the edge.
 */
inline __device__ void makeXDiffusionChoices(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex, const unsigned int blockXSize, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(window, sitesWindow, choices, latticeIndex, windowIndex-1, windowIndex, windowIndex+1, MPD_X_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowIndex >= blockXSize)
    {
        unsigned int apronLatticeIndex = (latticeXIndex>=blockXSize)?latticeIndex-blockXSize:latticeIndex+(latticeXSizeC-blockXSize);
        unsigned int apronWindowIndex = windowIndex-blockXSize;
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex, apronWindowIndex, apronWindowIndex+1, MPD_X_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeXIndex<latticeXSizeC-blockXSize)?latticeIndex+blockXSize:latticeIndex-latticeXSizeC+blockXSize;
        unsigned int apronWindowIndex = windowIndex+blockXSize;
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex-1, apronWindowIndex, apronWindowIndex, MPD_X_WINDOW_SIZE, timestepHash);
    }
}

inline __device__ void makeYDiffusionChoices(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeYIndex, unsigned int windowIndex, const unsigned int windowYIndex, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(window, sitesWindow, choices, latticeIndex, windowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowYIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?latticeIndex-(latticeXSizeC*MPD_APRON_SIZE):latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)+(latticeXYSizeC);
        unsigned int apronWindowIndex = windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex, apronWindowIndex, apronWindowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowYIndex >= TUNE_MPD_Y_BLOCK_Y_SIZE)
    {
        unsigned int apronLatticeIndex = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?latticeIndex+(latticeXSizeC*MPD_APRON_SIZE):latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)-(latticeXYSizeC);
        unsigned int apronWindowIndex = windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex,  apronWindowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex, MPD_Y_WINDOW_SIZE, timestepHash);
    }
}

inline __device__ void makeZDiffusionChoices(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(window, sitesWindow, choices, latticeIndex, windowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowZIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE):latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)+(global_latticeXYZSizeC);
        unsigned int apronWindowIndex = windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex, apronWindowIndex, apronWindowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        unsigned int apronLatticeIndex = (latticeZIndex<global_latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE):latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)-(global_latticeXYZSizeC);
        unsigned int apronWindowIndex = windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex, MPD_Z_WINDOW_SIZE, timestepHash);

    }
}

/**
 * Copies a segment of a lattice from device memory to shared memory.
 *
 * NOTE: For multi-GPU execution, the latticeIndex is the local memory index.  It is the host's responsibility to translate it to a global index
 */
inline __device__ void performPropagation(unsigned int * __restrict__ lattice, const unsigned int * __restrict__ window, const unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int windowIndexMinus, const unsigned int windowIndex, const unsigned int windowIndexPlus, const unsigned int windowSize, unsigned int * __restrict__ siteOverflowList)
{
    int nextParticle=0;
    unsigned int newParticlesInt[MPD_WORDS_PER_SITE*3];
    unsigned char* newParticles = (unsigned char*)newParticlesInt;

    // Initialize the particles that will be copied back to the lattice to zero.
	for(int i=0; i<MPD_WORDS_PER_SITE; i++)
		newParticlesInt[i] = 0;


	// Pointers to fudge
	const unsigned char* particles;
	unsigned char* particleChoices;
	const unsigned int *win_ptr;
	const unsigned int *choice_ptr;

	// Propagate particles staying at here
	win_ptr = window;
	choice_ptr = choices;
    #pragma unroll
	for(int w=0; w < MPD_WORDS_PER_SITE; w++, win_ptr = win_ptr+windowSize, choice_ptr = choice_ptr+windowSize)
	{
		particles = (const unsigned char*)&win_ptr[windowIndex];
		particleChoices = (unsigned char*)&choice_ptr[windowIndex];
		for (int i=0; i<4; i++)
			if (particleChoices[i] == MPD_MOVE_STAY)
				newParticles[nextParticle++] = particles[i];
	}

    // Propagate particles moving in from the minus site.
	win_ptr = window;
	choice_ptr = choices;
    #pragma unroll
	for(int w=0; w < MPD_WORDS_PER_SITE; w++, win_ptr = win_ptr+windowSize, choice_ptr = choice_ptr+windowSize)
	{
		particles = (const unsigned char*)&win_ptr[windowIndexMinus];
		particleChoices = (unsigned char*)&choice_ptr[windowIndexMinus];
		for (int i=0; i<4; i++)
			if (particleChoices[i] == MPD_MOVE_PLUS)
				newParticles[nextParticle++] = particles[i];
	}


    // Propagate particles moving in from the plus site.
	win_ptr = window;
	choice_ptr = choices;
    #pragma unroll
	for(int w=0; w < MPD_WORDS_PER_SITE; w++, win_ptr = win_ptr+windowSize, choice_ptr = choice_ptr+windowSize)
	{
		particles = (const unsigned char*)&win_ptr[windowIndexPlus];
		particleChoices = (unsigned char*)&choice_ptr[windowIndexPlus];
		for (int i=0; i<4; i++)
			if (particleChoices[i] == MPD_MOVE_MINUS)
				newParticles[nextParticle++] = particles[i];
	}

	for(int w=0; w < MPD_WORDS_PER_SITE; w++, lattice += latticeXYZSizeC)
		lattice[latticeIndex] = newParticlesInt[w];

    // Move any leftover particles to the overflow list.
    for (int i=MPD_PARTICLES_PER_SITE; i<nextParticle; i++)
    {
        int exceptionIndex = atomicAdd(siteOverflowList, 1);
        if (exceptionIndex < TUNE_MPD_MAX_PARTICLE_OVERFLOWS)
        {
            siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
            siteOverflowList[(exceptionIndex*2)+2]=newParticles[i];
        }
    }
}
