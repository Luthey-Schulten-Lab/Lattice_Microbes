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
    return TC[sourceSite*numberSiteTypesC*numberSpeciesC + destSite*numberSpeciesC + (particleType-1)];
}

/**
 * Gets the random choices for a one-dimension diffusion move.
 *
 * @return  For each bit-packed each particle: 0 = no particle, 1 = stay in place, 2 = move in the positive direction,
 *          4 = move in the minus direction.
 */
inline __device__ void makeDiffusionChoices(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, unsigned int * __restrict__ choices, const unsigned int latticeIndexG, const unsigned int windowIndexMinus, const unsigned int windowIndex, const unsigned int windowIndexPlus, const unsigned int windowSize, const unsigned long long timestepHash)
{
    #if MPD_WORDS_PER_SITE >= 2
    const unsigned int* window2 = window+windowSize;
    unsigned int* choices2 = choices+windowSize;
    #endif

    // Set the default choice to stay.
    choices[windowIndex]=0;
    #if MPD_WORDS_PER_SITE >= 2
    choices2[windowIndex]=0;
    #endif

    // If there are no particles, we are done.
    if (window[windowIndex] > 0)
    {
        // Get the particles as a byte pointer.
        const unsigned char* particles = (const unsigned char*)&window[windowIndex];
        unsigned char* particleChoices = (unsigned char*)&choices[windowIndex];

        const unsigned char site = sitesWindow[windowIndex];
        const unsigned char siteMinus = sitesWindow[windowIndexMinus];
        const unsigned char sitePlus = sitesWindow[windowIndexPlus];

        // Make choices for all of the particles in the site.
        for (int i=0; i<4; i++)
        {
            if (particles[i] > 0)
            {
                // Get the probability of moving plus and minus.
                float probMinus=lookupTransitionProbability(particles[i], site, siteMinus);
                float probPlus=lookupTransitionProbability(particles[i], site, sitePlus);

                // Get the random value and see which direction we should move in.
                float randomValue = getRandomHashFloat(latticeIndexG, MPD_PARTICLE_COUNT_BITS, i, timestepHash);
                particleChoices[i] = (randomValue < probMinus)?(MPD_MOVE_MINUS):(MPD_MOVE_STAY);
                particleChoices[i] = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_PLUS):(particleChoices[i]);
            }
        }

        // Make choices for all of the particles in the second site.
        #if MPD_WORDS_PER_SITE >= 2
        if (window2[windowIndex] > 0)
        {
            particles = (const unsigned char*)&window2[windowIndex];
            particleChoices = (unsigned char*)&choices2[windowIndex];
            for (int i=0; i<4; i++)
            {
                if (particles[i] > 0)
                {
                    // Get the probability of moving plus and minus.
                    float probMinus=lookupTransitionProbability(particles[i], site, siteMinus);
                    float probPlus=lookupTransitionProbability(particles[i], site, sitePlus);

                    // Get the random value and see which direction we should move in.
                    float randomValue = getRandomHashFloat(latticeIndexG, MPD_PARTICLE_COUNT_BITS, i+4, timestepHash);
                    particleChoices[i] = (randomValue < probMinus)?(MPD_MOVE_MINUS):(MPD_MOVE_STAY);
                    particleChoices[i] = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_PLUS):(particleChoices[i]);
                }
            }
        }
        #endif
    }
}

/**
 * Makes the X diffusion choices for all particles in loaded sites further than one from the edge.
 */
inline __device__ void makeXDiffusionChoices(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex, const unsigned int blockXSize, const unsigned long long timestepHash, const unsigned int globalLatticeIndex)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(window, sitesWindow, choices, globalLatticeIndex, windowIndex-1, windowIndex, windowIndex+1, MPD_X_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowIndex >= blockXSize)
    {
        unsigned int apronLatticeIndex = (latticeXIndex>=blockXSize)?latticeIndex-blockXSize:latticeIndex+(latticeXSizeC-blockXSize);
        unsigned int apronWindowIndex = windowIndex-blockXSize;
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex-1, apronWindowIndex, apronWindowIndex+1, MPD_X_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeXIndex<latticeXSizeC-blockXSize)?latticeIndex+blockXSize:latticeIndex-latticeXSizeC+blockXSize;
        unsigned int apronWindowIndex = windowIndex+blockXSize;
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex-1, apronWindowIndex, apronWindowIndex+1, MPD_X_WINDOW_SIZE, timestepHash);
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
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowYIndex >= TUNE_MPD_Y_BLOCK_Y_SIZE)
    {
        unsigned int apronLatticeIndex = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?latticeIndex+(latticeXSizeC*MPD_APRON_SIZE):latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)-(latticeXYSizeC);
        unsigned int apronWindowIndex = windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex,  apronWindowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, timestepHash);
    }
}

inline __device__ void makeZDiffusionChoices(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(window, sitesWindow, choices, latticeIndex, windowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowZIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE):latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)+(latticeXYZSizeC);
        unsigned int apronWindowIndex = windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        unsigned int apronLatticeIndex = (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE):latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)-(latticeXYZSizeC);
        unsigned int apronWindowIndex = windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoices(window, sitesWindow, choices, apronLatticeIndex, apronWindowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, timestepHash);

    }
}

/**
 * Copies a segment of a lattice from device memory to shared memory.
 */
inline __device__ void performPropagation(unsigned int * __restrict__ lattice, const unsigned int * __restrict__ window, const unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int windowIndexMinus, const unsigned int windowIndex, const unsigned int windowIndexPlus, const unsigned int windowSize, unsigned int * __restrict__ siteOverflowList)
{
    //Create any needed pointer aliases.
    #if MPD_WORDS_PER_SITE >= 2
    unsigned int * lattice2 = lattice+latticeXYZSizeC;
    const unsigned int * window2 = window+windowSize;
    const unsigned int * choices2 = choices+windowSize;
    #endif

    int nextParticle=0;
    unsigned int newParticlesInt[MPD_WORDS_PER_SITE*3];
    unsigned char* newParticles = (unsigned char*)newParticlesInt;

    // Initialize the particles that will be copied back to the lattice to zero.
    newParticlesInt[0] = 0;
    #if MPD_WORDS_PER_SITE >= 2
    newParticlesInt[1] = 0;
    #endif

    // Propagate particles staying at this site.
    const unsigned char* particles = (const unsigned char*)&window[windowIndex];
    unsigned char* particleChoices = (unsigned char*)&choices[windowIndex];

    //DEBUG
    //for (int i=0; i<4; i++) newParticles[nextParticle++] = (byte)(100.0f*TC[0]);
    //particles = (const unsigned char*)&window2[windowIndex];
    //particleChoices = (unsigned char*)&choices2[windowIndex];
    //#pragma unroll
    //for (int i=0; i<4; i++) newParticles[nextParticle++] = (byte)(100.0f*TC[1]);
    //DEBUG

    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_STAY) newParticles[nextParticle++] = particles[i];
    #if MPD_WORDS_PER_SITE >= 2
    if (particles[3])
    {
        particles = (const unsigned char*)&window2[windowIndex];
        particleChoices = (unsigned char*)&choices2[windowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_STAY) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving in from the minus site.
    particles = (const unsigned char*)&window[windowIndexMinus];
    particleChoices = (unsigned char*)&choices[windowIndexMinus];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_PLUS) newParticles[nextParticle++] = particles[i];
    #if MPD_WORDS_PER_SITE >= 2
    if (particles[3])
    {
        particles = (const unsigned char*)&window2[windowIndexMinus];
        particleChoices = (unsigned char*)&choices2[windowIndexMinus];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_PLUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving in from the plus site.
    particles = (const unsigned char*)&window[windowIndexPlus];
    particleChoices = (unsigned char*)&choices[windowIndexPlus];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_MINUS) newParticles[nextParticle++] = particles[i];
    #if MPD_WORDS_PER_SITE >= 2
    if (particles[3])
    {
        particles = (const unsigned char*)&window2[windowIndexPlus];
        particleChoices = (unsigned char*)&choices2[windowIndexPlus];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_MINUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Copy the new particles to the next lattice.
    lattice[latticeIndex] = newParticlesInt[0];
    #if MPD_WORDS_PER_SITE >= 2
    lattice2[latticeIndex] = newParticlesInt[1];
    #endif

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
