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

#if !defined MPD_MAX_PARTICLE_OVERFLOWS
#error "Must define the maximum size of overflow particle list."
#endif

#define MPD_PARTICLES_PER_SITE  LS_WORDS_PER_SITE*4
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
#define MPD_MOVE_X_PLUS         2
#define MPD_MOVE_X_MINUS        3
#define MPD_MOVE_Y_PLUS         4
#define MPD_MOVE_Y_MINUS        5
#define MPD_MOVE_Z_PLUS         6
#define MPD_MOVE_Z_MINUS        7

#if __CUDA_ARCH__ >= 120
#define WARP_COND(x) __any(x)
#else
#define WARP_COND(x) (x)
#endif

__device__ float lookupTransitionProbability(const unsigned int particleType, const unsigned char sourceSite, const unsigned char destSite)
{
    if (destSite == 0xFF) return 0.0;
    return 0.25;
}

/**
 * Gets the random choices for a one-dimension diffusion move.
 *
 * @return  For each bit-packed each particle: 0 = no particle, 1 = stay in place, 2 = move in the positive direction,
 *          4 = move in the minus direction.
 */
__device__ void makeDiffusionChoices(const unsigned int * __restrict__ window, unsigned int * __restrict__ sharedChoices, const unsigned int latticeIndex, const unsigned int windowIndex, const unsigned long long timestepHash)
{
    // Create the choices.
    unsigned int choicesInt[LS_WORDS_PER_SITE];
    unsigned char * choices = (unsigned char *)choicesInt;

    // Get the particles.
    unsigned int particlesInt[LS_WORDS_PER_SITE];
    unsigned char * particles = (unsigned char *)particlesInt;
    particlesInt[0] = window[windowIndex];
    #if LS_WORDS_PER_SITE >= 2
    particlesInt[1] = window[windowIndex+LS_XYZ_WINDOW_SIZE];
    #endif

    // If there are no particles, we are done.
    if (WARP_COND(particles[0]))
    {
        const unsigned char site = 0; //(particles&MPD_SITE_MASK)>>MPD_SITE_SHIFT;
        const unsigned char siteMinus = 0; //(window[windowIndexMinus]&MPD_SITE_MASK)>>MPD_SITE_SHIFT;
        const unsigned char sitePlus = 0; //(window[windowIndexPlus]&MPD_SITE_MASK)>>MPD_SITE_SHIFT;

        // Make choices for all of the particles in the site.
        for (int i=0; i<4; i++)
        {
            // Get the probability of moving plus and minus.
            float probMinus=lookupTransitionProbability(particles[i], site, siteMinus);
            float probPlus=lookupTransitionProbability(particles[i], site, sitePlus);

            // Get the random value and see which direction we should move in.
            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, i, timestepHash);
            choices[i] = (particles[i] > 0)?(MPD_MOVE_STAY):(MPD_MOVE_NONE);
            choices[i] = (randomValue < probMinus)?(MPD_MOVE_X_MINUS):(MPD_MOVE_STAY);
            choices[i] = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_X_PLUS):(MPD_MOVE_STAY);
        }
        if (WARP_COND(particles[4]))
        {
            for (int i=4; i<MPD_PARTICLES_PER_SITE; i++)
            {
                // Get the probability of moving plus and minus.
                float probMinus=lookupTransitionProbability(particles[i], site, siteMinus);
                float probPlus=lookupTransitionProbability(particles[i], site, sitePlus);

                // Get the random value and see which direction we should move in.
                float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, i, timestepHash);
                choices[i] = (particles[i] > 0)?(MPD_MOVE_STAY):(MPD_MOVE_NONE);
                choices[i] = (randomValue < probMinus)?(MPD_MOVE_X_MINUS):(MPD_MOVE_STAY);
                choices[i] = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_X_PLUS):(MPD_MOVE_STAY);
            }
        }
    }

    // Save the choices.
    sharedChoices[windowIndex] = choicesInt[0];
    #if LS_WORDS_PER_SITE >= 2
    sharedChoices[windowIndex+LS_XYZ_WINDOW_SIZE] = choicesInt[1];
    #endif
}

/**
 * Copies a segment of a lattice from device memory to shared memory.
 */
__device__ void performPropagation(unsigned int * __restrict__ lattice, const unsigned int * __restrict__ window, const unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned int windowIndex, unsigned int * __restrict__ siteOverflowList)
{
    //Create any needed pointer aliases.
    #if LS_WORDS_PER_SITE >= 2
    unsigned int * lattice2 = lattice+latticeXYZSize;
    const unsigned int * window2 = window+LS_XYZ_WINDOW_SIZE;
    const unsigned int * choices2 = choices+LS_XYZ_WINDOW_SIZE;
    #endif

    int nextParticle=0;
    #if __CUDA_ARCH__ >= 200
    unsigned int newParticlesInt[LS_WORDS_PER_SITE*7];
    #else
    __shared__ unsigned int newParticlesInt[LS_WORDS_PER_SITE*7];
    #endif
    unsigned char* newParticles = (unsigned char*)newParticlesInt;

    // Initialize the particles that will be copied back to the lattice to zero.
    newParticlesInt[0] = 0;
    #if LS_WORDS_PER_SITE >= 2
    newParticlesInt[1] = 0;
    #endif

    // Propagate particles staying at this site.
    unsigned int checkingWindowIndex = windowIndex;
    const unsigned char * particles = (const unsigned char*)&window[checkingWindowIndex];
    unsigned char * particleChoices = (unsigned char*)&choices[checkingWindowIndex];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_STAY) newParticles[nextParticle++] = particles[i];
    #if LS_WORDS_PER_SITE >= 2
    if (WARP_COND(particles[3]))
    {
        particles = (const unsigned char*)&window2[checkingWindowIndex];
        particleChoices = (unsigned char*)&choices2[checkingWindowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_STAY) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving into this site from the minus x direction.
    checkingWindowIndex = windowIndex-1;
    particles = (const unsigned char*)&window[checkingWindowIndex];
    particleChoices = (unsigned char*)&choices[checkingWindowIndex];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_X_PLUS) newParticles[nextParticle++] = particles[i];
    #if LS_WORDS_PER_SITE >= 2
    if (WARP_COND(particles[3]))
    {
        particles = (const unsigned char*)&window2[checkingWindowIndex];
        particleChoices = (unsigned char*)&choices2[checkingWindowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_X_PLUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving into this site from the plus x direction.
    checkingWindowIndex = windowIndex+1;
    particles = (const unsigned char*)&window[checkingWindowIndex];
    particleChoices = (unsigned char*)&choices[checkingWindowIndex];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_X_MINUS) newParticles[nextParticle++] = particles[i];
    #if LS_WORDS_PER_SITE >= 2
    if (WARP_COND(particles[3]))
    {
        particles = (const unsigned char*)&window2[checkingWindowIndex];
        particleChoices = (unsigned char*)&choices2[checkingWindowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_X_MINUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving into this site from the minus y direction.
    checkingWindowIndex = windowIndex-LS_XYZ_WINDOW_X_SIZE;
    particles = (const unsigned char*)&window[checkingWindowIndex];
    particleChoices = (unsigned char*)&choices[checkingWindowIndex];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Y_PLUS) newParticles[nextParticle++] = particles[i];
    #if LS_WORDS_PER_SITE >= 2
    if (WARP_COND(particles[3]))
    {
        particles = (const unsigned char*)&window2[checkingWindowIndex];
        particleChoices = (unsigned char*)&choices2[checkingWindowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Y_PLUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving into this site from the plus y direction.
    checkingWindowIndex = windowIndex+LS_XYZ_WINDOW_X_SIZE;
    particles = (const unsigned char*)&window[checkingWindowIndex];
    particleChoices = (unsigned char*)&choices[checkingWindowIndex];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Y_MINUS) newParticles[nextParticle++] = particles[i];
    #if LS_WORDS_PER_SITE >= 2
    if (WARP_COND(particles[3]))
    {
        particles = (const unsigned char*)&window2[checkingWindowIndex];
        particleChoices = (unsigned char*)&choices2[checkingWindowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Y_MINUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving into this site from the minus z direction.
    checkingWindowIndex = windowIndex-LS_XYZ_WINDOW_X_SIZE*LS_XYZ_WINDOW_Y_SIZE;
    particles = (const unsigned char*)&window[checkingWindowIndex];
    particleChoices = (unsigned char*)&choices[checkingWindowIndex];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Z_PLUS) newParticles[nextParticle++] = particles[i];
    #if LS_WORDS_PER_SITE >= 2
    if (WARP_COND(particles[3]))
    {
        particles = (const unsigned char*)&window2[checkingWindowIndex];
        particleChoices = (unsigned char*)&choices2[checkingWindowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Z_PLUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Propagate particles moving into this site from the plus z direction.
    checkingWindowIndex = windowIndex+LS_XYZ_WINDOW_X_SIZE*LS_XYZ_WINDOW_Y_SIZE;
    particles = (const unsigned char*)&window[checkingWindowIndex];
    particleChoices = (unsigned char*)&choices[checkingWindowIndex];
    #pragma unroll
    for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Z_MINUS) newParticles[nextParticle++] = particles[i];
    #if LS_WORDS_PER_SITE >= 2
    if (WARP_COND(particles[3]))
    {
        particles = (const unsigned char*)&window2[checkingWindowIndex];
        particleChoices = (unsigned char*)&choices2[checkingWindowIndex];
        #pragma unroll
        for (int i=0; i<4; i++) if (particleChoices[i] == MPD_MOVE_Z_MINUS) newParticles[nextParticle++] = particles[i];
    }
    #endif

    // Copy the new particles to the lattice.
    lattice[latticeIndex] = newParticlesInt[0];
    #if LS_WORDS_PER_SITE >= 2
    lattice2[latticeIndex] = newParticlesInt[1];
    #endif

    // Move any leftover particles to the overflow list.
    for (int i=MPD_PARTICLES_PER_SITE; i<nextParticle; i++)
    {
        int exceptionIndex = atomicAdd(siteOverflowList, 1);
        if (exceptionIndex < MPD_MAX_PARTICLE_OVERFLOWS)
        {
            siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
            siteOverflowList[(exceptionIndex*2)+2]=newParticles[i];
        }
    }
}

