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

#define MPD_PARTICLE_COUNT      7
#define MPD_PARTICLES_MASK      0x0FFFFFFF
#define MPD_PARTICLE_MASK       0x0000000F
#define MPD_PARTICLE_1_MASK     0x0000000F
#define MPD_PARTICLE_2_MASK     0x000000F0
#define MPD_PARTICLE_3_MASK     0x00000F00
#define MPD_PARTICLE_4_MASK     0x0000F000
#define MPD_PARTICLE_5_MASK     0x000F0000
#define MPD_PARTICLE_6_MASK     0x00F00000
#define MPD_PARTICLE_7_MASK     0x0F000000
#define MPD_PARTICLE_SHIFT      4
#define MPD_PARTICLE_1_SHIFT    0
#define MPD_PARTICLE_2_SHIFT    4
#define MPD_PARTICLE_3_SHIFT    8
#define MPD_PARTICLE_4_SHIFT    12
#define MPD_PARTICLE_5_SHIFT    16
#define MPD_PARTICLE_6_SHIFT    20
#define MPD_PARTICLE_7_SHIFT    24
#define MPD_PARTICLE_MAX_SHIFT  MPD_PARTICLE_7_SHIFT

#define MPD_SITE_MASK           0x70000000
#define MPD_SITE_SHIFT          28
#define MPD_SITE_MAX            0x00000006
#define MPD_SITE_BOUNDARY       0x00000007
#define MPD_SITE_COUNT          8

#define MPD_X_BLOCK_MAX_X_SIZE          128

#define MPD_Y_BLOCK_X_SIZE              16
#define MPD_Y_BLOCK_X_SIZE_SHIFT_MULT   4
#define MPD_Y_BLOCK_Y_SIZE              8
#define MPD_Y_BLOCK_Y_SIZE_SHIFT_MULT   3

#define MPD_Z_BLOCK_X_SIZE              16
#define MPD_Z_BLOCK_X_SIZE_SHIFT_MULT   4
#define MPD_Z_BLOCK_Z_SIZE              8
#define MPD_Z_BLOCK_Z_SIZE_SHIFT_MULT   3

#define MPD_WINDOW_APRON_SIZE           2
#define MPD_X_WINDOW_SIZE               (MPD_WINDOW_APRON_SIZE+MPD_X_BLOCK_MAX_X_SIZE+MPD_WINDOW_APRON_SIZE)
#define MPD_Y_WINDOW_SIZE               (MPD_Y_BLOCK_X_SIZE*(MPD_WINDOW_APRON_SIZE+MPD_Y_BLOCK_Y_SIZE+MPD_WINDOW_APRON_SIZE))
#define MPD_Z_WINDOW_SIZE               (MPD_Z_BLOCK_X_SIZE*(MPD_WINDOW_APRON_SIZE+MPD_Z_BLOCK_Z_SIZE+MPD_WINDOW_APRON_SIZE))

#define MPD_PARTICLE_COUNT_BITS 3
#define MPD_PARTICLE_1_STAY     0x00000001
#define MPD_PARTICLE_1_PLUS     0x00000002
#define MPD_PARTICLE_1_MINUS    0x00000004
#define MPD_PARTICLE_2_STAY     0x00000010
#define MPD_PARTICLE_2_PLUS     0x00000020
#define MPD_PARTICLE_2_MINUS    0x00000040
#define MPD_PARTICLE_3_STAY     0x00000100
#define MPD_PARTICLE_3_PLUS     0x00000200
#define MPD_PARTICLE_3_MINUS    0x00000400
#define MPD_PARTICLE_4_STAY     0x00001000
#define MPD_PARTICLE_4_PLUS     0x00002000
#define MPD_PARTICLE_4_MINUS    0x00004000
#define MPD_PARTICLE_5_STAY     0x00010000
#define MPD_PARTICLE_5_PLUS     0x00020000
#define MPD_PARTICLE_5_MINUS    0x00040000
#define MPD_PARTICLE_6_STAY     0x00100000
#define MPD_PARTICLE_6_PLUS     0x00200000
#define MPD_PARTICLE_6_MINUS    0x00400000
#define MPD_PARTICLE_7_STAY     0x01000000
#define MPD_PARTICLE_7_PLUS     0x02000000
#define MPD_PARTICLE_7_MINUS    0x04000000

#define LKCUDA_OVERFLOW_LIST_ENTRIES                 1+2*LKCUDA_MAX_PARTICLE_OVERFLOWS


__device__ float lookupTransitionProbability(const unsigned int particleType, const unsigned char sourceSite, const unsigned char destSite)
{
    return 0.25;
}

/**
 * Gets the random choices for a one-dimension diffusion move.
 *
 * @return  For each bit-packed each particle: 0 = no particle, 1 = stay in place, 2 = move in the positive direction,
 *          4 = move in the minus direction.
 */
__device__ void makeDiffusionChoices(const unsigned int latticeIndex, const unsigned int* latticeSegment, const unsigned int latticeSegmentIndexMinus, const unsigned int latticeSegmentIndex, const unsigned int latticeSegmentIndexPlus, unsigned int* latticeSegmentChoices, const unsigned long long timestepHash)
{
    // Set the initial choices to be no particle.
    unsigned int choices = 0x00000000;

    // Get the particles.
    unsigned int particles = latticeSegment[latticeSegmentIndex];

    // If we don't have any particle, don't do anything.
    if (particles&MPD_PARTICLES_MASK)
    {
        unsigned char site = (particles&MPD_SITE_MASK)>>MPD_SITE_SHIFT;
        unsigned char sitePlus = (latticeSegment[latticeSegmentIndexPlus]&MPD_SITE_MASK)>>MPD_SITE_SHIFT;
        unsigned char siteMinus = (latticeSegment[latticeSegmentIndexMinus]&MPD_SITE_MASK)>>MPD_SITE_SHIFT;

        // Make the choice for the first particle.
        if (particles&MPD_PARTICLE_1_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_1_MASK);
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            // Get the random value.
            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 0, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_1_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_1_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_1_STAY):(0);
        }

        // Make the choice for the second particle.
        if (particles&MPD_PARTICLE_2_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_2_MASK)>>MPD_PARTICLE_2_SHIFT;
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 1, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_2_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_2_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_2_STAY):(0);
        }

#if MPD_PARTICLE_COUNT >= 3
        // Make the choice for the third particle.
        if (particles&MPD_PARTICLE_3_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_3_MASK)>>MPD_PARTICLE_3_SHIFT;
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 2, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_3_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_3_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_3_STAY):(0);
        }
#endif
#if MPD_PARTICLE_COUNT >= 4
        // Make the choice for the fourth particle.
        if (particles&MPD_PARTICLE_4_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_4_MASK)>>MPD_PARTICLE_4_SHIFT;
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 3, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_4_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_4_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_4_STAY):(0);
        }
#endif
#if MPD_PARTICLE_COUNT >= 5
        // Make the choice for the fourth particle.
        if (particles&MPD_PARTICLE_5_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_5_MASK)>>MPD_PARTICLE_5_SHIFT;
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 4, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_5_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_5_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_5_STAY):(0);
        }
#endif
#if MPD_PARTICLE_COUNT >= 6
        // Make the choice for the fourth particle.
        if (particles&MPD_PARTICLE_6_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_6_MASK)>>MPD_PARTICLE_6_SHIFT;
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 5, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_6_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_6_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_6_STAY):(0);
        }
#endif
#if MPD_PARTICLE_COUNT >= 7
        // Make the choice for the fourth particle.
        if (particles&MPD_PARTICLE_7_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_7_MASK)>>MPD_PARTICLE_7_SHIFT;
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 6, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_7_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_7_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_7_STAY):(0);
        }
#endif
#if MPD_PARTICLE_COUNT >= 8
        // Make the choice for the fourth particle.
        if (particles&MPD_PARTICLE_8_MASK)
        {
            unsigned int particle=(particles&MPD_PARTICLE_8_MASK)>>MPD_PARTICLE_8_SHIFT;
            float probPlus=lookupTransitionProbability(particle, site, sitePlus);
            float probMinus=probPlus+lookupTransitionProbability(particle, site, siteMinus);

            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 7, timestepHash);
            choices |= (randomValue <= probPlus)?(MPD_PARTICLE_8_PLUS):(0);
            choices |= (randomValue > probPlus && randomValue <= probMinus)?(MPD_PARTICLE_8_MINUS):(0);
            choices |= (randomValue > probMinus)?(MPD_PARTICLE_8_STAY):(0);
        }
#endif
    }

    latticeSegmentChoices[latticeSegmentIndex] = choices;
}

/**
 * Copies a segment of a lattice from device memory to shared memory.
 */
__device__ void performPropagation(unsigned int* outLattice, const unsigned int latticeIndex, const unsigned int* latticeSegment, const unsigned int latticeSegmentIndexMinus, const unsigned int latticeSegmentIndex, const unsigned int latticeSegmentIndexPlus, const unsigned int* choices, unsigned int* siteOverflowList)
{
    // Get the particles and their chocies.
    unsigned int particles = latticeSegment[latticeSegmentIndex];
    unsigned int particlesChoices = choices[latticeSegmentIndex];

    // Preserve the site type, but none of the particles.
    unsigned int newParticles = particles&MPD_SITE_MASK;
    unsigned int newParticleShift=0;

    // Find particles that are not moving.
    if (particlesChoices&MPD_PARTICLE_1_STAY)
    {
        newParticles |=  (particles&MPD_PARTICLE_1_MASK)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
    if (particlesChoices&MPD_PARTICLE_2_STAY)
    {
        newParticles |=  ((particles&MPD_PARTICLE_2_MASK)>>MPD_PARTICLE_2_SHIFT)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
#if MPD_PARTICLE_COUNT >= 3
    if (particlesChoices&MPD_PARTICLE_3_STAY)
    {
        newParticles |=  ((particles&MPD_PARTICLE_3_MASK)>>MPD_PARTICLE_3_SHIFT)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
#endif
#if MPD_PARTICLE_COUNT >= 4

    if (particlesChoices&MPD_PARTICLE_4_STAY)
    {
        newParticles |=  ((particles&MPD_PARTICLE_4_MASK)>>MPD_PARTICLE_4_SHIFT)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
#endif
#if MPD_PARTICLE_COUNT >= 5

    if (particlesChoices&MPD_PARTICLE_5_STAY)
    {
        newParticles |=  ((particles&MPD_PARTICLE_5_MASK)>>MPD_PARTICLE_5_SHIFT)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
#endif
#if MPD_PARTICLE_COUNT >= 6

    if (particlesChoices&MPD_PARTICLE_6_STAY)
    {
        newParticles |=  ((particles&MPD_PARTICLE_6_MASK)>>MPD_PARTICLE_6_SHIFT)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
#endif
#if MPD_PARTICLE_COUNT >= 7

    if (particlesChoices&MPD_PARTICLE_7_STAY)
    {
        newParticles |=  ((particles&MPD_PARTICLE_7_MASK)>>MPD_PARTICLE_7_SHIFT)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
#endif
#if MPD_PARTICLE_COUNT >= 8

    if (particlesChoices&MPD_PARTICLE_8_STAY)
    {
        newParticles |=  ((particles&MPD_PARTICLE_8_MASK)>>MPD_PARTICLE_8_SHIFT)<<newParticleShift;
        newParticleShift += MPD_PARTICLE_SHIFT;
    }
#endif

    // Find the particles that are moving plus from the site to the minus.
    particles = latticeSegment[latticeSegmentIndexMinus];
    particlesChoices = choices[latticeSegmentIndexMinus];
    if (particlesChoices&MPD_PARTICLE_1_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  (particles&MPD_PARTICLE_1_MASK)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=particles&MPD_PARTICLE_1_MASK;
            }
        }
    }
    if (particlesChoices&MPD_PARTICLE_2_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_2_MASK)>>MPD_PARTICLE_2_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_2_MASK)>>MPD_PARTICLE_2_SHIFT;
            }
        }
    }
#if MPD_PARTICLE_COUNT >= 3
    if (particlesChoices&MPD_PARTICLE_3_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_3_MASK)>>MPD_PARTICLE_3_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_3_MASK)>>MPD_PARTICLE_3_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 4
    if (particlesChoices&MPD_PARTICLE_4_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_4_MASK)>>MPD_PARTICLE_4_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_4_MASK)>>MPD_PARTICLE_4_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 5
    if (particlesChoices&MPD_PARTICLE_5_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_5_MASK)>>MPD_PARTICLE_5_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_5_MASK)>>MPD_PARTICLE_5_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 6
    if (particlesChoices&MPD_PARTICLE_6_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_6_MASK)>>MPD_PARTICLE_6_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_6_MASK)>>MPD_PARTICLE_6_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 7
    if (particlesChoices&MPD_PARTICLE_7_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_7_MASK)>>MPD_PARTICLE_7_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_7_MASK)>>MPD_PARTICLE_7_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 8
    if (particlesChoices&MPD_PARTICLE_8_PLUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_8_MASK)>>MPD_PARTICLE_8_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_8_MASK)>>MPD_PARTICLE_8_SHIFT;
            }
        }
    }
#endif

    // Find the particles that are moving minus from the site to the plus.
    particles = latticeSegment[latticeSegmentIndexPlus];
    particlesChoices = choices[latticeSegmentIndexPlus];
    if (particlesChoices&MPD_PARTICLE_1_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  (particles&MPD_PARTICLE_1_MASK)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=particles&MPD_PARTICLE_1_MASK;
            }
        }
    }
    if (particlesChoices&MPD_PARTICLE_2_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_2_MASK)>>MPD_PARTICLE_2_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_2_MASK)>>MPD_PARTICLE_2_SHIFT;
            }
        }
    }
#if MPD_PARTICLE_COUNT >= 3
    if (particlesChoices&MPD_PARTICLE_3_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_3_MASK)>>MPD_PARTICLE_3_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_3_MASK)>>MPD_PARTICLE_3_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 4
    if (particlesChoices&MPD_PARTICLE_4_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_4_MASK)>>MPD_PARTICLE_4_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_4_MASK)>>MPD_PARTICLE_4_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 5
    if (particlesChoices&MPD_PARTICLE_5_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_5_MASK)>>MPD_PARTICLE_5_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_5_MASK)>>MPD_PARTICLE_5_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 6
    if (particlesChoices&MPD_PARTICLE_6_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_6_MASK)>>MPD_PARTICLE_6_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_6_MASK)>>MPD_PARTICLE_6_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 7
    if (particlesChoices&MPD_PARTICLE_7_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_7_MASK)>>MPD_PARTICLE_7_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_7_MASK)>>MPD_PARTICLE_7_SHIFT;
            }
        }
    }
#endif
#if MPD_PARTICLE_COUNT >= 8
    if (particlesChoices&MPD_PARTICLE_8_MINUS)
    {
        if (newParticleShift <= MPD_PARTICLE_MAX_SHIFT)
        {
            newParticles |=  ((particles&MPD_PARTICLE_8_MASK)>>MPD_PARTICLE_8_SHIFT)<<newParticleShift;
            newParticleShift += MPD_PARTICLE_SHIFT;
        }
        else
        {
            int exceptionIndex = atomicAdd(siteOverflowList, 1);
            if (exceptionIndex < LKCUDA_MAX_PARTICLE_OVERFLOWS)
            {
                siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                siteOverflowList[(exceptionIndex*2)+2]=(particles&MPD_PARTICLE_8_MASK)>>MPD_PARTICLE_8_SHIFT;
            }
        }
    }
#endif

    // Set the new site.
    outLattice[latticeIndex] = newParticles;

    // Wait for the whole region to be copied.
    __syncthreads();
}

/**
 * Copies a segment of a lattice from device memory to shared memory.
 */
__device__ void copyXWindowFromLattice(const unsigned int* inLattice, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int latticeXSize, unsigned int* latticeSegment, const unsigned int latticeSegmentIndex, const unsigned int blockXSize, const unsigned long long timestepHash)
{
    // Load the block.
    latticeSegment[latticeSegmentIndex] = inLattice[latticeIndex];

    // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
    if (latticeSegmentIndex >= blockXSize) latticeSegment[latticeSegmentIndex-blockXSize] = (latticeXIndex>=blockXSize)?inLattice[latticeIndex-blockXSize]:inLattice[latticeIndex-blockXSize+latticeXSize];

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
    if (latticeSegmentIndex < 2*MPD_WINDOW_APRON_SIZE) latticeSegment[latticeSegmentIndex+blockXSize] = (latticeXIndex<latticeXSize-blockXSize)?inLattice[latticeIndex+blockXSize]:inLattice[latticeIndex+blockXSize-latticeXSize];

    // Wait for the whole region to be copied.
    __syncthreads();
}

/**
 * Makes the X diffusion choices for all particles in loaded sites further than one from the edge.
 */
__device__ void makeXDiffusionChoices(const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int latticeXSize, const unsigned int* latticeSegment, const unsigned int latticeSegmentIndex, unsigned int* choices, const unsigned int blockXSize, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(latticeIndex, latticeSegment, latticeSegmentIndex-1, latticeSegmentIndex, latticeSegmentIndex+1, choices, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (latticeSegmentIndex >= (blockXSize+1))
    {
        unsigned int apronLatticeIndex = (latticeXIndex>=blockXSize)?latticeIndex-blockXSize:latticeIndex-blockXSize+latticeXSize;
        unsigned int apronLatticeSegmentIndex = latticeSegmentIndex-blockXSize;
        makeDiffusionChoices(apronLatticeIndex, latticeSegment, apronLatticeSegmentIndex-1, apronLatticeSegmentIndex, apronLatticeSegmentIndex+1, choices, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (latticeSegmentIndex < (2*MPD_WINDOW_APRON_SIZE)-1)
    {
        unsigned int apronLatticeIndex = (latticeXIndex<latticeXSize-blockXSize)?latticeIndex+blockXSize:latticeIndex+blockXSize-latticeXSize;
        unsigned int apronLatticeSegmentIndex = latticeSegmentIndex+blockXSize;
        makeDiffusionChoices(apronLatticeIndex, latticeSegment, apronLatticeSegmentIndex-1, apronLatticeSegmentIndex, apronLatticeSegmentIndex+1, choices, timestepHash);
    }

    // Wait for the whole segment to be processed.
    __syncthreads();
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
__device__ void copyYWindowFromLattice(const unsigned int* inLattice, const unsigned int latticeIndex, const unsigned int latticeYIndex, const unsigned int latticeXSize, const unsigned int latticeYSize, unsigned int* latticeSegment, const unsigned int latticeSegmentIndex, const unsigned int latticeSegmentYIndex, const unsigned long long timestepHash)
{
    // Load the block.
    latticeSegment[latticeSegmentIndex] = inLattice[latticeIndex];

    // If this thread is one that needs to load part of the leading apron, load it (unless this is in the first block).
    if (latticeSegmentYIndex < 2*MPD_WINDOW_APRON_SIZE) latticeSegment[latticeSegmentIndex-(MPD_Y_BLOCK_X_SIZE*MPD_WINDOW_APRON_SIZE)] = (latticeYIndex>=MPD_Y_BLOCK_Y_SIZE)?inLattice[latticeIndex-(latticeXSize*MPD_WINDOW_APRON_SIZE)]:inLattice[latticeIndex-(latticeXSize*MPD_WINDOW_APRON_SIZE)+(latticeXSize*latticeYSize)];

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is in the last block).
    if (latticeSegmentYIndex >= MPD_Y_BLOCK_Y_SIZE) latticeSegment[latticeSegmentIndex+(MPD_Y_BLOCK_X_SIZE*MPD_WINDOW_APRON_SIZE)] = (latticeYIndex<latticeYSize-MPD_Y_BLOCK_Y_SIZE)?inLattice[latticeIndex+(latticeXSize*MPD_WINDOW_APRON_SIZE)]:inLattice[latticeIndex+(latticeXSize*MPD_WINDOW_APRON_SIZE)-(latticeXSize*latticeYSize)];

    // Wait for the whole region to be copied.
    __syncthreads();
}

__device__ void makeYDiffusionChoices(const unsigned int latticeIndex, const unsigned int latticeYIndex, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int* latticeSegment, unsigned int latticeSegmentIndex, const unsigned int latticeSegmentYIndex, unsigned int* choices, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(latticeIndex, latticeSegment, latticeSegmentIndex-MPD_Y_BLOCK_X_SIZE, latticeSegmentIndex, latticeSegmentIndex+MPD_Y_BLOCK_X_SIZE, choices, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (latticeSegmentYIndex < (2*MPD_WINDOW_APRON_SIZE)-1)
    {
        unsigned int apronLatticeIndex = (latticeYIndex>=MPD_Y_BLOCK_Y_SIZE)?latticeIndex-(latticeXSize*(MPD_WINDOW_APRON_SIZE-1)):latticeIndex-(latticeXSize*(MPD_WINDOW_APRON_SIZE-1))+(latticeXSize*latticeYSize);
        unsigned int apronLatticeSegmentIndex = latticeSegmentIndex-(MPD_Y_BLOCK_X_SIZE*(MPD_WINDOW_APRON_SIZE-1));
        makeDiffusionChoices(apronLatticeIndex, latticeSegment, apronLatticeSegmentIndex-MPD_Y_BLOCK_X_SIZE, apronLatticeSegmentIndex, apronLatticeSegmentIndex+MPD_Y_BLOCK_X_SIZE, choices, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (latticeSegmentYIndex >= MPD_Y_BLOCK_Y_SIZE+1)
    {
        unsigned int apronLatticeIndex = (latticeYIndex<latticeYSize-MPD_Y_BLOCK_Y_SIZE)?latticeIndex+(latticeXSize*(MPD_WINDOW_APRON_SIZE-1)):latticeIndex+(latticeXSize*(MPD_WINDOW_APRON_SIZE-1))-(latticeXSize*latticeYSize);
        unsigned int apronLatticeSegmentIndex = latticeSegmentIndex+(MPD_Y_BLOCK_X_SIZE*(MPD_WINDOW_APRON_SIZE-1));
        makeDiffusionChoices(apronLatticeIndex, latticeSegment, apronLatticeSegmentIndex-MPD_Y_BLOCK_X_SIZE, apronLatticeSegmentIndex, apronLatticeSegmentIndex+MPD_Y_BLOCK_X_SIZE, choices, timestepHash);
    }

    // Wait for the whole segment to be processed.
    __syncthreads();
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
__device__ void copyZWindowFromLattice(const unsigned int* inLattice, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int latticeXYSize, const unsigned int latticeZSize, unsigned int* latticeSegment, const unsigned int latticeSegmentIndex, const unsigned int latticeSegmentZIndex, const unsigned long long timestepHash)
{
    // Load the block.
    latticeSegment[latticeSegmentIndex] = inLattice[latticeIndex];

    // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
    if (latticeSegmentZIndex < 2*MPD_WINDOW_APRON_SIZE) latticeSegment[latticeSegmentIndex-(MPD_Z_BLOCK_X_SIZE*MPD_WINDOW_APRON_SIZE)] = (latticeZIndex>=MPD_Z_BLOCK_Z_SIZE)?inLattice[latticeIndex-(latticeXYSize*MPD_WINDOW_APRON_SIZE)]:inLattice[latticeIndex-(latticeXYSize*MPD_WINDOW_APRON_SIZE)+(latticeXYSize*latticeZSize)];

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
    if (latticeSegmentZIndex >= MPD_Z_BLOCK_Z_SIZE) latticeSegment[latticeSegmentIndex+(MPD_Z_BLOCK_X_SIZE*MPD_WINDOW_APRON_SIZE)] = (latticeZIndex<latticeZSize-MPD_Z_BLOCK_Z_SIZE)?inLattice[latticeIndex+(latticeXYSize*MPD_WINDOW_APRON_SIZE)]:inLattice[latticeIndex+(latticeXYSize*MPD_WINDOW_APRON_SIZE)-(latticeXYSize*latticeZSize)];

    // Wait for the whole region to be copied.
    __syncthreads();
}

__device__ void makeZDiffusionChoices(const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int latticeXYSize, const unsigned int latticeZSize, const unsigned int* latticeSegment, const unsigned int latticeSegmentIndex, const unsigned int latticeSegmentZIndex, unsigned int* choices, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoices(latticeIndex, latticeSegment, latticeSegmentIndex-MPD_Z_BLOCK_X_SIZE, latticeSegmentIndex, latticeSegmentIndex+MPD_Z_BLOCK_X_SIZE, choices, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (latticeSegmentZIndex < (2*MPD_WINDOW_APRON_SIZE)-1)
    {
        unsigned int apronLatticeIndex = (latticeZIndex>=MPD_Z_BLOCK_Z_SIZE)?latticeIndex-(latticeXYSize*(MPD_WINDOW_APRON_SIZE-1)):latticeIndex-(latticeXYSize*(MPD_WINDOW_APRON_SIZE-1))+(latticeXYSize*latticeZSize);
        unsigned int apronLatticeSegmentIndex = latticeSegmentIndex-(MPD_Z_BLOCK_X_SIZE*(MPD_WINDOW_APRON_SIZE-1));
        makeDiffusionChoices(apronLatticeIndex, latticeSegment, apronLatticeSegmentIndex-MPD_Z_BLOCK_X_SIZE, apronLatticeSegmentIndex, apronLatticeSegmentIndex+MPD_Z_BLOCK_X_SIZE, choices, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (latticeSegmentZIndex >= MPD_Z_BLOCK_Z_SIZE+1)
    {
        unsigned int apronLatticeIndex = (latticeZIndex<latticeZSize-MPD_Z_BLOCK_Z_SIZE)?latticeIndex+(latticeXYSize*(MPD_WINDOW_APRON_SIZE-1)):latticeIndex+(latticeXYSize*(MPD_WINDOW_APRON_SIZE-1))-(latticeXYSize*latticeZSize);
        unsigned int apronLatticeSegmentIndex = latticeSegmentIndex+(MPD_Z_BLOCK_X_SIZE*(MPD_WINDOW_APRON_SIZE-1));
        makeDiffusionChoices(apronLatticeIndex, latticeSegment, apronLatticeSegmentIndex-MPD_Z_BLOCK_X_SIZE, apronLatticeSegmentIndex, apronLatticeSegmentIndex+MPD_Z_BLOCK_X_SIZE, choices, timestepHash);

    }

    // Wait for the whole segment to be processed.
    __syncthreads();
}
