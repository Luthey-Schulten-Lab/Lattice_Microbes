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

/*
Required flags:
MPD_WORDS_PER_SITE
MPD_APRON_SIZE
TUNE_MPD_X_BLOCK_MAX_X_SIZE
TUNE_MPD_Y_BLOCK_X_SIZE
TUNE_MPD_Y_BLOCK_Y_SIZE
TUNE_MPD_Z_BLOCK_X_SIZE
TUNE_MPD_Z_BLOCK_Z_SIZE

Optional flags:
MPD_BOUNDARY_PERIODIC        - defined if the lattice should laod with periodic boundary conditions
MPD_BOUNDARY_PARTICLE_VALUE  - the value to laod into lattice boundary sites (default 0)
MPD_BOUNDARY_SITE_VALUE      - the value to load into lattice boundary sites (default 0xFF)
TUNE_MPD_PACKED_SITES        - defined if objects are stored in lattice sites in such a way that
                              if one object position in a site is empty all higher words are empty
MPD_PACKED_LAST_OBJECT_MASK  - or mask for the last object in a word
MPD_LATTICE_SIZE_POWER_TWO   - defined if lattice x, y, and z sizes are enforced to be a power of two
*/

#if !defined MPD_WORDS_PER_SITE
#error "Must define the number of 32-bit words per lattice site."
#endif

#if !defined MPD_APRON_SIZE
#error "Must define the number of apron sites for the lattice window."
#endif

#if !defined TUNE_MPD_X_BLOCK_MAX_X_SIZE
#error "Must define the maximum x dimensions of a block for x window copying."
#endif

#if !defined TUNE_MPD_Y_BLOCK_X_SIZE || !defined TUNE_MPD_Y_BLOCK_Y_SIZE
#error "Must define the x and y dimensions of a block for y window copying."
#endif

#if !defined TUNE_MPD_Z_BLOCK_X_SIZE || !defined TUNE_MPD_Z_BLOCK_Z_SIZE
#error "Must define the x and z dimensions of a block for z window copying."
#endif

#if !defined MPD_BOUNDARY_PARTICLE_VALUE
#define MPD_BOUNDARY_PARTICLE_VALUE 0
#endif

#if !defined MPD_BOUNDARY_SITE_VALUE
#define MPD_BOUNDARY_SITE_VALUE 0xFF
#endif

#define MPD_X_WINDOW_SIZE               (MPD_APRON_SIZE+TUNE_MPD_X_BLOCK_MAX_X_SIZE+MPD_APRON_SIZE)
#define MPD_Y_WINDOW_SIZE               (TUNE_MPD_Y_BLOCK_X_SIZE*(MPD_APRON_SIZE+TUNE_MPD_Y_BLOCK_Y_SIZE+MPD_APRON_SIZE))
#define MPD_Z_WINDOW_SIZE               (TUNE_MPD_Z_BLOCK_X_SIZE*(MPD_APRON_SIZE+TUNE_MPD_Z_BLOCK_Z_SIZE+MPD_APRON_SIZE))


/**
 * Calculates the x, y, and z indices of the current thread block. When the lattice size is not
 * guaranteed to be a power of two, the integer divide is only performed once in the first thread
 * for optimal performance.
 */
inline __device__ void calculateBlockPosition(unsigned int * bx, unsigned int * by, unsigned int * bz, unsigned int gridXSize)
{
    #if defined MPD_LATTICE_SIZE_POWER_TWO
    #error "Unimplemented."
    //const unsigned int bx = blockIdx.x&xBlockMaskMod;
    //const unsigned int by = blockIdx.x>>yBlockShiftMult;
    //const unsigned int bz = blockIdx.y;
    #else
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
        *by = blockIdx.x/gridXSize;
        *bx = blockIdx.x-gridXSize*(*by);
        *bz = blockIdx.y;
    }
    __syncthreads();
    #endif
}

#define _MPD_COND_LOAD_SITE(window, window2, windowIndex, value) window2[windowIndex] = value;

/**
 * Copies a segment of a lattice from device memory to shared memory.
 */
inline __device__ void copyXWindowFromLattice(const unsigned int bx, const unsigned int * lattice, unsigned int * window, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex)
{
    if (latticeXIndex < latticeXSizeC)
    {
        // Load the block.
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
            window[windowIndex+windowOffset] = lattice[latticeIndex+latticeOffset];

        #if MPD_APRON_SIZE > 0
        // Calculate the effective thread block width, accounting for the fact that the block might not be aligned with the end of the lattice.
        int threadBlockWidth = ((bx+1)*blockDim.x <= latticeXSizeC)?(blockDim.x):(latticeXSizeC-(bx*blockDim.x));

        // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
        if (windowIndex >= threadBlockWidth)
        {
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
                #if defined MPD_BOUNDARY_PERIODIC
                window[windowIndex+windowOffset-threadBlockWidth] = (latticeXIndex>=threadBlockWidth)?lattice[latticeIndex+latticeOffset-threadBlockWidth]:lattice[latticeIndex+latticeOffset-threadBlockWidth+latticeXSizeC];
                #else
                window[windowIndex+windowOffset-threadBlockWidth] = (latticeXIndex>=threadBlockWidth)?lattice[latticeIndex+latticeOffset-threadBlockWidth]:MPD_BOUNDARY_PARTICLE_VALUE;
                #endif
        }

        // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
        if (windowIndex < 2*MPD_APRON_SIZE)
        {
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
                #if defined MPD_BOUNDARY_PERIODIC
                window[windowIndex+windowOffset+threadBlockWidth] = (latticeXIndex<latticeXSizeC-threadBlockWidth)?lattice[latticeIndex+latticeOffset+threadBlockWidth]:lattice[latticeIndex+latticeOffset+threadBlockWidth-latticeXSizeC];
                #else
                window[windowIndex+windowOffset+threadBlockWidth] = (latticeXIndex<latticeXSizeC-threadBlockWidth)?lattice[latticeIndex+latticeOffset+threadBlockWidth]:MPD_BOUNDARY_PARTICLE_VALUE;
                #endif
        }
        #endif
    }
}

inline __device__ void copyXWindowFromBuffer(const unsigned int bx, const unsigned int * buffer, unsigned int * window, const unsigned int bufIdx, const unsigned int latticeXIndex, const unsigned int windowIndex, const unsigned int packSize)
{
    if (latticeXIndex < latticeXSizeC)
    {
        // Load the block.
        for (uint w=0, packOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, packOffset+=packSize, windowOffset+=MPD_X_WINDOW_SIZE)
			window[windowIndex+windowOffset] = buffer[bufIdx+packOffset];

        #if MPD_APRON_SIZE > 0
        // Calculate the effective thread block width, accounting for the fact that the block might not be aligned with the end of the lattice.
        int threadBlockWidth = ((bx+1)*blockDim.x <= latticeXSizeC)?(blockDim.x):(latticeXSizeC-(bx*blockDim.x));

        // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
        if (windowIndex >= threadBlockWidth)
        {
			for (uint w=0, packOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, packOffset+=packSize, windowOffset+=MPD_X_WINDOW_SIZE)
                #if defined MPD_BOUNDARY_PERIODIC
                window[windowIndex+windowOffset-threadBlockWidth] = (latticeXIndex>=threadBlockWidth)?buffer[bufIdx+packOffset-threadBlockWidth]:buffer[bufIdx+packOffset-threadBlockWidth+latticeXSizeC];
                #else
				#pragma error "oops not supported"
                #endif
        }

        // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
        if (windowIndex < 2*MPD_APRON_SIZE)
        {
			for (uint w=0, packOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, packOffset+=packSize, windowOffset+=MPD_X_WINDOW_SIZE)
                #if defined MPD_BOUNDARY_PERIODIC
                window[windowIndex+windowOffset+threadBlockWidth] = (latticeXIndex<latticeXSizeC-threadBlockWidth)?buffer[bufIdx+packOffset+threadBlockWidth]:buffer[bufIdx+packOffset+threadBlockWidth-latticeXSizeC];
                #else
				#pragma error "oops not supported"
                #endif
        }
        #endif
    }
}

inline __device__ void copyXWindowFromSites(const unsigned int bx, const uint8_t * sites, uint8_t * sitesWindow, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex)
{
    if (latticeXIndex < latticeXSizeC)
    {
        // Load the block.
        sitesWindow[windowIndex] = sites[latticeIndex];

        #if MPD_APRON_SIZE > 0
        // Calculate the effective thread block width, accounting for the fact that the block might not be aligned with the end of the lattice.
        int threadBlockWidth = ((bx+1)*blockDim.x <= latticeXSizeC)?(blockDim.x):(latticeXSizeC-(bx*blockDim.x));

        // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
        if (windowIndex >= threadBlockWidth)
        {
            #if defined MPD_BOUNDARY_PERIODIC
            sitesWindow[windowIndex-threadBlockWidth] = (latticeXIndex>=threadBlockWidth)?sites[latticeIndex-threadBlockWidth]:sites[latticeIndex-threadBlockWidth+latticeXSizeC];
            #else
            sitesWindow[windowIndex-threadBlockWidth] = (latticeXIndex>=threadBlockWidth)?sites[latticeIndex-threadBlockWidth]:MPD_BOUNDARY_SITE_VALUE;
            #endif
        }

        // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
        if (windowIndex < 2*MPD_APRON_SIZE)
        {
            #if defined MPD_BOUNDARY_PERIODIC
            sitesWindow[windowIndex+threadBlockWidth] = (latticeXIndex<latticeXSizeC-threadBlockWidth)?sites[latticeIndex+threadBlockWidth]:sites[latticeIndex+threadBlockWidth-latticeXSizeC];
            #else
            sitesWindow[windowIndex+threadBlockWidth] = (latticeXIndex<latticeXSizeC-threadBlockWidth)?sites[latticeIndex+threadBlockWidth]:MPD_BOUNDARY_SITE_VALUE;
            #endif
        }
        #endif
    }
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
inline __device__ void copyYWindowFromLattice(const unsigned int* lattice, unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeYIndex, const unsigned int windowIndex, const unsigned int windowYIndex)
{
    // Load the block.
    for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
        window[windowIndex+windowOffset] = lattice[latticeIndex+latticeOffset];

    #if MPD_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it (unless this is in the first block).
    if (windowYIndex < 2*MPD_APRON_SIZE)
    {
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
            #if defined MPD_BOUNDARY_PERIODIC
            window[windowIndex+windowOffset-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?lattice[latticeIndex+latticeOffset-(latticeXSizeC*MPD_APRON_SIZE)]:lattice[latticeIndex+latticeOffset-(latticeXSizeC*MPD_APRON_SIZE)+(latticeXSizeC*latticeYSizeC)];
            #else
            window[windowIndex+windowOffset-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?lattice[latticeIndex+latticeOffset-(latticeXSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
    }

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is in the last block).
    if (windowYIndex >= TUNE_MPD_Y_BLOCK_Y_SIZE)
    {
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
            #if defined MPD_BOUNDARY_PERIODIC
            window[windowIndex+windowOffset+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?lattice[latticeIndex+latticeOffset+(latticeXSizeC*MPD_APRON_SIZE)]:lattice[latticeIndex+latticeOffset+(latticeXSizeC*MPD_APRON_SIZE)-(latticeXSizeC*latticeYSizeC)];
            #else
            window[windowIndex+windowOffset+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?lattice[latticeIndex+latticeOffset+(latticeXSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
    }
    #endif
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
inline __device__ void copyYWindowFromSites(const uint8_t * sites, uint8_t * sitesWindow, const unsigned int latticeIndex, const unsigned int latticeYIndex, const unsigned int windowIndex, const unsigned int windowYIndex)
{
    // Load the block.
    sitesWindow[windowIndex] = sites[latticeIndex];

    #if MPD_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it (unless this is in the first block).
    if (windowYIndex < 2*MPD_APRON_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        sitesWindow[windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?sites[latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)]:sites[latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)+(latticeXSizeC*latticeYSizeC)];
        #else
        sitesWindow[windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?sites[latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_SITE_VALUE;
        #endif
    }

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is in the last block).
    if (windowYIndex >= TUNE_MPD_Y_BLOCK_Y_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        sitesWindow[windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?sites[latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)]:sites[latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)-(latticeXSizeC*latticeYSizeC)];
        #else
        sitesWindow[windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?sites[latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_SITE_VALUE;
        #endif
    }
    #endif
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
inline __device__ void copyZWindowFromLattice(const unsigned int* lattice, unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex)
{
    // Load the block.
    for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
        window[windowIndex+windowOffset] = lattice[latticeIndex+latticeOffset];

    #if MPD_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
    if (windowZIndex < 2*MPD_APRON_SIZE)
    {
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
            #if defined MPD_BOUNDARY_PERIODIC
            window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?lattice[latticeIndex+latticeOffset-(latticeXYSizeC*MPD_APRON_SIZE)]:lattice[latticeIndex+latticeOffset-(latticeXYSizeC*MPD_APRON_SIZE)+(latticeXYZSizeC)];
            #else
            window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?lattice[latticeIndex+latticeOffset-(latticeXYSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
    }

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
            #if defined MPD_BOUNDARY_PERIODIC
            window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?lattice[latticeIndex+latticeOffset+(latticeXYSizeC*MPD_APRON_SIZE)]:lattice[latticeIndex+latticeOffset+(latticeXYSizeC*MPD_APRON_SIZE)-(latticeXYZSizeC)];
            #else
            window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?lattice[latticeIndex+latticeOffset+(latticeXYSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
    }
    #endif
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 * The Z Window is different in the multi-GPU case.  For periodic boundaries, the lattice is extended to include
 * the periodic portion from the remote GPU.  There is no need to wrap around.  On the non-periodic side, the indexing
 * has started not at the edge since the z kernel is not evaluated on the apron sites.
 * 
 */
inline __device__ void copyZWindowFromLattice_MGPU(const unsigned int* lattice, unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex)
{
    // Load the block.
    for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
        window[windowIndex+windowOffset] = lattice[latticeIndex+latticeOffset];

    #if MPD_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
    if (windowZIndex < 2*MPD_APRON_SIZE)
    {
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
		{
            #if defined MPD_BOUNDARY_PERIODIC
			window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = lattice[latticeIndex+latticeOffset-(latticeXYSizeC*MPD_APRON_SIZE)];
            #else
            window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE) ? 
				lattice[latticeIndex+latticeOffset-(latticeXYSizeC*MPD_APRON_SIZE)] :
				MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
		}
    }

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
            #if defined MPD_BOUNDARY_PERIODIC
            window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = lattice[latticeIndex+latticeOffset+(latticeXYSizeC*MPD_APRON_SIZE)];
            #else
            window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex<global_latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE) ? 
				lattice[latticeIndex+latticeOffset+(latticeXYSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
    }
    #endif
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
inline __device__ void copyZWindowFromSites(const uint8_t * sites, uint8_t * sitesWindow, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex)
{
    // Load the block.
    sitesWindow[windowIndex] = sites[latticeIndex];

    #if MPD_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
    if (windowZIndex < 2*MPD_APRON_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        sitesWindow[windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?sites[latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)]:sites[latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)+(latticeXYZSizeC)];
        #else
        sitesWindow[windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?sites[latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_SITE_VALUE;
        #endif
    }

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        sitesWindow[windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?sites[latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)]:sites[latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)-(latticeXYZSizeC)];
        #else
        sitesWindow[windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?sites[latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_SITE_VALUE;
        #endif
    }
    #endif
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
inline __device__ void copyZWindowFromSites_MGPU(const uint8_t * sites, uint8_t * sitesWindow, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex)
{
    // Load the block.
    sitesWindow[windowIndex] = sites[latticeIndex];

    #if MPD_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
    if (windowZIndex < 2*MPD_APRON_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        sitesWindow[windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = sites[latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)];
        #else
        sitesWindow[windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?sites[latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_SITE_VALUE;
        #endif
    }

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        sitesWindow[windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = sites[latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)];
        #else
        sitesWindow[windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (latticeZIndex<global_latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?sites[latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)]:MPD_BOUNDARY_SITE_VALUE;
        #endif
    }
    #endif
}

inline __device__ void copyXWindowToLattice(unsigned int* lattice, const unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex)
{
    if (latticeXIndex < latticeXSizeC)
    {
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
            lattice[latticeIndex+latticeOffset] = window[windowIndex+windowOffset];
    }
}


inline __device__ void copyYWindowToLattice(unsigned int* lattice, const unsigned int* window, const unsigned int latticeIndex, const unsigned int windowIndex)
{
    for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
        lattice[latticeIndex+latticeXYZSizeC] = window[windowIndex+windowOffset];
}


inline __device__ void copyZWindowToLattice(unsigned int* lattice, const unsigned int* window, const unsigned int latticeIndex, const unsigned int windowIndex)
{
    for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
        lattice[latticeIndex+latticeXYZSizeC] = window[windowIndex+windowOffset];
}




