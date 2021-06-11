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
LS_WORDS_PER_SITE
LS_APRON_SIZE
LS_XY_BLOCK_X_SIZE
LS_XY_BLOCK_Y_SIZE
LS_Z_BLOCK_Z_SIZE

Optional flags:
LS_BOUNDARY_PERIODIC        - defined if the lattice should laod with periodic boundary conditions
LS_BOUNDARY_VALUE           - the value to laod into lattice boundary sites (default 0)
LS_PACKED_SITES             - defined if objects are stored in lattice sites in such a way that
                              if one object position in a site is empty all higher words are empty
LS_PACKED_LAST_OBJECT_MASK  - or mask for the last object in a word
LS_LATTICE_SIZE_POWER_TWO   - defined if lattice x, y, and z sizes are enforced to be a power of two
*/

#if !defined LS_WORDS_PER_SITE
#error "Must define the number of 32-bit words per lattice site."
#endif

#if !defined LS_APRON_SIZE
#error "Must define the number of apron sites for the lattice window."
#endif

#if !defined LS_XY_BLOCK_X_SIZE || !defined LS_XY_BLOCK_Y_SIZE || !defined LS_Z_BLOCK_X_SIZE || !defined LS_Z_BLOCK_Z_SIZE
#error "Must define the x, y, and z dimensions of a block for xy and z window copying."
#endif

#if !defined LS_BOUNDARY_VALUE
#define LS_BOUNDARY_VALUE 0
#endif

#if defined LS_PACKED_SITES && !defined LS_PACKED_LAST_OBJECT_MASK
#error "Must set the last object mask when using packed sites."
#endif

#define LS_XY_WINDOW_X_SIZE             (LS_XY_BLOCK_X_SIZE+2*LS_APRON_SIZE)
#define LS_XY_WINDOW_Y_SIZE             (LS_XY_BLOCK_Y_SIZE+2*LS_APRON_SIZE)
#define LS_XY_WINDOW_SIZE               (LS_XY_WINDOW_X_SIZE*LS_XY_WINDOW_Y_SIZE)
#define LS_Z_WINDOW_SIZE                (LS_Z_BLOCK_X_SIZE*(LS_Z_BLOCK_Z_SIZE+2*LS_APRON_SIZE))

inline bool calculateXYLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int apronSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    if (latticeXSize%32 != 0 || latticeXSize%blockXSize != 0 || latticeYSize%blockYSize != 0)
        return false;

    // Calculate the grid size.
    *gridXSize = latticeXSize/blockXSize;
    gridSize->x = (*gridXSize)*(latticeYSize/blockYSize);
    gridSize->y = latticeZSize;
    gridSize->z = 1;
    threadBlockSize->x = blockXSize;
    threadBlockSize->y = blockYSize+2*apronSize;
    threadBlockSize->z = 1;
    return true;
}

/**
 * Calculates the x, y, and z indices of the current thread block. When the lattice size is not
 * guaranteed to be a power of two, the integer divide is only performed once in the first thread
 * for optimal performance.
 */
inline __device__ void calculateXYBlockIndices(unsigned int * bx, unsigned int * by, unsigned int * bz, unsigned int * gx, unsigned int * gy, unsigned int gridXSize)
{
    if (threadIdx.x == 0 && threadIdx.y == 0)
    {
        *by = blockIdx.x/gridXSize;
        *bx = blockIdx.x-gridXSize*(*by);
        *bz = blockIdx.y;
        *gx = gridXSize;
        *gy = gridDim.x/gridXSize;
    }
    __syncthreads();
}

inline __device__ void calculateXYThreadIndices(const unsigned int bx, const unsigned int by, const unsigned int bz, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYSize, int * blockXIndex, int * blockYIndex, int * latticeXIndex, int * latticeYIndex, unsigned int * latticeIndex, unsigned int * windowXIndex, unsigned int * windowYIndex, unsigned int * windowIndex)
{
     // Figure out the index of this thread in the block.
     *blockXIndex = threadIdx.x;
     *blockYIndex = threadIdx.y-LS_APRON_SIZE;

     // Figure out the index of this thread in the lattice.
     *latticeXIndex = (bx*LS_XY_BLOCK_X_SIZE) + *blockXIndex;
     *latticeYIndex = (by*LS_XY_BLOCK_Y_SIZE) + *blockYIndex;
     int latticeYWrap=0;
     if (*latticeYIndex < 0)
         latticeYWrap = latticeYSize;
     else if (*latticeYIndex >= latticeYSize)
         latticeYWrap = -latticeYSize;
     *latticeIndex = (bz*latticeXYSize) + (unsigned int)(((*latticeYIndex+latticeYWrap)*latticeXSize) + *latticeXIndex);

     // Figure out the index of this thread in the window.
     *windowXIndex =  threadIdx.x+LS_APRON_SIZE;
     *windowYIndex =  threadIdx.y;
     *windowIndex = (*windowYIndex*LS_XY_WINDOW_X_SIZE) + *windowXIndex;
}

// On architectures that support warp voting and if sites are packed, skip loading the next word if no threads in the warp had a completely filled previous word.
#if defined LS_PACKED_SITES && __CUDA_ARCH__ >= 120
#define _LS_COND_LOAD_SITE(window, window2, windowIndex, value)\
        if (__any(window[windowIndex]&LS_PACKED_LAST_OBJECT_MASK)) window2[windowIndex] = value;\
        else window2[windowIndex] = 0;

#else
#define _LS_COND_LOAD_SITE(window, window2, windowIndex, value) window2[windowIndex] = value;
#endif

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
inline __device__ void copyXYWindowFromLattice(const unsigned int* lattice, unsigned int* window, const unsigned int latticeIndex, const int latticeYIndex, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeXYZSize, const unsigned int windowIndex, const unsigned int windowXIndex, const unsigned int bx, const unsigned int gx)
{
    //Create any needed pointer aliases.
    #if LS_WORDS_PER_SITE >= 2
    const unsigned int* lattice2 = lattice+latticeXYZSize;
    unsigned int* window2 = window+LS_XY_WINDOW_SIZE;
    #endif

    #if !defined LS_BOUNDARY_PERIODIC
    int isYBoundary = (latticeYIndex < 0 || latticeYIndex >= latticeYSize);
    #endif

    #if LS_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it.
    if (windowXIndex >= LS_XY_BLOCK_X_SIZE)
    {
        #if defined LS_BOUNDARY_PERIODIC
        window[windowIndex-LS_XY_BLOCK_X_SIZE] = (bx>0)?lattice[latticeIndex-LS_XY_BLOCK_X_SIZE]:lattice[latticeIndex-LS_XY_BLOCK_X_SIZE+latticeXSize];
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex-LS_XY_BLOCK_X_SIZE, ((bx>0)?lattice2[latticeIndex-LS_XY_BLOCK_X_SIZE]:lattice2[latticeIndex-LS_XY_BLOCK_X_SIZE+latticeXSize]));
        #endif
        #else
        int isBoundary = (bx == 0 || isYBoundary);
        window[windowIndex-LS_XY_BLOCK_X_SIZE] = (!isBoundary)?lattice[latticeIndex-LS_XY_BLOCK_X_SIZE]:LS_BOUNDARY_VALUE;
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex-LS_XY_BLOCK_X_SIZE, ((!isBoundary)?lattice2[latticeIndex-LS_XY_BLOCK_X_SIZE]:LS_BOUNDARY_VALUE));
        #endif
        #endif
    }
    #endif

    // Load the sites.
    #if defined LS_BOUNDARY_PERIODIC
    window[windowIndex] = lattice[latticeIndex];
    #if LS_WORDS_PER_SITE >= 2
    _LS_COND_LOAD_SITE(window, window2, windowIndex, lattice2[latticeIndex]);
    #endif
    #else
    window[windowIndex] = (!isYBoundary)?lattice[latticeIndex]:LS_BOUNDARY_VALUE;
    #if LS_WORDS_PER_SITE >= 2
    _LS_COND_LOAD_SITE(window, window2, windowIndex, ((!isYBoundary)?lattice2[latticeIndex]:LS_BOUNDARY_VALUE));
    #endif
    #endif

    // Load the trailing apron.
    #if LS_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it.
    if (windowXIndex < 2*LS_APRON_SIZE)
    {
        #if defined LS_BOUNDARY_PERIODIC
        window[windowIndex+LS_XY_BLOCK_X_SIZE] = (bx<gx-1)?lattice[latticeIndex+LS_XY_BLOCK_X_SIZE]:lattice[latticeIndex+LS_XY_BLOCK_X_SIZE-latticeXSize];
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex+LS_XY_BLOCK_X_SIZE, ((bx<gx-1)?lattice2[latticeIndex+LS_XY_BLOCK_X_SIZE]:lattice2[latticeIndex+LS_XY_BLOCK_X_SIZE-latticeXSize]));
        #endif
        #else
        int isBoundary = (bx == gx-1 || isYBoundary);
        window[windowIndex+LS_XY_BLOCK_X_SIZE] = (!isBoundary)?lattice[latticeIndex+LS_XY_BLOCK_X_SIZE]:LS_BOUNDARY_VALUE;
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex+LS_XY_BLOCK_X_SIZE, ((!isBoundary)?lattice2[latticeIndex+LS_XY_BLOCK_X_SIZE]:LS_BOUNDARY_VALUE));
        #endif
        #endif
    }
    #endif
}

inline __device__ void copyXYWindowToLattice(unsigned int* lattice, const unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeXYZSize, const unsigned int windowIndex, const int blockYIndex)
{
    if (blockYIndex >= 0 && blockYIndex < LS_XY_BLOCK_Y_SIZE)
    {
        lattice[latticeIndex] = window[windowIndex];
        #if LS_WORDS_PER_SITE >= 2
        lattice[latticeIndex+latticeXYZSize] = window[windowIndex+LS_XY_WINDOW_SIZE];
        #endif
    }
}

/**
 * Gets the launch parameters for launching a z diffusion kernel.
 */
inline bool calculateZLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    if (latticeXSize%32 != 0 || latticeXSize%blockXSize != 0 || latticeZSize%blockZSize != 0)
        return false;

    *gridXSize = latticeXSize/blockXSize;
    gridSize->x = (*gridXSize)*(latticeYSize);
    gridSize->y = latticeZSize/blockZSize;
    gridSize->z = 1;
    threadBlockSize->x = blockXSize;
    threadBlockSize->y = 1;
    threadBlockSize->z = blockZSize;
    return true;
}

/**
 * Calculates the x, y, and z indices of the current thread block. When the lattice size is not
 * guaranteed to be a power of two, the integer divide is only performed once in the first thread
 * for optimal performance.
 */
inline __device__ void calculateZBlockIndices(unsigned int * bx, unsigned int * by, unsigned int * bz, unsigned int gridXSize)
{
    if (threadIdx.x == 0 && threadIdx.z == 0)
    {
        *by = blockIdx.x/gridXSize;
        *bx = blockIdx.x-gridXSize*(*by);
        *bz = blockIdx.y;
    }
    __syncthreads();
}

/**
 * Copies a window of a lattice from device memory to shared memory.
 */
inline __device__ void copyZWindowFromLattice(const unsigned int* lattice, unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned int windowIndex, const unsigned int windowZIndex)
{
    // Load the block.
    window[windowIndex] = lattice[latticeIndex];

    // Load the next part of the site, if it exists.
    #if LS_WORDS_PER_SITE >= 2
    const unsigned int* lattice2 = lattice+latticeXYZSize;
    unsigned int* window2 = window+LS_Z_WINDOW_SIZE;
    _LS_COND_LOAD_SITE(window, window2, windowIndex, lattice2[latticeIndex]);
    #endif

    #if LS_APRON_SIZE > 0
    // If this thread is one that needs to load part of the leading apron, load it (unless this is the lattice start).
    if (windowZIndex < 2*LS_APRON_SIZE)
    {
        #if defined LS_BOUNDARY_PERIODIC
        window[windowIndex-(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE)] = (latticeZIndex>=LS_Z_BLOCK_Z_SIZE)?lattice[latticeIndex-(latticeXYSize*LS_APRON_SIZE)]:lattice[latticeIndex-(latticeXYSize*LS_APRON_SIZE)+(latticeXYZSize)];
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex-(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE), ((latticeZIndex>=LS_Z_BLOCK_Z_SIZE)?lattice2[latticeIndex-(latticeXYSize*LS_APRON_SIZE)]:lattice2[latticeIndex-(latticeXYSize*LS_APRON_SIZE)+(latticeXYZSize)]));
        #endif
        #else
        window[windowIndex-(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE)] = (latticeZIndex>=LS_Z_BLOCK_Z_SIZE)?lattice[latticeIndex-(latticeXYSize*LS_APRON_SIZE)]:LS_BOUNDARY_VALUE;
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex-(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE), ((latticeZIndex>=LS_Z_BLOCK_Z_SIZE)?lattice2[latticeIndex-(latticeXYSize*LS_APRON_SIZE)]:LS_BOUNDARY_VALUE));
        #endif
        #endif
    }

    // If this thread is one that needs to load part of the trailing apron, load it (unless this is the lattice end).
    if (windowZIndex >= LS_Z_BLOCK_Z_SIZE)
    {
        #if defined LS_BOUNDARY_PERIODIC
        window[windowIndex+(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE)] = (latticeZIndex<latticeZSize-LS_Z_BLOCK_Z_SIZE)?lattice[latticeIndex+(latticeXYSize*LS_APRON_SIZE)]:lattice[latticeIndex+(latticeXYSize*LS_APRON_SIZE)-(latticeXYZSize)];
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex+(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE), ((latticeZIndex<latticeZSize-LS_Z_BLOCK_Z_SIZE)?lattice2[latticeIndex+(latticeXYSize*LS_APRON_SIZE)]:lattice2[latticeIndex+(latticeXYSize*LS_APRON_SIZE)-(latticeXYZSize)]));
        #endif
        #else
        window[windowIndex+(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE)] = (latticeZIndex<latticeZSize-LS_Z_BLOCK_Z_SIZE)?lattice[latticeIndex+(latticeXYSize*LS_APRON_SIZE)]:LS_BOUNDARY_VALUE;
        #if LS_WORDS_PER_SITE >= 2
        _LS_COND_LOAD_SITE(window, window2, windowIndex+(LS_Z_BLOCK_X_SIZE*LS_APRON_SIZE), ((latticeZIndex<latticeZSize-LS_Z_BLOCK_Z_SIZE)?lattice2[latticeIndex+(latticeXYSize*LS_APRON_SIZE)]:LS_BOUNDARY_VALUE));
        #endif
        #endif
    }
    #endif
}

inline __device__ void copyZWindowToLattice(unsigned int* lattice, const unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeXYZSize, const unsigned int windowIndex)
{
    lattice[latticeIndex] = window[windowIndex];
    #if LS_WORDS_PER_SITE >= 2
    lattice[latticeIndex+latticeXYZSize] = window[windowIndex+LS_Z_WINDOW_SIZE];
    #endif
}

