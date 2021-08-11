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
LS_WINDOW_APRON_SIZE
LS_XYZ_BLOCK_X_SIZE
LS_XYZ_BLOCK_Y_SIZE
LS_XYZ_BLOCK_Z_SIZE

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

#if !defined LS_XYZ_WINDOW_X_SIZE || !defined LS_XYZ_WINDOW_Y_SIZE || !defined LS_XYZ_WINDOW_Z_SIZE
#error "Must define the x, y, and z dimensions of the window for xyz window copying."
#endif

#define LS_XYZ_X_THREADS LS_XYZ_WINDOW_X_SIZE
#define LS_XYZ_Y_THREADS LS_XYZ_WINDOW_Y_SIZE
#if !defined LS_XYZ_Z_THREADS
#define LS_XYZ_Z_THREADS LS_XYZ_WINDOW_Z_SIZE
#endif

#if LS_XYZ_WINDOW_Z_SIZE%LS_XYZ_Z_THREADS != 0
#error "The z window size must be evenly divisible by the number of z threads."
#else
#define LS_XYZ_Z_LOOPS LS_XYZ_WINDOW_Z_SIZE/LS_XYZ_Z_THREADS
#endif

#if !defined LS_BOUNDARY_VALUE
#define LS_BOUNDARY_VALUE 0
#endif

#if defined LS_PACKED_SITES && !defined LS_PACKED_LAST_OBJECT_MASK
#error "Must set the last object mask when using packed sites."
#endif

#define LS_XYZ_BLOCK_X_SIZE               (LS_XYZ_WINDOW_X_SIZE-(2*LS_APRON_SIZE))
#define LS_XYZ_BLOCK_Y_SIZE               (LS_XYZ_WINDOW_Y_SIZE-(2*LS_APRON_SIZE))
#define LS_XYZ_BLOCK_Z_SIZE               (LS_XYZ_WINDOW_Z_SIZE-(2*LS_APRON_SIZE))
#define LS_XYZ_BLOCK_SIZE                 (LS_XYZ_BLOCK_X_SIZE*LS_XYZ_BLOCK_Y_SIZE*LS_XYZ_BLOCK_Z_SIZE)
#define LS_XYZ_WINDOW_SIZE                (LS_XYZ_WINDOW_X_SIZE*LS_XYZ_WINDOW_Y_SIZE*LS_XYZ_WINDOW_Z_SIZE)



inline bool calculateLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, const unsigned int xThreads, const unsigned int yThreads, const unsigned int zThreads)
{
    if (latticeXSize%blockXSize != 0 || latticeYSize%blockYSize != 0 || latticeZSize%blockZSize != 0)
        return false;

    // Calculate the grid size.
    unsigned int gx = latticeXSize/blockXSize;
    unsigned int gy = latticeYSize/blockYSize;
    unsigned int gz = latticeZSize/blockZSize;
    *gridXSize = gx;
    gridSize->x = gx*gy;
    gridSize->y = gz;

    // Calculate the thread block size.
    threadBlockSize->x = xThreads;
    threadBlockSize->y = yThreads;
    threadBlockSize->z = zThreads;
    return true;
}

/**
 * Calculates the x, y, and z indices of the current thread block. When the lattice size is not
 * guaranteed to be a power of two, the integer divide is only performed once in the first thread
 * for optimal performance.
 */
inline __device__ void calculateBlockIndices(unsigned int * bx, unsigned int * by, unsigned int * bz, unsigned int gridXSize)
{
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
        *by = blockIdx.x/gridXSize;
        *bx = blockIdx.x-gridXSize*(*by);
        *bz = blockIdx.y;
    }
    __syncthreads();
}

inline __device__ void calculateThreadIndices(const unsigned int bx, const unsigned int by, const unsigned int bz, const unsigned int latticeXSize, const unsigned int latticeXYSize, int * blockXIndex, int * blockYIndex, int * blockZIndex, int * latticeXIndex, int * latticeYIndex, int * latticeZIndex, int * latticeIndex, unsigned int * windowIndex)
{
    // Figure out the index of this thread in the block.
    *blockXIndex =  threadIdx.x-LS_APRON_SIZE;
    *blockYIndex =  threadIdx.y-LS_APRON_SIZE;
    *blockZIndex =  threadIdx.z-LS_APRON_SIZE;

    // Figure out the index of this thread in the lattice.
    *latticeXIndex = (bx*LS_XYZ_BLOCK_X_SIZE) + *blockXIndex;
    *latticeYIndex = (by*LS_XYZ_BLOCK_Y_SIZE) + *blockYIndex;
    *latticeZIndex = (bz*LS_XYZ_BLOCK_Z_SIZE) + *blockZIndex;
    *latticeIndex = (*latticeZIndex*latticeXYSize) + (*latticeYIndex*latticeXSize) + *latticeXIndex;

    // Figure out the index of this thread in the window.
    *windowIndex = (threadIdx.z*LS_XYZ_WINDOW_X_SIZE*LS_XYZ_WINDOW_Y_SIZE) + (threadIdx.y*LS_XYZ_WINDOW_X_SIZE) + threadIdx.x;
}



// On architectures that support warp voting and if sites are packed, skip loading the next word if no threads in the warp had a completely filled previous word.
#if defined LS_PACKED_SITES && __CUDA_ARCH__ >= 120

#define _LS_COND_LOAD_SITE(window, window2, windowIndex, value)\
        if (__any(window[windowIndex]&LS_PACKED_LAST_OBJECT_MASK)) window2[windowIndex] = value;\
        else window2[windowIndex] = 0;

#else
#define _LS_COND_LOAD_SITE(window, window2, windowIndex, value) window2[windowIndex] = value;
#endif

#define LS_Z_LOOP_Z_INDEX(baseIndex,loopCounter) (baseIndex+(loopCounter*LS_XYZ_Z_THREADS))
#define LS_Z_LOOP_LATTICE_INDEX(baseIndex,loopCounter,latticeXYSize) (baseIndex+(loopCounter*LS_XYZ_Z_THREADS*latticeXYSize))
#define LS_Z_LOOP_WINDOW_INDEX(baseIndex,loopCounter) (baseIndex+(loopCounter*LS_XYZ_Z_THREADS*LS_XYZ_WINDOW_X_SIZE*LS_XYZ_WINDOW_Y_SIZE))


/**
 * Copies a segment of a lattice from device memory to shared memory.
 */
inline __device__ void copyXYZWindowFromLattice(const unsigned int* lattice, unsigned int* window, const int latticeIndex, const int latticeXIndex, const int latticeYIndex, const int latticeZIndex, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned int windowIndex)
{
    //Create any needed pointer aliases.
    #if LS_WORDS_PER_SITE >= 2
    const unsigned int* lattice2 = lattice+latticeXYZSize;
    unsigned int* window2 = window+LS_XYZ_WINDOW_SIZE;
    #endif

    // Load the window while looping through the z planes.
    for (int i=0; i<LS_XYZ_Z_LOOPS; i++)
    {
        int loopLatticeZIndex = LS_Z_LOOP_Z_INDEX(latticeZIndex,i);
        int loopLatticeIndex = LS_Z_LOOP_LATTICE_INDEX(latticeIndex,i,latticeXYSize);
        unsigned int loopWindowIndex = LS_Z_LOOP_WINDOW_INDEX(windowIndex,i);

        #if LS_APRON_SIZE > 0
        // See if we are inside the lattice.
        if (latticeXIndex >= 0 && latticeXIndex < latticeXSize && latticeYIndex >= 0 && latticeYIndex < latticeYSize && loopLatticeZIndex >= 0 && loopLatticeZIndex < latticeZSize)
        {
        #endif
            // Load the first word.
            window[loopWindowIndex] = lattice[loopLatticeIndex];

            // Load the second word, if it exists.
            #if LS_WORDS_PER_SITE >= 2
            _LS_COND_LOAD_SITE(window, window2, loopWindowIndex, lattice2[loopLatticeIndex]);
            #endif

        #if LS_APRON_SIZE > 0
        }
        // Otherwise, fill in the lattice boundary.
        else
        {
            #if defined LS_BOUNDARY_PERIODIC
            #error "Unimplemented."
            #else
            window[loopWindowIndex] = LS_BOUNDARY_VALUE;
            #endif
        }
        #endif
    }
}

inline __device__ void copyXYZWindowToLattice(unsigned int* lattice, const unsigned int* window, const int latticeIndex, const unsigned int latticeXYSize, const unsigned int latticeXYZSize, const unsigned int windowIndex, const int blockXIndex, const int blockYIndex, const int blockZIndex)
{
    //Create any needed pointer aliases.
    #if LS_WORDS_PER_SITE >= 2
    unsigned int* lattice2 = lattice+latticeXYZSize;
    const unsigned int* window2 = window+LS_XYZ_WINDOW_SIZE;
    #endif

    // Save the block while looping through the z planes.
    for (int i=0; i<LS_XYZ_Z_LOOPS; i++)
    {
        // Figure out the new indices for this z loop.
        int loopBlockZIndex = LS_Z_LOOP_Z_INDEX(blockZIndex,i);
        int loopLatticeIndex = LS_Z_LOOP_LATTICE_INDEX(latticeIndex,i,latticeXYSize);
        unsigned int loopWindowIndex = LS_Z_LOOP_WINDOW_INDEX(windowIndex,i);


        // See if we are inside the block.
        #if LS_APRON_SIZE > 0
        if (blockXIndex >= 0 && blockXIndex < LS_XYZ_BLOCK_X_SIZE && blockYIndex >= 0 && blockYIndex < LS_XYZ_BLOCK_Y_SIZE && loopBlockZIndex >= 0 && loopBlockZIndex < LS_XYZ_BLOCK_Z_SIZE)
        {
        #endif
            // Save the first word.
            lattice[loopLatticeIndex] = window[loopWindowIndex];

            // Save the second word, if it exists.
            #if LS_WORDS_PER_SITE >= 2
            lattice2[loopLatticeIndex] = window2[loopWindowIndex];
            #endif

        #if LS_APRON_SIZE > 0
        }
        #endif

    }
}

