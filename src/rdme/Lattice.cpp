/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
 * Copyright 2012 Roberts Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 * 
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
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
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
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


#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include "config.h"
#include "core/Types.h"
#include "core/Exceptions.h"
#include "rdme/Lattice.h"

namespace lm {
namespace rdme {

// Function to return the compiled-in value for MPD_LATTICE_MAX_OCCUPANCY.  Used 
// by pyLM in order to request the appropriate density.
unsigned int getCompiledLatticeMaxOccupancy()
{
	return MPD_LATTICE_MAX_OCCUPANCY;
}

Lattice::Lattice(lattice_coord_t size, si_dist_t spacing)
:size(size),numberSites(size.x*size.y*size.z),spacing(spacing)
{}

Lattice::Lattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing)
:size(xSize,ySize,zSize),numberSites(size.x*size.y*size.z),spacing(spacing)
{
}

Lattice::~Lattice() {}

lattice_coord_t Lattice::getSize() const
{
    return size;
}

lattice_size_t Lattice::getXSize() const
{
	return size.x;
}

lattice_size_t Lattice::getYSize() const
{
	return size.y;
}

lattice_size_t Lattice::getZSize() const
{
	return size.z;
}

lattice_size_t Lattice::getNumberSites() const
{
	return numberSites;
}

si_dist_t Lattice::getSpacing() const
{
	return spacing;
}

std::vector<lattice_coord_t> Lattice::getNearbySites(lattice_size_t xc, lattice_size_t yc, lattice_size_t zc, uint minDistance, uint maxDistance)
{
	std::vector<lattice_coord_t> sites;
	double minDistanceSquared=(double)(minDistance*minDistance);
	double maxDistanceSquared=(double)(maxDistance*maxDistance);
	
	for (lattice_size_t x = (xc>=maxDistance)?(xc-maxDistance):0; x <= ((xc<=size.x-1-maxDistance)?(xc+maxDistance):(size.x-1)); x++)
	{
		for (lattice_size_t y = (yc>=maxDistance)?(yc-maxDistance):0; y <= ((yc<=size.y-1-maxDistance)?(yc+maxDistance):(size.y-1)); y++)
		{
			for (lattice_size_t z = (zc>=maxDistance)?(zc-maxDistance):0; z <= ((zc<=size.z-1-maxDistance)?(zc+maxDistance):(size.z-1)); z++)
			{
				double rSquared = pow((double)x-(double)xc,2.0)+pow((double)y-(double)yc,2.0)+pow((double)z-(double)zc,2.0);
				if (rSquared <= ((double)maxDistanceSquared) && (rSquared > minDistanceSquared || minDistance == maxDistance))
				{
					sites.push_back(lattice_coord_t(x,y,z));
				}
			}
		}
	}
	
	return sites;
}

void Lattice::removeAllParticles()
{
    for (lattice_size_t z=0; z<size.z; z++)
        for (lattice_size_t y=0; y<size.y; y++)
            for (lattice_size_t x=0; x<size.x; x++)
                removeParticles(x,y,z);
}

void Lattice::rowMajorByteSerialize(void * bufferPtr, void * latticePtr, size_t bufferSize)
{
    Lattice * lattice = (Lattice *)latticePtr;
    byte * buffer = (byte *)bufferPtr;
    lattice_coord_t size = lattice->size;
    if (bufferSize != size.x*size.y*size.z*lattice->getMaxOccupancy()) throw lm::InvalidArgException(__PRETTY_FUNCTION__, "the buffer size was not equal to the lattice data size", bufferSize, size.x*size.y*size.z*lattice->getMaxOccupancy());
    for (lattice_size_t x=0, index=0; x<size.x; x++)
        for (lattice_size_t y=0; y<size.y; y++)
            for (lattice_size_t z=0; z<size.z; z++)
                for (site_size_t p=0; p<lattice->getMaxOccupancy(); p++, index++)
                    buffer[index] = lattice->getParticle(x,y,z,p);
}

void Lattice::rowMajorIntSerialize(void * bufferPtr, void * latticePtr, size_t bufferSize)
{
    Lattice * lattice = (Lattice *)latticePtr;
    uint32_t * buffer = (uint32_t *)bufferPtr;
    lattice_coord_t size = lattice->size;
    if (bufferSize != size.x*size.y*size.z*lattice->getMaxOccupancy()) throw lm::InvalidArgException(__PRETTY_FUNCTION__, "the buffer size was not equal to the lattice data size", bufferSize, size.x*size.y*size.z*lattice->getMaxOccupancy());
    for (lattice_size_t x=0, index=0; x<size.x; x++)
        for (lattice_size_t y=0; y<size.y; y++)
            for (lattice_size_t z=0; z<size.z; z++)
                for (site_size_t p=0; p<lattice->getMaxOccupancy(); p++, index++)
                    buffer[index] = lattice->getParticle(x,y,z,p);
}

void Lattice::rowMajorByteSerializeSites(void * bufferPtr, void * latticePtr, size_t bufferSize)
{
    Lattice * lattice = (Lattice *)latticePtr;
    byte * buffer = (byte *)bufferPtr;
    lattice_coord_t size = lattice->size;
    if (bufferSize != size.x*size.y*size.z) throw lm::InvalidArgException(__PRETTY_FUNCTION__, "the buffer size was not equal to the lattice site data size");
    for (lattice_size_t x=0, index=0; x<size.x; x++)
        for (lattice_size_t y=0; y<size.y; y++)
            for (lattice_size_t z=0; z<size.z; z++, index++)
				buffer[index] = (byte)lattice->getSiteType(x,y,z);
}

void Lattice::print() const
{
    for (lattice_size_t z=0; z<size.z; z++)
    {
        for (lattice_size_t x=0; x<size.x; x++)
        {
            for (lattice_size_t y=0; y<size.y; y++)
            {
                printf("(%2d)",getSiteType(x,y,z));
                for (site_size_t p=0; p<getMaxOccupancy(); p++)
                {
                    printf("%2d%c",getParticle(x,y,z,p),p<getMaxOccupancy()-1?',':' ');
                }
            }
            printf("\n");
        }
        printf("---------------\n");
    }
}

}
}
