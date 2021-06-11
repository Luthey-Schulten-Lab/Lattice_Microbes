/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
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

#include <map>
#include <vector>
#include <climits>
#include <cstring>
#include "config.h"
#include "core/Types.h"
#include "core/Exceptions.h"
#include "rdme/Lattice.h"
#include "rdme/IntLattice.h"

namespace lm {
namespace rdme {

site_t IntLattice::getMaxSiteType() const {return 255;}
particle_t IntLattice::getMaxParticle() const {return UINT_MAX - 1;}
site_size_t IntLattice::getMaxOccupancy() const {return PARTICLES_PER_WORD*wordsPerSite;}

IntLattice::IntLattice(lattice_coord_t size, si_dist_t spacing, uint particlesPerSite)
:Lattice(size,spacing),wordsPerSite(particlesPerSite/PARTICLES_PER_WORD),particles(NULL),siteTypes(NULL)
{
    // Make sure that the memory layout is compatible with our CUDA assumptions.
    if (sizeof(uint8_t)*CHAR_BIT != 8) throw Exception("a byte lattice can only be created on architectures with an 8-bit char type.");
    if (sizeof(uint32_t)*CHAR_BIT != 32) throw Exception("a byte lattice can only be created on architectures with a 32-bit int type.");
    if (wordsPerSite*PARTICLES_PER_WORD != particlesPerSite) throw InvalidArgException("particlesPerSite", "must be evenly divisible by the number of particles per word");
    allocateMemory();
}

IntLattice::IntLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing, uint particlesPerSite)
:Lattice(xSize,ySize,zSize,spacing),wordsPerSite(particlesPerSite/PARTICLES_PER_WORD),particles(NULL),siteTypes(NULL)
{
    if (sizeof(uint8_t)*CHAR_BIT != 8) throw Exception("a byte lattice can only be created on architectures with an 8-bit char type.");
    if (sizeof(uint32_t)*CHAR_BIT != 32) throw Exception("a byte lattice can only be created on architectures with a 32-bit int type.");
    if (wordsPerSite*PARTICLES_PER_WORD != particlesPerSite) throw InvalidArgException("particlesPerSite", "must be evenly divisible by the number of particles per word");
    allocateMemory();
}

IntLattice::~IntLattice()
{
    deallocateMemory();
}

void IntLattice::getNeighboringSites(lattice_size_t index, lattice_size_t * neighboringIndices)
{
    //         y+1
    //         |  z+1
    //         | /
    //   x-1---i---x+1 
    //        /|
    //     z-1 |
    //         y-1
	lattice_size_t z = index/(size.x*size.y);
	lattice_size_t xy = index-(z*size.x*size.y);
	lattice_size_t y = xy/size.x;
	lattice_size_t x = xy-(y*size.x);

	lattice_size_t xySize = size.x*size.y;
	neighboringIndices[0] = (x>0)?(index-1):(index-1+size.x);
	neighboringIndices[1] = (x<(size.x-1))?(index+1):(index+1-size.x);
	neighboringIndices[2] = (y>0)?(index-size.x):(index-size.x+xySize);
	neighboringIndices[3] = (y<(size.y-1))?(index+size.x):(index+size.x-xySize);
	neighboringIndices[4] = (z>0)?(index-xySize):(index-xySize+numberSites);
	neighboringIndices[5] = (z<(size.z-1))?(index+xySize):(index+xySize-numberSites);
}

size_t IntLattice::getLatticeMemorySize() const
{
	return (size_t(numberSites)*size_t(wordsPerSite)*sizeof(uint32_t));
}

void IntLattice::allocateMemory()
{
    // Allocate the data for the particles.
    particles = new uint32_t[numberSites*wordsPerSite];
    memset(particles, 0, numberSites*wordsPerSite*sizeof(uint32_t));

    // Allocate the data for the site types.
    siteTypes = new uint8_t[numberSites];
    memset(siteTypes, 0, numberSites*sizeof(uint8_t));
}

void IntLattice::deallocateMemory()
{
    // Free any allocated memory.
    if (particles != NULL)
    {
        delete[] particles;
        particles = NULL;
    }
    if (siteTypes != NULL)
    {
        delete[] siteTypes;
        siteTypes = NULL;
    }
}

site_t IntLattice::getSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z) const
{
	// Make sure the arguments are valid.
	if (x >= size.x || y >= size.y || z >= size.z)
		throw InvalidSiteException(x,y,z);
	lattice_size_t latticeIndex = x+(y*size.x)+(z*size.x*size.y);
	
	return (site_t)siteTypes[latticeIndex];
}

site_t IntLattice::getSiteType(lattice_size_t index) const
{
	// Make sure the arguments are valid.
	if (index >= numberSites)
		throw InvalidSiteException(index);
	return (site_t)siteTypes[index];
}

void IntLattice::setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t siteType)
{
    // Make sure the arguments are valid.
    if (x >= size.x || y >= size.y || z >= size.z)
        throw InvalidSiteException(x,y,z);
    lattice_size_t latticeIndex = x+(y*size.x)+(z*size.x*size.y);

    siteTypes[latticeIndex] = (uint8_t)siteType;
}

void IntLattice::setSiteType(lattice_size_t index, site_t siteType)
{
    // Make sure the arguments are valid.
	if (index >= numberSites)
		throw InvalidSiteException(index);
    siteTypes[index] = (uint8_t)siteType;
}

site_size_t IntLattice::getOccupancy(lattice_size_t x, lattice_size_t y, lattice_size_t z) const
{
    // Make sure the arguments are valid.
    if (x >= size.x || y >= size.y || z >= size.z)
        throw InvalidSiteException(x,y,z);
    lattice_size_t latticeIndex = x+(y*size.x)+(z*size.x*size.y);
	
	// Count the number of particles at the site.
    site_size_t occupancy = 0;
	for (uint wi=0; wi<wordsPerSite; wi++)
	{
		if(particles[latticeIndex + wi * numberSites] == 0)
			break;

		occupancy++;
	}
		
	return occupancy;
}

site_size_t IntLattice::getOccupancy(lattice_size_t index) const
{
    // Make sure the arguments are valid.
	if (index >= numberSites)
		throw InvalidSiteException(index);

	// Count the number of particles at the site.
    site_size_t occupancy = 0;
	for (uint wi=0; wi<wordsPerSite; wi++)
	{
		if(particles[index + wi * numberSites] == 0)
			break;

		occupancy++;
	}
	return occupancy;
}

/**
 * Gets the current state of a give site in the lattice.
 * 
 * @param x	The zero based x index of the site to retrieve.
 * @param y	The zero based y index of the site to retrieve.
 * @param z	The zero based z index of the site to retrieve.
 * @return	The value in the lattice at the specified site.
 */
particle_t IntLattice::getParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_size_t particleIndex) const
{
	// Make sure the arguments are valid.
    if (x >= size.x || y >= size.y || z >= size.z) {
        throw InvalidSiteException(x,y,z);
    }
	if (particleIndex >= getMaxOccupancy()) {
		throw InvalidParticleException(particleIndex);
    }
    lattice_size_t latticeIndex = x+(y*size.x)+(z*size.x*size.y);
	
	latticeIndex += numberSites * particleIndex;
    return particles[latticeIndex];
}

particle_t IntLattice::getParticle(lattice_size_t index, site_size_t particleIndex) const
{
	// Make sure the arguments are valid.
	if (index >= numberSites)
		throw InvalidSiteException(index);
	if (particleIndex >= getMaxOccupancy())
		throw InvalidParticleException(particleIndex);

	index += numberSites * particleIndex;
    return particles[index];
}

/**
 * Sets the current state of a give site in the lattice.
 * 
 * @param x	The zero based x index of the site to set.
 * @param y	The zero based y index of the site to set.
 * @param z	The zero based z index of the site to set.
 * @param 
 */
void IntLattice::addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle)
{
	// Make sure the arguments are valid.
    if (x >= size.x || y >= size.y || z >= size.z)
        throw InvalidSiteException(x,y,z);
    lattice_size_t latticeIndex = x+(y*size.x)+(z*size.x*size.y);

    // Go through the particles at the site.
    for (uint wi=0; wi<wordsPerSite; wi++)
    {
		if (particles[latticeIndex + wi*numberSites] == 0)
		{
			particles[latticeIndex + wi*numberSites] = particle;
			return;
		}
    }

    // The site must have been full.
    throw InvalidParticleException(getMaxOccupancy());
}

void IntLattice::addParticle(lattice_size_t index, particle_t particle)
{
	// Make sure the arguments are valid.
	if (index >= numberSites)
		throw InvalidSiteException(index);

    // Go through the particles at the site.
    for (uint wi=0; wi<wordsPerSite; wi++)
    {
		if (particles[index + wi*numberSites] == 0)
		{
			particles[index + wi*numberSites] = particle;
			return;
		}
    }

    // The site must have been full.
    throw InvalidParticleException(getMaxOccupancy());
}

void IntLattice::removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z)
{
	// Make sure the arguments are valid.
    if (x >= size.x || y >= size.y || z >= size.z)
        throw InvalidSiteException(x,y,z);
    lattice_size_t latticeIndex = x+(y*size.x)+(z*size.x*size.y);
	
    // Reset all of the words for this site.
    for (uint wi=0; wi<wordsPerSite; wi++)
    {
        particles[latticeIndex] = 0;
        latticeIndex += numberSites;
    }
}

void IntLattice::removeParticles(lattice_size_t index)
{
	// Make sure the arguments are valid.
	if (index >= numberSites)
		throw InvalidSiteException(index);

    // Reset all of the words for this site.
    for (uint wi=0; wi<wordsPerSite; wi++)
    {
        particles[index] = 0;
        index += numberSites;
    }
}

void IntLattice::removeAllParticles()
{
    memset(particles, 0, numberSites*wordsPerSite*sizeof(uint32_t));
}

void IntLattice::nativeSerialize(void * destBuffer, void * latticePtr, size_t bufferSize)
{
    IntLattice * lattice = (IntLattice *)latticePtr;
    if (bufferSize != lattice->numberSites*lattice->wordsPerSite*sizeof(uint32_t))
        throw lm::InvalidArgException("IntLattice::nativeSerialize bufferSize", "the buffer size was not equal to the lattice data size");
    memcpy(destBuffer, lattice->particles, bufferSize);
}

void IntLattice::copyNativeToRowMajor(void * destBuffer, void * sourceBuffer, lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, uint particlesPerSite, size_t bufferSize)
{
    if (bufferSize != xSize*ySize*zSize*particlesPerSite*sizeof(uint32_t)) throw lm::InvalidArgException("IntLattice::copyNativeToRowMajor bufferSize", "the buffer size was not equal to the lattice data size");
    uint wordsPerSite = particlesPerSite;
    if (wordsPerSite*PARTICLES_PER_WORD != particlesPerSite) throw InvalidArgException("particlesPerSite", "must be evenly divisible by the number of particles per word");

    // Cast the buffers as appropriate.
    uint32_t * sourceParticles = (uint32_t *)sourceBuffer;
    uint32_t * destParticles = (uint32_t *)destBuffer;


    // Walk through the source buffer and copy the first word of particles.
    uint latticeIndex=0;
    for (uint w=0; w<wordsPerSite; w++)
    {
        for (uint z=0; z<zSize; z++)
        {
            for (uint y=0; y<ySize; y++)
            {
                for (uint x=0; x<xSize; x++, latticeIndex++)
                {
                    uint destIndex = x*ySize*zSize*particlesPerSite + y*zSize*particlesPerSite + z*particlesPerSite;
                    destParticles[destIndex+w*PARTICLES_PER_WORD] = sourceParticles[latticeIndex];
                }
            }
        }
    }
}

void IntLattice::copyRowMajorToNative(void * destBuffer, void * sourceBuffer, lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, uint particlesPerSite, size_t bufferSize)
{
	// buffer is of ints
    if (bufferSize != xSize*ySize*zSize*particlesPerSite) throw lm::InvalidArgException("IntLattice::copyRowMajorToNative bufferSize", "the buffer size was not equal to the lattice data size");
    uint wordsPerSite = particlesPerSite;
    if (wordsPerSite*PARTICLES_PER_WORD != particlesPerSite) throw InvalidArgException("particlesPerSite", "must be evenly divisible by the number of particles per word");

    // Cast the buffers as appropriate.
    uint32_t * destParticles = (uint32_t *)destBuffer;
    uint32_t * sourceParticles = (uint32_t *)sourceBuffer;
  

    // Walk through the dest buffer and copy into it the particles one word at a time.
    uint destIndex=0;
    for (uint w=0; w<wordsPerSite; w++)
    {
        for (uint z=0; z<zSize; z++)
        {
            for (uint y=0; y<ySize; y++)
            {
                for (uint x=0; x<xSize; x++, destIndex++)
                {
                    uint sourceIndex = x*ySize*zSize*particlesPerSite + y*zSize*particlesPerSite + z*particlesPerSite + w*PARTICLES_PER_WORD;
                    destParticles[destIndex] = sourceParticles[sourceIndex];
                }
            }
        }
    }
}

void IntLattice::copySitesRowMajorByteToNative(void * destBuffer, void * sourceBuffer, lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, size_t bufferSize)
{
    if (bufferSize != xSize*ySize*zSize) throw lm::InvalidArgException("IntLattice::copySitesRowMajorByteToNative bufferSize", "the buffer size was not equal to the lattice data size");

    // Cast the buffers as appropriate.
    uint8_t * destSites = (uint8_t *)destBuffer;
    uint8_t * sourceSites = (uint8_t *)sourceBuffer;


    // Walk through the dest buffer and copy into it sites.
    uint destIndex=0;
	for (uint z=0; z<zSize; z++)
	{
		for (uint y=0; y<ySize; y++)
		{
			for (uint x=0; x<xSize; x++, destIndex++)
			{
				uint sourceIndex = x*ySize*zSize + y*zSize + z;
				destSites[destIndex] = sourceSites[sourceIndex];
			}
		}
	}
}

// ONLY DEFINED BECAUSE IT IS PURE-VIRTUAL IN Lattice
void IntLattice::setFromRowMajorByteData(void * buffer, size_t bufferSize)
{
	throw("dont use this function");
    //copyRowMajorByteToNative(particles, buffer, size.x, size.y, size.z, wordsPerSite*PARTICLES_PER_WORD, bufferSize);
}
void IntLattice::setFromRowMajorData(void * buffer, size_t bufferSize)
{
    copyRowMajorToNative(particles, buffer, size.x, size.y, size.z, wordsPerSite*PARTICLES_PER_WORD, bufferSize);
}

void IntLattice::setSitesFromRowMajorByteData(void * buffer, size_t bufferSize)
{
    copySitesRowMajorByteToNative(siteTypes, buffer, size.x, size.y, size.z, bufferSize);
}

std::map<particle_t,uint> IntLattice::getParticleCounts()
{
    std::map<particle_t,uint> particleCountMap;
    for (lattice_size_t index=0; index<numberSites*wordsPerSite; index++)
    {
		unsigned int particle = particles[index];
		if(particle > 0)
		{
			if (particleCountMap.count(particle) == 0)
				particleCountMap[particle] = 1;
			else
				particleCountMap[particle]++;
        }
    }
    return particleCountMap;
}

std::vector<particle_loc_t> IntLattice::findParticles(particle_t minParticleType, particle_t maxParticleType)
{
    std::vector<particle_loc_t> ret;


    // Walk through the lattice.
    uint index=0;
    for (uint w=0; w<wordsPerSite; w++)
    {
        for (uint z=0; z<size.z; z++)
        {
            for (uint y=0; y<size.y; y++)
            {
                for (uint x=0; x<size.x; x++, index++)
                {
					if (particles[index] >= minParticleType && particles[index] <= maxParticleType)
						ret.push_back(particle_loc_t(particles[index], x, y, z, w*PARTICLES_PER_WORD));
                }
            }
        }
    }
    return ret;
}

uint32_t *IntLattice::getParticlesMemory()
{
	return particles;
}

uint8_t *IntLattice::getSitesMemory()
{
	return siteTypes;
}

void IntLattice::getSiteLatticeView(uint8_t **siteLattice, int *Nz, int *Ny, int *Nx) {
    *Nx = size.x;
    *Ny = size.y;
    *Nz = size.z;
    *siteLattice = siteTypes;
}

void IntLattice::getParticleLatticeView(uint32_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np) {
    *Nx = size.x;
    *Ny = size.y;
    *Nz = size.z;
    *Np = PARTICLES_PER_WORD; 
    *Nw = wordsPerSite;
    *particleLattice = (uint32_t*) particles;
}


}
}

