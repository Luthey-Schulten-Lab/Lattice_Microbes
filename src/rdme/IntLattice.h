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

#ifndef LM_RDME_INTLATTICE_H_
#define LM_RDME_INTLATTICE_H_

#include <map>
#include "core/Types.h"
#include "rdme/Lattice.h"

namespace lm {
namespace rdme {

/// @class IntLattice
/// @brief A Lattice that is based on one particle per word, with sites strided per particle
class IntLattice : public Lattice
{
	// TODO: total hack.  Needed to access the actual data structures
	// Should modify CudaByteLattice to do what I need.
	friend class MGPUMpdRdmeSolver;
	friend class MPIMpdRdmeSolver;
        friend class MGPUIntMpdRdmeSolver;

public:
    static void nativeSerialize(void * destBuffer, void * lattice, size_t latticeSize);
    static void copyNativeToRowMajor(void * destBuffer, void * sourceBuffer, lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, uint particlesPerSite, size_t bufferSize);
    static void copyRowMajorToNative(void * destBuffer, void * sourceBuffer, lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, uint particlesPerSite, size_t bufferSize);
    static void copySitesRowMajorByteToNative(void * destBuffer, void * sourceBuffer, lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, size_t bufferSize);

public:
    // Lattice limits.
    virtual site_t getMaxSiteType() const;
    virtual particle_t getMaxParticle() const;
    virtual site_size_t getMaxOccupancy() const;

public:
    IntLattice(lattice_coord_t size, si_dist_t spacing, uint particlesPerSite);
    IntLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing, uint particlesPerSite);
	virtual ~IntLattice();
	
	virtual void getNeighboringSites(lattice_size_t index, lattice_size_t * neighboringIndices);

	// Lattice site methods.
	virtual site_t getSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z) const;
	virtual site_t getSiteType(lattice_size_t subvolume) const;
	virtual void setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t siteType);
	virtual void setSiteType(lattice_size_t subvolume, site_t site);
	
	// Particle methods.
	virtual site_size_t getOccupancy(lattice_size_t x, lattice_size_t y, lattice_size_t z) const;
	virtual site_size_t getOccupancy(lattice_size_t subvolume) const;
	virtual particle_t getParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_size_t particleIndex) const;
	virtual particle_t getParticle(lattice_size_t subvolume, site_size_t particleIndex) const;
	virtual void addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle);
	virtual void addParticle(lattice_size_t subvolume, particle_t particle);
	virtual void removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z);
    virtual void removeParticles(lattice_size_t subvolume);
    virtual void removeAllParticles();
	
	// Particle searching.
    /*virtual particle_loc_t findParticle(particle_t particle);
    virtual particle_loc_t findNextParticle(particle_loc_t previousParticle);
    virtual particle_loc_t findNearbyParticle(particle_loc_t particle)=0;
	*/

    virtual std::map<particle_t,uint> getParticleCounts();
    virtual std::vector<particle_loc_t> findParticles(particle_t minParticleType, particle_t maxParticleType);
	
    // Methods to set the data directly.
    virtual void setFromRowMajorByteData(void * buffer, size_t bufferSize);
    virtual void setFromRowMajorData(void * buffer, size_t bufferSize);
    virtual void setSitesFromRowMajorByteData(void * buffer, size_t bufferSize);

    virtual void getParticleLatticeView(uint32_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np);
    virtual void getSiteLatticeView(uint8_t **siteLattice, int *Nz, int *Ny, int *Nx);
	
	virtual size_t getLatticeMemorySize() const;

protected:
	virtual void allocateMemory();
	virtual void deallocateMemory();
	virtual uint32_t* getParticlesMemory();
	virtual uint8_t* getSitesMemory();

protected:
    uint wordsPerSite;      // Number of bytes per site
    uint32_t * particles;   // Array of lists of particles at each site
    uint8_t * siteTypes;    // Array of site types representing the lattice

private:
    static const uint PARTICLES_PER_WORD = 1;
};

}
}

#endif
