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
 * Overflow algorithm in RDME solvers and CPU assignment (2012)
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

#ifndef LM_RDME_LATTICE_H_
#define LM_RDME_LATTICE_H_

//using namespace std; //TODO: workaround for cuda compiler bug, remove this line when fixed.
#include <iostream>
#include <vector>
#include <map>
#include "core/Types.h"
#include "core/Exceptions.h"

// Type to store a lattice index.
typedef uint32_t            lattice_size_t;

/// @struct lattice_coord_t
/// @brief Type to store a lattice coordinate.
struct lattice_coord_t {
    /// @brief Create a lattice coordinate
    /// @param x Lattice x point
    /// @param y Lattice y point
    /// @param z Lattice z point
    lattice_coord_t(lattice_size_t x=0, lattice_size_t y=0, lattice_size_t z=0):x(x),y(y),z(z){}
    lattice_size_t x;
    lattice_size_t y;
    lattice_size_t z;
};

// Type to store a particle index in a site.
typedef uint32_t            site_size_t;

// Type to store a lattice site type.
typedef uint32_t            site_t;

// Type to store a lattice particle type.
typedef uint32_t            particle_t;

/// @struct particle_loc_t
/// @brief Type to store a particle and position.
struct particle_loc_t {
    /// @brief Create a particle location
    /// @param p Particle type at the location
    /// @param x Lattice x point
    /// @param y Lattice y point
    /// @param z Lattice z point
    /// @param index Lattice index
    particle_loc_t(particle_t p=0, lattice_size_t x=0, lattice_size_t y=0, lattice_size_t z=0, site_size_t index=0):p(p),x(x),y(y),z(z),index(index){}
    particle_t p;
    
    // Location in the grid and index
    lattice_size_t x;
    lattice_size_t y;
    lattice_size_t z;
    site_size_t index;
};

namespace lm {
namespace rdme {

unsigned int getCompiledLatticeMaxOccupancy();

/// @class InvalidSiteException
/// @brief InvalidArgException for when an index into the lattice is out of bound or does not exist.
class InvalidSiteException : public InvalidArgException
{
public:
    /// @brief Create an exception based on the lattice location
    /// @param x Lattice x point
    /// @param y Lattice y point
    /// @param z Lattice z Point
    InvalidSiteException(lattice_size_t x, lattice_size_t y, lattice_size_t z) : InvalidArgException("lattice site index") {}
    /// @brief Create an exception based on the lattice index
    /// @param index Lattice index
    InvalidSiteException(lattice_size_t index) : InvalidArgException("lattice site index") {}
};

/// @class InvalidParticleException
/// @brief InvalidArgException for when a particle is 
class InvalidParticleException : public InvalidArgException
{
public:
    /// @brief Create an eception based on the particle index
    /// @param particleIndex Index of the particle
    InvalidParticleException(site_size_t particleIndex) : InvalidArgException("particleIndex", "Invalid particle index", particleIndex) {}
};

    
/// @class Lattice
/// @brief Base class for lattice type objects
class Lattice
{
public:
    static void rowMajorByteSerialize(void * destBuffer, void * lattice, size_t bufferSize);
    static void rowMajorIntSerialize(void * destBuffer, void * lattice, size_t bufferSize);
    static void rowMajorByteSerializeSites(void * destBuffer, void * lattice, size_t bufferSize);

public:
    // Lattice limits.
    /// @brief Get the maximum number of site types possible in the lattice
    virtual site_t getMaxSiteType() const =0;
    /// @brief Get the maximum number of particle types possible in the lattice
    virtual particle_t getMaxParticle() const =0;
    /// @brief Get the maximum number of particles that can live in a site
    virtual site_size_t getMaxOccupancy() const =0;
	
public:
    /// @brief Create a Lattice object
    /// @param size Size of the lattice as (x,y,z)
    /// @param spacing Physical spacing between lattice sites
    Lattice(lattice_coord_t size, si_dist_t spacing);
    /// @brief Create a Lattice object
    /// @param xSize Number of sites in x
    /// @param ySize Number of sites in y
    /// @param zSize Number of sites in z
    /// @param spacing Physical spacing between lattice sites
	Lattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing);
    /// @brief Destroy the Lattice object
	virtual ~Lattice();
    
    // Size related methods
    /// @brief Get size of the Lattice
	virtual lattice_coord_t getSize() const;
    /// @brief Get x dimension of the Lattice
	virtual lattice_size_t getXSize() const;
    /// @brief Get y dimension of the Lattice
	virtual lattice_size_t getYSize() const;
    /// @brief Get z dimension of the Lattice
	virtual lattice_size_t getZSize() const;
    /// @brief Get total number of sites in the Lattice
	virtual lattice_size_t getNumberSites() const;
    /// @brief Get spacing between lattice sites
	virtual si_dist_t getSpacing() const;
	
    
    /// @brief Get the sites that are neighbor to the indicated site
    /// @param index Index of the site for which to get neighbors
    /// @param neighboringIndicies An array to hold the indicies of the neighbor sites
	virtual void getNeighboringSites(lattice_size_t index, lattice_size_t * neighboringIndices)=0;

    
	// Lattice site methods.
    /// @brief Get the site type at the specified location
	virtual site_t getSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z) const = 0;
    /// @brief Get the site type at the specified location
	virtual site_t getSiteType(lattice_size_t index) const = 0;
    /// @brief Set the site type at the specified location
	virtual void setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t site) = 0;
    /// @brief Set the site type at the specified location
	virtual void setSiteType(lattice_size_t index, site_t site) = 0;
    /// @brief Get a list of sites near the specified site within a certain distance
    /// @param xc X location of the center site
    /// @param yc Y location of the center site
    /// @param zc Z location of the center site
    /// @param minDistance Minimum distance for halo of sites
    /// @param maxDistance Maximum distance for halo of sites
	std::vector<lattice_coord_t> getNearbySites(lattice_size_t xc, lattice_size_t yc, lattice_size_t zc, uint minDistance, uint maxDistance);

//	Need to remove this from the base class, can't assume a datatype for particlelattice
//        virtual void getParticleLatticeView(uint8_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np) = 0;
//        virtual void getSiteLatticeView(uint8_t **siteLattice, int *Nz, int *Ny, int *Nx) = 0;

    
	// Particle methods.
    /// @brief Get the number of particles in the specified lattice site
	virtual site_size_t getOccupancy(lattice_size_t x, lattice_size_t y, lattice_size_t z) const = 0;
    /// @brief Get the number of particles in the specified lattice site
	virtual site_size_t getOccupancy(lattice_size_t index) const = 0;
    
    /// @brief Get the particle at the specified site with at the specified number in the particle list
	virtual particle_t getParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_size_t particleIndex) const = 0;
    /// @brief Get the particle at the specified site with at the specified number in the particle list
	virtual particle_t getParticle(lattice_size_t index, site_size_t particleIndex) const = 0;
    
    /// @brief Add a particle to the specified site
	virtual void addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle) = 0;
    /// @brief Add a particle to the specified site
	virtual void addParticle(lattice_size_t index, particle_t particle) = 0;
    /// @brief Remove a particle to the specified site
    virtual void removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z) = 0;
    /// @brief Remove a particle to the specified site
    virtual void removeParticles(lattice_size_t index) = 0;
    /// @brief Empty all particles from the specified site
	virtual void removeAllParticles();
	
	/**
	 * Particle searching methods.
	 */

	// Search for a particle from the start of the lattice.
	/*virtual particle_loc_t findParticle(particle_t particle)=0;
	virtual particle_loc_t findNextParticle(particle_loc_t previousParticle)=0;
	
	// Search for a particle starting from a given location.
	virtual particle_loc_t findNearbyParticle(particle_loc_t particle)=0;
	*/
	
	// Counts all of the particles on the lattice.
    /// @brief Get the number of each particle type in the lattice
	virtual std::map<particle_t,uint> getParticleCounts()=0;
	
	//Finds all of the particle of a given range of types.
    /// @brief Get the number of the specified particles types in the lattice
	virtual std::vector<particle_loc_t> findParticles(particle_t minParticleType, particle_t maxParticleType)=0;

    /// @brief Print the lattice to the console
	virtual void print() const;

	// Methods to set the data directly.
	virtual void setFromRowMajorByteData(void * buffer, size_t bufferSize)=0;
	virtual void setSitesFromRowMajorByteData(void * buffer, size_t bufferSize)=0;

	virtual size_t getLatticeMemorySize() const = 0;

protected:
	lattice_coord_t size;           // Dimension of lattice as (x, y, z) 
	lattice_size_t numberSites;     // Number of total sites (e.g. x*y*z)
    si_dist_t spacing;              // Physical spacing of the lattice (i.e. in units of nm, um, mm, etc.)
};


}
}

#endif
