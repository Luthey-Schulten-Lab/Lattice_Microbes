/*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
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
 * Author(s): Andrew Magis, Elijah Roberts
 */

#ifndef LM_BUILDER_LATTICEBUILDER_H_
#define LM_BUILDER_LATTICEBUILDER_H_

#include <map>
#include <utility>
#include <vector>
#include "core/Types.h"
#include "builder/Shape.h"
#include "rdme/Lattice.h"
#include "rng/RandomGenerator.h"

using std::vector;
using lm::rng::RandomGenerator;

namespace lm {

namespace io {
class SpatialModel;
}

namespace rdme {
class Lattice;
}

namespace builder {

/// @class LatticeBuilder
/// @brief A class that defines regions of a lattice based on a set of geometries defined by shapes.  It also allows packing different types of particles into different regions. 
class LatticeBuilder
{
protected:
    /// @struct ParticlePlacement
    /// @brief A representation of a placement operation that needs to occur (i.e. where to place (siteType), what to place (particleType) and how many (count))
	struct ParticlePlacement
	{
    	particle_t particleType;
    	site_t siteType;
    	uint count;
        /// @brief Create a ParticlePlacement placement operation
        /// @param particleType The type of particle to place
        /// @param siteType The type of site in which to place the particles
        /// @param count The number of particles to place in that site type
    	ParticlePlacement(particle_t particleType, site_t siteType, uint count):particleType(particleType),siteType(siteType),count(count) {}
	};

public:
    /// @brief Create a Lattice Builder
    /// @param xLen Length of the domain along x-axis
    /// @param yLen Length of the domain along y-axis
    /// @param zLen Length of the domain along z-axis
    /// @param collisionGridSpacing The spacing for collision objects
    /// @param seedTop High 32 bits of the seed (allows a constant seed for debugging)
    /// @param seedBottom Low 32 bits of the seed
    LatticeBuilder(si_dist_t xLen, si_dist_t yLen, si_dist_t zLen, si_dist_t collisionGridSpacing, uint32_t seedTop, uint32_t seedBottom=0);
    virtual ~LatticeBuilder();

    /// @brief Add a region to the lattice
    /// @param shape A Shape object to add as a region
    virtual void addRegion(Shape * shape);
    /// @brief Add an shape to the lattice
    /// @param shape A Shape object to add as a region
    /// @return true if the object to place does not intersect another object
    virtual bool placeObject(Shape * shape);
    /// @brief Remove the shape from the lattice
    /// @brief s Shape to remove
    virtual void removeObject(Shape * s);
    
    /// @brief Place a sphere in the lattice (for obstacles)
    /// @param center The center point of sphere obstacle
    /// @param radius Radius of the sphere obstacle
    /// @param type The type of site in which to place sphere
    /// @return true if the sphere did not intersect
    virtual bool placeSphere(point center, si_dist_t radius, site_t type);
    /// @brief Remove a sphere in the lattice (for obstacles)
    /// @param center The center point of sphere obstacle
    /// @param radius Radius of the sphere obstacle
    /// @param type The type of site in which to place sphere
    virtual void removeSphere(point center, si_dist_t radius, site_t type);
    /// @brief Place a sphere randomly in the lattice (for obstacles)
    /// @param radius Radius of the sphere obstacle
    /// @param type The type of site in which to place sphere
    /// @param region The region in which to place obstacle randomly
    /// @return number of times a placement operation occured
    virtual uint placeRandomSphere(si_dist_t radius, site_t type, site_t region);
    /// @brief Place many spheres randomly in the lattice (for obstacles)
    /// @param radius Radius of the sphere obstacle
    /// @param type The type of site in which to place sphere
    /// @param region The region in which to place obstacle randomly
    /// @return number of times a placement operation occured
    virtual void placeRandomSpheres(uint count, si_dist_t radius, site_t type, site_t region);
    /// @brief Fill a region with random spheres to a specified volume fraction
    /// @param volumeFraction Total fraction of volume that should be filled with spheres
    /// @param radius Radius of spheres
    /// @param type The type of site to fill (i.e. the type of site to exclude other objects from)
    /// @param region The region of the lattice to place spheres into
    virtual void fillWithRandomSpheres(double volumeFraction, si_dist_t radius, site_t type, site_t region);
    /// @brief Gets a spatial model of the lattice for interface with python.  NOTE: this operation clears the object passed in from python
    /// @param spatialModel An object of a spatial model for interaction in python or HDF5.  The model will be filled with the current lattice
    virtual void getSpatialModel(lm::io::SpatialModel * spatialModel);

    
    /// @brief Add particles of a given type
    /// @param particleType The type of particles to randomly place in the lattice
    /// @param siteType Type of lattice site into which to place
    /// @param count Number of particles to place
    virtual void addParticles(particle_t particleType, site_t siteType, uint count);

    /// @brief Discretizes the regions to a square lattice
    /// @param lattice Lattice object into which to place
    /// @param obstacleSiteType An identifier for obstacle sites in the lattice
    /// @param fractionObstacleSitesOccupied Percentage of obstacle sites to be filled
    virtual void discretizeTo(lm::rdme::Lattice * lattice, site_t obstacleSiteType, double fractionObstacleSitesOccupied);

protected:
    /// @brief Discretizes the obstacles to a square lattice
    /// @param lattice Lattice object into which to place
    /// @param obstacleSiteType An identifier for obstacle sites in the lattice
    /// @param fractionObstacleSitesOccupied Percentage of obstacle sites to be filled
    virtual void discretizeObstaclesTo(lm::rdme::Lattice * lattice, site_t obstacleSiteType, double fractionObstacleSitesOccupied);

protected:
    // Lattice dimensions
    si_dist_t xLen, yLen, zLen;
    
    // Regions of the lattice (i.e. cell membrane, cytosol, extracellular space, etc.)
    std::vector<site_t> definedRegions;
    std::map<site_t,std::vector<Shape *> > regionShapes;
    std::map<site_t,bounding_box> regionBounds;
    std::vector<Shape *> objects;
    
    // Particle placement operations or in other words how many particles of each type that should be placed and where
    std::vector<ParticlePlacement> particles;
    
    // Collision objects, collision grid size and percent filled
    Shape *** collisionGrid;
    uint * collisionGridSize;
    uint * collisionGridOccupancy;
    si_dist_t collisionGridSpacing;
    double recipCollisionGridSpacing;
    uint collisionGridXSize, collisionGridYSize, collisionGridZSize;
    
    // Random number generator associated with this lattice
    RandomGenerator * rng;

};

}
}

#endif
