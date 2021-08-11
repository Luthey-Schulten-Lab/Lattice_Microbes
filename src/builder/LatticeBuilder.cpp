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

#include <list>
#include <cmath>
#include "config.h"
#include "core/Exceptions.h"
#include "core/Math.h"
#include "core/Print.h"
#include "builder/Cone.h"
#include "builder/Capsule.h"
#include "builder/CapsuleShell.h"
#include "builder/Cuboid.h"
#include "builder/Hemisphere.h"
#include "builder/LatticeBuilder.h"
#include "builder/Shape.h"
#include "builder/Sphere.h"
#include "builder/Difference.h"
#include "builder/Union.h"
#include "builder/UnionSet.h"
#include "builder/Intersection.h"
#include "builder/Torus.h"
#include "builder/Ellipse.h"
#include "builder/Cylinder.h"
#include "SpatialModel.pb.h"
#include "io/SimulationFile.h"
#include "rdme/Lattice.h"
#include "rng/RandomGenerator.h"
#include "rng/XORShift.h"
#ifdef OPT_CUDA
#include "rng/XORWow.h"
#include "cuda/lm_cuda.h"
#endif

using lm::io::SpatialModel;

namespace lm {
namespace builder {

LatticeBuilder::LatticeBuilder(si_dist_t xLen, si_dist_t yLen, si_dist_t zLen, si_dist_t collisionGridSpacing, uint32_t seedTop, uint32_t seedBottom)
:xLen(xLen),yLen(yLen),zLen(zLen),collisionGrid(NULL),collisionGridSize(NULL),collisionGridOccupancy(NULL),collisionGridSpacing(collisionGridSpacing),recipCollisionGridSpacing(1.0/collisionGridSpacing),rng(NULL)
{
    #ifdef OPT_CUDA
    // Create the cuda based rng, if we are using cuda and we have a cuda device assigned.
    //if (resources->cudaDevices.size() > 0)
   //{
		try
		{
			rng = new lm::rng::XORWow(0, seedTop, seedBottom, lm::rng::RandomGenerator::UNIFORM);
        Print::printf(Print::DEBUG, "Seeding xorwow rng with top word %u and bottom word %u", (unsigned int)(rng->getSeed()>>32), (unsigned int)(rng->getSeed()&0xFFFFFFFFLL));
		}
		catch(lm::CUDAException e)
		{
			Print::printf(Print::WARNING, "Initializing GPU random number generator failed: %s.  Falling back to CPU random number generator.", e.what());
			rng=NULL;
		}
    //}
    #endif

    if (rng == NULL)
    {
        rng = new lm::rng::XORShift(seedTop, seedBottom);
        Print::printf(Print::DEBUG, "Seeding xorshift rng with top word %u and bottom word %u", (unsigned int)(rng->getSeed()>>32), (unsigned int)(rng->getSeed()&0xFFFFFFFFLL));
    }

    // Allocate the collision detection grid.
    collisionGridXSize = (uint)lround(floor(xLen/collisionGridSpacing)+1);
    collisionGridYSize = (uint)lround(floor(yLen/collisionGridSpacing)+1);
    collisionGridZSize = (uint)lround(floor(zLen/collisionGridSpacing)+1);
    Print::printf(Print::DEBUG, "Creating collision grid of size %d %d %d", collisionGridXSize, collisionGridYSize, collisionGridZSize);
    collisionGrid = new Shape **[collisionGridXSize*collisionGridYSize*collisionGridZSize];
    collisionGridOccupancy = new uint[collisionGridXSize*collisionGridYSize*collisionGridZSize];
    collisionGridSize = new uint[collisionGridXSize*collisionGridYSize*collisionGridZSize];
    for (uint i=0; i<collisionGridXSize*collisionGridYSize*collisionGridZSize; i++)
    {
        collisionGrid[i] = new Shape*[2];
        collisionGridOccupancy[i] = 0U;
        collisionGridSize[i] = 2U;
    }
}

LatticeBuilder::~LatticeBuilder()
{
    // Delete the rng.
    if (rng != NULL) { delete rng; rng = NULL;}
    if (collisionGrid != NULL)
    {
        for (uint i=0; i<collisionGridXSize*collisionGridYSize*collisionGridZSize; i++)
        {
            if (collisionGrid[i] != NULL) {delete [] collisionGrid[i]; collisionGrid[i] = NULL;}
        }
        delete [] collisionGrid; collisionGrid = NULL;
    }
    if (collisionGridOccupancy != NULL) {delete collisionGridOccupancy; collisionGridOccupancy = NULL;}
    if (collisionGridSize != NULL) {delete collisionGridSize; collisionGridSize = NULL;}
}

void LatticeBuilder::addRegion(Shape * shape)
{
	// Make sure this region is in our defined list.
	bool found=false;
	for (uint i=0; i<definedRegions.size(); i++)
	{
		if (definedRegions[i] == shape->getType())
		{
			found = true;
			break;
		}
	}
	if (!found) definedRegions.push_back(shape->getType());

    // Add this region to the list.
    regionShapes[shape->getType()].push_back(shape);

    // Update the region's bounding box.
    if (regionBounds.count(shape->getType()) == 0)
        regionBounds[shape->getType()] = shape->getBoundingBox();
    else
        regionBounds[shape->getType()] = regionBounds[shape->getType()].joinWith(shape->getBoundingBox());

}

bool LatticeBuilder::placeObject(Shape * s)
{
    // Go through all of the grid volumes that this shape's bounding box intersects with.
    uint zMin = (uint)lround(floor(s->boundingBox.min.z*recipCollisionGridSpacing));
    uint zMax = (uint)lround(floor(s->boundingBox.max.z*recipCollisionGridSpacing));
    uint yMin = (uint)lround(floor(s->boundingBox.min.y*recipCollisionGridSpacing));
    uint yMax = (uint)lround(floor(s->boundingBox.max.y*recipCollisionGridSpacing));
    uint xMin = (uint)lround(floor(s->boundingBox.min.x*recipCollisionGridSpacing));
    uint xMax = (uint)lround(floor(s->boundingBox.max.x*recipCollisionGridSpacing));
    for (uint z=zMin; z<=zMax; z++)
    {
        for (uint y=yMin; y<=yMax; y++)
        {
            uint index=xMin+(y*collisionGridXSize)+(z*collisionGridXSize*collisionGridYSize);
            uint maxIndex = index + xMax - xMin;
            for (; index<=maxIndex; index++)
            {
                Print::printf(Print::VERBOSE_DEBUG, "Checking box at index %d",index);

                // Go through all of the shapes in this volume and see if we have a collision.
                uint occupancy = collisionGridOccupancy[index];
                for (uint i=0; i<occupancy; i++)
                {
                    Shape * shapeToCheck = collisionGrid[index][i];
                    if (shapeToCheck->intersects(s))
                    {
                        Print::printf(Print::VERBOSE_DEBUG, "Found another object intersecting.");
                        return false;
                    }
                }
            }
        }
    }

    // We need to add the shape to each grid volume that this shape's bounding box intersects with.
    for (uint z=zMin; z<=zMax; z++)
    {
        for (uint y=yMin; y<=yMax; y++)
        {
            uint index=xMin+(y*collisionGridXSize)+(z*collisionGridXSize*collisionGridYSize);
            for (uint x=xMin; x<=xMax; x++, index++)
            {
                // If there is no room in the list, increase its size.
                if (collisionGridOccupancy[index] >= collisionGridSize[index])
                {
                    uint currentSize = collisionGridSize[index];
                    Shape ** newList = new Shape*[currentSize*2];
                    memcpy(newList, collisionGrid[index], sizeof(Shape*)*currentSize);
                    collisionGridSize[index]=currentSize*2;
                    delete [] collisionGrid[index];
                    collisionGrid[index] = newList;
                }

                // Add the new shape to the list.
                int i=collisionGridOccupancy[index]++;
                collisionGrid[index][i] = s;
            }
        }
    }

    // Add the shape to the list of unique objects.
    objects.push_back(s);

    return true;
}

void LatticeBuilder::removeObject(Shape * s)
{
	// Add all of the obstacles.
    for (std::vector<Shape *>::iterator it=objects.begin(); it<objects.end(); it++)
    {
        Shape * shape = *it;
        if (s->getShapeType() == Shape::SPHERE && shape->getShapeType() == Shape::SPHERE)
        {
        	Sphere * target = (Sphere *)s;
            Sphere * sphere = (Sphere *)shape;
            if (target->getType() == sphere->getType() && target->getCenter().x == sphere->getCenter().x && target->getCenter().y == sphere->getCenter().y && target->getCenter().z == sphere->getCenter().z && target->getRadius() == sphere->getRadius())
            {
            	objects.erase(it);
            	return;
            }
        }
        else
        {
            throw Exception("Unsupported shape type.");
        }
    }
}

bool LatticeBuilder::placeSphere(point center, si_dist_t radius, site_t type)
{
	Sphere * s = new Sphere(center,radius,type);
	return placeObject(s);
}

void LatticeBuilder::removeSphere(point center, si_dist_t radius, site_t type)
{
	Sphere * s = new Sphere(center,radius,type);
	removeObject(s);
	delete s;
}

uint LatticeBuilder::placeRandomSphere(si_dist_t radius, site_t type, site_t region)
{
    uint numberAttempts=0;
    double randomValues[3];
    Sphere * s = new Sphere(0.0,radius,type);

    // Get the bounding box for this region.
    bounding_box b = regionBounds[region];

    // Randomly place the sphere until we find a location for it.
    while (true)
    {
        // Get a random location in the region's bounding box.
        numberAttempts++;
        rng->getRandomDoubles(randomValues, 3);
        point center(b.min.x+(((si_dist_t)randomValues[0])*(b.max.x-b.min.x)), b.min.y+(((si_dist_t)randomValues[1])*(b.max.y-b.min.y)), b.min.z+(((si_dist_t)randomValues[2])*(b.max.z-b.min.z)));
        s->setCenter(center);

        // Make sure the sphere is within one of the region's shapes.
        for (std::vector<Shape *>::iterator it=regionShapes[region].begin(); it<regionShapes[region].end(); it++)
        {
            if ((*it)->contains(s))
            {
            	// Make sure the object does not overlap with regions with a higher type value.
            	bool overlaps = false;
            	for (uint i=0; i<definedRegions.size() && overlaps == false; i++)
            	{
            		if (definedRegions[i] > region)
            		{
            			for (std::vector<Shape *>::iterator it2=regionShapes[definedRegions[i]].begin(); it2<regionShapes[definedRegions[i]].end(); it2++)
            			{
            				if ((*it2)->contains(s)) // TODO: this needs to be changes to intersects() after those methods are all implemented.
            				{
            					overlaps = true;
            					break;
            				}
            			}
            		}
            	}

                //  If it didn't overlap, try to place the object.
                if (!overlaps && placeObject(s))
                {
					Print::printf(Print::VERBOSE_DEBUG, "Placed sphere in region %d in %d attempts",region,numberAttempts);
					return numberAttempts;
                }
            }
        }
    }
}

void LatticeBuilder::placeRandomSpheres(uint count, si_dist_t radius, site_t type, site_t region)
{
    uint totalAttempts=0;
    for (uint i=0; i<count; i++)
    {
        totalAttempts+=placeRandomSphere(radius, type, region);
    }
    Print::printf(Print::DEBUG, "Placed %d spheres of type %d in region %d in %d total attempts",count,type,region,totalAttempts);
}

void LatticeBuilder::fillWithRandomSpheres(double volumeFraction, si_dist_t radius, site_t type, site_t region)
{
    double totalVolume = 0.0;
    std::vector<Shape *> shapes = regionShapes[region];
    for (std::vector<Shape *>::iterator it = shapes.begin(); it != shapes.end(); it++)
    {
        Shape * shape = *it;
        totalVolume += shape->getVolume();
    }

    uint numberSpheres = floor((totalVolume*volumeFraction)/((4.0/3.0)*PI*radius*radius*radius));
    Print::printf(Print::DEBUG, "Placing %d spheres of radius %0.2e in volume of %e m^3", numberSpheres, radius, totalVolume);
    placeRandomSpheres(numberSpheres, radius, type, region);

    if (Print::ifPrinted(9))  {
       uint totalGridObjects = 0;
       for (uint i=0; i<collisionGridXSize*collisionGridYSize*collisionGridZSize; i++)
	   totalGridObjects += collisionGridOccupancy[i];
       Print::printf(Print::DEBUG, "Number of grid volumes was %d with an mean of %0.2f objects each",collisionGridXSize*collisionGridYSize*collisionGridZSize,((double)totalGridObjects)/((double)collisionGridXSize*collisionGridYSize*collisionGridZSize));
    }
}

void LatticeBuilder::getSpatialModel(SpatialModel * spatialModel)
{
	spatialModel->Clear();

	// Set all of the regions.
	for (std::map<site_t,std::vector<Shape *> >::iterator it = regionShapes.begin(); it != regionShapes.end(); it++)
	{
		std::vector<Shape *> shapes = it->second;
		for (std::vector<Shape *>::iterator it2 = shapes.begin(); it2 != shapes.end(); it2++)
		{
			Shape * shape = *it2;
			spatialModel->add_region();
			uint i = spatialModel->region_size()-1;
			if (shape->getShapeType() == Shape::SPHERE) 
            {
				Sphere *sphere = (Sphere *)shape;
				spatialModel->mutable_region(i)->set_shape(sphere->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(sphere->getType());
				spatialModel->mutable_region(i)->add_shape_parameter(sphere->getCenter().x);
				spatialModel->mutable_region(i)->add_shape_parameter(sphere->getCenter().y);
				spatialModel->mutable_region(i)->add_shape_parameter(sphere->getCenter().z);
				spatialModel->mutable_region(i)->add_shape_parameter(sphere->getRadius());;
            }
            else if (shape->getShapeType() == Shape::HEMISPHERE)
            {
                Hemisphere * hemisphere = (Hemisphere *)shape;
                spatialModel->mutable_region(i)->set_shape(hemisphere->getShapeType());
                spatialModel->mutable_region(i)->set_site_type(hemisphere->getType());
                spatialModel->mutable_region(i)->add_shape_parameter(hemisphere->getCenter().x);
                spatialModel->mutable_region(i)->add_shape_parameter(hemisphere->getCenter().y);
                spatialModel->mutable_region(i)->add_shape_parameter(hemisphere->getCenter().z);
                spatialModel->mutable_region(i)->add_shape_parameter(hemisphere->getOrientation().x);
                spatialModel->mutable_region(i)->add_shape_parameter(hemisphere->getOrientation().y);
                spatialModel->mutable_region(i)->add_shape_parameter(hemisphere->getOrientation().z);
                spatialModel->mutable_region(i)->add_shape_parameter(hemisphere->getRadius());
            }
			else if (shape->getShapeType() == Shape::CAPSULE)
			{
				Capsule * capsule = (Capsule *)shape;
				spatialModel->mutable_region(i)->set_shape(capsule->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(capsule->getType());
				spatialModel->mutable_region(i)->add_shape_parameter(capsule->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(capsule->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(capsule->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(capsule->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(capsule->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(capsule->getP2().z);
				spatialModel->mutable_region(i)->add_shape_parameter(capsule->getRadius());
			}
			else if (shape->getShapeType() == Shape::CAPSULE_SHELL)
			{
				CapsuleShell * capsuleShell = (CapsuleShell *)shape;
				spatialModel->mutable_region(i)->set_shape(capsuleShell->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(capsuleShell->getType());
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getP2().z);
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getOuterRadius());
				spatialModel->mutable_region(i)->add_shape_parameter(capsuleShell->getInnerRadius());
			}
			else if (shape->getShapeType() == Shape::CUBOID)
			{
				Cuboid * cuboid = (Cuboid *)shape;
				spatialModel->mutable_region(i)->set_shape(cuboid->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(cuboid->getType());
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().z);
			}
			else if (shape->getShapeType() == Shape::UNION)
			{
				Union * uni = (Union *)shape;
				spatialModel->mutable_region(i)->set_shape(uni->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(uni->getType());
                // TODO: add union shape parameters
                /*
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().z);
                 */
			}
			else if (shape->getShapeType() == Shape::UNIONSET)
			{
				UnionSet * uni = (UnionSet *)shape;
				spatialModel->mutable_region(i)->set_shape(uni->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(uni->getType());
                // TODO: add union shape parameters
                /*
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().z);
                 */
			}
            else if (shape->getShapeType() == Shape::INTERSECTION)
			{
				Intersection * interObj = (Intersection *)shape;
				spatialModel->mutable_region(i)->set_shape(interObj->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(interObj->getType());
                // TODO: add intersection shape parameters
                /*
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cuboid->getP2().z);
                */
			}
            else if (shape->getShapeType() == Shape::DIFFERENCE)
			{
				Difference * diffObj = (Difference *)shape;
				spatialModel->mutable_region(i)->set_shape(diffObj->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(diffObj->getType());
                // TODO: add difference shape parameters
                /*
				spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().z);
                 */
			}
            else if (shape->getShapeType() == Shape::TORUS)
			{
				Torus * torObj = (Torus *)shape;
				spatialModel->mutable_region(i)->set_shape(torObj->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(torObj->getType());
                // TODO: add difference shape parameters
                /*
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().x);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().y);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().z);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().x);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().y);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().z);
                 */
			}
            else if (shape->getShapeType() == Shape::ELLIPSE)
			{
				Ellipse * ellObj = (Ellipse *)shape;
				spatialModel->mutable_region(i)->set_shape(ellObj->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(ellObj->getType());
                // TODO: add difference shape parameters
                /*
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().x);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().y);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP1().z);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().x);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().y);
                 spatialModel->mutable_region(i)->add_shape_parameter(diffObj->getP2().z);
                 */
			}
            else if (shape->getShapeType() == Shape::CYLINDER)
			{
				Cylinder * cylObj = (Cylinder *)shape;
				spatialModel->mutable_region(i)->set_shape(cylObj->getShapeType());
				spatialModel->mutable_region(i)->set_site_type(cylObj->getType());
				spatialModel->mutable_region(i)->add_shape_parameter(cylObj->getP1().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cylObj->getP1().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cylObj->getP1().z);
				spatialModel->mutable_region(i)->add_shape_parameter(cylObj->getP2().x);
				spatialModel->mutable_region(i)->add_shape_parameter(cylObj->getP2().y);
				spatialModel->mutable_region(i)->add_shape_parameter(cylObj->getP2().z);
				spatialModel->mutable_region(i)->add_shape_parameter(cylObj->getRadius());
			}
            else if (shape->getShapeType() == Shape::CONE)
            {
                Cone * coneObj = (Cone *)shape;
                spatialModel->mutable_region(i)->set_shape(coneObj->getShapeType());
                spatialModel->mutable_region(i)->set_site_type(coneObj->getType());
                spatialModel->mutable_region(i)->add_shape_parameter(coneObj->getCenter().x);
                spatialModel->mutable_region(i)->add_shape_parameter(coneObj->getCenter().y);
                spatialModel->mutable_region(i)->add_shape_parameter(coneObj->getCenter().z);
                spatialModel->mutable_region(i)->add_shape_parameter(coneObj->getRadius());
                spatialModel->mutable_region(i)->add_shape_parameter(coneObj->getHeight());
            }
			else
			{
				throw Exception("Unsupported shape type.");
			}
		}
	}

	// Add all of the obstacles.
    for (std::vector<Shape *>::iterator it=objects.begin(); it<objects.end(); it++)
    {
        Shape * shape = *it;
        spatialModel->add_obstacle();
        uint i = spatialModel->obstacle_size()-1;
        if (shape->getShapeType() == Shape::SPHERE)
        {
            Sphere * sphere = (Sphere *)shape;
            spatialModel->mutable_obstacle(i)->set_shape(sphere->getShapeType());
            spatialModel->mutable_obstacle(i)->set_site_type(sphere->getType());
            spatialModel->mutable_obstacle(i)->add_shape_parameter(sphere->getCenter().x);
            spatialModel->mutable_obstacle(i)->add_shape_parameter(sphere->getCenter().y);
            spatialModel->mutable_obstacle(i)->add_shape_parameter(sphere->getCenter().z);
            spatialModel->mutable_obstacle(i)->add_shape_parameter(sphere->getRadius());
        }
        else
        {
            throw Exception("Unsupported shape type.");
        }
    }
}

void LatticeBuilder::addParticles(particle_t particleType, site_t siteType, uint count)
{
	particles.push_back(ParticlePlacement(particleType, siteType, count));
}

void LatticeBuilder::discretizeTo(lm::rdme::Lattice * lattice, site_t obstacleSiteType, double fractionObstacleSitesOccupied)
{
	// Go through each region and discretize it to the lattice.
	for (uint i=0; i<definedRegions.size(); i++)
	{
		site_t region = definedRegions[i];

		// Go through each shape making up the region.
		for (std::vector<Shape *>::iterator it=regionShapes[region].begin(); it<regionShapes[region].end(); it++)
		{
			Shape * shape = *it;
			shape->discretizeTo(lattice);
		}
	}

	// Add the obstacles.
	this->discretizeObstaclesTo(lattice, obstacleSiteType, fractionObstacleSitesOccupied);

	// Add all of the particles to the lattice.
	for (uint i=0; i<particles.size(); i++)
	{
		ParticlePlacement p = particles[i];
		uint numberAttempts=0;
		uint placed=0;
		while (placed < p.count)
		{
			numberAttempts++;
			double randomValues[3];
			rng->getRandomDoubles(randomValues, 3);
			lattice_size_t x = (uint)floor(randomValues[0]*(double)lattice->getXSize());
			if (x == lattice->getXSize()) x--;
			lattice_size_t y = (uint)floor(randomValues[1]*(double)lattice->getYSize());
			if (y == lattice->getYSize()) y--;
			lattice_size_t z = (uint)floor(randomValues[2]*(double)lattice->getZSize());
			if (z == lattice->getZSize()) z--;
			if (lattice->getSiteType(x, y, z) == p.siteType && lattice->getOccupancy(x,y,z) < lattice->getMaxOccupancy())
			{
				lattice->addParticle(x, y, z, p.particleType+1);
				placed++;
			}
		}
		Print::printf(Print::DEBUG, "Placed %d particles of type %d in sites of type %d in %d attempts",p.count,p.particleType,p.siteType,numberAttempts);
	}
}

void LatticeBuilder::discretizeObstaclesTo(lm::rdme::Lattice * lattice, site_t obstacleSiteType, double fractionObstacleSitesOccupied)
{
	// Track the occupied volume of each subvolume.
	double * occupiedVolume = new double[lattice->getNumberSites()];
	for (uint i=0; i<lattice->getNumberSites(); i++) occupiedVolume[i] = 0.0;

	// Go through each obstacle and assign its volume to the appropriate subvolume.
    for (std::vector<Shape *>::iterator it=objects.begin(); it<objects.end(); it++)
    {
        Shape * shape = *it;
        if (shape->getShapeType() == Shape::SPHERE)
        {
            Sphere * sphere = (Sphere *)shape;

            // If the sphere is much bigger than the subvolume size, occupy all of the sites with which the sphere overlaps.
            if (sphere->getRadius() > lattice->getSpacing())
            {
            	point c=sphere->getCenter();
            	si_dist_t r=sphere->getRadius();
            	si_dist_t siteSubvolume = lattice->getSpacing()*lattice->getSpacing()*lattice->getSpacing();

        		// Get the bounds of the circle.
        		lattice_size_t x1 = (lattice_size_t)floor((c.x-r)/lattice->getSpacing());
        		lattice_size_t x2 = ((lattice_size_t)floor((c.x+r)/lattice->getSpacing()))+1;
        		lattice_size_t y1 = (lattice_size_t)floor((c.y-r)/lattice->getSpacing());
        		lattice_size_t y2 = ((lattice_size_t)floor((c.y+r)/lattice->getSpacing()))+1;
        		lattice_size_t z1 = (lattice_size_t)floor((c.z-r)/lattice->getSpacing());
        		lattice_size_t z2 = ((lattice_size_t)floor((c.z+r)/lattice->getSpacing()))+1;

        		// Fill in the circle.
        		for (lattice_size_t x = x1; x<=x2; x++)
        			for (lattice_size_t y = y1; y<=y2; y++)
            			for (lattice_size_t z = z1; z<=z2; z++)
						{
							point lc((((double)x)+0.5)*lattice->getSpacing()-c.x, (((double)y)+0.5)*lattice->getSpacing()-c.y, (((double)z)+0.5)*lattice->getSpacing()-c.z);
							if (lc.x*lc.x+lc.y*lc.y+lc.z*lc.z < r*r)
								occupiedVolume[x*lattice->getYSize()*lattice->getZSize() + y*lattice->getZSize() + z] += siteSubvolume;
						}
            }

            // Otherwise, the sphere is on the order of the subvolume size, try to distribute its volume to nearby subvolumes.
            else
            {
				// Figure out which half lattice site the subvolume is nearest to.
				uint hx = lround(sphere->getCenter().x/(lattice->getSpacing()/2));
				uint hy = lround(sphere->getCenter().y/(lattice->getSpacing()/2));
				uint hz = lround(sphere->getCenter().z/(lattice->getSpacing()/2));

				uint xsv[2];
				uint ysv[2];
				uint zsv[2];

				// Figure out which subvolumes to assign the density.
				if (hx%2 == 1)
				{
					// x falls in the center of a subvolume.
					xsv[0] = hx/2;
					xsv[1] = xsv[0];
				}
				else
				{
					// x falls on the edge of two subvolumes.
					xsv[0] = hx/2;
					xsv[1] = xsv[0]-1;
				}
				if (hy%2 == 1)
				{
					// y falls in the center of a subvolume.
					ysv[0] = hy/2;
					ysv[1] = ysv[0];
				}
				else
				{
					// y falls on the edge of two subvolumes.
					ysv[0] = hy/2;
					ysv[1] = ysv[0]-1;
				}
				if (hz%2 == 1)
				{
					// z falls in the center of a subvolume.
					zsv[0] = hz/2;
					zsv[1] = zsv[0];
				}
				else
				{
					// z falls on the edge of two subvolumes.
					zsv[0] = hz/2;
					zsv[1] = zsv[0]-1;
				}

				// Puts the density in the subvolumes.
				double eighthVolume = sphere->getVolume()/8;
				for (uint i=0; i<2; i++)
					for (uint j=0; j<2; j++)
						for (uint k=0; k<2; k++)
						{
							occupiedVolume[xsv[i]*lattice->getYSize()*lattice->getZSize() + ysv[j]*lattice->getZSize() + zsv[k]] += eighthVolume;
						}
            }
        }
        else
        {
            throw Exception("Unsupported shape type.");
        }
    }

    // Create a sorted list of all of the occupancy values.
    std::list<double> occupancyValues;
    for (uint x=0, index=0; x<lattice->getXSize(); x++)
    	for (uint y=0; y<lattice->getYSize(); y++)
    		for (uint z=0; z<lattice->getZSize(); z++, index++)
    			if (occupiedVolume[index] > 0.0)
    				occupancyValues.push_back(occupiedVolume[index]);
    occupancyValues.sort();
    occupancyValues.reverse();

    // Get the occupancy cutoff.
    uint maxOccupiedSites = lround(fractionObstacleSitesOccupied*((double)lattice->getNumberSites()));
    double cutoff = 0.0;
    uint i=0;
    for (std::list<double>::iterator it = occupancyValues.begin(); it != occupancyValues.end(); it++, i++)
    {
    	double occupancy = *it;
    	if (i >= maxOccupiedSites)
    	{
    		cutoff = occupancy;
    		break;
    	}
    }

    // Go through the subvolumes and mark any subvolumes that are over the threshold as occupied.
    uint numberOccupiedSites=0;
    for (uint x=0, index=0; x<lattice->getXSize(); x++)
    	for (uint y=0; y<lattice->getYSize(); y++)
    		for (uint z=0; z<lattice->getZSize(); z++, index++)
    		{
    			if (occupiedVolume[index] > 0.0 && occupiedVolume[index] >= cutoff)
    			{
    				lattice->setSiteType(x, y, z, obstacleSiteType);
    				numberOccupiedSites++;
    			}
    		}

    // Free the resources.
    delete [] occupiedVolume;

    Print::printf(Print::DEBUG, "Marked %d sites as occupied by obstacles",numberOccupiedSites);
}

}
}

