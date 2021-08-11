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
 * Author(s): Elijah Roberts
 */

#include "config.h"
#include <cmath>
#include "builder/Shape.h"
#include "builder/Sphere.h"
#include "builder/CapsuleShell.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {

CapsuleShell::CapsuleShell(point _p1, point _p2, si_dist_t _innerRadius, si_dist_t _outerRadius, site_t _type)
:Shape(CAPSULE_SHELL
      ,CapsuleShell::calcBoundingBox(_p1,_p2,_outerRadius)
      ,_type
      ,vector(_p1,_p2).unitize()
     ,vector::findPerpendicularVector(vector(_p1,_p2)).unitize())
,p1(_p1)
,p2(_p2)
,innerRadius(_innerRadius)
,outerRadius(_outerRadius)
,length(_p1.distance(_p2)+_outerRadius+_outerRadius)
{
}

bounding_box CapsuleShell::calcBoundingBox(point p1, point p2, si_dist_t radius)
{
	bounding_box bbox;
    double adj = max(radius, p1.distance(p2)+2.0*radius);
	bbox.min.x = (p1.x<p2.x)?(p1.x-adj):(p2.x-adj);
	bbox.max.x = (p1.x>p2.x)?(p1.x+adj):(p2.x+adj);
	bbox.min.y = (p1.y<p2.y)?(p1.y-adj):(p2.y-adj);
	bbox.max.y = (p1.y>p2.y)?(p1.y+adj):(p2.y+adj);
	bbox.min.z = (p1.z<p2.z)?(p1.z-adj):(p2.z-adj);
	bbox.max.z = (p1.z>p2.z)?(p1.z+adj):(p2.z+adj);
	return bbox;
}

CapsuleShell::~CapsuleShell()
{
}

bool CapsuleShell::intersects(Shape * query)
{
    // FIXME: NEEDS TO BE IMPLEMENTED //
    if (!boundingBoxesIntersect(query)) return false;

    if (query->getShapeType() == SPHERE)
    {
    	return false;
    }
    return false;
}

bool CapsuleShell::contains(point query)
{
    // Check if point is inside
    si_dist_t oradius2 = outerRadius*outerRadius;
    si_dist_t iradius2 = innerRadius*innerRadius;
    // 1) Compute projection onto normal to get length
    double proj = vector(p1, query).dot(at);
    // 2) Check if "below" base or "above" apex
    if(proj < 0.0 || proj > p1.distance(p2)) {
        if (query.distanceSquared(p1) <= oradius2 && query.distanceSquared(p1) >= iradius2) return true;
        if (query.distanceSquared(p2) <= oradius2 && query.distanceSquared(p2) >= iradius2) return true;
        return false;
    }
    // 3) Find scaled center (e.g. how far up the center spine)
    vector h_prime = at.scale(proj);
    vector c_prime = vector(p1.plus(h_prime.toPoint()));
    // 4) compute distance to scaled center
    if(c_prime.toPoint().distanceSquared(query) <= oradius2 && c_prime.toPoint().distanceSquared(query) >= iradius2)
        return true;

    return false;
}

bool CapsuleShell::contains(Shape * query)
{
	if (!boundingBoxesIntersect(query)) return false;


    if (query->getShapeType() == SPHERE)
    {
        Sphere * querySphere = (Sphere *)query;

        // Transform the sphere as if the capsule was oriented in the +z direction with p1 at 0,0,r and p2 at 0,0,l-r.
        point tc = querySphere->getCenter();
        si_dist_t tr=querySphere->getRadius();

        // Assume for now that the capsule is oriented with its long axis along the z axis and p1.z = r and p2.z = l-r.

        // Translate the sphere in the x-y plane.
        tc.x -= p1.x;
        tc.y -= p1.y;
        tc.z -= (p1.z-outerRadius);

        // See if the sphere is closer to the lower cap, the cylinder, or the upper cap.
        if (tc.z < outerRadius)
        {
        	double rSquared = tc.x*tc.x + tc.y*tc.y + (tc.z-(outerRadius))*(tc.z-(outerRadius));
            if (rSquared > (outerRadius-tr)*(outerRadius-tr) || rSquared < (innerRadius+tr)*(innerRadius+tr)) return false;
            return true;
        }
        else if (tc.z > length-outerRadius)
        {
        	double rSquared = tc.x*tc.x+tc.y*tc.y+(tc.z-(length-outerRadius))*(tc.z-(length-outerRadius));
        	if (rSquared > (outerRadius-tr)*(outerRadius-tr) || rSquared < (innerRadius+tr)*(innerRadius+tr)) return false;
            return true;
        }
        else
        {
            // Make sure the sphere does not fall outside the cylinder.
        	double rSquared = tc.x*tc.x+tc.y*tc.y;
            if (rSquared > (outerRadius-tr)*(outerRadius-tr) || rSquared < (innerRadius+tr)*(innerRadius+tr)) return false;
            return true;
        }
    }
    return false;
}

double CapsuleShell::getVolume()
{
	double vouter = ((4.0/3.0)*PI*outerRadius*outerRadius*outerRadius)+(PI*outerRadius*outerRadius*(length-outerRadius-outerRadius));
	double vinner = ((4.0/3.0)*PI*innerRadius*innerRadius*innerRadius)+(PI*innerRadius*innerRadius*(length-innerRadius-innerRadius));
    return vouter-vinner;
}

//void CapsuleShell::discretizeTo(lm::rdme::Lattice * lattice)
//{
//	// Assume that the capsule is oriented in the +z direction.
//
//	// Get the starting and ending z planes.
//	lattice_size_t z1 = (lattice_size_t)floor((p1.z-outerRadius+lattice->getSpacing()/2)/lattice->getSpacing());
//	lattice_size_t z2 = ((lattice_size_t)floor((p2.z+outerRadius-lattice->getSpacing()/2)/lattice->getSpacing()))+1;
//
//	// The first z plane is a solid circle.
//	{
//		lattice_size_t z = z1;
//
//		// Get the distance to p1.
//		si_dist_t b = p1.z-(((double)z)+0.5)*lattice->getSpacing();
//
//		// Get the radius of the circle in this plane.
//		si_dist_t pr = sqrt(outerRadius*outerRadius-b*b);
//
//		// Get the bounds of the circle in this plane.
//		lattice_size_t x1 = (lattice_size_t)floor((p1.x-pr)/lattice->getSpacing());
//		lattice_size_t x2 = ((lattice_size_t)floor((p1.x+pr)/lattice->getSpacing()))+1;
//		lattice_size_t y1 = (lattice_size_t)floor((p1.y-pr)/lattice->getSpacing());
//		lattice_size_t y2 = ((lattice_size_t)floor((p1.y+pr)/lattice->getSpacing()))+1;
//
//		// Fill in the circle.
//		for (lattice_size_t x = x1; x<=x2; x++)
//			for (lattice_size_t y = y1; y<=y2; y++)
//			{
//				point c((((double)x)+0.5)*lattice->getSpacing()-p1.x, (((double)y)+0.5)*lattice->getSpacing()-p1.y, 0.0);
//				if (c.x*c.x+c.y*c.y < pr*pr) lattice->setSiteType(x, y, z, type);
//			}
//	}
//
//	// Go through all of the intermediate z planes and fill them with rings.
//	for (lattice_size_t z=z1+1; z<z2-1; z++)
//	{
//		// See which area of the capsule we are in.
//		si_dist_t subvolumeZ = (((double)z)+0.5)*lattice->getSpacing();
//		si_dist_t pr1=0.0, pr2=0.0;
//		if (subvolumeZ <= p1.z)
//		{
//			si_dist_t b = p1.z-subvolumeZ;
//			pr1 = sqrt(outerRadius*outerRadius-b*b);
//			if (b<innerRadius)
//				pr2 = sqrt(innerRadius*innerRadius-b*b);
//			else
//				pr2 = 0.0;
//		}
//		else if (subvolumeZ > p2.z)
//		{
//			si_dist_t b = subvolumeZ-p2.z;
//			pr1 = sqrt(outerRadius*outerRadius-b*b);
//			if (b<innerRadius)
//				pr2 = sqrt(innerRadius*innerRadius-b*b);
//			else
//				pr2 = 0.0;
//		}
//		else
//		{
//			pr1 = outerRadius;
//			pr2 = innerRadius;
//		}
//
//		// Get the bounds of the circle in this plane.
//		lattice_size_t x1 = (lattice_size_t)floor((p1.x-pr1)/lattice->getSpacing());
//		lattice_size_t x2 = ((lattice_size_t)floor((p1.x+pr1)/lattice->getSpacing()))+1;
//		lattice_size_t y1 = (lattice_size_t)floor((p1.y-pr1)/lattice->getSpacing());
//		lattice_size_t y2 = ((lattice_size_t)floor((p1.y+pr1)/lattice->getSpacing()))+1;
//
//		// Fill in the ring.
//		for (lattice_size_t x = x1; x<=x2; x++)
//			for (lattice_size_t y = y1; y<=y2; y++)
//			{
//				point c((((double)x)+0.5)*lattice->getSpacing()-p1.x, (((double)y)+0.5)*lattice->getSpacing()-p1.y, 0.0);
//				if (c.x*c.x+c.y*c.y <= pr1*pr1 && c.x*c.x+c.y*c.y >= pr2*pr2) lattice->setSiteType(x, y, z, type);
//			}
//
//		// Make sure the ring is fully connected in x and y.
//	}
//
//	// The final z plane is also a solid circle.
//	{
//		lattice_size_t z = z2;
//
//		// Get the distance to p1.
//		si_dist_t b = ((((double)z)+0.5)*lattice->getSpacing())-p2.z;
//
//		// Get the radius of the circle in this plane.
//		si_dist_t pr = sqrt(outerRadius*outerRadius-b*b);
//
//		// Get the bounds of the circle in this plane.
//		lattice_size_t x1 = (lattice_size_t)floor((p2.x-pr)/lattice->getSpacing());
//		lattice_size_t x2 = ((lattice_size_t)floor((p2.x+pr)/lattice->getSpacing()))+1;
//		lattice_size_t y1 = (lattice_size_t)floor((p2.y-pr)/lattice->getSpacing());
//		lattice_size_t y2 = ((lattice_size_t)floor((p2.y+pr)/lattice->getSpacing()))+1;
//
//		// Fill in the circle.
//		for (lattice_size_t x = x1; x<=x2; x++)
//			for (lattice_size_t y = y1; y<=y2; y++)
//			{
//				point c((((double)x)+0.5)*lattice->getSpacing()-p2.x, (((double)y)+0.5)*lattice->getSpacing()-p2.y, 0.0);
//				if (c.x*c.x+c.y*c.y < pr*pr) lattice->setSiteType(x, y, z, type);
//			}
//	}
//}

}
}
