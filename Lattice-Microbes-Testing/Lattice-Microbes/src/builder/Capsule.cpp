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
#include "builder/Shape.h"
#include "builder/Sphere.h"
#include "builder/Capsule.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {


Capsule::Capsule(point _p1, point _p2, si_dist_t _radius, site_t _type)
:Shape(CAPSULE
      ,Capsule::calcBoundingBox(_p1,_p2,_radius)
      ,_type
      ,vector(_p1,_p2).unitize()
      ,vector::findPerpendicularVector(vector(_p1,_p2)).unitize())
,p1(_p1)
,p2(_p2)
,radius(_radius)
,length(_p1.distance(_p2)+_radius+_radius)
{
}

bounding_box Capsule::calcBoundingBox(point p1, point p2, si_dist_t radius)
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

Capsule::~Capsule()
{
}

bool Capsule::intersects(Shape * query)
{
    // FIXME: NEEDS TO BE IMPLEMENTED //
    if (!boundingBoxesIntersect(query)) return false;

    if (query->getShapeType() == SPHERE)
    {
    	return true;
    }
    return false;
}

bool Capsule::contains(point query)
{
    // Check if point is inside
    si_dist_t radius2 = radius*radius;
    if (query.distanceSquared(p1) <= radius2) // Lower cap region
    {
        return true;
    }
    else if (query.distanceSquared(p2) <= radius2) // Upper cap region
    {
        return true;
    }
    else // Cylinder region
    {
        // 1) Compute projection onto normal to get length
        double proj = vector(p1, query).dot(at);
        // 2) Check if "below" base or "above" apex
        if(proj < 0.0 || proj > p1.distance(p2))
            return false;
        // 3) Find scaled center (e.g. how far up the center spine)
        vector h_prime = at.scale(proj);
        vector c_prime = vector(p1.plus(h_prime.toPoint()));
        // 4) compute distance to scaled center
        if(c_prime.toPoint().distanceSquared(query) <= radius2)
            return true;
    }

    return false;
}

bool Capsule::contains(Shape * query)
{
	if (!boundingBoxesIntersect(query)) return false;


    if (query->getShapeType() == SPHERE)
    {
        Sphere * querySphere = (Sphere *)query;

        // Transform the sphere as if the capsule was oriented in the +z direction with p1 at 0,0,r and p2 at 0,0,l-r.
        point tc = querySphere->getCenter();
        si_dist_t tr=querySphere->getRadius();

        // Assume for now that the capsule is oriented with its long axis along the z axis.

        // Translate the sphere in the x-y plane.
        tc.x -= p1.x;
        tc.y -= p1.y;
        tc.z -= (p1.z-radius);

        // See if the sphere is closer to the lower cap, the cylinder, or the upper cap.
        if (tc.z < radius)
        {
            if (tc.x*tc.x+tc.y*tc.y+(tc.z-(radius))*(tc.z-(radius)) > (radius-tr)*(radius-tr)) return false;
            return true;
        }
        else if (tc.z > length-radius)
        {
            if (tc.x*tc.x+tc.y*tc.y+(tc.z-(length-radius))*(tc.z-(length-radius)) > (radius-tr)*(radius-tr)) return false;
            return true;
        }
        else
        {
            // Make sure the sphere does not fall outside the cylinder.
            if (tc.x*tc.x+tc.y*tc.y > (radius-tr)*(radius-tr)) return false;
            return true;
        }
    }
    return false;
}

double Capsule::getVolume()
{
    return ((4.0/3.0)*PI*radius*radius*radius)+(PI*radius*radius*(length-radius-radius));
}


}
}
