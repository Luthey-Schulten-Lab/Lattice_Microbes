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

#include "config.h"
#include "builder/Shape.h"
#include "builder/Hemisphere.h"
#include "builder/Sphere.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {

Hemisphere::Hemisphere(point center, si_dist_t radius, vector orientation, site_t type)
:Shape(HEMISPHERE,bounding_box(point(center.x-radius, center.y-radius, center.z-radius),point(center.x+radius, center.y+radius, center.z+radius)), type, orientation, vector::findPerpendicularVector(orientation)),center(center),radius(radius),orientation(orientation)
{
}

Hemisphere::~Hemisphere()
{
}

bool Hemisphere::intersects(Shape * query)
{
    // FIXME: NEEDS TO BE IMPLEMENTED //
    if (!boundingBoxesIntersect(query)) return false;

    //Sphere * querySphere = dynamic_cast<Sphere *>(query);
    Sphere * querySphere = (Sphere *)query;
    point c = querySphere->getCenter();
    si_dist_t r = querySphere->getRadius();

    // If the distance between the two centers is less than the sum of their radii, these two spheres intersect
    return (center.distanceSquared(c) <= (radius+r)*(radius+r));
}

bool Hemisphere::contains(point query)
{
    // FIXME: NEEDS TO BE IMPLEMENTED //
    if( center.distanceSquared(query) <= radius*radius ) {
        vector vec = vector(center, query);
        if( orientation.dot(vec) >= 0) {
            return true;
        }
    }

	return false;
}

bool Hemisphere::contains(Shape * query)
{
    if (!boundingBoxesIntersect(query)) return false;

    if (query->getShapeType() == SPHERE)
    {
        Sphere * querySphere = (Sphere *)query;

        // Transform the sphere as if the hemisphere was oriented in the +z direction and centered at 0,0,0.
        point c = querySphere->getCenter();
        si_dist_t x=c.x-center.x;
        si_dist_t y=c.y-center.y;
        si_dist_t z=c.z-center.z;
        si_dist_t r=querySphere->getRadius();

        // Make sure none of the sphere falls below the xy plane.
        if (z < r) return false;

        // Make sure the sphere is not closer than r to the hemisphere surface.
        if ((x*x+y*y+z*z) > (radius-r)*(radius-r)) return false;

        return true;
    }
    return false;
}

double Hemisphere::getVolume()
{
    return (2.0/3.0)*PI*radius*radius*radius;
}

}
}
