/*
 * University of Illinois Open Source License
 * Copyright 2014-2018 Luthey-Schulten Group,
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
 * Author(s): Joseph R. Peterson
 */

#include "config.h"
#include "builder/Mesh.h"
#include "builder/Shape.h"
#include "builder/Intersection.h"
#include "rdme/Lattice.h"

namespace lm {
    namespace builder {
        
        Mesh::Mesh(Shape *s1, Shape *s2, site_t type)
        :Shape(MESH,bounding_box(point(center.x-radius, center.y-radius, center.z-radius),point(center.x+radius, center.y+radius, center.z+radius)), type),center(center),radius(radius)
        {
        }
        
        Sphere::~Sphere()
        {
        }
        
        bool Sphere::intersects(Shape * query)
        {
            // FIXME: NEEDS TO BE IMPLEMENTED //
            if (!boundingBoxesIntersect(query)) return false;
            
            if (query->getShapeType() == SPHERE)
            {
                //Sphere * querySphere = dynamic_cast<Sphere *>(query);
                Sphere * querySphere = (Sphere *)query;
                
                // If the distance between the two centers is less than the sum of their radii, these two spheres intersect
                return (center.distanceSquared(querySphere->center) <= (radius+querySphere->radius)*(radius+querySphere->radius));
            }
            else
            {
                return query->intersects(this);
            }
        }
        
        bool Sphere::contains(point query)
        {
            
            if ((center.x-query.x)*(center.x-query.x)+
                (center.y-query.y)*(center.y-query.y)+
                (center.z-query.z)*(center.z-query.z) < radius*radius )
                return true;
            
            return false;
        }
        
        bool Sphere::contains(Shape * query)
        {
            if (!boundingBoxesIntersect(query)) return false;
            
            
            if (query->getShapeType() == SPHERE)
            {
                Sphere * querySphere = (Sphere *)query;
                
                // Transform the sphere as if the capsule was oriented in the +z direction with p1 at 0,0,r and p2 at 0,0,l-r.
                point tc = querySphere->getCenter();
                si_dist_t tr=querySphere->getRadius();
                
                si_dist_t dx = center.x - tc.x;
                si_dist_t dy = center.y - tc.y;
                si_dist_t dz = center.z - tc.z;
                si_dist_t dr = radius - tr;
                
                if (dx*dx+dy*dy+dz*dz < dr*dr)
                    return true;
                
            }
            return false;
            
        }
        
        void Sphere::setCenter(point center)
        {
            this->center = center;
            boundingBox.min = point(center.x-radius,
                                    center.y-radius,
                                    center.z-radius);
            boundingBox.max = point(center.x+radius,
                                    center.y+radius,
                                    center.z+radius);
        }
        
        double Sphere::getVolume()
        {
            return (4.0/3.0)*PI*radius*radius*radius;
        }
        
    }
}
