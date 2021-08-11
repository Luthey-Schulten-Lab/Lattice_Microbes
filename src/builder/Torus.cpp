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
 * Author(s): Joseph R. Peterson
 */

#include "config.h"
#include "builder/Shape.h"
#include "builder/Torus.h"
#include "rdme/Lattice.h"

namespace lm {
    namespace builder {
        
        Torus::Torus(point center, si_dist_t r1, si_dist_t r2, site_t type, vector orientation, vector up)
        :Shape(TORUS,bounding_box(point(center.x-(r1 + r2), center.y-(r1 + r2), center.z-r2),point(center.x+(r1 + r2), center.y+(r1 + r2), center.z+r2)), type, orientation, up),center(center),radius1(r1),radius2(r2)
        {
			// Compute a looser bounding box
			// TODO: Compute the tightest bounding box
			double minx = center.x - (r1+r2);
			double miny = center.y - (r1+r2);
			double minz = center.z - (r1+r2);
			minx = max(minx,0.0);
			miny = max(miny,0.0);
			minz = max(minz,0.0);

			double maxx = center.x + (r1+r2);
			double maxy = center.y + (r1+r2);
			double maxz = center.z + (r1+r2);

			this->boundingBox = bounding_box(point(minx,miny,minz), point(maxx,maxy,maxz));
        }
        
        Torus::~Torus()
        {
        }
        
        bool Torus::intersects(Shape * query)
        {
            if (!boundingBoxesIntersect(query)) return false;
            
            // FIXME: NEEDS TO BE IMPLEMENTED //
            if (query->getShapeType() == TORUS)
            {
                return true;
            }
            else
            {
                return query->intersects(this);
            }
        }
        
        bool Torus::contains(point query)
        {
			// Rotate the point
			//P'= R.T.P
			point q2(query);

			// Translate to origin
			q2.x -= center.x;
			q2.y -= center.y;
			q2.z -= center.z;

			// Do rotation
			vector query2 = rotation.mult(vector(q2)).toPoint();

            
            if (pow((query2.x*query2.x+
                     query2.y*query2.y+
                     query2.z*query2.z) -
                    (radius1*radius1+radius2*radius2),   2.0) -
                4.0*radius1*radius2*(radius2*radius2-(query2.z)*(query2.z)) < 0.0)
                return true;
            
            return false;
        }
        
        bool Torus::contains(Shape * query)
        {
            if (!boundingBoxesIntersect(query)) return false;
            
            // TODO: Correctly specify shape intersections
            if (query->getShapeType() == TORUS)
            {
                return true;
                
            }
            return false;
            
        }
        
        void Torus::setCenter(point center)
        {
            this->center = center;
            boundingBox.min = point(center.x-radius1,
                                    center.y-radius1,
                                    center.z-radius1);
            boundingBox.max = point(center.x+radius1,
                                    center.y+radius1,
                                    center.z+radius1);
        }
        
        double Torus::getVolume()
        {
            return 2.0*PI*PI*radius1*radius2*radius2;
        }
        
    }
}
