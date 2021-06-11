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
#include "builder/Shape.h"
#include "builder/Ellipse.h"
#include "rdme/Lattice.h"

namespace lm {
    namespace builder {
        
        Ellipse::Ellipse(point _center, si_dist_t _r1, si_dist_t _r2, si_dist_t _r3, site_t type, vector orientation1, vector orientation2)
        :Shape(ELLIPSE,
               bounding_box(point(_center.x - max(max(_r1,_r2),_r3),
                                  _center.y - max(max(_r1,_r2),_r3),
                                  _center.z - max(max(_r1,_r2),_r3)),
                            point(_center.x + max(max(_r1,_r2),_r3),
                                  _center.y + max(max(_r1,_r2),_r3),
                                  _center.z + max(max(_r1,_r2),_r3))),
               type,
               orientation1,
               orientation2)
        {
            r1 = _r1;
            r2 = _r2;
            r3 = _r3;
            center = _center;

			// TODO: Compute a tighter bounding box
        }
        
        Ellipse::~Ellipse()
        {
        }
        
        bool Ellipse::intersects(Shape * query)
        {
            if (!boundingBoxesIntersect(query)) return false;
            
            // TODO:
            return true;
        }
        
        bool Ellipse::contains(point query)
        {
			// P' = R.T.P
			point q2(query);
			
			// Do translation
			q2.x -= center.x;
			q2.y -= center.y;
			q2.z -= center.z;
			
			// Do rotation
			vector query2 = rotation.mult(vector(q2)).toPoint();

			// TODO: Take care of Gimbal Lock

            if ((query2.x)*(query2.x)/(r1*r1) +
                (query2.y)*(query2.y)/(r2*r2) +
                (query2.z)*(query2.z)/(r3*r3) <= 1.0) {
                return true;
            }
            
            return false;
        }
        
        bool Ellipse::contains(Shape * query)
        {
            if (!boundingBoxesIntersect(query)) return false;

            // TODO: Contains
            return true;
            
        }
        
        void Ellipse::setCenter(point center)
        {
            this->center = center;
            double maxRadius = max(r1,max(r2,r3));
            boundingBox.min = point(center.x-maxRadius,
                                    center.y-maxRadius,
                                    center.z-maxRadius);
            boundingBox.max = point(center.x+maxRadius,
                                    center.y+maxRadius,
                                    center.z+maxRadius);
        }
        
        double Ellipse::getVolume()
        {
            return 4.0/3.0*PI*r1*r2*r3;
        }
        
    }
}
