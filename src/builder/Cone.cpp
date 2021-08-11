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
#include "builder/Cone.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {

    Cone::Cone(point _center, si_dist_t _radius, si_dist_t _height, site_t _type, vector _normal)
    :Shape(CONE
          ,bounding_box(point(_center.x-2.0*::max(_radius,_height), _center.y-2.0*::max(_radius,_height), _center.z-2.0*::max(_radius,_height))
                       ,point(_center.x+2.0*::max(_radius,_height), _center.y+2.0*::max(_radius,_height), _center.z+2.0*::max(_radius,_height)))
          ,_type)
    ,center(_center)
    ,radius(_radius)
    ,height(_height)
    ,normal(_normal.unitize()) 
    {
    }

    Cone::~Cone() {
    }
    
    bool Cone::intersects(Shape * query) {
        // TODO: implement
        return false;
    }

    bool Cone::contains(point query) {
       // 1) Compute projection onto normal to get length
       double proj = vector(center, query).dot(normal);
       // 2) Check if "below" base or "above" apex
       if(proj < 0.0 || proj > height)
           return false;
       // 3) Find scaled center (e.g. how far up the center spine)
       vector h_prime = normal.scale(proj);
       vector c_prime = vector(center.plus(h_prime.toPoint()));
       // 4) Find half angle of cone
       double angle = atan(radius/height);
       // 5) Find scaled radius
       double r_prime = (height-h_prime.length())*tan(angle);
       // 6) compute distance to scaled center
       if(c_prime.toPoint().distanceSquared(query) <= r_prime*r_prime)
           return true;

       return false;
    }

    bool Cone::contains(Shape * query) {
        // TODO: implement
        return false;
    }

    double Cone::getVolume() {
       return 1.0/3.0*(PI*radius*radius)*height;
    }

}
}

