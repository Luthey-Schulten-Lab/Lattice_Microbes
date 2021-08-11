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
#include "builder/Difference.h"
#include "rdme/Lattice.h"

namespace lm {
    namespace builder {
        
        Difference::Difference(Shape *s1, Shape *s2, site_t type, bool symmetric)
        :Shape(DIFFERENCE,
               s1->getBoundingBox().joinWith(s2->getBoundingBox()),
               type)
        {
            shape1 = s1;
            shape2 = s2;
            hasBeenIntegrated = false;
            symmetry = false;
        }
        
        Difference::~Difference()
        {
        }
        
        bool Difference::intersects(Shape * query)
        {
            if (!boundingBoxesIntersect(query)) return false;
            
            bool inter1 = shape1->intersects(query);
            bool inter2 = shape2->intersects(query);
            // Cannot intersect both
            if(!symmetry) {
                if(inter1 && !inter2)
                    return true;
            } else {
                if((inter1 && !inter2) ||
                   (!inter1 && inter2))
                    return true;
            }
            
            return false;
        }
        
        bool Difference::contains(point query)
        {
            bool c1 = shape1->contains(query);
            bool c2 = shape2->contains(query);
            if(!symmetry) {
                // Cannot be contained in both
                if(c1 && !c2) {
                    return true;
                }
            } else {
                if((c1 && !c2) ||
                   (!c1 && c2))
                    return true;
            }
            
            return false;
        }
        
        bool Difference::contains(Shape * query)
        {
            if (!boundingBoxesIntersect(query)) return false;
            
            bool c1 = shape1->contains(query);
            bool c2 = shape2->contains(query);
            if(!symmetry) {
                // Cannot be contained in both
                if(c1 && !c2) {
                    return true;
                }
            } else {
                if((c1 && !c2) ||
                   (!c1 && c2))
                    return true;
            }
            
            return false;
            
        }
        
        double Difference::getVolume(bool reintegrate)
        {
            // This uses a Monte-Carlo integrator provided by the Shape base class
            if(!hasBeenIntegrated && !reintegrate) {
                hasBeenIntegrated = true;
                return (storedVolume = this->integrateToPercentThreshold(1.0e-6));
            } else
                return storedVolume;
        }
        
    }
}
