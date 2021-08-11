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

#ifndef LM_BUILDER_TORUS_H_
#define LM_BUILDER_TORUS_H_

#include "builder/Shape.h"
#include "rdme/Lattice.h"

namespace lm {
    namespace builder {
        
        /// @class Torus
        /// @brief A Shape that represents a Torus
        class Torus : public Shape
        {
        public:
            /// @brief Create a Torus
            /// @param center Point center of the circle of the slice plane through the sphere
            /// @param r1 Large radius of the Torus
            /// @param r2 Small radius of the Torus
            /// @param type The type of the sites within the sphere
			/// @param orientation The direction vector around which the ring is made, default: (0,0,1)
			/// @param up A vector to define the y axis of the object's coordinate system , default: (0,1,0)
            Torus(point center, si_dist_t r1, si_dist_t r2, site_t type, vector orientation = vector(0.0,0.0,1.0), vector up = vector(0.0,1.0,0.0));
            /// @brief Destroy the Torus
            virtual ~Torus();
            
            /// @brief Check if two shapes intersect
            /// @param query The other shape to check
            /// @return true/false
            virtual bool intersects(Shape * query);
            /// @brief Determine if the shape contains the specified point
            /// @param query Point to test
            /// @return true/false
            virtual bool contains(point query);
            /// @brief Determine if the shape contains the specified shape
            /// @param query Shape to test
            /// @return true/false
            virtual bool contains(Shape * query);
            /// @brief Set the center of the sphere
            /// @param center Point of the center
            virtual void setCenter(point center);
            /// @brief Get the center of the sphere
            virtual point getCenter() {return center;}
            /// @brief Get the large radius of the sphere
            virtual si_dist_t getRadius1() {return radius1;}
            /// @brief Get the small radius of the sphere
            virtual si_dist_t getRadius2() {return radius2;}
            /// @brief Get the volume bounded by the sphere
            virtual double getVolume();
            
        protected:
            point center;
            si_dist_t radius1; // big radius
            si_dist_t radius2; // small radius
            
        };
        
    }
}

#endif

