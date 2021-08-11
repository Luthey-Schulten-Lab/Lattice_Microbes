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

#ifndef LM_BUILDER_DIFFERENCE_H_
#define LM_BUILDER_DIFFERENCE_H_

namespace lm {
    namespace builder {
        
        /// @class Difference
        /// @brief A Shape that represents a Difference between the first and second object
        class Difference : public Shape
        {
        public:
            /// @brief Create a Difference
            /// @param s1 The first shape to Difference
            /// @param s2 The second shape to Difference
            /// @param type The type of the sites within the difference
            /// @param symmetric Determine if the difference is symmetric.  If false, difference is the 1st shape minus the second shape
            Difference(Shape *s1, Shape *s2, site_t type, bool symmetric = false);
            /// @brief Destroy the Sphere
            virtual ~Difference();
            
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
            /// @brief Get the volume bounded by the sphere
            virtual double getVolume(bool reintegrate = false);
            virtual double getVolume() {
                return this->getVolume(false);
            }
            
        protected:
            Shape *shape1;
            Shape *shape2;
            bool hasBeenIntegrated;
            bool symmetry;
            double storedVolume;
            
        };
        
    }
}

#endif /* defined(LM_BUILDER_DIFFERENCE) */
