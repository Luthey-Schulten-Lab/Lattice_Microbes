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

#ifndef LM_BUILDER_CONE_H_
#define LM_BUILDER_CONE_H_

#include "builder/Shape.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {

/// @class Cone
/// @brief A Shape that represents a Cone
class Cone : public Shape
{
public:
    /// @brief Create a Cone
    /// @param center Point center of the circle of the base
    /// @param radius Radius of the base
    /// @param height Height of the cone
    /// @param type The type of the sites within the cone
    /// @param normal Normal to the center of the cone base
    Cone(point center, si_dist_t radius, si_dist_t height, site_t type, vector normal = vector(1.0,0.0,0.0));
    /// @brief Destroy the Sphere
    virtual ~Cone();
    
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
    /// @brief Get the center of the cone
    virtual point getCenter() {return center;}
    /// @brief Get the radius of the cone
    virtual si_dist_t getRadius() {return radius;}
    /// @brief Get the height of the cone
    virtual si_dist_t getHeight() {return height;}
    /// @brief Get the volume bounded by the cone
    virtual double getVolume();

protected:
    point center;
    si_dist_t radius;
    si_dist_t height;
    vector normal;

};


}
}

#endif
