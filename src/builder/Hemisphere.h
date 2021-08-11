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

#ifndef LM_BUILDER_HEMISPHERE_H_
#define LM_BUILDER_HEMISPHERE_H_

#include "builder/Shape.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {

/// @class Hemisphere
/// @brief A hemisphere Shape
class Hemisphere : public Shape
{
public:
    /// @brief Create a Hemisphere
    /// @param center Point center of the circle of the slice plane through the Hemisphere
    /// @param radius Radius of the Hemisphere
    /// @param orientation Orientation normal to the center point of the center point in the direction of the curved part of the hemisphere
    /// @param type The type of the sites within the hemisphere
    Hemisphere(point center, si_dist_t radius, vector orientation, site_t type);
    /// @brief Destroy the Hemisphere
    virtual ~Hemisphere();
    
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
    /// @brief Get the center of the Hemisphere
    virtual point getCenter() {return center;}
    /// @brief Get the radius of the Hemisphere
    virtual si_dist_t getRadius() {return radius;}
    /// @brief Get the orientation vector of the Hemisphere
    virtual vector getOrientation() {return orientation;}
    /// @brief Get the volume contained by the Hemisphere
    virtual double getVolume();

protected:
    point center;
    si_dist_t radius;
    vector orientation;
};

}
}

#endif

