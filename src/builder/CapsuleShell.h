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
 * Author(s): Elijah Roberts
 */

#ifndef LM_BUILDER_CAPSULESHELL_H_
#define LM_BUILDER_CAPSULESHELL_H_

#include "builder/Shape.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {

/// @class CapsuleShell
/// @brief A capsule shell Shape with a hollow inside
class CapsuleShell : public Shape
{
public:
    /// @brief Create a Capsule
    /// @param p1 Point of first end of cylinder
    /// @param p2 Point of second end of cylinder
    /// @param innerRadius Inner radius of the capsule shell
    /// @param outerRadius Outer radius of the capsule shell
    /// @param type Type of the sites in the capsule shell
	CapsuleShell(point p1, point p2, si_dist_t innerRadius, si_dist_t outerRadius, site_t type);
    virtual ~CapsuleShell();
    
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
    
    /// @brief Get point of first end of cylinder
    virtual point getP1() {return p1;}
    /// @brief Get point of second end of cylinder
    virtual point getP2() {return p2;}
    /// @brief Get the inner radius of the cylinder
    virtual si_dist_t getInnerRadius() {return innerRadius;}
    /// @brief Get the outer radius of the cylinder
    virtual si_dist_t getOuterRadius() {return outerRadius;}
    /// @brief Get the total enclosed voluem of the CapsuleShell
    virtual double getVolume();

    /// @brief Discretize the object to the specified lattice
    /// @param lattice Lattice in which to discretize the shape
    //virtual void discretizeTo(lm::rdme::Lattice * lattice);

private:
    static bounding_box calcBoundingBox(point p1, point p2, si_dist_t outerRadius);

protected:
    point p1, p2;
    si_dist_t innerRadius, outerRadius;
    si_dist_t length;
};

}
}

#endif

