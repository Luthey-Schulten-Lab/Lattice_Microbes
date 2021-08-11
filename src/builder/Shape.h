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

#ifndef LM_BUILDER_SHAPE_H_
#define LM_BUILDER_SHAPE_H_

#include <cmath>
#include <cstdlib>
#include <ctime>
#include "core/Math.h"
#include "core/Types.h"
#include "rdme/Lattice.h"

namespace lm {

namespace rdme {
class Lattice;
}

namespace builder {

// Forward declarations
struct point;
struct vector;
struct matrix;



/// @struct point
/// @brief Type to store a position in space.
typedef struct point {
    /// @brief Create a point
    /// @param x X location
    /// @param y Y location
    /// @param z Z location
    point(si_dist_t x=0.0, si_dist_t y=0.0, si_dist_t z=0.0):x(x),y(y),z(z){}
    si_dist_t x;
    si_dist_t y;
    si_dist_t z;

	/// @brief Create a point from another point
	/// @param p2 The point to copy
	point(const point &p2):x(p2.x),y(p2.y),z(p2.z) {}

    /// @brief Determine the distance squared between two points
    /// @param p2 Second point
    /// @return scalar distance squared
    si_dist_t distanceSquared(const point & p2)
    {
        si_dist_t dx = p2.x - x;
        si_dist_t dy = p2.y - y;
        si_dist_t dz = p2.z - z;
        return (dx*dx + dy*dy + dz*dz);
    }

	/// @brief Subtract the vector from the point
	/// @param v The vector to subtract
	/// @return A new point at that location
	point minus(const vector &r);

	/// @brief Add the vector from the point
	/// @param v The vector to add
	/// @return A new point at that location
	point plus(const vector &r);

    /// @brief Determine the distance to another point
    /// @param p2 Second point
    /// @return scalar distance
    si_dist_t distance(const point & p2) {return sqrt(distanceSquared(p2));}
} point;

/// @struct bounding_box
/// @brief The bounds for a shape represented as a rectangular box
typedef struct bounding_box {
    /// @brief Create a bounding box for coordinates
    /// @param x1 Bottom x
    /// @param y1 Bottom y
    /// @param z1 Bottom z
    /// @param x2 Top x
    /// @param y2 Top y
    /// @param z2 Top z
    bounding_box(si_dist_t x1=0.0, si_dist_t y1=0.0, si_dist_t z1=0.0, si_dist_t x2=0.0, si_dist_t y2=0.0, si_dist_t z2=0.0):min(x1,y1,z1),max(x2,y2,z2){}
    /// @brief Create a bounding box from points
    /// @param min Minimum coordinate
    /// @param max Maximum coordinate
    bounding_box(point min, point max):min(min),max(max){}
    
    // Actual points of the bounding box
    point min, max;

    /// @brief Return a bounding box spanning the two bounding boxes
    /// @param j Other bounding box
    /// @return new bounding box
    bounding_box joinWith(bounding_box j)
    {
        return bounding_box(::min(j.min.x,min.x),::min(j.min.y,min.y),::min(j.min.z,min.z),::max(j.max.x,max.x),::max(j.max.y,max.y),::max(j.max.z,max.z));
    }
    
    /// @brief Returns the bounding box volume
    /// @return volume
    double volume() { return fabs((max.x - min.x)*(max.y - min.y)*(max.z - min.z)); }
    
    /// @brief Retun a random point inside the bounding box
    /// @return point
    point randomPointInside() {
        double rx = (float)rand() / (float)RAND_MAX;
        double ry = (float)rand() / (float)RAND_MAX;
        double rz = (float)rand() / (float)RAND_MAX;
        return point(rx * (max.x - min.x) + min.x,
                     ry * (max.y - min.y) + min.y,
                     rz * (max.z - min.z) + min.z);
    }
} bounding_box;

/// @struct vector
/// @brief A vector which points in a direction
typedef struct vector {
//////////////////
// Constructors //
//////////////////
    /// @brief Create a vector
    /// @param x Distance in x
    /// @param y Distance in y
    /// @param z Distance in z
    vector(si_dist_t x=0.0, si_dist_t y=0.0, si_dist_t z=0.0):x(x),y(y),z(z){}
   
    /// @brief Construct a vector pointing from one point to another
    /// @param from The point at the root of the vector
    /// @param to The point at the end of the vector
    vector(const point &from, const point &to) {
      x = to.x-from.x;
      y = to.y-from.y;
      z = to.z-from.z;
    }

	/// @brief Construct a vector from another vector
	/// @param v The vector to copy
	vector(const vector &v):x(v.x),y(v.y),z(v.z) {}

	/// @brief Construct a vector from a point
	/// @param p The point to cast
	vector(const point &p):x(p.x),y(p.y),z(p.z) {}

	/// @brief convert the vector to a point
	point toPoint() {
		return point(x,y,z);
	}

    static vector findPerpendicularVector(vector v) {
      if ( fabs(v.z) >= 1.0e-12 && fabs(-v.x - v.y) >= 1.0e-12)
        return vector(v.z,v.z, -v.x-v.y);
      else
        return vector(-v.y-v.z,v.x,v.x);
    }
////////////////////// 
// Unary Operations //
////////////////////// 
    /// @brief Get the vector length
    si_dist_t length();

    /// @brief Get a unit vector pointing in the same direction as this vector
    vector unitize();
/////////////////////// 
// Binary Operations //
/////////////////////// 
    /// @brief Compute the dot product between vectors
    /// @param r The right hand vector
    si_dist_t dot(const vector &r);

    /// @brief Compute the cross product between vectors
    /// @param r The right hand vector
    vector cross(const vector &r);

    /// @brief Scale the vector by a constant
    /// @param r The scalar
    vector scale(const si_dist_t r);

    /// @brief Multiply with another vector (taking this vector to be a row vector) to form a matrix
    /// @param r The column vector
    matrix mult(const vector &r);

    /// @brief Multiply with a matrix to form a new vector (assuming this vector is a row vector)
    /// @param r The matrix
    vector mult(const matrix &r);

/////////////// 
// Variables //
/////////////// 
    // coordinates
    si_dist_t x;
    si_dist_t y;
    si_dist_t z;
} vector;



/// @struct matrix
/// @brief A matrix used for rotations
typedef struct matrix {
//////////////////
// Constructors //
//////////////////
    /// @brief Create a matrix with the speficied elements
    matrix(si_dist_t m11 = 0.0, si_dist_t m12 = 0.0,si_dist_t m13 = 0.0,
           si_dist_t m21 = 0.0, si_dist_t m22 = 0.0,si_dist_t m23 = 0.0,
           si_dist_t m31 = 0.0, si_dist_t m32 = 0.0,si_dist_t m33 = 0.0)
     : m11(m11), m12(m12), m13(m13),
       m21(m21), m22(m22), m23(m23),
       m31(m31), m32(m32), m33(m33)
    {
    }

    /// @brief Get a forward rotation matrix from angles
    /// @param phi The angle around x (in radians)
    /// @param theta The angle around y (in radians)
    /// @param psi The angle around z (in radians)
    static matrix eulerMatrixFromAngles(si_dist_t phi, si_dist_t theta, si_dist_t psi) {
      matrix m( cos(theta)*cos(psi), cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi), sin(phi)*sin(psi) - cos(phi)*sin(theta)*cos(psi),
               -cos(theta)*sin(psi), cos(phi)*cos(psi) - sin(phi)*sin(theta)*sin(psi), sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi),
                sin(theta),         -sin(phi)*cos(theta),                              cos(phi)*cos(theta));

      return m;
    }

    /// @brief Get an identity matrix
    static matrix Identity() {
      return matrix(1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0);
    }

//////////////////////
// Unary Operations //
//////////////////////
    matrix transpose();

    si_dist_t determinant();

    si_dist_t trace();

///////////////////////
// Binary Operations //
///////////////////////
    /// @brief Multiply by a vector
    vector mult(const vector &r);

    /// @brief Multiply by a matrix
    matrix mult(const matrix &r);

///////////////
// Variables //
///////////////
    // Elements
    si_dist_t m11;
    si_dist_t m12;
    si_dist_t m13;
    si_dist_t m21;
    si_dist_t m22;
    si_dist_t m23;
    si_dist_t m31;
    si_dist_t m32;
    si_dist_t m33;
} matrix;


/// @class Shape
/// @brief Abstract base class for all the shapes in Lattice Microbes simulation builder
class Shape
{
public:
    /// @enum ShapeType
    /// @brief Possible shape types that can be used in Lattice Microbes
    enum ShapeType
    {
       SPHERE           =  1,
       HEMISPHERE       =  2,
       CYLINDER         =  3,
       CAPSULE          =  4,
       CUBOID           =  5,
       CAPSULE_SHELL    =  6,
       UNION            =  7,
       DIFFERENCE       =  8,
       INTERSECTION     =  9,
       MESH             = 10,
       TORUS            = 11,
       ELLIPSE          = 12,
       CONE             = 13,
       UNIONSET         = 14
    };

public:
    /// @brief Create a Shape
    /// @param shapeType The type of shape should be of type ShapeType
    /// @param boundingBox The extents of the object used for fast collision/contains checking
    /// @param type The type of site that the object represents
    /// @param at A vector describing what direction the object is pointing "at"; default: (0,0,1)
    /// @param up A vector describing what direction is "up" relative to "at"; default: (0,1,0)
    Shape(ShapeType shapeType, bounding_box boundingBox, site_t type, vector at = vector(0.0,0.0,1.0), vector up = vector(0.0,1.0,0.0));
    /// @brief Destroy the Shape
    virtual ~Shape();
    
    /// @brief Checks if another shape's bounding box interstects with this shape's bounding box
    /// @param query The other shape to test
    /// @return true/false
    virtual bool boundingBoxesIntersect(Shape * query);
    /// @brief Check if two shapes intersect
    /// @param query The other shape to check
    /// @return true/false
    virtual bool intersects(Shape * query) = 0;
    /// @brief Determine if the shape contains the specified point
    /// @param query Point to test
    /// @return true/false
    virtual bool contains(point query) = 0;
    /// @brief Determine if the shape contains the specified shape
    /// @param query Shape to test
    /// @return true/false
    virtual bool contains(Shape * query) = 0;
    
    /// @brief Get the bounding box
    virtual bounding_box getBoundingBox() {return boundingBox;}
    /// @brief Get the site type associated with the shape
    virtual site_t getType() {return type;}
    /// @brief Get the shape type
    virtual ShapeType getShapeType() {return shapeType;}
    /// @brief Get the total internal volume of the shape
    virtual double getVolume() = 0;

    /// @brief Discretize the object to the specified lattice
    /// @param lattice Lattice in which to discretize the shape
    virtual void discretizeTo(lm::rdme::Lattice * lattice);

protected:
    // Functions used for Monte-Carlo integration of volume
    double integrateToPercentThreshold(double percent); // Integrate until the fluctuation in value is less than
    double integrateToCount(int count); // Integrate a static number of steps
    
    virtual point minBounds(point p1, point p2) const {return point(min(p1.x,p2.x),min(p1.y,p2.y),min(p1.z,p2.z));}
    virtual point maxBounds(point p1, point p2) const {return point(max(p1.x,p2.x),max(p1.y,p2.y),max(p1.z,p2.z));}
    
    ShapeType shapeType;
    bounding_box boundingBox;
    site_t type;
   
    // An orientation for the shape
    vector at; // A vector pointing to where the object is heading
    vector up; // A vector representing "up"
	matrix rotation; // The rotation matrix based on Euler angles that converts a point to the same axes as (at, up)
 
    friend class LatticeBuilder;
};

}
}

#endif

