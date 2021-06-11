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
 * Author(s): Andrew Magis, Elijah Roberts, Joseph R. Peterson
 */

#include "config.h"
#include "core/Print.h"
#include "builder/Shape.h"
#include "rdme/Lattice.h"


namespace lm {
namespace builder {


point point::minus(const vector &r) {
  return point(x-r.x, y-r.y, z-r.z);
}

point point::plus(const vector &r) {
  return point(x+r.x, y+r.y, z+r.z);
}


si_dist_t vector::length() {
   return sqrtf(x*x + y*y + z*z);
}

vector vector::unitize() {
   si_dist_t l = length();
   if(l > 0)
     return vector(x/l, y/l, z/l);
   else
     return vector(0.0, 0.0, 0.0);
}


si_dist_t vector::dot(const vector &r) {
   return x*r.x + y*r.y + z*r.z;
}

vector vector::cross(const vector &r) {
   return vector(y*r.z - z*r.y,
                 z*r.x - x*r.z,
                 x*r.y - y*r.x);
}

vector vector::scale(const si_dist_t r) {
  return vector(r*x, r*y, r*z);
}

matrix vector::mult(const vector &r) {
  return matrix(x*r.x, x*r.y, x*r.z,
                y*r.x, y*r.y, y*r.z,
                z*r.x, z*r.y, z*r.z);
}

vector vector::mult(const matrix &r) {
   return vector(x*r.m11 + y*r.m21 + z*r.m31,
                 y*r.m12 + y*r.m22 + z*r.m32,
                 z*r.m13 + y*r.m23 + z*r.m33);
}

matrix matrix::transpose() {
  return matrix(m11, m21, m31,
                m12, m22, m32,
                m13, m23, m33);
}

si_dist_t matrix::determinant() {
   return m11*m22*m33 + m12*m23*m31 + m13*m21*m32 - (m13*m22*m31 + m12*m21*m33 + m11*m23*m32);
}

si_dist_t matrix::trace() {
   return m11*m11 + m22*m22 + m33*m33;
}


vector matrix::mult(const vector &r) {
  return vector(m11*r.x + m12*r.y + m13*r.z,
                m21*r.x + m22*r.y + m23*r.z,
                m31*r.x + m32*r.y + m33*r.z);
}

matrix matrix::mult(const matrix &r) {
  return matrix(m11*r.m11+m12*r.m21+m13*r.m31, m11*r.m12+m12*r.m22+m13*r.m32, m11*r.m13+m12*r.m23+m13*r.m33,
                m21*r.m11+m22*r.m21+m23*r.m31, m21*r.m12+m22*r.m22+m23*r.m32, m21*r.m13+m22*r.m23+m23*r.m33,
                m31*r.m11+m32*r.m21+m33*r.m31, m31*r.m12+m32*r.m22+m33*r.m32, m31*r.m13+m32*r.m23+m33*r.m33);
}




// Shape Class //
Shape::Shape(ShapeType shapeType, bounding_box boundingBox, site_t type, vector at, vector up)
:shapeType(shapeType),boundingBox(boundingBox),type(type),at(at),up(up)
{
  // Determine angles between at and up
  double phi   = 0.0; // Rotation about x
  double theta = 0.0; // Rotation about y
  double psi   = 0.0; // Rotation about z

  // Determine unitized vectors
  vector Xvec = at.unitize();     // X in the rotated frame
  vector Yvec = up.unitize();     // Y in the rotated frame
  vector Zvec = Xvec.cross(Yvec); // Z in the rotated frame

//  std::cout << "X ("<< Xvec.x << "," << Xvec.y << "," << Xvec.z << ")" << std::endl;
//  std::cout << "Y ("<< Yvec.x << "," << Yvec.y << "," << Yvec.z << ")" << std::endl;
//  std::cout << "Z ("<< Zvec.x << "," << Zvec.y << "," << Zvec.z << ")" << std::endl << std::endl;


  phi   = atan2(Zvec.x, Zvec.y);
  theta = acos(Zvec.z);
  psi   = atan2(Yvec.z, Xvec.z);

  Print::printf(Print::DEBUG, "phi:%.4g\ttheta:%.4g\tpsi:%.4g",phi,theta,psi);

  // Determine the rotation matrix
  rotation = matrix::eulerMatrixFromAngles(phi, theta, psi);
}

Shape::~Shape()
{
}

bool Shape::boundingBoxesIntersect(Shape * query)
{
    bool xOverlap = (query->boundingBox.max.x >= boundingBox.min.x && query->boundingBox.min.x <= boundingBox.max.x);
    bool yOverlap = (query->boundingBox.max.y >= boundingBox.min.y && query->boundingBox.min.y <= boundingBox.max.y);
    bool zOverlap = (query->boundingBox.max.z >= boundingBox.min.z && query->boundingBox.min.z <= boundingBox.max.z);
    return xOverlap && yOverlap && zOverlap;
}

void Shape::discretizeTo(lm::rdme::Lattice * lattice)
{
	// Get the starting and ending subvolumes coordinates.
	lattice_size_t x1 = (lattice_size_t) max(0.0,floor(boundingBox.min.x/lattice->getSpacing()));
	lattice_size_t x2 = ((lattice_size_t)floor(boundingBox.max.x/lattice->getSpacing()))+1;
	lattice_size_t y1 = (lattice_size_t) max(0.0,floor(boundingBox.min.y/lattice->getSpacing()));
	lattice_size_t y2 = ((lattice_size_t)floor(boundingBox.max.y/lattice->getSpacing()))+1;
	lattice_size_t z1 = (lattice_size_t) max(0.0,floor(boundingBox.min.z/lattice->getSpacing()));
	lattice_size_t z2 = ((lattice_size_t)floor(boundingBox.max.z/lattice->getSpacing()))+1;
    
    // Constrain to the RDME domain
    x1 = ::max((lattice_size_t)0, x1);
    y1 = ::max((lattice_size_t)0, y1);
    z1 = ::max((lattice_size_t)0, z1);
    x2 = ::min(lattice->getXSize()-1, x2);
    y2 = ::min(lattice->getYSize()-1, y2);
    z2 = ::min(lattice->getZSize()-1, z2);

	// Go through each subvolume that contains a portion of this shape.
	for (lattice_size_t x = x1; x<=x2; x++)
		for (lattice_size_t y = y1; y<=y2; y++)
			for (lattice_size_t z = z1; z<=z2; z++)
			{
				// If the center of this subvolume is in the cuboid, mark it as part of the shape.
				point c((((double)x)+0.5)*lattice->getSpacing(), (((double)y)+0.5)*lattice->getSpacing(), (((double)z)+0.5)*lattice->getSpacing());
				if (contains(c)) lattice->setSiteType(x, y, z, type);
			}
}

    
// Functions used for Monte-Carlo integration of volume
double Shape::integrateToPercentThreshold(double percent) {
    double bbVol = this->boundingBox.volume();
    srand(time(NULL));
    
    double totalPoints = 0.0;
    double pointsInside = 0.0;
    double volumeThreshold = percent/100.0 * bbVol;
    double delta = 66e6;
    double lastChangedValue = 66e6;
    while(delta > volumeThreshold) {
        totalPoints += 1.0;
        if(this->contains(boundingBox.randomPointInside())) {
            pointsInside += 1.0;
            delta = fabs(lastChangedValue - (pointsInside / totalPoints) * bbVol);
            lastChangedValue = (pointsInside / totalPoints) * bbVol;
        }
    }
    
    
    return (pointsInside / totalPoints) * bbVol;
}
double Shape::integrateToCount(int count = 1000000) {
    double bbVol = this->boundingBox.volume();
    srand(time(NULL));
    
    // Count fraction of random points that end up inside
    double pointsInside = 0.0;
    for(int i = 0; i < count; i++) {
        if(this->contains(boundingBox.randomPointInside()))
            pointsInside += 1.0;
    }
    
    // Return fraction of volume that is inside points
    return pointsInside/count * bbVol;
}
    
    
    
}
}
