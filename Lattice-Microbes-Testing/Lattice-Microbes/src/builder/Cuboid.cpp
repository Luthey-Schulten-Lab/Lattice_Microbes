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

#include "config.h"
#include "builder/Shape.h"
#include "builder/Sphere.h"
#include "builder/Cuboid.h"
#include "rdme/Lattice.h"

namespace lm {
namespace builder {

Cuboid::Cuboid(point p1, point p2, site_t type)
:Shape(CUBOID,bounding_box(minBounds(p1,p2), maxBounds(p1,p2)), type, vector(1.0,0.0,0.0), vector(0.0,1.0,0.0)),p1(minBounds(p1,p2)),p2(maxBounds(p1,p2))
{
	width  = boundingBox.max.x - boundingBox.min.x;
	height = boundingBox.max.y - boundingBox.min.x;
	depth  = boundingBox.max.z - boundingBox.min.z;

	cen = point((boundingBox.max.x - boundingBox.min.x)/2.0 + boundingBox.min.x, 
                (boundingBox.max.y - boundingBox.min.y)/2.0 + boundingBox.min.y,
                (boundingBox.max.z - boundingBox.min.z)/2.0 + boundingBox.min.z);
}

Cuboid::Cuboid(point center, si_dist_t w, si_dist_t h, si_dist_t d, site_t type, vector at, vector up) 
:Shape(CUBOID,bounding_box(point(center.x-max(max(w,h),d), center.y-max(max(w,h),d), center.z-max(max(w,h),d)), point(center.x+max(max(w,h),d), center.y+max(max(w,h),d), center.z+max(max(w,h),d))), type, at, up)
{
	// Compute half distances
	si_dist_t w2 = w/2.0;
	si_dist_t h2 = h/2.0;
	si_dist_t d2 = d/2.0;
	width  = w;
	height = h;
	depth  = d;
	cen = center;


    // up (height)
	// ^   ^ cross (height)
    // |  / 
	// | /
    // |/
    // -------> at (width)
	//  P7 _______ P8
	//    /      /|
	// P3/______/ |P4
	//   |   C  | |
	//   | P5   | /P6
	//   |______|/
	//   P1     P2
	vector upU   = up.unitize();
	vector atU   = at.unitize();
	vector cross = atU.cross(upU).unitize();

	// Compute the lower and upper points
	point corners[8];

	// front face
	corners[0] = center.minus(atU.scale(w2)).minus(upU.scale(h2)).minus(cross.scale(d2));
	corners[1] = center.plus( atU.scale(w2)).minus(upU.scale(h2)).minus(cross.scale(d2));
	corners[2] = center.minus(atU.scale(w2)).plus( upU.scale(h2)).minus(cross.scale(d2));
	corners[3] = center.plus( atU.scale(w2)).plus( upU.scale(h2)).minus(cross.scale(d2));
	// back face
	corners[4] = center.minus(atU.scale(w2)).minus(upU.scale(h2)).plus(cross.scale(d2));
	corners[5] = center.plus( atU.scale(w2)).minus(upU.scale(h2)).plus(cross.scale(d2));
	corners[6] = center.minus(atU.scale(w2)).plus( upU.scale(h2)).plus(cross.scale(d2));
	corners[7] = center.plus( atU.scale(w2)).plus( upU.scale(h2)).plus(cross.scale(d2));
	
	// Compute bounding box limits
	double minX = 666e6;
	double minY = 666e6;
	double minZ = 666e6;
	double maxX = -666e6;
	double maxY = -666e6;
	double maxZ = -666e6;

	// TODO: Fix this so it is a tighter box
	for(int i = 0; i < 8; i++) {
		// Find minimum points
		// minX
		if(corners[i].x < minX)
			minX = corners[i].x;
		if(corners[i].y < minX)
			minX = corners[i].y;
		if(corners[i].z < minX)
			minX = corners[i].z;
		// minY
		if(corners[i].y < minY)
			minY = corners[i].y;
		if(corners[i].z < minY)
			minY = corners[i].z;
		if(corners[i].x < minY)
			minY = corners[i].x;
		// minZ
		if(corners[i].z < minZ)
			minZ = corners[i].z;
		if(corners[i].x < minZ)
			minZ = corners[i].x;
		if(corners[i].y < minZ)
			minZ = corners[i].y;

		// Find maximum points
		// maxX
		if(corners[i].x > maxX)
			maxX = corners[i].x;
		if(corners[i].y > maxX)
			maxX = corners[i].y;
		if(corners[i].z > maxX)
			maxX = corners[i].z;
		// maxY
		if(corners[i].y > maxY)
			maxY = corners[i].y;
		if(corners[i].z > maxY)
			maxY = corners[i].z;
		if(corners[i].x > maxY)
			maxY = corners[i].x;
		// maxZ
		if(corners[i].z > maxZ)
			maxZ = corners[i].z;
		if(corners[i].x > maxZ)
			maxZ = corners[i].x;
		if(corners[i].y > maxZ)
			maxZ = corners[i].y;
	}

	p1 = corners[0];
	p2 = corners[7];

	// Update the bounding box
	minX = max(0.0,minX);
	minY = max(0.0,minY);
	minZ = max(0.0,minZ);
	this->boundingBox = bounding_box(point(minX, minY, minZ), point(maxX, maxY, maxZ));
}

Cuboid::~Cuboid()
{
}

bool Cuboid::intersects(Shape * query)
{
    // FIXME: NEEDS TO BE IMPLEMENTED //
    if (!boundingBoxesIntersect(query)) return false;

    if (query->getShapeType() == SPHERE)
    {
        return false;
    }

    return false;
}

bool Cuboid::contains(point query)
{
    // P' = R.T.P
    point q2(query);

    // Do translation
    q2.x -= cen.x;
    q2.y -= cen.y;
    q2.z -= cen.z;

    // Do rotation
    vector query2 = rotation.mult(vector(q2)).toPoint();

    // TODO: Take care of Gimbal Lock

	return (query2.x >= -width/2.0  && query2.x <= width/2.0 &&
            query2.y >= -height/2.0 && query2.y <= height/2.0 &&
            query2.z >= -depth/2.0  && query2.z <= depth/2.0);
}

bool Cuboid::contains(Shape * query)
{
    if (!boundingBoxesIntersect(query)) return false;

    if (query->getShapeType() == SPHERE)
    {
        Sphere * querySphere = (Sphere *)query;
        point c = querySphere->getCenter();
        si_dist_t r=querySphere->getRadius();

        if (c.x >= p1.x+r && c.x <= p2.x-r &&
            c.y >= p1.y+r && c.y <= p2.y-r &&
            c.z >= p1.z+r && c.z <= p2.z-r) return true;
        return false;
    }

    return false;
}

double Cuboid::getVolume()
{
    return fabs(width*height*depth);
}

}
}
