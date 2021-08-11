/*
 * University of Illinois Open Source License
 * Copyright 2011 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
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

#include <string>
#include <map>
#include <vector>
#include "TestingHelper.h"
#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/rdme/Lattice.h"
#include "lm/rdme/ByteLattice.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

class ByteLatticeTester : public lm::rdme::ByteLattice
{
public:
    ByteLatticeTester(lattice_coord_t sizeA, si_dist_t spacingA, uint particlesPerSiteA):ByteLattice(sizeA,spacingA,particlesPerSiteA),sizeP(&size),numberSitesP(&numberSites),spacingP(&spacing),wordsPerSiteP(&wordsPerSite),particlesP(&particles),siteTypesP(&siteTypes){}
    ByteLatticeTester(lattice_size_t xSizeA, lattice_size_t ySizeA, lattice_size_t zSizeA, si_dist_t spacingA, uint particlesPerSiteA):ByteLattice(xSizeA,ySizeA,zSizeA,spacingA,particlesPerSiteA),sizeP(&size),numberSitesP(&numberSites),spacingP(&spacing),wordsPerSiteP(&wordsPerSite),particlesP(&particles),siteTypesP(&siteTypes){}
    lattice_coord_t * sizeP;
    lattice_size_t * numberSitesP;
    si_dist_t * spacingP;
    uint * wordsPerSiteP;
    uint32_t ** particlesP;
    uint8_t ** siteTypesP;
    virtual void deallocateMemory() throw(std::bad_alloc) {ByteLattice::deallocateMemory();}
};

BOOST_AUTO_TEST_SUITE(LatticeTest)

BOOST_AUTO_TEST_CASE(constructAndDestruct)
{
    ByteLatticeTester * t = NULL;

    try
    {
    	new ByteLatticeTester(32,64,128,10.0,8);
    }
    catch (lm::Exception & e)
    {
    	printf("%s\n", e.what());
    }
    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(32,64,128,10.0,8));
    BOOST_CHECK_EQUAL(t->sizeP->x, 32);
    BOOST_CHECK_EQUAL(t->sizeP->y, 64);
    BOOST_CHECK_EQUAL(t->sizeP->z, 128);
    BOOST_CHECK_EQUAL(*t->numberSitesP, 32*64*128);
    BOOST_CHECK_CLOSE(*t->spacingP, 10.0, 0.00001);
    BOOST_CHECK_EQUAL(*t->wordsPerSiteP, 2);
    BOOST_CHECK(*t->particlesP != NULL);
    BOOST_CHECK(*t->siteTypesP != NULL);
    if (t != NULL) delete t; t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(32,64,128),10.0,12));
    BOOST_CHECK_EQUAL(t->sizeP->x, 32);
    BOOST_CHECK_EQUAL(t->sizeP->y, 64);
    BOOST_CHECK_EQUAL(t->sizeP->z, 128);
    BOOST_CHECK_EQUAL(*t->numberSitesP, 32*64*128);
    BOOST_CHECK_CLOSE(*t->spacingP, 10.0, 0.00001);
    BOOST_CHECK_EQUAL(*t->wordsPerSiteP, 3);
    BOOST_CHECK(*t->particlesP != NULL);
    BOOST_CHECK(*t->siteTypesP != NULL);
    if (t != NULL) delete t; t = NULL;

    // Make sure the delete memory function works.
    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(32,64,128),10.0,8));
    BOOST_REQUIRE_NO_THROW(t->deallocateMemory());
    BOOST_CHECK(*t->particlesP == NULL);
    BOOST_CHECK(*t->siteTypesP == NULL);
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(getProperties)
{
    ByteLatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(32,64,128),10.0,8));
    BOOST_CHECK_EQUAL(t->getMaxSiteType(), 255);
    BOOST_CHECK_EQUAL(t->getMaxParticle(), 255);
    BOOST_CHECK_EQUAL(t->getMaxOccupancy(), 8);
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(getSetSiteType)
{
    uint xSize=32, ySize=64, zSize=128;
    ByteLatticeTester * t = NULL;

    // Test setting and getting sites.
    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(xSize,ySize,zSize),10.0,8));
    for (uint z=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++)
            {
                BOOST_CHECK_EQUAL(t->getSiteType(x,y,z), 0);
            }
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                BOOST_CHECK_NO_THROW(t->setSiteType(x,y,z,i%255));
            }
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                BOOST_CHECK_EQUAL(t->getSiteType(x,y,z), i%255);
            }
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(getParticles)
{
    uint xSize=32, ySize=64, zSize=128;
    ByteLatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(xSize,ySize,zSize),10.0,8));

    // Check for exception for an invalid site.
    BOOST_CHECK_THROW(t->getParticle(xSize,0,0,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(0,ySize,0,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(0,0,zSize,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(xSize,ySize,zSize,0), lm::rdme::InvalidSiteException);

    // Check for exception for an invalid particle.
    BOOST_CHECK_THROW(t->getParticle(0,0,0,8), lm::rdme::InvalidParticleException);

    // Test getting particles.
    for (uint z=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++)
            {
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,0), 0);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,1), 0);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,2), 0);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,3), 0);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,4), 0);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,5), 0);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,6), 0);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,7), 0);
            }
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                (*t->particlesP)[i] = i%255;
                (*t->particlesP)[i] |= (i%192)<<8;
                (*t->particlesP)[i] |= (i%128)<<16;
                (*t->particlesP)[i] |= (i%96)<<24;
                (*t->particlesP)[i+*(t->numberSitesP)] = i%77;
                (*t->particlesP)[i+*(t->numberSitesP)] |= (i%66)<<8;
                (*t->particlesP)[i+*(t->numberSitesP)] |= (i%55)<<16;
                (*t->particlesP)[i+*(t->numberSitesP)] |= (i%44)<<24;
            }
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,0), i%255);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,1), i%192);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,2), i%128);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,3), i%96);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,4), i%77);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,5), i%66);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,6), i%55);
                BOOST_CHECK_EQUAL(t->getParticle(x,y,z,7), i%44);
            }
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(addParticles)
{
    uint xSize=32, ySize=64, zSize=128;
    ByteLatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(xSize,ySize,zSize),10.0,8));

    // Check for exception for an invalid site.
    BOOST_CHECK_THROW(t->getParticle(xSize,0,0,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(0,ySize,0,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(0,0,zSize,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(xSize,ySize,zSize,0), lm::rdme::InvalidSiteException);

    // Test adding particles.
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
                for (uint j=0; j<i%9; j++)
                {
                    BOOST_CHECK_NO_THROW(t->addParticle(x,y,z,(i%100)+j));
                }
    for (uint z=0, i=0; z<1; z++)
        for (uint y=0; y<1; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                for (uint j=0; j<i%9; j++)
                {
                    BOOST_CHECK_EQUAL(t->getParticle(x,y,z,j), (i%100)+j);
                }
            }

    // Test adding a particle over the limit.
    BOOST_CHECK_THROW(t->addParticle(8,0,0,1), lm::rdme::InvalidParticleException);

    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(getOccupancy)
{
    uint xSize=32, ySize=64, zSize=128;
    ByteLatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(xSize,ySize,zSize),10.0,8));

    // Check for exception for an invalid site.
    BOOST_CHECK_THROW(t->getParticle(xSize,0,0,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(0,ySize,0,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(0,0,zSize,0), lm::rdme::InvalidSiteException);
    BOOST_CHECK_THROW(t->getParticle(xSize,ySize,zSize,0), lm::rdme::InvalidSiteException);

    // Test adding particles.
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
                for (uint j=0; j<i%9; j++)
                {
                    BOOST_CHECK_NO_THROW(t->addParticle(x,y,z,(i%100)+j));
                }
    for (uint z=0, i=0; z<1; z++)
        for (uint y=0; y<1; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                BOOST_CHECK_EQUAL(t->getOccupancy(x,y,z), i%9);
            }

    if (t != NULL) delete t; t = NULL;
}


BOOST_AUTO_TEST_CASE(removeParticles)
{
    uint xSize=32, ySize=64, zSize=128;
    ByteLatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(xSize,ySize,zSize),10.0,8));

    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
                for (uint j=0; j<i%9; j++)
                {
                    BOOST_CHECK_NO_THROW(t->addParticle(x,y,z,(i%100)+j));
                }
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                BOOST_CHECK_NO_THROW(t->removeParticles(x,y,z));
                BOOST_CHECK_EQUAL(t->getOccupancy(x,y,z), 0);
            }
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(removeAllParticles)
{
    uint xSize=32, ySize=64, zSize=128;
    ByteLatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(xSize,ySize,zSize),10.0,8));

    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
                for (uint j=0; j<i%9; j++)
                {
                    BOOST_CHECK_NO_THROW(t->addParticle(x,y,z,(i%100)+j));
                }
    t->removeAllParticles();
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
                BOOST_CHECK_EQUAL(t->getOccupancy(x,y,z), 0);
            }
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(getNeighboringSites)
{
    uint xSize=32, ySize=64, zSize=128;
    ByteLatticeTester * t = NULL;

    lattice_size_t neighbors[6];
    BOOST_REQUIRE_NO_THROW(t=new ByteLatticeTester(lattice_coord_t(xSize,ySize,zSize),10.0,8));
    for (uint z=0, i=0; z<zSize; z++)
        for (uint y=0; y<ySize; y++)
            for (uint x=0; x<xSize; x++, i++)
            {
            	t->getNeighboringSites(i, neighbors);
            	if (x > 0) BOOST_REQUIRE_EQUAL(neighbors[0], i-1);
            	else BOOST_REQUIRE_EQUAL(neighbors[0], i-1+xSize);
            	if (x < xSize-1) BOOST_REQUIRE_EQUAL(neighbors[1], i+1);
            	else BOOST_REQUIRE_EQUAL(neighbors[1], i+1-xSize);

            	if (y > 0) BOOST_REQUIRE_EQUAL(neighbors[2], i-xSize);
            	else BOOST_REQUIRE_EQUAL(neighbors[2], i-xSize+(xSize*ySize));
            	if (y < ySize-1) BOOST_REQUIRE_EQUAL(neighbors[3], i+xSize);
            	else BOOST_REQUIRE_EQUAL(neighbors[3], i+xSize-(xSize*ySize));

            	if (z > 0) BOOST_REQUIRE_EQUAL(neighbors[4], i-xSize*ySize);
            	else BOOST_REQUIRE_EQUAL(neighbors[4], i-xSize*ySize+(xSize*ySize*zSize));
            	if (z < zSize-1) BOOST_REQUIRE_EQUAL(neighbors[5], i+xSize*ySize);
            	else BOOST_REQUIRE_EQUAL(neighbors[5], i+xSize*ySize-(xSize*ySize*zSize));
            }

    if (t != NULL) delete t; t = NULL;
}



BOOST_AUTO_TEST_SUITE_END()
