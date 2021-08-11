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
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>


class LatticeTester : public lm::rdme::Lattice
{
public:
    LatticeTester(lattice_coord_t sizeA, si_dist_t spacingA):Lattice(sizeA,spacingA),sizeP(&size),numberSitesP(&numberSites),spacingP(&spacing){}
    LatticeTester(lattice_size_t xSizeA, lattice_size_t ySizeA, lattice_size_t zSizeA, si_dist_t spacingA):Lattice(xSizeA,ySizeA,zSizeA,spacingA),sizeP(&size),numberSitesP(&numberSites),spacingP(&spacing){}
    lattice_coord_t * sizeP;
    lattice_size_t * numberSitesP;
    si_dist_t * spacingP;

    virtual site_t getMaxSiteType() const {return 0;}
    virtual particle_t getMaxParticle() const {return 0;}
    virtual site_size_t getMaxOccupancy() const {return 0;}
    virtual site_t getSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z) const throw(lm::rdme::InvalidSiteException) {return 0;}
    virtual void setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t site) throw(lm::rdme::InvalidSiteException) {}
    virtual site_size_t getOccupancy(lattice_size_t x, lattice_size_t y, lattice_size_t z) const throw(lm::rdme::InvalidSiteException) {return 0;}
    virtual particle_t getParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_size_t particleIndex) const throw(lm::rdme::InvalidSiteException,lm::rdme::InvalidParticleException) {return 0;}
    virtual void addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle) throw(lm::rdme::InvalidSiteException,lm::rdme::InvalidParticleException) {}
    virtual void removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z) throw(lm::rdme::InvalidSiteException) {}

};

BOOST_AUTO_TEST_SUITE(LatticeTest)

BOOST_AUTO_TEST_CASE(CreateLattice)
{
    LatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new LatticeTester(32,64,128,10.0));
    BOOST_CHECK_EQUAL(t->sizeP->x, 32);
    BOOST_CHECK_EQUAL(t->sizeP->y, 64);
    BOOST_CHECK_EQUAL(t->sizeP->z, 128);
    BOOST_CHECK_EQUAL(*t->numberSitesP, 32*64*128);
    BOOST_CHECK_CLOSE(*t->spacingP, 10.0, 0.00001);
    if (t != NULL) delete t; t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new LatticeTester(lattice_coord_t(32,64,128),10.0));
    BOOST_CHECK_EQUAL(t->sizeP->x, 32);
    BOOST_CHECK_EQUAL(t->sizeP->y, 64);
    BOOST_CHECK_EQUAL(t->sizeP->z, 128);
    BOOST_CHECK_EQUAL(*t->numberSitesP, 32*64*128);
    BOOST_CHECK_CLOSE(*t->spacingP, 10.0, 0.00001);
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(getProperties)
{
    LatticeTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new LatticeTester(lattice_coord_t(32,64,128),10.0));
    BOOST_CHECK_EQUAL(t->getXSize(), 32);
    BOOST_CHECK_EQUAL(t->getYSize(), 64);
    BOOST_CHECK_EQUAL(t->getZSize(), 128);
    BOOST_CHECK_EQUAL(t->getSize().x, 32);
    BOOST_CHECK_EQUAL(t->getSize().y, 64);
    BOOST_CHECK_EQUAL(t->getSize().z, 128);
    BOOST_CHECK_EQUAL(t->getNumberSites(), 32*64*128);
    BOOST_CHECK_CLOSE(t->getSpacing(), 10.0, 0.00001);
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_SUITE_END()
