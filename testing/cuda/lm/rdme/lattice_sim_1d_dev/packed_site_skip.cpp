/*
 * University of Illinois Open Source License
 * Copyright 2010 Luthey-Schulten Group,
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

#include "TestingHelper.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

void cu_PackedSiteSkipX(unsigned int * inLattice, unsigned int * outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize);
void cu_PackedSiteSkipY(unsigned int * inLattice, unsigned int * outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize);
void cu_PackedSiteSkipZ(unsigned int * inLattice, unsigned int * outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize);

BOOST_AUTO_TEST_SUITE(PackedSiteSkip)

BOOST_AUTO_TEST_CASE(PackedSiteSkipX)
{
    unsigned int latticeXSize=128, latticeYSize=96, latticeZSize=64;
    unsigned int latticeXYZSize=latticeXSize*latticeYSize*latticeZSize;
    unsigned int * inLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];

    // Make sure a lattice with only three byte filled does not copy the second word.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = x | (y<<8) | (z<<16);
                inLattice[i+latticeXYZSize] = z | (x<<8) | (y<<16);
            }
    cu_PackedSiteSkipX(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], 0 );
            }

    // Make sure a lattice with four bytes filled does copy the second word.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                if (x < 64)
                    inLattice[i] = x | (y<<8) | (z<<16);
                else
                    inLattice[i] = x | (y<<8) | (z<<16) | (1<<24);
                inLattice[i+latticeXYZSize] = z | (x<<8) | (y<<16);
            }
    cu_PackedSiteSkipX(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                if (x < 64)
                    BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], 0 );
                else
                    BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], inLattice[i+latticeXYZSize] );
            }

    delete[] outLattice;
    delete[] inLattice;
}

BOOST_AUTO_TEST_CASE(PackedSiteSkipY)
{
    unsigned int latticeXSize=128, latticeYSize=96, latticeZSize=64;
    unsigned int latticeXYZSize=latticeXSize*latticeYSize*latticeZSize;
    unsigned int * inLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];

    // Make sure a lattice with only three byte filled does not copy the second word.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = x | (y<<8) | (z<<16);
                inLattice[i+latticeXYZSize] = z | (x<<8) | (y<<16);
            }
    cu_PackedSiteSkipY(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], 0 );
            }

    // Make sure a lattice with four bytes filled does copy the second word.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                if (x < 64)
                    inLattice[i] = x | (y<<8) | (z<<16);
                else
                    inLattice[i] = x | (y<<8) | (z<<16) | (1<<24);
                inLattice[i+latticeXYZSize] = z | (x<<8) | (y<<16);
            }
    cu_PackedSiteSkipY(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                if (x < 64)
                    BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], 0 );
                else
                    BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], inLattice[i+latticeXYZSize] );
            }

    delete[] outLattice;
    delete[] inLattice;
}

BOOST_AUTO_TEST_CASE(PackedSiteSkipZ)
{
    unsigned int latticeXSize=128, latticeYSize=96, latticeZSize=64;
    unsigned int latticeXYZSize=latticeXSize*latticeYSize*latticeZSize;
    unsigned int * inLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];

    // Make sure a lattice with only three byte filled does not copy the second word.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = x | (y<<8) | (z<<16);
                inLattice[i+latticeXYZSize] = z | (x<<8) | (y<<16);
            }
    cu_PackedSiteSkipZ(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], 0 );
            }

    // Make sure a lattice with four bytes filled does copy the second word.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                if (x < 64)
                    inLattice[i] = x | (y<<8) | (z<<16);
                else
                    inLattice[i] = x | (y<<8) | (z<<16) | (1<<24);
                inLattice[i+latticeXYZSize] = z | (x<<8) | (y<<16);
            }
    cu_PackedSiteSkipZ(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                if (x < 64)
                    BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], 0 );
                else
                    BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], inLattice[i+latticeXYZSize] );
            }

    delete[] outLattice;
    delete[] inLattice;
}

BOOST_AUTO_TEST_SUITE_END()
