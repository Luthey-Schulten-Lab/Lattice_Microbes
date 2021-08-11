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

void cu_IndexXYWindowSites(unsigned int * inLattice, unsigned int * outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize);
void cu_IndexXYWindowAprons(unsigned int * inLattice, unsigned int * outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize);
void cu_CopyXYWindowSites(unsigned int * inLattice, unsigned int * outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize);
void cu_CopyXYWindowAprons(int mode, unsigned int * inLattice, unsigned int * outLattice, unsigned int latticeXSize, unsigned int latticeYSize, unsigned int latticeZSize);

BOOST_AUTO_TEST_SUITE(CopyXYWindow)

BOOST_AUTO_TEST_CASE(IndexXYWindowSites)
{
    unsigned int latticeXSize=128, latticeYSize=96, latticeZSize=64;
    unsigned int latticeXYZSize=latticeXSize*latticeYSize*latticeZSize;
    unsigned int * inLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];

    // Make sure we get the right indices.
    cu_IndexXYWindowSites(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_MESSAGE( outLattice[i] == ((x) | (y<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice[i] << " expecting " << ((x) | (y<<10) | (z<<20)) );
                BOOST_REQUIRE_MESSAGE( outLattice[i+latticeXYZSize] == ((x%16)+2 | (((y%8)+2)<<10) | ((((y%8)+2)*20+((x%16)+2))<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice[i+latticeXYZSize] << " expecting " << (z | ((x+2)<<10) | (y<<20)) );
            }

    delete[] outLattice;
    delete[] inLattice;
}

BOOST_AUTO_TEST_CASE(IndexXYWindowAprons)
{
    unsigned int latticeXSize=128, latticeYSize=96, latticeZSize=64;
    unsigned int latticeXYZSize=latticeXSize*latticeYSize*latticeZSize;
    unsigned int * inLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];

    // Make sure we get the right indices.
    cu_IndexXYWindowAprons(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;

                if (y < 8 && y%8 < 2)
                {
                    BOOST_REQUIRE_MESSAGE( outLattice[i] == ((x) | ((y-2)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice[i] << " expecting " << ((x) | ((y-2)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice[i+latticeXYZSize] == z*latticeXSize*latticeYSize + (y-2+latticeYSize)*latticeXSize + x, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice[i+latticeXYZSize] << " expecting " << z*latticeXSize*latticeYSize + (y-2+latticeYSize)*latticeXSize + x );
                }
                else if (y >= 8 && y%8 < 2)
                {
                    BOOST_REQUIRE_MESSAGE( outLattice[i] == ((x) | ((y-2)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice[i] << " expecting " << ((x) | ((y-2)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice[i+latticeXYZSize] == z*latticeXSize*latticeYSize + (y-2)*latticeXSize + x, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice[i+latticeXYZSize] << " expecting " << z*latticeXSize*latticeYSize + (y-2)*latticeXSize + x );
                }
                else if (y < latticeYSize-8 && y%8 >= 8-2)
                {
                    BOOST_REQUIRE_MESSAGE( outLattice[i] == ((x) | ((y+2)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice[i] << " expecting " << ((x) | ((y+2)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice[i+latticeXYZSize] == z*latticeXSize*latticeYSize + (y+2)*latticeXSize + x, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice[i+latticeXYZSize] << " expecting " << z*latticeXSize*latticeYSize + (y+2)*latticeXSize + x );
                }
                else if (y >= latticeYSize-8 && y%8 >= 8-2)
                {
                    BOOST_REQUIRE_MESSAGE( outLattice[i] == ((x) | ((y+2)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice[i] << " expecting " << ((x) | ((y+2)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice[i+latticeXYZSize] == z*latticeXSize*latticeYSize + (y+2-latticeYSize)*latticeXSize + x, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice[i+latticeXYZSize] << " expecting " << z*latticeXSize*latticeYSize + (y+2-latticeYSize)*latticeXSize + x );
                }
                else
                {
                    BOOST_REQUIRE_MESSAGE( outLattice[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice[i+latticeXYZSize] << " expecting " << 0 );
                }
            }

    delete[] outLattice;
    delete[] inLattice;
}

/*
    unsigned int * outLattice = outLattice1;
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0],outLattice[1],outLattice[2],outLattice[3],outLattice[4],outLattice[5],outLattice[6],outLattice[7],outLattice[8],outLattice[9],outLattice[10],outLattice[11],outLattice[12],outLattice[13],outLattice[14],outLattice[15],outLattice[16]);
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+latticeXSize],outLattice[1+latticeXSize],outLattice[2+latticeXSize],outLattice[3+latticeXSize],outLattice[4+latticeXSize],outLattice[5+latticeXSize],outLattice[6+latticeXSize],outLattice[7+latticeXSize],outLattice[8+latticeXSize],outLattice[9+latticeXSize],outLattice[10+latticeXSize],outLattice[11+latticeXSize],outLattice[12+latticeXSize],outLattice[13+latticeXSize],outLattice[14+latticeXSize],outLattice[15+latticeXSize],outLattice[16+latticeXSize]);
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+2*latticeXSize],outLattice[1+2*latticeXSize],outLattice[2+2*latticeXSize],outLattice[3+2*latticeXSize],outLattice[4+2*latticeXSize],outLattice[5+2*latticeXSize],outLattice[6+2*latticeXSize],outLattice[7+2*latticeXSize],outLattice[8+2*latticeXSize],outLattice[9+2*latticeXSize],outLattice[10+2*latticeXSize],outLattice[11+2*latticeXSize],outLattice[12+2*latticeXSize],outLattice[13+2*latticeXSize],outLattice[14+2*latticeXSize],outLattice[15+2*latticeXSize],outLattice[16+2*latticeXSize]);
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+3*latticeXSize],outLattice[1+3*latticeXSize],outLattice[2+3*latticeXSize],outLattice[3+3*latticeXSize],outLattice[4+3*latticeXSize],outLattice[5+3*latticeXSize],outLattice[6+3*latticeXSize],outLattice[7+3*latticeXSize],outLattice[8+3*latticeXSize],outLattice[9+3*latticeXSize],outLattice[10+3*latticeXSize],outLattice[11+3*latticeXSize],outLattice[12+3*latticeXSize],outLattice[13+3*latticeXSize],outLattice[14+3*latticeXSize],outLattice[15+3*latticeXSize],outLattice[16+3*latticeXSize]);
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+4*latticeXSize],outLattice[1+4*latticeXSize],outLattice[2+4*latticeXSize],outLattice[3+4*latticeXSize],outLattice[4+4*latticeXSize],outLattice[5+4*latticeXSize],outLattice[6+4*latticeXSize],outLattice[7+4*latticeXSize],outLattice[8+4*latticeXSize],outLattice[9+4*latticeXSize],outLattice[10+4*latticeXSize],outLattice[11+4*latticeXSize],outLattice[12+4*latticeXSize],outLattice[13+4*latticeXSize],outLattice[14+4*latticeXSize],outLattice[15+4*latticeXSize],outLattice[16+4*latticeXSize]);
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+5*latticeXSize],outLattice[1+5*latticeXSize],outLattice[2+5*latticeXSize],outLattice[3+5*latticeXSize],outLattice[4+5*latticeXSize],outLattice[5+5*latticeXSize],outLattice[6+5*latticeXSize],outLattice[7+5*latticeXSize],outLattice[8+5*latticeXSize],outLattice[9+5*latticeXSize],outLattice[10+5*latticeXSize],outLattice[11+5*latticeXSize],outLattice[12+5*latticeXSize],outLattice[13+5*latticeXSize],outLattice[14+5*latticeXSize],outLattice[15+5*latticeXSize],outLattice[16+5*latticeXSize]);
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+6*latticeXSize],outLattice[1+6*latticeXSize],outLattice[2+6*latticeXSize],outLattice[3+6*latticeXSize],outLattice[4+6*latticeXSize],outLattice[5+6*latticeXSize],outLattice[6+6*latticeXSize],outLattice[7+6*latticeXSize],outLattice[8+6*latticeXSize],outLattice[9+6*latticeXSize],outLattice[10+6*latticeXSize],outLattice[11+6*latticeXSize],outLattice[12+6*latticeXSize],outLattice[13+6*latticeXSize],outLattice[14+6*latticeXSize],outLattice[15+6*latticeXSize],outLattice[16+6*latticeXSize]);
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+7*latticeXSize],outLattice[1+7*latticeXSize],outLattice[2+7*latticeXSize],outLattice[3+7*latticeXSize],outLattice[4+7*latticeXSize],outLattice[5+7*latticeXSize],outLattice[6+7*latticeXSize],outLattice[7+7*latticeXSize],outLattice[8+7*latticeXSize],outLattice[9+7*latticeXSize],outLattice[10+7*latticeXSize],outLattice[11+7*latticeXSize],outLattice[12+7*latticeXSize],outLattice[13+7*latticeXSize],outLattice[14+7*latticeXSize],outLattice[15+7*latticeXSize],outLattice[16+7*latticeXSize]);
    printf("-----------------------------------------------------\n");
    printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u |  %u\n",outLattice[0+8*latticeXSize],outLattice[1+8*latticeXSize],outLattice[2+8*latticeXSize],outLattice[3+8*latticeXSize],outLattice[4+8*latticeXSize],outLattice[5+8*latticeXSize],outLattice[6+8*latticeXSize],outLattice[7+8*latticeXSize],outLattice[8+8*latticeXSize],outLattice[9+8*latticeXSize],outLattice[10+8*latticeXSize],outLattice[11+8*latticeXSize],outLattice[12+8*latticeXSize],outLattice[13+8*latticeXSize],outLattice[14+8*latticeXSize],outLattice[15+8*latticeXSize],outLattice[16+8*latticeXSize]);
*/

BOOST_AUTO_TEST_CASE(CopyXYWindowSites)
{
    unsigned int latticeXSize=128, latticeYSize=96, latticeZSize=64;
    unsigned int latticeXYZSize=latticeXSize*latticeYSize*latticeZSize;
    unsigned int * inLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];

    // Test a simple lattice with ones in the first word.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = 1;
                inLattice[i+latticeXYZSize] = 0;
            }
    cu_CopyXYWindowSites(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], inLattice[i+latticeXYZSize] );
            }

    // Test a simple lattice with ones in both words.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = 1;
                inLattice[i+latticeXYZSize] = 1;
            }
    cu_CopyXYWindowSites(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], inLattice[i+latticeXYZSize] );
            }

    // Test a lattice with the x,y,z in the first word and z,x,y in the second.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = x | (y<<8) | (z<<16);
                inLattice[i+latticeXYZSize] = z | (x<<8) | (y<<16);
            }
    cu_CopyXYWindowSites(inLattice, outLattice, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                BOOST_REQUIRE_EQUAL( outLattice[i], inLattice[i] );
                BOOST_REQUIRE_EQUAL( outLattice[i+latticeXYZSize], inLattice[i+latticeXYZSize] );
            }

    delete[] outLattice;
    delete[] inLattice;
}

BOOST_AUTO_TEST_CASE(CopyXYWindowAprons)
{
    unsigned int latticeXSize=384, latticeYSize=96, latticeZSize=64;
    unsigned int latticeXYZSize=latticeXSize*latticeYSize*latticeZSize;
    unsigned int apronSize=2;
    unsigned int * inLattice = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice1 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice2 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice3 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice4 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice5 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice6 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice7 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];
    unsigned int * outLattice8 = new unsigned int[latticeXSize*latticeYSize*latticeZSize*2];

    // Test a lattice with the 0 in the first word and 0 in the second to make sure we get the lattice boundaries correct.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = 0;
                inLattice[i+latticeXYZSize] = 0;
            }
    cu_CopyXYWindowAprons(1, inLattice, outLattice1, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(2, inLattice, outLattice2, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(3, inLattice, outLattice3, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(4, inLattice, outLattice4, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(5, inLattice, outLattice5, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(6, inLattice, outLattice6, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(7, inLattice, outLattice7, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(8, inLattice, outLattice8, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int bx = x/16;
                unsigned int by = y/8;
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;

                //-x,-y
                if ((bx == 0 || by == 0) && x%16 < apronSize && y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice1[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice1[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice1[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice1[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice1[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice1[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice1[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice1[i+latticeXYZSize] << " expecting " << 0 );
                }

                //-y
                if ((by == 0) && y < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice2[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice2[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice2[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice2[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice2[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice2[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice2[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice2[i+latticeXYZSize] << " expecting " << 0 );
                }

                //-x,-y
                if ((bx == ((latticeXSize/16)-1) || by == 0) && x%16 >= (16-apronSize) && y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice3[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice3[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice3[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice3[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice3[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice3[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice3[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice3[i+latticeXYZSize] << " expecting " << 0 );
                }

                //-x
                if ((bx == 0) && x < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice4[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice4[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice4[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice4[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice4[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice4[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice4[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice4[i+latticeXYZSize] << " expecting " << 0 );
                }

                //+x
                if ((bx == ((latticeXSize/16)-1)) && x%16 >= (16-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice5[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice5[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice5[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice5[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice5[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice5[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice5[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice5[i+latticeXYZSize] << " expecting " << 0 );
                }

                //-x,+y
                if ((bx == 0 || by == ((latticeYSize/8)-1)) && x%16 < apronSize && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice6[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice6[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice6[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice6[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice6[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice6[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice6[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice6[i+latticeXYZSize] << " expecting " << 0 );
                }

                //+y
                if ((by == ((latticeYSize/8)-1)) && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice7[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice7[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice7[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice7[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice7[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice7[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice7[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice7[i+latticeXYZSize] << " expecting " << 0 );
                }

                //+x,+y
                if ((bx == ((latticeXSize/16)-1) || by == ((latticeYSize/8)-1)) && x%16 >= (16-apronSize) && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice8[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice8[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice8[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice8[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice8[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice8[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice8[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice8[i+latticeXYZSize] << " expecting " << 0 );
                }
            }

    // Test a lattice with the x,y,z in the first word and z,x,y in the second.
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;
                inLattice[i] = x | (y<<10) | (z<<20);
                inLattice[i+latticeXYZSize] = z | (x<<10) | (y<<20);
            }
    cu_CopyXYWindowAprons(1, inLattice, outLattice1, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(2, inLattice, outLattice2, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(3, inLattice, outLattice3, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(4, inLattice, outLattice4, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(5, inLattice, outLattice5, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(6, inLattice, outLattice6, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(7, inLattice, outLattice7, latticeXSize, latticeYSize, latticeZSize);
    cu_CopyXYWindowAprons(8, inLattice, outLattice8, latticeXSize, latticeYSize, latticeZSize);
    for (unsigned int z=0; z<latticeZSize; z++)
        for (unsigned int y=0; y<latticeYSize; y++)
            for (unsigned int x=0; x<latticeXSize; x++)
            {
                unsigned int bx = x/16;
                unsigned int by = y/8;
                unsigned int i = z*latticeXSize*latticeYSize + y*latticeXSize + x;

                //-x,-y
                if ((bx == 0 || by == 0) && x%16 < apronSize && y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice1[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice1[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice1[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice1[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (x%16 < apronSize && y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice1[i] == ((x-apronSize) | ((y-apronSize)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice1[i] << " expecting " << ((x-apronSize) | ((y-apronSize)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice1[i+latticeXYZSize] == (z | ((x-apronSize)<<10) | ((y-apronSize)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice1[i+latticeXYZSize] << " expecting " << (z | ((x-apronSize)<<10) | ((y-apronSize)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice1[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice1[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice1[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice1[i+latticeXYZSize] << " expecting " << 0 );
                }

                //0,-y
                if ((by == 0) && y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice2[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice2[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice2[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice2[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice2[i] == ((x) | ((y-apronSize)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice2[i] << " expecting " << ((x-apronSize) | ((y)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice2[i+latticeXYZSize] == (z | ((x)<<10) | ((y-apronSize)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice2[i+latticeXYZSize] << " expecting " << (z | ((x)<<10) | ((y-apronSize)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice2[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice2[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice2[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice2[i+latticeXYZSize] << " expecting " << 0 );
                }

                //+x,-y
                if ((bx == ((latticeXSize/16)-1) || by == 0) && x%16 >= (16-apronSize) && y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice3[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice3[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice3[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice3[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (x%16 >= (16-apronSize) && y%8 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice3[i] == ((x+apronSize) | ((y-apronSize)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice3[i] << " expecting " << ((x+apronSize) | ((y-apronSize)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice3[i+latticeXYZSize] == (z | ((x+apronSize)<<10) | ((y-apronSize)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice3[i+latticeXYZSize] << " expecting " << (z | ((x+apronSize)<<10) | ((y-apronSize)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice3[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice3[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice3[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice3[i+latticeXYZSize] << " expecting " << 0 );
                }

                //-x
                if ((bx == 0) && x%16 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice4[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice4[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice4[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice4[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (x%16 < apronSize)
                {
                    BOOST_CHECK_MESSAGE( outLattice4[i] == ((x-apronSize) | ((y)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice4[i] << " expecting " << ((x-apronSize) | ((y)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice4[i+latticeXYZSize] == (z | ((x-apronSize)<<10) | ((y)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice4[i+latticeXYZSize] << " expecting " << (z | ((x-apronSize)<<10) | ((y)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice4[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice4[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice4[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice4[i+latticeXYZSize] << " expecting " << 0 );
                }

                //+x
                if ((bx == ((latticeXSize/16)-1)) && x%16 >= (16-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice5[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice5[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice5[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice5[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (x%16 >= (16-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice5[i] == ((x+apronSize) | ((y)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice5[i] << " expecting " << ((x+apronSize) | ((y)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice5[i+latticeXYZSize] == (z | ((x+apronSize)<<10) | ((y)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice5[i+latticeXYZSize] << " expecting " << (z | ((x+apronSize)<<10) | ((y)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice5[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice5[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice5[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice5[i+latticeXYZSize] << " expecting " << 0 );
                }

                //-x,+y
                if ((bx == 0 || by == ((latticeYSize/8)-1)) && x%16 < apronSize && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice6[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice6[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice6[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice6[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (x%16 < apronSize && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice6[i] == ((x-apronSize) | ((y+apronSize)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice6[i] << " expecting " << ((x-apronSize) | ((y+apronSize)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice6[i+latticeXYZSize] == (z | ((x-apronSize)<<10) | ((y+apronSize)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice6[i+latticeXYZSize] << " expecting " << (z | ((x-apronSize)<<10) | ((y+apronSize)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice6[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice6[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice6[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice6[i+latticeXYZSize] << " expecting " << 0 );
                }

                //+y
                if ((by == ((latticeYSize/8)-1)) && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice7[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice7[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice7[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice7[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice7[i] == ((x) | ((y+apronSize)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice7[i] << " expecting " << ((x) | ((y+apronSize)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice7[i+latticeXYZSize] == (z | ((x)<<10) | ((y+apronSize)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice7[i+latticeXYZSize] << " expecting " << (z | ((x)<<10) | ((y+apronSize)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice7[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice7[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice7[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice7[i+latticeXYZSize] << " expecting " << 0 );
                }

                //+x,+y
                if ((bx == ((latticeXSize/16)-1) || by == ((latticeYSize/8)-1)) && x%16 >= (16-apronSize) && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice8[i] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice8[i] << " expecting " << 0xFFEEDDCC );
                    BOOST_REQUIRE_MESSAGE( outLattice8[i+latticeXYZSize] == 0xFFEEDDCC, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice8[i+latticeXYZSize] << " expecting " << 0xFFEEDDCC );
                }
                else if (x%16 >= (16-apronSize) && y%8 >= (8-apronSize))
                {
                    BOOST_CHECK_MESSAGE( outLattice8[i] == ((x+apronSize) | ((y+apronSize)<<10) | (z<<20)), "Site " << x << "," << y << "," << z << " word 1 was " << outLattice8[i] << " expecting " << ((x+apronSize) | ((y+apronSize)<<10) | (z<<20)) );
                    BOOST_REQUIRE_MESSAGE( outLattice8[i+latticeXYZSize] == (z | ((x+apronSize)<<10) | ((y+apronSize)<<20)), "Site " << x << "," << y << "," << z << " word 2 was " << outLattice8[i+latticeXYZSize] << " expecting " << (z | ((x+apronSize)<<10) | ((y+apronSize)<<20)) );
                }
                else
                {
                    BOOST_CHECK_MESSAGE( outLattice8[i] == 0, "Site " << x << "," << y << "," << z << " word 1 was " << outLattice8[i] << " expecting " << 0 );
                    BOOST_REQUIRE_MESSAGE( outLattice8[i+latticeXYZSize] == 0, "Site " << x << "," << y << "," << z << " word 2 was " << outLattice8[i+latticeXYZSize] << " expecting " << 0 );
                }
            }

    delete[] outLattice8;
    delete[] outLattice7;
    delete[] outLattice6;
    delete[] outLattice5;
    delete[] outLattice4;
    delete[] outLattice3;
    delete[] outLattice2;
    delete[] outLattice1;
    delete[] inLattice;
}
BOOST_AUTO_TEST_SUITE_END()
