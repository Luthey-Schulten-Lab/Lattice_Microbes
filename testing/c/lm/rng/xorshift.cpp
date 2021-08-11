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

#include <iostream>
#include <fstream>
#include <vector>
#include "TestingHelper.h"
#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/rng/XORShift.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

using std::ofstream;
using std::vector;
using lm::rng::XORShift;

class XORShiftTester : public XORShift
{
public:
    XORShiftTester(uint32_t seedTopA, uint32_t seedBottomA)
    :XORShift(seedTopA,seedBottomA),stateP(&state)
    {}
    unsigned long long * stateP;
};

BOOST_AUTO_TEST_SUITE(XORShiftTest)

BOOST_AUTO_TEST_CASE(Constructor)
{
    XORShiftTester * t = NULL;

    // Test a simple constructor.
    BOOST_REQUIRE_NO_THROW(t=new XORShiftTester(1, 2011));
    BOOST_CHECK_EQUAL(*(t->stateP), 1);
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}

/*BOOST_AUTO_TEST_CASE(GetRandomDouble)
{
    XORWowTester * t = NULL;
    double d;
    void * p1;
    void * p2;

    // Test a simple series of calls to getRandomDouble.
    BOOST_REQUIRE_NO_THROW(t=new XORWowTester(0, 1, 2011));
    BOOST_CHECK_EQUAL(*(t->nextValueP), *(t->numberValuesP));

    p1=*(t->randomValuesP);
    p2=*(t->nextRandomValuesP);

    // Check the first set.
    d = t->getRandomDouble();
    BOOST_CHECK_GE(d, 0.0);
    BOOST_CHECK_LE(d, 1.0);
    BOOST_CHECK_EQUAL(*(t->nextValueP), 1);
    BOOST_CHECK_EQUAL(*(t->randomValuesP), p2);
    BOOST_CHECK_EQUAL(*(t->nextRandomValuesP), p1);
    for (uint i=1; i<*(t->numberValuesP); i++)
    {
        d = t->getRandomDouble();
        if (i < 10) printf("%d: %e\n",i,d);
        BOOST_CHECK_GE(d, 0.0);
        BOOST_CHECK_LE(d, 1.0);
        BOOST_CHECK_EQUAL(*(t->nextValueP), i+1);
    }

    // Check the start of the second set.
    BOOST_CHECK_EQUAL(*(t->randomValuesP), p2);
    BOOST_CHECK_EQUAL(*(t->nextRandomValuesP), p1);
    d = t->getRandomDouble();
    BOOST_CHECK_GE(d, 0.0);
    BOOST_CHECK_LE(d, 1.0);
    BOOST_CHECK_EQUAL(*(t->nextValueP), 1);
    BOOST_CHECK_EQUAL(*(t->randomValuesP), p1);
    BOOST_CHECK_EQUAL(*(t->nextRandomValuesP), p2);

    // Check the second set.
    for (uint i=1; i<*(t->numberValuesP); i++)
    {
        d = t->getRandomDouble();
        BOOST_CHECK_GE(d, 0.0);
        BOOST_CHECK_LE(d, 1.0);
        BOOST_CHECK_EQUAL(*(t->nextValueP), i+1);
    }
    BOOST_CHECK_EQUAL(*(t->randomValuesP), p1);
    BOOST_CHECK_EQUAL(*(t->nextRandomValuesP), p2);
    d = t->getRandomDouble();
    BOOST_CHECK_GE(d, 0.0);
    BOOST_CHECK_LE(d, 1.0);
    BOOST_CHECK_EQUAL(*(t->nextValueP), 1);
    BOOST_CHECK_EQUAL(*(t->randomValuesP), p2);
    BOOST_CHECK_EQUAL(*(t->nextRandomValuesP), p1);

    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}

BOOST_AUTO_TEST_CASE(GetExpRandomDouble)
{
    XORWowTester * t = NULL;
    double d;
    void * p1;
    void * p2;

    // Test a simple series of calls to getRandomDouble.
    BOOST_REQUIRE_NO_THROW(t=new XORWowTester(0, 1, 2011));
    BOOST_CHECK_EQUAL(*(t->nextValueP), *(t->numberValuesP));

    p1=*(t->expRandomValuesP);
    p2=*(t->nextExpRandomValuesP);

    // Check the first set.
    d = t->getExpRandomDouble();
    BOOST_CHECK_GT(d, 0.0);
    BOOST_CHECK_LE(d, 50.0);
    BOOST_CHECK_EQUAL(*(t->nextValueP), 1);
    BOOST_CHECK_EQUAL(*(t->expRandomValuesP), p2);
    BOOST_CHECK_EQUAL(*(t->nextExpRandomValuesP), p1);
    for (uint i=1; i<*(t->numberValuesP); i++)
    {
        d = t->getExpRandomDouble();
        if (i < 10) printf("%d: %e\n",i,d);
        BOOST_CHECK_GT(d, 0.0);
        BOOST_CHECK_LE(d, 50.0);
        BOOST_CHECK_EQUAL(*(t->nextValueP), i+1);
    }

    // Check the start of the second set.
    BOOST_CHECK_EQUAL(*(t->expRandomValuesP), p2);
    BOOST_CHECK_EQUAL(*(t->nextExpRandomValuesP), p1);
    d = t->getExpRandomDouble();
    BOOST_CHECK_GT(d, 0.0);
    BOOST_CHECK_LE(d, 50.0);
    BOOST_CHECK_EQUAL(*(t->nextValueP), 1);
    BOOST_CHECK_EQUAL(*(t->expRandomValuesP), p1);
    BOOST_CHECK_EQUAL(*(t->nextExpRandomValuesP), p2);

    // Check the second set.
    for (uint i=1; i<*(t->numberValuesP); i++)
    {
        d = t->getExpRandomDouble();
        BOOST_CHECK_GT(d, 0.0);
        BOOST_CHECK_LE(d, 50.0);
        BOOST_CHECK_EQUAL(*(t->nextValueP), i+1);
    }
    BOOST_CHECK_EQUAL(*(t->expRandomValuesP), p1);
    BOOST_CHECK_EQUAL(*(t->nextExpRandomValuesP), p2);
    d = t->getExpRandomDouble();
    BOOST_CHECK_GT(d, 0.0);
    BOOST_CHECK_LE(d, 50.0);
    BOOST_CHECK_EQUAL(*(t->nextValueP), 1);
    BOOST_CHECK_EQUAL(*(t->expRandomValuesP), p2);
    BOOST_CHECK_EQUAL(*(t->nextExpRandomValuesP), p1);

    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}*/

/**
 * To check in Matlab:
    data=hdrload('/tmp/lm_rng_xorshift_gaussian1.dat',0,0);
    [N,X]=hist(data,100);
    D=N./(sum(N)*(X(2)-X(1)));
    subplot(2,1,1);
    plot(X,D);
    hold('on');
    plot(X,normpdf(X,0,1),'r');
    hold('off');
    subplot(2,1,2);
    semilogy(X,D);
    hold('on');
    semilogy(X,normpdf(X,0,1),'r');
    hold('off');
 *
 */
BOOST_AUTO_TEST_CASE(GetNormRandomDouble)
{
    XORShiftTester * t = NULL;
    ofstream f;

    // Test a simple constructor.
    f.open((TMP_DIR+"/lm_rng_xorshift_gaussian1.dat").c_str(), std::ios::out|std::ios::binary);
    BOOST_REQUIRE_NO_THROW(t=new XORShiftTester(1, 2011));
    for (int i=0; i<1000000; i++)
    {
        f << t->getNormRandomDouble() << "\n";
    }
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;

    f.close();

}

/**
 * To check in Matlab:
    data=hdrload('/tmp/lm_rng_xorshift_gaussian2.dat',0,0);
    [N,X]=hist(data,100);
    D=N./(sum(N)*(X(2)-X(1)));
    subplot(2,1,1);
    plot(X,D);
    hold('on');
    plot(X,normpdf(X,0,1),'r');
    hold('off');
    subplot(2,1,2);
    semilogy(X,D);
    hold('on');
    semilogy(X,normpdf(X,0,1),'r');
    hold('off');
 *
 */
BOOST_AUTO_TEST_CASE(GetNormRandomDoubles)
{
    XORShiftTester * t = NULL;
    ofstream f;
    double * d;

    // Test a simple constructor.
    d = new double[1000000];
    f.open((TMP_DIR+"/lm_rng_xorshift_gaussian2.dat").c_str(), std::ios::out|std::ios::binary);
    BOOST_REQUIRE_NO_THROW(t=new XORShiftTester(1, 2011));
    t->getNormRandomDoubles(d, 1000000);
    for (int i=0; i<1000000; i++)
    {
        f << d[i] << "\n";
    }
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    f.close();
    delete [] d;
}

BOOST_AUTO_TEST_SUITE_END()
