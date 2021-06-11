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
#include <string>
#include "TestingHelper.h"
#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/cme/FluctuatingNRSolver.h"
#include "lm/main/ResourceAllocator.h"
#include "lm/rng/RandomGenerator.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

using std::map;
using std::string;
using lm::main::ResourceAllocator;
using lm::rng::RandomGenerator;

class FluctuatingNRSolverTester : public lm::cme::FluctuatingNRSolver
{
public:
    FluctuatingNRSolverTester():FluctuatingNRSolver(),replicateP(&replicate),parametersP(&parameters),resourcesP(&resources),rngP(&rng),\
                             numberSpeciesP(&numberSpecies),numberReactionsP(&numberReactions),initialSpeciesCountsP(&initialSpeciesCounts),\
                             speciesCountsP(&speciesCounts),SP(&S),DP(&D),propensityFunctionsP(&propensityFunctions),propensityFunctionArgsP(&propensityFunctionArgs),\
                             numberDependentSpeciesP(&numberDependentSpecies),dependentSpeciesP(&dependentSpecies),dependentSpeciesChangeP(&dependentSpeciesChange),numberDependentReactionsP(&numberDependentReactions),dependentReactionsP(&dependentReactions){}
    unsigned int * replicateP;
    map<string,string> ** parametersP;
    ResourceAllocator::ComputeResources ** resourcesP;
    RandomGenerator ** rngP;
    uint * numberSpeciesP;
    uint * numberReactionsP;
    uint ** initialSpeciesCountsP;
    uint ** speciesCountsP;
    int ** SP;
    uint ** DP;
    void *** propensityFunctionsP;
    void *** propensityFunctionArgsP;
    uint **numberDependentSpeciesP;
    uint *** dependentSpeciesP;
    int *** dependentSpeciesChangeP;
    uint **numberDependentReactionsP;
    uint *** dependentReactionsP;
    void destroyModelT() {destroyModel();}

    void * ZerothOrderPropensityFunction() {return (void *)&zerothOrderPropensity;}
    void * FirstOrderPropensityFunction() {return (void *)&firstOrderPropensity;}
    void * OUFirstOrderPropensityFunction() {return (void *)&ouFirstOrderPropensity;}
    double ZerothOrderPropensityArgs_k(void * ptr) {return ((ZerothOrderPropensityArgs *)ptr)->k;}
    uint FirstOrderPropensityArgs_si(void * ptr) {return ((FirstOrderPropensityArgs *)ptr)->si;}
    double FirstOrderPropensityArgs_k(void * ptr) {return ((FirstOrderPropensityArgs *)ptr)->k;}

    uint OUFirstOrderPropensityArgs_si(void * ptr) {return ((OUFirstOrderPropensityArgs *)ptr)->si;}
    double OUFirstOrderPropensityArgs_k(void * ptr) {return ((OUFirstOrderPropensityArgs *)ptr)->k;}
    uint OUFirstOrderPropensityArgs_oui(void * ptr) {return ((OUFirstOrderPropensityArgs *)ptr)->oui;}
    double OUFirstOrderPropensityArgs_noiseVariance(void * ptr) {return ((OUFirstOrderPropensityArgs *)ptr)->noiseVariance;}
    double OUFirstOrderPropensityArgs_noiseTau(void * ptr) {return ((OUFirstOrderPropensityArgs *)ptr)->noiseTau;}
    void * OUFirstOrderPropensityArgs_rng(void * ptr) {return (void *)((OUFirstOrderPropensityArgs *)ptr)->rng;}
};

BOOST_AUTO_TEST_SUITE(FluctuatingNRSolverTest)

BOOST_AUTO_TEST_CASE(Constructor)
{
    FluctuatingNRSolverTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new FluctuatingNRSolverTester());
    BOOST_CHECK_EQUAL(*(t->replicateP), (uint)-1);
    BOOST_CHECK_EQUAL(*(t->parametersP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->resourcesP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->rngP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), 0U);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), 0U);
    BOOST_CHECK_EQUAL(*(t->propensityFunctionsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->propensityFunctionArgsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->numberDependentSpeciesP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->dependentSpeciesP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->dependentSpeciesChangeP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->numberDependentReactionsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->dependentReactionsP), (void *)NULL);
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(BuildModel)
{
    FluctuatingNRSolverTester * t = NULL;

    // Make sure we can build a model with no noise terms.
    {
    uint ns         = 4;
    uint nr         = 6;
    uint c[]        = {  1,   0,   0,   0};
    uint o[]        = {  9,   9,   1,   1,   1,   1};
    double k[]      = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,
                         1,  -1,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,
                         0,   0,   0,   1,   0,  -1};
    uint D[]        = {  1,   0,   0,   0,   0,   0,
                         0,   1,   1,   0,   0,   0,
                         0,   0,   0,   1,   1,   0,
                         1,   1,   0,   0,   0,   1};
    double nvar[]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double ntau[]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double recalc   = 0.001;
    BOOST_REQUIRE_NO_THROW(t=new FluctuatingNRSolverTester());
    BOOST_REQUIRE_NO_THROW(t->buildModel(ns, nr, c, o, k, S, D, nvar, ntau, recalc));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], c[i]);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->SP))[i], S[i]);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->DP))[i], D[i]);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    for (uint i=2; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionsP))[i], (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionArgsP), (void *)NULL);
    for (uint i=2; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionArgsP))[i], (void *)NULL);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[2]), 1U);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[3]), 2U);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[4]), 2U);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[5]), 3U);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[2]), 3.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[3]), 4.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[4]), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[5]), 6.0, 1e-9);
    BOOST_CHECK_EQUAL((*(t->numberDependentSpeciesP))[0], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentSpeciesP))[1], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentSpeciesP))[2], 1U);
    BOOST_CHECK_EQUAL((*(t->numberDependentSpeciesP))[3], 1U);
    BOOST_CHECK_EQUAL((*(t->numberDependentSpeciesP))[4], 1U);
    BOOST_CHECK_EQUAL((*(t->numberDependentSpeciesP))[5], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[0][0], 0U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[0][1], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[1][0], 0U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[1][1], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[2][0], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[3][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[4][0], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesP))[5][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[0][0], -1);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[0][1], 1);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[1][0], 1);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[1][1], -1);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[2][0], 1);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[3][0], 1);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[4][0], -1);
    BOOST_CHECK_EQUAL((*(t->dependentSpeciesChangeP))[5][0], -1);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[0], 3U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[1], 3U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[2], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[3], 3U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[4], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[5], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][0], 0U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][1], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][2], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][0], 0U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][1], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][2], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[2][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[2][1], 4U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[3][0], 0U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[3][1], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[3][2], 5U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[4][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[4][1], 4U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[5][0], 0U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[5][1], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[5][2], 5U);
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    }

    // Build a model with a single noisy reaction.
    {
    uint ns         = 1;
    uint nr         = 1;
    uint nfr        = 1;
    uint c[]        = {  1};
    uint o[]        = {  1};
    double k[]      = {3.0};
    int S[]         = { -1};
    uint D[]        = {  1};
    double nvar[]   = {2.0};
    double ntau[]   = {10.0};
    double recalc   = 0.001;
    uint dependencyIndex[] = {0};
    uint noiseReactionIndex[] = {0};
    BOOST_REQUIRE_NO_THROW(t=new FluctuatingNRSolverTester());
    BOOST_REQUIRE_NO_THROW(t->buildModel(ns, nr, c, o, k, S, D, nvar, ntau, recalc));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns+nfr);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr+nfr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], c[i]);
    for (uint i=ns; i<ns+nfr; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], 0U);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=ns; i<ns+nfr; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=0; i<ns+nfr; i++)
        for (uint j=0; j<nr+nfr; j++)
            if (i<ns && j<nr)
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], S[(i*(nr))+j]);
            else if (i-ns == j-nr)
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], 1);
            else
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], 0);
    for (uint i=0; i<ns+nfr; i++)
        for (uint j=0; j<nr+nfr; j++)
            if (i<ns && j<nr)
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], D[(i*(nr))+j]);
            else if (i>=ns && j<nr && j==noiseReactionIndex[i-ns])
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], 1U);
            else
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], 0U);
    for (uint j=0, noiseIndex=0; j<nr+nfr; j++)
    {
        if (j<nr && nvar[j]>0)
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->OUFirstOrderPropensityFunction());
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[j]), dependencyIndex[j]);
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_oui((*(t->propensityFunctionArgsP))[j]), ns+noiseIndex);
            BOOST_CHECK_CLOSE(t->OUFirstOrderPropensityArgs_noiseVariance((*(t->propensityFunctionArgsP))[j]), nvar[noiseReactionIndex[noiseIndex]], 1e-9);
            BOOST_CHECK_CLOSE(t->OUFirstOrderPropensityArgs_noiseTau((*(t->propensityFunctionArgsP))[j]), ntau[noiseReactionIndex[noiseIndex]], 1e-9);
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_rng((*(t->propensityFunctionArgsP))[j]), *(t->rngP));
            noiseIndex++;
        }
        else if (j<nr)
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->FirstOrderPropensityFunction());
            BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[j]), dependencyIndex[j]);
            BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[j]), k[j], 1e-9);
        }
        else
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->ZerothOrderPropensityFunction());
            BOOST_CHECK_CLOSE(t->ZerothOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[j]), 1/(ntau[noiseReactionIndex[j-nr]]*recalc), 1e-9);
        }
    }
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    }

    // Build a model with two noisy reactions.
    {
    uint ns         = 1;
    uint nr         = 2;
    uint nfr        = 2;
    uint c[]        = {   1};
    uint o[]        = {   1,   1};
    double k[]      = { 3.0, 5.6};
    int S[]         = {  -1,   2};
    uint D[]        = {   1,   1};
    double nvar[]   = { 2.0, 0.1};
    double ntau[]   = {10.0, 1e5};
    double recalc   = 0.001;
    uint dependencyIndex[] = {0, 0};
    uint noiseReactionIndex[] = {0, 1};
    BOOST_REQUIRE_NO_THROW(t=new FluctuatingNRSolverTester());
    BOOST_REQUIRE_NO_THROW(t->buildModel(ns, nr, c, o, k, S, D, nvar, ntau, recalc));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns+nfr);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr+nfr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], c[i]);
    for (uint i=ns; i<ns+nfr; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], 0U);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=ns; i<ns+nfr; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=0; i<ns+nfr; i++)
        for (uint j=0; j<nr+nfr; j++)
            if (i<ns && j<nr)
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], S[(i*(nr))+j]);
            else if (i-ns == j-nr)
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], 1);
            else
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], 0);
    for (uint i=0; i<ns+nfr; i++)
        for (uint j=0; j<nr+nfr; j++)
            if (i<ns && j<nr)
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], D[(i*(nr))+j]);
            else if (i>=ns && j<nr && j==noiseReactionIndex[i-ns])
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], 1U);
            else
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], 0U);
    for (uint j=0, noiseIndex=0; j<nr+nfr; j++)
    {
        if (j<nr && nvar[j]>0)
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->OUFirstOrderPropensityFunction());
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[j]), dependencyIndex[j]);
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_oui((*(t->propensityFunctionArgsP))[j]), ns+noiseIndex);
            BOOST_CHECK_CLOSE(t->OUFirstOrderPropensityArgs_noiseVariance((*(t->propensityFunctionArgsP))[j]), nvar[noiseReactionIndex[noiseIndex]], 1e-9);
            BOOST_CHECK_CLOSE(t->OUFirstOrderPropensityArgs_noiseTau((*(t->propensityFunctionArgsP))[j]), ntau[noiseReactionIndex[noiseIndex]], 1e-9);
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_rng((*(t->propensityFunctionArgsP))[j]), *(t->rngP));
            noiseIndex++;
        }
        else if (j<nr)
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->FirstOrderPropensityFunction());
            BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[j]), dependencyIndex[j]);
            BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[j]), k[j], 1e-9);
        }
        else
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->ZerothOrderPropensityFunction());
            BOOST_CHECK_CLOSE(t->ZerothOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[j]), 1/(ntau[noiseReactionIndex[j-nr]]*recalc), 1e-9);
        }
    }
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    }

    // Build a model with noisy and non-noisy reactions.
    {
    uint ns         = 2;
    uint nr         = 5;
    uint nfr        = 2;
    uint c[]        = {   1,  10};
    uint o[]        = {   1,   1,   1,   1,   1};
    double k[]      = { 3.0, 5.6, 1e1, 1e2, 1e3};
    int S[]         = {  -1,   2,   0,   0,  -2,
                          0,  -1,   1,  -1,   1};
    uint D[]        = {   1,   0,   0,   0,   0,
                          0,   1,   1,   1,   1};
    double nvar[]   = { 0.0, 0.1, 0.0, 1.0, 0.0};
    double ntau[]   = { 0.0, 1e1, 0.0, 1e2, 0.0};
    double recalc   = 0.001;
    uint dependencyIndex[] = {0, 1, 1, 1, 1};
    uint noiseReactionIndex[] = {1, 3};
    BOOST_REQUIRE_NO_THROW(t=new FluctuatingNRSolverTester());
    BOOST_REQUIRE_NO_THROW(t->buildModel(ns, nr, c, o, k, S, D, nvar, ntau, recalc));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns+nfr);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr+nfr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], c[i]);
    for (uint i=ns; i<ns+nfr; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], 0U);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=ns; i<ns+nfr; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=0; i<ns+nfr; i++)
        for (uint j=0; j<nr+nfr; j++)
            if (i<ns && j<nr)
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], S[(i*(nr))+j]);
            else if (i-ns == j-nr)
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], 1);
            else
                BOOST_CHECK_EQUAL((*(t->SP))[(i*(nr+nfr))+j], 0);
    for (uint i=0; i<ns+nfr; i++)
        for (uint j=0; j<nr+nfr; j++)
            if (i<ns && j<nr)
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], D[(i*(nr))+j]);
            else if (i>=ns && j<nr && j==noiseReactionIndex[i-ns])
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], 1U);
            else
                BOOST_CHECK_EQUAL((*(t->DP))[(i*(nr+nfr))+j], 0U);
    for (uint j=0, noiseIndex=0; j<nr+nfr; j++)
    {
        if (j<nr && nvar[j]>0)
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->OUFirstOrderPropensityFunction());
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[j]), dependencyIndex[j]);
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_oui((*(t->propensityFunctionArgsP))[j]), ns+noiseIndex);
            BOOST_CHECK_CLOSE(t->OUFirstOrderPropensityArgs_noiseVariance((*(t->propensityFunctionArgsP))[j]), nvar[noiseReactionIndex[noiseIndex]], 1e-9);
            BOOST_CHECK_CLOSE(t->OUFirstOrderPropensityArgs_noiseTau((*(t->propensityFunctionArgsP))[j]), ntau[noiseReactionIndex[noiseIndex]], 1e-9);
            BOOST_CHECK_EQUAL(t->OUFirstOrderPropensityArgs_rng((*(t->propensityFunctionArgsP))[j]), *(t->rngP));
            noiseIndex++;
        }
        else if (j<nr)
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->FirstOrderPropensityFunction());
            BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[j]), dependencyIndex[j]);
            BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[j]), k[j], 1e-9);
        }
        else
        {
            BOOST_CHECK_EQUAL((*(t->propensityFunctionsP))[j], t->ZerothOrderPropensityFunction());
            BOOST_CHECK_CLOSE(t->ZerothOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[j]), 1/(ntau[noiseReactionIndex[j-nr]]*recalc), 1e-9);
        }
    }
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    }

/*
    for (uint i=0; i<ns+nfr; i++)
    {
        for (uint j=0; j<nr+nfr; j++)
            printf("%2d ",(*(t->SP))[(i*(nr+nfr))+j]);
        printf("\n");
    }
    printf("---\n");
    for (uint i=0; i<ns+nfr; i++)
    {
        for (uint j=0; j<nr+nfr; j++)
            printf("%2d ",(*(t->DP))[(i*(nr+nfr))+j]);
        printf("\n");
    }
*/
}

BOOST_AUTO_TEST_CASE(DestroyModel)
{
    FluctuatingNRSolverTester * t = NULL;

    uint ns         = 4;
    uint nr         = 6;
    uint c[]        = {  1,   0,   0,   0};
    uint o[]        = {  1,   1,   1,   1,   1,   1};
    double k[]      = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,
                         1,  -1,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,
                         0,   0,   0,   1,   0,  -1};
    uint D[]        = {  1,   0,   0,   0,   0,   0,
                         0,   1,   1,   0,   0,   0,
                         0,   0,   0,   1,   1,   0,
                         0,   0,   0,   0,   0,   1};
    BOOST_REQUIRE_NO_THROW(t=new FluctuatingNRSolverTester());
    BOOST_REQUIRE_NO_THROW(t->buildModel(ns, nr, c, o, k, S, D));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionArgsP), (void *)NULL);

    BOOST_REQUIRE_NO_THROW(t->destroyModelT());
    BOOST_CHECK_EQUAL(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->SP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->DP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->propensityFunctionsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->propensityFunctionArgsP), (void *)NULL);

    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}


BOOST_AUTO_TEST_SUITE_END()

