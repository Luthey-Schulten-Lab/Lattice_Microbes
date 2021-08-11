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
#include "lm/cme/GillespieDSolver.h"
#include "lm/io/ReactionModel.pb.h"
#include "lm/main/ResourceAllocator.h"
#include "lm/rng/RandomGenerator.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

using std::map;
using std::string;
using lm::io::ReactionModel;
using lm::main::ResourceAllocator;
using lm::rng::RandomGenerator;

class GillespieDSolverTester : public lm::cme::GillespieDSolver
{
public:
    GillespieDSolverTester():GillespieDSolver(),replicateP(&replicate),parametersP(&parameters),resourcesP(&resources),rngP(&rng),\
                             numberSpeciesP(&numberSpecies),numberReactionsP(&numberReactions),initialSpeciesCountsP(&initialSpeciesCounts),\
                             speciesCountsP(&speciesCounts),propensitiesP(&propensities),SP(&S),DP(&D),propensityFunctionsP(&propensityFunctions),propensityFunctionArgsP(&propensityFunctionArgs),\
                             numberDependentSpeciesP(&numberDependentSpecies),dependentSpeciesP(&dependentSpecies),dependentSpeciesChangeP(&dependentSpeciesChange),numberDependentReactionsP(&numberDependentReactions),dependentReactionsP(&dependentReactions){}
    unsigned int * replicateP;
    map<string,string> ** parametersP;
    ResourceAllocator::ComputeResources ** resourcesP;
    RandomGenerator ** rngP;
    uint * numberSpeciesP;
    uint * numberReactionsP;
    uint ** initialSpeciesCountsP;
    uint ** speciesCountsP;
    double ** propensitiesP;
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

    double ZerothOrderPropensityArgs_k(void * ptr) {return ((ZerothOrderPropensityArgs *)ptr)->k;}
    uint FirstOrderPropensityArgs_si(void * ptr) {return ((FirstOrderPropensityArgs *)ptr)->si;}
    double FirstOrderPropensityArgs_k(void * ptr) {return ((FirstOrderPropensityArgs *)ptr)->k;}
    uint SecondOrderPropensityArgs_s1i(void * ptr) {return ((SecondOrderPropensityArgs *)ptr)->s1i;}
    uint SecondOrderPropensityArgs_s2i(void * ptr) {return ((SecondOrderPropensityArgs *)ptr)->s2i;}
    double SecondOrderPropensityArgs_k(void * ptr) {return ((SecondOrderPropensityArgs *)ptr)->k;}
    uint SecondOrderSelfPropensityArgs_si(void * ptr) {return ((SecondOrderSelfPropensityArgs *)ptr)->si;}
    double SecondOrderSelfPropensityArgs_k(void * ptr) {return ((SecondOrderSelfPropensityArgs *)ptr)->k;}
};

BOOST_AUTO_TEST_SUITE(GillespieDSolverTest)

BOOST_AUTO_TEST_CASE(Constructor)
{
    GillespieDSolverTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new GillespieDSolverTester());
    BOOST_CHECK_EQUAL(*(t->replicateP), (uint)-1);
    BOOST_CHECK_EQUAL(*(t->parametersP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->resourcesP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->rngP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), 0);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), 0);
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
    GillespieDSolverTester * t = NULL;

    // Build a model with an valid dependency matrix.
    {
    uint ns         = 4;
    uint nr         = 6;
    uint c[]        = {  100, 200, 300, 400};
    uint o[]        = {  0,   1,   2,   3,   1,   2};
    double k[]      = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,
                         1,  -1,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,
                         0,   0,   0,   1,   0,  -1};
    uint D[]        = {  0,   1,   0,   0,   0,   1,
                         0,   0,   1,   0,   0,   0,
                         0,   0,   1,   0,   1,   0,
                         0,   0,   0,   1,   0,   1};
    BOOST_REQUIRE_NO_THROW(t=new GillespieDSolverTester());
    BOOST_REQUIRE_NO_THROW(t->buildModel(ns, nr, c, o, k, S, D));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensitiesP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], c[i]);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=0; i<nr; i++) BOOST_CHECK_EQUAL((*(t->propensitiesP))[i], 0.0);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->SP))[i], S[i]);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->DP))[i], D[i]);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    for (uint i=0; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionsP))[i], (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionArgsP), (void *)NULL);
    for (uint i=0; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionArgsP))[i], (void *)NULL);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[1]), 0U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s1i((*(t->propensityFunctionArgsP))[2]), 1U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s2i((*(t->propensityFunctionArgsP))[2]), 2U);
    BOOST_CHECK_EQUAL(t->SecondOrderSelfPropensityArgs_si((*(t->propensityFunctionArgsP))[3]), 3U);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[4]), 2U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s1i((*(t->propensityFunctionArgsP))[5]), 0U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s2i((*(t->propensityFunctionArgsP))[5]), 3U);
    BOOST_CHECK_CLOSE(t->ZerothOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[0]), 1.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[1]), 2.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[2]), 3.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderSelfPropensityArgs_k((*(t->propensityFunctionArgsP))[3]), 4.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[4]), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[5]), 6.0, 1e-9);
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
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[3], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[4], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[5], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][0], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][1], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][2], 5U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][0], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][1], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][2], 5U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[2][0], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[2][1], 4U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[3][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[3][1], 5U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[4][0], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[4][1], 4U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[5][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[5][1], 5U);
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    }
}

BOOST_AUTO_TEST_CASE(SetReactionModel)
{
    GillespieDSolverTester * t = NULL;
    lm::io::ReactionModel * m = NULL;

    // Build a model from a ReactionModel structure.
    {
    uint ns         = 4;
    uint nr         = 6;
    uint c[]        = {  100, 200, 300, 400};
    uint o[]        = {  0,   1,   2,   3,   1,   2};
    double k[]      = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,
                         1,  -1,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,
                         0,   0,   0,   1,   0,  -1};
    uint D[]        = {  0,   1,   0,   0,   0,   1,
                         0,   0,   1,   0,   0,   0,
                         0,   0,   1,   0,   1,   0,
                         0,   0,   0,   1,   0,   1};
    m = new ReactionModel();
    m->set_number_species(ns);
    m->set_number_reactions(nr);
    for (uint i=0; i<ns; i++) m->add_initial_species_count(c[i]);
    for (uint i=0; i<nr; i++)
    {
        m->add_reaction();
        m->mutable_reaction(i)->set_type(o[i]);
        m->mutable_reaction(i)->add_rate_constant(k[i]);
    }
    for (uint i=0; i<ns*nr; i++) m->add_stoichiometric_matrix(S[i]);
    for (uint i=0; i<ns*nr; i++) m->add_dependency_matrix(D[i]);
    BOOST_REQUIRE_NO_THROW(t=new GillespieDSolverTester());
    BOOST_REQUIRE_NO_THROW(t->setReactionModel(m));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_REQUIRE_NE(*(t->propensitiesP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], c[i]);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=0; i<nr; i++) BOOST_CHECK_EQUAL((*(t->propensitiesP))[i], 0.0);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->SP))[i], S[i]);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->DP))[i], D[i]);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    for (uint i=0; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionsP))[i], (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionArgsP), (void *)NULL);
    for (uint i=0; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionArgsP))[i], (void *)NULL);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[1]), 0U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s1i((*(t->propensityFunctionArgsP))[2]), 1U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s2i((*(t->propensityFunctionArgsP))[2]), 2U);
    BOOST_CHECK_EQUAL(t->SecondOrderSelfPropensityArgs_si((*(t->propensityFunctionArgsP))[3]), 3U);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[4]), 2U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s1i((*(t->propensityFunctionArgsP))[5]), 0U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s2i((*(t->propensityFunctionArgsP))[5]), 3U);
    BOOST_CHECK_CLOSE(t->ZerothOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[0]), 1.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[1]), 2.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[2]), 3.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderSelfPropensityArgs_k((*(t->propensityFunctionArgsP))[3]), 4.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[4]), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[5]), 6.0, 1e-9);
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
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[3], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[4], 2U);
    BOOST_CHECK_EQUAL((*(t->numberDependentReactionsP))[5], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][0], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][1], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[0][2], 5U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][0], 1U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][1], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[1][2], 5U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[2][0], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[2][1], 4U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[3][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[3][1], 5U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[4][0], 2U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[4][1], 4U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[5][0], 3U);
    BOOST_CHECK_EQUAL((*(t->dependentReactionsP))[5][1], 5U);
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    if (m != NULL) BOOST_REQUIRE_NO_THROW(delete m); m = NULL;
    }

    // Build a model from a ReactionModel structure with multiple k values per reaction.
    {
    uint ns         = 4;
    uint nr         = 6;
    uint c[]        = {  100, 200, 300, 400};
    uint o[]        = {  0,   1,   2,   3,   1,   2};
    double k[]      = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int S[]         = { -1,   1,   0,   0,   0,   0,
                         1,  -1,   0,   0,   0,   0,
                         0,   0,   1,   0,  -1,   0,
                         0,   0,   0,   1,   0,  -1};
    uint D[]        = {  0,   1,   0,   0,   0,   1,
                         0,   0,   1,   0,   0,   0,
                         0,   0,   1,   0,   1,   0,
                         0,   0,   0,   1,   0,   1};
    m = new ReactionModel();
    m->set_number_species(ns);
    m->set_number_reactions(nr);
    for (uint i=0; i<ns; i++) m->add_initial_species_count(c[i]);
    for (uint i=0; i<nr; i++)
    {
        m->add_reaction();
        m->mutable_reaction(i)->set_type(o[i]);
        for (uint j=0; j<o[i]+1; j++)
            m->mutable_reaction(i)->add_rate_constant(k[i]);
    }
    for (uint i=0; i<ns*nr; i++) m->add_stoichiometric_matrix(S[i]);
    for (uint i=0; i<ns*nr; i++) m->add_dependency_matrix(D[i]);
    BOOST_REQUIRE_NO_THROW(t=new GillespieDSolverTester());
    BOOST_REQUIRE_NO_THROW(t->setReactionModel(m));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensitiesP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->initialSpeciesCountsP))[i], c[i]);
    for (uint i=0; i<ns; i++) BOOST_CHECK_EQUAL((*(t->speciesCountsP))[i], 0U);
    for (uint i=0; i<nr; i++) BOOST_CHECK_EQUAL((*(t->propensitiesP))[i], 0.0);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->SP))[i], S[i]);
    for (uint i=0; i<ns*nr; i++) BOOST_CHECK_EQUAL((*(t->DP))[i], D[i]);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    for (uint i=0; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionsP))[i], (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionArgsP), (void *)NULL);
    for (uint i=0; i<nr; i++) BOOST_CHECK_NE((*(t->propensityFunctionArgsP))[i], (void *)NULL);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[1]), 0U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s1i((*(t->propensityFunctionArgsP))[2]), 1U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s2i((*(t->propensityFunctionArgsP))[2]), 2U);
    BOOST_CHECK_EQUAL(t->SecondOrderSelfPropensityArgs_si((*(t->propensityFunctionArgsP))[3]), 3U);
    BOOST_CHECK_EQUAL(t->FirstOrderPropensityArgs_si((*(t->propensityFunctionArgsP))[4]), 2U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s1i((*(t->propensityFunctionArgsP))[5]), 0U);
    BOOST_CHECK_EQUAL(t->SecondOrderPropensityArgs_s2i((*(t->propensityFunctionArgsP))[5]), 3U);
    BOOST_CHECK_CLOSE(t->ZerothOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[0]), 1.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[1]), 2.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[2]), 3.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderSelfPropensityArgs_k((*(t->propensityFunctionArgsP))[3]), 4.0, 1e-9);
    BOOST_CHECK_CLOSE(t->FirstOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[4]), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(t->SecondOrderPropensityArgs_k((*(t->propensityFunctionArgsP))[5]), 6.0, 1e-9);
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
    if (m != NULL) BOOST_REQUIRE_NO_THROW(delete m); m = NULL;
    }
}


BOOST_AUTO_TEST_CASE(DestroyModel)
{
    GillespieDSolverTester * t = NULL;

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
    BOOST_REQUIRE_NO_THROW(t=new GillespieDSolverTester());
    BOOST_REQUIRE_NO_THROW(t->buildModel(ns, nr, c, o, k, S, D));
    BOOST_CHECK_EQUAL(*(t->numberSpeciesP), ns);
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), nr);
    BOOST_CHECK_NE(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensitiesP), (void *)NULL);
    BOOST_CHECK_NE(*(t->SP), (void *)NULL);
    BOOST_CHECK_NE(*(t->DP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->propensityFunctionArgsP), (void *)NULL);

    BOOST_REQUIRE_NO_THROW(t->destroyModelT());
    BOOST_CHECK_EQUAL(*(t->initialSpeciesCountsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->speciesCountsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->propensitiesP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->SP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->DP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->propensityFunctionsP), (void *)NULL);
    BOOST_CHECK_EQUAL(*(t->propensityFunctionArgsP), (void *)NULL);

    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}


BOOST_AUTO_TEST_SUITE_END()

//virtual void setModel(uint numberSpeciesA, uint numberReactionsA, uint * initialSpeciesCountsA, uint * oA, double * kA, int * SA, uint * CA);
