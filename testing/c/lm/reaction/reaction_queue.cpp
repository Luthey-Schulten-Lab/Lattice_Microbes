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

#include "TestingHelper.h"
#include "lm/Types.h"
#include "lm/reaction/ReactionQueue.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

class ReactionQueueTester : public lm::reaction::ReactionQueue
{
public:
    ReactionQueueTester(uint numberReactionsA):ReactionQueue(numberReactionsA),numberReactionsP(&numberReactions),reactionEventsP(&reactionEvents),reactionQueueP(&reactionQueue),reactionPositionsP(&reactionPositions){}

    uint *numberReactionsP;
    ReactionEvent ** reactionEventsP;
    uint ** reactionQueueP;
    uint ** reactionPositionsP;
};

BOOST_AUTO_TEST_SUITE(GillespieDSolverTest)

BOOST_AUTO_TEST_CASE(Constructor)
{
    ReactionQueueTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ReactionQueueTester(6));
    BOOST_CHECK_EQUAL(*(t->numberReactionsP), 6);
    BOOST_CHECK_NE(*(t->reactionEventsP), (void *)NULL);
    BOOST_CHECK_NE(*(t->reactionQueueP), (void *)NULL);
    BOOST_CHECK_NE(*(t->reactionPositionsP), (void *)NULL);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[0], 0);
    for (uint i=0; i<6; i++)
    {
        BOOST_CHECK_EQUAL((*(t->reactionEventsP))[i].time, INFINITY);
        BOOST_CHECK_EQUAL((*(t->reactionEventsP))[i].propensity, 0.0);
        BOOST_CHECK_EQUAL((*(t->reactionQueueP))[i+1], i);
        BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[i], i+1);
    }

    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}

BOOST_AUTO_TEST_CASE(GetNextReaction)
{
    ReactionQueueTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ReactionQueueTester(6));
    BOOST_CHECK_EQUAL(t->getNextReaction(), 0);
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}

BOOST_AUTO_TEST_CASE(GetReactionEvent)
{
    ReactionQueueTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ReactionQueueTester(6));
    for (uint i=0; i<6; i++)
    {
        BOOST_CHECK_EQUAL(t->getReactionEvent(i).time, INFINITY);
        BOOST_CHECK_EQUAL(t->getReactionEvent(i).propensity, 0.0);
    }
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}


BOOST_AUTO_TEST_CASE(UpdateReactionEvent)
{
    ReactionQueueTester * t = NULL;

    BOOST_REQUIRE_NO_THROW(t=new ReactionQueueTester(6));

    // Update reaction 0, keeping it in place.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(0, 15.765, 8.999));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 15.765, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 8.999, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 0);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 1);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[2].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[2].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 2);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 4);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[5].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[5].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 5);


    // Update reaction 5, moving it tow the top.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(5, 1.786, 32.78));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 15.765, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 8.999, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 0);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 1);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[2].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[2].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 2);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Update reaction 2 without changing its queue position.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(2, 100.01, 1e4));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 15.765, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 8.999, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 0);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 100.01, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e4, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 2);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Update reaction 0, moving it down one level.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(0, 1752.0, 0.567));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 1752.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 0.567, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 0);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[1].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 100.01, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e4, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 2);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[3].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[4].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Set the other reaction keeping them in place.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(1, 5.5, 1.01));
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(3, 6789.0, 1.02));
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(4, 1234.0, 1.03));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 1752.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 0.567, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 0);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].time, 5.5, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].propensity, 1.01, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 100.01, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e4, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 2);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].time, 6789.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].propensity, 1.02, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].time, 1234.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].propensity, 1.03, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Update reaction 3 moving it up one level.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(3, 3.4, 1.09));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 1752.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 0.567, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 0);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].time, 5.5, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].propensity, 1.01, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 100.01, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e4, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 2);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].time, 3.4, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].propensity, 1.09, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 3);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].time, 1234.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].propensity, 1.03, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Update reaction 3 moving it up down level to the first child.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(3, 2500.07, 1.09));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 1752.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 0.567, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 0);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].time, 5.5, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].propensity, 1.01, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 100.01, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e4, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 2);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].time, 2500.07, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].propensity, 1.09, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].time, 1234.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].propensity, 1.03, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Update reaction 1 moving it down a level to the second child.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(1, 1782.0, 11.0));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 1752.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 0.567, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 0);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].time, 1782.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].propensity, 11.0, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 100.01, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e4, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 2);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].time, 2500.07, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].propensity, 1.09, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].time, 1234.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].propensity, 1.03, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Update reaction 2 moving it down a level to the only child.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(2, 2011, 1e-6));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 1752.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 0.567, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 0);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].time, 1782.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].propensity, 11.0, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 2011, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e-6, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 2);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].time, 2500.07, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].propensity, 1.09, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].time, 1234.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].propensity, 1.03, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 4);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].time, 1.786, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[5].propensity, 32.78, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 5);

    // Update reaction 5 to Infinity, moving it to the bottom.
    BOOST_REQUIRE_NO_THROW(t->updateReactionEvent(5, INFINITY, 0.0));

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].time, 1752.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[0].propensity, 0.567, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[0], 3);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[3], 0);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].time, 1782.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[1].propensity, 11.0, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[1], 2);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[2], 1);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].time, 2011, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[2].propensity, 1e-6, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[2], 6);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[6], 2);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].time, 2500.07, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[3].propensity, 1.09, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[3], 4);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[4], 3);

    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].time, 1234.0, 1e-9);
    BOOST_CHECK_CLOSE((*(t->reactionEventsP))[4].propensity, 1.03, 1e-9);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[4], 1);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[1], 4);

    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[5].time, INFINITY);
    BOOST_CHECK_EQUAL((*(t->reactionEventsP))[5].propensity, 0.0);
    BOOST_CHECK_EQUAL((*(t->reactionPositionsP))[5], 5);
    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[5], 5);

    BOOST_CHECK_EQUAL((*(t->reactionQueueP))[0], 0);
    if (t != NULL) BOOST_REQUIRE_NO_THROW(delete t); t = NULL;
}


BOOST_AUTO_TEST_SUITE_END()
