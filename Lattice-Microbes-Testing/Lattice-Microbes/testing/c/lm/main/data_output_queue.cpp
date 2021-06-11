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

#include <queue>
#include "TestingHelper.h"
#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/main/DataOutputQueue.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

class DataOutputQueueTester : public lm::main::DataOutputQueue
{
public:
    DataOutputQueueTester():dataQueueP(&dataQueue){}
    std::queue<DataSet *> * dataQueueP;

    void popDataSet(byte ** data, size_t * size)
    {
        *data = dataQueueP->front()->data;
        *size = dataQueueP->front()->size;
        dataQueueP->pop();
    }
};

BOOST_AUTO_TEST_SUITE(DataOutputQueueTest)

BOOST_AUTO_TEST_CASE(PushDataSetBuffer)
{
    DataOutputQueueTester * t = NULL;
    byte * data = NULL;
    size_t dataSize = 0;

    // Test pushing a data set from a buffer.
    BOOST_REQUIRE_NO_THROW(t=new DataOutputQueueTester());
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 0);
    dataSize = 100;
    data = new byte[dataSize];
    for (uint i=0; i<dataSize; i++) data[i] = i;
    BOOST_CHECK_NO_THROW(t->pushDataSet(data, dataSize));
    if (data != NULL) delete[] data; data = NULL; dataSize = 0;
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 1);
    t->popDataSet(&data, &dataSize);
    BOOST_CHECK_EQUAL(dataSize, 100);
    for (uint i=0; i<dataSize; i++)
    {
        BOOST_CHECK_EQUAL(data[i], i);
    }
    if (t != NULL) delete t; t = NULL;

    // Test pushing mulitple data sets from a buffer.
    BOOST_REQUIRE_NO_THROW(t=new DataOutputQueueTester());
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 0);
    dataSize = 100;
    data = new byte[dataSize];
    for (uint i=0; i<dataSize; i++) data[i] = i;
    BOOST_CHECK_NO_THROW(t->pushDataSet(data, dataSize));
    if (data != NULL) delete[] data; data = NULL; dataSize = 0;
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 1);

    dataSize = 200;
    data = new byte[dataSize];
    for (uint i=0; i<dataSize; i++) data[i] = i+5;
    BOOST_CHECK_NO_THROW(t->pushDataSet(data, dataSize));
    if (data != NULL) delete[] data; data = NULL; dataSize = 0;
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 2);

    dataSize = 300;
    data = new byte[dataSize];
    for (uint i=0; i<dataSize; i++) data[i] = i+7;
    BOOST_CHECK_NO_THROW(t->pushDataSet(data, dataSize));
    if (data != NULL) delete[] data; data = NULL; dataSize = 0;
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 3);

    t->popDataSet(&data, &dataSize);
    BOOST_CHECK_EQUAL(dataSize, 100);
    for (uint i=0; i<dataSize; i++)
    {
        BOOST_CHECK_EQUAL(data[i], i);
    }

    t->popDataSet(&data, &dataSize);
    BOOST_CHECK_EQUAL(dataSize, 200);
    for (uint i=0; i<dataSize; i++)
    {
        BOOST_CHECK_EQUAL(data[i], i+5);
    }

    t->popDataSet(&data, &dataSize);
    BOOST_CHECK_EQUAL(dataSize, 300);
    for (uint i=0; i<dataSize; i++)
    {
        BOOST_CHECK_EQUAL(data[i], (i+7)%256);
    }


    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(PushDataSetMessage)
{
    DataOutputQueueTester * t = NULL;
    lm::io::SpeciesCounts * message = NULL;
    lm::io::SpeciesCounts * newMessage = NULL;
    size_t messageSize = 0;
    byte * payload = NULL;
    size_t payloadSize = 0;
    byte * data = NULL;
    size_t dataSize = 0;

    // Test pushing a data set from a buffer.
    BOOST_REQUIRE_NO_THROW(t=new DataOutputQueueTester());
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 0);
    message = new lm::io::SpeciesCounts();
    message->set_number_species(10);
    message->set_number_entries(33);
    for (int i=0; i<10*33; i++)
    {
        message->add_species_count(i);
    }
    for (int i=0; i<3; i++)
    {
        message->add_time(i);
    }
    messageSize = message->ByteSize();
    payloadSize = 300;
    payload = new byte[payloadSize];
    for (uint i=0; i<payloadSize; i++) payload[i] = i+7;
    BOOST_CHECK_NO_THROW(t->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, 117, message, payload, payloadSize));
    BOOST_CHECK_EQUAL(t->dataQueueP->size(), 1);
    t->popDataSet(&data, &dataSize);
    BOOST_CHECK_EQUAL(dataSize, sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize+sizeof(size_t)+payloadSize);
    BOOST_CHECK_EQUAL(*((uint *)&data[0]), 10);
    BOOST_CHECK_EQUAL(*((uint *)&data[sizeof(uint)]), 117);
    BOOST_CHECK_EQUAL(*((size_t *)&data[sizeof(uint)+sizeof(uint)]), messageSize);
    newMessage = new lm::io::SpeciesCounts();
    newMessage->ParseFromArray(&data[sizeof(uint)+sizeof(uint)+sizeof(size_t)], messageSize);
    BOOST_CHECK_EQUAL(newMessage->number_species(), message->number_species());
    BOOST_CHECK_EQUAL(newMessage->number_entries(), message->number_entries());
    BOOST_CHECK_EQUAL(newMessage->species_count_size(), message->species_count_size());
    for (int i=0; i<newMessage->species_count_size(); i++)
    {
        BOOST_CHECK_EQUAL(newMessage->species_count(i), message->species_count(i));
    }
    BOOST_CHECK_EQUAL(newMessage->time_size(), message->time_size());
    for (int i=0; i<newMessage->time_size(); i++)
    {
        BOOST_CHECK_EQUAL(newMessage->time(i), message->time(i));
    }
    BOOST_CHECK_EQUAL(*((size_t *)&data[sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize]), payloadSize);
    for (uint i=0; i<payloadSize; i++)
    {
        BOOST_CHECK_EQUAL(data[sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize+sizeof(size_t)+(size_t)i], payload[(size_t)i]);
    }
    if (message != NULL) delete message; message = NULL;
    if (newMessage != NULL) delete newMessage; newMessage = NULL;
    if (payload != NULL) delete[] payload; payload = NULL; payloadSize = 0;
    if (data != NULL) delete[] data; data = NULL; dataSize = 0;
    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_SUITE_END()
