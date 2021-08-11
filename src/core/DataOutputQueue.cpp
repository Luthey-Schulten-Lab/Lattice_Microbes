/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
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
#include <google/protobuf/stubs/common.h>
#include "config.h"
#include "core/Exceptions.h"
#include "core/Types.h"
#include "core/DataOutputQueue.h"

namespace lm {
namespace main {

DataOutputQueue * DataOutputQueue::instance = NULL;

void DataOutputQueue::setInstance(DataOutputQueue * instance)
{
    DataOutputQueue::instance = instance;
}

DataOutputQueue * DataOutputQueue::getInstance()
{
    return instance;
}

DataOutputQueue::DataOutputQueue()
{
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_init(&dataMutex, NULL));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_init(&returnMutex, NULL));
    PTHREAD_EXCEPTION_CHECK(pthread_cond_init(&returnAvailable, NULL));
}

DataOutputQueue::~DataOutputQueue()
{
    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_mutex_destroy(&dataMutex));
    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_mutex_destroy(&returnMutex));
    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_cond_destroy(&returnAvailable));
}

void DataOutputQueue::pushDataSet(const H5MetaData md, uint replicate, void * payload)
{
    size_t buffer_sz = sizeof(md)+md.payloadSize;
    DataSet* dataSet = new DataSet(buffer_sz);
    H5MetaData *md_ptr = reinterpret_cast<H5MetaData *>(dataSet->data);
    unsigned char *data_ptr = dataSet->data + sizeof(md);
    memcpy(md_ptr, &md, sizeof(md));
    memcpy(data_ptr, payload, md.payloadSize);
    md_ptr->type = ARBITRARY_H5;
    md_ptr->replicate = replicate;
    pushDataSet(dataSet);
}



void DataOutputQueue::pushDataSet(uint type, uint replicate, ::google::protobuf::Message * message, void * payload, size_t payloadSize, void (*payloadSerializer)(void *, void *, size_t))
{
    // Get the size of the message.
    size_t messageSize = 0;
    if (message != NULL) messageSize = message->ByteSize();

    // Allocate a new data set buffer for the data.
    DataSet * dataSet = new DataSet(sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize+sizeof(size_t)+payloadSize);

    // Serialize the data into the buffer.
    byte * buffer = dataSet->data;
    memcpy(buffer, &type, sizeof(uint));
    buffer += sizeof(uint);
    memcpy(buffer, &replicate, sizeof(uint));
    buffer += sizeof(uint);

    // Serialize the message.
    memcpy(buffer, &messageSize, sizeof(size_t));
    buffer += sizeof(size_t);
    if (messageSize > 0)
    {
        message->SerializeToArray(buffer, messageSize);
        buffer += messageSize;
    }

    // Serialize the payload.
    memcpy(buffer, &payloadSize, sizeof(size_t));
    buffer += sizeof(size_t);
    if (payloadSize > 0)
    {
        if (payloadSerializer == NULL)
            memcpy(buffer, payload, payloadSize);
        else
            payloadSerializer(buffer, payload, payloadSize);
    }

    // Add the data set to the queue.
    pushDataSet(dataSet);
}

void DataOutputQueue::pushDataSet(void * data, size_t dataSize)
{
    pushDataSet(new DataSet(data, dataSize));
}

void DataOutputQueue::pushDataSet(DataSet * dataSet)
{
    //// BEGIN CRITICAL SECTION: dataMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&dataMutex));
    if (dataSet != NULL) dataQueue.push(dataSet);
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&dataMutex));
    //// END CRITICAL SECTION: dataMutex
}

}
}

