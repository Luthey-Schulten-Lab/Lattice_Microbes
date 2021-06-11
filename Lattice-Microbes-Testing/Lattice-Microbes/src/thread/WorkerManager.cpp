/*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
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

#include <list>
#include "config.h"
#include "core/Print.h"
#include "thread/Thread.h"
#include "thread/Worker.h"
#include "thread/WorkerManager.h"

namespace lm {
namespace thread {

WorkerManager WorkerManager::instance;

WorkerManager * WorkerManager::getInstance()
{
    return &instance;
}

WorkerManager::WorkerManager()
{
    // Create the mutex.
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_init(&mutex, NULL));
}

WorkerManager::~WorkerManager()
{
    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_mutex_destroy(&mutex));
}

void WorkerManager::addWorker(Worker * worker)
{
    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));
    workers.push_back(worker);
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex
}

void WorkerManager::removeWorker(Worker * worker)
{
    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));
    workers.remove(worker);
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex
}

void WorkerManager::abortWorkers()
{
    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));

    // Abort the worker threads.
    for (list<Worker *>::iterator it=workers.begin(); it != workers.end(); it++)
    {
        (*it)->abort();
    }

    // Wait for the worker threads to stop.
    for (list<Worker *>::iterator it=workers.begin(); it != workers.end(); it++)
    {
        (*it)->stop();
    }

    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex
}

void WorkerManager::checkpointWorkers()
{
    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));

    // Tell the worker threads to checkpoint.
    for (list<Worker *>::iterator it=workers.begin(); it != workers.end(); it++)
    {
        (*it)->checkpoint();
    }

    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex
}

void WorkerManager::stopWorkers()
{
    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));

    // Wait for the worker threads to stop.
    for (list<Worker *>::iterator it=workers.begin(); it != workers.end(); it++)
    {
        Print::printf(Print::DEBUG, "Stopping worker thread: %u", (*it)->getId());
        (*it)->stop();
        Print::printf(Print::DEBUG, "Stopped worker thread: %u", (*it)->getId());
    }

    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex
}

}
}
