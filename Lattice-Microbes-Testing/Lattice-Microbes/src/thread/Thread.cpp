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

#include "config.h"
#include <pthread.h>
#include "core/Print.h"
#include "thread/Thread.h"

namespace lm {
namespace thread {

Thread::Thread()
:threadId(0),running(false),cpuNumber(-1)
{
    // Create the control mutex.
	pthread_mutexattr_t attr;
	PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_init(&attr));
	PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_NORMAL));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_init(&controlMutex, &attr));
    PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_destroy(&attr));
}

Thread::~Thread()
{
    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_mutex_destroy(&controlMutex));
}

void Thread::setAffinity(int cpuNumber)
{
    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));

    // Set the assigned cpu.
	this->cpuNumber = cpuNumber;

	// If we are already running, reset the processor affinity.
    if (running)
    {
		#if defined(LINUX) && !defined(ARM)
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(cpuNumber, &cpuset);
		if (pthread_setaffinity_np(threadId, sizeof(cpu_set_t), &cpuset) != 0)
			Print::printf(Print::WARNING, "Could not bind thread %u to CPU core %d", threadId, cpuNumber);
		#endif
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex
}

void Thread::start()
{
    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    if (!running)
    {
        running=true;
        pthread_attr_t attr;
        PTHREAD_EXCEPTION_CHECK(pthread_attr_init(&attr));
        PTHREAD_EXCEPTION_CHECK(pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE));
        PTHREAD_EXCEPTION_CHECK(pthread_create(&threadId, &attr, &Thread::start_thread, this));
        PTHREAD_EXCEPTION_CHECK(pthread_attr_destroy(&attr));

    	// Set the processor affinity, if we have a cpu assigned.
        if (cpuNumber >= 0)
        {
			#if defined(LINUX) && !defined(ARM)
			cpu_set_t cpuset;
			CPU_ZERO(&cpuset);
			CPU_SET(cpuNumber, &cpuset);
			if (pthread_setaffinity_np(threadId, sizeof(cpu_set_t), &cpuset) != 0)
				Print::printf(Print::WARNING, "Could not bind thread %u to CPU core %d", threadId, cpuNumber);
			#endif
        }
        Print::printf(Print::DEBUG, "Started thread %u.", threadId);
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex
}

void * Thread::start_thread(void * obj)
{
    unsigned long long ret = (reinterpret_cast<Thread *>(obj))->run();
    pthread_exit((void *)ret);
}

void Thread::stop()
{
    bool waitForThread=false;

    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    if (running)
    {
        Print::printf(Print::DEBUG, "Stopping thread %u.", threadId);
        running = false;
        waitForThread = true;
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex

    // Wait for the thread to stop, if we really stopped it.
    if (waitForThread)
    {
        // Try to wake up the thread in case it is sleeping.
        wake();

        // Join with the thread.
        void * ret;
        PTHREAD_EXCEPTION_CHECK(pthread_join(threadId, &ret));
        Print::printf(Print::DEBUG, "Thread %u stopped.", threadId);
    }
}

}
}
