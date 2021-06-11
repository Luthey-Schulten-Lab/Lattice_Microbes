/*
 * University of Illinois Open Source License
 * Copyright 2012-2018 Luthey-Schulten Group,
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

#include <cerrno>
#include <ctime>
#include "config.h"
#if defined(MACOSX)
#include <sys/time.h>
#endif
#include <pthread.h>
#include "core/Print.h"
#include "core/CheckpointSignaler.h"
#include "thread/Thread.h"
#include "thread/Worker.h"
#include "thread/WorkerManager.h"

using lm::thread::PthreadException;
using lm::thread::WorkerManager;

namespace lm {
namespace main {

CheckpointSignaler::CheckpointSignaler()
: checkpointInterval(0)
{
	nextCheckpoint.tv_sec = 0;
	nextCheckpoint.tv_nsec = 0;
	PTHREAD_EXCEPTION_CHECK(pthread_cond_init(&controlChange, NULL));
}

CheckpointSignaler::~CheckpointSignaler()
{
    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_cond_destroy(&controlChange));
}

void CheckpointSignaler::startCheckpointing(time_t checkpointInterval)
{
	Print::printf(Print::INFO, "Creating a checkpoint file every %d seconds.", checkpointInterval);

    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
	this->checkpointInterval = checkpointInterval;
	setNextCheckpointTime();
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex

	wake();
}

void CheckpointSignaler::stopCheckpointing()
{
    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
	this->checkpointInterval = 0;
	setNextCheckpointTime();
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex

	wake();
}

void CheckpointSignaler::wake()
{
    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&controlChange));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex
}

int CheckpointSignaler::run()
{
    try
    {
        Print::printf(Print::DEBUG, "Checkpoint signaler thread running.");

        bool looping = true;
        while (looping)
        {
        	bool doCheckpoint = false;

            //// BEGIN CRITICAL SECTION: controlMutex
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));

            // See if we are checkpointing.
            if (checkpointInterval > 0)
            {
				// Sleep until our next checkpoint interval or we receive a control signal.
            	Print::printf(Print::DEBUG, "Checkpoint signaler thread waiting until the next checkpoint time.");
				int ret=pthread_cond_timedwait(&controlChange, &controlMutex, &nextCheckpoint);
				if (ret != 0 && ret != ETIMEDOUT) throw lm::thread::PthreadException(ret,__FILE__,__LINE__);
            }
            else
            {
                // Otherwise we are not checkpointing, so just sleep until our next control change.
            	Print::printf(Print::DEBUG, "Checkpoint signaler thread waiting for a control change.");
            	PTHREAD_EXCEPTION_CHECK(pthread_cond_wait(&controlChange, &controlMutex));
            }

            Print::printf(Print::DEBUG, "Checkpoint signaler thread awoken.");

            // If we should abort or have been stopped, stop looping.
            if (aborted || !running) looping=false;

            // Otherwise, see if we need to perform a checkpoint.
            else if (checkpointInterval > 0)
            {
            	// Get the current time.
				#if defined(LINUX)
				struct timespec now;
				clock_gettime(CLOCK_REALTIME, &now);
				#elif defined(MACOSX)
				struct timeval now;
				gettimeofday(&now, NULL);
				#endif

            	// See if we are past the checkpoint time.
				if (now.tv_sec >= nextCheckpoint.tv_sec)
				{
					// Mark that we should perform the checkpoint once we are out of the critical section.
					doCheckpoint = true;
				}
            }

			PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
			//// END CRITICAL SECTION: controlMutex

			// See if we should perform the checkpoint.
			if (doCheckpoint)
			{
				// Perform the checkpoint.
				Print::printf(Print::DEBUG, "Performing checkpoint.");
				WorkerManager::getInstance()->checkpointWorkers();

				// Update the next checkpoint time.
				setNextCheckpointTime();
			}
        }
        Print::printf(Print::DEBUG, "Checkpoint signaler thread finished.");
        return 0;
    }
    catch (lm::thread::PthreadException & e)
    {
        Print::printf(Print::FATAL, "Pthread exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (lm::Exception & e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (std::exception & e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown Exception during execution (%s:%d)", __FILE__, __LINE__);
    }

    return -1;
}

void CheckpointSignaler::setNextCheckpointTime()
{
	// Set the next checkpoint time.
	#if defined(LINUX)
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);
	nextCheckpoint.tv_sec = now.tv_sec+checkpointInterval;
	#elif defined(MACOSX)
	struct timeval now;
	gettimeofday(&now, NULL);
	nextCheckpoint.tv_sec = now.tv_sec+checkpointInterval;
	#endif
}


}
}
