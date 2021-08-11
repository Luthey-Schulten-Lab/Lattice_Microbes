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

#include <csignal>
#include "config.h"
#include "core/Print.h"
#include "core/LocalDataOutputWorker.h"
#include "core/Globals.h"
#include "core/SignalHandler.h"
#include "thread/Thread.h"
#include "thread/WorkerManager.h"

using lm::thread::PthreadException;
using lm::thread::WorkerManager;

namespace lm {
namespace main {


SignalHandler::SignalHandler()
{
    // Block signals from interrupting the calling thread.
    sigemptyset(&signalMask);
    sigaddset(&signalMask, SIGINT);
    sigaddset(&signalMask, SIGTERM);
    sigaddset(&signalMask, SIGUSR1);
    sigaddset(&signalMask, SIGUSR2);
    PTHREAD_EXCEPTION_CHECK(pthread_sigmask(SIG_BLOCK, &signalMask, NULL));
}

SignalHandler::~SignalHandler()
{
}

void SignalHandler::wake()
{
    PTHREAD_EXCEPTION_CHECK(pthread_kill(threadId, SIGUSR1));
}

int SignalHandler::run()
{
    try
    {
        Print::printf(Print::DEBUG, "Signal handler thread running.");

        bool looping = true;
        while (looping)
        {
            // Wait for a signal.
            int signalCaught;
            if (sigwait(&signalMask, &signalCaught) != 0) throw lm::Exception("Error waiting for signal.");

            Print::printf(Print::DEBUG, "Received signal %d, control signal is %d.", signalCaught,SIGUSR1);

            // If this was our control signal, check for control changes.
            if (signalCaught == SIGUSR1)
            {
                //// BEGIN CRITICAL SECTION: controlMutex
                PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));

                // If we should abort, stop looping.
                if (aborted || !running) looping=false;

                PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
                //// END CRITICAL SECTION: controlMutex
            }

            // Otherwise, see what we need to do.
            if ((signalCaught == SIGUSR1 && looping == true) || signalCaught == SIGUSR2 || signalCaught == SIGINT || signalCaught == SIGTERM)
            {
                // Set the global abort flag.
                Print::printf(Print::WARNING, "Global abort signaled.");
                globalAbort = true;
            }
        }
        Print::printf(Print::DEBUG, "Signal handler thread finished.");
        return 0;
    }
    catch (lm::thread::PthreadException e)
    {
        Print::printf(Print::FATAL, "Pthread exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (lm::Exception e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (std::exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown Exception during execution (%s:%d)", __FILE__, __LINE__);
    }

    return -1;
}

}
}
