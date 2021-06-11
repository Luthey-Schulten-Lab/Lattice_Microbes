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

#ifndef LM_THREAD_THREAD_H_
#define LM_THREAD_THREAD_H_

#include <cstdlib>
#include <pthread.h>
#include "core/Exceptions.h"

namespace lm {
namespace thread {

/// @class PthreadException
/// @brief An Exception class for handling pthread exceptions.
class PthreadException : public Exception
{
public:
    PthreadException(int error, const char * file, const int line) : Exception("Error in pthread library", error, file, line) {}
};

/// @class Thread
/// @brief A base class wrapping pthreads.
class Thread
{
public:
    /// @brief Creates a pthread locking mechanism and initializes "Thread"
    Thread();
    /// @brief Destory the Thread
    virtual ~Thread();
    
    /// @brief If no thread exists, creates a new thread and begins execution
    virtual void start();
    /// @brief Joins the thread with the parent waiting if necessary
    virtual void stop();
    /// @brief Wakes a sleeping thead
    virtual void wake()=0;
    /// @brief Returns the pthread based id for the Thread
    /// @return pthread assigned id
    virtual pthread_t getId() {return threadId;}
    /// @brief Binds the thread to a CPU core
    /// @param The cpu core for which to bind the thread
    virtual void setAffinity(int cpuNumber);

protected:
    virtual int run()=0;

private:
    static void * start_thread(void * obj);

protected:
    pthread_mutex_t controlMutex;   // thread locking mechanism
    pthread_t threadId;             // thread ID
    volatile bool running;          // thread status
    int cpuNumber;                  // which cpu thread is running on
};

}
}


/// @def PTHREAD_EXCEPTION_CHECK
/// @brief Exception wrapping for the pthread api.
#define PTHREAD_EXCEPTION_CHECK(pthread_call) {int _pthread_ret_=pthread_call; if (_pthread_ret_ != 0) throw lm::thread::PthreadException(_pthread_ret_,__FILE__,__LINE__);}
#define PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_call) {int _pthread_ret_=pthread_call; if (_pthread_ret_ != 0) {fprintf(stderr, "%s:%d: Pthread error, returned %d\n", __FILE__, __LINE__, _pthread_ret_); abort();}}

#endif
