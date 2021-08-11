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

#ifndef LM_THREAD_WORKER_MANAGER_H_
#define LM_THREAD_WORKER_MANAGER_H_

#include <list>
#include "thread/Thread.h"
#include "thread/Worker.h"

using std::list;

namespace lm {
namespace thread {

/// @class WorkerManager
/// @brief A singleton manager that creates and manages workers that run on a processing node
class WorkerManager
{
public:
    /// @brief Get the global worker manager
    /// @return globalWorkerManager Global thread manager instance handle
    static WorkerManager * getInstance();

private:
    static WorkerManager instance;  // Global thread manager

public:
    /// @brief Create the WorkerManager
    WorkerManager();
    /// @brief Destroy the WorkerManager
    ~WorkerManager();

    /// @brief Adds a worker to the manager
    /// @param worker Worker to add
    void addWorker(Worker * worker);
    /// @brief Removes a worker from the manager
    /// @param worker Worker handle for which to remove
    void removeWorker(Worker * worker);
    /// @brief Aborts all the worker threads
    void abortWorkers();
    /// @brief Causes all the worker threads to checkpoint (checkpointing is currently unimplemented)
    void checkpointWorkers();
    /// @brief Stops all the worker threads (e.g. merges them with the master)
    void stopWorkers();

private:
    pthread_mutex_t mutex;      // Main mutex for performing thread operations
    list<Worker *> workers;     // List of all workers
};

}
}


#endif
