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

#ifndef LM_MAIN_LOCALDATAOUTPUTWORKER
#define LM_MAIN_LOCALDATAOUTPUTWORKER

#include <queue>
#include <cstring>
#include "Types.h"
#include "io/SimulationFile.h"
#include "core/DataOutputQueue.h"
#include "thread/Thread.h"
#include "thread/Worker.h"

using std::queue;
using lm::thread::PthreadException;
using lm::thread::Worker;

using lm::io::hdf5::SimulationFile;

namespace lm {
namespace main {

/// @class LocalDataOutputWorker
/// @brief A class that handles output on a particular core
class LocalDataOutputWorker : public Worker, public DataOutputQueue
{
public:
    /// @brief Create a new LocalDataOuputWorker
    /// @param file The file in which to save the current data stream
    LocalDataOutputWorker(SimulationFile * file);
    virtual ~LocalDataOutputWorker();

    /// @brief Push a DataSet for writing
    /// @param dataSet DataSet to add to the queue for writing
    virtual void pushDataSet(DataSet * dataSet);

    /// @brief Wake the thread from sleeping state
    virtual void wake();
    /// @brief Abort the thread
    virtual void abort();
    /// @brief Tell the thread to checkpoint at the next available times
    virtual void checkpoint();

protected:
    /// @brief Start running the thread
    /// @return 0 on success, -1 on failure
    virtual int run();

private:

    SimulationFile * file;
    pthread_cond_t dataAvailable;
    bool shouldCheckpoint;
    bool shouldAbort;
};

}
}


#endif
