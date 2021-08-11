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

#ifndef LM_MAIN_MPIREMOTEDATAOUTPUTQUEUE
#define LM_MAIN_MPIREMOTEDATAOUTPUTQUEUE

#include <queue>
#include "thread/Thread.h"
#include "core/DataOutputQueue.h"

using std::queue;
using lm::thread::PthreadException;

namespace lm {
namespace main {

/// @class MPIRemoteDataOutputQueue
/// @brief A DataOutputQueue that works for MPI jobs across various different nodes that outputs only on the "master" thread
class MPIRemoteDataOutputQueue : public DataOutputQueue
{
public:
    /// @brief Create a new MPIRemoteDataOutputQueue
    MPIRemoteDataOutputQueue() {}
    virtual ~MPIRemoteDataOutputQueue() {}

    /// @brief Tell if the queue is empty
    /// @return true if the queue is empty
    virtual bool isEmpty();
    /// @brief Pop data from queue and write to the buffer
    /// @param buffer The data buffer into which the next DataSet object should be written
    /// @param bufferSize The size of the buffer 
    /// @return Number of bytes written to the buffer
    virtual size_t popDataSetIntoBuffer(void * buffer, size_t bufferSize);
};

}
}


#endif
