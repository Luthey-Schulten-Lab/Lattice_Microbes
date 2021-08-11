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

#ifndef LM_MAIN_CHECKPOINTSIGNALER
#define LM_MAIN_CHECKPOINTSIGNALER

#include <ctime>
#include <pthread.h>
#include "thread/Thread.h"
#include "thread/Worker.h"

using lm::thread::PthreadException;
using lm::thread::Worker;


namespace lm {
namespace main {

/// @class CheckpointSignaler
/// @brief A type of worker thread that checkpoints at a specified interval
class CheckpointSignaler : public Worker
{
public:
    /// @brief Creates a new CheckpointSignaler Worker thread
	CheckpointSignaler();
    virtual ~CheckpointSignaler();
    
    /// @brief Tells the thread to start checkpointing every checkpointInterval seconds
    /// @param checkpointInterval How often in seconds to checkpoint
    virtual void startCheckpointing(time_t checkpointInterval);
    /// @brief Tell the thread to stop checkpointing
    virtual void stopCheckpointing();
    /// @brief Wake the thread if inactive
    virtual void wake();

private:
    pthread_cond_t controlChange;
    time_t checkpointInterval;
    struct timespec nextCheckpoint;

protected:
    /// @brief A method to actually run a timestep while checkpointing every "checkpointInterval" time increment
    /// @return 0 on success, -1 on failure
    virtual int run();
    /// @brief Set the next checkpoint time for this thread
    void setNextCheckpointTime();
};

}
}


#endif
