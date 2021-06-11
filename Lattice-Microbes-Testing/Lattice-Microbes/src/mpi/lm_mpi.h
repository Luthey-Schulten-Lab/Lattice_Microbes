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

#ifndef LM_MPI_H_
#define LM_MPI_H_

#include <mpi.h>
#include "core/Exceptions.h"

void MPIErrorHandler(MPI_Comm *, int *rc, ...);

namespace lm {

/// @class MPIException
/// @brief An Exception thrown when an MPI call fails.
class MPIException : public Exception
{
public:
    int errorCode;              // Error code thrown by MPI
    
    /// @brief Create the MPIException
    /// @param error MPI error code
    MPIException(int error);
};

/// @class MPI
/// @brief Handles the MPI capabilities and properties of the simulation.
class MPI
{
public:
    static int version;             // MPI main version
    static int subversion;          // MPI minor version
    static int threadSupport;       // MPI level of thread support
    static int worldSize;           // Size of MPI groups in the simulation
    static int worldRank;           // Number of MPI processes in the simulation
    static const int MASTER=0;

    // MPI messages.
    static const int MSG_RUN_SIMULATION         = 1;
    static const int MSG_SIMULATION_FINISHED    = 2;
    static const int MSG_OUTPUT_DATA_STATIC     = 10;
    static const int MSG_EXIT                   = 99;

    static const int OUTPUT_DATA_STATIC_MAX_SIZE    = 10*1024*1024;

    /// @brief Initialize the MPI runtime
    /// @param argc Command line number of arguments
    /// @param argv Command line array of arguments
    static void init(int argc, char** argv);
    /// @brief Print the capabilities of the current processing device
    static void printCapabilities();
    /// @brief Close and cleanup the MPI runtime
    static void finalize();
};

}

/// @def MPI_EXCEPTION_CHECK
/// @brief A macro to encapsulate calling an MPI function and resulting error handling
#define MPI_EXCEPTION_CHECK(mpi_call) {int _mpi_ret_=mpi_call; if (_mpi_ret_ != MPI_SUCCESS) throw lm::MPIException(_mpi_ret_);}


#endif /*LM_MPI_H_*/
