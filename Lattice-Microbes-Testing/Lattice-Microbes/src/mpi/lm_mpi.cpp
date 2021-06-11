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

#include <cstdio>
#include <iostream>
#include <mpi.h>
#include "config.h"
#include "core/Exceptions.h"
#include "mpi/lm_mpi.h"

namespace lm {

MPIException::MPIException(int error)
:Exception("MPI Exception"),errorCode(error)
{
    char errorMessage[MPI_MAX_ERROR_STRING];
    int errorMessageLength=0;
    MPI_Error_string(error, errorMessage, &errorMessageLength);
    snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s", errorMessage);
}

int MPI::version = -1;
int MPI::subversion = -1;
int MPI::threadSupport = -1;
int MPI::worldSize = -1;
int MPI::worldRank = -1;

void MPI::init(int argc, char** argv)
{
    // Get the MPI version info.
    MPI_EXCEPTION_CHECK(MPI_Get_version(&version, &subversion));

    // Initialize the library.
    MPI_EXCEPTION_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadSupport));

    // Set the error handler to return errors for the world.
    MPI_EXCEPTION_CHECK(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN));

    // Get the size of the world.
    MPI_EXCEPTION_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &worldSize));

    // Get our rank.
    MPI_EXCEPTION_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &worldRank));
}

void MPI::printCapabilities()
{
    std::cout << "MPI version " << lm::MPI::version << "." << lm::MPI::subversion << " with thread support ";
    if (lm::MPI::threadSupport == MPI_THREAD_SINGLE)
        std::cout << "MPI_THREAD_SINGLE." << std::endl;
    else if (lm::MPI::threadSupport == MPI_THREAD_FUNNELED)
        std::cout << "MPI_THREAD_FUNNELED." << std::endl;
    else if (lm::MPI::threadSupport == MPI_THREAD_SERIALIZED)
        std::cout << "MPI_THREAD_SERIALIZED." << std::endl;
    else if (lm::MPI::threadSupport == MPI_THREAD_MULTIPLE)
        std::cout << "MPI_THREAD_MULTIPLE." << std::endl;
    else
        std::cout << "UNKNOWN." << std::endl;
    std::cout << "Running with " << worldSize << " process(es)." << std::endl;
    std::cout << std::endl;
}


void MPI::finalize()
{
    // Close own the MPI library.
    MPI_EXCEPTION_CHECK(MPI_Finalize());
}

}

