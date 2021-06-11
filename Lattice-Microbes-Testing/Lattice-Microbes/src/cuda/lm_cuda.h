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

#ifndef LM_CUDAHELPER_H_
#define LM_CUDAHELPER_H_

#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include "core/Exceptions.h"

namespace lm
{

/**
 * CUDA runtime exception.
 */
class CUDAException : public Exception
{
public:
	CUDAException(cudaError_t error, const char * file, const int line) : Exception("CUDA error",cudaGetErrorString(error),file,line) {}
};

/**
 * CUDA driver exception.
 */
class CUDADriverException : public Exception
{
public:
	CUDADriverException(CUresult error, const char * file, const int line) : Exception("CUDA driver error", error,file,line) {}
};

/**
 * Class for accessing CUDA functions.
 */
class CUDA
{
public:
	static CUresult _cuda_dret_;
	
	static int getNumberDevices();
	static void setCurrentDevice(int device);
	static int getCurrentDevice();
	static size_t getFreeMemory(int device);
	static void getComputeCapabilities(int device, int * major, int * minor);
	static void printCapabilities(int device);
    static std::string getCapabilitiesString(int device);
};

}

// Macros for wrapping CUDA calls within an exception handler.
#define CUDA_EXCEPTION_CHECK(cuda_call) {if ((cuda_call) != cudaSuccess) throw lm::CUDAException(cudaGetLastError(),__FILE__,__LINE__);}
#define CUDA_EXCEPTION_DRIVER_CHECK(cuda_call) {if ((::lm::CUDA::_cuda_dret_=cuda_call) != CUDA_SUCCESS) throw lm::CUDADriverException(::lm::CUDA::_cuda_dret_,__FILE__,__LINE__);}
#define CUDA_EXCEPTION_EXECUTE(cuda_exec) {cuda_exec; if (cudaPeekAtLastError() != cudaSuccess) throw lm::CUDAException(cudaGetLastError(),__FILE__,__LINE__);}
#define CUDA_EXCEPTION_CHECK_NOTHROW(cuda_call) do { \
    if ((cuda_call) != cudaSuccess) { \
        fprintf(stderr, "Cuda error in destructor: %s @ %s:%d\n", cudaGetErrorString(cudaGetLastError()), __FILE__, __LINE__); \
        std::terminate(); \
    } \
} while (0)

#endif /*LM_CUDAHELPER_H_*/
