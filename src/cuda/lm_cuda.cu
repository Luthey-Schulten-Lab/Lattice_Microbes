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

#include <iostream>
#include <sstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "config.h"
#include "core/Exceptions.h"
#include "cuda/lm_cuda.h"

namespace lm
{

/*
 * Static members.
 */
CUresult CUDA::_cuda_dret_ = CUDA_SUCCESS;

int CUDA::getNumberDevices()
{
    // Get the number of CUDA devices found.
    int numberDevices;
    CUDA_EXCEPTION_CHECK(cudaGetDeviceCount(&numberDevices));
	
	return numberDevices;
}

void CUDA::setCurrentDevice(int device)
{
	if (device < 0 || device >= getNumberDevices())
	{
		throw InvalidArgException("The specified device was not valid.");
	}
	else
	{
		// Set that we are active on the CUDA device.
		CUDA_EXCEPTION_CHECK(cudaSetDevice(device));
	}
}

int CUDA::getCurrentDevice()
{
	int device;
	CUDA_EXCEPTION_CHECK(cudaGetDevice(&device));
	return device;
}

size_t CUDA::getFreeMemory(int device)
{
	// Make sure it is a valid device.
	if (device < 0 || device >= getNumberDevices())
	{
		throw InvalidArgException("The specified device was not valid.");
	}

    size_t free=0, total=0;
	cudaSetDevice(device);
    CUDA_EXCEPTION_CHECK(cudaMemGetInfo(&free, &total));
    return (size_t)free;
}

void CUDA::getComputeCapabilities(int device, int * major, int * minor)
{
    // Make sure it is a valid device.
    if (device < 0 || device >= getNumberDevices())
        throw InvalidArgException("The specified device was not valid.");
    if (major == NULL || minor == NULL)
        throw InvalidArgException("The specified return pointers were not valid.");

    // Get the device properties.
    cudaDeviceProp p;
    CUDA_EXCEPTION_CHECK(cudaGetDeviceProperties(&p, device));
    *major = p.major;
    *minor = p.minor;
}

void CUDA::printCapabilities(int device)
{
	// Make sure it is a valid device.
	if (device < 0 || device >= getNumberDevices())
		throw InvalidArgException("The specified device was not valid.");
	
    // Get the device properties.
    cudaDeviceProp p;
    CUDA_EXCEPTION_CHECK(cudaGetDeviceProperties(&p, device));
    
    // Print the device capabilities.
    std::cout << "CUDA Device(" << device << "): " << p.name << " v" << p.major << "." << p.minor << std::endl;
    std::cout << "  Device memory:             " << p.totalGlobalMem << " bytes" << std::endl;
    std::cout << "  Free memory:               " << CUDA::getFreeMemory(device) << " bytes" << std::endl;
    std::cout << "  Constant memory:           " << p.totalConstMem << " bytes" << std::endl;
    std::cout << "  Shared memory per block:   " << p.sharedMemPerBlock << " bytes" << std::endl;
    std::cout << "  Registers per block:       " << p.regsPerBlock << std::endl;
    std::cout << "  Warp size:                 " << p.warpSize << std::endl;
    std::cout << "  Maximum threads per block: " << p.maxThreadsPerBlock << std::endl;
    std::cout << "  Maximum block size:        " << p.maxThreadsDim[0] << "x" << p.maxThreadsDim[1] << "x" << p.maxThreadsDim[2] << std::endl;
    std::cout << "  Maximum grid size:         " << p.maxGridSize[0] << "x" << p.maxGridSize[1] << "x" << p.maxGridSize[2] << std::endl;
    std::cout << "  Maximum memory pitch:      " << p.memPitch << " bytes" << std::endl;
    std::cout << "  Texture alignment:         " << p.textureAlignment << " bytes" << std::endl;
    std::cout << "  Clock rate:                " << p.clockRate << " kilohertz" << std::endl;
    std::cout << std::endl;
}

std::string CUDA::getCapabilitiesString(int device)
{
    // Make sure it is a valid device.
    if (device < 0 || device >= getNumberDevices())
        throw InvalidArgException("The specified device was not valid.");

    // Get the device properties.
    cudaDeviceProp p;
    CUDA_EXCEPTION_CHECK(cudaGetDeviceProperties(&p, device));

    // Get a string version of the properties.
    std::ostringstream out (std::ostringstream::out);
    out << device << ":n=\"" << p.name \
        << "\",v=" << p.major << "." << p.minor \
        << ",mt=" <<  p.totalGlobalMem \
        << ",mf=" << CUDA::getFreeMemory(device) \
        << ",mc=" << p.totalConstMem \
        << ",ms=" << p.sharedMemPerBlock \
        << ",r=" << p.regsPerBlock \
        << ",w=" << p.warpSize \
        << ",c=" << p.clockRate;
    return out.str();
}

}
