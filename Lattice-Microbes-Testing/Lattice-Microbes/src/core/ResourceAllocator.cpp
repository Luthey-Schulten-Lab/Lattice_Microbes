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

#include <iostream>
#include <pthread.h>
#include <sstream>
#include "config.h"
#include "core/Exceptions.h"
#include "core/Math.h"
#include "core/Types.h"
#include "core/ResourceAllocator.h"
#include "thread/Thread.h"

using lm::thread::PthreadException;

namespace lm {
namespace main {

ResourceAllocator::ResourceAllocator(int processNumber, int numberCpuCores, float cpuCoresPerReplicate)
:processNumber(processNumber),numberCpuCores(numberCpuCores),reservedCpuCores(0),cpuSlots(NULL),cudaSlots(NULL)
{
    initialize(cpuCoresPerReplicate, 0.0);
}

ResourceAllocator::ResourceAllocator(int processNumber, int numberCpuCores, float cpuCoresPerReplicate, vector<int> cudaDevices, float cudaDevicesPerReplicate)
:processNumber(processNumber),numberCpuCores(numberCpuCores),reservedCpuCores(0),cudaDevices(cudaDevices),cpuSlots(NULL),cudaSlots(NULL)
{
    initialize(cpuCoresPerReplicate, cudaDevicesPerReplicate);
}

void ResourceAllocator::initialize(float cpuCoresPerReplicate, float cudaDevicesPerReplicate)
{
    // Create the mutex.
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_init(&mutex, NULL));

    // Figure out how many cpu slots we have and how many we need per replicate.
    if (cpuCoresPerReplicate == 0.0)
    {
        cpuSlotsPerCore = 0;
        cpuSlotsPerReplicate = 0;
    }
    else if (lroundf(floorf(cpuCoresPerReplicate)) < 1)
    {
        cpuSlotsPerCore = lroundf(1.0f/cpuCoresPerReplicate);
        cpuSlotsPerReplicate = 1;
    }
    else
    {
        cpuSlotsPerCore = 1;
        cpuSlotsPerReplicate = lroundf(cpuCoresPerReplicate);
    }

    // Figure out how many cuda slots we have and how many we need per replicate.
    if (cudaDevicesPerReplicate == 0.0)
    {
        cudaSlotsPerDevice = 0;
        cudaSlotsPerReplicate = 0;
    }
    else if (lroundf(floorf(cudaDevicesPerReplicate)) < 1)
    {
        cudaSlotsPerDevice = lroundf(1.0f/cudaDevicesPerReplicate);
        cudaSlotsPerReplicate = 1;
    }
    else
    {
        cudaSlotsPerDevice = 1;
        cudaSlotsPerReplicate = lroundf(cudaDevicesPerReplicate);
    }

    // Create the thread assignment table.

    // Create the cpu assignment table.
    cpuSlots = new int*[numberCpuCores];
    for (int i=0; i<numberCpuCores; i++)
    {
        cpuSlots[i] = new int[cpuSlotsPerCore];
        for (int j=0; j<cpuSlotsPerCore; j++)
            cpuSlots[i][j] = -1;
    }

    // Create the cuda assignment table.
    cudaSlots = new int*[cudaDevices.size()];
    for (int i=0; i<(int)cudaDevices.size(); i++)
    {
        cudaSlots[i] = new int[cudaSlotsPerDevice];
        for (int j=0; j<cudaSlotsPerDevice; j++)
            cudaSlots[i][j] = -1;
    }
}

ResourceAllocator::~ResourceAllocator()
{
    // Deallocate the cpu assignment table.
    if (cpuSlots != NULL)
    {
        for (int i=0; i<numberCpuCores; i++)
            delete[] cpuSlots[i];
        delete[] cpuSlots;
    }
    cpuSlots = NULL;

    // Deallocate the cuda assignment table.
    if (cudaSlots != NULL)
    {
        for (int i=0; i<(int)cudaDevices.size(); i++)
            delete[] cudaSlots[i];
        delete[] cudaSlots;
    }
    cudaSlots = NULL;

    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_mutex_destroy(&mutex));
}

int ResourceAllocator::getMaxSimultaneousReplicates()
{
    // Figure out how many replicates we can run simultaneously.
    if (cpuSlotsPerReplicate > 0 && cudaSlotsPerReplicate == 0)
        return cpuSlotsPerCore*(numberCpuCores-reservedCpuCores)/cpuSlotsPerReplicate;
    else if (cpuSlotsPerReplicate == 0 && cudaSlotsPerReplicate > 0)
        return cudaSlotsPerDevice*cudaDevices.size()/cudaSlotsPerReplicate;
    else if (cpuSlotsPerReplicate > 0 && cudaSlotsPerReplicate > 0)
    	return min((uint)(cpuSlotsPerCore*(numberCpuCores-reservedCpuCores)/cpuSlotsPerReplicate), (uint)(cudaSlotsPerDevice*cudaDevices.size()/cudaSlotsPerReplicate));

    return 0;
}

ResourceAllocator::ComputeResources ResourceAllocator::assignReplicate(int replicate)
{
    ComputeResources allocatedResources;
    allocatedResources.processNumber = processNumber;

    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));

    // Walk through the cpu table and find free resources.
    for (int j=0; j<cpuSlotsPerCore && (int)allocatedResources.cpuCores.size()<cpuSlotsPerReplicate; j++)
    {
        for (int i=0; i<numberCpuCores && (int)allocatedResources.cpuCores.size()<cpuSlotsPerReplicate; i++)
        {
            if (cpuSlots[i][j] == -1)
            {
                allocatedResources.cpuCores.push_back(i);
                cpuSlots[i][j] = replicate;
            }
        }
    }

    // Walk through the gpu table and find free resources.
    for (int j=0; j<cudaSlotsPerDevice && (int)allocatedResources.cudaDevices.size()<cudaSlotsPerReplicate; j++)
    {
        for (int i=0; i<(int)cudaDevices.size() && (int)allocatedResources.cudaDevices.size()<cudaSlotsPerReplicate; i++)
        {
            if (cudaSlots[i][j] == -1)
            {
                allocatedResources.cudaDevices.push_back(cudaDevices[i]);
                cudaSlots[i][j] = replicate;
            }
        }
    }

    // Make sure we found enough resources.
    if ((int)allocatedResources.cpuCores.size() != cpuSlotsPerReplicate || (int)allocatedResources.cudaDevices.size() != cudaSlotsPerReplicate)
    {
        PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
        removeReplicate(replicate);
        throw lm::Exception("Could not find enough resources to assign the replicate.");
    }

    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex

    return allocatedResources;
}

int ResourceAllocator::reserveCpuCore()
{
    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));

    // Walk through the CPU table and find a core that has no slots in use.
    int reservedReplicate=-1;
    for (int i=0; i<numberCpuCores; i++)
    {
    	bool allSlotsFree=true;
        for (int j=0; j<cpuSlotsPerCore; j++)
        {
        	if (cpuSlots[i][j] != -1)
        	{
        		allSlotsFree = false;
        		break;
        	}
        }

        // If this core had no slots in use, reserve it.
        if (allSlotsFree)
        {
        	reservedCpuCores++;
        	reservedReplicate = i;
        	for (int j=0; j<cpuSlotsPerCore; j++)
        		cpuSlots[i][j] = -2;
        	break;
        }
    }

    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex

    // Return the replicate or throw an exception if we couldn't reserve one.
    if (reservedReplicate >= 0)
    	return reservedReplicate;
    else
    	throw lm::Exception("No free cpu cores available to reserve.");
}

void ResourceAllocator::removeReplicate(int replicate)
{
    //// BEGIN CRITICAL SECTION: mutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&mutex));

    // Walk through the cpu table and remove any resources assigned to the replicate.
    for (int i=0; i<numberCpuCores; i++)
    {
        for (int j=0; j<cpuSlotsPerCore; j++)
        {
            if (cpuSlots[i][j] == replicate)
            {
                cpuSlots[i][j] = -1;
            }
        }
    }

    // Walk through the gpu table and remove any resources assigned to the replicate.
    for (int i=0; i<(int)cudaDevices.size(); i++)
    {
        for (int j=0; j<cudaSlotsPerDevice; j++)
        {
            if (cudaSlots[i][j] == replicate)
            {
                cudaSlots[i][j] = -1;
            }
        }
    }

    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&mutex));
    //// END CRITICAL SECTION: mutex
}

string ResourceAllocator::ComputeResources::toString()
{
    std::ostringstream out(std::ostringstream::out);
    for (int i=0; i<(int)cpuCores.size(); i++)
    {
        if (i==0) out << "cpu_id=";
        else out << ",";
        out << cpuCores[i];
    }
    for (int i=0; i<(int)cudaDevices.size(); i++)
    {
        if (i==0) out << "/cuda_id=";
        else out << ",";
        out << cudaDevices[i];
    }
    return out.str();
}

}
}
