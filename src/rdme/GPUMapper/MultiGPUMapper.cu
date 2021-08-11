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
 * Author(s): Mike Hallock
 */

#include <cuda.h>
#include <cuda_runtime.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "config.h"

#include "SegmentDescriptor.h"
#include "MultiGPUMapper.h"

#include <pthread.h>

#ifdef MPD_NUMA_SUPPORT
#include <numa.h>
#endif

#include "cuda/lm_cuda.h"

#include "core/Print.h"
using lm::Print;

#include <vector>
#include <string>
extern std::vector<int> cudaDevices;
extern std::vector<int> numaNodes;
extern bool mgpu_disablePeering;

MultiGPUMapper::MultiGPUMapper(dim3 ldim, size_t cellsize, int apron=1, int overlap=0, int ngpus=0, int* devices=NULL, int pages=1)
:lattice_dim(ldim), cellsize(cellsize), apron(apron), overlap(overlap), device_id(NULL), pagecount(pages)
{
	if(! ngpus)
	{
		cudaGetDeviceCount(&num_gpus);
	}
	else
	{
		num_gpus=ngpus;
	}

	device_id=new int[num_gpus];
	pthread_key_create(&affinity,NULL);

	for(int g=0; g<num_gpus; g++)
	{
		Print::printf(Print::DEBUG, "[mgpu] device list %d: %d", g, devices[g]);
        device_id[g]=devices[g];
	}
	
	descriptor=new SegmentDescriptor_s*[num_gpus];
	lb_weights=new float[num_gpus];
	lb_cost=new int[num_gpus];
	device_memory=new size_t[num_gpus];
	for(int i=0; i<num_gpus; i++)
	{
		//lm::CUDA::printCapabilities(device_id[i]);
		lb_weights[i]=1;
		lb_cost[i]=100;
		descriptor[i]=new SegmentDescriptor_s;
		device_memory[i]=lm::CUDA::getFreeMemory(device_id[i]);
		//Print::printf(Print::DEBUG, "[mgpu] Device %d free mem %llu\n", i, device_memory[i]);
	}

	// Uncomment for testing to artifically constrain memory
	// device_memory[0]=1024*1024*6;
}

MultiGPUMapper::~MultiGPUMapper()
{
	if(device_id)
		delete device_id;

	pthread_key_delete(affinity);
}

int MultiGPUMapper::get_num_gpus()
{
	return num_gpus;
}

bool MultiGPUMapper::use(int gpu)
{
	if(gpu < 0 || gpu >= num_gpus)
		return false;

	cudaError_t err=cudaSetDevice(device_id[gpu]);

	return (err == cudaSuccess);
}

int MultiGPUMapper::get_overlap()
{
	return overlap;
}

int MultiGPUMapper::get_apron()
{
	return apron;
}

void MultiGPUMapper::set_affinity(int gpu)
{
	pthread_setspecific(affinity, reinterpret_cast<void *>(gpu));
}

int MultiGPUMapper::get_affinity()
{
	return (size_t)pthread_getspecific(affinity);
}

dim3 MultiGPUMapper::get_lattice_dim()
{
	return lattice_dim;
}

SegmentDescriptor_s* MultiGPUMapper::getSegmentDescriptor(int gpu)
{
	if(gpu < 0 || gpu >= num_gpus)
		return NULL;

	return descriptor[gpu];
}

void MultiGPUMapper::build_descriptor(int gpu, dim3 ldim, int3 goffset, dim3 active, dim3 loffset)
{
	SegmentDescriptor_s *seg=descriptor[gpu];
	seg->local_dimensions=ldim;
	seg->global_offset=goffset;
	seg->active_dimensions=active;
	seg->active_offset=loffset;

	Print::printf(Print::DEBUG, "* Descriptor %d: "
                           " local dim:     %d x %d x %d\n"
           "*                active dim:    %d x %d x %d\n"
           "*                global offset: %d x %d x %d\n"
           "*                active offset: %d x %d x %d\n",gpu,
			ldim.x, ldim.y, ldim.z,
			active.x, active.y, active.z,
			goffset.x, goffset.y, goffset.z,
			loffset.x, loffset.y, loffset.z
	);
}

bool MultiGPUMapper::enable_peer_access(int src, int dst)
{
	if(dst < 0 || dst >= num_gpus)
		return false;

	if(mgpu_disablePeering)
		return false;

	bool is_peered=false;
	int local=device_id[src];
	int peer=device_id[dst];
	if(use(src))
	{
		std::string msg;
		int can_access=0;
		cudaDeviceCanAccessPeer(&can_access, local, peer);
		if(can_access)
		{
			switch(cudaDeviceEnablePeerAccess(peer,0))
			{
				case cudaSuccess:
				case cudaErrorPeerAccessAlreadyEnabled:
				
					msg="Peer access enabled";
					cudaGetLastError(); // clear out potential already-enabled error
					is_peered=true;
					break;

				default:
					msg="Peer access setup FAILED";
			}
		}
		else
		{
			msg="NOTICE: peer access not available";
		}

		Print::printf(Print::DEBUG, "%s from device %d to %d (logical %d->%d)", msg.c_str(), local, peer, src, dst);
	}

	return is_peered;
}

bool MultiGPUMapper::numa_bind_thread(int gpu)
{
#ifdef MPD_NUMA_SUPPORT
	if(gpu >= numaNodes.size())
		return false;
#ifdef DEBUG
	Print::printf(Print::DEBUG_VERBOSE, "Binding gpu thread %d to node %d", gpu, numaNodes[gpu]);
#endif

	nodemask_t nm;
	nodemask_zero(&nm);
	nodemask_set(&nm, numaNodes[gpu]);
	numa_bind(&nm);
	return true;
#else
	Print::printf(Print::DEBUG, "Built without -DMPD_NUMA_SUPPORT, cannot bind thread");
	return false;
#endif
}

void MultiGPUMapper::record_execution_cost(int gpu, int etime)
{
	lb_cost[gpu]=etime;	
}

bool MultiGPUMapper::rebalance()
{
	if(! determine_load_balance())
		return false;

	printf("\n*** Rebalance *** \n");
	initialize();
	printf("***************** \n\n");

	return true;
}

void MultiGPUMapper::initialize_gpu(int gpu)
{
	if(! use(gpu))
		throw "Failed to use gpu!";

}

size_t MultiGPUMapper::get_global_size()
{
	return DIMSIZE(lattice_dim) * cellsize;
}
