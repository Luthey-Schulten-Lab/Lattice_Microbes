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

#ifndef __MULTIGPUMAPPER__
#define __MULTIGPUMAPPER__

#define POSDIM(_p, _d) ((_p).x + (_p).y * (_d).x + (_p).z * (_d).x * (_d).y)
#define DIMSIZE(_d) ((_d).x*(_d).y*(_d).z)

#include <pthread.h>
#include "SegmentDescriptor.h"

class MultiGPUMapper
{
	protected:
	int num_gpus;
	int *device_id;
	dim3 lattice_dim;
	int overlap;
	size_t cellsize;
	pthread_key_t affinity;
	int apron;
	SegmentDescriptor_s **descriptor;
	float *lb_weights;
	int *lb_cost;
	size_t *device_memory;
	int pagecount;

	virtual void initialize()=0;
	bool enable_peer_access(int src, int dst);
	void build_descriptor(int gpu, dim3 ldim, int3 goffset, dim3 active, dim3 loffset);
	void compute_balances();

	public:

	MultiGPUMapper(dim3 ldim, size_t cellsize, int apron, int overlap, int num_gpus, int* devices, int pages);
	virtual ~MultiGPUMapper();

	int get_num_gpus();
	bool use(int gpu);

	int get_overlap();

	int get_apron();

	void set_affinity(int);
	int get_affinity();

	dim3 get_lattice_dim();
	SegmentDescriptor_s* getSegmentDescriptor(int gpu);
	size_t get_global_size();

	void record_execution_cost(int, int);
	bool rebalance();
	bool numa_bind_thread(int);

	// virtual interface
	virtual dim3 get_global_dim(int gpu)=0;
	virtual dim3 get_local_dim(int gpu)=0;
	virtual int3 get_global_offset(int gpu)=0;
	virtual size_t get_local_size(int gpu)=0;
	virtual size_t get_authority_size(int gpu)=0;
	virtual ssize_t get_global_input_offset(int gpu)=0;
	virtual size_t get_global_output_offset(int gpu)=0;
	virtual size_t get_authority_offset(int gpu)=0;
	virtual void stage_in(int gpu, void *dptr, void *hptr)=0;
	virtual void stage_in_sites(int gpu, void *dptr, void *hptr)=0;
	virtual void stage_out(int gpu, void *hptr, void *dptr)=0;
	virtual void publish_state(int gpu, void *dptr, int timestamp)=0;
	virtual void refresh(int gpu, void *dptr, int timestamp)=0;

	virtual void schedule_send(int gpu, void *dptr, int timestamp, int neighbor, cudaStream_t stream)=0;
	virtual void schedule_recv(int gpu, void *dptr, int timestamp, int neighbor, cudaStream_t stream)=0;

	virtual int map_index_to_gpu(size_t index)=0;
	
	virtual void initialize_gpu(int gpu);

	virtual bool determine_load_balance()=0;

};

#define check_error() ({ cudaError_t err=cudaGetLastError(); if(err!=cudaSuccess) { printf("Cuda error %s:%d: %s\n", __FILE__, __LINE__,cudaGetErrorString(err)); exit(1); } })


#endif 
