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

#ifndef __ZDIVMULTIGPUMAPPER__
#define __ZDIVMULTIGPUMAPPER__

#include "MultiGPUMapper.h"

#define NUM_SHADOWS 2
#define Z_SHADOW_TOP 0
#define Z_SHADOW_BOTTOM 1

#define MCLKR_MAX_NBUFFERS 2

struct neighbor_buffer
{
	char *buffer[MCLKR_MAX_NBUFFERS];
	size_t size;
};

struct gpu_info
{
	dim3 global_pos;		// Position in the global lattice that this gpu is associated with
	dim3 global_dim;		// Dimensions in the global lattice that this is responsible for

	dim3 local_dim;			// Local GPU dimensions, including all padding
	dim3 local_authority;	// Position in local GPU that begins authoritative data

	int3 dependent_pos;		// Position in global lattice where we can get local_dim data
							// NOTE: Can be negative for periodic conditions

	dim3 overlap_dim;
	dim3 overlap_send[NUM_SHADOWS];
	dim3 overlap_recv[NUM_SHADOWS];
	int neighbor[NUM_SHADOWS];
	neighbor_buffer *read_buffer[NUM_SHADOWS];
	neighbor_buffer *write_buffer[NUM_SHADOWS];

	int lb_chunks;			// Number of lattice chunks assigned for load balancing
	int lb_max_chunks;
	float lb_imbalance;

	char *tmp_buffer[MCLKR_MAX_NBUFFERS];
};

	
class ZDivMultiGPUMapper : public MultiGPUMapper
{
	private:
		gpu_info *info;
		bool periodic_z;
	
	protected:
		virtual void initialize();
		//void copy_from_neighbor(int gpu, int neighbor, void *dptr, int key);
		//void copy_to_neighbor(int gpu, int neighbor, int key);
		void stage_in_real(int gpu, void *dptr, void *hptr, unsigned int element_size);
		void stage_out_real(int gpu, void *hptr, void *dptr);

		bool determine_load_balance();
		int alter_chunks(int gpu, int count);

	public:
		ZDivMultiGPUMapper(dim3, size_t, int, int, int, int* gl=NULL, bool pz=false, int pages=1);
		ZDivMultiGPUMapper(int x, int y, int z, size_t, int, int, int, int* gl=NULL, bool pz=false, int pages=1);
		~ZDivMultiGPUMapper();
	
		virtual void initialize_gpu(int gpu);
		dim3 get_global_dim(int gpu);
		dim3 get_local_dim(int gpu);
		int3 get_global_offset(int gpu);
		size_t get_local_size(int gpu);
		size_t get_authority_size(int gpu);
		ssize_t get_global_input_offset(int gpu);
		size_t get_global_output_offset(int gpu);
		size_t get_authority_offset(int gpu);

		void stage_in(int gpu, void *dptr, void *hptr);
		void stage_in_sites(int gpu, void *dptr, void *hptr);
		void stage_out(int gpu, void *hptr, void *dptr);

		void publish_state(int gpu, int key, cudaStream_t top, cudaStream_t bot, void *dptr = NULL);
		void refresh(int gpu, void *dptr, int key);
		virtual void schedule_send(int gpu, void *dptr, int timestamp, int neighbor, cudaStream_t stream);
		virtual void schedule_recv(int gpu, void *dptr, int timestamp, int neighbor, cudaStream_t stream);

		int map_index_to_gpu(size_t index);

		gpu_info* getinfo(int gpu);
		unsigned int* gettbuf(int gpu, int key, int neighbor);
		unsigned int* getrbuf(int gpu, int key, int neighbor);
};

#endif
