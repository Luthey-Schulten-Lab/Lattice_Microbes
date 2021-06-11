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

#ifndef __ZDIVMMPIGPUMAPPER__
#define __ZDIVMMPIGPUMAPPER__

#include "SegmentDescriptor.h"

#define NUM_SHADOWS 2
#define Z_SHADOW_TOP 0
#define Z_SHADOW_BOTTOM 1

#define MCLKR_MAX_NBUFFERS 2


struct neighbor_buffer
{
	char *buffer[MCLKR_MAX_NBUFFERS];
	size_t size;
};

class ZDivMPIGPUMapper 
{
	private:
		bool periodic_z;

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
	neighbor_buffer *comm_buffer[NUM_SHADOWS];

	int rank, world_size, host_gpu;
//mgpumapper imports
	        dim3 lattice_dim;
        int overlap, apron;
        size_t cellsize;
      int pagecount;

	size_t local_size, global_size;
	size_t authority_size, authority_offset;
	ssize_t global_input_offset;
	size_t global_output_offset;

	
	protected:
		void stage_in_real(void *dptr, void *hptr, unsigned int element_size);
		void stage_out_real(void *hptr, void *dptr);

	public:
		ZDivMPIGPUMapper(int x, int y, int z, size_t, int, int, int, bool pz=false, int pages=1);
		~ZDivMPIGPUMapper();
	
		virtual SegmentDescriptor_s* initialize();
		virtual void initialize_gpu();

		dim3 get_global_dim();
		dim3 get_local_dim();
		int3 get_global_offset();
		size_t get_local_size();
		size_t get_authority_size();
		ssize_t get_global_input_offset();
		size_t get_global_output_offset();
		size_t get_authority_offset();

		void stage_in(void *dptr, void *hptr);
		void stage_in_sites(void *dptr, void *hptr);
		void stage_out(void *hptr, void *dptr);

/*
		void publish_state(int gpu, void *dptr, int key);
		void refresh(int gpu, void *dptr, int key);
		virtual void schedule_send(int gpu, void *dptr, int timestamp, int neighbor, cudaStream_t stream);
		virtual void schedule_recv(int gpu, void *dptr, int timestamp, int neighbor, cudaStream_t stream);

		int map_index_to_gpu(size_t index);
*/

	void communicate_edges(void *dptr, int key);
	void recv_edge(int key, int neigh);
	void send_edge(int key, int neigh);
	void copy_edge_to_device(void *dptr, int key, int neighbor, cudaStream_t s);
	void copy_edge_to_host(void *dptr, int key, int neighbor, cudaStream_t s);
	void stage_in_from_slave(int rank, void *hptr);
	void stage_out_to_master(void *hptr);

};

#endif
