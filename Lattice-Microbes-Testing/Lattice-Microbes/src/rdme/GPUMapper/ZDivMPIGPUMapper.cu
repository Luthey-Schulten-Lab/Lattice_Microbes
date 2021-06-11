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
#include <assert.h>

#include "config.h"
#include "lptf/Profile.h"

#include "MultiGPUMapper.h"
#include "ZDivMPIGPUMapper.h"

#include "core/Print.h"
using lm::Print;

#include <mpi.h>

ZDivMPIGPUMapper::ZDivMPIGPUMapper(int x, int y, int z, size_t _cellsize, int _apron, int _overlap, int gpu, bool pz, int pages) 
{
	if(MPI_Comm_size(MPI_COMM_WORLD, &world_size)!=MPI_SUCCESS)
		exit(1);
	if(MPI_Comm_rank(MPI_COMM_WORLD, &rank)!=MPI_SUCCESS)
		exit(1);

	apron = _apron;
	overlap = _overlap;
	periodic_z = pz;
	assert(apron == 2);
	assert(pz == true);
	assert(overlap == 0);

	lattice_dim=dim3(x,y,z);
	pagecount=pages;
	cellsize=_cellsize;
	host_gpu=gpu;
}

ZDivMPIGPUMapper::~ZDivMPIGPUMapper()
{
}

SegmentDescriptor_s *ZDivMPIGPUMapper::initialize()
{

/*{
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
        sleep(5);
}*/

	// Split up lattice as equally as possible
	// must be in "chunks" of the z block size
	int total_chunks = lattice_dim.z / TUNE_MPD_Z_BLOCK_Z_SIZE;
	int slice = total_chunks / world_size;
	//printf("rank %d total %d slice %d\n", rank, total_chunks, slice);

	// determine the number of chunks for this rank
	int csize,zdim;
	if(rank == world_size-1)
		csize=total_chunks - (slice*rank);
	else
		csize=slice;
	zdim = csize*TUNE_MPD_Z_BLOCK_Z_SIZE;

	// starting z index for this rank
	int zstart = slice*TUNE_MPD_Z_BLOCK_Z_SIZE*rank;

	//printf("rank %d csize %d zdim %d zstart %d\n", rank, csize, zdim, zstart);

	int top_overlap = ((periodic_z || rank > 0) ? overlap : 0);
	int bottom_overlap = ((periodic_z || rank < world_size-1) ? overlap : 0);

	int top_apron = ((periodic_z || rank > 0) ? apron : 0);
	int bottom_apron = ((periodic_z || rank < world_size-1) ? apron : 0);

	int top_shadow = top_overlap+top_apron;
	int bottom_shadow = bottom_overlap+bottom_apron;

	// Store global position and dimension that this rank processes
	global_pos=dim3(0,0, zstart);
	global_dim=dim3(lattice_dim.x, lattice_dim.y, zdim);

	// Store the dimensions of the local lattice
	local_dim=dim3(lattice_dim.x, lattice_dim.y, zdim+top_shadow+bottom_shadow);
	
	// Local authority is where in the local lattice authoritative data begins
	local_authority=dim3(0,0, top_shadow);

	// Dependent position is where the global authority, less shadow region
	// NOTE: can be negative in periodic cases!
	dependent_pos.x=0;
	dependent_pos.y=0;
	dependent_pos.z=zstart - top_shadow;

	// dimensions of the shadow regions
	overlap_dim=dim3(lattice_dim.x, lattice_dim.y, overlap+apron);

	// Where my neighbor reads from me
	overlap_send[Z_SHADOW_TOP]=dim3(0,0,overlap+apron);
	overlap_send[Z_SHADOW_BOTTOM]=dim3(0,0, zdim-(overlap+apron) + top_shadow);

	// Where I receive
	overlap_recv[Z_SHADOW_TOP]=dim3(0,0,0);
	overlap_recv[Z_SHADOW_BOTTOM]=dim3(0,0, zdim + top_shadow);

	neighbor[Z_SHADOW_TOP]=rank-1;
	neighbor[Z_SHADOW_BOTTOM]=rank+1;

	// If Z-Periodic, set top and bottom to be neighbors
	if(periodic_z)
	{
		if(rank==0)
		{
			neighbor[Z_SHADOW_TOP]=world_size-1;
		}
		else if(rank==world_size-1)
		{
			neighbor[Z_SHADOW_BOTTOM]=0;
		}
	}

	SegmentDescriptor_s *seg=new SegmentDescriptor_s;
	seg->local_dimensions=local_dim;
	seg->global_offset=dependent_pos;
	seg->active_dimensions=dim3(lattice_dim.x, lattice_dim.y, (zdim+top_overlap+bottom_overlap));
	seg->active_offset=dim3(0, 0, top_apron);

	Print::printf(Print::DEBUG, "* Descriptor %d: "
                           " local dim:     %d x %d x %d\n"
           "*                active dim:    %d x %d x %d\n"
           "*                global offset: %d x %d x %d\n"
           "*                active offset: %d x %d x %d\n",
			rank,
			local_dim.x, local_dim.y, local_dim.z,
			seg->active_dimensions.x, seg->active_dimensions.y, seg->active_dimensions.z,
			seg->global_offset.x, seg->global_offset.y, seg->global_offset.z,
			seg->active_offset.x, seg->active_offset.y, seg->active_offset.z
	);

	// Other stored things
	local_size = DIMSIZE(local_dim)*cellsize;
	global_size = DIMSIZE(lattice_dim) * cellsize;
	authority_size = DIMSIZE(global_dim)*cellsize;
	authority_offset = POSDIM(local_authority, local_dim);
	global_input_offset = dependent_pos.x + dependent_pos.y * (int)global_dim.x + dependent_pos.z * (int)global_dim.x * (int)global_dim.y;
	global_output_offset=POSDIM(global_pos, global_dim);


	return seg;
}

void ZDivMPIGPUMapper::initialize_gpu()
{
	cudaError_t err=cudaSetDevice(host_gpu);
	if (err != cudaSuccess)
		throw("Failed to use device");

	comm_buffer[Z_SHADOW_TOP]=NULL;
	comm_buffer[Z_SHADOW_BOTTOM]=NULL;

	// Allocate memory for top neighbor
	size_t ss=DIMSIZE(overlap_dim)*cellsize;

	int neigh = neighbor[Z_SHADOW_TOP];
	if(neigh >= 0 && neigh < world_size)	
	{
		neighbor_buffer *nb=new neighbor_buffer;
		nb->size=ss;

			cudaHostAlloc(&(nb->buffer[0]), ss*pagecount, cudaHostAllocPortable);
			cudaHostAlloc(&(nb->buffer[1]), ss*pagecount, cudaHostAllocPortable);

		comm_buffer[Z_SHADOW_TOP]=nb;
	}

	// Peer with bottom neighbor
	neigh = neighbor[Z_SHADOW_BOTTOM];
	if(neigh >= 0 && neigh < world_size)	
	{
		neighbor_buffer *nb=new neighbor_buffer;
		nb->size=ss;

		cudaHostAlloc(&(nb->buffer[0]), ss*pagecount, cudaHostAllocPortable);
		cudaHostAlloc(&(nb->buffer[1]), ss*pagecount, cudaHostAllocPortable);

		comm_buffer[Z_SHADOW_BOTTOM]=nb;
	}

	check_error();

	// Attach MPI Buffer

	size_t buffer_size = (ss*2+MPI_BSEND_OVERHEAD) * 10;
	void *buffer=malloc(buffer_size);
	MPI_Buffer_attach(buffer, buffer_size);
}


dim3 ZDivMPIGPUMapper::get_global_dim()
{
	return global_dim;
}

dim3 ZDivMPIGPUMapper::get_local_dim()
{
	return local_dim;
}

/*
int3 ZDivMPIGPUMapper::get_global_offset(int gpu)
{
	return info[gpu].dependent_pos;
}
*/

size_t ZDivMPIGPUMapper::get_local_size()
{
	return local_size;
}

/*
size_t ZDivMPIGPUMapper::get_authority_size(int gpu)
{
	return DIMSIZE(info[gpu].global_dim)*cellsize;
}

ssize_t ZDivMPIGPUMapper::get_global_input_offset(int gpu)
{
	const int3 dp=info[gpu].dependent_pos;
	const dim3 gd=info[gpu].global_dim;

	return dp.x + dp.y * (int)gd.x + dp.z * (int)gd.x * (int)gd.y;
}

size_t ZDivMPIGPUMapper::get_global_output_offset(int gpu)
{
	return POSDIM(info[gpu].global_pos, info[gpu].global_dim);
}

size_t ZDivMPIGPUMapper::get_authority_offset(int gpu)
{
	return POSDIM(info[gpu].local_authority, info[gpu].local_dim);
}

*/

void ZDivMPIGPUMapper::stage_in(void *dptr, void *hptr)
{
	// Strided on the device by local_size
	// Strided on the host by global_size
	for(int p=0; p<pagecount; p++)
	{
		char *dpage = (char*)dptr + p*local_size;
		char *hpage = (char*)hptr + p*global_size;
		stage_in_real(dpage, hpage, cellsize);
	}
}

void ZDivMPIGPUMapper::stage_in_sites(void *dptr, void *hptr)
{
	stage_in_real(dptr, hptr, 1);
}

// Copy from host memory to device memory
// TODO dont need edge data, will get that later
void ZDivMPIGPUMapper::stage_in_real(void *dptr, void *hptr, unsigned int element_size) 
{
	PROF_BEGIN(PROF_H2D);
	char *src=(char*)hptr;
	char *dst=(char*)dptr;
	ssize_t offset=global_input_offset*element_size;
	size_t localsize=DIMSIZE(local_dim)*element_size;

	size_t latsize=DIMSIZE(lattice_dim)*element_size;

	if(offset < 0)
	{
		// Read adj bytes from the end of the lattice
		size_t adj=abs(offset);
		cudaMemcpyAsync(dst, src+latsize-adj, adj, cudaMemcpyDefault);
		check_error();
	
		dst+=adj;
		offset=0;
		localsize-=adj;
	}
	else if((offset+localsize) >=latsize)
	{
		// read adj bytes from the beginning of the lattice
		size_t adj=(offset+localsize)-latsize;
		cudaMemcpyAsync(dst+localsize-adj, src, adj, cudaMemcpyDefault);
		check_error();
		localsize-=adj;
	}
		
	src += offset;
	cudaMemcpy(dst, src, localsize, cudaMemcpyDefault);
	check_error();
	PROF_END(PROF_H2D);
}

void ZDivMPIGPUMapper::stage_out(void *hptr, void *dptr) 
{
	for(int p=0; p<pagecount; p++)
	{
		char *dpage = (char*)dptr + p*local_size;
		char *hpage = (char*)hptr + p*global_size;
		stage_out_real(hpage, dpage);
	}
}

void ZDivMPIGPUMapper::stage_out_real(void *hptr, void *dptr) 
{
	PROF_BEGIN(PROF_D2H);
	char *src=(char*)dptr;
	src += authority_offset * cellsize;

	char *dst=(char*)hptr;
	dst += global_output_offset * cellsize;

	cudaMemcpy(dst, src, authority_size, cudaMemcpyDefault);
	check_error();
	PROF_END(PROF_D2H);
}

/*
void ZDivMPIGPUMapper::publish_state(int gpu, void *dptr, int timestep)
{
	schedule_send(dptr, timestep, Z_SHADOW_TOP, 0);
	schedule_send(dptr, timestep, Z_SHADOW_BOTTOM, 0);
	if(cudaDeviceSynchronize() != cudaSuccess)
		throw("cuda error");
}

void ZDivMPIGPUMapper::schedule_recv(int gpu, void *dptr, int key, int neighbor, cudaStream_t s)
{
	int other=info[gpu].neighbor[neighbor];
	if(other < 0 || other >= num_gpus)
		return;	

	key=key&1;

	neighbor_buffer *nb=info[gpu].read_buffer[neighbor];
	
	char *lbuf=(char*)dptr;
	size_t lofs=POSDIM(info[gpu].overlap_recv[neighbor], info[gpu].local_dim);
	lbuf+=lofs*cellsize;

	PROF_CUDA_BEGIN(PROF_MRBASE+neighbor, s);
	char *nbuf=nb->buffer[key];
	for(int p=0; p<pagecount; p++)
		cudaMemcpyAsync(lbuf+p*local_size,
			nbuf+p*(nb->size),
			nb->size, cudaMemcpyDefault,s);
	PROF_CUDA_END(PROF_MRBASE+neighbor, s);
}

void ZDivMPIGPUMapper::schedule_send(int gpu, void *dptr, int key, int neighbor, cudaStream_t s)
{
	int other=info[gpu].neighbor[neighbor];
	if(other < 0 || other >= num_gpus)
		return;	

	size_t ofs;
	char *buf;

	key=(key&1)^1;
	
	neighbor_buffer *nb=info[gpu].write_buffer[neighbor];

	ofs=POSDIM(info[gpu].overlap_send[neighbor], info[gpu].local_dim);
	buf=(char*)dptr+ofs*cellsize;

	PROF_CUDA_BEGIN(PROF_MSBASE+neighbor, s);
	char *nbuf=nb->buffer[key];
	for(int p=0; p<pagecount; p++)
		cudaMemcpyAsync(nbuf+p*(nb->size),
			buf+p*local_size,
			nb->size, cudaMemcpyDefault,s);
	PROF_CUDA_END(PROF_MSBASE+neighbor, s);
}


void ZDivMPIGPUMapper::refresh(int gpu, void *dptr, int timestep)
{
	schedule_recv(gpu, dptr, timestep, Z_SHADOW_TOP, 0);
	schedule_recv(gpu, dptr, timestep, Z_SHADOW_BOTTOM, 0);
	if(cudaDeviceSynchronize() != cudaSuccess)
		throw("cuda error");
}

*/

void ZDivMPIGPUMapper::copy_edge_to_host(void *dptr, int key, int neighbor, cudaStream_t s)
{
	neighbor_buffer *nb=comm_buffer[neighbor];
	if(nb == NULL)
		return;	

	size_t ofs;
	char *buf;

	key=(key&1)^1;
	//printf("rank %d edge to host neigh %d key %d\n", rank, neighbor, key);

	ofs=POSDIM(overlap_send[neighbor], local_dim);
	buf=(char*)dptr+ofs*cellsize;

	PROF_CUDA_BEGIN(PROF_MSBASE+neighbor, s);
	char *nbuf=nb->buffer[key];
	for(int p=0; p<pagecount; p++)
		cudaMemcpyAsync(nbuf+p*(nb->size),
			buf+p*local_size,
			nb->size, cudaMemcpyDefault,s);
	PROF_CUDA_END(PROF_MSBASE+neighbor, s);
}

void ZDivMPIGPUMapper::copy_edge_to_device(void *dptr, int key, int neighbor, cudaStream_t s)
{
	neighbor_buffer *nb=comm_buffer[neighbor];
	if(nb == NULL)
		return;	

	key=key&1;
	//printf("rank %d edge to dev neigh %d key %d\n", rank, neighbor, key);

	char *lbuf=(char*)dptr;
	size_t lofs=POSDIM(overlap_recv[neighbor], local_dim);
	lbuf+=lofs*cellsize;

	PROF_CUDA_BEGIN(PROF_MRBASE+neighbor, s);
	char *nbuf=nb->buffer[key];
	for(int p=0; p<pagecount; p++)
		cudaMemcpyAsync(lbuf+p*local_size,
			nbuf+p*(nb->size),
			nb->size, cudaMemcpyDefault,s);
	PROF_CUDA_END(PROF_MRBASE+neighbor, s);
}

void ZDivMPIGPUMapper::send_edge(int key, int neigh)
{
	neighbor_buffer *nb=comm_buffer[neigh];
	if(nb == NULL)
		return;	

	key=(key&1)^1;

	// printf("rank %d edge send neigh %d key %d\n", rank, neigh, key);

	int retval=MPI_Bsend(nb->buffer[key],
		nb->size*pagecount,
		MPI_BYTE,
		neighbor[neigh],
		1001, MPI_COMM_WORLD);		// TODO define this tag number somewhere
	if(retval != MPI_SUCCESS)
	{
		char error[MPI_MAX_ERROR_STRING];
		int errlen;
		MPI_Error_string(retval, error, &errlen);
		printf("MPI_Bsend failed on rank %d key %d neigh %d with error %d: %s\n", rank, key, neigh,retval, error);
		MPI_Abort(MPI_COMM_WORLD, retval);
	}
}

void ZDivMPIGPUMapper::recv_edge(int key, int neigh)
{
	neighbor_buffer *nb=comm_buffer[neigh];
	if(nb == NULL)
		return;	

	key=key&1;

	// printf("rank %d edge recv neigh %d key %d\n", rank, neigh, key);

    MPI_Recv(nb->buffer[key],
        nb->size*pagecount,
        MPI_BYTE,
        neighbor[neigh],
        1001, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void ZDivMPIGPUMapper::communicate_edges(void *dptr, int key)
{
	//printf("** comm edge rank %d\n", rank);
	copy_edge_to_host(dptr, key, Z_SHADOW_TOP, 0);
	copy_edge_to_host(dptr, key, Z_SHADOW_BOTTOM, 0);
	cudaDeviceSynchronize();

	send_edge(key, Z_SHADOW_TOP);
	send_edge(key, Z_SHADOW_BOTTOM);
	
	recv_edge(key, Z_SHADOW_BOTTOM);
	recv_edge(key, Z_SHADOW_TOP);

	copy_edge_to_device(dptr, key, Z_SHADOW_TOP, 0);
	copy_edge_to_device(dptr, key, Z_SHADOW_BOTTOM, 0);
	cudaDeviceSynchronize();
	//printf("** end comm edge rank %d\n", rank);
}

void ZDivMPIGPUMapper::stage_in_from_slave(int srcrank, void *hptr)
{
	size_t offset=0;
	size_t size=0;

	for(int p=0; p<pagecount; p++)
	{
		MPI_Recv(&offset, 1, MPI_UINT64_T, srcrank, 1002, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&size, 1, MPI_UINT64_T, srcrank, 1002, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		char *h_dst = (char*)hptr + offset;

	//printf("from rank %d got offset %d and size %d writing to base %x loc %x\n", srcrank, offset, size, hptr, h_dst);

		MPI_Recv(h_dst, size, MPI_BYTE, srcrank, 1002, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

void ZDivMPIGPUMapper::stage_out_to_master(void *hptr)
{
	for(int p=0; p<pagecount; p++)
	{
		size_t offset = p * global_size + global_output_offset * cellsize;
		char *h_src = (char*)hptr + offset;
		MPI_Send(&offset, 1, MPI_UINT64_T, 0, 1002, MPI_COMM_WORLD);
		MPI_Send(&authority_size, 1, MPI_UINT64_T, 0, 1002, MPI_COMM_WORLD);
		MPI_Send(h_src, authority_size, MPI_BYTE, 0, 1002, MPI_COMM_WORLD);
		//printf("from rank %d sent offset %d and size %d\n", rank, offset, authority_size);
	}
}

