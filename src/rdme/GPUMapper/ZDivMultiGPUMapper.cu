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
#include "ZDivMultiGPUMapper.h"

#include "core/Print.h"
using lm::Print;

#define LB_MAX_MEM_THRESHOLD 0.9f
#define LB_IMBALANCE_DAMPENER 0.40f

ZDivMultiGPUMapper::ZDivMultiGPUMapper(dim3 l, size_t m, int a, int c, int n, int *gl, bool pz, int pages)
: MultiGPUMapper(l, m, a, c, n, gl, pages), info(NULL), periodic_z(pz)
{
	info=new gpu_info[num_gpus];

	size_t total_mem=0;
	int fair_chunks=int(round((float)l.z/TUNE_MPD_Z_BLOCK_Z_SIZE)/num_gpus);
	for(int i=0; i<num_gpus; i++)
	{
		total_mem += device_memory[i];
		size_t max_mem=(device_memory[i]*LB_MAX_MEM_THRESHOLD)/2;
		info[i].lb_max_chunks=max_mem/(l.x*l.y*TUNE_MPD_Z_BLOCK_Z_SIZE*sizeof(float));
		info[i].lb_chunks=(info[i].lb_max_chunks < fair_chunks ? info[i].lb_max_chunks : fair_chunks);
		info[i].lb_imbalance=0;
	}

	if(total_mem < l.x*l.y*l.z*sizeof(float)*2)
	{
		throw("Insufficient aggregate memory");
	}

	determine_load_balance();
	initialize();
}

ZDivMultiGPUMapper::ZDivMultiGPUMapper(int x, int y, int z, size_t cellsize, int apron, int overlap, int ngpu, int *gl, bool pz, int pages)
: MultiGPUMapper(dim3(x,y,z), cellsize, apron, overlap, ngpu, gl, pages), info(NULL), periodic_z(pz)
{
	info=new gpu_info[num_gpus];

	size_t total_mem=0;
	int fair_chunks=int(round((float)z/TUNE_MPD_Z_BLOCK_Z_SIZE)/num_gpus);
	for(int i=0; i<num_gpus; i++)
	{
		total_mem += device_memory[i];
		size_t max_mem=(device_memory[i]*LB_MAX_MEM_THRESHOLD)/2;
		info[i].lb_max_chunks=max_mem/(x*y*TUNE_MPD_Z_BLOCK_Z_SIZE*sizeof(float));
		info[i].lb_chunks=(info[i].lb_max_chunks < fair_chunks ? info[i].lb_max_chunks : fair_chunks);
		info[i].lb_imbalance=0;
	}

	if(total_mem < x*y*z*sizeof(float)*2*pages)
	{
		throw("Insufficient aggregate memory");
	}

	determine_load_balance();
	initialize();
}

ZDivMultiGPUMapper::~ZDivMultiGPUMapper()
{
}

void ZDivMultiGPUMapper::initialize()
{
	int zstart=0;

	for(int i=0; i<num_gpus; i++)
	{
		gpu_info g=info[i];
		
		// Determine the length of the lattice as dictated by load balancing
		int zdim=g.lb_chunks*TUNE_MPD_Z_BLOCK_Z_SIZE;
	
		assert(zdim <= lattice_dim.z);

		int top_overlap = ((periodic_z || i > 0) ? overlap : 0);
		int bottom_overlap = ((periodic_z || i < num_gpus-1) ? overlap : 0);

		int top_apron = ((periodic_z || i > 0) ? apron : 0);
		int bottom_apron = ((periodic_z || i < num_gpus-1) ? apron : 0);

		int top_shadow = top_overlap+top_apron;
		int bottom_shadow = bottom_overlap+bottom_apron;

		// Store global position and dimension that this segment represents
		g.global_pos=dim3(0,0,zstart);
		g.global_dim=dim3(lattice_dim.x, lattice_dim.y, zdim);

		// Store the dimensions of the local lattice
		g.local_dim=dim3(lattice_dim.x, lattice_dim.y, zdim+top_shadow+bottom_shadow);
		
		// Local authority is where in the local lattice authoritative data begins
		g.local_authority=dim3(0,0, top_shadow);

		// Dependent position is where the global authority, less shadow region
		// NOTE: can be negative in periodic cases!
		g.dependent_pos.x=0;
		g.dependent_pos.y=0;
		g.dependent_pos.z=zstart - top_shadow;

		// dimensions of the shadow regions
		g.overlap_dim=dim3(lattice_dim.x, lattice_dim.y, overlap+apron);

		// Where my neighbor reads from me
		g.overlap_send[Z_SHADOW_TOP]=dim3(0,0,overlap+apron);
		g.overlap_send[Z_SHADOW_BOTTOM]=dim3(0,0, zdim-(overlap+apron) + top_shadow);
		
		// Where I receive
		g.overlap_recv[Z_SHADOW_TOP]=dim3(0,0,0);
		g.overlap_recv[Z_SHADOW_BOTTOM]=dim3(0,0, zdim + top_shadow);

		g.neighbor[Z_SHADOW_TOP]=i-1;
		g.neighbor[Z_SHADOW_BOTTOM]=i+1;

		// If Z-Periodic, set top and bottom to be neighbors
		if(periodic_z)
		{
			if(i==0)
			{
				g.neighbor[Z_SHADOW_TOP]=num_gpus-1;
			}
			if(i==num_gpus-1)
			{
				g.neighbor[Z_SHADOW_BOTTOM]=0;
			}
		}

		info[i]=g;

		build_descriptor(i, g.local_dim, g.dependent_pos,
			dim3(lattice_dim.x, lattice_dim.y, (zdim+top_overlap+bottom_overlap)),
			dim3(0, 0, top_apron));

		zstart+=zdim;
	}
}

void ZDivMultiGPUMapper::initialize_gpu(int gpu)
{
	MultiGPUMapper::initialize_gpu(gpu);

	// we assume peering is always commutative; so if we can peer with a
	// neighbor, they can peer with us.  If we can peer, then allocate memory
	// here so that they can write to it during publish.  Otherwise, make a
	// host bounce buffer.

	size_t ss=DIMSIZE(info[gpu].overlap_dim)*cellsize;

	// Peer with top neighbor
	int neighbor = info[gpu].neighbor[Z_SHADOW_TOP];
	
	if(neighbor >= 0 && neighbor < num_gpus)	
	{
		neighbor_buffer *nb=new neighbor_buffer;
		nb->size=ss;

		if(enable_peer_access(gpu, neighbor))
		{
			cudaMalloc(&(nb->buffer[0]), ss*pagecount);
			cudaMalloc(&(nb->buffer[1]), ss*pagecount);
			cudaMemset(nb->buffer[0], 255, ss*pagecount);
			cudaMemset(nb->buffer[1], 255, ss*pagecount);
		}
		else
		{
			cudaHostAlloc(&(nb->buffer[0]), ss*pagecount, cudaHostAllocPortable);
			cudaHostAlloc(&(nb->buffer[1]), ss*pagecount, cudaHostAllocPortable);
		}

		info[gpu].read_buffer[Z_SHADOW_TOP]=nb;
		info[neighbor].write_buffer[Z_SHADOW_BOTTOM]=nb;

		// make local temp buffer to collect pages before P2P
		cudaMalloc(&(info[gpu].tmp_buffer[0]), ss*pagecount);
	}

	// Peer with bottom neighbor
	neighbor = info[gpu].neighbor[Z_SHADOW_BOTTOM];
	
	if(neighbor >= 0 && neighbor < num_gpus)	
	{
		neighbor_buffer *nb=new neighbor_buffer;
		nb->size=ss;

		if(enable_peer_access(gpu, neighbor))
		{
			cudaMalloc(&(nb->buffer[0]), ss*pagecount);
			cudaMalloc(&(nb->buffer[1]), ss*pagecount);
			cudaMemset(nb->buffer[0], 255, ss*pagecount);
			cudaMemset(nb->buffer[1], 255, ss*pagecount);
		}
		else
		{
			cudaHostAlloc(&(nb->buffer[0]), ss*pagecount, cudaHostAllocPortable);
			cudaHostAlloc(&(nb->buffer[1]), ss*pagecount, cudaHostAllocPortable);
		}

		info[gpu].read_buffer[Z_SHADOW_BOTTOM]=nb;
		info[neighbor].write_buffer[Z_SHADOW_TOP]=nb;

		// make local temp buffer to collect pages before P2P
		cudaMalloc(&(info[gpu].tmp_buffer[1]), ss*pagecount);
	}

	check_error();
}

dim3 ZDivMultiGPUMapper::get_global_dim(int gpu)
{
	return info[gpu].global_dim;
}

dim3 ZDivMultiGPUMapper::get_local_dim(int gpu)
{
	return info[gpu].local_dim;
}

int3 ZDivMultiGPUMapper::get_global_offset(int gpu)
{
	return info[gpu].dependent_pos;
}

size_t ZDivMultiGPUMapper::get_local_size(int gpu)
{
	return DIMSIZE(info[gpu].local_dim)*cellsize;
}

size_t ZDivMultiGPUMapper::get_authority_size(int gpu)
{
	return DIMSIZE(info[gpu].global_dim)*cellsize;
}

ssize_t ZDivMultiGPUMapper::get_global_input_offset(int gpu)
{
	const int3 dp=info[gpu].dependent_pos;
	const dim3 gd=info[gpu].global_dim;

	return dp.x + dp.y * (int)gd.x + dp.z * (int)gd.x * (int)gd.y;
}

size_t ZDivMultiGPUMapper::get_global_output_offset(int gpu)
{
	return POSDIM(info[gpu].global_pos, info[gpu].global_dim);
}

size_t ZDivMultiGPUMapper::get_authority_offset(int gpu)
{
	return POSDIM(info[gpu].local_authority, info[gpu].local_dim);
}

void ZDivMultiGPUMapper::stage_in(int gpu, void *dptr, void *hptr)
{
	for (int p = 0; p < pagecount; p++)
	{
		char *dpage = (char*)dptr + p*get_local_size(gpu);
		char *hpage = (char*)hptr + p*get_global_size();
		stage_in_real(gpu, dpage, hpage, cellsize);
	}
}

void ZDivMultiGPUMapper::stage_in_sites(int gpu, void *dptr, void *hptr)
{
	stage_in_real(gpu, dptr, hptr, 1);
}

void ZDivMultiGPUMapper::stage_in_real(int gpu, void *dptr, void *hptr, unsigned int element_size) 
{
	PROF_BEGIN(PROF_H2D);
	char *src = (char*)hptr;
	char *dst = (char*)dptr;

	ssize_t offset = get_global_input_offset(gpu) * element_size;
	size_t localsize = DIMSIZE(info[gpu].local_dim) * element_size;
	size_t latsize = DIMSIZE(lattice_dim) * element_size;

	if (offset < 0)
	{
		// Read adj bytes from the end of the lattice
		size_t adj = abs(offset);
		cudaMemcpy(dst, src + latsize - adj, adj, cudaMemcpyDefault);
		check_error();
	
		dst += adj;
		offset = 0;
		localsize -= adj;
	}
	else if ((offset + localsize) >= latsize)
	{
		// read adj bytes from the beginning of the lattice
		size_t adj = (offset + localsize) - latsize;
		cudaMemcpy(dst + localsize - adj, src, adj, cudaMemcpyDefault);
		check_error();
		localsize -= adj;
	}
		
	src += offset;
	cudaMemcpy(dst, src, localsize, cudaMemcpyDefault);
	check_error();
	PROF_END(PROF_H2D);
}

void ZDivMultiGPUMapper::stage_out(int gpu, void *hptr, void *dptr) 
{
	for (int p = 0; p < pagecount; p++)
	{
		char *dpage = (char*)dptr + p*get_local_size(gpu);
		char *hpage = (char*)hptr + p*get_global_size();
		stage_out_real(gpu, hpage, dpage);
	}
}

void ZDivMultiGPUMapper::stage_out_real(int gpu, void *hptr, void *dptr) 
{
	PROF_BEGIN(PROF_D2H);
	char *src = (char*)dptr;
	src += get_authority_offset(gpu) * cellsize;

	char *dst = (char*)hptr;
	dst += get_global_output_offset(gpu) * cellsize;

	cudaMemcpy(dst, src, get_authority_size(gpu), cudaMemcpyDefault);
	check_error();
	PROF_END(PROF_D2H);
}

void ZDivMultiGPUMapper::publish_state(int gpu, int timestep,
                                       cudaStream_t top, cudaStream_t bot,
									   void *dptr)
{
	schedule_send(gpu, dptr, timestep, Z_SHADOW_TOP,    top);
	schedule_send(gpu, dptr, timestep, Z_SHADOW_BOTTOM, bot);
	cudaDeviceSynchronize();
}

void ZDivMultiGPUMapper::schedule_recv(int gpu, void *dptr, int key, int neighbor, cudaStream_t s)
{
	int other = info[gpu].neighbor[neighbor];
	if (other < 0 || other >= num_gpus)
		return;	

	key = key & 1;

	neighbor_buffer *nb=info[gpu].read_buffer[neighbor];
	
	char *lbuf = (char*)dptr;
	size_t lofs = POSDIM(info[gpu].overlap_recv[neighbor], info[gpu].local_dim);
	lbuf += lofs*cellsize;

	PROF_CUDA_BEGIN(PROF_MRBASE + neighbor, s);
	char *nbuf = nb->buffer[key];
	for (int p = 0; p < pagecount; p++)
		cudaMemcpyAsync(lbuf + p*get_local_size(gpu),
	                    nbuf + p*(nb->size),
						nb->size, cudaMemcpyDefault, s);
	PROF_CUDA_END(PROF_MRBASE + neighbor, s);
}

void ZDivMultiGPUMapper::schedule_send(int gpu, void *dptr, int key, int neighbor, cudaStream_t s)
{
	int other           = info[gpu].neighbor[neighbor];
	char *tbuf          = info[gpu].tmp_buffer[neighbor];
	neighbor_buffer *nb = info[gpu].write_buffer[neighbor];

	key = (key & 1) ^ 1;
	char *nbuf = nb->buffer[key];

	if (dptr)
	{
		size_t ofs;
		char *buf;

		ofs = POSDIM(info[gpu].overlap_send[neighbor], info[gpu].local_dim);
		buf = (char*)dptr + ofs*cellsize;

		for (int p = 0; p < pagecount; p++)
			cudaMemcpy(tbuf + p*(nb->size),
		               buf + p*get_local_size(gpu),
					   nb->size, cudaMemcpyDefault);
	}

	cudaMemcpyPeerAsync(nbuf, other, tbuf, gpu, nb->size * pagecount, s);
}

gpu_info* ZDivMultiGPUMapper::getinfo(int gpu)
{
	return &(info[gpu]);
}

unsigned int* ZDivMultiGPUMapper::gettbuf(int gpu, int key, int neighbor)
{
	return (unsigned int*)info[gpu].tmp_buffer[neighbor];
}

unsigned int* ZDivMultiGPUMapper::getrbuf(int gpu, int key, int neighbor)
{
	key = key & 1;

	neighbor_buffer *nb = info[gpu].read_buffer[neighbor];
	return (unsigned int*)nb->buffer[key];
}

/*
void ZDivMultiGPUMapper::copy_to_neighbor(int gpu, int neighbor, int key)
{
	int other=info[gpu].neighbor[neighbor];
	if(other < 0 || other >= num_gpus)
		return;	

	size_t ofs;
	char *buf;
	
	neighbor_buffer *nb=info[gpu].write_buffer[neighbor];

	ofs=POSDIM(info[gpu].overlap_send[neighbor], info[gpu].local_dim);
	buf=(char*)info[gpu].current_state+ofs*cellsize;

	PROF_BEGIN(PROF_D2D);
	cudaMemcpy(nb->buffer[key], buf, nb->size, cudaMemcpyDefault);
	check_error();
	PROF_END(PROF_D2D);
}

void ZDivMultiGPUMapper::copy_from_neighbor(int gpu, int neighbor, void *dptr, int key) 
{
	int other=info[gpu].neighbor[neighbor];
	if(other < 0 || other >= num_gpus)
		return;	

	neighbor_buffer *nb=info[gpu].read_buffer[neighbor];
	
	char *lbuf=(char*)dptr;
	size_t lofs=POSDIM(info[gpu].overlap_recv[neighbor], info[gpu].local_dim);
	lbuf+=lofs*cellsize;

	//printf("gpu %d dir %d from buf %X\n",gpu,neighbor,nptr);
	PROF_BEGIN(PROF_D2D);
	cudaMemcpy(lbuf, nb->buffer[key], nb->size, cudaMemcpyDefault);
	PROF_END(PROF_D2D);
}
*/

void ZDivMultiGPUMapper::refresh(int gpu, void *dptr, int timestep)
{
	schedule_recv(gpu, dptr, timestep, Z_SHADOW_TOP, 0);
	schedule_recv(gpu, dptr, timestep, Z_SHADOW_BOTTOM, 0);
	if(cudaDeviceSynchronize() != cudaSuccess)
		throw("cuda error");
}

int ZDivMultiGPUMapper::map_index_to_gpu(size_t index)
{
	size_t plane=lattice_dim.x * lattice_dim.y;
	size_t zplane = index / plane;

	int i;
	for(i=0; i<num_gpus; i++)
	{
		if(zplane < info[i].global_pos.z)
			break;
	}
			
	return i-1;
}

bool ZDivMultiGPUMapper::determine_load_balance()
{
	bool changed=false;

	// Determine the average runtime among 'active' gpus
	// That is, ones that are not capacity constrained.
	unsigned int sum=0;
	unsigned int active=0;

	for(int i=0; i<num_gpus; i++)
	{
		if(info[i].lb_chunks < info[i].lb_max_chunks)
		{
			active++;
			sum+=lb_cost[i];
		}
	}
	
	float avg=(float)sum/active;
	//printf("* sum = %d, active = %d, avg = %f\n", sum, active, avg);

	// Determine balance for non-constrained gpus
	int zblocks=lattice_dim.z/TUNE_MPD_Z_BLOCK_Z_SIZE;
	float chunk_pct=1.0f/(float)zblocks;
	float equalpoint=1.0f/active;

	for(int i=0; i<num_gpus; i++)
	{
		// skip constrained gpus
		if(info[i].lb_chunks == info[i].lb_max_chunks)
			continue;

		// Calculate percent differnence from the average
		// and scale it against the equal point
		float ib_acc=((lb_cost[i]-avg)/avg)*equalpoint;
		//printf("%d cost %d ca %f eq %f ib now %f ib acc %f\n",
		//	i, lb_cost[i], lb_cost[i]/avg, equalpoint, ib_acc, info[i].lb_imbalance);

		// combine this imbalance with the cumulative 
		// effects of the past, dampened by a factor
		float imbalance=ib_acc+info[i].lb_imbalance*LB_IMBALANCE_DAMPENER;
		info[i].lb_imbalance=imbalance;

		// Compute the magnitude of the number of chunks to adjust
		int magnitude=(int)floor(fabs(imbalance)/chunk_pct);

		if(magnitude > 0)
		{
			int sign=(imbalance > 0 ? -1 : 1);
			int adjust=magnitude*sign;
			alter_chunks(i, adjust);
			changed=true;
			info[i].lb_imbalance=0;
			printf("[%d] %s %d blocks\n",i,(adjust>0 ? "added" : "shaved"),magnitude);
		}
	}
	
	// Check to make sure everything is accounted for
	int blocks_consumed=0;
	for(int i=0; i<num_gpus; i++)
    {
		blocks_consumed+=info[i].lb_chunks;
	}

	if(blocks_consumed != zblocks)
	{
        for(int i=0; i<num_gpus; i++)
        {
            alter_chunks(i, 1);
            blocks_consumed++;
            if(blocks_consumed == zblocks) break;
        }
        
/*
		// Attempt to equally distribute any discrepency
		int correction=(int)round((zblocks-blocks_consumed)/active);
		if(correction > 0)
		{
			blocks_consumed=0;
			for(int i=0; i<num_gpus; i++)
			{
				if(info[i].lb_chunks == info[i].lb_max_chunks)
				{
					blocks_consumed+=info[i].lb_chunks;
					continue;
				}

				blocks_consumed+=alter_chunks(i, correction);
				printf("[%d] augmented %d blocks\n",i,correction);
			}
		}

		// If still uneven, find the most imbalanced and apply difference
		if(blocks_consumed != zblocks)
		{

			int most_off=0;
			float most_off_amt=0.0f;
			for(int i=0; i<num_gpus; i++)
			{
				if(info[i].lb_chunks == info[i].lb_max_chunks)
					continue;

				if(fabs(info[i].lb_imbalance) >= most_off_amt)
				{
					most_off_amt=fabs(info[i].lb_imbalance);
					most_off=i;
				}
			}
			
			alter_chunks(most_off, zblocks-blocks_consumed);
			printf("[%d] forced %d blocks\n",most_off,zblocks-blocks_consumed);
			info[most_off].lb_imbalance=0;
		}

*/
		changed=true;
	}

	Print::printf(Print::DEBUG, "LB State: ");
	for(int i=0; i<num_gpus; i++)
		Print::printf(Print::DEBUG, "%+.02f%% (n=%d) \t", info[i].lb_imbalance*100, info[i].lb_chunks);
	Print::printf(Print::DEBUG, "\n");
		
	return changed;
}

int ZDivMultiGPUMapper::alter_chunks(int gpu, int count)
{
	count+=info[gpu].lb_chunks;

	if(count < 1)
		count=1;

	//printf("[%d] req of %d, max is %d\n", gpu, count, max_chunks);
	count = (info[gpu].lb_max_chunks < count ? info[gpu].lb_max_chunks : count);
	info[gpu].lb_chunks=count;

	return count;
}
			
		
		

