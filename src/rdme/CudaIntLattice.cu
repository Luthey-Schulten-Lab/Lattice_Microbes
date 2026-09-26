#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
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
 * Author(s): Elijah Roberts, Ron Acda
 *   (Ron Acda: using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
 */

#include "config.h"
#include "core/Types.h"
#include "core/Exceptions.h"
#include "rdme/CudaIntLattice.h"
#include "rdme/Lattice.h"

namespace lm {
namespace rdme {

CudaIntLattice::CudaIntLattice(lattice_coord_t size, si_dist_t latticeSpacing, uint particlesPerSite)
:IntLattice(size,latticeSpacing,particlesPerSite),cudaParticlesCurrent(0),cudaParticlesSize(0),cudaSiteTypesSize(0),cudaSiteTypes(NULL),isGPUMemorySynched(false),hostRegistered(false)
{
    // Initialize the pointers.
    cudaParticles[0] = NULL;
    cudaParticles[1] = NULL;

    // Make sure the lattice dimensions are divisible by 32.
    if (size.x%32 != 0 || size.y%32 != 0 || size.z%32 != 0) throw InvalidArgException("size","each dimension of a CUDA lattice must be divisible by 32");
    allocateCudaMemory();
}

CudaIntLattice::CudaIntLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t latticeSpacing, uint particlesPerSite)
:IntLattice(xSize,ySize,zSize,latticeSpacing,particlesPerSite),cudaParticlesCurrent(0),cudaParticlesSize(0),cudaSiteTypesSize(0),cudaSiteTypes(NULL),isGPUMemorySynched(false),hostRegistered(false)
{
    // Initialize the pointers.
    cudaParticles[0] = NULL;
    cudaParticles[1] = NULL;

    // Make sure the lattice dimensions are divisible by 32.
    if (size.x%32 != 0 || size.y%32 != 0 || size.z%32 != 0) throw InvalidArgException("size","each dimension of a CUDA lattice must be divisible by 32");
    allocateCudaMemory();
}

CudaIntLattice::~CudaIntLattice()
{
    deallocateCudaMemory();
}

void CudaIntLattice::allocateCudaMemory()
{
    // Allocate memory on the CUDA device.
    cudaParticlesSize=numberSites*wordsPerSite*sizeof(uint32_t);
    CUDA_EXCEPTION_CHECK(cudaMalloc(&cudaParticles[0], cudaParticlesSize)); //TODO: track memory usage.
    CUDA_EXCEPTION_CHECK(cudaMalloc(&cudaParticles[1], cudaParticlesSize)); //TODO: track memory usage.
    // both buffers start zeroed (the int solver's per-site occupancy arrays start at 0 for both).
    CUDA_EXCEPTION_CHECK(cudaMemset(cudaParticles[0], 0, cudaParticlesSize));
    CUDA_EXCEPTION_CHECK(cudaMemset(cudaParticles[1], 0, cudaParticlesSize));
    cudaSiteTypesSize=numberSites*sizeof(uint8_t);
    CUDA_EXCEPTION_CHECK(cudaMalloc(&cudaSiteTypes, cudaSiteTypesSize)); //TODO: track memory usage.

    // page-lock the host lattice buffers (allocated once by IntLattice, never reallocated) so the
    // per-hook copyFromGPU/copyToGPU run at pinned-memory bandwidth. Same bytes are copied; if the driver
    // refuses the registration the copies simply stay pageable.
    hostRegistered = (cudaHostRegister(particles, cudaParticlesSize, cudaHostRegisterPortable) == cudaSuccess);
    if (hostRegistered && cudaHostRegister(siteTypes, cudaSiteTypesSize, cudaHostRegisterPortable) != cudaSuccess)
    {
        cudaHostUnregister(particles); hostRegistered = false;
    }
    if (!hostRegistered) cudaGetLastError();   // clear the sticky error of a refused registration
}

void CudaIntLattice::deallocateCudaMemory()
{
    if (wcmDevBuf) { cudaFree(wcmDevBuf); wcmDevBuf = nullptr; wcmDevBufSize = 0; }
    if (hostRegistered)
    {
        cudaHostUnregister(particles);
        cudaHostUnregister(siteTypes);
        hostRegistered = false;
    }
    // If we have any allocated device memory, free it.
    if (cudaParticles[0] != NULL)
    {
        CUDA_EXCEPTION_CHECK(cudaFree(cudaParticles[0])); //TODO: track memory usage.
        cudaParticles[0] = NULL;
    }
    if (cudaParticles[1] != NULL)
    {
        CUDA_EXCEPTION_CHECK(cudaFree(cudaParticles[1])); //TODO: track memory usage.
        cudaParticles[1] = NULL;
    }
    cudaParticlesSize = 0;
    if (cudaSiteTypes != NULL)
    {
        CUDA_EXCEPTION_CHECK(cudaFree(cudaSiteTypes)); //TODO: track memory usage.
        cudaSiteTypes = NULL;
        cudaSiteTypesSize = 0;
    }
}

// lazy particle download ------------------------------------------------------------------------------------------------------------------
__global__ void wcm_ribo_mask_kernel(const unsigned int *p, unsigned int n, unsigned int words, unsigned int ridx, unsigned char *mask)
{
    const unsigned int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    unsigned char m = 0;
    for (unsigned int w = 0; w < words; w++) { const unsigned int v = p[(size_t)w * n + k]; if (v == 0) break; if (v == ridx) { m = 1; break; } }
    mask[k] = m;
}
__global__ void wcm_gather_kernel(const unsigned int *p, unsigned int n, unsigned int words, const int *sites, int nsites, unsigned int *out)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= nsites) return;
    const size_t k = (size_t)sites[j];
    for (unsigned int w = 0; w < words; w++) out[(size_t)w * nsites + j] = p[(size_t)w * n + k];
}
static void *wcm_scratch(void *&buf, size_t &cap, size_t need)
{
    if (need > cap) { if (buf) cudaFree(buf); CUDA_EXCEPTION_CHECK(cudaMalloc(&buf, need)); cap = need; }
    return buf;
}
void CudaIntLattice::wcmRiboSiteMask(unsigned char *mask, int n, unsigned int ridx)
{
    if ((size_t)n != (size_t)numberSites) throw InvalidArgException("mask", "wcmRiboSiteMask: mask length must be the number of sites");
    unsigned char *d = (unsigned char *)wcm_scratch(wcmDevBuf, wcmDevBufSize, (size_t)n);
    wcm_ribo_mask_kernel<<<(n + 255) / 256, 256>>>((const unsigned int *)cudaParticles[cudaParticlesCurrent], (unsigned int)n, (unsigned int)wordsPerSite, ridx, d);
    CUDA_EXCEPTION_CHECK(cudaMemcpy(mask, d, (size_t)n, cudaMemcpyDeviceToHost));
}
void CudaIntLattice::wcmGatherSlots(int *sites, int nsites, unsigned int *slots, int nslots)
{
    if ((size_t)nslots != (size_t)nsites * wordsPerSite) throw InvalidArgException("slots", "wcmGatherSlots: slots length must be nsites x words per site");
    if (nsites == 0) return;
    for (int j = 0; j < nsites; j++) if (sites[j] < 0 || (size_t)sites[j] >= (size_t)numberSites) throw InvalidArgException("sites", "wcmGatherSlots: site index out of range");
    const size_t sb = (size_t)nsites * sizeof(int), ob = (size_t)nslots * sizeof(unsigned int);
    char *d = (char *)wcm_scratch(wcmDevBuf, wcmDevBufSize, sb + ob);
    CUDA_EXCEPTION_CHECK(cudaMemcpy(d, sites, sb, cudaMemcpyHostToDevice));
    wcm_gather_kernel<<<(nsites + 255) / 256, 256>>>((const unsigned int *)cudaParticles[cudaParticlesCurrent], (unsigned int)numberSites, (unsigned int)wordsPerSite,
                                                      (const int *)d, nsites, (unsigned int *)(d + sb));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(slots, d + sb, ob, cudaMemcpyDeviceToHost));
}
void CudaIntLattice::copyFromGPULazy()
{
    CUDA_EXCEPTION_CHECK(cudaMemcpy(siteTypes, cudaSiteTypes, cudaSiteTypesSize, cudaMemcpyDeviceToHost));
    wcmParticlesStale = true;
    isGPUMemorySynched = true;
}
void CudaIntLattice::wcm_sync_host_particles() const
{
    if (!wcmParticlesStale) return;
    CUDA_EXCEPTION_CHECK(cudaMemcpy(particles, cudaParticles[cudaParticlesCurrent], cudaParticlesSize, cudaMemcpyDeviceToHost));
    wcmParticlesStale = false;
}
bool CudaIntLattice::wcmCopyToGPU()
{
    if (isGPUMemorySynched) return false;
    if (wcmParticlesStale)
    {   // the host never read the particles, so nothing can have changed them: the device copy is current
        CUDA_EXCEPTION_CHECK(cudaMemcpy(cudaSiteTypes, siteTypes, cudaSiteTypesSize, cudaMemcpyHostToDevice));
        isGPUMemorySynched = true;
        return false;
    }
    copyToGPU();
    return true;
}
// ----------------------------------------------------------------------------------------------------------------------------

void CudaIntLattice::copyToGPU()
{
	if (!isGPUMemorySynched && wcmParticlesStale)   // never upload host particles that were not downloaded
	{
		CUDA_EXCEPTION_CHECK(cudaMemcpy(cudaSiteTypes, siteTypes, cudaSiteTypesSize, cudaMemcpyHostToDevice));
		isGPUMemorySynched = true;
		return;
	}
	if (!isGPUMemorySynched)
	{
		CUDA_EXCEPTION_CHECK(cudaMemcpy(cudaParticles[cudaParticlesCurrent], particles, cudaParticlesSize, cudaMemcpyHostToDevice));
        CUDA_EXCEPTION_CHECK(cudaMemcpy(cudaSiteTypes, siteTypes, cudaSiteTypesSize, cudaMemcpyHostToDevice));
		isGPUMemorySynched = true;
	}
}

void CudaIntLattice::copySiteTypesToGPU()
{
    static const bool verify = std::getenv("WCM_SITE_ONLY_UPLOAD_VERIFY") != NULL;
    if (verify && !wcmParticlesStale)
    {
        static std::vector<unsigned char> dev;
        static long checked = 0, bad = 0;
        dev.resize(cudaParticlesSize);
        CUDA_EXCEPTION_CHECK(cudaMemcpy(dev.data(), cudaParticles[cudaParticlesCurrent], cudaParticlesSize, cudaMemcpyDeviceToHost));
        const bool same = std::memcmp(dev.data(), particles, cudaParticlesSize) == 0;
        checked++; if (!same) bad++;
        if (!same || checked % 1000 == 1)
            printf("site-types-only upload check: site-types-only upload %ld: host particles %s the GPU's (%ld of %ld differ)\n",
                   checked, same ? "IDENTICAL to" : "DIFFERENT from", bad, checked);
    }
    CUDA_EXCEPTION_CHECK(cudaMemcpy(cudaSiteTypes, siteTypes, cudaSiteTypesSize, cudaMemcpyHostToDevice));
    isGPUMemorySynched = true;
}

void CudaIntLattice::copyFromGPU()
{
	wcmParticlesStale = false;
	CUDA_EXCEPTION_CHECK(cudaMemcpy(particles, cudaParticles[cudaParticlesCurrent], cudaParticlesSize, cudaMemcpyDeviceToHost));
    CUDA_EXCEPTION_CHECK(cudaMemcpy(siteTypes, cudaSiteTypes, cudaSiteTypesSize, cudaMemcpyDeviceToHost));
	isGPUMemorySynched = true;
}

void * CudaIntLattice::getGPUMemorySrc()
{
    return cudaParticles[cudaParticlesCurrent];
}

void * CudaIntLattice::getGPUMemoryDest()
{
    return cudaParticles[cudaParticlesCurrent==0?1:0];
}

void CudaIntLattice::swapSrcDest()
{
    cudaParticlesCurrent = cudaParticlesCurrent==0?1:0;
}

void * CudaIntLattice::getGPUMemorySiteTypes()
{
    return cudaSiteTypes;
}

void CudaIntLattice::setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t site) 
{
    IntLattice::setSiteType(x,y,z,site);
    isGPUMemorySynched = false;
}

void CudaIntLattice::addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle) 
{
    IntLattice::addParticle(x,y,z,particle);
	isGPUMemorySynched = false;
}

void CudaIntLattice::removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z) 
{
    IntLattice::removeParticles(x,y,z);
    isGPUMemorySynched = false;
}

void CudaIntLattice::setSiteType(lattice_size_t index, site_t site) 
{
    IntLattice::setSiteType(index,site);
    isGPUMemorySynched = false;
}

void CudaIntLattice::addParticle(lattice_size_t index, particle_t particle) 
{
    IntLattice::addParticle(index,particle);
	isGPUMemorySynched = false;
}

void CudaIntLattice::removeParticles(lattice_size_t index) 
{
    IntLattice::removeParticles(index);
    isGPUMemorySynched = false;
}

void CudaIntLattice::removeAllParticles()
{
    IntLattice::removeAllParticles();
	isGPUMemorySynched = false;
}

void CudaIntLattice::setFromRowMajorByteData(void * buffer, size_t bufferSize)
{
    IntLattice::setFromRowMajorByteData(buffer, bufferSize);
    isGPUMemorySynched = false;
}

}
}
