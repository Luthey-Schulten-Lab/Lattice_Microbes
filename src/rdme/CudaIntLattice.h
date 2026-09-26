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

#ifndef LM_RDME_CUDAINTLATTICE_H_
#define LM_RDME_CUDAINTLATTICE_H_

#include <vector>
#include <map>
#include "core/Exceptions.h"
#include "cuda/lm_cuda.h"
#include "rdme/IntLattice.h"
#include "rdme/Lattice.h"

namespace lm {
namespace rdme {

class CudaIntLattice : public IntLattice
{
public:
    CudaIntLattice(lattice_coord_t size, si_dist_t spacing, uint particlesPerSite);
    CudaIntLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing, uint particlesPerSite);
	virtual ~CudaIntLattice();
	
    virtual void copyToGPU();
    virtual void copyFromGPU();
    // upload the site types only (the hook left the particles unchanged); WCM_SITE_ONLY_UPLOAD_VERIFY=1 first compares the
    // host particles with the GPU's, byte for byte, and reports any difference
    virtual void copySiteTypesToGPU();
    // lazy particle download. copyFromGPULazy() downloads the site types only and marks the host particles stale;
    // any host particle access (IntLattice accessors, views, getParticlesMemory) then downloads them first. While they are
    // stale copyToGPU() uploads the site types only (the device particles are the current ones) and returns false.
    virtual void copyFromGPULazy();
    virtual void wcm_sync_host_particles() const;
    bool wcmHostParticlesStale() const { return wcmParticlesStale; }
    bool wcmCopyToGPU();   // copyToGPU(); true when the particles were uploaded
    // the ribosome hook's two particle reads, on the device (no host particles needed):
    //   wcmRiboSiteMask: mask[site] = 1 where a slot holds particle type ridx (slots packed from 0: stops at the first empty)
    //   wcmGatherSlots:  slots[w * nsites + j] = particle word w at site index sites[j]
    void wcmRiboSiteMask(unsigned char *mask, int n, unsigned int ridx);
    void wcmGatherSlots(int *sites, int nsites, unsigned int *slots, int nslots);
    virtual void * getGPUMemorySrc();
    virtual void * getGPUMemoryDest();
    virtual void swapSrcDest();
    virtual void * getGPUMemorySiteTypes();

	// Override methods that can cause the GPU memory to become stale.
	virtual void setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t site);
	virtual void setSiteType(lattice_size_t index, site_t site);
	virtual void addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle);
	virtual void addParticle(lattice_size_t index, particle_t particle);
    virtual void removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z);
    virtual void removeParticles(lattice_size_t index);
	virtual void removeAllParticles();

    // Methods to set the data directly.
    virtual void setFromRowMajorByteData(void * buffer, size_t bufferSize);
	
protected:
    virtual void allocateCudaMemory();
    virtual void deallocateCudaMemory();
	
protected:
    uint cudaParticlesCurrent;
    size_t cudaParticlesSize;
    void * cudaParticles[2];
    size_t cudaSiteTypesSize;
    void * cudaSiteTypes;
    bool isGPUMemorySynched;
    bool hostRegistered;      // host buffers page-locked with cudaHostRegister
    mutable bool wcmParticlesStale = false;
    void *wcmDevBuf = nullptr; size_t wcmDevBufSize = 0;   // device scratch for the mask / gather
};

}
}

#endif
