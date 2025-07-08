/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
 * Copyright 2012 Roberts Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 * 
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
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
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts, Zane Thornburg
 */

#ifndef LM_RDME_MPDRDMESOLVER_H_
#define LM_RDME_MPDRDMESOLVER_H_

#include "cuda/lm_cuda.h"
#include "rdme/RDMESolver.h"
#include "rdme/CudaByteLattice.h"

#define OVERFLOW_MODE_CLASSIC 0
#define OVERFLOW_MODE_RELAXED 1

using lm::rdme::RDMESolver;
using lm::rdme::Lattice;

namespace lm {

namespace io{
class Lattice;
class SpeciesCounts;
}
namespace rdme {

class MpdRdmeSolver : public RDMESolver
{
private:
    using CMESolver::hookSimulation;
    using CMESolver::onBeginTrajectory;
    using CMESolver::onEndTrajectory;
    using RDMESolver::buildDiffusionModel;
    std::vector<unsigned int> maxSiteCounts, maxParticleCounts;
    std::vector<unsigned int> currentMaxSiteCounts, currentMaxParticleCounts;

public:
    MpdRdmeSolver();
    virtual ~MpdRdmeSolver();
    virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
    virtual bool needsReactionModel() {return true;}
    virtual bool needsDiffusionModel()  {return true;}
    virtual void buildModel(const uint numberSpeciesA, const uint numberReactionsA, const uint * initialSpeciesCountsA, const uint * reactionTypeA, const double * kA, const int * SA, const uint * DA, const uint kCols=1);
    virtual void buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData=true);
    virtual void generateTrajectory();
    virtual void setReactionRate(unsigned int rxid, float rate);

protected:
	virtual int hookSimulation(double time, CudaByteLattice *lattice);
	virtual int onBeginTrajectory(CudaByteLattice *lattice);
	virtual int onEndTrajectory(CudaByteLattice *lattice);
	virtual int onWriteLattice(double time, CudaByteLattice *lattice);


    virtual void allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing);
    virtual void writeLatticeData(double time, CudaByteLattice * lattice, lm::io::Lattice * latticeDataSet);
	virtual void writeLatticeSites(double time, CudaByteLattice * lattice);
    virtual void initMaxCounts(CudaByteLattice *lattice);
    virtual void writeMaxCounts();
    virtual void recordSpeciesCounts(double time, CudaByteLattice * lattice, lm::io::SpeciesCounts * speciesCountsDataSet);
    virtual void writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet);
    virtual void hookCheckSimulation(double time, CudaByteLattice * lattice);
    virtual void runTimestep(CudaByteLattice * lattice, uint32_t timestep);
    virtual uint64_t getTimestepSeed(uint32_t timestep, uint32_t substep);
    virtual void computePropensities();
    virtual void copyModelsToDevice();

    #ifdef MPD_CUDA_3D_GRID_LAUNCH
    virtual void calculateXLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateYLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateZLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateReactionLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    #else
    virtual void calculateXLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateYLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateZLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    virtual void calculateReactionLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize);
    #endif

protected:
    uint32_t seed;
    double tau;
    bool reactionModelModified;

    void *       cudaOverflowList;
    cudaStream_t cudaStream;

    uint32_t overflowTimesteps;
    uint32_t overflowListUses;
    int      overflow_handling;

    // Stored model parameters for const memroy
    float* model_reactionRates;
    size_t zeroOrderSize, firstOrderSize, secondOrderSize;
    float *zeroOrder, *firstOrder, *secondOrder;

#ifdef MPD_GLOBAL_S_MATRIX
	uint8_t *RLG;	// Device global memory pointer for RL matrix
	int8_t  *SG;	// Device global memory pointer for S matrix
#endif

#ifdef MPD_GLOBAL_T_MATRIX
	float *TG;
#endif

#ifdef MPD_GLOBAL_R_MATRIX
    float        *reactionRatesG;
    unsigned int *reactionOrdersG;
    unsigned int *reactionSitesG;
    unsigned int *D1G;
    unsigned int *D2G;
#endif

    float *propZeroOrder;
	float *propFirstOrder;
	float *propSecondOrder;
};


namespace mpdrdme_dev {

#ifdef MPD_FREAKYFAST
#ifdef MPD_GLOBAL_R_MATRIX
    __global__ void precomp_reaction_kernel(const unsigned int* inLattice,
                                            const uint8_t * inSites,
                                            unsigned int* outLattice,
                                            const unsigned long long timestepHash,
                                            unsigned int* siteOverflowList,
                                            const __restrict__ int8_t *SG,
                                            const __restrict__ uint8_t *RLG,
                                            const unsigned int* __restrict__ reactionOrdersG,
                                            const unsigned int* __restrict__ reactionSitesG,
                                            const unsigned int* __restrict__ D1G,
                                            const unsigned int* __restrict__ D2G,
                                            const float* reactionRatesG,
                                            const float* __restrict__ qp0,
                                            const float* __restrict__ qp1,
                                            const float* __restrict__ qp2);
#else
    __global__ void precomp_reaction_kernel(const unsigned int* inLattice,
                                            const uint8_t * inSites,
                                            unsigned int* outLattice,
                                            const unsigned long long timestepHash,
                                            unsigned int* siteOverflowList,
                                            const __restrict__ int8_t *SG,
                                            const __restrict__ uint8_t *RLG,
                                            const float* __restrict__ qp0,
                                            const float* __restrict__ qp1,
                                            const float* __restrict__ qp2);
#endif
#endif

#ifdef MPD_CUDA_3D_GRID_LAUNCH
    __global__ void mpd_x_kernel(const unsigned int* inLattice,
                                 const uint8_t * inSites,
                                 unsigned int* outLattice,
                                 const unsigned long long timestepHash,
                                 unsigned int* siteOverflowList);
    __global__ void mpd_y_kernel(const unsigned int* inLattice,
                                 const uint8_t * inSites,
                                 unsigned int* outLattice,
                                 const unsigned long long timestepHash,
                                 unsigned int* siteOverflowList);
    __global__ void mpd_z_kernel(const unsigned int* inLattice,
                                 const uint8_t * inSites,
                                 unsigned int* outLattice,
                                 const unsigned long long timestepHash,
                                 unsigned int* siteOverflowList);
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
    __global__ void reaction_kernel(const unsigned int* inLattice,
                                    const uint8_t * inSites,
                                    unsigned int* outLattice,
                                    const unsigned long long timestepHash,
                                    unsigned int* siteOverflowList,
                                    const int8_t* __restrict__ SG,
                                    const uint8_t* __restrict__ RLG,
                                    const unsigned int* __restrict__ reactionOrdersG,
                                    const unsigned int* __restrict__ reactionSitesG,
                                    const unsigned int* __restrict__ D1G,
                                    const unsigned int* __restrict__ D2G,
                                    const float* __restrict__ reactionRatesG);
#else
    __global__ void reaction_kernel(const unsigned int* inLattice,
                                    const uint8_t * inSites,
                                    unsigned int* outLattice,
                                    const unsigned long long timestepHash,
                                    unsigned int* siteOverflowList,
                                    const int8_t* __restrict__ SG,
                                    const uint8_t* __restrict__ RLG);
#endif
#else
    __global__ void reaction_kernel(const unsigned int* inLattice,
                                    const uint8_t * inSites,
                                    unsigned int* outLattice,
                                    const unsigned long long timestepHash,
                                    unsigned int* siteOverflowList);
#endif

#else
    __global__ void mpd_x_kernel(const unsigned int* inLattice,
                                 const uint8_t * inSites,
                                 unsigned int* outLattice,
                                 const unsigned int gridXSize,
                                 const unsigned long long timestepHash,
                                 unsigned int* siteOverflowList);
    __global__ void mpd_y_kernel(const unsigned int* inLattice,
                                 const uint8_t * inSites,
                                 unsigned int* outLattice,
                                 const unsigned int gridXSize,
                                 const unsigned long long timestepHash,
                                 unsigned int* siteOverflowList);
    __global__ void mpd_z_kernel(const unsigned int* inLattice,
                                 const uint8_t * inSites,
                                 unsigned int* outLattice,
                                 const unsigned int gridXSize,
                                 const unsigned long long timestepHash,
                                 unsigned int* siteOverflowList);
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
    __global__ void reaction_kernel(const unsigned int* inLattice,
                                    const uint8_t * inSites,
                                    unsigned int* outLattice,
                                    const unsigned int gridXSize,
                                    const unsigned long long timestepHash,
                                    unsigned int* siteOverflowList,
                                    const int8_t* __restrict__ SG,
                                    const uint8_t* __restrict__ RLG,
                                    const unsigned int* __restrict__ reactionOrdersG,
                                    const unsigned int* __restrict__ reactionSitesG,
                                    const unsigned int* __restrict__ D1G,
                                    const unsigned int* __restrict__ D2G,
                                    const float* __restrict__ reactionRatesG);
#else
    __global__ void reaction_kernel(const unsigned int* inLattice,
                                    const uint8_t * inSites,
                                    unsigned int* outLattice,
                                    const unsigned int gridXSize,
                                    const unsigned long long timestepHash,
                                    unsigned int* siteOverflowList,
                                    const int8_t* __restrict__ SG,
                                    const uint8_t* __restrict__ RLG);
#endif
#else
    __global__ void reaction_kernel(const unsigned int* inLattice,
                                    const uint8_t * inSites,
                                    unsigned int* outLattice,
                                    const unsigned int gridXSize,
                                    const unsigned long long timestepHash,
                                    unsigned int* siteOverflowList);
#endif
#endif

__global__ void correct_overflows(unsigned int* inLattice, unsigned int* siteOverflowList);
__global__ void correct_overflowsEV(unsigned int* inLattice, const float* evData, unsigned int* siteOverflowList);
__global__ void sanity_check(const unsigned int* L1, const unsigned int* L2);

}
}
}

#endif
