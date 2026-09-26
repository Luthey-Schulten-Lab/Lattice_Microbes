/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
 * Copyright 2012 Roberts Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.uiuc.edu/~schulten
 * 
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
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
 * Author(s): Elijah Roberts, Zane Thornburg, Ron Acda
 *   (Ron Acda: using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
 */

#include <map>
#include <algorithm>
#include <string>
#include <cstdlib>
#include "config.h"
#if defined(MACOSX)
#include <mach/mach_time.h>
#elif defined(LINUX)
#include <time.h>
#endif
#include "cuda/lm_cuda.h"
#include <climits>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "core/Math.h"
#include "core/Print.h"
#include "cme/CMESolver.h"
#include "DiffusionModel.pb.h"
#include "Lattice.pb.h"
#include "SpeciesCounts.pb.h"
#include "core/DataOutputQueue.h"
#include "core/ResourceAllocator.h"
#include "rdme/ByteLattice.h"
#include "rdme/IntLattice.h"
#include "rdme/CudaIntLattice.h"
#include "rdme/IntMpdRdmeSolver.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include "core/Timer.h"

#define MPD_WORDS_PER_SITE    MPD_LATTICE_MAX_OCCUPANCY // 16 
#define MPD_APRON_SIZE        1
#ifndef WCM_RXN_CHUNK
#define WCM_RXN_CHUNK 8
#endif

#include "cuda/constant.cuh"
#if defined(MPD_FREAKYFAST) && defined(MPD_GLOBAL_S_MATRIX) && defined(MPD_GLOBAL_R_MATRIX)
#define WCM_FUSE_RXN_AVAILABLE 1   // the fused kernel can take the reaction step (precomp_reaction_kernel path)
#endif

namespace lm {
namespace rdme {
namespace intmpdrdme_dev {
#include "rdme/dev/xor_random_dev.cu"
#include "rdme/dev/lattice_sim_1d_dev.cu"
#include "rdme/dev/word_diffusion_1d_dev.cu"
#include "rdme/dev/word_reaction_dev.cu"
#include "rdme/dev/word_occ_dev.cu"
}}}

namespace lm { namespace rdme { namespace intmpdrdme_dev {
// occupancy of a buffer just uploaded from the host = highest nonzero slot + 1 per site.
__global__ void recount_occupancy_kernel(const unsigned int* lattice, uint8_t* occ, const unsigned int numberSites)
{
    const unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= numberSites) return;
    unsigned int hi = 0;
    for (unsigned int w=0; w<MPD_WORDS_PER_SITE; w++)
        if (lattice[i + w*numberSites] != 0) hi = w+1;
    occ[i] = (uint8_t)hi;
}

// bounding box (inclusive x0,x1,y0,y1,z0,z1) of every site that holds or can gain a particle: sites of a live type
// (some species can diffuse into it, or it has zero-order reactions; types >= 32 always count as live) and sites occupied in
// either particle buffer. Outside it both buffers are empty and stay empty (nothing can move into or appear in a dead-type site),
// so a diffusion or reaction block whose own sites all lie outside it would only rewrite zeros over zeros: it exits at once.
__device__ int wcm_boxD[6];
// reaction lists built from the model (IntMpdRdmeSolver::buildModel). Per particle value v (species+1), the reactions
// with v as a reactant (D1 or D2), ascending; the zero-order reactions; per reaction its reactants and products as
// (species << 8 | count), ascending species. wcm_rxnListsOn = 0: the original determineReactionIndex / evaluateReaction.
__device__ const unsigned int * wcm_specRxnOff; __device__ const unsigned int * wcm_specRxn;
__device__ const unsigned int * wcm_zeroRxn;    __device__ unsigned int wcm_nZeroRxn;
__device__ const unsigned int * wcm_reacOff;    __device__ const unsigned int * wcm_reac;
__device__ const unsigned int * wcm_prodOff;    __device__ const unsigned int * wcm_prod;
__device__ int wcm_rxnListsOn;

__global__ void wcm_box_init_kernel(int* box, int x0, int x1, int y0, int y1, int z0, int z1)
{
    box[0] = x0; box[1] = x1; box[2] = y0; box[3] = y1; box[4] = z0; box[5] = z1;
}

// warp-wide integer min / max: the sm_80+ reduction instructions, a shuffle tree on older GPUs (same result)
__device__ __forceinline__ int wcm_warp_min(int v) {
#if __CUDA_ARCH__ >= 800
    return __reduce_min_sync(0xffffffffu, v);
#else
    for (int o = 16; o > 0; o >>= 1) v = min(v, __shfl_xor_sync(0xffffffffu, v, o));
    return v;
#endif
}
__device__ __forceinline__ int wcm_warp_max(int v) {
#if __CUDA_ARCH__ >= 800
    return __reduce_max_sync(0xffffffffu, v);
#else
    for (int o = 16; o > 0; o >>= 1) v = max(v, __shfl_xor_sync(0xffffffffu, v, o));
    return v;
#endif
}

__global__ void __launch_bounds__(256) wcm_box_kernel(const uint8_t* __restrict__ sites, const uint8_t* __restrict__ occA, const uint8_t* __restrict__ occB, const unsigned int liveMask, const unsigned int numberSites, const unsigned int X, const unsigned int XY, int* box)
{
    const unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    bool live = false; int x = 0, y = 0, z = 0;
    if (i < numberSites)
    {
        const unsigned int t = sites[i];
        live = (t >= 32u) || ((liveMask >> t) & 1u) || occA[i] != 0 || occB[i] != 0;
        x = (int)(i % X); y = (int)((i % XY) / X); z = (int)(i / XY);
    }
    int v[6] = { live ? x : INT_MAX, live ? x : -1, live ? y : INT_MAX, live ? y : -1, live ? z : INT_MAX, live ? z : -1 };
    #pragma unroll
    for (int k = 0; k < 6; k++) v[k] = (k & 1) ? wcm_warp_max(v[k]) : wcm_warp_min(v[k]);
    __shared__ int sv[8][6];
    const int warp = threadIdx.x >> 5, lane = threadIdx.x & 31;
    if (lane == 0) for (int k = 0; k < 6; k++) sv[warp][k] = v[k];
    __syncthreads();
    if (threadIdx.x == 0)
    {
        for (int w = 1; w < (int)(blockDim.x >> 5); w++)
            for (int k = 0; k < 6; k++) v[k] = (k & 1) ? max(v[k], sv[w][k]) : min(v[k], sv[w][k]);
        if (v[1] >= 0)
        {
            atomicMin(&box[0], v[0]); atomicMax(&box[1], v[1]); atomicMin(&box[2], v[2]);
            atomicMax(&box[3], v[3]); atomicMin(&box[4], v[4]); atomicMax(&box[5], v[5]);
        }
    }
}

// true when none of the block's own sites [x0,x1] x [y0,y1] x [z0,z1] lies in the box
inline __device__ bool wcm_block_outside_box(const int x0, const int x1, const int y0, const int y1, const int z0, const int z1)
{
    return x1 < wcm_boxD[0] || x0 > wcm_boxD[1] || y1 < wcm_boxD[2] || y0 > wcm_boxD[3] || z1 < wcm_boxD[4] || z0 > wcm_boxD[5];
}
}}}

// One occupancy array per device particle buffer, found by the buffer's address (the lattice swaps src/dest every pass).
static void *   wcm_occ_buf[2] = {NULL, NULL};
static uint8_t * wcm_occ_arr[2] = {NULL, NULL};
static const void * wcm_occ_owner = NULL;   // the lattice the arrays belong to (a new solver/lattice in the same process re-initialises)
static uint8_t * wcm_occ(lm::rdme::CudaIntLattice * lattice, void * buf)
{
    if (wcm_occ_owner != (const void *)lattice)
    {
        for (int b = 0; b < 2; b++) if (wcm_occ_arr[b] != NULL) { cudaFree(wcm_occ_arr[b]); wcm_occ_arr[b] = NULL; }
        wcm_occ_owner = (const void *)lattice;
        wcm_occ_buf[0] = lattice->getGPUMemorySrc(); wcm_occ_buf[1] = lattice->getGPUMemoryDest();
        const size_t n = lattice->getNumberSites();
        for (int b = 0; b < 2; b++)
        {
            CUDA_EXCEPTION_CHECK(cudaMalloc(&wcm_occ_arr[b], n));
            CUDA_EXCEPTION_CHECK(cudaMemset(wcm_occ_arr[b], 0, n));   // both particle buffers start zeroed (CudaIntLattice)
        }
    }
    if (buf == wcm_occ_buf[0]) return wcm_occ_arr[0];
    if (buf == wcm_occ_buf[1]) return wcm_occ_arr[1];
    throw lm::Exception("unknown particle buffer");
}
static void wcm_recount_occupancy(lm::rdme::CudaIntLattice * lattice, cudaStream_t stream)
{
    const unsigned int n = (unsigned int)lattice->getNumberSites();
    CUDA_EXCEPTION_EXECUTE((lm::rdme::intmpdrdme_dev::recount_occupancy_kernel<<<(n+255)/256, 256, 0, stream>>>((const unsigned int *)lattice->getGPUMemorySrc(), wcm_occ(lattice, lattice->getGPUMemorySrc()), n)));
}

// site types into which some species can diffuse (bit t), from the transition table; zero-order reactions are added
// from the current propensities when the box is computed. WCM_RDME_BOX_OFF=1: the box is the whole lattice.
static unsigned int wcm_diffLiveTypes = 0xFFFFFFFFu;
// the fused kernel's tile for the next steps, chosen from the box of the previous wcm_update_box: its copy into
// pinned host memory was queued at the end of the previous batch (or single step) and is complete after that batch's stream
// synchronisation. the first tile in the order 24x6x6, 17x7x7, 16x4x4 whose tiles over the box all fit on the GPU at
// once (tiles <= SMs x resident blocks per SM from the occupancy calculator; 16x4x4 when none fits). Measured on 4DWCM cell-cycle states
// (fused kernel us, 24x6x6 / 17x7x7 / 16x4x4): t=900 box 47^3 23.5 / 25.4 / 26.8; t=1800-2400 box 49^3 27.3-28.2 / 25.3-26.1 /
// 27.5-28.3; t>=2700 box >= 51^3 33.5-39.0 / 31.5-38.9 / 27.2-28.6. All tiles give the same lattice.
static int  wcm_fuse_auto = 0;
static int * wcm_box_host = NULL;
static long wcm_tiles_over_box(const lm::rdme::CudaIntLattice * lattice, const int cx, const int cy, const int cz)
{
    const auto size = lattice->getSize();
    const int L[3] = { (int)size.x, (int)size.y, (int)size.z }, C[3] = { cx, cy, cz };
    long tiles = 1;
    for (int d = 0; d < 3; d++)
    {
        const int lo = std::max(wcm_box_host[2*d], 0), hi = std::min(wcm_box_host[2*d+1], L[d]-1);
        tiles *= (hi < lo) ? 0 : (hi - lo) / C[d] + 1;
    }
    return tiles;
}
namespace lm { namespace rdme { template<int CX, int CY, int CZ, int NT, int MINB> long wcm_resident_tiles(); } }   // defined with the fused launcher
static void wcm_pick_tile(const lm::rdme::CudaIntLattice * lattice)
{
    static int nsm = 0, changes = 0;
    static long slots0 = 0, slots3 = 0;
    if (nsm == 0)
    {
        int dev = 0; CUDA_EXCEPTION_CHECK(cudaGetDevice(&dev)); CUDA_EXCEPTION_CHECK(cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev));
        slots0 = nsm * lm::rdme::wcm_resident_tiles<24,6,6,1024,1>(); slots3 = nsm * lm::rdme::wcm_resident_tiles<17,7,7,1024,1>();
    }
    if (wcm_box_host == NULL) return;   // no box yet: keep the default
    const long t0 = wcm_tiles_over_box(lattice, 24, 6, 6), t3 = wcm_tiles_over_box(lattice, 17, 7, 7);
    const int v = (t0 <= slots0) ? 0 : ((t3 <= slots3) ? 3 : 2);
    if (v != wcm_fuse_auto && changes < 20)
    {
        changes++;
        lm::Print::printf(lm::Print::INFO, "fused tile %s (tiles 24x6x6 %ld of %ld resident, 17x7x7 %ld of %ld; box x %d-%d y %d-%d z %d-%d)", v == 0 ? "24x6x6" : (v == 3 ? "17x7x7" : "16x4x4"),
            t0, slots0, t3, slots3, wcm_box_host[0], wcm_box_host[1], wcm_box_host[2], wcm_box_host[3], wcm_box_host[4], wcm_box_host[5]);
    }
    wcm_fuse_auto = v;
}
// The copy is queued after the batch's last kernel (wcm_queue_box_copy, just before the batch's stream synchronisation), not
// between the box kernel and the first step; the box does not change within a batch. The very first box is read at once.
static int * wcm_box_dev = NULL;
static void wcm_copy_box_to_host(const lm::rdme::CudaIntLattice * lattice, int * box, cudaStream_t stream)
{
    wcm_box_dev = box;
    if (wcm_box_host == NULL)
    {
        CUDA_EXCEPTION_CHECK(cudaMallocHost(&wcm_box_host, 6*sizeof(int)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(wcm_box_host, box, 6*sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(stream)); wcm_pick_tile(lattice);
    }
}
static void wcm_queue_box_copy(cudaStream_t stream)
{
    if (wcm_box_host != NULL && wcm_box_dev != NULL)
        CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(wcm_box_host, wcm_box_dev, 6*sizeof(int), cudaMemcpyDeviceToHost, stream));
}

static void wcm_update_box(lm::rdme::CudaIntLattice * lattice, cudaStream_t stream, const float * zeroOrder, size_t zeroOrderSize)
{
    int * box = NULL;
    CUDA_EXCEPTION_CHECK(cudaGetSymbolAddress((void **)&box, lm::rdme::intmpdrdme_dev::wcm_boxD));
    const auto size = lattice->getSize();
    wcm_pick_tile(lattice);
    static const bool off = getenv("WCM_RDME_BOX_OFF") != NULL;
    if (off)
    {
        CUDA_EXCEPTION_EXECUTE((lm::rdme::intmpdrdme_dev::wcm_box_init_kernel<<<1, 1, 0, stream>>>(box, 0, (int)size.x-1, 0, (int)size.y-1, 0, (int)size.z-1)));
        wcm_copy_box_to_host(lattice, box, stream);
        return;
    }
    unsigned int live = wcm_diffLiveTypes;
    for (size_t t = 0; t < zeroOrderSize && t < 32; t++) if (zeroOrder[t] != 0.0f) live |= (1u << t);
    const unsigned int n = (unsigned int)lattice->getNumberSites();
    CUDA_EXCEPTION_EXECUTE((lm::rdme::intmpdrdme_dev::wcm_box_init_kernel<<<1, 1, 0, stream>>>(box, INT_MAX, -1, INT_MAX, -1, INT_MAX, -1)));
    CUDA_EXCEPTION_EXECUTE((lm::rdme::intmpdrdme_dev::wcm_box_kernel<<<(n+255)/256, 256, 0, stream>>>((const uint8_t *)lattice->getGPUMemorySiteTypes(),
        wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()), live, n, (unsigned int)size.x, (unsigned int)(size.x*size.y), box)));
    wcm_copy_box_to_host(lattice, box, stream);
}

// per-species reactant flags (IntMpdRdmeSolver::computePropensities)
static uint8_t * wcm_reactiveG = NULL;

// device copy of the source particle buffer at the start of a batch, and the batch length
static void *  wcm_snap = NULL;
static size_t  wcm_snap_bytes = 0;
static uint32_t wcm_batch_steps()
{
    static int k = -1;
    if (k < 0) { const char * e = getenv("WCM_RDME_BATCH"); k = (e && atoi(e) > 0) ? atoi(e) : 32; }
    return (uint32_t)k;
}

extern bool globalAbort;

using std::map;
using lm::io::DiffusionModel;
using lm::rdme::Lattice;
using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {

IntMpdRdmeSolver::IntMpdRdmeSolver()
:RDMESolver(lm::rng::RandomGenerator::NONE),
 seed(0), tau(0.0),
 cudaOverflowList(NULL), cudaStream(0),
 overflowTimesteps(0), overflowListUses(0),
 model_reactionRates(NULL),
 zeroOrder(NULL), firstOrder(NULL), secondOrder(NULL)
{
}

void IntMpdRdmeSolver::initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources)
{
	RDMESolver::initialize(replicate, parameters, resources);

    // Figure out the random seed.
    uint32_t seedTop=(unsigned int)atoi((*parameters)["seed"].c_str());
    if (seedTop == 0)
    {
        #if defined(MACOSX)
        seedTop = (uint32_t)mach_absolute_time();
        #elif defined(LINUX)
        struct timespec seed_timespec;
        if (clock_gettime(CLOCK_REALTIME, &seed_timespec) != 0) throw lm::Exception("Error getting time to use for random seed.");
        seedTop = seed_timespec.tv_nsec;
        #endif
    }
    seed = (seedTop<<16)|(replicate&0x0000FFFF);

#ifdef MPD_MAPPED_OVERFLOWS
    // mapped, pinned host memory as in MpdRdmeSolver: the host reads the list after the stream sync
    // instead of a synchronous cudaMemcpy of it every timestep.
    CUDA_EXCEPTION_CHECK(cudaHostAlloc(&cudaOverflowList, MPD_OVERFLOW_LIST_SIZE, cudaHostAllocPortable|cudaHostAllocMapped));
    memset(cudaOverflowList, 0, MPD_OVERFLOW_LIST_SIZE);
#else
    // Allocate memory on the device for the exception list.
    CUDA_EXCEPTION_CHECK(cudaMalloc(&cudaOverflowList, MPD_OVERFLOW_LIST_SIZE)); //TODO: track memory usage.
    CUDA_EXCEPTION_CHECK(cudaMemset(cudaOverflowList, 0, MPD_OVERFLOW_LIST_SIZE));
#endif

    // Create a stream for synchronizing the events.
    CUDA_EXCEPTION_CHECK(cudaStreamCreate(&cudaStream));
}

IntMpdRdmeSolver::~IntMpdRdmeSolver()
{
    // Free any device memory.
    if (cudaOverflowList != NULL)
    {
#ifdef MPD_MAPPED_OVERFLOWS
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFreeHost(cudaOverflowList));
#else
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaFree(cudaOverflowList)); //TODO: track memory usage.
#endif
        cudaOverflowList = NULL;
    }

    // If we have created a stream, destroy it.
    if (cudaStream != 0)
    {
        CUDA_EXCEPTION_CHECK_NOTHROW(cudaStreamDestroy(cudaStream));
        cudaStream = NULL;
    }

    // Free allocated memory in IntMpdRdmeSolver
    if (model_reactionRates)  {delete [] model_reactionRates;}

    if (zeroOrder)   {delete [] zeroOrder;}
	if (firstOrder)  {delete [] firstOrder;}
	if (secondOrder) {delete [] secondOrder;}
}

void IntMpdRdmeSolver::allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing)
{
	if(bytes_per_particle != 4)
	{
		Print::printf(Print::ERROR, "IntMpdRdmeSolver works only with lattices of 4 bytes per particle, not %d", bytes_per_particle);
		throw Exception("incorrect lattice data type");
	}

	if(particlesPerSite != MPD_LATTICE_MAX_OCCUPANCY)
	{
		Print::printf(Print::ERROR, "requested allocation for %d particles per site is not %d", particlesPerSite, MPD_LATTICE_MAX_OCCUPANCY);
		throw Exception("incorrect particle density");
	}

    lattice = (Lattice *)new CudaIntLattice(latticeXSize, latticeYSize, latticeZSize, latticeSpacing, particlesPerSite);
}

void IntMpdRdmeSolver::buildModel(const uint numberSpeciesA,
                                  const uint numberReactionsA,
                                  const uint * initialSpeciesCountsA,
                                  const uint * reactionTypesA,
                                  const double * KA,
                                  const int * SA,
                                  const uint * DA,
                                  const uint kCols)
{
    CMESolver::buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypesA, KA, SA, DA, kCols);

    // Get the time step.
    tau=atof((*parameters)["timestep"].c_str());
    if (tau <= 0.0) throw InvalidArgException("timestep", "A positive timestep must be specified for the solver.");

    // Make sure we can support the reaction model.
    if (numberReactions > MPD_MAX_REACTION_TABLE_ENTRIES) throw Exception("The number of reaction table entries exceeds the maximum supported by the solver.");
#ifndef MPD_GLOBAL_S_MATRIX
    if (numberSpecies*numberReactions > MPD_MAX_S_MATRIX_ENTRIES) throw Exception("The number of S matrix entries exceeds the maximum supported by the solver.");
#endif

    // Setup the cuda reaction model.
    unsigned int * reactionOrders = new unsigned int[numberReactions];
    unsigned int * reactionSites = new unsigned int[numberReactions];
    unsigned int * D1 = new unsigned int[numberReactions];
    unsigned int * D2 = new unsigned int[numberReactions];

    for (uint i=0; i<numberReactions; i++)
    {
	    if(reactionTypes[i] == ZerothOrderPropensityArgs::REACTION_TYPE)
        {
    	    reactionOrders[i] = MPD_ZERO_ORDER_REACTION;
    		reactionSites[i] = 0;
    		D1[i] = 0; 
    		D2[i] = 0;
	    }
    	else if (reactionTypes[i] == FirstOrderPropensityArgs::REACTION_TYPE)
    	{
    		reactionOrders[i] = MPD_FIRST_ORDER_REACTION;
    		reactionSites[i] = 0;
    		D1[i] = ((FirstOrderPropensityArgs *)propensityFunctionArgs[i])->si+1;
    		D2[i] = 0;
    	}
    	else if (reactionTypes[i] == SecondOrderPropensityArgs::REACTION_TYPE)
    	{
    		reactionOrders[i] = MPD_SECOND_ORDER_REACTION;
    		reactionSites[i] = 0;
    		D1[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s1i+1;
    		D2[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s2i+1;
    	}
    	else if (reactionTypes[i] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
    	{
    		reactionOrders[i] = MPD_SECOND_ORDER_SELF_REACTION;
    		reactionSites[i] = 0;
    		D1[i] = ((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i])->si+1;
    		D2[i] = 0;
    	}
    	else
    	{
    		throw InvalidArgException("reactionTypeA", "the reaction type was not supported by the solver", reactionTypes[i]);
    	}
    }

    // Setup the cuda S matrix.
	// Transpose S on the device because we will want to read along the numSpecies axis
	// Before A1A2A3A4A5 B1B2B3B4B5 C1C2C3C4C5 ...
	// After A1B1C1 A2B2C2 ...
/*
    int8_t * tmpS = new int8_t[numberSpecies*numberReactions];
    for (uint i=0; i<numberSpecies*numberReactions; i++)
    {
    	tmpS[i] = S[i];
    }
*/
    int8_t * tmpS = new int8_t[numberSpecies*numberReactions];
	for(uint rx = 0; rx < numberReactions; rx++)
	{
		for (uint p=0; p<numberSpecies; p++)
		{
			//tmpS[numberSpecies*rx + p] = S[numberSpecies*rx + p];
			tmpS[rx * numberSpecies + p] = S[numberReactions*p + rx];
		}
	}

    // Copy the reaction model and S matrix to constant memory on the GPU.
#ifdef MPD_GLOBAL_R_MATRIX
    // R matrix put in global memory
    //cudaMalloc(&numberReactionsG, sizeof(unsigned int));
    //cudaMemcpy(numberReactionsG, &numberReactions, sizeof(unsigned int), cudaMemcpyHostToDevice);
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberReactionsC, &numberReactions, sizeof(unsigned int)));
    cudaMalloc(&reactionOrdersG, numberReactions*sizeof(unsigned int));
    cudaMemcpy(reactionOrdersG, reactionOrders, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&reactionSitesG, numberReactions*sizeof(unsigned int));
    cudaMemcpy(reactionSitesG, reactionSites, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&D1G, numberReactions*sizeof(unsigned int));
    cudaMemcpy(D1G, D1, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc(&D2G, numberReactions*sizeof(unsigned int));
    cudaMemcpy(D2G, D2, numberReactions*sizeof(unsigned int), cudaMemcpyHostToDevice);
#else
    // R matrix stored in constant memory
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberReactionsC, &numberReactions, sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionOrdersC, reactionOrders, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionSitesC, reactionSites, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(D1C, D1, numberReactions*sizeof(unsigned int)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(D2C, D2, numberReactions*sizeof(unsigned int)));
#endif

    // reaction lists for the reacting-site path of precomp_reaction_kernel (see wcm_determineReactionIndex /
    // wcm_evaluateReaction). WCM_RDME_RXN_LISTS_OFF=1, or a reaction with more than 8 reactant species or an entry beyond
    // 255 / species beyond 2^24: the original functions.
    {
        std::vector<std::vector<unsigned int>> bySpecies(numberSpecies + 2);
        std::vector<unsigned int> zeroList, reacOff(1, 0), reac, prodOff(1, 0), prod;
        bool ok = getenv("WCM_RDME_RXN_LISTS_OFF") == NULL && numberSpecies < (1u << 24);
        // a second-order reaction (two different reactants) has a nonzero propensity only with both present, so it is
        // listed under the reactant that takes part in fewer reactions (ties: the lower value); the hub species (ribosomes,
        // polymerases: up to 1012 reactions) keep short lists.
        std::vector<unsigned int> memberships(numberSpecies + 2, 0);
        for (uint rx=0; rx<numberReactions; rx++)
            if (reactionOrders[rx] != MPD_ZERO_ORDER_REACTION)
            {
                if (D1[rx] > 0 && D1[rx] <= numberSpecies) memberships[D1[rx]]++;
                if (D2[rx] > 0 && D2[rx] <= numberSpecies && D2[rx] != D1[rx]) memberships[D2[rx]]++;
            }
        const bool rarer = getenv("WCM_RDME_RXN_LISTS_BOTH") == NULL;
        for (uint rx=0; rx<numberReactions; rx++)
        {
            if (reactionOrders[rx] == MPD_ZERO_ORDER_REACTION) zeroList.push_back(rx);
            else if (rarer && reactionOrders[rx] == MPD_SECOND_ORDER_REACTION && D1[rx] > 0 && D2[rx] > 0 && D1[rx] <= numberSpecies && D2[rx] <= numberSpecies && D1[rx] != D2[rx])
            {
                const unsigned int a = D1[rx], b = D2[rx];
                const bool pick_a = (memberships[a] < memberships[b]) || (memberships[a] == memberships[b] && a < b);
                bySpecies[pick_a ? a : b].push_back(rx);
            }
            else
            {
                if (D1[rx] > 0 && D1[rx] <= numberSpecies) bySpecies[D1[rx]].push_back(rx);
                if (D2[rx] > 0 && D2[rx] <= numberSpecies && D2[rx] != D1[rx]) bySpecies[D2[rx]].push_back(rx);
            }
            int nneg = 0;
            for (uint sp=0; sp<numberSpecies; sp++)
            {
                const int v = S[numberReactions*sp + rx];
                if (v < 0) { reac.push_back((sp << 8) | (unsigned int)(-v)); nneg++; if (-v > 255) ok = false; }
                else if (v > 0) { prod.push_back((sp << 8) | (unsigned int)v); if (v > 255) ok = false; }
            }
            if (nneg > 8) ok = false;
            reacOff.push_back((unsigned int)reac.size()); prodOff.push_back((unsigned int)prod.size());
        }
        std::vector<unsigned int> specOff(numberSpecies + 2, 0), specList;
        for (uint v=0; v<=numberSpecies; v++) { specOff[v] = (unsigned int)specList.size(); specList.insert(specList.end(), bySpecies[v].begin(), bySpecies[v].end()); }
        specOff[numberSpecies+1] = (unsigned int)specList.size();
        auto up = [](const std::vector<unsigned int> &h) -> unsigned int * {
            unsigned int *d = NULL;
            CUDA_EXCEPTION_CHECK(cudaMalloc(&d, std::max<size_t>(1, h.size()) * sizeof(unsigned int)));
            if (!h.empty()) CUDA_EXCEPTION_CHECK(cudaMemcpy(d, h.data(), h.size() * sizeof(unsigned int), cudaMemcpyHostToDevice));
            return d;
        };
        const unsigned int *d_specOff = up(specOff), *d_spec = up(specList), *d_zero = up(zeroList), *d_reacOff = up(reacOff), *d_reac = up(reac), *d_prodOff = up(prodOff), *d_prod = up(prod);
        const unsigned int nz = (unsigned int)zeroList.size();
        const int on = ok ? 1 : 0;
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_specRxnOff, &d_specOff, sizeof(d_specOff)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_specRxn, &d_spec, sizeof(d_spec)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_zeroRxn, &d_zero, sizeof(d_zero)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_nZeroRxn, &nz, sizeof(nz)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_reacOff, &d_reacOff, sizeof(d_reacOff)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_reac, &d_reac, sizeof(d_reac)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_prodOff, &d_prodOff, sizeof(d_prodOff)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_prod, &d_prod, sizeof(d_prod)));
        CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(intmpdrdme_dev::wcm_rxnListsOn, &on, sizeof(on)));
        Print::printf(Print::INFO, "reaction lists %s (%u species-reaction entries, %u zero-order, %zu reactant / %zu product entries)", on ? "on" : "off", (unsigned int)specList.size(), nz, reac.size(), prod.size());
    }

#ifdef MPD_GLOBAL_S_MATRIX
	// If S matrix stored in global memory, allocate space and perform copy
	cudaMalloc(&SG, numberSpecies*numberReactions * sizeof(int8_t));
	cudaMemcpy(SG, tmpS, numberSpecies*numberReactions * sizeof(int8_t), cudaMemcpyHostToDevice);
#else
	// S matrix is in constant memory
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(SC, tmpS, numberSpecies*numberReactions*sizeof(int8_t)));
#endif

    // Free any temporary resources.
    delete [] reactionSites;
    delete [] D1;
    delete [] D2;
    delete [] reactionOrders;
    delete [] tmpS;
}

void IntMpdRdmeSolver::buildDiffusionModel(const uint numberSiteTypesA,
                                           const double * DFA,
                                           const uint * RLA,
                                           lattice_size_t latticeXSize,
                                           lattice_size_t latticeYSize,
                                           lattice_size_t latticeZSize,
                                           site_size_t particlesPerSite,
                                           const unsigned int bytes_per_particle,
                                           si_dist_t latticeSpacing,
                                           const uint8_t * latticeData,
                                           const uint8_t * latticeSitesData,
                                           bool rowMajorData)
{
    RDMESolver::buildDiffusionModel(numberSiteTypesA, DFA, RLA,
                                    latticeXSize, latticeYSize, latticeZSize,
                                    particlesPerSite, bytes_per_particle,
                                    latticeSpacing, latticeData, latticeSitesData, rowMajorData);

    // Get the time step.
    tau=atof((*parameters)["timestep"].c_str());
    if (tau <= 0.0) throw InvalidArgException("timestep", "A positive timestep must be specified for the solver.");

    // Setup the cuda transition matrix.
    const size_t DFmatrixSize = numberSpecies*numberSiteTypes*numberSiteTypes;
    if (DFmatrixSize > MPD_MAX_TRANSITION_TABLE_ENTRIES) throw Exception("The number of transition table entries exceeds the maximum supported by the solver.");
    #ifndef MPD_GLOBAL_S_MATRIX
    if (numberReactions*numberSiteTypes > MPD_MAX_RL_MATRIX_ENTRIES) throw Exception("The number of RL matrix entries exceeds the maximum supported by the solver.");
    #endif

    #ifndef MPD_GLOBAL_T_MATRIX
    if (DFmatrixSize > MPD_MAX_TRANSITION_TABLE_ENTRIES) throw Exception("The number of transition table entries exceeds the maximum supported by the solver.");
    #endif

    // Calculate the probability from the diffusion coefficient and the lattice properties.
    /*
     * p0 = probability of staying at the site, q = probability of moving in plus or minus direction
     *
     * D=(1-p0)*lambda^2/2*tau
     * q=(1-p0)/2
     * D=2q*lambda^2/2*tau
     * q=D*tau/lambda^2
     */
    float * T = new float[DFmatrixSize];
    for (uint i=0; i<DFmatrixSize; i++)
    {
        float q=(float)(DF[i]*tau/pow(latticeSpacing,2));
        if (q > 0.50f) throw InvalidArgException("D", "The specified diffusion coefficient is too high for the diffusion model.");
        T[i] = q;
    }

    // Setup the cuda reaction location matrix.
    uint8_t * tmpRL = new uint8_t[numberReactions*numberSiteTypes];
    for (uint i=0; i<numberReactions*numberSiteTypes; i++)
    {
    	tmpRL[i] = RL[i];
    }

    #ifdef MPD_GLOBAL_T_MATRIX
	CUDA_EXCEPTION_CHECK(cudaMalloc(&TG, DFmatrixSize*sizeof(float)));
	CUDA_EXCEPTION_CHECK(cudaMemcpy(TG, T, DFmatrixSize*sizeof(float), cudaMemcpyHostToDevice));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(TC, &TG, sizeof(float*)));
    #else
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(TC, T, DFmatrixSize*sizeof(float)));
    #endif

    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberSpeciesC, &numberSpecies, sizeof(numberSpeciesC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(numberSiteTypesC, &numberSiteTypes, sizeof(numberSiteTypesC)));
    const unsigned int latticeXYSize = latticeXSize*latticeYSize;
    const unsigned int latticeXYZSize = latticeXSize*latticeYSize*latticeZSize;
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXSizeC, &latticeXSize, sizeof(latticeYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeYSizeC, &latticeYSize, sizeof(latticeYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeZSizeC, &latticeZSize, sizeof(latticeZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(global_latticeZSizeC, &latticeZSize, sizeof(latticeZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXYSizeC, &latticeXYSize, sizeof(latticeXYSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(latticeXYZSizeC, &latticeXYZSize, sizeof(latticeXYZSizeC)));
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(global_latticeXYZSizeC, &latticeXYZSize, sizeof(latticeXYZSizeC)));
    #ifdef MPD_GLOBAL_S_MATRIX
	// Store RL in global memory too, since I'm going to assume if S is too big, then RL is too.
	cudaMalloc(&RLG, numberReactions*numberSiteTypes * sizeof(uint8_t));
	cudaMemcpy(RLG, tmpRL, numberReactions*numberSiteTypes * sizeof(uint8_t), cudaMemcpyHostToDevice);
    #else
	// RL is stored in constant memory
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(RLC, tmpRL, numberReactions*numberSiteTypes*sizeof(uint8_t)));
    #endif
    delete [] tmpRL;
    // site type t is live for diffusion when some species has a nonzero transition probability into it from any type
    wcm_diffLiveTypes = 0;
    for (uint t=0; t<numberSiteTypes && t<32; t++)
        for (uint src=0; src<numberSiteTypes; src++)
            for (uint sp=0; sp<numberSpecies; sp++)
                if (T[src*numberSiteTypes*numberSpecies + t*numberSpecies + sp] != 0.0f) wcm_diffLiveTypes |= (1u << t);
    delete [] T;

    // Set the cuda reaction model rates now that we have the subvolume size.
    model_reactionRates = new float[numberReactions];

    // Set up pre-configured propensity matrices
	zeroOrderSize=numberSiteTypes;
	firstOrderSize=numberSpecies*numberSiteTypes;
	secondOrderSize=numberSpecies*numberSpecies*numberSiteTypes;
	zeroOrder=new float[zeroOrderSize];
	firstOrder=new float[firstOrderSize];
	secondOrder=new float[secondOrderSize];

    #ifdef MPD_GLOBAL_R_MATRIX
    cudaMalloc(&reactionRatesG, numberReactions*sizeof(float));
    #endif
		
	cudaMalloc(&propZeroOrder, zeroOrderSize*sizeof(float));
	cudaMalloc(&propFirstOrder, firstOrderSize*sizeof(float));
	cudaMalloc(&propSecondOrder, secondOrderSize*sizeof(float));

    computePropensities();
    copyModelsToDevice();
	reactionModelModified = false;
}

void IntMpdRdmeSolver::setLatticeData(const uint8_t* latticeData)
{
	lattice_coord_t s = lattice->getSize();
	site_size_t p = lattice->getMaxOccupancy();

	// reinterpret the input as a 32-bit array
	uint32_t *data = (uint32_t*)latticeData;

	// Set the lattice data.
	for (uint i=0, index=0; i<s.x; i++)
	{
		for (uint j=0; j<s.y; j++)
		{
			for (uint k=0; k<s.z; k++)
			{
				for (uint l=0; l<p; l++, index++)
				{
					if (data[index] != 0)
					{
						if (data[index] > numberSpecies) throw InvalidArgException("latticeData", "an invalid species was found",data[index]);
						lattice->addParticle(i,j,k,data[index]);
					}
				}
			}
		}
	}
}

void IntMpdRdmeSolver::setReactionRate(unsigned int rxid, float rate)
{
    if (reactionTypes[rxid] ==  ZerothOrderPropensityArgs::REACTION_TYPE)
	{
		((ZerothOrderPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
	else if (reactionTypes[rxid] == FirstOrderPropensityArgs::REACTION_TYPE)
	{
		((FirstOrderPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
	else if (reactionTypes[rxid] == SecondOrderPropensityArgs::REACTION_TYPE)
	{
		((SecondOrderPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
	else if (reactionTypes[rxid] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
	{
		((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[rxid])->k = rate;
	}
    else
    {
    	throw InvalidArgException("reactionTypeA", "the reaction type was not supported by the solver", reactionTypes[rxid]);
    }

	// Flag to trigger re-computation of rdme propensities
	reactionModelModified = true;
}

void IntMpdRdmeSolver::computePropensities()
{
	unsigned int latticeXSize = lattice->getXSize();
	unsigned int latticeYSize = lattice->getYSize();
	unsigned int latticeZSize = lattice->getZSize();
    for (uint i=0; i<numberReactions; i++)
    {
    	if (reactionTypes[i] == ZerothOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((ZerothOrderPropensityArgs *)propensityFunctionArgs[i])->k*tau/(latticeXSize*latticeYSize*latticeZSize);
    	}
    	else if (reactionTypes[i] == FirstOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((FirstOrderPropensityArgs *)propensityFunctionArgs[i])->k*tau;
    	}
    	else if (reactionTypes[i] == SecondOrderPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->k*tau*latticeXSize*latticeYSize*latticeZSize;
    	}
    	else if (reactionTypes[i] == SecondOrderSelfPropensityArgs::REACTION_TYPE)
    	{
    		model_reactionRates[i] = ((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i])->k*tau*latticeXSize*latticeYSize*latticeZSize;
    	}
    	else
    	{
    		throw InvalidArgException("reactionTypeA", "the reaction type was not supported by the solver", reactionTypes[i]);
    	}
    }

    // species that are a reactant of a first- or second-order reaction (the only ones with nonzero qp1/qp2 terms)
    {
        std::vector<uint8_t> reactive(numberSpecies, 0);
        for (uint i=0; i<numberReactions; i++)
        {
            if (reactionTypes[i] == FirstOrderPropensityArgs::REACTION_TYPE) reactive[((FirstOrderPropensityArgs *)propensityFunctionArgs[i])->si] = 1;
            else if (reactionTypes[i] == SecondOrderPropensityArgs::REACTION_TYPE) { reactive[((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s1i] = 1; reactive[((SecondOrderPropensityArgs *)propensityFunctionArgs[i])->s2i] = 1; }
            else if (reactionTypes[i] == SecondOrderSelfPropensityArgs::REACTION_TYPE) reactive[((SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i])->si] = 1;
        }
        if (wcm_reactiveG == NULL) CUDA_EXCEPTION_CHECK(cudaMalloc(&wcm_reactiveG, numberSpecies > 0 ? numberSpecies : 1));
        if (numberSpecies > 0) CUDA_EXCEPTION_CHECK(cudaMemcpy(wcm_reactiveG, reactive.data(), numberSpecies, cudaMemcpyHostToDevice));
    }
    float scale=latticeXSize*latticeYSize*latticeZSize;
	for(uint site=0; site<numberSiteTypes; site++)
	{
		uint o1=site*numberSpecies;
		uint o2=site*numberSpecies*numberSpecies;

		zeroOrder[site]=0.0f;
		for (uint i=0; i<numberSpecies; i++)
		{
			firstOrder[o1 + i]=0.0f;
			for (uint j=0; j<numberSpecies; j++)
				secondOrder[o2 + i*numberSpecies + j]=0.0f;	
		}

		for (uint i=0; i<numberReactions; i++)
		{
			if(! RL[i*numberSiteTypes + site])
				continue;

			switch(reactionTypes[i])
			{
				case ZerothOrderPropensityArgs::REACTION_TYPE:
				{
					ZerothOrderPropensityArgs *rx=(ZerothOrderPropensityArgs *)propensityFunctionArgs[i];
					zeroOrder[site]+=rx->k*tau/scale;
				}   break;

				case FirstOrderPropensityArgs::REACTION_TYPE:
				{
				    FirstOrderPropensityArgs *rx=(FirstOrderPropensityArgs *)propensityFunctionArgs[i];
				    firstOrder[o1 + rx->si]+=rx->k*tau;
				}   break;
				
				case SecondOrderPropensityArgs::REACTION_TYPE:
				{
				    SecondOrderPropensityArgs *rx=(SecondOrderPropensityArgs *)propensityFunctionArgs[i];
				    secondOrder[o2 + (rx->s1i) * numberSpecies + (rx->s2i)]=rx->k*tau*scale;
				    secondOrder[o2 + (rx->s2i) * numberSpecies + (rx->s1i)]=rx->k*tau*scale;
				}   break; 

				case SecondOrderSelfPropensityArgs::REACTION_TYPE:
				{
				    SecondOrderSelfPropensityArgs *rx=(SecondOrderSelfPropensityArgs *)propensityFunctionArgs[i];
				    secondOrder[o2 + (rx->si) * numberSpecies + (rx->si)]+=rx->k*tau*scale*2;
				}
			}
		}
	}
}

void IntMpdRdmeSolver::copyModelsToDevice()
{
#ifdef MPD_GLOBAL_R_MATRIX
    cudaMemcpy(reactionRatesG, model_reactionRates, numberReactions*sizeof(float), cudaMemcpyHostToDevice);
#else
    CUDA_EXCEPTION_CHECK(cudaMemcpyToSymbol(reactionRatesC, model_reactionRates, numberReactions*sizeof(float)));
#endif

	cudaMemcpy(propZeroOrder, zeroOrder, zeroOrderSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(propFirstOrder, firstOrder, firstOrderSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(propSecondOrder, secondOrder, secondOrderSize*sizeof(float), cudaMemcpyHostToDevice);
}

void IntMpdRdmeSolver::generateTrajectory()
{
    // Shadow the lattice member as a cuda lattice.
    CudaIntLattice * lattice = (CudaIntLattice *)this->lattice;

    // Get the interval for writing species counts and lattices.
    uint32_t speciesCountsWriteInterval=atol((*parameters)["writeInterval"].c_str());
    uint32_t nextSpeciesCountsWriteTime = speciesCountsWriteInterval;
    lm::io::SpeciesCounts speciesCountsDataSet;
    speciesCountsDataSet.set_number_species(numberSpeciesToTrack);
    speciesCountsDataSet.set_number_entries(0);
    uint32_t latticeWriteInterval=atol((*parameters)["latticeWriteInterval"].c_str());
    uint32_t nextLatticeWriteTime = latticeWriteInterval;
    lm::io::Lattice latticeDataSet;

    // Get the simulation time limit.
    double maxTime=atof((*parameters)["maxTime"].c_str());

    Print::printf(Print::INFO,
                  "Running mpd rdme simulation with %d species, %d reactions, %d site types for %e s with tau %e. Writing species at %e and lattice at %e intervals",
                  numberSpecies, numberReactions, numberSiteTypes,
                  maxTime, tau,
                  speciesCountsWriteInterval, latticeWriteInterval);

    // Set the initial time.
	double time = 0.0;
	uint32_t current_timestep=0;

    bool hookEnabled=false;
    uint32_t nextHookTime=0;
    uint32_t hookInterval=0;

    // Find out at what interval to hook simulations
    if((*parameters)["hookInterval"] != "")
    {
        hookInterval=atol((*parameters)["hookInterval"].c_str());
        hookEnabled=true;
        nextHookTime=hookInterval;
    }

    // Find out at what interval to write status messages
    double printPerfInterval = 60;  
    if((*parameters)["perfPrintInterval"] != "")
    {
        printPerfInterval=atof((*parameters)["perfPrintInterval"].c_str());
    }
         

    // Call beginning hook. Any modifications to the lattice will be
    // accounted for because we have not yet copied to device memory
    // thus we do not need to check return value
    // onBeginTrajectory(lattice);

    // Synchronize the cuda memory.
    lattice->copyToGPU(); wcm_recount_occupancy(lattice, cudaStream);

    // Record the initial species counts.
    recordSpeciesCounts(time, lattice, &speciesCountsDataSet);

    // Write the initial lattice.
    writeLatticeData(time, lattice, &latticeDataSet);

    // initialize max counts from initial particle/site lattice
    // initMaxCounts(lattice);

    Timer timer;
    timer.tick();
    double lastT     = 0;
    int    lastSteps = 0;

    // Perform an initial hook check
    hookCheckSimulation(time, lattice);

    // Loop until we have finished the simulation
    uint32_t wcm_lastK = 1;   // timesteps run by the previous iteration
    while (time < maxTime)
    {

		if(globalAbort)
		{
			printf("Global abort: terminating solver\n");
			break;
		}

        lastT     += timer.tock();
        lastSteps += wcm_lastK;

        if (lastT >= printPerfInterval)
        {
            double stepTime = lastT/lastSteps;
            double completionTime = stepTime*(maxTime-time)/tau;
            std::string units;
            if (completionTime > 60*60*24*365) {
                units = "weeks";
                completionTime /= 60*60*24*365;
            } else if (completionTime > 60*60*24*30) {
                units = "months";
                completionTime /= 60*60*24*30;
            } else if (completionTime > 60*60*24*7) {
                units = "weeks";
                completionTime /= 60*60*24*7;
            } else if (completionTime > 60*60*24) {
                units = "days";
                completionTime /= 60*60*24;
            } else if (completionTime > 60*60) {
                units = "hours";
                completionTime /= 60*60;
            } else if (completionTime > 60) {
                units = "minutes";
                completionTime /= 60;
            } else {
                units = "seconds";
            }

            Print::printf(Print::INFO, "Average walltime per timestep: %.2f ms. Progress: %.4fs/%.4fs (% .3g%% done / %.2g %s walltime remaining)",
                                       1000.0*stepTime, time, maxTime, 100.0*time/maxTime, completionTime, units.c_str());

            lastT     = 0;
            lastSteps = 0;
        }

        timer.tick();

        // Run the next timestep(s). up to wcm_batch_steps() per host synchronisation, never past the next lattice
        // write, species-count write or hook step (so every check below fires on the same step as before) nor past maxTime.
        uint32_t K = wcm_batch_steps();
        if (nextLatticeWriteTime > current_timestep) K = std::min<uint32_t>(K, nextLatticeWriteTime - current_timestep); else K = 1;
        if (nextSpeciesCountsWriteTime > current_timestep) K = std::min<uint32_t>(K, nextSpeciesCountsWriteTime - current_timestep); else K = 1;
        if (hookEnabled) { if (nextHookTime > current_timestep) K = std::min<uint32_t>(K, nextHookTime - current_timestep); else K = 1; }
        for (uint32_t j = 1; j < K; j++) if ((current_timestep + j)*tau >= maxTime) { K = j; break; }
        wcm_runTimesteps(lattice, current_timestep, K);
        current_timestep += K;
        wcm_lastK = K;

        // Update the time.
	    time = current_timestep*tau;

        // See if we need to write out the any data.
        if (current_timestep >= nextLatticeWriteTime
            || current_timestep >= nextSpeciesCountsWriteTime
            || (hookEnabled && current_timestep >= nextHookTime))
        {
            reactionModelModified = false;

            // Synchronize the lattice. a step that only runs the hook (no lattice or species-count output due) takes
            // the site types now and the particles only when the hook first reads them (WCM_LAZY_DOWNLOAD_OFF=1: always both now).
            static const bool offLazy = std::getenv("WCM_LAZY_DOWNLOAD_OFF") != NULL;
            const bool hook_only = hookEnabled && current_timestep >= nextHookTime &&
                                   current_timestep < nextLatticeWriteTime && current_timestep < nextSpeciesCountsWriteTime;
            if (hook_only && !offLazy) lattice->copyFromGPULazy();
            else lattice->copyFromGPU();

	        Print::printf(Print::INFO, "Time is %.14f", time);

	        // Check if we need to execute the hook
            if (hookEnabled && current_timestep >= nextHookTime)
            {
		        Print::printf(Print::INFO, "Hook time is %.14f, in steps is %d", time, nextHookTime);
   	            hookCheckSimulation(time, lattice);
                nextHookTime += hookInterval;
		        Print::printf(Print::INFO, "Next hook time is %d", nextHookTime);
            }

            // See if we need to write the lattice.
	        if (current_timestep >= nextLatticeWriteTime)
            {
                PROF_BEGIN(PROF_SERIALIZE_LATTICE);
		        Print::printf(Print::INFO, "Lattice write time is %.14f, in steps is %.d", time, nextLatticeWriteTime);
                writeLatticeData(time, lattice, &latticeDataSet);
                nextLatticeWriteTime += latticeWriteInterval;
		        Print::printf(Print::INFO, "Next lattice write time is %.d", nextLatticeWriteTime);
                PROF_END(PROF_SERIALIZE_LATTICE);
            }

            // See if we need to write the species counts.
            if (current_timestep >= nextSpeciesCountsWriteTime)
            {
                PROF_BEGIN(PROF_DETERMINE_COUNTS);
                recordSpeciesCounts(time, lattice, &speciesCountsDataSet);
                nextSpeciesCountsWriteTime += speciesCountsWriteInterval;
                PROF_END(PROF_DETERMINE_COUNTS);

                // See if we have accumulated enough species counts to send.
                if (speciesCountsDataSet.number_entries() >= TUNE_SPECIES_COUNTS_BUFFER_SIZE)
                {
                    PROF_BEGIN(PROF_SERIALIZE_COUNTS);
                    writeSpeciesCounts(&speciesCountsDataSet);
                    PROF_END(PROF_SERIALIZE_COUNTS);
                }
            }
        }
    }

    lattice->copyFromGPU();

    // Write any remaining species counts.
    writeSpeciesCounts(&speciesCountsDataSet);
}

int IntMpdRdmeSolver::hookSimulation(double time, CudaIntLattice *lattice)
{                                                       
    // Overload this function in derivative classes
	// Return 0 if the lattice state is unchanged
	// Return 1 if the lattice state has been modified, 
    //          and it needs to be copied back to the GPU
	// Return 2 if lattice sites have changed, and should
	//			be copied back to the GPU *and* be recorded
	//			in the output file           
    return 0;                                             
}

void IntMpdRdmeSolver::writeLatticeData(double time, CudaIntLattice * lattice, lm::io::Lattice * latticeDataSet)
{
    Print::printf(Print::DEBUG, "Writing lattice at %e s", time);

    // Record the lattice data.
    latticeDataSet->Clear();
    latticeDataSet->set_lattice_x_size(lattice->getSize().x);
    latticeDataSet->set_lattice_y_size(lattice->getSize().y);
    latticeDataSet->set_lattice_z_size(lattice->getSize().z);
    latticeDataSet->set_particles_per_site(lattice->getMaxOccupancy());
    latticeDataSet->set_time(time);

    // Push it to the output queue.
    size_t payloadSize = lattice->getSize().x*lattice->getSize().y*lattice->getSize().z*lattice->getMaxOccupancy()*sizeof(uint32_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::INT_LATTICE, replicate, latticeDataSet, lattice, payloadSize, &lm::rdme::IntLattice::nativeSerialize);
}

void IntMpdRdmeSolver::writeLatticeSites(double time, CudaIntLattice * lattice)
{
    Print::printf(Print::DEBUG, "Writing lattice sites at %e s", time);

	lm::io::Lattice latticeDataSet;
    // Record the lattice data.
    latticeDataSet.Clear();
    latticeDataSet.set_lattice_x_size(lattice->getSize().x);
    latticeDataSet.set_lattice_y_size(lattice->getSize().y);
    latticeDataSet.set_lattice_z_size(lattice->getSize().z);
    latticeDataSet.set_particles_per_site(lattice->getMaxOccupancy());
    latticeDataSet.set_time(time);

    // Push it to the output queue.
    size_t payloadSize = lattice->getSize().x*lattice->getSize().y*lattice->getSize().z*sizeof(uint8_t);
    lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SITE_LATTICE, replicate, &latticeDataSet, lattice, payloadSize, &lm::rdme::ByteLattice::nativeSerializeSites);
}

void IntMpdRdmeSolver::recordSpeciesCounts(double time, CudaIntLattice * lattice, lm::io::SpeciesCounts * speciesCountsDataSet)
{
    std::map<particle_t,uint> particleCounts = lattice->getParticleCounts();
    speciesCountsDataSet->set_number_entries(speciesCountsDataSet->number_entries()+1);
    speciesCountsDataSet->add_time(time);
    for (particle_t p=0; p<numberSpeciesToTrack; p++)
    {
        speciesCountsDataSet->add_species_count((particleCounts.count(p+1)>0)?particleCounts[p+1]:0);
    }
}

void IntMpdRdmeSolver::writeSpeciesCounts(lm::io::SpeciesCounts * speciesCountsDataSet)
{
    if (speciesCountsDataSet->number_entries() > 0)
    {
        // Push it to the output queue.
        lm::main::DataOutputQueue::getInstance()->pushDataSet(lm::main::DataOutputQueue::SPECIES_COUNTS, replicate, speciesCountsDataSet);

        // Reset the data set.
        speciesCountsDataSet->Clear();
        speciesCountsDataSet->set_number_species(numberSpeciesToTrack);
        speciesCountsDataSet->set_number_entries(0);
    }
}

void IntMpdRdmeSolver::hookCheckSimulation(double time, CudaIntLattice * lattice)
{
    switch(hookSimulation(time, lattice))
    {
        case 0:
            break; 

        // particles the hook never read are still the device's (wcmCopyToGPU uploads the site types only); the
        // occupancy is recounted as before whenever the host particles were current
        case 1:
            lattice->wcmCopyToGPU(); if (!lattice->wcmHostParticlesStale()) wcm_recount_occupancy(lattice, cudaStream);
            break;

        case 2:
            lattice->wcmCopyToGPU(); if (!lattice->wcmHostParticlesStale()) wcm_recount_occupancy(lattice, cudaStream);
            writeLatticeSites(time, lattice);
            break;

        case 4:   // the hook changed site types only (particles and hence occupancy unchanged)
            lattice->copySiteTypesToGPU();
            break;

        case 3:
            printf("hook return value is 3, force to stop.\n");
            writeLatticeSites(time, lattice);
            return;

        default:
            throw("Unknown hook return value");
    }

    if(reactionModelModified)
		computePropensities();
}

uint64_t IntMpdRdmeSolver::getTimestepSeed(uint32_t timestep, uint32_t substep)
{
    uint64_t timestepHash = (((((uint64_t)seed)<<30)+timestep)<<2)+substep;
    timestepHash = timestepHash * 3202034522624059733ULL + 4354685564936845319ULL;
    timestepHash ^= timestepHash >> 20; timestepHash ^= timestepHash << 41; timestepHash ^= timestepHash >> 5;
    timestepHash *= 7664345821815920749ULL;
    return timestepHash;
}

// run K timesteps with one host synchronisation.
// Each timestep used to end with cudaStreamSynchronize + a host check of the overflow list (~8 us of idle GPU per ~117 us step
// on the capture lattice; overflows are rare: <= 6 in 20,000 steps). The K steps are queued back to back after a device copy of
// the source particle buffer; if the overflow list is empty after the batch the result is exactly what K single steps give
// (same kernels, same per-timestep seeds, and without overflows the per-step host path does nothing). Otherwise the lattice is
// restored from the copy (occupancy recounted as after any copyToGPU), the list cleared, and the K steps rerun one at a time
// through runTimestep, which handles each overflow on the host as before. WCM_RDME_BATCH=<n> sets K (default 32, 1 = off).

void IntMpdRdmeSolver::wcm_runTimesteps(CudaIntLattice * lattice, uint32_t first, uint32_t K)
{
    if (K <= 1) { runTimestep(lattice, first); return; }
    wcm_update_box(lattice, cudaStream, zeroOrder, zeroOrderSize);
    const size_t bytes = (size_t)lattice->getNumberSites() * lattice->getMaxOccupancy() * sizeof(uint32_t);
    if (wcm_snap_bytes < bytes)
    {
        if (wcm_snap != NULL) cudaFree(wcm_snap);
        CUDA_EXCEPTION_CHECK(cudaMalloc(&wcm_snap, bytes));
        wcm_snap_bytes = bytes;
    }
    void * src0 = lattice->getGPUMemorySrc();
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(wcm_snap, src0, bytes, cudaMemcpyDeviceToDevice, cudaStream));
    for (uint32_t k = 0; k < K; k++) wcm_launchTimestep(lattice, first + k);
    wcm_queue_box_copy(cudaStream);
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(cudaStream));
#ifndef MPD_MAPPED_OVERFLOWS
    uint32_t numberExceptions = 0;
    CUDA_EXCEPTION_CHECK(cudaMemcpy(&numberExceptions, cudaOverflowList, sizeof(uint32_t), cudaMemcpyDeviceToHost));
#else
    uint32_t numberExceptions = ((uint32_t*)cudaOverflowList)[0];
#endif
    if (numberExceptions == 0)
    {
        // what the K per-step overflow checks would have counted (no uses)
        for (uint32_t k = 0; k < K; k++)
        {
            overflowTimesteps++;
            if (overflowTimesteps >= 1000)
            {
                if (overflowListUses > 10)
                    Print::printf(Print::WARNING, "%d uses of the particle overflow list in the last 1000 timesteps, performance may be degraded.", overflowListUses);
                overflowTimesteps = 0;
                overflowListUses = 0;
            }
        }
        return;
    }
    // an overflow happened somewhere in the batch: back to the batch start, then one step at a time
    Print::printf(Print::DEBUG, "overflow in batch at timestep %u..%u, rerunning step by step", first, first + K - 1);
    if (lattice->getGPUMemorySrc() != src0) lattice->swapSrcDest();
    CUDA_EXCEPTION_CHECK(cudaMemcpyAsync(src0, wcm_snap, bytes, cudaMemcpyDeviceToDevice, cudaStream));
    wcm_recount_occupancy(lattice, cudaStream);
    CUDA_EXCEPTION_CHECK(cudaMemsetAsync(cudaOverflowList, 0, MPD_OVERFLOW_LIST_SIZE, cudaStream));
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(cudaStream));
    for (uint32_t k = 0; k < K; k++) runTimestep(lattice, first + k);
}

// one fused x/y/z diffusion kernel per timestep (wcm_mpd_xyz_kernel) instead of mpd_x/y/z_kernel, same lattice.
// WCM_FUSE_OFF=1 keeps the three kernels; WCM_FUSE_TILE=<v> fixes the tile core (0: 24x6x6 / 1024 threads; 1: 32x8x4 / 1024;
// 2: 16x4x4 / 256; 3: 17x7x7 / 1024); by default (or =auto) the solver picks 0, 3 or 2 per batch from the box (wcm_pick_tile). The path is used only when every lattice dimension is at least the tile plus its apron
// and particle values fit in 16 bits (numberSpecies < 65535); otherwise the three kernels run.
// the fused kernel also takes the reaction check of precomp_reaction_kernel (propensity sums over a block-wide list of
// the nonzero terms, checkForReaction); the few sites that react are listed and finished by wcm_rxn_tail_kernel.
// WCM_FUSE_RXN_OFF=1: the separate precomp_reaction_kernel as in 126.
namespace intmpdrdme_dev {
namespace wcm_fuse {
template<int CX, int CY, int CZ> struct Tile
{
    static constexpr int IX = CX+2, IY = CY+2, IZ = CZ+2, NI = IX*IY*IZ, NCORE = CX*CY*CZ;
    static_assert(NI < 4096, "work items hold the site index in 12 bits");
    static_assert(NI < (1 << 24), "reaction term descriptors hold the site index in 24 bits");
    static constexpr int CAP = 2*NI;   // particle work list of a pass (a larger pass falls back to the per-site loop)
    // dynamic shared memory: uint16 A[16*NI], uint16 B[16*NI], uint32 ch[NI], uint32 li[NI], int wsum[64], uint16 items[CAP],
    // uint8 oA[NI], oB[NI], st[NI]
    static constexpr size_t SMEM = (size_t)NI*(2*MPD_WORDS_PER_SITE*2 + 4 + 4 + 3) + 64*4 + (size_t)CAP*2;
};
// what precomp_reaction_kernel receives besides the lattice (on == false: no reaction step in the fused kernel)
struct Rxn
{
    bool on; unsigned long long hash;
    const int8_t * SG; const uint8_t * RLG; const unsigned int * reactionOrdersG; const unsigned int * reactionSitesG;
    const unsigned int * D1G; const unsigned int * D2G; const float * reactionRatesG; const float * qp0; const float * qp1; const float * qp2;
    uint2 * list; unsigned int * ctr;
    const uint8_t * reactive;   // per species: 1 if it is a reactant of some first- or second-order reaction (else all its qp1/qp2 terms are 0)   // reacting sites (lattice index, total propensity bits); ctr[0] = count, ctr[1] = tail blocks done
};
}
template<int CX, int CY, int CZ, int NT, int MINB, bool RXN>
__global__ void __launch_bounds__(NT, MINB) wcm_mpd_xyz_kernel(const unsigned int* __restrict__ inLattice, const uint8_t * __restrict__ inSites, unsigned int* __restrict__ outLattice, const unsigned long long hashX, const unsigned long long hashY, const unsigned long long hashZ, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc, const bool wcmCompact, const wcm_fuse::Rxn rxn);
#ifdef WCM_FUSE_RXN_AVAILABLE
__global__ void __launch_bounds__(128) wcm_rxn_tail_kernel(unsigned int* lattice, const uint8_t * __restrict__ inSites, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, uint8_t* __restrict__ occ, const uint2* __restrict__ list, unsigned int* ctr);
#endif
}
static int wcm_fuse_variant()
{
    static int v = -2;
    if (v == -2)
    {
        const char * off = getenv("WCM_FUSE_OFF");
        const char * t = getenv("WCM_FUSE_TILE");
        v = (off != NULL && atoi(off) != 0) ? -1 : ((t != NULL && strcmp(t, "auto") != 0) ? atoi(t) : 100);
        lm::Print::printf(lm::Print::INFO, "fused x/y/z diffusion kernel %s (tile variant %d%s)", v >= 0 ? "on" : "off", v, v == 100 ? " = auto" : "");
    }
    return (v == 100) ? wcm_fuse_auto : v;
}
static bool wcm_fuse_compact()
{
    static const bool v = getenv("WCM_FUSE_NOCOMPACT") == NULL;   // =1: per-site choice loop only (same results; for verification)
    return v;
}
template<int CX, int CY, int CZ, int NT, int MINB>
static bool wcm_launch_fused(lm::rdme::CudaIntLattice * lattice, cudaStream_t stream, uint64_t hx, uint64_t hy, uint64_t hz, void * overflowList, unsigned int numberSpecies, const intmpdrdme_dev::wcm_fuse::Rxn & rxn)
{
    const auto size = lattice->getSize();
    const unsigned int X = size.x, Y = size.y, Z = size.z;
    if (X < CX+2 || Y < CY+2 || Z < CZ+2 || numberSpecies >= 65535u || MPD_WORDS_PER_SITE > 16) return false;
    constexpr size_t smem = intmpdrdme_dev::wcm_fuse::Tile<CX,CY,CZ>::SMEM;
    static bool attr = false;
    if (!attr)
    {
        CUDA_EXCEPTION_CHECK(cudaFuncSetAttribute(intmpdrdme_dev::wcm_mpd_xyz_kernel<CX,CY,CZ,NT,MINB,false>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smem));
        CUDA_EXCEPTION_CHECK(cudaFuncSetAttribute(intmpdrdme_dev::wcm_mpd_xyz_kernel<CX,CY,CZ,NT,MINB,true>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smem));
        attr = true;
    }
    dim3 grid((X+CX-1)/CX, (Y+CY-1)/CY, (Z+CZ-1)/CZ), block(NT, 1, 1);
    if (rxn.on)
    {    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::wcm_mpd_xyz_kernel<CX,CY,CZ,NT,MINB,true><<<grid, block, smem, stream>>>((const unsigned int *)lattice->getGPUMemorySrc(), (const uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), hx, hy, hz, (unsigned int *)overflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()), wcm_fuse_compact(), rxn))); }
    else
    {    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::wcm_mpd_xyz_kernel<CX,CY,CZ,NT,MINB,false><<<grid, block, smem, stream>>>((const unsigned int *)lattice->getGPUMemorySrc(), (const uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), hx, hy, hz, (unsigned int *)overflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()), wcm_fuse_compact(), rxn))); }
    return true;
}

// resident blocks per SM of a fused-kernel instance (the reaction instance, the one normally launched)
template<int CX, int CY, int CZ, int NT, int MINB>
long wcm_resident_tiles()
{
    constexpr size_t smem = intmpdrdme_dev::wcm_fuse::Tile<CX,CY,CZ>::SMEM;
    CUDA_EXCEPTION_CHECK(cudaFuncSetAttribute(intmpdrdme_dev::wcm_mpd_xyz_kernel<CX,CY,CZ,NT,MINB,true>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smem));
    int n = 0;
    CUDA_EXCEPTION_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n, intmpdrdme_dev::wcm_mpd_xyz_kernel<CX,CY,CZ,NT,MINB,true>, NT, smem));
    return n > 0 ? n : 1;
}

// the kernel launches of one timestep (x, y, z diffusion + reaction), without the synchronisation and overflow
// handling, so several timesteps can be queued back to back (wcm_runTimesteps).
void IntMpdRdmeSolver::wcm_launchTimestep(CudaIntLattice * lattice, uint32_t timestep)
{

    if(reactionModelModified)
        copyModelsToDevice();

    // Calculate some properties of the lattice.
    lattice_coord_t size = lattice->getSize();
    const unsigned int latticeXSize = size.x;
    const unsigned int latticeYSize = size.y;
    const unsigned int latticeZSize = size.z;

    dim3 gridSize, threadBlockSize;

    // Execute the kernel for the x direction.
    PROF_CUDA_START(cudaStream);

    // the three diffusion passes in one kernel when the tile fits the lattice
    bool wcm_fused = false;
    // the reaction step inside the fused kernel too (WCM_FUSE_RXN_OFF=1: separate precomp_reaction_kernel as in 126)
    intmpdrdme_dev::wcm_fuse::Rxn rxn = {};
    #ifdef WCM_FUSE_RXN_AVAILABLE
    {
        static const bool rxnOff = getenv("WCM_FUSE_RXN_OFF") != NULL && atoi(getenv("WCM_FUSE_RXN_OFF")) != 0;
        static bool said = false;
        if (!said) { Print::printf(Print::INFO, "reaction step in the fused kernel %s", rxnOff ? "off" : "on"); said = true; }
        if (!rxnOff && numberReactions > 0)
        {
            rxn.on = true; rxn.hash = getTimestepSeed(timestep,3);
            rxn.SG = SG; rxn.RLG = RLG; rxn.reactionOrdersG = reactionOrdersG; rxn.reactionSitesG = reactionSitesG;
            rxn.D1G = D1G; rxn.D2G = D2G; rxn.reactionRatesG = reactionRatesG;
            rxn.qp0 = propZeroOrder; rxn.qp1 = propFirstOrder; rxn.qp2 = propSecondOrder;
            static uint2 * list = NULL; static unsigned int * ctr = NULL; static size_t listN = 0;
            const size_t nSites = lattice->getNumberSites();
            if (listN < nSites)
            {
                if (list != NULL) { cudaFree(list); cudaFree(ctr); }
                CUDA_EXCEPTION_CHECK(cudaMalloc(&list, nSites*sizeof(uint2)));
                CUDA_EXCEPTION_CHECK(cudaMalloc(&ctr, 2*sizeof(unsigned int)));
                CUDA_EXCEPTION_CHECK(cudaMemset(ctr, 0, 2*sizeof(unsigned int)));
                listN = nSites;
            }
            rxn.list = list; rxn.ctr = ctr; rxn.reactive = wcm_reactiveG;
        }
    }
    #endif
    {
        const uint64_t hx = getTimestepSeed(timestep,0), hy = getTimestepSeed(timestep,1), hz = getTimestepSeed(timestep,2);
        switch (wcm_fuse_variant())
        {
        case 0: wcm_fused = wcm_launch_fused<24,6,6,1024,1>(lattice, cudaStream, hx, hy, hz, cudaOverflowList, numberSpecies, rxn); break;
        case 1: wcm_fused = wcm_launch_fused<32,8,4,1024,1>(lattice, cudaStream, hx, hy, hz, cudaOverflowList, numberSpecies, rxn); break;
        case 2: wcm_fused = wcm_launch_fused<16,4,4,256,4>(lattice, cudaStream, hx, hy, hz, cudaOverflowList, numberSpecies, rxn); break;
        case 3: wcm_fused = wcm_launch_fused<17,7,7,1024,1>(lattice, cudaStream, hx, hy, hz, cudaOverflowList, numberSpecies, rxn); break;
        default: break;
        }
    }
    if (wcm_fused)
        lattice->swapSrcDest();
    else
    {
    PROF_CUDA_BEGIN(PROF_MPD_X_DIFFUSION,cudaStream);
    #ifdef MPD_CUDA_3D_GRID_LAUNCH
    calculateXLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::mpd_x_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,0), (unsigned int*)cudaOverflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()))));
    #else
    unsigned int gridXSize;
    calculateXLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_X_BLOCK_MAX_X_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::mpd_x_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), gridXSize, getTimestepSeed(timestep,0), (unsigned int*)cudaOverflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()))));
    #endif
    PROF_CUDA_END(PROF_MPD_X_DIFFUSION,cudaStream);
    lattice->swapSrcDest();

    // Execute the kernel for the y direction.
    PROF_CUDA_BEGIN(PROF_MPD_Y_DIFFUSION,cudaStream);
    #ifdef MPD_CUDA_3D_GRID_LAUNCH
    calculateYLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_Y_BLOCK_X_SIZE, TUNE_MPD_Y_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::mpd_y_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,1), (unsigned int*)cudaOverflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()))));
    #else
    calculateYLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_Y_BLOCK_X_SIZE, TUNE_MPD_Y_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::mpd_y_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), gridXSize, getTimestepSeed(timestep,1), (unsigned int*)cudaOverflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()))));
    #endif
    PROF_CUDA_END(PROF_MPD_Y_DIFFUSION,cudaStream);
    lattice->swapSrcDest();

    // Execute the kernel for the z direction.
    PROF_CUDA_BEGIN(PROF_MPD_Z_DIFFUSION,cudaStream);
    #ifdef MPD_CUDA_3D_GRID_LAUNCH
    calculateZLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_Z_BLOCK_X_SIZE, TUNE_MPD_Z_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::mpd_z_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), getTimestepSeed(timestep,2), (unsigned int*)cudaOverflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()))));
    #else
    calculateZLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_Z_BLOCK_X_SIZE, TUNE_MPD_Z_BLOCK_Z_SIZE, latticeXSize, latticeYSize, latticeZSize);
    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::mpd_z_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemoryDest(), gridXSize, getTimestepSeed(timestep,2), (unsigned int*)cudaOverflowList, wcm_occ(lattice, lattice->getGPUMemorySrc()), wcm_occ(lattice, lattice->getGPUMemoryDest()))));
    #endif
    PROF_CUDA_END(PROF_MPD_Z_DIFFUSION,cudaStream);
    lattice->swapSrcDest();
    }

    #ifdef WCM_FUSE_RXN_AVAILABLE
    if (wcm_fused && rxn.on)
        CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::wcm_rxn_tail_kernel<<<64, 128, 0, cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), rxn.hash, (unsigned int*)cudaOverflowList, SG, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG, wcm_occ(lattice, lattice->getGPUMemorySrc()), rxn.list, rxn.ctr)));
    #endif
    if (numberReactions > 0 && !(wcm_fused && rxn.on))
    {
        // Execute the kernel for the reaction, this kernel updates the lattice in-place, so only the src pointer is passed.
        PROF_CUDA_BEGIN(PROF_MPD_REACTION,cudaStream);
        #ifdef MPD_CUDA_3D_GRID_LAUNCH
        calculateReactionLaunchParameters(&gridSize, &threadBlockSize, TUNE_MPD_REACTION_BLOCK_X_SIZE, TUNE_MPD_REACTION_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
            #ifdef MPD_GLOBAL_S_MATRIX
            #ifdef MPD_GLOBAL_R_MATRIX
                    #ifdef MPD_FREAKYFAST
                    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::precomp_reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG, propZeroOrder, propFirstOrder, propSecondOrder, wcm_occ(lattice, lattice->getGPUMemorySrc()))));
                    #else
                    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG)));
                    #endif
            #else
                #ifdef MPD_FREAKYFAST
                CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::precomp_reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG, propZeroOrder, propFirstOrder, propSecondOrder, wcm_occ(lattice, lattice->getGPUMemorySrc()))));
                #else
                    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG)));
                #endif
            #endif
            #else
                    CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList)));
            #endif
                    #else
                    calculateReactionLaunchParameters(&gridXSize, &gridSize, &threadBlockSize, TUNE_MPD_REACTION_BLOCK_X_SIZE, TUNE_MPD_REACTION_BLOCK_Y_SIZE, latticeXSize, latticeYSize, latticeZSize);
        #ifdef MPD_GLOBAL_S_MATRIX
                CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), gridXSize, getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList, SG, RLG)));
        #else
                CUDA_EXCEPTION_EXECUTE((intmpdrdme_dev::reaction_kernel<<<gridSize,threadBlockSize,0,cudaStream>>>((unsigned int *)lattice->getGPUMemorySrc(), (uint8_t *)lattice->getGPUMemorySiteTypes(), (unsigned int *)lattice->getGPUMemorySrc(), gridXSize, getTimestepSeed(timestep,3), (unsigned int*)cudaOverflowList)));
        #endif
        #endif
        PROF_CUDA_END(PROF_MPD_REACTION,cudaStream);
    }

}

void IntMpdRdmeSolver::runTimestep(CudaIntLattice * lattice, uint32_t timestep)
{
    PROF_BEGIN(PROF_MPD_TIMESTEP);
    wcm_update_box(lattice, cudaStream, zeroOrder, zeroOrderSize);   // the host may have changed the lattice since
    wcm_launchTimestep(lattice, timestep);
    wcm_queue_box_copy(cudaStream);

    // Wait for the kernels to complete.
    PROF_BEGIN(PROF_MPD_SYNCHRONIZE);
    CUDA_EXCEPTION_CHECK(cudaStreamSynchronize(cudaStream));
    PROF_END(PROF_MPD_SYNCHRONIZE);

    // Handle any particle overflows.
    PROF_BEGIN(PROF_MPD_OVERFLOW);
    
    overflowTimesteps++;
#ifndef MPD_MAPPED_OVERFLOWS
    uint32_t overflowList[1+2*TUNE_MPD_MAX_PARTICLE_OVERFLOWS];
    CUDA_EXCEPTION_CHECK(cudaMemcpy(overflowList, cudaOverflowList, MPD_OVERFLOW_LIST_SIZE, cudaMemcpyDeviceToHost));
#else
    uint32_t *overflowList = (uint32_t*)cudaOverflowList;
#endif
    uint numberExceptions = overflowList[0];
    if (numberExceptions > 0)
    {
        Print::printf(Print::DEBUG, "%d overflows", numberExceptions);
        static const bool wcm_ovflog = getenv("WCM_OVERFLOW_LOG") != NULL;   // overflow log for verification runs
        if (wcm_ovflog) Print::printf(Print::INFO, "overflow: %u particle(s) at timestep %u", numberExceptions, timestep);
        
        // Make sure we did not exceed the overflow buffer.
        if (numberExceptions > TUNE_MPD_MAX_PARTICLE_OVERFLOWS)
            throw Exception("Too many particle overflows for the available buffer", numberExceptions);
            
        // Synchronize the lattice.
        lattice->copyFromGPU();
        
        // Go through each exception.
        for (uint i=0; i<numberExceptions; i++)
        {
            // Extract the index and particle type.
            lattice_size_t latticeIndex = overflowList[(i*2)+1];
            particle_t particle = overflowList[(i*2)+2];
            
            // Get the x, y, and z coordiantes.
            lattice_size_t x = latticeIndex%lattice->getXSize();
            lattice_size_t y = (latticeIndex/lattice->getXSize())%lattice->getYSize();
            lattice_size_t z = latticeIndex/(lattice->getXSize()*lattice->getYSize());
            
            // Put the particles back into a nearby lattice site.
            bool replacedParticle = false;
            for (uint searchRadius=0; !replacedParticle && searchRadius <= TUNE_MPD_MAX_OVERFLOW_REPLACEMENT_DIST; searchRadius++)
            {
                // Get the nearby sites.
                std::vector<lattice_coord_t> sites = lattice->getNearbySites(x,y,z,(searchRadius>0)?searchRadius-1:0,searchRadius);
                
                // TODO: Shuffle the sites.
                
                // Try to find one that in not fully occupied and of the same type.
                for (std::vector<lattice_coord_t>::iterator it=sites.begin(); it<sites.end(); it++)
                {
                    lattice_coord_t site = *it;
                    if (lattice->getOccupancy(site.x,site.y,site.z) < lattice->getMaxOccupancy() && lattice->getSiteType(site.x,site.y,site.z) == lattice->getSiteType(x,y,z))
                    {
                        lattice->addParticle(site.x, site.y, site.z, particle);
                        replacedParticle = true;
                        Print::printf(Print::VERBOSE_DEBUG, "Handled overflow of particle %d at site %d,%d,%d type=%d occ=%d by placing at site %d,%d,%d type=%d newocc=%d dist=%0.2f", particle, x, y, z, lattice->getSiteType(x,y,z), lattice->getOccupancy(x,y,z), site.x, site.y, site.z, lattice->getSiteType(site.x,site.y,site.z), lattice->getOccupancy(site.x,site.y,site.z), sqrt(pow((double)x-(double)site.x,2.0)+pow((double)y-(double)site.y,2.0)+pow((double)z-(double)site.z,2.0)));
                        break;
                    }
                }
            }
            
            // If we were not able to fix the exception, throw an error.
            if (!replacedParticle)
                throw Exception("Unable to find an available site to handle a particle overflow.");
        }
        
            
        // Copy the changes back to the GPU.
        lattice->copyToGPU(); wcm_recount_occupancy(lattice, cudaStream);
            
        // Reset the overflow list.
        CUDA_EXCEPTION_CHECK(cudaMemset(cudaOverflowList, 0, MPD_OVERFLOW_LIST_SIZE));
        
        // Track that we used the overflow list.
        overflowListUses++;
    }
    
    // If the overflow lsit is being used too often, print a warning.
    if (overflowTimesteps >= 1000)
    {
        if (overflowListUses > 10)
            Print::printf(Print::WARNING, "%d uses of the particle overflow list in the last 1000 timesteps, performance may be degraded.", overflowListUses);
        overflowTimesteps = 0;
        overflowListUses = 0;
    }
    PROF_END(PROF_MPD_OVERFLOW);

    PROF_CUDA_FINISH(cudaStream);
    PROF_END(PROF_MPD_TIMESTEP);
}

#ifdef MPD_CUDA_3D_GRID_LAUNCH
/**
 * Gets the launch parameters for launching an x diffusion kernel.
 */
void IntMpdRdmeSolver::calculateXLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    unsigned int xBlockXSize = min(maxXBlockSize,latticeXSize);
    unsigned int gridXSize = latticeXSize/xBlockXSize;
    if (gridXSize*xBlockXSize != latticeXSize)
	{
		// Find the largest number of warps that is divisible
		unsigned int tryx=32;
		while(tryx < maxXBlockSize)
		{
			if (latticeXSize % tryx == 0)
				xBlockXSize = tryx;
			tryx +=32;
		}
		gridXSize = latticeXSize/xBlockXSize;
	}
			
    (*gridSize).x = gridXSize;
    (*gridSize).y = latticeYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = xBlockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a y diffusion kernel.
 */
void IntMpdRdmeSolver::calculateYLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize/blockYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a z diffusion kernel.
 */
void IntMpdRdmeSolver::calculateZLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize;
    (*gridSize).z = latticeZSize/blockZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = blockZSize;
}

/**
 * Gets the launch parameters for launching a y diffusion kernel.
 */
void IntMpdRdmeSolver::calculateReactionLaunchParameters(dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    (*gridSize).x = latticeXSize/blockXSize;
    (*gridSize).y = latticeYSize/blockYSize;
    (*gridSize).z = latticeZSize;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

#else
/**
 * Gets the launch parameters for launching an x diffusion kernel.
 */
void IntMpdRdmeSolver::calculateXLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int maxXBlockSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    unsigned int xBlockXSize = min(maxXBlockSize,latticeXSize);
    *gridXSize = latticeXSize/xBlockXSize;
    if (gridXSize*xBlockXSize != latticeXSize)
	{
		// Find the largest number of warps that is divisible
		unsigned int tryx=32;
		while(tryx < maxXBlockSize)
		{
			if (latticeXSize % tryx == 0)
				xBlockXSize = tryx;
			tryx +=32;
		}
		gridXSize = latticeXSize/xBlockXSize;
	}

    (*gridSize).x = (*gridXSize)*latticeYSize;
    (*gridSize).y = latticeZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = xBlockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a y diffusion kernel.
 */
void IntMpdRdmeSolver::calculateYLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    *gridXSize = latticeXSize/blockXSize;
    (*gridSize).x = (*gridXSize)*(latticeYSize/blockYSize);
    (*gridSize).y = latticeZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}

/**
 * Gets the launch parameters for launching a z diffusion kernel.
 */
void IntMpdRdmeSolver::calculateZLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockZSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    *gridXSize = latticeXSize/blockXSize;
    (*gridSize).x = (*gridXSize)*(latticeYSize);
    (*gridSize).y = latticeZSize/blockZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = 1;
    (*threadBlockSize).z = blockZSize;
}

/**
 * Gets the launch parameters for launching a reaction diffusion kernel.
 */
void IntMpdRdmeSolver::calculateReactionLaunchParameters(unsigned int * gridXSize, dim3 * gridSize, dim3 * threadBlockSize, const unsigned int blockXSize, const unsigned int blockYSize, const unsigned int latticeXSize, const unsigned int latticeYSize, const unsigned int latticeZSize)
{
    *gridXSize = latticeXSize/blockXSize;
    (*gridSize).x = (*gridXSize)*(latticeYSize/blockYSize);
    (*gridSize).y = latticeZSize;
    (*gridSize).z = 1;
    (*threadBlockSize).x = blockXSize;
    (*threadBlockSize).y = blockYSize;
    (*threadBlockSize).z = 1;
}
#endif


namespace intmpdrdme_dev {
/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
#ifdef MPD_CUDA_3D_GRID_LAUNCH
__global__ void __launch_bounds__(TUNE_MPD_X_BLOCK_MAX_X_SIZE,1) mpd_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc)
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
__global__ void __launch_bounds__(TUNE_MPD_X_BLOCK_MAX_X_SIZE,1) mpd_x_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    if (wcm_block_outside_box(bx*blockDim.x, bx*blockDim.x+blockDim.x-1, by, by, bz, bz)) return;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeXIndex = (bx*blockDim.x) + threadIdx.x;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (by*latticeXSizeC) + latticeXIndex;
    unsigned int windowIndex = threadIdx.x+MPD_APRON_SIZE;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment.
    __shared__ unsigned int window[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_X_WINDOW_SIZE];
    __shared__ uint8_t occWin[MPD_X_WINDOW_SIZE];   // occupancy of each window site

    // Copy the x window from device memory into shared memory.
    // the destination site's old occupancy, the particle window and the site types are all loaded before a single
    // barrier; the empty-window test uses the OR of the occupancies each thread stored into occWin.
    const unsigned int oldO = outOcc[latticeIndex];
    const unsigned int wcm_or = copyXWindowFromLatticeOcc2(bx, inLattice, inOcc, window, occWin, latticeIndex, latticeXIndex, windowIndex);
    copyXWindowFromSites(bx, inSites, sitesWindow, latticeIndex, latticeXIndex, windowIndex);
    // a window with no particle at all (most blocks outside the cell) moves nothing: the destination sites become
    // empty, exactly as performPropagationOcc2 leaves them with no incoming particle (slots below the old occupancy zeroed,
    // occupancy 0); no choices or random draws are needed. Otherwise the kernel continues unchanged.
    if (!__syncthreads_or(wcm_or))
    {
        for (unsigned int w=0; w<oldO; w++) outLattice[latticeIndex + w*latticeXYZSizeC] = 0;
        outOcc[latticeIndex] = 0;
        return;
    }

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ uint8_t choices[MPD_X_WINDOW_SIZE*MPD_WORDS_PER_SITE];   // one byte per choice (values 0..3)

    // Make the choices.
    makeXDiffusionChoicesOcc(window, sitesWindow, occWin, choices, latticeIndex, latticeXIndex, windowIndex, blockDim.x, timestepHash);
    __syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Propagate the choices to the new lattice segment.
    performPropagationOcc2(outLattice, outOcc, window, occWin, choices, latticeIndex, windowIndex-1, windowIndex, windowIndex+1, MPD_X_WINDOW_SIZE, siteOverflowList, oldO);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
#ifdef MPD_CUDA_3D_GRID_LAUNCH
__global__ void __launch_bounds__(TUNE_MPD_Y_BLOCK_X_SIZE*TUNE_MPD_Y_BLOCK_Y_SIZE,4) mpd_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc)
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
__global__ void __launch_bounds__(TUNE_MPD_Y_BLOCK_X_SIZE*TUNE_MPD_Y_BLOCK_Y_SIZE,4) mpd_y_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    if (wcm_block_outside_box(bx*blockDim.x, bx*blockDim.x+blockDim.x-1, by*blockDim.y, by*blockDim.y+blockDim.y-1, bz, bz)) return;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (latticeYIndex*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;
    unsigned int windowYIndex = threadIdx.y+MPD_APRON_SIZE;
    unsigned int windowIndex = (windowYIndex*blockDim.x) + threadIdx.x;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[MPD_Y_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_Y_WINDOW_SIZE];
    __shared__ uint8_t occWin[MPD_Y_WINDOW_SIZE];   // occupancy of each window site

    // Copy the x window from device memory into shared memory.
    // the destination site's old occupancy, the particle window and the site types are all loaded before a single
    // barrier; the empty-window test uses the OR of the occupancies each thread stored into occWin.
    const unsigned int oldO = outOcc[latticeIndex];
    const unsigned int wcm_or = copyYWindowFromLatticeOcc2(inLattice, inOcc, window, occWin, latticeIndex, latticeYIndex, windowIndex, windowYIndex);
    copyYWindowFromSites(inSites, sitesWindow, latticeIndex, latticeYIndex, windowIndex, windowYIndex);
    // a window with no particle at all (most blocks outside the cell) moves nothing: the destination sites become
    // empty, exactly as performPropagationOcc2 leaves them with no incoming particle (slots below the old occupancy zeroed,
    // occupancy 0); no choices or random draws are needed. Otherwise the kernel continues unchanged.
    if (!__syncthreads_or(wcm_or))
    {
        for (unsigned int w=0; w<oldO; w++) outLattice[latticeIndex + w*latticeXYZSizeC] = 0;
        outOcc[latticeIndex] = 0;
        return;
    }

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ uint8_t choices[MPD_Y_WINDOW_SIZE*MPD_WORDS_PER_SITE];   // one byte per choice (values 0..3)

    // Make the choices.
    makeYDiffusionChoicesOcc(window, sitesWindow, occWin, choices, latticeIndex, latticeYIndex, windowIndex, windowYIndex, timestepHash);
    __syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagationOcc2(outLattice, outOcc, window, occWin, choices, latticeIndex, windowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, siteOverflowList, oldO);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
#ifdef MPD_CUDA_3D_GRID_LAUNCH
__global__ void __launch_bounds__(TUNE_MPD_Z_BLOCK_X_SIZE*TUNE_MPD_Z_BLOCK_Z_SIZE,1) mpd_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc)
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
__global__ void __launch_bounds__(TUNE_MPD_Z_BLOCK_X_SIZE*TUNE_MPD_Z_BLOCK_Z_SIZE,1) mpd_z_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc)
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    if (wcm_block_outside_box(bx*blockDim.x, bx*blockDim.x+blockDim.x-1, by, by, bz*blockDim.z, bz*blockDim.z+blockDim.z-1)) return;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeZIndex = (bz*blockDim.z) + threadIdx.z;
    unsigned int latticeIndex = (latticeZIndex*latticeXYSizeC) + (by*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;
    unsigned int windowZIndex = threadIdx.z+MPD_APRON_SIZE;
    unsigned int windowIndex = (windowZIndex*blockDim.x) + threadIdx.x;

    ///////////////////////////////////////////
    // Load the lattice into shared memory. //
    ///////////////////////////////////////////

    // Shared memory to store the lattice segment. Each lattice site has four particles, eight bits for each particle.
    __shared__ unsigned int window[MPD_Z_WINDOW_SIZE*MPD_WORDS_PER_SITE];
    __shared__ uint8_t sitesWindow[MPD_Z_WINDOW_SIZE];
    __shared__ uint8_t occWin[MPD_Z_WINDOW_SIZE];   // occupancy of each window site

    // Copy the x window from device memory into shared memory.
    // the destination site's old occupancy, the particle window and the site types are all loaded before a single
    // barrier; the empty-window test uses the OR of the occupancies each thread stored into occWin.
    const unsigned int oldO = outOcc[latticeIndex];
    const unsigned int wcm_or = copyZWindowFromLatticeOcc2(inLattice, inOcc, window, occWin, latticeIndex, latticeZIndex, windowIndex, windowZIndex);
    copyZWindowFromSites(inSites, sitesWindow, latticeIndex, latticeZIndex, windowIndex, windowZIndex);
    // a window with no particle at all (most blocks outside the cell) moves nothing: the destination sites become
    // empty, exactly as performPropagationOcc2 leaves them with no incoming particle (slots below the old occupancy zeroed,
    // occupancy 0); no choices or random draws are needed. Otherwise the kernel continues unchanged.
    if (!__syncthreads_or(wcm_or))
    {
        for (unsigned int w=0; w<oldO; w++) outLattice[latticeIndex + w*latticeXYZSizeC] = 0;
        outOcc[latticeIndex] = 0;
        return;
    }

    ////////////////////////////////////////
    // Make the choice for each particle. //
    ////////////////////////////////////////

    __shared__ uint8_t choices[MPD_Z_WINDOW_SIZE*MPD_WORDS_PER_SITE];   // one byte per choice (values 0..3)

    // Make the choices.
    makeZDiffusionChoicesOcc(window, sitesWindow, occWin, choices, latticeIndex, latticeZIndex, windowIndex, windowZIndex, timestepHash);
    __syncthreads();

    //////////////////////////////////////////////////////////
    // Create version of the lattice at the next time step. //
    //////////////////////////////////////////////////////////

    // Progate the choices to the new lattice segment.
    performPropagationOcc2(outLattice, outOcc, window, occWin, choices, latticeIndex, windowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, siteOverflowList, oldO);
}


/**
 * x, y and z multiparticle diffusion of one timestep in one kernel. A block takes a CX x CY x CZ core of sites and loads
 * the core plus a one-site apron on every side, (CX+2) x (CY+2) x (CZ+2) sites, into shared memory once. The three passes run on
 * that tile in the order of the three kernels, each reading the previous pass's result: the x pass computes the x core over the
 * whole y/z extent of the tile, the y pass the x and y core over the whole z extent, and the z pass writes the core straight to
 * the destination buffer. Each site of a pass is evaluated with the same function of the same inputs as in mpd_{x,y,z}_kernel: the
 * choice of slot w of site L is getRandomHashFloat(L, MPD_PARTICLE_COUNT_BITS, w, hash of that pass) against
 * lookupTransitionProbability with the site's type and its neighbour types along the pass axis, and arrivals are stored in the
 * same order (stay, from minus, from plus; ascending slots), the 17th and later going to the overflow list. Apron sites, which
 * belong to another block's core, are recomputed here with identical results; only the block that owns a site writes it and
 * records its overflows (a tile moved back at the lattice's high end owns only its part past the previous tile). At the tile edge
 * along a pass axis the missing outer neighbour type is replaced by the site's own type, as the kernels do for their outermost
 * apron sites: only the move into the tile is used there, and with q <= 0.5 (enforced in buildModel) the minus/plus decisions of
 * the choice formula do not depend on the other direction's probability. Particle values are held as 16 bits (the host enables
 * the path only when numberSpecies < 65535) and a site's choices as 2 bits per slot.
 */

// choices of the sites [lx0,lx1) x [ly0,ly1) x [lz0,lz1) of buffer buf along the axis with stride S (axis coordinate range [0,L))
template<int IX, int IY, int IZ, int AXIS, int NT>
__device__ __forceinline__ void wcm_fuse_choices(const uint16_t * __restrict__ buf, const uint8_t * __restrict__ occ, const uint8_t * __restrict__ st, const unsigned int * __restrict__ li, unsigned int * __restrict__ ch, const int lx0, const int lx1, const int ly0, const int ly1, const int lz0, const int lz1, const unsigned long long timestepHash)
{
    constexpr int NI = IX*IY*IZ;
    constexpr int S = (AXIS == 0) ? 1 : ((AXIS == 1) ? IX : IX*IY);
    constexpr int L = (AXIS == 0) ? IX : ((AXIS == 1) ? IY : IZ);
    constexpr int J = (NI + NT - 1) / NT;
    const int nx = lx1-lx0, ny = ly1-ly0, n = nx*ny*(lz1-lz0);
    // the thread's sites are handled together: their slot-0 transition probabilities are fetched before any is used (the
    // lookups are global-memory reads; most occupied sites hold one particle); further slots follow site by site
    int idx[J]; unsigned int oc[J], p0[J]; unsigned char sT[J], sM[J], sP[J]; float pm[J], pp[J];
    #pragma unroll
    for (int j = 0; j < J; j++)
    {
        const int k = threadIdx.x + j*NT;
        oc[j] = 0; p0[j] = 0; pm[j] = 0.0f; pp[j] = 0.0f; idx[j] = 0; sT[j] = sM[j] = sP[j] = 0;
        if (k < n)
        {
            const int lx = lx0 + k % nx, ly = ly0 + (k / nx) % ny, lz = lz0 + k / (nx*ny);
            const int i = (lz*IY + ly)*IX + lx;
            idx[j] = i; oc[j] = occ[i];
            if (oc[j] > 0)
            {
                const int a = (AXIS == 0) ? lx : ((AXIS == 1) ? ly : lz);
                sT[j] = st[i]; sM[j] = (a > 0) ? st[i-S] : sT[j]; sP[j] = (a < L-1) ? st[i+S] : sT[j];
                p0[j] = buf[i];
                if (p0[j] > 0) { pm[j] = lookupTransitionProbability(p0[j], sT[j], sM[j]); pp[j] = lookupTransitionProbability(p0[j], sT[j], sP[j]); }
            }
        }
    }
    #pragma unroll
    for (int j = 0; j < J; j++)
    {
        if (oc[j] == 0) continue;
        const int i = idx[j];
        const unsigned int latticeIndex = li[i];
        unsigned int c = 0;
        if (p0[j] > 0)
        {
            float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, 0u, timestepHash);
            unsigned int cc = (randomValue < pm[j])?(MPD_MOVE_MINUS):(MPD_MOVE_STAY);
            cc = (randomValue >= 0.5f && randomValue < (pp[j]+0.5f))?(MPD_MOVE_PLUS):(cc);
            c = cc;
        }
        for (unsigned int w = 1; w < oc[j]; w++)
        {
            const unsigned int particle = buf[i + w*NI];
            if (particle > 0)
            {
                float probMinus=lookupTransitionProbability(particle, sT[j], sM[j]);
                float probPlus=lookupTransitionProbability(particle, sT[j], sP[j]);
                float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, w, timestepHash);
                unsigned int cc = (randomValue < probMinus)?(MPD_MOVE_MINUS):(MPD_MOVE_STAY);
                cc = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_PLUS):(cc);
                c |= cc << (2*w);
            }
        }
        ch[i] = c;
    }
}

// the same choices as wcm_fuse_choices, evaluated over a compacted list of the occupied (site, slot) pairs of the
// region (a block-wide prefix sum of the occupancies) so that the threads of a warp all work on particles; each pair ORs its
// 2-bit choice into its site's word. A particle with both probabilities 0 stays without a draw (the formula gives STAY for any
// random value in [0,1)). Returns false (nothing done) when the region holds more than CAP particles.
template<int IX, int IY, int IZ, int AXIS, int NT, int CAP>
__device__ __forceinline__ bool wcm_fuse_choices_c(const uint16_t * __restrict__ buf, const uint8_t * __restrict__ occ, const uint8_t * __restrict__ st, const unsigned int * __restrict__ li, unsigned int * __restrict__ ch, int * __restrict__ wsum, uint16_t * __restrict__ items, const int lx0, const int lx1, const int ly0, const int ly1, const int lz0, const int lz1, const unsigned long long timestepHash)
{
    constexpr int NI = IX*IY*IZ;
    constexpr int S = (AXIS == 0) ? 1 : ((AXIS == 1) ? IX : IX*IY);
    constexpr int L = (AXIS == 0) ? IX : ((AXIS == 1) ? IY : IZ);
    constexpr int J = (NI + NT - 1) / NT;
    constexpr int NW = NT / 32;
    const int nx = lx1-lx0, ny = ly1-ly0, n = nx*ny*(lz1-lz0);
    // this thread's sites: consecutive region indices k = J*threadIdx.x + j (keeps the work list in site order)
    int idx[J]; unsigned int oc[J]; int mine = 0;
    #pragma unroll
    for (int j = 0; j < J; j++)
    {
        const int k = J*threadIdx.x + j;
        idx[j] = 0; oc[j] = 0;
        if (k < n)
        {
            const int lx = lx0 + k % nx, ly = ly0 + (k / nx) % ny, lz = lz0 + k / (nx*ny);
            idx[j] = (lz*IY + ly)*IX + lx; oc[j] = occ[idx[j]]; mine += oc[j];
            if (oc[j] > 0) ch[idx[j]] = 0;
        }
    }
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    int incl = mine;
    #pragma unroll
    for (int d = 1; d < 32; d <<= 1) { const int t = __shfl_up_sync(0xffffffffu, incl, d); if (lane >= d) incl += t; }
    if (lane == 31) wsum[warp] = incl;
    __syncthreads();
    if (warp == 0)
    {
        int v = (lane < NW) ? wsum[lane] : 0;
        #pragma unroll
        for (int d = 1; d < 32; d <<= 1) { const int t = __shfl_up_sync(0xffffffffu, v, d); if (lane >= d) v += t; }
        if (lane < NW) wsum[32 + lane] = v;
    }
    __syncthreads();
    const int total = wsum[32 + NW - 1];
    if (total > CAP) return false;   // block-uniform
    int base = incl - mine + ((warp > 0) ? wsum[32 + warp - 1] : 0);
    #pragma unroll
    for (int j = 0; j < J; j++)
        for (unsigned int w = 0; w < oc[j]; w++) items[base++] = (uint16_t)((idx[j] << 4) | w);
    __syncthreads();
    for (int t = threadIdx.x; t < total; t += NT)
    {
        const unsigned int it = items[t];
        const int i = (int)(it >> 4); const unsigned int w = it & 15u;
        const unsigned int particle = buf[i + w*NI];
        if (particle == 0) continue;
        const int a = (AXIS == 0) ? (i % IX) : ((AXIS == 1) ? ((i / IX) % IY) : (i / (IX*IY)));
        const unsigned char site = st[i];
        const unsigned char siteMinus = (a > 0) ? st[i-S] : site;
        const unsigned char sitePlus = (a < L-1) ? st[i+S] : site;
        float probMinus=lookupTransitionProbability(particle, site, siteMinus);
        float probPlus=lookupTransitionProbability(particle, site, sitePlus);
        unsigned int cc = MPD_MOVE_STAY;
        if (probMinus != 0.0f || probPlus != 0.0f)
        {
            float randomValue = getRandomHashFloat(li[i], MPD_PARTICLE_COUNT_BITS, w, timestepHash);
            cc = (randomValue < probMinus)?(MPD_MOVE_MINUS):(MPD_MOVE_STAY);
            cc = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_PLUS):(cc);
        }
        atomicOr(&ch[i], cc << (2*w));
    }
    return true;
}

// the slots w < o of a site whose 2-bit choice in c equals m, as a mask of bit 2w (ascending w = slot order)
__device__ __forceinline__ unsigned int wcm_fuse_match(const unsigned int c, const unsigned int o, const unsigned int m)
{
    unsigned int e = ~(c ^ (m * 0x55555555u));
    e &= (e >> 1) & 0x55555555u;
    return (o >= 16u) ? e : (e & ((1u << (2*o)) - 1u));
}

// propagation of the sites [lx0,lx1) x [ly0,ly1) x [lz0,lz1) along the axis with stride S, from (in, oin) into (out, oout); a
// site records its overflows only if it is in the block's core ([1,I-1) on the two other axes, own[] = the axes to test)
template<int IX, int IY, int IZ, int AXIS>
__device__ __forceinline__ void wcm_fuse_propagate(const uint16_t * __restrict__ in, const uint8_t * __restrict__ oin, const unsigned int * __restrict__ ch, const unsigned int * __restrict__ li, uint16_t * __restrict__ out, uint8_t * __restrict__ oout, const int lx0, const int lx1, const int ly0, const int ly1, const int lz0, const int lz1, unsigned int * __restrict__ siteOverflowList, const int ox0, const int ox1, const int oy0, const int oy1, const int oz0, const int oz1)
{
    constexpr int NI = IX*IY*IZ;
    constexpr int S = (AXIS == 0) ? 1 : ((AXIS == 1) ? IX : IX*IY);
    const int nx = lx1-lx0, ny = ly1-ly0, n = nx*ny*(lz1-lz0);
    for (int k = threadIdx.x; k < n; k += blockDim.x)
    {
        const int lx = lx0 + k % nx, ly = ly0 + (k / nx) % ny, lz = lz0 + k / (nx*ny);
        const int i = (lz*IY + ly)*IX + lx;
        const bool owner = lx-1 >= ox0 && lx-1 < ox1 && ly-1 >= oy0 && ly-1 < oy1 && lz-1 >= oz0 && lz-1 < oz1;
        unsigned int nn = 0;
        #define WCM_FPUT(p) { const unsigned int wcm_p = (p); \
            if (nn < MPD_WORDS_PER_SITE) out[i + nn*NI] = (uint16_t)wcm_p; \
            else if (owner) { int exceptionIndex = atomicAdd(siteOverflowList, 1); \
                   if (exceptionIndex < TUNE_MPD_MAX_PARTICLE_OVERFLOWS) { siteOverflowList[(exceptionIndex*2)+1]=li[i]; siteOverflowList[(exceptionIndex*2)+2]=wcm_p; } } \
            nn++; }
        if (oin[i] > 0)   { unsigned int m = wcm_fuse_match(ch[i],   oin[i],   MPD_MOVE_STAY);  while (m) { const unsigned int w = (__ffs(m)-1) >> 1; m &= m-1; WCM_FPUT(in[i + w*NI]) } }
        if (oin[i-S] > 0) { unsigned int m = wcm_fuse_match(ch[i-S], oin[i-S], MPD_MOVE_PLUS);  while (m) { const unsigned int w = (__ffs(m)-1) >> 1; m &= m-1; WCM_FPUT(in[i-S + w*NI]) } }
        if (oin[i+S] > 0) { unsigned int m = wcm_fuse_match(ch[i+S], oin[i+S], MPD_MOVE_MINUS); while (m) { const unsigned int w = (__ffs(m)-1) >> 1; m &= m-1; WCM_FPUT(in[i+S + w*NI]) } }
        #undef WCM_FPUT
        oout[i] = (uint8_t)((nn < MPD_WORDS_PER_SITE) ? nn : (unsigned int)MPD_WORDS_PER_SITE);
    }
}

// a site whose reaction check fired (checkForReaction with the substep-3 hash on its total propensity) is listed for
// wcm_rxn_tail_kernel, which chooses and applies the reaction exactly as precomp_reaction_kernel does.
__device__ __forceinline__ void wcm_fuse_list_reaction(const unsigned int latticeIndex, const float total, const wcm_fuse::Rxn & rxn)
{
    const unsigned int e = atomicAdd(&rxn.ctr[0], 1u);
    rxn.list[e] = make_uint2(latticeIndex, __float_as_uint(total));
}

__device__ __forceinline__ float wcm_fuse_propensity(const uint16_t * __restrict__ buf, const int i, const int NI, const unsigned int o, const uint8_t siteType, const wcm_fuse::Rxn & rxn)
{
    float totalReactionPropensity = read_element(rxn.qp0, siteType);
    const unsigned int base1 = siteType * numberSpeciesC;
    const unsigned int base2 = siteType * numberSpeciesC * numberSpeciesC;
    unsigned int ta = 0, tb = 0;
    while (ta < o)
    {
        float v[WCM_RXN_CHUNK]; bool ok[WCM_RXN_CHUNK];
        #pragma unroll
        for (int k=0; k<WCM_RXN_CHUNK; k++)
        {
            ok[k] = (ta < o);
            v[k] = 0.0f;
            if (ok[k])
            {
                const unsigned int qa = buf[i + ta*NI];
                v[k] = (tb == ta) ? read_element(rxn.qp1, base1 + (qa-1))
                                  : read_element(rxn.qp2, base2 + (qa-1)*numberSpeciesC + ((unsigned int)buf[i + tb*NI]-1));
                if (++tb >= o) { ta++; tb = ta; }
            }
        }
        #pragma unroll
        for (int k=0; k<WCM_RXN_CHUNK; k++) if (ok[k]) totalReactionPropensity += v[k];
    }
    return totalReactionPropensity;
}

// one tile, core origin (x0, y0, z0); the block writes (and records the overflows of) only the core sites it owns,
// [ox0,ox1) x [oy0,oy1) x [oz0,oz1) in core coordinates
template<int CX, int CY, int CZ, int NT, bool RXN>
__device__ __forceinline__ void wcm_fuse_tile(const int x0, const int y0, const int z0, const int ox0, const int ox1, const int oy0, const int oy1, const int oz0, const int oz1, const unsigned int* __restrict__ inLattice, const uint8_t * __restrict__ inSites, unsigned int* __restrict__ outLattice, const unsigned long long hashX, const unsigned long long hashY, const unsigned long long hashZ, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc, const bool wcmCompact, const wcm_fuse::Rxn rxn)
{
    typedef wcm_fuse::Tile<CX,CY,CZ> T;
    constexpr int IX = T::IX, IY = T::IY, IZ = T::IZ, NI = T::NI, NCORE = T::NCORE;
    constexpr int CPT = (NCORE + NT - 1) / NT;

    extern __shared__ uint4 wcm_fuse_smem[];
    uint16_t * A = (uint16_t *)wcm_fuse_smem;
    uint16_t * B = A + MPD_WORDS_PER_SITE*NI;
    unsigned int * ch = (unsigned int *)(B + MPD_WORDS_PER_SITE*NI);
    unsigned int * li = ch + NI;
    int * wsum = (int *)(li + NI);
    uint16_t * items = (uint16_t *)(wsum + 64);
    uint8_t * oA = (uint8_t *)(items + T::CAP);
    uint8_t * oB = oA + NI;
    uint8_t * st = oB + NI;

    const unsigned int X = latticeXSizeC, Y = latticeYSizeC, Z = latticeZSizeC, XY = latticeXYSizeC, XYZ = latticeXYZSizeC;

    // the destination's old occupancy of this thread's core sites (as the kernels load it at their start)
    unsigned int coreLi[CPT], oldO[CPT];
    #pragma unroll
    for (int c = 0; c < CPT; c++)
    {
        const int k = threadIdx.x + c*NT;
        coreLi[c] = 0xFFFFFFFFu; oldO[c] = 0;
        const int cx = k % CX, cy = (k / CX) % CY, cz = k / (CX*CY);
        if (k < NCORE && cx >= ox0 && cx < ox1 && cy >= oy0 && cy < oy1 && cz >= oz0 && cz < oz1)
        {
            coreLi[c] = (unsigned int)(z0+cz)*XY + (unsigned int)(y0+cy)*X + (unsigned int)(x0+cx);
            oldO[c] = outOcc[coreLi[c]];
        }
    }

    // load the tile (periodic lattice, as MPD_BOUNDARY_PERIODIC in the window copies); the thread's sites are loaded together:
    // occupancies and site types first, then slot 0 of each, then the further slots
    constexpr int JL = (NI + NT - 1) / NT;
    unsigned int orv = 0;
    {
        unsigned int l[JL], o[JL], t[JL], v0[JL];
        #pragma unroll
        for (int j = 0; j < JL; j++)
        {
            const int i = threadIdx.x + j*NT;
            l[j] = 0; o[j] = 0; t[j] = 0;
            if (i < NI)
            {
                const int lx = i % IX, ly = (i / IX) % IY, lz = i / (IX*IY);
                int gx = x0 + lx - 1, gy = y0 + ly - 1, gz = z0 + lz - 1;
                gx += (gx < 0) ? (int)X : 0; gx -= (gx >= (int)X) ? (int)X : 0;
                gy += (gy < 0) ? (int)Y : 0; gy -= (gy >= (int)Y) ? (int)Y : 0;
                gz += (gz < 0) ? (int)Z : 0; gz -= (gz >= (int)Z) ? (int)Z : 0;
                l[j] = (unsigned int)gz*XY + (unsigned int)gy*X + (unsigned int)gx;
                o[j] = inOcc[l[j]]; t[j] = inSites[l[j]];
            }
        }
        #pragma unroll
        for (int j = 0; j < JL; j++) v0[j] = (o[j] > 0) ? inLattice[l[j]] : 0u;
        #pragma unroll
        for (int j = 0; j < JL; j++)
        {
            const int i = threadIdx.x + j*NT;
            if (i >= NI) continue;
            li[i] = l[j]; oA[i] = (uint8_t)o[j]; st[i] = (uint8_t)t[j]; orv |= o[j];
            if (o[j] > 0) A[i] = (uint16_t)v0[j];
            for (unsigned int w0 = 1; w0 < o[j]; w0 += 8)
            {
                unsigned int v[8];
                #pragma unroll
                for (unsigned int k = 0; k < 8; k++) v[k] = (w0+k < o[j]) ? inLattice[l[j] + (w0+k)*XYZ] : 0u;
                #pragma unroll
                for (unsigned int k = 0; k < 8; k++) if (w0+k < o[j]) A[i + (w0+k)*NI] = (uint16_t)v[k];
            }
        }
    }
    // nothing in the tile: the core becomes empty
    if (!__syncthreads_or(orv))
    {
        #pragma unroll
        for (int c = 0; c < CPT; c++)
        {
            if (coreLi[c] == 0xFFFFFFFFu) continue;
            unsigned int o = 0;
            if (RXN)
            {
                // an empty site can still react through a zero-order reaction (propensity qp0 of its type)
                const int k = threadIdx.x + c*NT;
                const int i = ((1 + k / (CX*CY))*IY + 1 + (k / CX) % CY)*IX + 1 + k % CX;
                const uint8_t siteType = st[i];
                const float total = read_element(rxn.qp0, siteType);
                if (total != 0.0f && checkForReaction(coreLi[c], calculateReactionProbability(total), rxn.hash))
                    wcm_fuse_list_reaction(coreLi[c], total, rxn);
            }
            for (unsigned int w=o; w<oldO[c]; w++) outLattice[coreLi[c] + w*XYZ] = 0;
            outOcc[coreLi[c]] = (uint8_t)o;
        }
        return;
    }

    // x pass: A -> B on x core, all y, z
    if (!(wcmCompact && wcm_fuse_choices_c<IX,IY,IZ,0,NT,T::CAP>(A, oA, st, li, ch, wsum, items, 0, IX, 0, IY, 0, IZ, hashX))) { __syncthreads(); wcm_fuse_choices<IX,IY,IZ,0,NT>(A, oA, st, li, ch, 0, IX, 0, IY, 0, IZ, hashX); }
    __syncthreads();
    wcm_fuse_propagate<IX,IY,IZ,0>(A, oA, ch, li, B, oB, 1, IX-1, 0, IY, 0, IZ, siteOverflowList, ox0, ox1, oy0, oy1, oz0, oz1);
    __syncthreads();
    // y pass: B -> A on x core, y core, all z
    if (!(wcmCompact && wcm_fuse_choices_c<IX,IY,IZ,1,NT,T::CAP>(B, oB, st, li, ch, wsum, items, 1, IX-1, 0, IY, 0, IZ, hashY))) { __syncthreads(); wcm_fuse_choices<IX,IY,IZ,1,NT>(B, oB, st, li, ch, 1, IX-1, 0, IY, 0, IZ, hashY); }
    __syncthreads();
    wcm_fuse_propagate<IX,IY,IZ,1>(B, oB, ch, li, A, oA, 1, IX-1, 1, IY-1, 0, IZ, siteOverflowList, ox0, ox1, oy0, oy1, oz0, oz1);
    __syncthreads();
    // z pass: A -> destination lattice on the core
    if (!(wcmCompact && wcm_fuse_choices_c<IX,IY,IZ,2,NT,T::CAP>(A, oA, st, li, ch, wsum, items, 1, IX-1, 1, IY-1, 0, IZ, hashZ))) { __syncthreads(); wcm_fuse_choices<IX,IY,IZ,2,NT>(A, oA, st, li, ch, 1, IX-1, 1, IY-1, 0, IZ, hashZ); }
    __syncthreads();
    if (RXN)
    {
        // z pass into B (core sites), then the reaction check of each owned core site, then one store of the result
        // (a site whose check fires is listed and rewritten by wcm_rxn_tail_kernel).
        // The propensity terms of all core sites (per site: for each particle ta its first-order term, then its pair terms with
        // every later particle tb, the original summation order) are listed in the free A buffer, fetched in parallel by the whole block,
        // and each site then adds its own terms in that order onto qp0 of its type: the same float sum as precomp_reaction_kernel.
        // A tile with more terms than the buffer holds sums site by site (wcm_fuse_propensity).
        wcm_fuse_propagate<IX,IY,IZ,2>(A, oA, ch, li, B, oB, 1, IX-1, 1, IY-1, 1, IZ-1, siteOverflowList, ox0, ox1, oy0, oy1, oz0, oz1);
        __syncthreads();
        constexpr int TCAP = (MPD_WORDS_PER_SITE*NI*2) / 4;
        unsigned int * terms = (unsigned int *)A;
        int ci[CPT]; unsigned int co[CPT], rm[CPT]; int nt[CPT]; int mine = 0;
        #pragma unroll
        for (int c = 0; c < CPT; c++)
        {
            const int k = threadIdx.x + c*NT;
            ci[c] = ((1 + k / (CX*CY))*IY + 1 + (k / CX) % CY)*IX + 1 + k % CX;
            co[c] = (coreLi[c] == 0xFFFFFFFFu) ? 0u : (unsigned int)oB[ci[c]];
            // only the terms of reactant species are listed: every other term is exactly 0 (qp1/qp2 entries are sums over the
            // reactions with that reactant), and adding +-0 to the running sum (which starts at qp0 >= +0) leaves it unchanged
            rm[c] = 0;
            for (unsigned int w = 0; w < co[c]; w++) if (rxn.reactive[(unsigned int)B[ci[c] + w*NI] - 1]) rm[c] |= 1u << w;
            const unsigned int r = __popc(rm[c]);
            nt[c] = (int)(r*(r+1)/2); mine += nt[c];
        }
        const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
        constexpr int NW = NT / 32;
        int incl = mine;
        #pragma unroll
        for (int d = 1; d < 32; d <<= 1) { const int t = __shfl_up_sync(0xffffffffu, incl, d); if (lane >= d) incl += t; }
        if (lane == 31) wsum[warp] = incl;
        __syncthreads();
        if (warp == 0)
        {
            int v = (lane < NW) ? wsum[lane] : 0;
            #pragma unroll
            for (int d = 1; d < 32; d <<= 1) { const int t = __shfl_up_sync(0xffffffffu, v, d); if (lane >= d) v += t; }
            if (lane < NW) wsum[32 + lane] = v;
        }
        __syncthreads();
        const int nterms = wsum[32 + NW - 1];
        const bool listed = (nterms <= TCAP);   // block-uniform
        int off[CPT];
        {
            int base = incl - mine + ((warp > 0) ? wsum[32 + warp - 1] : 0);
            #pragma unroll
            for (int c = 0; c < CPT; c++)
            {
                off[c] = base;
                if (listed)
                    for (unsigned int ma = rm[c]; ma; ma &= ma-1)
                    {
                        const unsigned int ta = __ffs(ma) - 1;
                        for (unsigned int mb = ma; mb; mb &= mb-1) terms[base++] = ((unsigned int)ci[c] << 8) | (ta << 4) | (unsigned int)(__ffs(mb) - 1);
                    }
                else base += nt[c];
            }
        }
        __syncthreads();
        if (listed)
        {
            for (int t = threadIdx.x; t < nterms; t += NT)
            {
                const unsigned int d = terms[t];
                const int i = (int)(d >> 8); const unsigned int ta = (d >> 4) & 15u, tb = d & 15u;
                const unsigned int siteType = st[i];
                const unsigned int qa = B[i + ta*NI];
                const float v = (tb == ta) ? read_element(rxn.qp1, siteType*numberSpeciesC + (qa-1))
                                           : read_element(rxn.qp2, siteType*numberSpeciesC*numberSpeciesC + (qa-1)*numberSpeciesC + ((unsigned int)B[i + tb*NI]-1));
                terms[t] = __float_as_uint(v);
            }
        }
        __syncthreads();
        #pragma unroll
        for (int c = 0; c < CPT; c++)
        {
            if (coreLi[c] == 0xFFFFFFFFu) continue;
            const int i = ci[c];
            const unsigned int latticeIndex = coreLi[c];
            const uint8_t siteType = st[i];
            unsigned int o = co[c];
            float total;
            if (listed)
            {
                total = read_element(rxn.qp0, siteType);
                for (int t = off[c]; t < off[c] + nt[c]; t++) total += __uint_as_float(terms[t]);
            }
            else total = wcm_fuse_propensity(B, i, NI, o, siteType, rxn);
            if (total != 0.0f && checkForReaction(latticeIndex, calculateReactionProbability(total), rxn.hash))
                wcm_fuse_list_reaction(latticeIndex, total, rxn);
            for (unsigned int w=0; w<o; w++) outLattice[latticeIndex + w*XYZ] = B[i + w*NI];
            for (unsigned int w=o; w < oldO[c]; w++) outLattice[latticeIndex + w*XYZ] = 0;
            outOcc[latticeIndex] = (uint8_t)o;
        }
        return;
    }
    #pragma unroll
    for (int c = 0; c < CPT; c++)
    {
        const int k = threadIdx.x + c*NT;
        if (coreLi[c] == 0xFFFFFFFFu) continue;
        const int lx = 1 + k % CX, ly = 1 + (k / CX) % CY, lz = 1 + k / (CX*CY);
        const int i = (lz*IY + ly)*IX + lx;
        constexpr int S = IX*IY;
        const unsigned int latticeIndex = coreLi[c];
        unsigned int n = 0;
        #define WCM_GPUT(p) { const unsigned int wcm_p = (p); \
            if (n < MPD_WORDS_PER_SITE) outLattice[latticeIndex + n*XYZ] = wcm_p; \
            else { int exceptionIndex = atomicAdd(siteOverflowList, 1); \
                   if (exceptionIndex < TUNE_MPD_MAX_PARTICLE_OVERFLOWS) { siteOverflowList[(exceptionIndex*2)+1]=latticeIndex; siteOverflowList[(exceptionIndex*2)+2]=wcm_p; } } \
            n++; }
        if (oA[i] > 0)   { unsigned int m = wcm_fuse_match(ch[i],   oA[i],   MPD_MOVE_STAY);  while (m) { const unsigned int w = (__ffs(m)-1) >> 1; m &= m-1; WCM_GPUT(A[i + w*NI]) } }
        if (oA[i-S] > 0) { unsigned int m = wcm_fuse_match(ch[i-S], oA[i-S], MPD_MOVE_PLUS);  while (m) { const unsigned int w = (__ffs(m)-1) >> 1; m &= m-1; WCM_GPUT(A[i-S + w*NI]) } }
        if (oA[i+S] > 0) { unsigned int m = wcm_fuse_match(ch[i+S], oA[i+S], MPD_MOVE_MINUS); while (m) { const unsigned int w = (__ffs(m)-1) >> 1; m &= m-1; WCM_GPUT(A[i+S + w*NI]) } }
        #undef WCM_GPUT
        const unsigned int newO = (n < MPD_WORDS_PER_SITE) ? n : (unsigned int)MPD_WORDS_PER_SITE;
        for (unsigned int w=newO; w < oldO[c]; w++) outLattice[latticeIndex + w*XYZ] = 0;
        outOcc[latticeIndex] = (uint8_t)newO;
    }
}

template<int CX, int CY, int CZ, int NT, int MINB, bool RXN>
__global__ void __launch_bounds__(NT, MINB) wcm_mpd_xyz_kernel(const unsigned int* __restrict__ inLattice, const uint8_t * __restrict__ inSites, unsigned int* __restrict__ outLattice, const unsigned long long hashX, const unsigned long long hashY, const unsigned long long hashZ, unsigned int* siteOverflowList, const uint8_t* __restrict__ inOcc, uint8_t* __restrict__ outOcc, const bool wcmCompact, const wcm_fuse::Rxn rxn)
{
    // the tiles start at the box's low corner (clamped to the lattice) and cover the box; a tile that would cross the
    // lattice's high end is moved back to end there and owns only its sites past the previous tile
    const int b0[3] = { wcm_boxD[0], wcm_boxD[2], wcm_boxD[4] }, b1[3] = { wcm_boxD[1], wcm_boxD[3], wcm_boxD[5] };
    const int C[3] = { CX, CY, CZ }, Ls[3] = { (int)latticeXSizeC, (int)latticeYSizeC, (int)latticeZSizeC };
    const int bi[3] = { (int)blockIdx.x, (int)blockIdx.y, (int)blockIdx.z };
    int org[3], own0[3], own1[3];
    #pragma unroll
    for (int d = 0; d < 3; d++)
    {
        const int lo = max(b0[d], 0), hi = min(b1[d], Ls[d]-1);
        if (hi < lo) return;   // empty box
        if (bi[d] > (hi - lo) / C[d]) return;   // beyond the box
        const int ns = lo + bi[d]*C[d];
        org[d] = min(ns, Ls[d] - C[d]);
        own0[d] = ns - org[d]; own1[d] = min(ns + C[d], Ls[d]) - org[d];
    }
    wcm_fuse_tile<CX,CY,CZ,NT,RXN>(org[0], org[1], org[2], own0[0], own1[0], own0[1], own1[1], own0[2], own1[2], inLattice, inSites, outLattice, hashX, hashY, hashZ, siteOverflowList, inOcc, outOcc, wcmCompact, rxn);
}

/**
 * Multiparticle diffusion performed by copying the lattice section to shared memory, making a choice for each lattice
 * site, storing the new lattice into shared memory, and then updating the global lattice.
 */
#if defined(MPD_GLOBAL_S_MATRIX) && defined(MPD_GLOBAL_R_MATRIX)
// determineReactionIndex over the reactions that can have a nonzero propensity at this site only (zero-order ones and
// those with a reactant present), visited in ascending index with the same propensity function and the same subtraction as the
// full scan: a reaction left out has propensity 0 (or -0), which the full scan skips as well, so the chosen index is the same.
__device__ __noinline__ unsigned int wcm_determineReactionIndex(const uint8_t siteType, const unsigned int * __restrict__ particles, const unsigned int latticeIndex, const float totalReactionPropensity, const unsigned long long timestepHash, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
{
    if (!wcm_rxnListsOn) return determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
    float randomPropensity = getRandomHashFloat(latticeIndex, 1, 1, timestepHash)*totalReactionPropensity;
    unsigned int reactionIndex = 0;
    unsigned int val[MPD_PARTICLES_PER_SITE], cur[MPD_PARTICLES_PER_SITE], end[MPD_PARTICLES_PER_SITE];
    int m = 0;
    for (int i=0; i<MPD_PARTICLES_PER_SITE; i++)
    {
        const unsigned int p = particles[i];
        if (p == 0) continue;
        bool seen = false;
        for (int k=0; k<m; k++) seen |= (val[k] == p);
        if (seen) continue;
        val[m] = p; cur[m] = wcm_specRxnOff[p]; end[m] = wcm_specRxnOff[p+1]; m++;
    }
    unsigned int z = 0;
    const unsigned int nz = wcm_nZeroRxn;
    while (true)
    {
        unsigned int best = 0xFFFFFFFFu;
        for (int k=0; k<m; k++) if (cur[k] < end[k]) best = min(best, wcm_specRxn[cur[k]]);
        if (z < nz) best = min(best, wcm_zeroRxn[z]);
        if (best == 0xFFFFFFFFu) break;
        for (int k=0; k<m; k++) if (cur[k] < end[k] && wcm_specRxn[cur[k]] == best) cur[k]++;
        if (z < nz && wcm_zeroRxn[z] == best) z++;
        float propensity = calculateReactionPropensity(siteType, particles, best, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
        if (propensity > 0.0f)
        {
            if (randomPropensity > 0.0f)
                reactionIndex = best;
            randomPropensity -= propensity;
        }
    }
    return reactionIndex;
}

// evaluateReaction without the per-event device malloc and the copy of the reaction's whole stoichiometry row: a
// particle is dropped while its species still has reactant count left (S < 0 in the copy), else kept, in slot order; then the
// products are added in ascending species order (the positive entries, which the drop loop never changes), with the same
// overflow handling; the remaining slots are cleared.
__device__ __noinline__ void wcm_evaluateReaction(const unsigned int latticeIndex, const uint8_t siteType, unsigned int * __restrict__ particles, const unsigned int reactionIndex, unsigned int * siteOverflowList, const int8_t* __restrict__ SG)
{
    if (!wcm_rxnListsOn) { evaluateReaction(latticeIndex, siteType, particles, reactionIndex, siteOverflowList, SG); return; }
    unsigned int rs[8]; int rc[8]; int nr = 0;
    for (unsigned int q = wcm_reacOff[reactionIndex]; q < wcm_reacOff[reactionIndex+1] && nr < 8; q++, nr++)
        { const unsigned int e = wcm_reac[q]; rs[nr] = e >> 8; rc[nr] = (int)(e & 0xFFu); }
    int nextParticle=0;
    for (uint i=0; i<MPD_PARTICLES_PER_SITE; i++)
    {
        const unsigned int particle = particles[i];
        if (particle > 0)
        {
            bool drop = false;
            for (int k=0; k<nr; k++) if (rs[k] == particle-1 && rc[k] > 0) { rc[k]--; drop = true; break; }
            if (!drop) particles[nextParticle++] = particle;
        }
    }
    for (unsigned int q = wcm_prodOff[reactionIndex]; q < wcm_prodOff[reactionIndex+1]; q++)
    {
        const unsigned int e = wcm_prod[q], sp = e >> 8, cnt = e & 0xFFu;
        for (uint j=0; j<cnt; j++)
        {
            if (nextParticle < MPD_PARTICLES_PER_SITE)
                particles[nextParticle++] = sp+1;
            else
            {
                int exceptionIndex = atomicAdd(siteOverflowList, 1);
                if (exceptionIndex < TUNE_MPD_MAX_PARTICLE_OVERFLOWS)
                {
                    siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
                    siteOverflowList[(exceptionIndex*2)+2]=sp+1;
                }
            }
        }
    }
    while (nextParticle < MPD_PARTICLES_PER_SITE)
        particles[nextParticle++] = 0;
}
#endif

#ifdef MPD_CUDA_3D_GRID_LAUNCH
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG)
#endif
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList)
#endif
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;
#else
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG)
#endif
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,1) reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned int gridXSize, const unsigned long long timestepHash, unsigned int* siteOverflowList)
#endif
{
    __shared__ unsigned int bx, by, bz;
    calculateBlockPosition(&bx, &by, &bz, gridXSize);
#endif

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (latticeYIndex*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;

    ///////////////////////////////////////////
    // Load the particles and site.          //
    ///////////////////////////////////////////

    unsigned int particles[MPD_WORDS_PER_SITE];
    for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
        particles[w] = inLattice[latticeIndex+latticeOffset];
    uint8_t siteType = inSites[latticeIndex];

    ////////////////////////////////////////
    // Perform the reactions.             //
    ////////////////////////////////////////

    // Calculate the kinetic rate for each reaction at this site.
    float totalReactionPropensity = 0.0f;
    for (int i=0; i<numberReactionsC; i++)
    {
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        totalReactionPropensity += calculateReactionPropensity(siteType, particles, i, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        totalReactionPropensity += calculateReactionPropensity(siteType, particles, i, RLG);
#endif
#else
        totalReactionPropensity += calculateReactionPropensity(siteType, particles, i);
#endif
    }

	// If propensity is zero, no reaction can occur.
	if(totalReactionPropensity == 0.0f)
		return;

    // See if a reaction occurred at the site.
    float reactionProbability = calculateReactionProbability(totalReactionPropensity);
    unsigned int reactionOccurred = checkForReaction(latticeIndex, reactionProbability, timestepHash);

    // If there was a reaction, process it.
    if (reactionOccurred)
    {
        // Figure out which reaction occurred.
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        unsigned int reactionIndex = determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash, RLG);
#endif
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash);
#endif

        // Construct the new site.
#ifdef MPD_GLOBAL_S_MATRIX
        evaluateReaction(latticeIndex, siteType, particles, reactionIndex, siteOverflowList, SG);
#else
        evaluateReaction(latticeIndex, siteType, particles, reactionIndex, siteOverflowList);
#endif

        // Copy the new particles back into the lattice.
        for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
             outLattice[latticeIndex+latticeOffset] = particles[w];
    }
}

#ifdef MPD_GLOBAL_R_MATRIX
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,4) precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, uint8_t* __restrict__ occ)
#else
__global__ void __launch_bounds__(TUNE_MPD_REACTION_BLOCK_X_SIZE*TUNE_MPD_REACTION_BLOCK_Y_SIZE,4) precomp_reaction_kernel(const unsigned int* inLattice, const uint8_t * inSites, unsigned int* outLattice, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const float* __restrict__ qp0, const float* __restrict__ qp1, const float* __restrict__ qp2, uint8_t* __restrict__ occ)
#endif
{
    unsigned int bx=blockIdx.x, by=blockIdx.y, bz=blockIdx.z;

    if (wcm_block_outside_box(bx*blockDim.x, bx*blockDim.x+blockDim.x-1, by*blockDim.y, by*blockDim.y+blockDim.y-1, bz, bz)) return;

    // Figure out the offset of this thread in the lattice and the lattice segment.
    unsigned int latticeYIndex = (by*blockDim.y) + threadIdx.y;
    unsigned int latticeIndex = (bz*latticeXYSizeC) + (latticeYIndex*latticeXSizeC) + (bx*blockDim.x) + threadIdx.x;

    ///////////////////////////////////////////
    // Load the particles and site.          //
    ///////////////////////////////////////////

    unsigned int particles[MPD_WORDS_PER_SITE];
    const unsigned int occHere = occ[latticeIndex];   // slots >= occ are zero
    for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
        WCM_LOAD(particles[w], inLattice[latticeIndex+latticeOffset], occHere, w)
    uint8_t siteType = inSites[latticeIndex];

    ////////////////////////////////////////
    // Perform the reactions.             //
    ////////////////////////////////////////

    // Calculate the kinetic rate for each reaction at this site.
    float totalReactionPropensity = read_element(qp0,siteType);
    //float totalReactionPropensity = qp0[siteType];
    // the propensity terms in their original order (for each nonzero particle: its first-order term, then its pair
    // term with every later nonzero particle; slots >= occHere are zero, 021 invariant) are fetched four at a time with the
    // loads issued first, then added in that same order: the float sum is the one of the former nested loops.
    {
        unsigned int q[MPD_WORDS_PER_SITE];
        unsigned int n = 0;
        for (uint i=0; i<occHere; i++) if (particles[i] > 0) q[n++] = particles[i];
        const unsigned int base1 = siteType * numberSpeciesC;
        const unsigned int base2 = siteType * numberSpeciesC * numberSpeciesC;
        unsigned int ta = 0, tb = 0;   // next term: tb == ta -> first-order term of q[ta]; tb > ta -> pair (q[ta], q[tb])
        // WCM_RXN_CHUNK terms per round of loads (was 4): a crowded site (8-16 particles, up to 136 terms) waited one
        // global round trip per 4 terms; the sum is still taken term by term in the original order.
        while (ta < n)
        {
            float v[WCM_RXN_CHUNK]; bool ok[WCM_RXN_CHUNK];
            #pragma unroll
            for (int k=0; k<WCM_RXN_CHUNK; k++)
            {
                ok[k] = (ta < n);
                v[k] = 0.0f;
                if (ok[k])
                {
                    v[k] = (tb == ta) ? read_element(qp1, base1 + (q[ta]-1))
                                      : read_element(qp2, base2 + (q[ta]-1)*numberSpeciesC + (q[tb]-1));
                    if (++tb >= n) { ta++; tb = ta; }
                }
            }
            #pragma unroll
            for (int k=0; k<WCM_RXN_CHUNK; k++) if (ok[k]) totalReactionPropensity += v[k];
        }
    }


	// If propensity is zero, no reaction can occur.
	if(totalReactionPropensity == 0.0f)
		return;

    // See if a reaction occurred at the site.
    float reactionProbability = calculateReactionProbability(totalReactionPropensity);
    unsigned int reactionOccurred = checkForReaction(latticeIndex, reactionProbability, timestepHash);

    // If there was a reaction, process it.
    if (reactionOccurred)
    {
        // Figure out which reaction occurred.
#ifdef MPD_GLOBAL_S_MATRIX
#ifdef MPD_GLOBAL_R_MATRIX
        unsigned int reactionIndex = wcm_determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash, RLG);
#endif
#else
        unsigned int reactionIndex = determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash);
#endif

        // Construct the new site.
#ifdef MPD_GLOBAL_S_MATRIX
        wcm_evaluateReaction(latticeIndex, siteType, particles, reactionIndex, siteOverflowList, SG);
#else
        evaluateReaction(latticeIndex, siteType, particles, reactionIndex, siteOverflowList);
#endif

        // Copy the new particles back into the lattice.
        for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
             outLattice[latticeIndex+latticeOffset] = particles[w];
        unsigned int occNew = 0;
        for (uint w=0; w<MPD_WORDS_PER_SITE; w++) if (particles[w] != 0) occNew = w+1;
        occ[latticeIndex] = (uint8_t)occNew;
    }
}

#ifdef WCM_FUSE_RXN_AVAILABLE
// the reacting sites listed by wcm_mpd_xyz_kernel<RXN>: the rest of precomp_reaction_kernel's reacting branch
// (wcm_determineReactionIndex, wcm_evaluateReaction, all slots stored, occupancy recounted) on the site's particles as the fused
// kernel stored them, with its total propensity. The sites are independent, so the list order does not matter; the reaction
// overflows are appended after all diffusion overflows of the step, as with the separate reaction kernel. The last block to
// finish clears the counters for the next step.
__global__ void __launch_bounds__(128) wcm_rxn_tail_kernel(unsigned int* lattice, const uint8_t * __restrict__ inSites, const unsigned long long timestepHash, unsigned int* siteOverflowList, const __restrict__ int8_t *SG, const __restrict__ uint8_t *RLG, const unsigned int* __restrict__ reactionOrdersG, const unsigned int* __restrict__ reactionSitesG, const unsigned int* __restrict__ D1G, const unsigned int* __restrict__ D2G, const float* __restrict__ reactionRatesG, uint8_t* __restrict__ occ, const uint2* __restrict__ list, unsigned int* ctr)
{
    const unsigned int n = *(volatile unsigned int *)&ctr[0];
    for (unsigned int e = blockIdx.x*blockDim.x + threadIdx.x; e < n; e += gridDim.x*blockDim.x)
    {
        const uint2 it = list[e];
        const unsigned int latticeIndex = it.x;
        const float totalReactionPropensity = __uint_as_float(it.y);
        unsigned int particles[MPD_WORDS_PER_SITE];
        const unsigned int occHere = occ[latticeIndex];
        for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
            WCM_LOAD(particles[w], lattice[latticeIndex+latticeOffset], occHere, w)
        const uint8_t siteType = inSites[latticeIndex];
        unsigned int reactionIndex = wcm_determineReactionIndex(siteType, particles, latticeIndex, totalReactionPropensity, timestepHash, RLG, reactionOrdersG, reactionSitesG, D1G, D2G, reactionRatesG);
        wcm_evaluateReaction(latticeIndex, siteType, particles, reactionIndex, siteOverflowList, SG);
        for (uint w=0, latticeOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC)
             lattice[latticeIndex+latticeOffset] = particles[w];
        unsigned int occNew = 0;
        for (uint w=0; w<MPD_WORDS_PER_SITE; w++) if (particles[w] != 0) occNew = w+1;
        occ[latticeIndex] = (uint8_t)occNew;
    }
    __syncthreads();
    if (threadIdx.x == 0)
    {
        __threadfence();
        if (atomicAdd(&ctr[1], 1u) == gridDim.x - 1) { ctr[0] = 0; ctr[1] = 0; }
    }
}
#endif
}

}
}