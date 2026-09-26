/*
 * Part of Lattice Microbes; same license as the rest of the repository (see the header of IntMpdRdmeSolver.cu).
 *
 * Author(s): Ron Acda
 *   (Ron Acda: using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
 */
/*
 * per-site occupancy for the int-lattice MPD-RDME kernels (IntMpdRdmeSolver only).
 *
 * Every device particle buffer b has an occupancy byte per site, occ_b[site]. Invariant: in buffer b every slot >= occ_b[site]
 * of that site is zero. Readers load only slots < occ from global memory and zero-fill the rest of the shared window, so the choice
 * and propagation code sees exactly the windows it saw before (same particles, same RNG draws per slot index). Writers store slot w
 * of a site only if w < max(new occupancy, old occupancy of the destination site) (above that the destination is already zero) and
 * record the new occupancy. Host uploads are followed by a recount of the uploaded buffer; both buffers start zeroed.
 * The copy functions mirror copy{X,Y,Z}WindowFromLattice of lattice_sim_1d_dev.cu line for line, indices included.
 */

#define WCM_LOAD(dst, src_ptr, o, w) { unsigned int _v = 0; if ((w) < (o)) _v = (src_ptr); (dst) = _v; }

inline __device__ void copyXWindowFromLatticeOcc(const unsigned int bx, const unsigned int * lattice, const uint8_t * occ, unsigned int * window, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex)
{
    if (latticeXIndex < latticeXSizeC)
    {
        const unsigned int o = occ[latticeIndex];
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
            WCM_LOAD(window[windowIndex+windowOffset], lattice[latticeIndex+latticeOffset], o, w)

        #if MPD_APRON_SIZE > 0
        int threadBlockWidth = ((bx+1)*blockDim.x <= latticeXSizeC)?(blockDim.x):(latticeXSizeC-(bx*blockDim.x));

        if (windowIndex >= threadBlockWidth)
        {
            #if defined MPD_BOUNDARY_PERIODIC
            const unsigned int ai = (latticeXIndex>=threadBlockWidth)?(latticeIndex-threadBlockWidth):(latticeIndex-threadBlockWidth+latticeXSizeC);
            const unsigned int ao = occ[ai];
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
                WCM_LOAD(window[windowIndex+windowOffset-threadBlockWidth], lattice[ai+latticeOffset], ao, w)
            #else
            if (latticeXIndex>=threadBlockWidth)
            {
                const unsigned int ai = latticeIndex-threadBlockWidth; const unsigned int ao = occ[ai];
                for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
                    WCM_LOAD(window[windowIndex+windowOffset-threadBlockWidth], lattice[ai+latticeOffset], ao, w)
            }
            else
                for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_X_WINDOW_SIZE)
                    window[windowIndex+windowOffset-threadBlockWidth] = MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
        }

        if (windowIndex < 2*MPD_APRON_SIZE)
        {
            #if defined MPD_BOUNDARY_PERIODIC
            const unsigned int ai = (latticeXIndex<latticeXSizeC-threadBlockWidth)?(latticeIndex+threadBlockWidth):(latticeIndex+threadBlockWidth-latticeXSizeC);
            const unsigned int ao = occ[ai];
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
                WCM_LOAD(window[windowIndex+windowOffset+threadBlockWidth], lattice[ai+latticeOffset], ao, w)
            #else
            if (latticeXIndex<latticeXSizeC-threadBlockWidth)
            {
                const unsigned int ai = latticeIndex+threadBlockWidth; const unsigned int ao = occ[ai];
                for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_X_WINDOW_SIZE)
                    WCM_LOAD(window[windowIndex+windowOffset+threadBlockWidth], lattice[ai+latticeOffset], ao, w)
            }
            else
                for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_X_WINDOW_SIZE)
                    window[windowIndex+windowOffset+threadBlockWidth] = MPD_BOUNDARY_PARTICLE_VALUE;
            #endif
        }
        #endif
    }
}

inline __device__ void copyYWindowFromLatticeOcc(const unsigned int* lattice, const uint8_t * occ, unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeYIndex, const unsigned int windowIndex, const unsigned int windowYIndex)
{
    const unsigned int o = occ[latticeIndex];
    for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
        WCM_LOAD(window[windowIndex+windowOffset], lattice[latticeIndex+latticeOffset], o, w)

    #if MPD_APRON_SIZE > 0
    if (windowYIndex < 2*MPD_APRON_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?(latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)):(latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)+(latticeXSizeC*latticeYSizeC));
        const unsigned int ao = occ[ai];
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
            WCM_LOAD(window[windowIndex+windowOffset-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        #else
        if (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)
        {
            const unsigned int ai = latticeIndex-(latticeXSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai];
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
                WCM_LOAD(window[windowIndex+windowOffset-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        }
        else
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Y_WINDOW_SIZE)
                window[windowIndex+windowOffset-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE;
        #endif
    }

    if (windowYIndex >= TUNE_MPD_Y_BLOCK_Y_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?(latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)):(latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)-(latticeXSizeC*latticeYSizeC));
        const unsigned int ao = occ[ai];
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
            WCM_LOAD(window[windowIndex+windowOffset+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        #else
        if (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)
        {
            const unsigned int ai = latticeIndex+(latticeXSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai];
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Y_WINDOW_SIZE)
                WCM_LOAD(window[windowIndex+windowOffset+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        }
        else
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Y_WINDOW_SIZE)
                window[windowIndex+windowOffset+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE;
        #endif
    }
    #endif
}

inline __device__ void copyZWindowFromLatticeOcc(const unsigned int* lattice, const uint8_t * occ, unsigned int* window, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex)
{
    const unsigned int o = occ[latticeIndex];
    for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
        WCM_LOAD(window[windowIndex+windowOffset], lattice[latticeIndex+latticeOffset], o, w)

    #if MPD_APRON_SIZE > 0
    if (windowZIndex < 2*MPD_APRON_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?(latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)):(latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)+(latticeXYZSizeC));
        const unsigned int ao = occ[ai];
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
            WCM_LOAD(window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        #else
        if (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)
        {
            const unsigned int ai = latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai];
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
                WCM_LOAD(window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        }
        else
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Z_WINDOW_SIZE)
                window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE;
        #endif
    }

    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?(latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)):(latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)-(latticeXYZSizeC));
        const unsigned int ao = occ[ai];
        for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
            WCM_LOAD(window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        #else
        if (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)
        {
            const unsigned int ai = latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai];
            for (uint w=0, latticeOffset=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, latticeOffset+=latticeXYZSizeC, windowOffset+=MPD_Z_WINDOW_SIZE)
                WCM_LOAD(window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)], lattice[ai+latticeOffset], ao, w)
        }
        else
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Z_WINDOW_SIZE)
                window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE;
        #endif
    }
    #endif
}

/* performPropagation of word_diffusion_1d_dev.cu with the occupancy write-back (same particle order, same overflow handling). */
inline __device__ void performPropagationOcc(unsigned int * __restrict__ lattice, uint8_t * __restrict__ dstOcc, const unsigned int * __restrict__ window, const unsigned int * __restrict__ choices, const unsigned int latticeIndex, const unsigned int windowIndexMinus, const unsigned int windowIndex, const unsigned int windowIndexPlus, const unsigned int windowSize, unsigned int * __restrict__ siteOverflowList)
{
    int nextParticle=0;
    unsigned int newParticles[MPD_WORDS_PER_SITE*3];
    for(int i=0; i<MPD_WORDS_PER_SITE*3; i++)
        newParticles[i] = 0;

    #pragma unroll
    for(int w=0; w < MPD_WORDS_PER_SITE; w++)
        if(choices[windowIndex + w*windowSize] == MPD_MOVE_STAY)
            newParticles[nextParticle++] = window[windowIndex + w*windowSize];
    #pragma unroll
    for(int w=0; w < MPD_WORDS_PER_SITE; w++)
        if(choices[windowIndexMinus + w*windowSize] == MPD_MOVE_PLUS)
            newParticles[nextParticle++] = window[windowIndexMinus + w*windowSize];
    #pragma unroll
    for(int w=0; w < MPD_WORDS_PER_SITE; w++)
        if(choices[windowIndexPlus + w*windowSize] == MPD_MOVE_MINUS)
            newParticles[nextParticle++] = window[windowIndexPlus + w*windowSize];

    const unsigned int oldO = dstOcc[latticeIndex];
    const unsigned int newO = (nextParticle < MPD_WORDS_PER_SITE) ? (unsigned int)nextParticle : (unsigned int)MPD_WORDS_PER_SITE;
    const unsigned int wOut = (newO > oldO) ? newO : oldO;
    for(unsigned int w=0; w < MPD_WORDS_PER_SITE; w++, lattice += latticeXYZSizeC)
        if (w < wOut)
            lattice[latticeIndex] = newParticles[w];
    dstOcc[latticeIndex] = (uint8_t)newO;

    for (int i=MPD_PARTICLES_PER_SITE; i<nextParticle; i++)
    {
        int exceptionIndex = atomicAdd(siteOverflowList, 1);
        if (exceptionIndex < TUNE_MPD_MAX_PARTICLE_OVERFLOWS)
        {
            siteOverflowList[(exceptionIndex*2)+1]=latticeIndex;
            siteOverflowList[(exceptionIndex*2)+2]=newParticles[i];
        }
    }
}



/* ===================================================================================================================
 * occupancy-bounded window, choices and propagation. The window copies also record each window site's occupancy
 * (occWin) and store only its occupied slots; makeDiffusionChoicesOcc and performPropagationOcc2 loop over a site's occupied
 * slots only. Slots >= occupancy hold no particle, got the 'no move' choice, and contributed nothing to the propagation; the
 * random draw of a particle is a hash of (site, slot), so every particle gets the same choice and lands in the same slot.
 * =================================================================================================================== */


inline __device__ void makeDiffusionChoicesOcc(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, const uint8_t * __restrict__ occWin, uint8_t * __restrict__ choices, const unsigned int latticeIndex, const unsigned int windowIndexMinus, const unsigned int windowIndex, const unsigned int windowIndexPlus, const unsigned int windowSize, const unsigned long long timestepHash)
{
	const unsigned char site = sitesWindow[windowIndex];
	const unsigned char siteMinus = sitesWindow[windowIndexMinus];
	const unsigned char sitePlus = sitesWindow[windowIndexPlus];

	const int occHere = occWin[windowIndex];
	for(int w=0; w < occHere; w++, window = window+windowSize, choices = choices+windowSize)
	{
		// Set the default choice to none.
		choices[windowIndex]=0;

		// If there are no particles, we are done.
		if (window[windowIndex] > 0)
		{
			const unsigned int particle = window[windowIndex];

			if (particle > 0)
			{
				// Get the probability of moving plus and minus.
				float probMinus=lookupTransitionProbability(particle, site, siteMinus);
				float probPlus=lookupTransitionProbability(particle, site, sitePlus);

				// Get the random value and see which direction we should move in.
				float randomValue = getRandomHashFloat(latticeIndex, MPD_PARTICLE_COUNT_BITS, w, timestepHash);
				uint8_t c = (randomValue < probMinus)?(MPD_MOVE_MINUS):(MPD_MOVE_STAY);
				choices[windowIndex] = (randomValue >= 0.5f && randomValue < (probPlus+0.5f))?(MPD_MOVE_PLUS):(c);
			}
		}
	}
}

inline __device__ void makeXDiffusionChoicesOcc(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, const uint8_t * __restrict__ occWin, uint8_t * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex, const unsigned int blockXSize, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, latticeIndex, windowIndex-1, windowIndex, windowIndex+1, MPD_X_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowIndex >= blockXSize)
    {
        unsigned int apronLatticeIndex = (latticeXIndex>=blockXSize)?latticeIndex-blockXSize:latticeIndex+(latticeXSizeC-blockXSize);
        unsigned int apronWindowIndex = windowIndex-blockXSize;
        makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, apronLatticeIndex, apronWindowIndex, apronWindowIndex, apronWindowIndex+1, MPD_X_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeXIndex<latticeXSizeC-blockXSize)?latticeIndex+blockXSize:latticeIndex-latticeXSizeC+blockXSize;
        unsigned int apronWindowIndex = windowIndex+blockXSize;
        makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, apronLatticeIndex, apronWindowIndex-1, apronWindowIndex, apronWindowIndex, MPD_X_WINDOW_SIZE, timestepHash);
    }
}

inline __device__ void makeYDiffusionChoicesOcc(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, const uint8_t * __restrict__ occWin, uint8_t * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeYIndex, unsigned int windowIndex, const unsigned int windowYIndex, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, latticeIndex, windowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowYIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?latticeIndex-(latticeXSizeC*MPD_APRON_SIZE):latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)+(latticeXYSizeC);
        unsigned int apronWindowIndex = windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, apronLatticeIndex, apronWindowIndex, apronWindowIndex, apronWindowIndex+TUNE_MPD_Y_BLOCK_X_SIZE, MPD_Y_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowYIndex >= TUNE_MPD_Y_BLOCK_Y_SIZE)
    {
        unsigned int apronLatticeIndex = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?latticeIndex+(latticeXSizeC*MPD_APRON_SIZE):latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)-(latticeXYSizeC);
        unsigned int apronWindowIndex = windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, apronLatticeIndex,  apronWindowIndex-TUNE_MPD_Y_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex, MPD_Y_WINDOW_SIZE, timestepHash);
    }
}

inline __device__ void makeZDiffusionChoicesOcc(const unsigned int * __restrict__ window, const uint8_t * __restrict__ sitesWindow, const uint8_t * __restrict__ occWin, uint8_t * __restrict__ choices, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex, const unsigned long long timestepHash)
{
    // Calculate the diffusion choices for the segment index.
    makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, latticeIndex, windowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, windowIndex, windowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, timestepHash);

    // If this thread is one that needs to calculate choices for the leading apron, calculate them.
    if (windowZIndex < (2*MPD_APRON_SIZE))
    {
        unsigned int apronLatticeIndex = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE):latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)+(global_latticeXYZSizeC);
        unsigned int apronWindowIndex = windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, apronLatticeIndex, apronWindowIndex, apronWindowIndex, apronWindowIndex+TUNE_MPD_Z_BLOCK_X_SIZE, MPD_Z_WINDOW_SIZE, timestepHash);
    }

    // If this thread is one that needs to calculate choices for the trailing apron, calculate them.
    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        unsigned int apronLatticeIndex = (latticeZIndex<global_latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE):latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)-(global_latticeXYZSizeC);
        unsigned int apronWindowIndex = windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE);
        makeDiffusionChoicesOcc(window, sitesWindow, occWin, choices, apronLatticeIndex, apronWindowIndex-TUNE_MPD_Z_BLOCK_X_SIZE, apronWindowIndex, apronWindowIndex, MPD_Z_WINDOW_SIZE, timestepHash);

    }
}

// copies a site's o occupied slots into the window with the loads issued eight at a time before their shared stores
// (a load-store loop waits one global round trip per slot); same window contents.
inline __device__ void wcm_copy_slots(const unsigned int * __restrict__ lattice, const unsigned int li, unsigned int * window, const unsigned int wi, const unsigned int wstride, const unsigned int o)
{
    for (unsigned int w0=0; w0<o; w0+=8)
    {
        unsigned int v[8];
        #pragma unroll
        for (unsigned int k=0; k<8; k++) v[k] = (w0+k < o) ? lattice[li + (w0+k)*latticeXYZSizeC] : 0u;
        #pragma unroll
        for (unsigned int k=0; k<8; k++) if (w0+k < o) window[wi + (w0+k)*wstride] = v[k];
    }
}

// returns the OR of every occupancy it stored into occWin (0 = the window holds no particle).
inline __device__ unsigned int copyXWindowFromLatticeOcc2(const unsigned int bx, const unsigned int * lattice, const uint8_t * occ, unsigned int * window, uint8_t * occWin, const unsigned int latticeIndex, const unsigned int latticeXIndex, const unsigned int windowIndex)
{
    unsigned int wcm_or = 0;
    if (latticeXIndex < latticeXSizeC)
    {
        const unsigned int o = occ[latticeIndex]; occWin[windowIndex] = (uint8_t)o; wcm_or |= o;
        wcm_copy_slots(lattice, latticeIndex, window, windowIndex, MPD_X_WINDOW_SIZE, o);

        #if MPD_APRON_SIZE > 0
        int threadBlockWidth = ((bx+1)*blockDim.x <= latticeXSizeC)?(blockDim.x):(latticeXSizeC-(bx*blockDim.x));

        if (windowIndex >= threadBlockWidth)
        {
            #if defined MPD_BOUNDARY_PERIODIC
            const unsigned int ai = (latticeXIndex>=threadBlockWidth)?(latticeIndex-threadBlockWidth):(latticeIndex-threadBlockWidth+latticeXSizeC);
            const unsigned int ao = occ[ai]; occWin[windowIndex-threadBlockWidth] = (uint8_t)ao; wcm_or |= ao;
            wcm_copy_slots(lattice, ai, window, windowIndex-threadBlockWidth, MPD_X_WINDOW_SIZE, ao);
            #else
            if (latticeXIndex>=threadBlockWidth)
            {
                const unsigned int ai = latticeIndex-threadBlockWidth; const unsigned int ao = occ[ai]; occWin[windowIndex-threadBlockWidth] = (uint8_t)ao; wcm_or |= ao;
                wcm_copy_slots(lattice, ai, window, windowIndex-threadBlockWidth, MPD_X_WINDOW_SIZE, ao);
            }
            else
                { occWin[windowIndex-threadBlockWidth] = (uint8_t)MPD_WORDS_PER_SITE; wcm_or |= MPD_WORDS_PER_SITE;
                for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_X_WINDOW_SIZE)
                    window[windowIndex+windowOffset-threadBlockWidth] = MPD_BOUNDARY_PARTICLE_VALUE; }
            #endif
        }

        if (windowIndex < 2*MPD_APRON_SIZE)
        {
            #if defined MPD_BOUNDARY_PERIODIC
            const unsigned int ai = (latticeXIndex<latticeXSizeC-threadBlockWidth)?(latticeIndex+threadBlockWidth):(latticeIndex+threadBlockWidth-latticeXSizeC);
            const unsigned int ao = occ[ai]; occWin[windowIndex+threadBlockWidth] = (uint8_t)ao; wcm_or |= ao;
            wcm_copy_slots(lattice, ai, window, windowIndex+threadBlockWidth, MPD_X_WINDOW_SIZE, ao);
            #else
            if (latticeXIndex<latticeXSizeC-threadBlockWidth)
            {
                const unsigned int ai = latticeIndex+threadBlockWidth; const unsigned int ao = occ[ai]; occWin[windowIndex+threadBlockWidth] = (uint8_t)ao; wcm_or |= ao;
                wcm_copy_slots(lattice, ai, window, windowIndex+threadBlockWidth, MPD_X_WINDOW_SIZE, ao);
            }
            else
                { occWin[windowIndex+threadBlockWidth] = (uint8_t)MPD_WORDS_PER_SITE; wcm_or |= MPD_WORDS_PER_SITE;
                for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_X_WINDOW_SIZE)
                    window[windowIndex+windowOffset+threadBlockWidth] = MPD_BOUNDARY_PARTICLE_VALUE; }
            #endif
        }
        #endif
    }
    return wcm_or;
}

// returns the OR of every occupancy it stored into occWin (0 = the window holds no particle).
inline __device__ unsigned int copyYWindowFromLatticeOcc2(const unsigned int* lattice, const uint8_t * occ, unsigned int* window, uint8_t * occWin, const unsigned int latticeIndex, const unsigned int latticeYIndex, const unsigned int windowIndex, const unsigned int windowYIndex)
{
    unsigned int wcm_or = 0;
    const unsigned int o = occ[latticeIndex]; occWin[windowIndex] = (uint8_t)o; wcm_or |= o;
    wcm_copy_slots(lattice, latticeIndex, window, windowIndex, MPD_Y_WINDOW_SIZE, o);

    #if MPD_APRON_SIZE > 0
    if (windowYIndex < 2*MPD_APRON_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)?(latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)):(latticeIndex-(latticeXSizeC*MPD_APRON_SIZE)+(latticeXSizeC*latticeYSizeC));
        const unsigned int ao = occ[ai]; occWin[windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
        wcm_copy_slots(lattice, ai, window, windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Y_WINDOW_SIZE, ao);
        #else
        if (latticeYIndex>=TUNE_MPD_Y_BLOCK_Y_SIZE)
        {
            const unsigned int ai = latticeIndex-(latticeXSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai]; occWin[windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
            wcm_copy_slots(lattice, ai, window, windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Y_WINDOW_SIZE, ao);
        }
        else
            { occWin[windowIndex-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)MPD_WORDS_PER_SITE; wcm_or |= MPD_WORDS_PER_SITE;
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Y_WINDOW_SIZE)
                window[windowIndex+windowOffset-(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE; }
        #endif
    }

    if (windowYIndex >= TUNE_MPD_Y_BLOCK_Y_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)?(latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)):(latticeIndex+(latticeXSizeC*MPD_APRON_SIZE)-(latticeXSizeC*latticeYSizeC));
        const unsigned int ao = occ[ai]; occWin[windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
        wcm_copy_slots(lattice, ai, window, windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Y_WINDOW_SIZE, ao);
        #else
        if (latticeYIndex<latticeYSizeC-TUNE_MPD_Y_BLOCK_Y_SIZE)
        {
            const unsigned int ai = latticeIndex+(latticeXSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai]; occWin[windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
            wcm_copy_slots(lattice, ai, window, windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Y_WINDOW_SIZE, ao);
        }
        else
            { occWin[windowIndex+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)MPD_WORDS_PER_SITE; wcm_or |= MPD_WORDS_PER_SITE;
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Y_WINDOW_SIZE)
                window[windowIndex+windowOffset+(TUNE_MPD_Y_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE; }
        #endif
    }
    #endif
    return wcm_or;
}

// returns the OR of every occupancy it stored into occWin (0 = the window holds no particle).
inline __device__ unsigned int copyZWindowFromLatticeOcc2(const unsigned int* lattice, const uint8_t * occ, unsigned int* window, uint8_t * occWin, const unsigned int latticeIndex, const unsigned int latticeZIndex, const unsigned int windowIndex, const unsigned int windowZIndex)
{
    unsigned int wcm_or = 0;
    const unsigned int o = occ[latticeIndex]; occWin[windowIndex] = (uint8_t)o; wcm_or |= o;
    wcm_copy_slots(lattice, latticeIndex, window, windowIndex, MPD_Z_WINDOW_SIZE, o);

    #if MPD_APRON_SIZE > 0
    if (windowZIndex < 2*MPD_APRON_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)?(latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)):(latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE)+(latticeXYZSizeC));
        const unsigned int ao = occ[ai]; occWin[windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
        wcm_copy_slots(lattice, ai, window, windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Z_WINDOW_SIZE, ao);
        #else
        if (latticeZIndex>=TUNE_MPD_Z_BLOCK_Z_SIZE)
        {
            const unsigned int ai = latticeIndex-(latticeXYSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai]; occWin[windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
            wcm_copy_slots(lattice, ai, window, windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Z_WINDOW_SIZE, ao);
        }
        else
            { occWin[windowIndex-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)MPD_WORDS_PER_SITE; wcm_or |= MPD_WORDS_PER_SITE;
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Z_WINDOW_SIZE)
                window[windowIndex+windowOffset-(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE; }
        #endif
    }

    if (windowZIndex >= TUNE_MPD_Z_BLOCK_Z_SIZE)
    {
        #if defined MPD_BOUNDARY_PERIODIC
        const unsigned int ai = (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)?(latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)):(latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE)-(latticeXYZSizeC));
        const unsigned int ao = occ[ai]; occWin[windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
        wcm_copy_slots(lattice, ai, window, windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Z_WINDOW_SIZE, ao);
        #else
        if (latticeZIndex<latticeZSizeC-TUNE_MPD_Z_BLOCK_Z_SIZE)
        {
            const unsigned int ai = latticeIndex+(latticeXYSizeC*MPD_APRON_SIZE); const unsigned int ao = occ[ai]; occWin[windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)ao; wcm_or |= ao;
            wcm_copy_slots(lattice, ai, window, windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE), MPD_Z_WINDOW_SIZE, ao);
        }
        else
            { occWin[windowIndex+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = (uint8_t)MPD_WORDS_PER_SITE; wcm_or |= MPD_WORDS_PER_SITE;
            for (uint w=0, windowOffset=0; w<MPD_WORDS_PER_SITE; w++, windowOffset+=MPD_Z_WINDOW_SIZE)
                window[windowIndex+windowOffset+(TUNE_MPD_Z_BLOCK_X_SIZE*MPD_APRON_SIZE)] = MPD_BOUNDARY_PARTICLE_VALUE; }
        #endif
    }
    #endif
    return wcm_or;
}

inline __device__ void performPropagationOcc2(unsigned int * __restrict__ lattice, uint8_t * __restrict__ dstOcc, const unsigned int * __restrict__ window, const uint8_t * __restrict__ occWin, const uint8_t * __restrict__ choices, const unsigned int latticeIndex, const unsigned int windowIndexMinus, const unsigned int windowIndex, const unsigned int windowIndexPlus, const unsigned int windowSize, unsigned int * __restrict__ siteOverflowList, const unsigned int oldO)
{
    // each arriving particle is stored straight into the next destination slot, in the order the former local
    // newParticles[48] array collected them (stay, from minus, from plus); the thread keeps no local-memory array. Slots from the
    // new occupancy up to the old one are zeroed, particles beyond MPD_PARTICLES_PER_SITE go to the overflow list in the same
    // order: the destination lattice, occupancy and overflow entries are the same as before.
    unsigned int n = 0;
    #define WCM_PUT(p) { const unsigned int wcm_p = (p); \
        if (n < MPD_WORDS_PER_SITE) lattice[latticeIndex + n*latticeXYZSizeC] = wcm_p; \
        else { int exceptionIndex = atomicAdd(siteOverflowList, 1); \
               if (exceptionIndex < TUNE_MPD_MAX_PARTICLE_OVERFLOWS) { siteOverflowList[(exceptionIndex*2)+1]=latticeIndex; siteOverflowList[(exceptionIndex*2)+2]=wcm_p; } } \
        n++; }
    for(int w=0, wn=occWin[windowIndex]; w < wn; w++)
        if(choices[windowIndex + w*windowSize] == MPD_MOVE_STAY) WCM_PUT(window[windowIndex + w*windowSize])
    for(int w=0, wn=occWin[windowIndexMinus]; w < wn; w++)
        if(choices[windowIndexMinus + w*windowSize] == MPD_MOVE_PLUS) WCM_PUT(window[windowIndexMinus + w*windowSize])
    for(int w=0, wn=occWin[windowIndexPlus]; w < wn; w++)
        if(choices[windowIndexPlus + w*windowSize] == MPD_MOVE_MINUS) WCM_PUT(window[windowIndexPlus + w*windowSize])
    #undef WCM_PUT

    // oldO (the destination site's previous occupancy) is loaded by the caller at kernel start.
    const unsigned int newO = (n < MPD_WORDS_PER_SITE) ? n : (unsigned int)MPD_WORDS_PER_SITE;
    for(unsigned int w=newO; w < oldO; w++)
        lattice[latticeIndex + w*latticeXYZSizeC] = 0;
    dstOcc[latticeIndex] = (uint8_t)newO;
}
