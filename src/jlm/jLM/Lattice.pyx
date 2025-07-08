# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
    
"""Cythonized lattice manipulation functions"""

import numpy as npy
cimport numpy as cnpy
from libc.stdlib cimport *
from posix.stdlib cimport *
from libc.time cimport *
from libc.stdint cimport *

ctypedef fused plattice_type:
    uint8_t
    uint32_t

cdef void latticeHistogram3(int[:,:,:] target, unsigned char[:,:,:,:] src, int[:] idxs) nogil:
    cdef int nx,ny,nz,np,ni,x,y,z,p,i
    nx = src.shape[0]
    ny = src.shape[1]
    nz = src.shape[2]
    np = src.shape[3]
    ni = idxs.shape[0]

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for p in range(np):
                    for i in range(ni):
                        if src[x,y,z,p] == idxs[i]:
                            target[x,y,z] += 1
                            break

cdef void latticeHistogram2(int[:,:] target, unsigned char[:,:,:,:] src, int[:] idxs, int intDim) nogil:
    cdef int nu,nv,nx,ny,nz,np,ni,x,y,z,p,i
    nx = src.shape[0]
    ny = src.shape[1]
    nz = src.shape[2]
    np = src.shape[3]
    ni = idxs.shape[0]

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for p in range(np):
                    for i in range(ni):
                        if src[x,y,z,p] == idxs[i]:
                            if intDim == 0:
                                target[y,z] += 1
                            elif intDim == 1:
                                target[z,x] += 1
                            else:
                                target[x,y] += 1
                            break

cdef void latticeHistogram1(int[:] target, unsigned char[:,:,:,:] src, int[:] idxs, int keepDim) nogil:
    cdef int nu,nv,nx,ny,nz,np,ni,x,y,z,p,i
    nx = src.shape[0]
    ny = src.shape[1]
    nz = src.shape[2]
    np = src.shape[3]
    ni = idxs.shape[0]

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for p in range(np):
                    for i in range(ni):
                        if src[x,y,z,p] == idxs[i]:
                            if keepDim == 0:
                                target[x] += 1
                            elif keepDim == 1:
                                target[y] += 1
                            else:
                                target[z] += 1
                            break

def latticeHistogram(unsigned char[:,:,:,:] src, int[:] idxs, axes=None, target=None):
    """Compute histogram of particle types.

    Args:
        src (:py:class:`numpy.ndarray`):
            Result array, contains the number of particles at each lattice 
            site in `idxs`. (shape=(nx,ny,nz,pps), dtype=uint8)
        idxs (:py:class:`numpy.ndarray`):
            List of particle types to bin (shape=(nparticles,), dtype=int32)
            
    Keyword Args:
        axes (list):
            Axis indexes to sum out (shape=(naxes,), dtype=int32)
        target:
            Output array. None or an :py:class:`numpy.ndarray` of the proper
            size

    Returns:
        :py:class:`numpy.ndarray`:
            New ndarray or reference to target
    """
    cdef int u,v,dim,histDim

    if axes is not None and len(axes)>0:
        if len(axes) == 1:
            dim = 2
            dims = [0,1,2]
            u = dims[(axes[0]+1)%3]
            v = dims[(axes[0]+2)%3]
            histDim = int(axes[0])
            targetShape = (src.shape[u],src.shape[v])
        elif len(axes) == 2:
            dim = 1
            u = next(iter(set([0,1,2])-set(axes)))
            histDim = int(u)
            targetShape = (src.shape[u],)
        elif len(axes) == 3:
            targetShape = (src.shape[0], src.shape[1], src.shape[2])
            dim = 3
        else:
            raise RuntimeError("Wrong axes args")
    else:
        targetShape = (src.shape[0], src.shape[1], src.shape[2])
        dim = 3

    if target is None:
        target = npy.zeros(targetShape,dtype='i4')

    if dim == 1:
        latticeHistogram1(target, src, idxs, histDim)
    elif dim==2:
        latticeHistogram2(target, src, idxs, histDim)
    else:
        latticeHistogram3(target, src, idxs)

    return target



def countParticlesInLattice(plattice_type[:,:,:,:,] pLattice, unsigned char[:,:,:] sLattice, int[:] spIdList, int[:] regIdList):
    """Count particles.

    Given a list of region ids and species ids, count all particles of those 
    types within those regions.

    Args:
        pLattice (:py:class:`numpy.ndarray`):
            Particle lattice (HDF5 format, shape=(nx,ny,nz,pps))
        sLattice (:py:class:`numpy.ndarray`):
            Site lattice (shape=(nx,ny,nz))
        spIdList (:py:class:`numpy.ndarray`):
            List of particle ids to count (shape=(nsps,))
        regIdList (:py:class:`numpy.ndarray`):
            List of region ids to consider (shape=(nreg,))
    Returns:
        int: Count
    """
    cdef int nx,ny,nz,np,ns,nr,x,y,z,p,s,r,ct
    nx = pLattice.shape[0]
    ny = pLattice.shape[1]
    nz = pLattice.shape[2]
    np = pLattice.shape[3]
    ns = spIdList.shape[0]
    nr = regIdList.shape[0]

    assert nx == sLattice.shape[0]
    assert ny == sLattice.shape[1]
    assert nz == sLattice.shape[2]

    ct = 0

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for r in range(nr):
                    if sLattice[x,y,z] == regIdList[r]:
                        for p in range(np):
                            for s in range(ns):
                                if pLattice[x,y,z,p] == spIdList[s]:
                                    ct += 1
    return ct

def hdf2lmRep(plattice_type[:,:,:,:,:] plattice1, plattice_type[:,:,:,:] plattice0, txTable=None):
    """Convert an HDF5 particle lattice representation to LM native format

    Args:
        plattice1 (:py:class:`numpy.ndarray`):
            Output LM format (dtype=uint8/uint32, shape=(nw,nz,ny,nx,np))
        plattice0 (:py:class:`numpy.ndarray`):
            Input HDF5 format (dtype=uint8/uint32, shape=(nx,ny,nz,nw*np))
    Keyword Args:
        txTable (:py:class:`numpy.ndarray`):
            Optional particle id translation table: lmNative = txTable[h5Native]
    """
    cdef int nx,ny,nz,np,x,y,z,p,w
    nw = plattice1.shape[0]
    nz = plattice1.shape[1]
    ny = plattice1.shape[2]
    nx = plattice1.shape[3]
    np = plattice1.shape[4]

    cdef unsigned char[:] table

    assert nw*np == plattice0.shape[3]
    assert nz == plattice0.shape[2]
    assert ny == plattice0.shape[1]
    assert nx == plattice0.shape[0]

    if txTable is None:
        for w in range(nw):
            for z in range(nz):
                for y in range(ny):
                    for x in range(nx):
                        for p in range(np):
                            plattice1[w,z,y,x,p] = plattice0[x,y,z,p+np*w]
    else:
        table = txTable
        for w in range(nw):
            for z in range(nz):
                for y in range(ny):
                    for x in range(nx):
                        for p in range(np):
                            plattice1[w,z,y,x,p] = table[plattice0[x,y,z,p+np*w]]


def latticeStatsAll_h5fmt(plattice_type[:,:,:,:] plattice, unsigned char[:,:,:] slattice):
    """Compute number of particles in a particular site type, and the number of those site types

    Args:
        plattice (:py:class:`numpy.ndarray`):
            Particle lattice (H5 format: shape=(nx,ny,nz,pps), dtype=uint8/uint32)
        slattice (:py:class:`numpy.ndarray`):
            Site lattice (HDF5 format: shape=(nx,ny,nz), dtype=uint8)

    Returns:
        (:py:class:`numpy.ndarray`,:py:class:`numpy.ndarray`): 
            Number of particles, number of sites
    """
    cdef int nx,ny,nz,np,x,y,z,p

    nx = plattice.shape[0]
    ny = plattice.shape[1]
    nz = plattice.shape[2]
    np = plattice.shape[3]

    assert nx == slattice.shape[0]
    assert ny == slattice.shape[1]
    assert nz == slattice.shape[2]

    cdef cnpy.ndarray nd_pCount, nd_sCount
    cdef int64_t[:,:] pCount
    cdef int64_t[:] sCount
    np_pCount = npy.zeros((16384,16384), dtype=npy.int64)
    np_sCount = npy.zeros(16384, dtype=npy.int64)
    pCount = np_pCount
    sCount = np_sCount

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                sCount[slattice[x,y,z]] += 1
                for p in range(np):
                    pCount[slattice[x,y,z], plattice[x,y,z,p]] += 1

    return np_pCount, np_sCount

def latticeStatsAll(plattice_type[:,:,:,:,:] plattice, unsigned char[:,:,:] slattice):
    """ Compute number of particles in a particular site type, and the number of those site types

    Args:
        plattice (:py:class:`numpy.ndarray`):
            Particle lattice (shape=(nw,nz,ny,nx,np), dtype=uint8/uint32)
        slattice (:py:class:`numpy.ndarray`):
            Site lattice (HDF5 format: shape=(nz,ny,nx), dtype=uint8)

    Returns:
        (:py:class:`numpy.ndarray`,:py:class:`numpy.ndarray`): 
            Number of particles, number of sites
    """
    cdef int nx,ny,nz,np,nsts,x,y,z,p,st
    nw = plattice.shape[0]
    nx = plattice.shape[1]
    ny = plattice.shape[2]
    nz = plattice.shape[3]
    np = plattice.shape[4]

    cdef cnpy.ndarray nd_pCount, nd_sCount
    cdef int64_t[:,:] pCount
    cdef int64_t[:] sCount
    np_pCount = npy.zeros((16384,16384), dtype=npy.int64)
    np_sCount = npy.zeros(16384, dtype=npy.int64)
    pCount = np_pCount
    sCount = np_sCount

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                sCount[slattice[x,y,z]] += 1
                for w in range(nw):
                    for p in range(np):
                        pCount[slattice[x,y,z], plattice[w,x,y,z,p]] += 1

    return np_pCount, np_sCount



def countSites(unsigned char[:,:,:] slattice):
    """Count number of subvolumes

    Args:
        slattice (:py:class:`numpy.ndarray`):
            Site lattice

    Returns:
        :py:class:`numpy.ndarray`:
            Ndarray of shape=(256,) counting the number of site of each type
    """
    cdef cnpy.ndarray np_count
    cdef int64_t[:] count
    cdef int nx,ny,nz,x,y,z
    nx = slattice.shape[0]
    ny = slattice.shape[1]
    nz = slattice.shape[2]

    np_count = npy.zeros(16384, dtype=npy.int64)
    count = np_count

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                count[slattice[x,y,z]] += 1

    return np_count


def appendParticle(plattice_type[:,:,:,:,:] pLattice, uint64_t x, uint64_t y, uint64_t z, uint64_t sp):
    """Add a particle to a lattice site

    Args:
        plattice (:py:class:`numpy.ndarray`):
            Particle lattice (H5 format: shape=(nw,nz,ny,nz,np), dtype=uint8/uint32)
        x (int):
            x coordinate
        y (int):
            y coordinate
        z (int):
            z coordinate
        sp (uint64):
            Species type

    Returns:
        int: 1 if added, 0 if site full
    """
    cdef uint64_t Nw, Np, w, p
    Nw = pLattice.shape[0]
    Np = pLattice.shape[4]
    for w in range(Nw):
        for p in range(Np):
            if pLattice[w,z,y,x,p] == 0:
                pLattice[w,z,y,x,p] = sp
                return 1
    else:
        return 0


# random int w/o modulo bias
cdef int64_t randIntLessThan(int64_t n):
    cdef int x = random()
    while x >= RAND_MAX - (RAND_MAX % n):
          x = random()
    return x % n

cdef uint32_t chooseParticle(double* rs, int n):
    cdef double acc = 0
    cdef double r = drand48()
    cdef int i = 0

    for i in range(n):
        acc += rs[i]
        if r < acc:
            return i

def populateLattice(pLattice, sLattice, particleDensity, exact=True, mask=None):
    r"""Add particles to lattice at random locations.

    The `particleDensityView` argument is a 
    :math:`N_{siteTypes}\times N_{particleTypes}` array of slot occupancy 
    probabilities. It describes the probability to find a particle type 
    (including empty==0) at a particular particle slot.

    Args:
        pLattice (:py:class:`numpy.ndarray`):
            Particle lattice (H5 format, shape=(nw,nz,ny,nx,np), dtype=uint8/uint32)
        sLattice (:py:class:`numpy.ndarray`):
            Site lattice (H5 format, shape=(nz,ny,nx), dtype=uint8)
        particleDensity (:py:class:`numpy.ndarray`):
            Particle probability, (shape=(Nsites,Nsps), dtype=float64)
    Keyword Args:
        exact (bool):
            If `True`, make a second pass that corrects the particle counts 
            exactly
        mask (:py:class:`numpy.ndarray`):
            Boolean mask to set region of lattice to populate, `None` 
            (default) populates the entire lattice
    """
    cdef cnpy.ndarray mk
    cdef int ex

    if mask is None:
        mk = npy.zeros_like(sLattice)+1
    elif mask.dtype is not npy.uint8:
        mk = mask.astype(npy.uint8)
    else:
        mk = mask

    if exact:
        ex = 1
    else:
        ex = 0

    if npy.issubdtype(pLattice.dtype, npy.uint32):
        populateLattice_[uint32_t](pLattice, sLattice, mk, particleDensity, ex)
    else:
        populateLattice_[uint8_t](pLattice, sLattice, mk, particleDensity, ex)


cdef void populateLattice_(plattice_type[:,:,:,:,:] pLattice, uint8_t[:,:,:] sLattice, uint8_t[:,:,:] mask,
                           double[:,:] particleDensityView, int exact):
    cdef int64_t sp,st
    cdef int64_t w,x,y,z,p
    cdef int64_t spIdx,taIdx,laIdx,ltIdx
    cdef int64_t lutOffset, lutCount, lutEmptyOffset, lutEmptyCount

    cdef int64_t Nsts = particleDensityView.shape[0]
    cdef int64_t Nsps = particleDensityView.shape[1]
    cdef int64_t Ntab = Nsps*Nsts

    cdef int64_t Nw = pLattice.shape[0]
    cdef int64_t Nz = pLattice.shape[1]
    cdef int64_t Ny = pLattice.shape[2]
    cdef int64_t Nx = pLattice.shape[3]
    cdef int64_t Np = pLattice.shape[4]
    cdef int64_t siteSize = Nw*Np

    cdef plattice_type* pLatticeFlat = <plattice_type*> &(pLattice[0,0,0,0,0])
    cdef double* particleDensity= <double*> &(particleDensityView[0,0])

    cdef uint32_t* siteTmp = <uint32_t*> calloc(siteSize,sizeof(uint32_t))
    cdef int64_t* realCounts = <int64_t*> calloc(Ntab,sizeof(int64_t))
    cdef int64_t* offset = <int64_t*> calloc(Ntab,sizeof(int64_t))
    cdef int64_t* cursor = <int64_t*> calloc(Ntab,sizeof(int64_t))
    cdef int64_t* lookup = <int64_t*> calloc(Nw*Nz*Ny*Nx*Np,sizeof(int64_t))
    cdef int64_t* countsDiff = <int64_t*> calloc(Ntab,sizeof(int64_t))

    cdef int64_t* counts = <int64_t*> calloc(Ntab,sizeof(int64_t))
    cdef int64_t* siteVolume = <int64_t*> calloc(Nsts,sizeof(int64_t))

    srand48(time(<time_t *>0))
    srandom(time(<time_t *>0))

    # place particles using probabilities

    for w in range(Nw):
        for z in range(Nz):
            for y in range(Ny):
                for x in range(Nx):
                    for p in range(Np):
                        st = sLattice[z,y,x]
                        if mask[z,y,x] != 0:
                            sp = chooseParticle(particleDensity+Nsps*st, Nsps)
                            pLattice[w,z,y,x,p] = sp
                            realCounts[sp+Nsps*st] += 1
                            siteVolume[st] += 1

    if exact == 1:
        for st in range(Nsts):
            for sp in range(Nsps):
                counts[sp+Nsps*st] = int(0.5+siteVolume[st]*particleDensity[sp+Nsps*st])

        # compute placement error
        for taIdx in range(Ntab):
            countsDiff[taIdx] = realCounts[taIdx] - counts[taIdx]


        # get offsets for each (species, site) into lookup table
        for taIdx in range(1, Ntab):
            offset[taIdx] = offset[taIdx-1]+realCounts[taIdx-1]
            cursor[taIdx] = offset[taIdx]

        # construct lookup table
        for w in range(Nw):
            for z in range(Nz):
                for y in range(Ny):
                    for x in range(Nx):
                        st = sLattice[z,y,x]
                        if mask[z,y,x] != 0:
                            for p in range(Np):
                                sp = pLattice[w,z,y,x,p]
                                taIdx = sp + Nsps*st
                                lookup[cursor[taIdx]] = p+Np*(x+Nx*(y+Ny*(z+Nz*w)))
                                cursor[taIdx] += 1

        # correct placement error
        for st in range(Nsts):
            lutEmptyOffset = offset[st*Nsps]     # species index 0
            lutEmptyCount = realCounts[st*Nsps] # species index 0

            for sp in range(1,Nsps): # ignore holes
                taIdx = sp+Nsps*st
                lutOffset = offset[taIdx]
                lutCount = realCounts[taIdx]

                if countsDiff[taIdx] > 0: # over count
                    while countsDiff[taIdx] > 0:
                        laIdx = -1
                        while laIdx < 0:
                            ltIdx = randIntLessThan(lutCount) + lutOffset
                            laIdx = lookup[ltIdx]

                        pLatticeFlat[laIdx] = 0 # delete the particle
                        countsDiff[taIdx] -= 1
                        lookup[ltIdx] = -1 # ignore lattice index later


                elif countsDiff[taIdx] < 0: # under count
                    while countsDiff[taIdx] != 0:
                        laIdx = -1
                        while laIdx < 0:
                            ltIdx = randIntLessThan(lutEmptyCount) + lutEmptyOffset
                            laIdx = lookup[ltIdx]

                        pLatticeFlat[laIdx] = sp # add a particle
                        countsDiff[taIdx] += 1
                        lookup[ltIdx] = -1 # ignore lattice index later

    # compact each site
    for z in range(Nz):
        for y in range(Ny):
            for x in range(Nx):
                if mask[z,y,x] != 0:
                    for w in range(Nw):
                        for p in range(Np):
                            siteTmp[p + Np*w] = pLattice[w,z,y,x,p]

                    p = 0
                    for spIdx in range(siteSize):
                        sp = siteTmp[spIdx]
                        if siteTmp[spIdx] != 0:
                            siteTmp[spIdx] = 0
                            siteTmp[p] = sp
                            p += 1

                    for w in range(Nw):
                        for p in range(Np):
                            pLattice[w,z,y,x,p] = siteTmp[p + Np*w]

    free(siteTmp)
    free(realCounts)
    free(offset)
    free(cursor)
    free(lookup)
    free(countsDiff)

cdef int makeQuads(int du, int u, int[:,:] mask, int[:,:] uFaceCount,
              int[:] quadsDir, int[:] quadsU, int[:] quadsV, int[:] quadsW, 
              int[:] quadsDV, int[:] quadsDW, int[:] quadsCCW):
    # build mesh greedily using lexicographic ordering
    cdef int nv, nw, v, w, v1, w1, lenW, lenV, nq
    nq = 0

    nv = mask.shape[0]
    nw = mask.shape[1]


    for v in range(nv):
        w = 0
        while w < nw:
            if mask[v,w]:
                for w1 in range(w+1,nw): # measure w direction
                    if not mask[v,w1]:
                        break
                else:
                    w1 = nw
                lenW = w1-w

                vdirDone = False
                for v1 in range(v+1,nv): # measure v direction 
                    for w1 in range(w,w+lenW): 
                        if not mask[v1,w1]:
                            vdirDone = True
                            break

                    if vdirDone:
                        break
                else:
                    v1 = nv

                lenV = v1-v

                quadsDir[nq] = du 
                quadsU[nq] = u+1 
                quadsV[nq] = v 
                quadsW[nq] = w 
                quadsDV[nq] = lenV 
                quadsDW[nq] = lenW
                quadsCCW[nq] = uFaceCount[v,w]%2 # first quad has uFaceCount=1, so vertex order is CW
                nq += 1

                # remove quad from mask
                for v1 in range(v,v+lenV):
                    for w1 in range(w,w+lenW):
                        mask[v1,w1] = 0

                w += lenW
            else:
                w += 1

    return nq

cdef int makeIdx(int[:] vert, int* nverts, int[:,:] verts, int[:,:,:] lookup):
    cdef int idx = lookup[vert[0],vert[1],vert[2]]
    if idx == -1:
        verts[nverts[0],:] = vert[:]
        idx = nverts[0]
        nverts[0] +=  1
        lookup[vert[0],vert[1],vert[2]] = idx

    return idx


def greedyMesh(unsigned char[:,:,:] binaryLattice):
    """Generate a simplified triangle mesh from a binary lattice.

    Args:
        binaryLattice (:py:class:`numpy.ndarray`):
            Site lattice. One if site occupied, zero otherwise. 
            (shape=(nx,ny,nz), dtype=uint8)

    Returns:
        (:py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`):
            Vertices (shape=(nverts,3), dtype=int32), and faces, each 
            element which indexes into verticies. (shape=(nfaces, 3), 
            dtype=uint8).

    Todo:
       The correct winding order is not always chosen.
    """
    cdef int du,dv,dw,u,v,w,nu,nv,nw,qu,qv,qw,ls,i,nq
    cdef int p0, p1
    cdef int[:] X, Q
    cdef int[:,:] mask, uFaceCount
    cdef int[:] quadsU, quadsV, quadsW, quadsDV, quadsDW, quadsCCW, quadsDir

    ls = binaryLattice.shape[0]*binaryLattice.shape[1]*binaryLattice.shape[2]

    X = npy.zeros(3,dtype=npy.int32)
    Q = npy.zeros(3,dtype=npy.int32)
    quadsDir = npy.zeros(ls,dtype=npy.int32)
    quadsU = npy.zeros(ls,dtype=npy.int32)
    quadsV = npy.zeros(ls,dtype=npy.int32)
    quadsW = npy.zeros(ls,dtype=npy.int32)
    quadsDV = npy.zeros(ls,dtype=npy.int32)
    quadsDW = npy.zeros(ls,dtype=npy.int32)
    quadsCCW= npy.zeros(ls,dtype=npy.int32)
    quadsU[:] = -1
    nq = 0

    # loop over dimensions
    for du in range(3):
        dv = (du+1)%3
        dw = (du+2)%3
        nu = binaryLattice.shape[du]
        nv = binaryLattice.shape[dv]
        nw = binaryLattice.shape[dw]

        mask = npy.zeros((nv,nw),dtype=npy.int32)
        uFaceCount = npy.zeros((nv,nw),dtype=npy.int32)

        # displacment to next plane
        Q[du] = 1
        Q[dv] = 0
        Q[dw] = 0

        # loop over planes orthogonal to du
        for u in range(-1,nu):
            X[du] = u

            # generate mask
            for v in range(nv):
                X[dv] = v
                for w in range(nw):
                    X[dw] = w

                    p0 = 0
                    p1 = 0

                    if 0 <= u:
                        p0 = binaryLattice[X[0], X[1], X[2]]

                    if u < nu-1:
                        p1 = binaryLattice[X[0]+Q[0], X[1]+Q[1], X[2]+Q[2]]

                    mask[v,w] = p0 != p1

                    uFaceCount[v,w] += mask[v,w]

            nq += makeQuads(du, u, mask, uFaceCount, quadsDir[nq:],
                            quadsU[nq:], quadsV[nq:], quadsW[nq:], 
                            quadsDV[nq:], quadsDW[nq:], quadsCCW[nq:])

    cdef int nverts, nfaces
    cdef int[:] lu_shape, A, B, C, D
    cdef int[:,:] verts, faces
    cdef int[:,:,:] lookup
    cdef cnpy.ndarray np_faces, np_verts

    np_faces = npy.zeros((2*nq,3), dtype=npy.int32)
    np_verts = npy.zeros((4*nq,3), dtype=npy.int32)
    faces=np_faces
    verts=np_verts
    
    A = npy.zeros(3, dtype=npy.int32)
    B = npy.zeros(3, dtype=npy.int32)
    C = npy.zeros(3, dtype=npy.int32)
    D = npy.zeros(3, dtype=npy.int32)
    lu_shape = npy.zeros(3,dtype=npy.int32)
    lu_shape[:] = -1

    for i in range(nq):
        du = quadsDir[i]
        dv = (du+1)%3
        dw = (du+2)%3
        lu_shape[du] = max(lu_shape[du], 1+quadsU[i])
        lu_shape[dv] = max(lu_shape[dv], 1+quadsV[i])
        lu_shape[dw] = max(lu_shape[dw], 1+quadsW[i])

    lookup = npy.zeros(lu_shape, dtype=npy.int32)
    lookup[:,:,:] = -1
    verts[:,:] = -1
    faces[:,:] = -1
    nverts = 0
    nfaces = 0

    for i in range(nq):
        du = quadsDir[i]
        dv = (du+1)%3
        dw = (du+2)%3
        A[du] = quadsU[i]
        A[dv] = quadsV[i]
        A[dw] = quadsW[i]

        B[du] = quadsU[i]
        B[dv] = quadsV[i] + quadsDV[i]
        B[dw] = quadsW[i]

        C[du] = quadsU[i]
        C[dv] = quadsV[i] + quadsDV[i]
        C[dw] = quadsW[i] + quadsDW[i]

        D[du] = quadsU[i]
        D[dv] = quadsV[i]
        D[dw] = quadsW[i] + quadsDW[i]

        if quadsCCW[i] == 1:
            faces[nfaces,0] = makeIdx(A, &nverts, verts, lookup)
            faces[nfaces,1] = makeIdx(C, &nverts, verts, lookup)
            faces[nfaces,2] = makeIdx(B, &nverts, verts, lookup)
            nfaces += 1

            faces[nfaces,0] = makeIdx(A, &nverts, verts, lookup)
            faces[nfaces,1] = makeIdx(D, &nverts, verts, lookup)
            faces[nfaces,2] = makeIdx(C, &nverts, verts, lookup)
            nfaces += 1
        else:
            faces[nfaces,0] = makeIdx(A, &nverts, verts, lookup)
            faces[nfaces,1] = makeIdx(B, &nverts, verts, lookup)
            faces[nfaces,2] = makeIdx(C, &nverts, verts, lookup)
            nfaces += 1

            faces[nfaces,0] = makeIdx(A, &nverts, verts, lookup)
            faces[nfaces,1] = makeIdx(C, &nverts, verts, lookup)
            faces[nfaces,2] = makeIdx(D, &nverts, verts, lookup)
            nfaces += 1

    return np_verts[:nverts], np_faces

def checkSite(unsigned char[:,:,:] sites,
              unsigned char[:] targets,
              long int[:,:,:] mask):
    """
    Generate a lattice set to 1 where the site type is in `targets`, 0 otherwise.

    :param sites:     Site lattice.
    :type  sites:     unsigned char[:,:,:]
    :param targets:   List of site types to flag.
    :type  targets:   unsigned char[:]
    :param mask:      Output lattice (preallocated)
    :type  mask:      long int[:,:,:]
    """

    cdef int ni = sites.shape[0]
    cdef int nj = sites.shape[1]
    cdef int nk = sites.shape[2]
    cdef int nt = targets.shape[0]
    cdef int st

    mask[:,:,:] = 0

    for t in range(nt):
        st = targets[t]
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    if sites[i,j,k]==st:
                        mask[i,j,k]=1

    cdef cnpy.ndarray np_faces, np_verts



def checkParticle(plattice_type[:,:,:,:] lattice,
                  plattice_type[:] targets,
                  long int[:,:,:] mask):
    """
    Generate a lattice set to 1 where if a particle in `targets` is present in at the lattice site, 0 otherwise.

    :param lattice:   Particle lattice.
    :type  lattice:   uint8/uint32[:,:,:,:]
    :param targets:   List of particle types to flag.
    :type  targets:   unsigned char[:]
    :param mask:      Output lattice (preallocated)
    :type  mask:      long int[:,:,:]
    """
    cdef int ni = lattice.shape[0]
    cdef int nj = lattice.shape[1]
    cdef int nk = lattice.shape[2]
    cdef int np = lattice.shape[3]
    cdef int nt = targets.shape[0]

    mask[:,:,:] = 0

    for t in range(nt):
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    for p in range(np):
                        if lattice[i,j,k,p]==targets[t]:
                            mask[i,j,k]+=1

