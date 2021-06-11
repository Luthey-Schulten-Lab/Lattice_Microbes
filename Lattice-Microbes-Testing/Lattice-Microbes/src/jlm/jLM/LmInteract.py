# 
# University of Illinois Open Source License
# Copyright 2016-2018 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Tyler M. Earnest
# 

"""SpatialModel mixins for interaction with LM"""
import h5py
import numpy as np
import itertools, os
from . import Lattice
from contextlib import contextmanager

class LmWriteMixin:
    @contextmanager
    def h5(self):
        with h5py.File(self.filename, "r+") as h5:
            yield h5


    def finalize(self):
        import lm
        """Write LM data file"""
        try:
            os.remove(self.filename)
        except FileNotFoundError:
            pass

        lm.SimulationFile.create(self.filename)
        lmFile = lm.SimulationFile(self.filename)

        self.Nrxn=len(self.reactionList)
        self.Nsps=len(self.speciesList)
        self.Nsts=len(self.regionList)

        self._placeAllParticles()

        rm = self._reactionModel()
        dm = self._diffusionModel()
        sm = self._spatialModel()

        lmFile.setReactionModel(rm)
        lmFile.setDiffusionModel(dm)
        lmFile.setSpatialModel(sm)
        lmFile.setDiffusionModelLattice(dm, self.lattice)
        lmFile.setParameter("writeInterval", "{:.10e}".format(self.speciesWriteInterval))
        lmFile.setParameter("latticeWriteInterval", "{:.10e}".format(self.latticeWriteInterval))
        lmFile.setParameter("perfPrintInterval", str(self.perfPrintInterval))
        if self.hookInterval is not None:
            lmFile.setParameter("hookInterval", "{:.10e}".format(self.hookInterval))
        else:
            lmFile.setParameter("hookInterval", "")
        lmFile.setParameter("maxTime", "{:.10e}".format(self.simulationTime))
        lmFile.setParameter("timestep", "{:.10e}".format(self.timestep))
        lmFile.setParameter("speciesNames", ",".join(s.name for s in sorted(self.speciesList,key=lambda x:x.idx)))
        lmFile.setParameter("siteTypeNames", ",".join(s.name for s in sorted(self.regionList,key=lambda x:x.idx)))
        lmFile.close()

        self._storeParamData()

    def copyParticleLattice(self, fname, diffusion_model=False, replicate=1, frame=-1):
        """Copy particle lattice from existing simulation

        Args:
            fname (str):
                LM simulation H5 file

        Keyword Args:
            diffusion_model (bool):
                Use the initial condition particle lattice
            replicate (int):
                LM simulation replicate
            frame (int):
                LM simulation particle lattice frame. Negative values
                index from end.
        """
        with h5py.File(fname, "r") as h5:
            if diffusion_model:
                plattice = h5['/Model/Diffusion/Lattice'][...]
            else:
                ts = h5['/Simulations/{:07d}/LatticeTimes'.format(replicate)][...]
                frame = (ts.size + frame)%ts.size
                plattice = h5['/Simulations/{:07d}/Lattice/{:010d}'.format(replicate, frame)][...]
            inSps = h5['Parameters'].attrs['speciesNames'].decode().split(',')

        spsTranslationMap = {inSps.index(sp.name)+1: sp.idx for sp in self.speciesList}
        txMap = np.zeros(16384, dtype=np.uint32)
        for iOld,iNew in spsTranslationMap.items():
            txMap[iOld] = iNew

        Lattice.hdf2lmRep(self.particleLattice, plattice, txMap)
        self._particlesPlaced = True

        spCts = np.zeros(len(self.speciesList), dtype=int)

        pCount, sCount = Lattice.latticeStatsAll(self.particleLattice, self.siteLattice)

        pCount = np.sum(pCount, axis=0)
        spCts[...] = pCount[1:len(self.speciesList)+1]

        self._spCts = spCts
        self._particlePlacement = []
        self._particleDistCount = []
        self._particleDistConc = []


    def _placeAllParticles(self):
        if self._particlesPlaced:
            return

        self.particleLattice[...] = 0
        particleDensity = np.zeros((len(self.regionList), len(self.speciesList)+1), dtype=np.double)
        siteCounts = Lattice.countSites(self.siteLattice)

        spCts = np.zeros(len(self.speciesList), dtype=int)

        for regid, spid, count in self._particleDistCount:
            spCts[spid-1] += count
            v = 1.0*count/siteCounts[regid]/self.pps
            particleDensity[regid,spid] += v

        for regid, spid, conc in self._particleDistConc:
            count = int(self.siteNAV*conc*siteCounts[regid])
            spCts[spid-1] += count
            v = 1.0*count/siteCounts[regid]/self.pps
            particleDensity[regid,spid] += v

        # particleDensity is the probability that a particle will be placed
        # we need to compute the probability to place empty particles
        particleDensity[:,0] = 1 - np.sum(particleDensity[:,1:], axis=1)
        if np.any(particleDensity < 0):
            raise RuntimeError("Lattice capacity exceeded")

        Lattice.populateLattice(self.particleLattice, self.siteLattice, particleDensity)

        for x,y,z,sp in self._particlePlacement:
            spCts[sp-1] += 1
            if Lattice.appendParticle(self.particleLattice, x, y, z, sp) == 0:
                raise RuntimeError("Failed to add {} at {}: site full.".format(self.speciesList[sp].name, (x,y,z)))
        self._spCts = spCts


    def _buildReactionMatricies(self):
        Nrxn, Nsps = self.Nrxn, self.Nsps

        # Build reaction matricies
        rxtypes = np.zeros(Nrxn, dtype=np.int)
        rxconst = np.zeros(Nrxn, dtype=np.double)
        stoich = np.zeros((Nrxn,Nsps), dtype=np.int)
        depmat = np.zeros((Nrxn,Nsps), dtype=np.int)
        rxnParamMatrix=np.zeros(Nrxn, dtype=np.int)

        for rx in self.reactionList:
            rxconst[rx.idx]=rx.rate._toLM()
            rxnParamMatrix[rx.idx]=rx.rate.idx
            if rx.order == 2: # Second Order Reaction
                if rx.reactants[0] != rx.reactants[1]:
                    rxtypes[rx.idx]=2
                    for sr in rx.reactants:
                        depmat[rx.idx,sr.idx-1]=1
                        stoich[rx.idx,sr.idx-1] -= 1
                else: # Second Order Self Reaction
                    rxtypes[rx.idx]=3
                    depmat[rx.idx,rx.reactants[0].idx-1]=1
                    stoich[rx.idx,rx.reactants[0].idx-1] -= 2
            elif rx.order == 0: # zeroth order
                rxtypes[rx.idx] = 0
            elif rx.order == 1: # First order
                    rxtypes[rx.idx]=1
                    depmat[rx.idx, rx.reactants[0].idx-1]=1
                    stoich[rx.idx, rx.reactants[0].idx-1] -= 1
            else:
                raise NotImplementedError("{} order reations".format(rx.order))

            for sp in rx.products:
                stoich[rx.idx, sp.idx-1] += 1

        return rxtypes, rxconst, stoich, depmat, rxnParamMatrix


    def _reactionModel(self):
        import lm
        rm=lm.ReactionModel()

        Nrxn, Nsps = self.Nrxn, self.Nsps

        rm.set_number_reactions(Nrxn)
        rm.set_number_species(Nsps)

        for s in self._spCts:
            rm.add_initial_species_count(int(s))

        rxtypes, rxconst, stoich, depmat, rxnParamMatrix = self._buildReactionMatricies()

        # Populate reaction object data
        for r in range(Nrxn):
            rm.add_reaction()
            rm.mutable_reaction(r).set_type(int(rxtypes[r]))
            rm.mutable_reaction(r).add_rate_constant(rxconst[r])

        for s in range(Nsps):
            for r in range(Nrxn):
                rm.add_stoichiometric_matrix(int(stoich[r,s]))

        for s in range(Nsps):
            for r in range(Nrxn):
                rm.add_dependency_matrix(int(depmat[r,s]))

        self._rxnParamMatrix=rxnParamMatrix

        return rm


    def _buildReactionLocationMatrix(self):
        RL=np.zeros((self.Nsts, self.Nrxn), dtype=np.int)
        for rxid, regid in self._reactionLocations:
            if regid is None:
                RL[:,rxid]=1
            else:
                RL[regid,rxid]=1
        return RL



    def _buildDiffusionMatrix(self):
        D=np.zeros((self.Nsts,self.Nsts,self.Nsps), dtype=np.double)
        diffParamMatrix=np.zeros((self.Nsts,self.Nsts,self.Nsps), dtype=np.int)
        D[:,:,:] = -1 # flag not visited

        for sp, rFrom, rTo, rateIdx in self._transitionRates:
            if sp is None:
                sp = slice(None)
            else:
                sp = sp - 1
            if rFrom is None:
                rFrom = slice(None)
            if rTo is None:
                rTo = slice(None)

            D[rFrom, rTo,sp] = self.diffRateList[rateIdx].value
            diffParamMatrix[rFrom, rTo, sp] = rateIdx

        if np.any(D<0):
            raise RuntimeError("Diffusion table not complete, ensure that all "
                               "combinations of (species, fromRegion, toRegion) are defined")

        self._diffParamMatrix = diffParamMatrix

        return D


    def _diffusionModel(self):
        import lm
        dm=lm.DiffusionModel()

        Nsps, Nsts, Nrxn = self.Nsps, self.Nsts, self.Nrxn

        dm.set_number_reactions(Nrxn)
        dm.set_number_species(Nsps)
        dm.set_number_site_types(Nsts)
        dm.set_lattice_x_size(self.shape[0])
        dm.set_lattice_y_size(self.shape[1])
        dm.set_lattice_z_size(self.shape[2])
        dm.set_lattice_spacing(self.latticeSpacing)
        dm.set_particles_per_site(self.pps)
        dm.set_bytes_per_particle(self.bytesPerParticle)

        D = self._buildDiffusionMatrix()
        RL = self._buildReactionLocationMatrix()

        for via in range(Nsts):
            for to in range(Nsts):
                for p in range(Nsps):
                    dm.add_diffusion_matrix(D[via,to,p])

        for r in range(Nrxn):
            for s in range(Nsts):
                dm.add_reaction_location_matrix(int(RL[s,r]))

        return dm

    def _spatialModel(self):
        import lm
        return lm.SpatialModel()

    def _storeParamData(self):
        # need to store mapping between parameter objects and ids

        with h5py.File(self.filename) as h5:
            diffConstNames = ",".join(x.name for x in self.diffRateList)
            h5['Parameters'].attrs['diffConstNames'] = np.string_(diffConstNames)
            h5.create_dataset("Parameters/DiffusionParameterMatrix", data=self._diffParamMatrix)
            diffParams = np.zeros(len(self.diffRateList), dtype=np.double)
            for k in self.diffRateList:
                diffParams[k.idx] = k.value
            h5.create_dataset("Parameters/DiffusionParameterValues", data=diffParams)

            rxnConstNames = ",".join(x.name for x in self.rxnRateList)
            h5['Parameters'].attrs['reactionConstNames'] = np.string_(rxnConstNames)
            h5.create_dataset("Parameters/ReactionParameterMatrix", data=self._rxnParamMatrix)
            rateParams = np.zeros(len(self.rxnRateList), dtype=np.double)
            for k in self.rxnRateList:
                rateParams[k.idx] = k.value
            h5.create_dataset("Parameters/ReactionParameterValues", data=rateParams)
            rxnConstOrder = np.array([r.order for r in self.rxnRateList],dtype=np.int)
            h5.create_dataset("Parameters/ReactionParameterOrder", data=rxnConstOrder)
            h5['Model'].attrs['name'] = np.string_(self.name)


class LmReadMixin:

    def _loadFromH5(self, replicate):
        self._initParameters()
        self._initSpecies()
        self._initRegions()
        self._initDiffParams()
        self._initRateParams()
        self._readDiffusionTable()
        if self.Nrxn > 0:
            self._readReactions()
        self._readReplicates()
        self._initLattice()
        if len(self.replicates) > 0:
            self.setReplicate(replicate)

    def _readReplicates(self):
        self.replicates = [int(x) for x in self.h5['Simulations']]

    def setReplicate(self, replicate=1):
        """Select replicate from HDF5 file

        Args:
            replicate (int):
                Replicate index (1-based)
        """
        self.replicate = replicate
        replicate = "{:07d}".format(replicate)
        self.trajData = self.h5['Simulations'][replicate]
        self._loadReplicate()
        return self

    def _initLattice(self):
        self.siteLattice[...] = self.h5['/Model/Diffusion/LatticeSites'][...].transpose(2,1,0)
        inputPl= np.array(self.h5['/Model/Diffusion/Lattice'], dtype=np.uint32)
        Lattice.hdf2lmRep(self.particleLattice, inputPl)

    def _loadReplicate(self):
         self.speciesCountTimes = self.trajData['SpeciesCountTimes'][...]

         cs = self.trajData['SpeciesCounts'][...]

         self.speciesCounts = dict()
         for sp in self.speciesList:
             self.speciesCounts[sp] = cs[:,sp.idx-1]

         self.latticeTimes = self.trajData['LatticeTimes'][...]

         last = self.latticeTimes.size-1

         inputPl= np.array(self.trajData['Lattice']['{:010d}'.format(last)], dtype=np.uint32)
         Lattice.hdf2lmRep(self.particleLattice, inputPl)


    def _initSpecies(self):
        speciesNames = self.h5['Parameters'].attrs['speciesNames'].decode().split(',')
        for s in speciesNames:
            self.speciesList.get(s)
        self.speciesList.freeze()
        assert self.Nsps == len(speciesNames)

    def _initRegions(self):
        try:
            regionNames = self.h5['Parameters'].attrs['siteTypeNames'].decode().split(',')
        except KeyError:
            regionNames = ["site{:03d}".format(f) for f in range(int(self.h5['Model/Diffusion'].attrs['numberSiteTypes']))]
        for r in regionNames:
            self.regionList.get(r)
        self.regionList.freeze()
        assert self.Nsts == len(regionNames)

    def _initParameters(self):
        self.timestep = float(self.h5['Parameters'].attrs['timestep'])
        self.latticeWriteInterval = float(self.h5['Parameters'].attrs['latticeWriteInterval'])
        self.speciesWriteInterval = float(self.h5['Parameters'].attrs['writeInterval'])
        p = self.h5['Parameters'].attrs.get('hookInterval', b"")
        if p == b'':
            self.hookInterval = 0.0
        else:
            self.hookInterval = float(p)
        self.simulationTime = float(self.h5['Parameters'].attrs['maxTime'])
        self.bytesPerParticle = int(self.h5['Model/Diffusion'].attrs['bytes_per_particle'])
        self.pps = int(self.h5['Model/Diffusion'].attrs['particlesPerSite'])
        self.Nsps = int(self.h5['Model/Diffusion'].attrs['numberSpecies'])
        self.Nrxn = int(self.h5['Model/Diffusion'].attrs['numberReactions'])
        self.Nsts = int(self.h5['Model/Diffusion'].attrs['numberSiteTypes'])


    def _initDiffParams(self):
        diffRateNames = self.h5['Parameters'].attrs['diffConstNames'].decode().split(',')
        diffParams = self.h5["Parameters/DiffusionParameterValues"][...]
        for i in range(len(diffRateNames)):
            self.diffusionConst(diffRateNames[i], diffParams[i])
        self.NdiffParams = len(self.diffRateList)

    def _readDiffusionTable(self):
        D = self.h5['Model/Diffusion/DiffusionMatrix'][...]
        assert D.shape == (self.Nsts, self.Nsts, self.Nsps)

        diffParamMap = self.h5["Parameters/DiffusionParameterMatrix"][...]

        for rFrom, rTo, sp in itertools.product(self.regionList, self.regionList, self.speciesList):
            v = D[rFrom.idx, rTo.idx, sp.idx-1]
            rate = self.diffRateList[diffParamMap[rFrom.idx,rTo.idx, sp.idx-1]]
            assert v==rate.value
            self.transitionRate(sp, rFrom, rTo, rate)


    def _initRateParams(self):
        if self.Nrxn > 0 :
            rxnOrder = self.h5["Parameters/ReactionParameterOrder"][...]
            rxnRateNames = self.h5['Parameters'].attrs['reactionConstNames'].decode().split(',')
            rxnParams = self.h5['Parameters/ReactionParameterValues'][...]
            for i in range(len(rxnRateNames)):
                self.rateConst(rxnRateNames[i], rxnParams[i], rxnOrder[i])
            self.NrxnParams = len(self.rxnRateList)
        else:
            self.NrxnParams = 0



    def _readReactions(self):
        rxnParamMatrix = self.h5['Parameters/ReactionParameterMatrix'][...]
        S = self.h5['Model/Reaction/StoichiometricMatrix'][...]
        rMat = self.h5['Model/Reaction/DependencyMatrix'][...]
        pMat = S+rMat

        RL = self.h5['Model/Diffusion/ReactionLocationMatrix'][...]
        rates = self.h5['Model/Reaction/ReactionRateConstants'][...]
        Nsps, Nrxn = S.shape
        _, Nsts = RL.shape
        assert Nsps == self.Nsps
        assert Nrxn == self.Nrxn
        assert Nsts == self.Nsts

        for xid in range(Nrxn):
            ps, rs, regions = [], [], []
            for sid in range(Nsps):
                ct = rMat[sid, xid]
                if ct > 0:
                    while ct != 0:
                        rs.append( self.speciesList[sid+1] )
                        ct -= 1
                ct = pMat[sid, xid]
                if ct > 0:
                    while ct != 0:
                        ps.append( self.speciesList[sid+1] )
                        ct -= 1
            for rid in range(Nsts):
                if RL[xid, rid] == 1:
                    regions.append( self.regionList[rid] )

            v = rates[xid, 0]
            rate = self.rxnRateList[rxnParamMatrix[xid]]
            assert rate._toLM() == v
            rxn = self.reaction(rs, ps, rate)

            for r in regions:
                self.assignReaction(rxn, r)

