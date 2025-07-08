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

"""Simulation analysis mixins for SpatialModel"""
import os.path
import numpy as np
import warnings
from h5py import File as H5File
from . import Lattice
from . import Types

class LatticeAnalysisMixin:

    def particleStatistics(self, particleLattice=None, siteLattice=None):
        """Compute properties of the particle lattice

        Dictionary keys

        countBySpeciesRegion
           Particle counts for each species type for each region.

        concBySpeciesRegion
           Concentration for each species type for each region.

        countByRegion
           Total particle counts for each region.

        concByRegion
           Total concentration for each region.
         
        countBySpecies
           Total particle counts for each species.

        concBySpecies
           Total concentration of each species, averaged over 
           simulation volume

        count
           Total particle count

        conc
           Total concentration, averaged over simulation volume

        vol
           Total simulation volume 

        siteCount
           Total number of subvolumes

        regionVol
           Volume of each region

        regionCounts
           Number of subvolumes in each region
          
        Returns: 
            dict:
               Particle lattice statistics
        """
        particleLattice =  self.particleLattice if particleLattice is None else particleLattice
        siteLattice =  self.siteLattice if siteLattice is None else siteLattice
        
        if particleLattice.ndim == 5:
            self._pCount, self._sCount = Lattice.latticeStatsAll(particleLattice, siteLattice)
        elif particleLattice.ndim == 4:
            self._pCount, self._sCount = Lattice.latticeStatsAll_h5fmt(particleLattice, siteLattice.T)
        else:
            raise ValueError("Wrong number of dimensions for particle lattice")

        regionCounts = {rt: self._sCount[rt.idx] for rt in self.regionList}
        regionVol = {rt: self.siteV*self._sCount[rt.idx] for rt in self.regionList}
        concScl = {rt: 1/self.NA/max(self.siteV,self.siteV*self._sCount[rt.idx]) for rt in self.regionList}
        siteCount = sum(regionCounts.values())
        vol = self.siteV*siteCount
        countBySpeciesRegion  = { sp:{ rt: self._pCount[rt.idx,sp.idx] for rt in self.regionList } for sp in self.speciesList}
        countBySpecies = { sp: np.sum(self._pCount[:,sp.idx]) for sp in self.speciesList }
        countByRegion = { rt: np.sum(self._pCount[rt.idx,1:]) for rt in self.regionList }
        count = np.sum(self._pCount[:,1:])
        conc = count/(vol*self.NA)
        concBySpecies = {sp: ct/(vol*self.NA) for sp,ct in countBySpecies.items()}
        concByRegion  = {rt: ct*concScl[rt] for rt,ct in countByRegion.items()}
        concBySpeciesRegion  = { sp:{ rt: countBySpeciesRegion[sp][rt]*concScl[rt] for rt in self.regionList } for sp in self.speciesList}

        return dict(countBySpeciesRegion=countBySpeciesRegion,
                    countByRegion=countByRegion,
                    countBySpecies=countBySpecies,
                    count=count,
                    siteCount=siteCount,
                    vol=vol,
                    regionVol=regionVol,
                    regionCounts=regionCounts,
                    concBySpeciesRegion=concBySpeciesRegion,
                    concByRegion=concByRegion,
                    concBySpecies=concBySpecies,
                    conc=conc)


class TrajAnalysisMixin:
    def _speciesIdParse(self, regex, species, startIndex=0):
        """
        Parse and convert species identifiers to a list of integer indices.

            This function handles different ways of specifying species:
            1. Using a regex pattern to match species names
            2. A single Species object
            3. A single species name as a string
            4. A list of Species objects
            5. A list of species names as strings
            6. A list of integer indices

            The function returns a list of integer indices representing the species,
            adjusted by the startIndex.

            Args:
                regex (str): A regular expression to match species names
                species (Species, str, list): Species specification
                startIndex (int): An offset to add to the returned indices

            Returns:
                list: A list of integer indices representing the specified species

            Raises:
                RuntimeError: If both or neither regex and species are provided,
                            or if the species specification is unrecognized
        """
        
        if (regex is None) == (species is None):
            raise RuntimeError("`regex` and `species` are mutually exclusive")
        elif regex is not None:
            return [x.idx-1 + startIndex for x in self.speciesList.matchRegex(regex)]
        else:
            if isinstance(species, Types.Species):
                return [species.idx - 1 + startIndex]
            elif isinstance(species, str):
                return [self.speciesList[species].idx - 1 + startIndex]
            elif hasattr(species, "__len__"):
                if isinstance(species[0], Types.Species):
                    return [x.idx-1 + startIndex for x in species]
                elif isinstance(species[0], str):
                    return [self.speciesList[s].idx - 1 + startIndex for s in species]
                elif isinstance(species[0], int):
                    return species

        raise RuntimeError("Unrecognized species specification")


    def _regionIdParse(self, regex, region):
        if (regex is None) == (region is None):
            raise RuntimeError("`regex` xor `region` must be defined")
        elif regex is not None:
            return [x.idx for x in self.regionList.matchRegex(regex)]
        else:
            if isinstance(region, Types.Region):
                return [region.idx]
            elif isinstance(region, str):
                return [self.regionList[region].idx]
            elif hasattr(region, "__len__"):
                if isinstance(region[0], Types.Region):
                    return [x.idx for x in region]
                elif isinstance(region[0], str):
                    return [self.speciesList[s].idx for s in region]
                elif isinstance(region[0], int):
                    return region

        raise RuntimeError("Unrecognized region specification")


    @staticmethod
    def _frameParse(ts, frameStart, frameEnd, timeStart, timeEnd):

        if (frameStart is not None and timeStart is not None):
            RuntimeError("Specify either frameStart or timeStart")
        if (frameEnd is not None and timeEnd is not None):
            RuntimeError("Specify either frameStart or timeStart")

        if timeStart is not None:
            for i,t in enumerate(ts):
                if t>timeStart:
                    frameStart = i-1
                    break

        if timeEnd is not None:
            for i,t in reversed(list(enumerate(ts))):
                if t<timeEnd:
                    frameEnd = i
                    break

        if frameStart is None:
            frameStart = 0

        if frameEnd is None:
            frameEnd = -1

        if frameStart < 0:
            frameStart += len(ts)

        if frameEnd < 0:
            frameEnd += len(ts)

        return frameStart, frameEnd+1


    def _replicateParse(self,replicate):
        if replicate == 'all':
            return list(self.replicates)
        if replicate is None:
            replicate = self.replicate
        if isinstance(replicate, int):
            return [ replicate ]
        return [int(r) for r in replicate]

    def getNumberTrajectory(self,regex=None, species=None, replicate=None,frameStart=None,frameEnd=None,timeStart=None,timeEnd=None):
        """Calculate particle number trajectories

        The time course of the particle numbers in the simulation can 
        be queried with this function.  Species can be selected 
        individually with the `species` key. If a list of species are 
        given instead the sum of the particle numbers for those species 
        types will be returned. A regular expression can also be used 
        to select the species to return. Species types or species names 
        in the form of a string are both acceptible to use in the 
        `species` key. The `regex` and `species` options are mutually 
        exclusive.  Instead of returning the entire trajectory, a
        segment can be selected either through the frame numbers 
        (`frameStart`, `frameEnd`) or through the simulation time 
        (`timeStart`, `timeEnd`).

        Args:
            regex (str):
                Regular expression matching species names. Multiple 
                matches will be summed.
            species ([:py:class:`~jLM.Types.Species`]):
                Species
            replicate (int):
                Replicate to return, default is current replicate.
            frameStart (int):
                Starting frame
            frameEnd (int):
                Ending frame
            timeStart (int):
                Starting time
            timeEnd (int):
                Ending Time

        Returns:
            (:py:class:`~numpy.ndarray(shape=(nt,), dtype=float)`, :py:class:`~numpy.ndarray(shape=(nt,), dtype=float)`):
                The evaluation times and particle counts
        """
        spIdList = self._speciesIdParse(regex, species, startIndex=0)
        replicate = self._replicateParse(replicate)
        if len(replicate) > 1:
            raise ValueError("only one replicate may be specified")

        replicate = "{:07d}".format(replicate[0])

        ts = self.h5['Simulations'][replicate]['SpeciesCountTimes'][...]

        frameStart, frameEnd = self._frameParse(ts, frameStart, frameEnd, timeStart, timeEnd)

        ts = ts[frameStart:frameEnd]
        traj = np.zeros(ts.shape)

        ys = self.h5['Simulations'][replicate]['SpeciesCounts'][frameStart:frameEnd,...]


        for i in spIdList:
            traj += ys[:,i]

        return ts, traj

    def getNumberTrajectoryFromRegion(self, spRegex=None, species=None, regRegex=None, region=None, replicate=None, 
                                            frameStart=None, frameEnd=None, frameDownScale= None, timeStart=None, timeEnd=None, timeDownScale=None):
        """Calculate particle number trajectories for specific regions

        The time course of the particle numbers in the simulation 
        within specific regions can be queried with this function.  
        Species can be selected individually with the species key. If 
        a list of species are given instead the sum of the particle 
        numbers for those species types will be returned. A regular 
        expression can also be used to select the species to return. 
        Species types or species names in the form of a string are 
        both acceptible to use in the `species` key. The `spRegex` and
        `species` options are mutually exclusive.  Instead of 
        returning the entire trajectory, a segment can be selected 
        either through the frame numbers (`frameStart`, `frameEnd`) or 
        through the simulation time (`timeStart`, `timeEnd`).If both given, the frame numbers will be used. 
        The regions to compute particle numbers over are selected similar 
        to the species through the options `regRegex` and `region`.

        Args:
            spRegex (str):
                Regular expression matching species names. Multiple 
                matches will be summed.
            species ([:py:class:`~jLM.Types.Species`]):
                Species
            regRegex (str):
                Regular expression matching region names. Multiple 
                matches will be summed.
            region ([:py:class:`~jLM.Types.Region`]):
                Region
            replicate (int):
                Replicate to return, default is current replicate.
            frameStart (int):
                Starting frame
            frameEnd (int):
                Ending frame
            timeStart (int):
                Starting time
            timeEnd (int):
                Ending Time

        Returns:
            (:py:class:`~numpy.ndarray(shape=(nt,), dtype=float)`, :py:class:`~numpy.ndarray(shape=(nt,), dtype=float)`):
                The evaluation times and particle counts
          
        Note:
            Since this calculation requires reading the entire 
            particle lattice trajectory, it can be slow.
        """
        spIdList = np.array(self._speciesIdParse(spRegex, species, startIndex=1), dtype=np.int32)
        regIdList = np.array(self._regionIdParse(regRegex, region), dtype=np.int32)
        replicate = self._replicateParse(replicate)
        if len(replicate) > 1:
            raise ValueError("only one replicate may be specified")
        
        replicate = "{:07d}".format(replicate[0])
        # if input is time, convert to frame
        if frameStart is None and timeStart is not None and timeEnd is not None:
            dt = float(self.h5['Parameters'].attrs['timestep'])
            frameStart = int(timeStart / dt)
            frameEnd = int(timeEnd / dt)
            # here we support downscale for time
            if timeDownScale is not None:
                frameDownScale = int(timeDownScale / dt)
        elif frameStart is None and timeStart is None:
            raise ValueError("Either frameStart, frameEnd or timeStart, timeEnd must be specified")    
        
        # Here we need to choose between frame or time steps
        if frameStart is not None and frameEnd is not None:
            # Check if the frameStart and frameEnd are within the bounds of the lattice times
            lattice_times_length = len(self.h5['Simulations'][replicate]['LatticeTimes'])
            if frameStart < 0:
                raise ValueError(f"frameStart ({frameStart}) is less than 0, which exceeds the lower boundary of the lattice times.")
            if frameEnd > lattice_times_length:
                raise ValueError(f"frameEnd ({frameEnd}) exceeds the upper boundary ({lattice_times_length}) of the lattice times.")
            
            # not exceed the bounds of the lattice times
            ts = self.h5['Simulations'][replicate]['LatticeTimes'][frameStart:frameEnd]
            if frameDownScale is not None and frameDownScale > 1 and frameDownScale < lattice_times_length:
                ts = ts[::frameDownScale]
            
        else:
            raise ValueError("Either you must give frameStart and frameEnd or timeStart and timeEnd")    
    
        # get the region types in the lattice, and species types in the lattice 
        max_regions = int(self.h5['Model/Diffusion'].attrs['numberSiteTypes'])
        max_species = int(self.h5['Model/Diffusion'].attrs['numberSpecies'])
        actual_regions = len(regIdList)
        actual_species = len(spIdList)

        cacheKey = "SiteParticleCount", replicate
        spc = self._cachedResult(cacheKey)
        
        if spc is None:
            sCount, spc = self._postprocess_siteparticlecounts(rep=replicate,
                                                               reg_num=actual_regions,
                                                               spe_num=actual_species,
                                                               max_regions=max_regions,
                                                               max_species=max_species)
            self._cachedResult(cacheKey, spc)
            self._cachedResult(("SiteCount", replicate), sCount)
        
        # cartesian product over sites and particles
        rs = np.tile(regIdList, len(spIdList))
        ss = np.repeat(spIdList, len(regIdList))
        if len(rs) == 1:
            traj = spc[:,rs[0],ss[0]]
            
        else:
            traj = np.sum(spc[:,np.tile(regIdList, len(spIdList)), np.repeat(spIdList, len(regIdList))], axis=-1)

        return ts, traj


    def getLatticeTrajectory(self, regex=None, species=None, integrated="", replicate=None, 
                             frameStart=None, frameEnd=None, timeStart=None, timeEnd=None):
        """Calculate the spatial distribution of particles versus time

        Species can be selected individually with the species key. If 
        a list of species are given instead the sum of the particle 
        numbers for those species types will be returned. A regular 
        expression can also be used to select the species to return.
        Species types or species names in the form of a string are 
        both acceptible to use in the `species` key. The `regex` and 
        `species` options are mutually exclusive.  Instead of 
        returning the entire trajectory, a segment can be selected 
        either through the frame numbers (`frameStart`, `frameEnd`) or
        through the simulation time (`timeStart`, `timeEnd`). If the 
        full 3-D lattice is not needed, any of the directions can be 
        integrated out.

        Args:
            regex (str):
                Regular expression matching species names. Multiple
                matches will be summed.
            species ([:py:class:`~jLM.Types.Species`]):
                Species
            integrated (str):
                Combination of 'x', 'y', 'z' specifying the directions 
                to integrate out
            replicate (int):
                Replicate to return, default is current replicate.
            frameStart (int):
                Starting frame
            frameEnd (int):
                Ending frame
            timeStart (int):
                Starting time
            timeEnd (int):
                Ending Time

        Returns:
            (:py:class:`~numpy.ndarray(shape=(nt,), dtype=float)`, :py:class:`~numpy.ndarray(dtype=float)`):
                The evaluation times and particle counts
          
        Note:
            Since this calculation requires reading the entire particle
            lattice trajectory, it can be slow.
        """
        spIdList = self._speciesIdParse(regex, species, startIndex=1)
        replicateList = self._replicateParse(replicate)
        if len(replicateList) > 1:
            raise ValueError("only one replicate may be specified")
        replicate = "{:07d}".format(replicateList[0])

        ts = self.h5['Simulations'][replicate]['LatticeTimes'][...]

        frameStart, frameEnd = self._frameParse(ts, frameStart, frameEnd, timeStart, timeEnd)

        def hdfLattice(x): 
            return self.h5['Simulations'][replicate]['Lattice']['{:010d}'.format(x)]

        hdfData = hdfLattice(0)
        latticeDims = hdfData.shape
            
        nFrames = frameEnd-frameStart

        shape = [nFrames]
        integratedDims = []
        for i,d in enumerate('xyz'):
            if d not in integrated:
                shape.append(latticeDims[i])
            else:
                integratedDims.append(i)

        shape = tuple(shape)
        
        integratedDims.append(3) # integrate over the occupancy
        integratedDims.reverse() # integrate higher number axes first to ensure that lower axes do not shift

        cacheKey = "getLatticeTrajectory", replicate, spIdList, [frameStart, frameEnd], integratedDims

        traj = self._cachedResult(cacheKey)
        if traj is None:
            traj=[]

            def latticeGen(f0,f1):
                hdfData = hdfLattice(0)
                buf = np.zeros(hdfData.shape, dtype=hdfData.dtype)
                for i in range(f0,f1):
                    hdfLattice(i).read_direct(buf)
                    yield buf

            for lattice in latticeGen(frameStart, frameEnd):
                ind = np.zeros(latticeDims)
                for spId in spIdList:
                    ind += lattice==spId
                for axis in integratedDims:
                    ind = np.sum(ind,axis=axis)
                traj.append(ind.T)

            traj = np.array(traj)
            self._cachedResult(cacheKey, traj)

        return ts[frameStart:frameEnd], traj

    def getLatticeHistogram(self, regex=None, species=None, integrated="xyz", replicate=None, 
                            frameStart=None, frameEnd=None, timeStart=None,timeEnd=None):
        """Calculate the spatial distribution of particles over an interval

        Species can be selected individually with the species key. If 
        a list of species are given instead the sum of the particle 
        numbers for those species types will be returned. A regular 
        expression can also be used to select the species to return. 
        Species types or species names in the form of a string are 
        both acceptible to use in the `species` key. The `regex` and 
        `species` options are mutually exclusive.  The interval of
        time in which the histogram is computed can be selected either
        through the frame numbers (`frameStart`, `frameEnd`) or 
        through the simulation time (`timeStart`, `timeEnd`). If the 
        full 3-D lattice is not needed, any of the directions can be 
        integrated out.  

        Args:
            regex (str):
                Regular expression matching species names. Multiple 
                matches will be summed.
            species ([:py:class:`~jLM.Types.Species`]):
                Species
            integrated (str):
                Combination of 'x', 'y', 'z' specifying the directions
                to integrate out
            replicate (int):
                Replicate to return, default is current replicate.
            frameStart (int):
                Starting frame
            frameEnd (int):
                Ending frame
            timeStart (int):
                Starting time
            timeEnd (int):
                Ending Time

        Returns:
            :py:class:`~numpy.ndarray(dtype=float)`:
                The average particle counts
          
        Note:
            Since this calculation requires reading the entire particle 
            lattice trajectory, it can be slow.
        """
        spIdList = np.array(self._speciesIdParse(regex, species, startIndex=1), dtype=np.int32)
        assert len(spIdList) > 0
        replicateList = self._replicateParse(replicate)

        def hdfLattice(r,f): 
            return self.h5['Simulations/{:07d}/Lattice/{:010d}'.format(r,f)]

        def hdfLatticeTimes(r): 
            return self.h5['Simulations/{:07d}/LatticeTimes'.format(r)]


        integratedDims = ["xyz".index(a) for a in integrated]

        hist = None
        for r in replicateList:
            hdfData = hdfLatticeTimes(r)
            ts = np.zeros(hdfData.shape, dtype=hdfData.dtype)
            hdfData.read_direct(ts)
            frameStart, frameEnd = self._frameParse(ts, frameStart, frameEnd, timeStart, timeEnd)

            cacheKey = "getLatticeHistogram", "{:07d}".format(r), spIdList, [frameStart, frameEnd], integratedDims
            cachedHist = self._cachedResult(cacheKey)
            if cachedHist is None:
                for f in range(frameStart, frameEnd):
                    lattice = hdfLattice(r,f)
                    if hist is None:
                        latticeTmp = np.zeros(lattice.shape, dtype=lattice.dtype)
                        lattice.read_direct(latticeTmp)
                        hist = Lattice.latticeHistogram(latticeTmp, spIdList, axes=integratedDims, target=None)
                    else:
                        lattice.read_direct(latticeTmp)
                        Lattice.latticeHistogram(latticeTmp, spIdList, axes=integratedDims, target=hist)
                self._cachedResult(cacheKey, hist)
            else:
                if hist is None:
                    hist = cachedHist
                else:
                    hist += cachedHist


        return hist.T/(len(replicateList)*(frameEnd-frameStart))

    def _postprocess_siteparticlecounts(self, rep, reg_num = 16384, spe_num = 16384, max_regions=16384, max_species=16384):
        '''
        This function is used to generate the site-particle counts for a given replicate.
        It is used to cache the results of the site-particle counts so that they can be used
        to generate the site-particle counts for a given replicate.
        return: 
            sCount: The site counts
            pCounts: The particle counts
        '''
        traj = self.h5
        siteLattice = traj['/Model/Diffusion/LatticeSites'][...]

        npl = sum(len(traj['/Simulations/'+rep+'/Lattice']) for rep in traj['/Simulations'])
        warnings.warn("Generating Site/Particle counts for {} frames".format(npl))

        particleLattice = None
        frames = sorted(traj['Simulations/'+rep+"/Lattice"])
        nframes = len(frames)

        pCounts = np.zeros((nframes, reg_num, spe_num), dtype=np.int32)

        for i,frame in enumerate(frames):
            if particleLattice is not None:
                traj['Simulations/'+rep+"/Lattice/"+frame].read_direct(particleLattice)
            else:
                particleLattice = traj['Simulations/'+rep+"/Lattice/"+frame][...]

            pCount, sCount = Lattice.latticeStatsAll_h5fmt(particleLattice, siteLattice)
            
            pCounts[i,:,:] = pCount[:reg_num,:spe_num]

        return sCount, pCounts

    def _cachedResult(self, key, val=None):
        call, rep, *rest = key
        base, _ = os.path.splitext(self.filename)
        cacheFile = base+"-jlmcache.h5"
        argKey = ":".join( ",".join(str(x) for x in xs) for xs in rest)
        if argKey:
            path = "/{}/{}/{}".format(rep, call, argKey)
        else:
            path = "/{}/{}".format(rep, call)

        with H5File(cacheFile, "a", libver='latest') as h5:
            if val is not None:
                if path in h5:
                    h5[path][...] = val
                else:
                    h5.create_dataset(path, data=val, compression='gzip', compression_opts=9)
                if not hasattr(self, "_postprocCache"):
                    self._postprocCache = dict()
                if rep not in self._postprocCache:
                    self._postprocCache[rep] = dict()
                if call not in self._postprocCache[rep]:
                    self._postprocCache[rep][call] = dict()
                self._postprocCache[rep][call][argKey] = val
            else:
                if path in h5:
                    return h5[path][...]