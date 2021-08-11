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
# Author(s): Tyler M. Earnest, Zane Thornburg
# 

"""Spatial model"""
import h5py

from . import BaseTypes as BT
from . import Types as T
from . import JupyterDisplay as JD
from . import LmInteract as LM
from . import Analysis as AN

def _idxOrNone(x):
    if x is None:
        return x
    else:
        return x.idx

def _copyDoc(method):
    def _doc(func):
        func.__doc__ = method.__doc__
        return func
    return _doc

class SpatialModel(AN.LatticeAnalysisMixin, JD.JupyterDisplayMixin):
    """Base RDME model """
    @property
    def simulationTime(self):
        """Total simulated time"""
        return getattr(self, "_simulationTime", None)

    @simulationTime.setter
    def simulationTime(self, val):
        setattr(self, "_simulationTime", val)

    @property
    def speciesWriteInterval(self):
        """Simulation time between writing the species counts to the HDF5 file"""
        return getattr(self, "_speciesWriteInterval", None)

    @speciesWriteInterval.setter
    def speciesWriteInterval(self, val):
        setattr(self, "_speciesWriteInterval", val)

    @property
    def latticeWriteInterval(self):
        """Simulation time between writing the lattice to the HDF5 file"""
        return getattr(self, "_latticeWriteInterval", None)

    @latticeWriteInterval.setter
    def latticeWriteInterval(self, val):
        setattr(self, "_latticeWriteInterval", val)

    @property
    def perfPrintInterval(self):
        """Real elapsed time between writing simulation progress to stdout"""
        if not hasattr(self, "_perfPrintInterval"):
            self._perfPrintInterval = 60
        return getattr(self, "_perfPrintInterval", None)

    @perfPrintInterval.setter
    def perfPrintInterval(self, val):
        setattr(self, "_perfPrintInterval", val)

    @property
    def hookInterval(self):
        """Simulation time between calls to hookSimulation"""
        return getattr(self, "_hookInterval", None)

    @hookInterval.setter
    def hookInterval(self, val):
        setattr(self, "_hookInterval", val)

    def __init__(self, name, filename, dimensions, latticeSpacing ):
        self.name = name                                   #: str: Simulation name

        self.filename = filename                           #: str: LM data filename

        self.NA = 6.02214085774e23                         #: float: Avogadro's number

        self._particlePlacement = [] # [ (x,y,z,spid), ... ]
        self._particleDistCount = [] # [ (regid,spid, count), ... ]
        self._particleDistConc = [] # [ (regid,spid, conc), ... ]
        self._transitionRates = [] # [ (sp, from, to, rateId), ...] where {sp,from,to} = None|id, None means all
        self._reactionLocations = [] # [ (rxid, regid), ...] 

        self.resizeLattice(dimensions, latticeSpacing)

    def resizeLattice(self, dimensions, latticeSpacing):
        import lm
        """Resize lattice

        Args:
            dimensions (int,int,int):
                New lattice dimensions
            latticeSpacing (float):
                Lattice spacing in meters
        """
        nx,ny,nz = dimensions
        self.siteV = 1000 * latticeSpacing**3                             #: float: Subvolume size in liters
        self.siteNAV = self.siteV * self.NA                               #: float: Conversion factor from concentration to site occupancy
        self.shape = nz,ny,nx                                             #: (int, int, int): Site lattice dimensions

        self.latticeSpacing = latticeSpacing                              #: float: Lattice spacing in meters
        self.pps = lm.getCompiledLatticeMaxOccupancy()                    #: int: Particles per site

        self.lattice = lm.IntLattice(nz,ny,nx, latticeSpacing, self.pps)  #: :py:class:`lm.ByteLattice`: LM lattice object
        self.siteLattice = self.lattice.getSiteLatticeView()              #: :py:class:`~numpy.ndarray`: NumPy view of site lattice
        self.particleLattice = self.lattice.getParticleLatticeView()      #: :py:class:`~numpy.ndarray`: NumPy view of particle lattice
        self.maxConcentration = self.pps/self.siteNAV                     #: float: Concentration of a fully packed lattice

    def placeNumber(self, sp, x, y, z, n):
        """Place a number of particles at a specific subvolume

        Args:
            sp (:py:class:`~jLM.Types.Species`):
                Species type
            x (int):
                x-coordinate
            y (int):
                y-coordinate
            z (int):
                z-coordinate
            n (int):
                Number of particles
        """
        x,y,z,n = map(int, (x,y,z,n))
        for _ in range(n):
            self._particlePlacement.append(( z, y, x, sp.idx))

    def distributeNumber(self, sp, reg, count):
        """Distribute a number of particles uniformly through a region

        Args:
            sp (:py:class:`~jLM.Types.Species`]):
                Species type
            reg (:py:class:`~jLM.Types.Region`]):
                Region type
            count (int):
                Number of particles
        """
        count = int(count)
        cs = self._particleDistCount
        cs = [x for x in cs if not (x[0]==reg.idx and x[1]==sp.idx)]
        cs.append((reg.idx, sp.idx, count))
        self._particleDistCount = cs

    def distributeConcentration(self, sp, reg, conc):
        """Distribute a concentration of particles uniformly through a region

        Args:
            sp (:py:class:`~jLM.Types.Species`]):
                Species type
            reg (:py:class:`~jLM.Types.Region`]):
                Region type
            conc (float):
                Concentration of particles
        """
        cs = self._particleDistConc
        cs = [x for x in cs if not (x[0]==reg.idx and x[1]==sp.idx)]
        cs.append((reg.idx, sp.idx, conc))
        self._particleDistConc = cs
        

    def transitionRate(self, sp, rFrom, rTo, rate, value=None):
        """Define the diffusive transition rate between regions

        This allows for fine grained detail over the diffusion matrix.
        The transition is defined for the species `sp` from region
        `rFrom` to `rTo`. If `None` is given for any of these
        parameters, then the entire axis of the matrix is affected,
        i.e. `None` is a wildcard. The rate is specified with the
        `rate` parameter, either as a DiffusionConst instance, or as
        a string to lookup the previously created DiffusionConst. If a
        string is provided and a DiffusionConst instance has not been
        created, providing the parameter `value` will create the new
        constant.

        Args:
            sp (:py:class:`~jLM.Types.Species`]): 
                Species type or None
            rFrom (:py:class:`~jLM.Types.Region`]):
                Region type or None
            rTo (:py:class:`~jLM.Types.Region`]):
                Region type or None
            rate (:py:class:`~jLM.Types.DiffusionConst`]):
                Diffusion rate
            value (float):
                Value of diffusion constant if new
        """
        if not isinstance(rate, T.DiffusionConst):
            if isinstance(rate, str) and isinstance(value, float):
                rate = self.diffusionConst(rate, value)
            else:
                raise TypeError("rate")
        self._transitionRates.append( (_idxOrNone(sp), _idxOrNone(rFrom), _idxOrNone(rTo), rate.idx) )

    def assignReaction(self, reaction, region):
        """Assign a reaction to a region

        Args:
            reaction (:py:class:`~jLM.Types.Reaction`):
                Reaction
            region (:py:class:`~jLM.Types.Region`):
                Target region
        """
        e = (reaction.idx, region.idx)
        if e not in self._reactionLocations:
            self._reactionLocations.append( (reaction.idx, region.idx) )

    def species(self, name, **kwargs):
        """Lookup/create species instance

        Args:
            name (str):
                Name of the species
            texRepr (str):
                TeX math mode representation of the species
            annotation (str):
                Description text

        Returns:
            :py:class:`~jLM.Types.Species`: 
                Species instance
        """
        return self.speciesList.get(name, **kwargs)

    def region(self, name, **kwargs):
        """Lookup/create region instance

        Args:
            name (str):
                Name of the region
            texRepr (str):
                TeX math mode representation of the region
            annotation (str):
                Description text

        Returns:
            :py:class:`~jLM.Types.Region`:
               Region instance
        """
        return self.regionList.get(name, **kwargs)

    def reaction(self, reactants, products, rate, value=None, regions=None, **kwargs):
        r"""Lookup/create reaction instance

        The product/reactants can be specified as a string, Species 
        instance, list of strings, list of Species instances, or the 
        empty list, denoting no species. The reaction rate can be 
        created in place if a string is given for `rate` and the value
        of the rate constant is specified by `value`.  Rates must be 
        provided in standard chemical kinetic units, e.g. for a reaction

        .. math:: \mathrm{A} + \mathrm{B} \overset{k}{\to} \mathrm{C},

        the rate coefficent, :math:`k` has units of 
        :math:`\mathrm{M}^{-1}\mathrm{s}^{-1}`


        Args:
            reactants ([:py:class:`~jLM.Types.Species`]):
                Species, or list of species acting as reactants
            products ([:py:class:`~jLM.Types.Species`]):
                Species, or list of species acting as products
            rate ([:py:class:`~jLM.Types.RateConst`]):
                Rate constant
            value (float):
                Optional value of rate constant
            regions (py:class:`~jLM.Types.Region`]):
                List of applicable regions

        Returns:
            :py:class:`~jLM.Types.Reaction`:
               Reaction instance
        """
        if not isinstance(rate, T.RateConst):
            if isinstance(rate, str) and isinstance(value, float):
                rate = self.rateConst(rate, value, len(reactants))
            else:
                raise TypeError("rate")

        rxn = self.reactionList.get(reactants, products, rate, **kwargs)
        if regions is not None:
            if not isinstance(regions, list):
                regions = [regions]
            for reg in regions:
                reg = self.region(reg) if isinstance(reg, str) else reg
                self.assignReaction(rxn, reg)
        return rxn


    def rateConst(self, rate, value, order, **kwargs):
        """Lookup/create reaction rate constant instance

        Args:
            rate (str):  
                Name of rate constant
            value (float): 
                Value of rate constant
            order (int): 
                Order of reaction
            texRepr (str): 
                TeX math mode representation of the rate constant
            annotation (str): 
                Description text

        Returns:
            :py:class:`~jLM.Types.RateConst`:
                RateConst instance
        """
        return self.rxnRateList.get(rate, value, order, **kwargs)

    def diffusionConst(self, rate, value, **kwargs):
        """Lookup/create diffusion constant instance

        Args:
            rate (str):
                Name of diffusion constant
            value (float):
                Value of diffusion constant
            texRepr (str):
                TeX math mode representation of the rate constant
            annotation (str):
                Description text

        Returns:
            :py:class:`~jLM.Types.DiffusionConst`:
                DiffusionConst instance
        """
        return self.diffRateList.get(rate, value, **kwargs)

    @property
    def diffusionZero(self):
        """Lookup/create a zero-valued diffusion constant instance

        Returns:
            :py:class:`~jLM.Types.DiffusionConst`:
                DiffusionConst instance
        """
        return self.diffRateList.get('0', 0.0, texRepr=r'D_\varnothing')

    def maxDiffusionRate(self, latticeSpacing=None, dt=None):
        """Compute max allowed diffusion constant for the simulation 

        Args:
            latticeSpacing (float):
                Lattice spacing in meters
            dt (float):
                Timestep in seconds

        Returns:
            float:
                Maximum diffusion constant in m^/s
        """
        if dt is None:
            dt = self.timestep
        if latticeSpacing is None:
            latticeSpacing = self.latticeSpacing
        
        return self.latticeSpacing**2/6/dt

    @property
    def diffusionFast(self):
        """Lookup/create a maximum-valued diffusion constant instance

        Returns:
            :py:class:`~jLM.Types.DiffusionConst`: 
                DiffusionConst instance
        """
        return self.diffRateList.get("fast", self.maxDiffusionRate(), texRepr=r'D_\infty')

    def setMaximumTimestep(self):
        """Set the simulation timestep using the fastest diffusion rate"""
        v = 0
        for r in self.diffRateList:
            v = max(v, r.value)

        if v <= 0:
            raise RuntimeError("Define diffusion constants prior to calling")
        else:
            return (self.latticeSpacing)**2 / (6*v)

    def run(self, solver=None, replicate=1, seed=None, cudaDevices=None, checkpointInterval=0):
        """Run the RDME simulation

        Args:
            solver (:py:class:`lm.RDMESolver`):
                Rdme solver
            replicate (int):
                Index of replicate
            seed (int):
                RNG seed
            cudaDevices ([int]):
                List of CUDA device indexes
            checkpointInterval (int):
                Number of seconds between checkpoints

        Returns:
            :py:class:`jLM.File`:
                Simulation result
        """
        import lm
        solver = solver or lm.MpdRdmeSolver()
        cudaDevices = cudaDevices or [0]
        if seed is not None:
            f = lm.simulationFile(self.filename)
            f.setParameter("seed", str(seed))
            f.close()
        lm.runSolver(self.filename, replicate, solver, cudaDevices, int(checkpointInterval))
        return File(self.filename)


class Sim(LM.LmWriteMixin, SpatialModel):
    """Define and run an RDME simulation"""
    @_copyDoc(AN.LatticeAnalysisMixin.particleStatistics)
    def particleStatistics(self, particleLattice=None, siteLattice=None):
        self._placeAllParticles()
        return super().particleStatistics(particleLattice,siteLattice)

    def __init__(self, name, filename, dimensions, latticeSpacing, regionName, dt=None):
        """Create new RDME object

        Args:
            name (str):
                Name of simulation
            filename (str):
                LM data filename
            dimensions ((nx,ny,nz)):
                Dimensions of lattice
            latticeSpacing (float):
                Lattice spacing in meters
            regionName (str):
                Name of the default region
            dt (float):
                Timestep
        """
        super().__init__(name, filename, dimensions, latticeSpacing)

        self.speciesList = BT.SimObjs(self, T.BuilderSpecies,idbase=1) #: :py:class:`~jLM.BaseTypes.SimObjs`: List of species types
        self.regionList = BT.SimObjs(self, T.BuilderRegion)            #: :py:class:`~jLM.BaseTypes.SimObjs`: List of regions
        self.reactionList = BT.SimObjs(self, T.BuilderReaction)        #: :py:class:`~jLM.BaseTypes.SimObjs`: List of reactions
        self.rxnRateList = BT.SimObjs(self, T.RateConst)               #: :py:class:`~jLM.BaseTypes.SimObjs`: List of reaction rates
        self.diffRateList = BT.SimObjs(self, T.DiffusionConst)         #: :py:class:`~jLM.BaseTypes.SimObjs`: List of diffusion constants
        self.sp = self.speciesList.getAutoNamespace()             #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient species access
        self.reg = self.regionList.getAutoNamespace()             #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient region access
        self.rc = self.rxnRateList.getAutoNamespace()             #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient reaction rate access
        self.dc = self.diffRateList.getAutoNamespace()            #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient diffusion constant access

        self._particlesPlaced = False
        self.filename = filename
        self.bytesPerParticle = 4
        self.region(regionName)
        if dt is not None:
            self.timestep = dt


class File(AN.TrajAnalysisMixin, JD.FileJupyterMixin, LM.LmReadMixin, SpatialModel):
    """Load a previously defined simulation"""
    def __init__(self, fname, replicate=1):
        """Load a RDME simulation file

        Args:
            fname (str):
                LM data filename
            replicate (int):
                Replicate to load initially

        Note:
            jLM includes some metadata in the Lattice Microbes HDF5 file which 
            is necessary to reload the model.
        """
        self.h5 = h5py.File(fname, "r")

        name = self.h5['Model'].attrs['name'].decode("ascii")
        naturalDims = (int(self.h5['Model/Diffusion'].attrs['latticeZSize']),
                       int(self.h5['Model/Diffusion'].attrs['latticeYSize']),
                       int(self.h5['Model/Diffusion'].attrs['latticeXSize']))
        latticeSpacing = float(self.h5['Model/Diffusion'].attrs['latticeSpacing'])

        super().__init__(name, fname, naturalDims, latticeSpacing)
        self.speciesList = BT.SimObjs(self, T.TrajSpecies,idbase=1) 
        self.regionList = BT.SimObjs(self, T.TrajRegion)
        self.sp = self.speciesList.getAutoNamespace()
        self.reg = self.regionList.getAutoNamespace()

        self.speciesList = BT.SimObjs(self, T.TrajSpecies,idbase=1) #: :py:class:`~jLM.BaseTypes.SimObjs`: List of species types
        self.regionList = BT.SimObjs(self, T.TrajRegion)            #: :py:class:`~jLM.BaseTypes.SimObjs`: List of regions
        self.reactionList = BT.SimObjs(self, T.Reaction)            #: :py:class:`~jLM.BaseTypes.SimObjs`: List of reactions
        self.rxnRateList = BT.SimObjs(self, T.RateConst)            #: :py:class:`~jLM.BaseTypes.SimObjs`: List of reaction rates
        self.diffRateList = BT.SimObjs(self, T.DiffusionConst)      #: :py:class:`~jLM.BaseTypes.SimObjs`: List of diffusion constants
        self.sp = self.speciesList.getAutoNamespace()          #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient species access
        self.reg = self.regionList.getAutoNamespace()          #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient region access
        self.rc = self.rxnRateList.getAutoNamespace()          #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient reaction rate access
        self.dc = self.diffRateList.getAutoNamespace()         #: :py:class:`~jLM.BaseTypes.Namespace`: Convienient diffusion constant access

        self.counts= dict()
        self._loadFromH5(replicate)
