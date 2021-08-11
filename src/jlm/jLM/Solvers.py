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

"""Generic RDME solvers"""
import numpy as np

from . import RDME
from . import Lattice
from . import Types as T

def makeSolver(rdmeSolver, customSolver):
    """Create custom solver

    Note:
        The user defined class should define a 
        `hookSimulation(self, t, lattice)` method where `t` is the current 
        simulation time and `lattice` is the current state of the 
        :py:class:`lm.CudaByteLattice`. 

    Warning:
        If `__init__` is overridden, `super().__init__` must be called first.

    Args:
        rdmeSolver (:py:class:`lm.RDMESolver`): 
            Base solver class
        customSolver (class):
            User defined solver class

    Returns:
        class:
            Composed class
    """
    name = customSolver.__name__ + "_" + rdmeSolver.__name__
    return type(name, (customSolver, rdmeSolver), {})

class ConstBoundaryConc:
    """Solver with fixed concentration boundary conditions

    The simulation hook runs a cythonized function which destructively sets
    the specified boundaries to a prescribed concentration. """
    def __init__(self, lmFile, exact=False):
        super(ConstBoundaryConc,self).__init__()
        self.exact=exact
        if isinstance(lmFile, (RDME.Sim, RDME.File)):
            self.rdme = lmFile
        else:
            self.rdme = RDME.File(lmFile)

        self.dists = []
        self.distTable = None
        self.bMask = np.zeros_like(self.rdme.siteLattice)
        self.bLattice = np.zeros_like(self.rdme.siteLattice)

    def setBoundary(self, species, concs, boundary):
        r"""Specify boundary conditions

        Specify a boundary condition using a species (or list of species), a 
        concentration (or list of concentrations; in mol/L), and a specification 
        of the boundary in terms of a binary lattice. Subsequent calls will add 
        new boundary conditions.  When the boundary region given by a subsequent 
        call overlaps with a previously described region, the later call will 
        override the earlier conditions.

        Args:
            species (:py:class:`~jLM.Types.Species`):
                Species at boundary
            concs (float):
                Concentrations
            boundary (:py:class:`~numpy.ndarray`):
                Binary lattice describing the boundary
        """
        if isinstance(species, str):
            species = [self.rdme.speciesList[species]]
            concs = [concs]
        elif isinstance(species, T.Species):
            concs = [concs]
            species = [species]
        elif isinstance(species[0], str):
            species = [self.rdme.speciesList[x] for x in species]
        elif isinstance(species[0], T.Species):
            pass
        else:
            raise TypeError()

        dist = np.zeros(256)
        for sp, conc in zip(species, concs):
            particle_number = self.rdme.siteNAV*conc
            dist[sp.idx] = particle_number
        dist[0] = self.rdme.pps - np.sum(dist)
        dist /= self.rdme.pps

        if dist[0] < 0:
            raise RuntimeError("Requested concentration exceeds site capacity")

        self.bMask[boundary>0] = 1
        self.bLattice[boundary>0] = len(self.dists)
        self.dists.append(dist)

    def hookSimulation(self, t, lattice):
        """Replaces all particles on boundary according to B.C.s"""
        if self.distTable is None:
            self.distTable = np.vstack(self.dists)
        particles = lattice.getParticleLatticeView()
        Lattice.populateLattice(particles, self.bLattice,  self.distTable, mask=self.bMask, exact=self.exact)
        return 1
