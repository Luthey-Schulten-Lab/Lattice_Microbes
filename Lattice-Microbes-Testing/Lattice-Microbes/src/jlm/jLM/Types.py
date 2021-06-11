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

"""Simulation objects. Species, Reaction, etc."""
from . import BaseTypes as BT
from . import ColorGen as CG
import numpy as np

def _rmTex(x):
    return r"\mathrm{" + x.replace("_", r"\_") + "}"

class Species(BT.SimObj):
    """Chemical species"""
    @classmethod
    def _unique_id(cls, name, **kwargs):
        return name

    def _TeX(self):
        return _rmTex(self.name)

    def _html(self):
        return "<span class='jLMnum'>"+ self.name + "</span>"

    def show(self):
        return self._sim.showSpecies(self)

class BuilderSpecies(Species):
    """Chemical species, in the context of simulation construction"""
    def placeParticle(self, x, y, z, count):
        """Place in a subvolume.
        
        Args:
            x (int):
                x-coordinate
            y (int):
                y-coordinate
            z (int):
                z-coordinate
            count (int):
                Number of particles
        """
        self._sim.placeNumber(self, x, y, z, count)
        return self

    def placeNumberInto(self, region, count):
        """Distribute a number within a region
        
        Args:
            region (:py:class:`~jLM.Types.Region`):
                Target region
            count (int):
                number
        """
        self._sim.distributeNumber(self, region, count)
        return self

    def placeConcentrationInto(self, region, conc):
        """Distribute a concentration within a region
        
        Args:
            region (:py:class:`~jLM.Types.Region`):
                Target region
            conc (float):
                Concentration in mol/L
        """
        self._sim.distributeConcentration(self, region, conc)
        return self

    def diffusionRate(self, region, rate, value=None):
        """Set the diffusion rate

        Args:
            region (:py:class:`~jLM.Types.Region`):
                Region
            rate (:py:class:`~jLM.Types.DiffusionConst`):
                Diffusion rate
            value (float):
                If str given for rate and the DiffusionConst does not exist, 
                create it with this value
        """
        self._sim.transitionRate(self, region, region, rate, value=value)
        return self

    def transitionRate(self, regionFrom, regionTo, rate, value=None):
        """Set the diffusion rate

        Args:
            regionFrom (:py:class:`~jLM.Types.Region`):
                Source region
            regionTo (:py:class:`~jLM.Types.Region`):
                Destination region
            rate (:py:class:`~jLM.Types.DiffusionConst`):
                Diffusion rate
            value (float):
                If str given for rate and the DiffusionConst does not exist,
                create it with this value
        """
        self._sim.transitionRate(self, regionFrom, regionTo, rate, value=value)
        return self


class TrajSpecies(Species):
    """Chemical species, in the context of trajectory analysis"""
    def getNumberTrajectory(self, **kwargs):
        """Get particle number trajectory for this species

        See :py:meth:`~jLM.RDME.File.getNumberTrajectory`
        """
        return self._sim.getNumberTrajectory(species=self, **kwargs)
    def getLatticeTrajectory(self, **kwargs):
        """Get lattice trajectory for this species

        See :py:meth:`~jLM.RDME.File.getLatticeTrajectory`
        """
        return self._sim.getLatticeTrajectory(species=self, **kwargs)
    def getLatticeHistogram(self, **kwargs):
        """Get lattice histogram for this species

        See :py:meth:`~jLM.RDME.File.getLatticeHistogram`
        """
        return self._sim.getLatticeHistogram(species=self, **kwargs)
    def getNumberTrajectoryFromRegion(self, **kwargs):
        """Get number trajectory for a region spec for this species

        See :py:meth:`~jLM.RDME.File.getNumberTrajectoryFromRegion`
        """
        return self._sim.getNumberTrajectoryFromRegion(species=self, **kwargs)


class Region(BT.SimObj):
    """Simulation compartment"""
    @classmethod
    def _unique_id(cls, name, **kwargs):
        return name

    @property
    def volume(self):
        return np.sum(self._sim.siteLattice==self.idx)*self._sim.latticeSpacing**3*1000

    @property
    def NAV(self):
        return self.volume*self._sim.NA

    def _floatFgColor(self):
        return CG.ColorSeq.floatFg(self.idx)

    def _hexFgColor(self):
        return CG.ColorSeq.strFg(self.idx)

    def _floatColor(self):
        return CG.ColorSeq.float(self.idx)

    def _hexColor(self):
        return CG.ColorSeq.str(self.idx)

    def _html(self):
        return r'<span class="jLMregion" style="color:{};background:{};">{}</span>'.format(self._hexFgColor(), self._hexColor(), self.name)

    def _TeX(self):
       return r"\texttt{" + self.name.replace("_", r"\textunderscore") + "}"

class BuilderRegion(Region):
    """Compartment, in the context of simulation construction"""
    def diffusionRate(self, rate, value=None):
        """Set diffusion rate for all particles

        Args:
            rate (:py:class:`~jLM.Types.DiffusionConst`):
                Diffusion rate
            value (float):
                If str given for rate and the DiffusionConst does not exist, 
                create it with this value
        """
        self._sim.transitionRate(None, self, self, rate, value=value)
        return self

    def transitionRateOut(self, rate, value=None):
        """Set diffusion rate for all particles leaving this region

        Args:
            rate (:py:class:`~jLM.Types.DiffusionConst`):
                Diffusion rate
            value (float):
                If str given for rate and the DiffusionConst does not exist,
                create it with this value
        """
        self._sim.transitionRate(None, self, None, rate, value=value)
        return self

    def transitionRateIn(self, rate, value=None):
        """Set diffusion rate for all particles entering this region

        Args:
            rate (:py:class:`~jLM.Types.DiffusionConst`):
                Diffusion rate
            value (float):
                If str given for rate and the DiffusionConst does not exist,
                create it with this value
        """
        self._sim.transitionRate(None, None, self, rate, value=value)
        return self

    def placeSpeciesNumber(self, sps, count):
        """Distribute a number of species
        
        Args:
            sps (:py:class:`~jLM.Types.Species`):
                Species type
            count (int):
                Number
        """
        self._sim.distributeNumber(sps, self, count)
        return self

    def placeSpeciesConcentration(self, sps, conc):
        """Distribute a concentration of species
        
        Args:
            sps (:py:class:`~jLM.Types.Species`):
                Species type
            conc (float):
                Concentration
        """
        self._sim.distributeConcentration(sps, self, conc)
        return self

    def addReaction(self, *args, **kwargs):
        """Add a reaction to this region

        Args:
            *args:
                Either a Reaction instance or the parameters needed to create a 
                reaction
            *kwargs:
                kwargs necessary to create new reaction
        """
        if isinstance(args[0], Reaction):
            reaction = args[0]
        else:
            reaction = self._sim.reaction(*args, **kwargs)

        self._sim.assignReaction(reaction, self)
        return self

class TrajRegion(Region):
    """Compartment, in the context of trajectory analysis"""
    def getNumberTrajectoryFromRegion(self, **kwargs):
        """Get number trajectory a species spec for this region

        See :py:meth:`~jLM.RDME.File.getNumberTrajectoryFromRegion`
        """
        return self._sim.getNumberTrajectoryFromRegion(self, region=self)

class DiffusionConst(BT.SimObj):
    """Diffusion constant"""
    def __float__(self):
        return float(self.value)

    @classmethod
    def _unique_id(cls, name, value, **kwargs):
        return name

    def _setup(self, name,  value):
        self.value = value
        self._unit = "m^{2}*s^{-1}"

    def _TeX(self):
       return "D_{" + _rmTex(self.name) + "}"


class RateConst(BT.SimObj):
    """Reaction rate constant"""
    def __float__(self):
        return float(self.value)

    @classmethod
    def _unique_id(cls, name, value, order, **kwargs):
        return name

    @property
    def stochasticRate(self):
        """Return the stochastic rate constant :math:`k_{det}*N_A*V_sv^{1-o}`"""
        return self.value*self._sim.siteNAV**self._concPower

    def _setup(self, name, value, order, stochasticRate=False):
        self.order = order
        self._concPower = 1-order

        if stochasticRate:
            self.value = value*self._sim.siteNAV**-self._concPower
        else:
            self.value = value
        concPowerStr = "" if self._concPower==1 else "^{"+str(self._concPower)+"}"
        concUnit = "" if self._concPower == 0 else "M" + concPowerStr + "*"
        self._unit = concUnit + "s^{-1}"

    def _TeX(self):
       return "k_{" + _rmTex(self.name) + "}"

    def _toLM(self):
        # XXX For some reason the rates are scaled by Nsites^(order-1) in 
        #   MpdRdmeSolver::buildDiffusionModel. This corrects for this
        return self.stochasticRate*float(self._sim.siteLattice.size)**self._concPower
    

class Reaction(BT.SimObj):
    """Chemical reaction"""
    @classmethod
    def _unique_id(cls, reactants, products, rate, **kwargs):
        return ("#".join(x.name for x in cls._parse(reactants)) +
                "@" + "#".join(x.name for x in cls._parse(products)) +
                "@" + rate.name)

    def _setup(self, reactants, products, rate):
        self.reactants = self._parse(reactants)
        self.products = self._parse(products)
        self.order = len(self.reactants)
        self.rate = rate

    @staticmethod
    def _parse(x):
        if x is None:
            return []
        elif isinstance(x, Species):
            return [x]
        elif all(isinstance(y, Species) for y in x):
            return x
        else:
            raise TypeError

    def _TeX(self):
        def side(ss):
            if len(ss) == 0:
               return r"\varnothing"
            else:
               return " + ".join(x._TeXMath() for x in ss)

        return side(self.reactants) + r"\overset{" + self.rate._TeXMath() + r"}{\longrightarrow}" + side(self.products)

class BuilderReaction(Reaction):
    """Chemical reaction, in the context of simulation construction"""
    def assignToRegion(self, region):
        """Assign to a region

        Args:
            region (:py:class:`~jLM.Types.Region`): 
                Target region
        """
        self._sim.assignReactionRegion(self, region)
        return self
