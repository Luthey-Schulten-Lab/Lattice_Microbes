# 
# University of Illinois Open Source License
# Copyright 2008-2018 Luthey-Schulten Group,
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
# Author(s): Michael J. Hallock and Joseph R. Peterson
# 
#

import sys
import os
import random
import pyLM
from pyLM.LMLogger import *


## @class ParticleGenerator
class ParticleGenerator:
	"""A particle generator that adds particles to a particular region of a simulation.

    Args:
        gen:
            A PDF set for x,y,z that on the domain [-1,1] and range [0,1) where x is the range of interest
        ptype:
            The name of the particle
    """
	def __init__(self, gen, ptype):
		self.generator=gen
		if len(self.generator) != 3:
			LMLogger.error("Must pass a list of generators: [xgen,ygen,zgen]")
			raise ValueError("Axial generator list expected")	
		self.particlename = ptype
		self.totalAdded = 0
		self.addedLast = 0

	def addParticlesToRegion(self,rdmesim, region, number, directions, spatial = None):
		"""Add a number of particles to a region based on Monte-Carlo sampling of the generator

        Args:
            rdmesim:
                An RDME simulation object that has already been "discretized"
            region:
                A region name of the simulation to add particles to
            number:
                The number of the particle to add
            directions:
                A list of the form (x,y,z) where a 0 indicates a uniform distribution in 
                that direction and a 1 indicates the particles are sampled from the generator in 
                that direction 
            spatial (OPTIONAL):
                A definition of a subsection of the domain that this generation should occur in 
                as a list of lists [[xl,yl,zl],[xh,yh,zh]], otherwise the generator is scaled 
                from (0,0,0) to domain extends (w,h,d)
        """
		# bulletproofing
		if not isinstance(rdmesim, pyLM.RDME.RDMESimulation):
			LMLogger.error("addParticlesToRegion requires that an RDMESimulation be specified, given type: %s", rdmesim.__class__.__name__)
			raise TypeError("Expected RDMESimulation")
		if rdmesim.hasBeenDiscretized == False:
			LMLogger.error("Must discretize region before adding particles from a generator")
			raise Exception("Discretization required")
		if number < 1:
			LMLogger.warning('Must add at least one particle, given: %d', number)
			return
		if len(directions) != 3:
			LMLogger.error("Must specify a direction of form (x,y,z), given: %s", directions)
			raise ValueError("3-tuple expected")	
	
		# Get the dimensions to normalize distribution to
		mins=3*[0.]
		maxs=[d for d in rdmesim.continuousDimensions]
		if spatial != None:
			if len(spatial) != 2:
				LMLogger.error("spatial argument must be list of list: [[x,y,z],[xh,yh,xh]]")
				raise ValueError("bounds expected")	
			if len(spatial[0]) !=3 or len(spatial[1]) !=3:
				LMLogger.error("spatial must be list of lists: [[x,y,z],[xh,yh,xh]]")
				raise ValueError("bounds expected")	
			mins=spatial[0]
			maxs=spatial[1]


		self.addedLast=0
		# Get the particle number
		ptype=rdmesim.particleMap[self.particlename]		
		stype = rdmesim.siteTypes[region]
		lat = rdmesim.lattice
		latspace=rdmesim.latticeSpacing
	
		dels = [maxs[i]-mins[i] for i in range(3)]

		# start adding particles
		while self.addedLast < number:
			# generate randoms
			xlocal=random.uniform(0.0,1.0)
			ylocal=random.uniform(0.0,1.0)
			zlocal=random.uniform(0.0,1.0)
			if directions[0] == 1:
				good=False
				while not good:
					rlocal = random.uniform(0,1.0)
					loc = self.generator[0](xlocal)
					LMLogger.info("x trying: %f with: %f"%(rlocal,loc))
					if rlocal < loc:
						good=True
					else:
						xlocal=random.uniform(0.0,1.0)
			if directions[1] == 1:
				good=False
				while not good:
					rlocal = random.uniform(0,1.0)
					loc = self.generator[1](ylocal)
					LMLogger.info("y trying: %f with: %f"%(rlocal,loc))
					if rlocal < loc:
						good=True
					else:
						ylocal=random.uniform(0.0,1.0)
			if directions[2] == 1:
				good=False
				while not good:
					rlocal = random.uniform(0,1.0)
					loc = self.generator[2](zlocal)
					LMLogger.info("z trying: %f with: %f"%(rlocal,loc))
					if rlocal < loc:
						good=True
					else:
						zlocal=random.uniform(0.0,1.0)
		
			LMLogger.info("trying (%f,%f,%f)"%(xlocal, ylocal, zlocal))
	
			# convert to ints
			xint = int(xlocal*dels[0]/latspace)
			yint = int(ylocal*dels[1]/latspace)
			zint = int(zlocal*dels[2]/latspace)
			LMLogger.info("loc (%d,%d,%d)"%(xint, yint, zint))

			if lat.getSiteType(xint,yint,zint) == stype:
				# TODO: Add bulletproofing if a exception of type (InvalidSiteException,InvalidParticleException) is raised
				LMLogger.info("Added a particle of type %s to (%d,%d,%d)", self.particlename, xint, yint, zint)
				lat.addParticle(xint,yint,zint, ptype)
				self.addedLast += 1

		# Count total number added
		self.totalAdded += self.addedLast
		rdmesim.customAddedParticleList.append((self.particlename,self.totalAdded))
		return


