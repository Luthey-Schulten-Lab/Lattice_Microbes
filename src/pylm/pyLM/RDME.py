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

import lm
import math
import os

from  .LMLogger import *
from  .ipyInterface import *
try:
	from tqdm import tqdm
except:
	def tqdm(x,ascii=False):
		return x

## @class RDMERegion
class RDMERegion:
	"""A class that represents a type region of an RDME simulation.

    Reactions may be specified that live within a region.  In addition
    a specie's diffusion constant is region dependent and diffusion
    between regions must be specified.
    For example: cytosol, membrane, extracellular, nucleoid, etc.

    Args:
        name:
            The name of the region (e.g. cytosol)
    """
	def __init__(self, name):
		self.name=name
		self.reactions=[]
		self.defaultDiffusionRate=0.0
		self.diffusionRate={}
	
	def addReaction(self, reactant, product, rate):
		"""Adds a 0th, 1st or 2nd order reaction that can occur in the region

        Reaction rates are specified as *stochastic* rates: i.e. in units of
        :math:`1/\mathrm{s}`.

        Args:
            reactant:
                A set of reactants either as a singleton or a list
            product:
                A set of products either as a singeton or a list
            rate:
                The *stochastic* rate of the reaction.

        Returns:
            RDMERegion
        """
		if rate <= 0.0:
			LMLogger.error("In RDME.addReaction, rate must be positive")
		self.reactions.append((reactant,product,rate))
		return self

	def setDefaultDiffusionRate(self, rate):
		"""Specifies the default diffusion rate of all particles in the region 

        Args:
            rate:
                The rate of diffusion in um^2/s or um/s for 3D or 2D diffusion

        Returns:
            RDMERegion
        """
		self.defaultDiffusionRate=rate
		return self

	def setDiffusionRate(self, species, rate):
		"""Specify the diffusion rate for a particular particle type

        Args:
            species:
                The particle type
            rate:
                The rate of diffusion in um^2/s or um/s for 3D or 2D diffusion

        Returns:
            RDMERegion
        """
		self.diffusionRate[species]=rate
		return self

	def getReactionCount(self):
		"""Return the number of reactions defined in this region

        Returns:
            Get the number of reactions in the region
        """
		return len(self.reactions)

	## Create a representation that is loadable in iPython
	# @return A string containing an HTML object that can be displayed in the Jupyter notebook
	def _repr_html_(self):
		# Write out major parameters
		s  = ""
		s += "<h1> Region: %s </h1>"%self.name
		s += "Default Diffusion Rate: %f<br/>"%(self.defaultDiffusionRate)
		s += "<h2> Reactions </h2>"
		s += writeTable(["Reaction","Rate","Units"], [getReactionString(rxn[0], rxn[1], rxn[2]) for rxn in self.reactions])
		s += "<h2> Diffusion </h2>"
		s += writeTable(["Specie","Diffusion Rate"], [(specie,rate) for specie, rate in self.diffusionRate.items()])
		return s


class RDMESimulation:
	"""A class that contains all regions, reactions, diffusions and rules for a RDME simulation

	Specify a cuboid region that represents the extents to the reaction region as well as the lattice spacing

    Args:
        dimensions:
            A list of [x,y,z]
        spacing:
            Lattice spacing
        name:
            The name of the RDME simulation; default: "unnamed"
        defaultRegion:
            The name of the region that is associated with the lattice sites before any other regions are added; default:"default"
    """
	def __init__(self, dimensions, spacing, name="unnamed", defaultRegion="default"):

		# site type dictionary (necessary?)
		self.siteTypes={}
		# particleMap
		self.particleMap={}
		self.customAddedParticleList=[] # Particles that were added at specific locations
		self.species_id=[]
		self.regions={}
		self.transitionRates=[]
		self.initial_counts={}
		self.parameters={}
		
		self.continousDimensions=dimensions
		self.latticeSpacing=spacing
		self.bytesPerParticle=1

		self.addRegion(defaultRegion)

		self.lm_builder=lm.LatticeBuilder(dimensions[0], dimensions[1], dimensions[2], spacing, 1, 0)
		extraCellular=lm.Cuboid(lm.point(0.0,0.0,0.0),lm.point(*self.continousDimensions),self.siteTypes[defaultRegion])
		self.lm_builder.addRegion(extraCellular)
		extraCellular.thisown=0	# PREVENT GC
		sm=lm.SpatialModel()
		self.lm_builder.getSpatialModel(sm)

		lx,ly,lz = [int(round(cDim/self.latticeSpacing)) for cDim in self.continousDimensions]
		#self.lattice = lm.ByteLattice(lx, ly, lz, spacing, lm.getCompiledLatticeMaxOccupancy())

		self.hasBeenDiscretized = False

		# Metadata for Jupyter
		self.name = name
		self.replicates = []
		self.filename = ""

	def addRegion(self, region):
		"""Add a region to the simulation

        Args:
            region:
                The region to add to the simulation

        Returns:
            region just added
        """
		if region not in self.siteTypes:
			typeid=len(self.siteTypes)
			self.siteTypes[region]=typeid
			self.regions[region]=RDMERegion(region)
		return self.regions[region]

	def addCuboidRegion(self, name,a,b):
		"""Add a cuboid to the builder

        Args:
            name:
                Name of the site type for this region
            a:
                tuple for the first corner in continous space
            b:
                tuple for the second corner in continous space
        """
		if name not in self.siteTypes:
			raise Exception("Cannot add cuboid of an undefined region (%s)" % name)
		pa=lm.point(*a)
		pb=lm.point(*b)
		shape=lm.Cuboid(pa, pb, self.siteTypes[name])
		self.lm_builder.addRegion(shape)
		shape.thisown=0

	def addShape(self, shape):
		"""Add a Shape to the builder

        Args:
            shape:
                The Shape to add to the builder
        """
		self.lm_builder.addRegion(shape)
		shape.thisown=0

	def modifyRegion(self, region):
		"""Return a pointer to a region so that it may be modified

        Args:
            region:
                Get a region that is attached to the simulation for modification

        Returns:
            The region to modify
        """
		return self.regions[region]

	def defineSpecies(self, species):
		"""Define a specie/s of particles that exist in the simulation

        Args:
            species:
                A list of species to add to the simulation

        Returns:
            The simulation object
        """
		for s in species:
			if s not in self.species_id:
				self.species_id.append(s)
				self.particleMap[s]=len(self.species_id)
				self.initial_counts[s]=0
		return self

	def siteVolume(self):
		"""Get the actual volume of a specific site in L

        Returns:
            The volume of the site in L
        """
		return self.latticeSpacing * self.latticeSpacing * self.latticeSpacing * 1000

	def buildCapsidCell(self, length, diameter, membraneThickness, points = False):
		"""Build a capsule based shell in this RDMESimulation centered within the simulation domain that includes a membrane and cytoplasm

        Args:
            length:
                The length of the capsule from one sphere origin to the other
                diameter The diameter of the cell
            membraneThickness:
                The thickness of the membrane
            points:
                OPTIONAL: List of lists containing the coordinates of the centers 
                of the sphereoids that cap the capsid cell, e.g. [[x1, y1, z1], 
                [x2, y2, z2]]. If not given, cell is centered in the volume and 
                aligned in the z-direction

        Returns:
            The simulation object
        """
		self.addRegion('cytoplasm')
		self.addRegion('membrane')
		
		r = diameter / 2.0

		if (points != False):
			
			#test if the points lie outside the box
			if (min(points[0][0], points[1][0]) - r < 0) or (min(points[0][1], points[1][1]) - r < 0) or (min(points[0][2], points[1][2]) - r < 0):
				LMLogger.info("cell is outside the lattice volume")
				return
			if (max(points[0][0], points[1][0]) + r > self.continuousDimensions[0]) or (max(points[0][1], points[1][1]) + r > self.continuousDimensions[1]) or (max(points[0][2], points[1][2])+ r > self.continuousDimensions[2]):
				LMLogger.info("cell is outside the lattice volume")
				return
			
			p1 = lm.point(*points[0])
			p2 = lm.point(*points[1])
			
			dx,dy,dz = [points[1][i] - points[0][i] for i in range(3)]
			l = math.sqrt(dx*dx + dy*dy + dz*dz)
			if (abs(length - l) / l > 1e-2):
				LMLogger.warning("Cell length given does not match length defined by points!  Using points!")
				length = l
			
			LMLogger.info("Creating capsid with length=%g, diameter=%g, and thickness=%g", length, diameter, membraneThickness)
			
			membrane = lm.CapsuleShell(p1, p2, r - membraneThickness, r, self.siteTypes['membrane'])
			cytoplasm = lm.Capsule(p1, p2, r - membraneThickness, self.siteTypes['cytoplasm'])
			
			self.lm_builder.addRegion(membrane)
			self.lm_builder.addRegion(cytoplasm)
			
			membrane.thisown = 0 # PREVENT GC
			cytoplasm.thisown = 0 # PREVENT GC
			
		else:
			x=self.continousDimensions[0]/2.0
			y=self.continousDimensions[1]/2.0
			z=(self.continousDimensions[2]-length)/2.0
			LMLogger.info("Creating capsid with length=%g, diameter=%g, and thickness=%g", length, diameter, membraneThickness)
			p1=lm.point(x,y,z+r)
			p2=lm.point(x,y,z+length-r)
			membrane =lm.CapsuleShell(p1, p2, r-membraneThickness, r, self.siteTypes['membrane'])
			cytoplasm=lm.Capsule     (p1, p2, r-membraneThickness, self.siteTypes['cytoplasm'])
			self.lm_builder.addRegion(membrane)
			self.lm_builder.addRegion(cytoplasm)
			membrane.thisown=0	# PREVENT GC
			cytoplasm.thisown=0	# PREVENT GC
		return self

	def buildSphericalCell(self, diameter, membraneThickness, point = False):
		"""Build a spherical based shell in this RDMESimulation centered within the simulation domain that includes a membrane and cytoplasm

        Args:
            diameter:
                The diameter of the cell
            membraneThickness:
                The thickness of the membrane
            point:
                The center of the spherical cell

        Returns:
            The simulation object
        """
		self.addRegion('cytoplasm')
		self.addRegion('membrane')

		r=diameter/2.0
		
		if (point != False):
			#test if the point is inside the sim volume:
			if (point[0] - r < 0) or (point[0] + r > self.continuousDimensions[0]) or (point[1] - r < 0) or (point[1] + r > self.continuousDimensions[1]) or (point[2] - r < 0) or (point[2] + r > self.continuousDimensions[2]):
				LMLogger.info("cell is outside the lattice volume")
				return
			
			p1 = lm.point(*point)
			
			LMLogger.info("Creating spheroid with diameter=%g, and thickness=%g", diameter, membraneThickness)
			
			membrane = lm.Sphere(p1, r, self.siteTypes['membrane'])
			cytoplasm = lm.Sphere(p1, r - membraneThickness, self.siteTypes['cytoplasm'])
			
			self.lm_builder.addRegion(membrane)
			self.lm_builder.addRegion(cytoplasm)
			membrane.thisown = 0 # PREVENT GC
			cytoplasm.thisown = 0 # PREVENT GC
		
		else:
			x,y,z = [cdim/2.0 for cdim in self.continousDimensions]
			LMLogger.info("Creating spherid with diameter=%g, and thickness=%g", diameter, membraneThickness)
			p1=lm.point(x,y,z)
			membrane =lm.Sphere(p1, r, self.siteTypes['membrane'])
			cytoplasm=lm.Sphere(p1, r-membraneThickness, self.siteTypes['cytoplasm'])
			self.lm_builder.addRegion(membrane)
			self.lm_builder.addRegion(cytoplasm)
			membrane.thisown = 0 # PREVENT GC
			cytoplasm.thisown = 0 # PREVENT GC
		return self

	def addParticles(self, species='unknown', region='default', count=1):
		"""Add a specified number of particles of the specified type to the  specified region

        Args:
            species:
                The species to add to the region
            region:
                The region to add particles to
            count:
                Number of particles to add (default: 1)

        Returns:
            The simulation object
        """
		try:
			particleNum=self.particleMap[species]
		except KeyError:
			particleNum=1

		siteTypeNum=self.siteTypes[region]

		self.initial_counts[species]+=count
		self.lm_builder.addParticles(particleNum-1, siteTypeNum, count)
		return self

	def packRegion(self, region, radius, percentage, obstacleID):
		"""Add nonmoving obstacles to a particular region

        Args:
            region:
                The name of the region in which to add particles to
            radius:
                The radius of the particles
            percentage:
                The percentage of the total region volume that should be packed
            obstacleID:
                an identifier for the obstacle

        Returns:
            The simulation object
        """
		siteTypeNum=self.siteTypes[region]
		
		self.lm_builder.fillWithRandomSpheres(percentage/100.0, radius, 255-obstacleID, siteTypeNum)
		return self


	def setTransitionRate(self, species, via, to, rate):
		"""Specify the diffusion rate between species; this is a one directional rate e.g. membrane->cytosol or extracellular->membrane

        Args:
            species:
                The specie that can transition between regions
            via:
                From this region
            to:
                To this region
            rate:
                Diffusion rate between regions

        Returns:
            The simulation object
        """
		self.transitionRates.append((species, via, to, rate))
		return self

	def setTwoWayTransitionRate(self, species, one, two, rate):
		"""Specify the diffusion rate between species; this is a two directional rate 

        Args:
            species:
                The specie that can transition between regions
            one:
                A region
            two:
                The other region
            rate:
                Diffusion rate between regions

        Returns:
            The simulation object
        """
		self.setTransitionRate(species, one, two, rate)
		self.setTransitionRate(species, two, one, rate)
		return self
	
	def buildDiffusionModel(self):
		"""Return the Lattice Microbes DiffusionModel object for fine-tuning

        Returns:
            Simulation diffusion model
        """
		dm=lm.DiffusionModel()
	
		numReactions=0
		for r in self.regions:
			numReactions += self.regions[r].getReactionCount()
		LMLogger.info("number of reactions = %d", numReactions)
		dm.set_number_reactions(numReactions)

		numSpecies=len(self.particleMap)
		LMLogger.info("number of species   = %d", numSpecies)
		dm.set_number_species(numSpecies)

		numSiteTypes=len(self.siteTypes)
		LMLogger.info("number of sitetypes   = %d", numSiteTypes)
		dm.set_number_site_types(numSiteTypes)
		
		lx,ly,lz = [int(round(cDim / self.latticeSpacing)) for cDim in self.continousDimensions]
		particlesPerSite = lm.getCompiledLatticeMaxOccupancy();
		LMLogger.info("Lattice is %d x %d x %d with %g nm spacing and %d particles per site", lx, ly, lz, self.latticeSpacing*1e9, particlesPerSite)
		dm.set_lattice_x_size(lx)
		dm.set_lattice_y_size(ly)
		dm.set_lattice_z_size(lz)
		dm.set_lattice_spacing(self.latticeSpacing)
		dm.set_particles_per_site(particlesPerSite)
		dm.set_bytes_per_particle(self.bytesPerParticle)

		D=[0] * numSpecies * numSiteTypes * numSiteTypes
		RL=[0] * numSiteTypes * numReactions
		rxnum=0

		# This is about to get messy fast.  Use a function to determine the 
		# proper offset in the diffusion matrix given the site types and species
		dix=lambda _t1,_t2,_p: _t1*(numSiteTypes*numSpecies)+_t2*(numSpecies)+_p

		# Get unique reactions and the locations where they occur
		uniqueReactions = []
		uniqueReactionsLoc = []
		for region in self.regions:
			for rx in self.regions[region].reactions:
				rxnHash = (rx[0],rx[1],rx[2])
				if rxnHash in uniqueReactions:
					idx = uniqueReactions.index(rxnHash)
					uniqueReactionsLoc[idx].append(region)
				else:
					uniqueReactions.append(rxnHash)
					uniqueReactionsLoc.append( [region] )

		for k in self.regions:
			r=self.regions[k]
			LMLogger.debug("D for %s (%d)",k,self.siteTypes[k])
			for s in range(numSpecies):
				LMLogger.debug("\tD[%d][%d][%d] ([%d]) = %g", self.siteTypes[k],self.siteTypes[k],s, dix(self.siteTypes[k],self.siteTypes[k],s), r.defaultDiffusionRate)
				#D[offset+s]=r.defaultDiffusionRate
				D[dix(self.siteTypes[k],self.siteTypes[k],s)]=r.defaultDiffusionRate
			for s in r.diffusionRate:
				LMLogger.debug("\tD[%d][%d][%d] ([%d]) = %g", self.siteTypes[k],self.siteTypes[k],self.particleMap[s]-1, dix(self.siteTypes[k],self.siteTypes[k],self.particleMap[s]-1), r.diffusionRate[s])
				D[dix(self.siteTypes[k], self.siteTypes[k], (self.particleMap[s]-1))]=r.diffusionRate[s]

		for i, rx in enumerate(uniqueReactions):
			for region in uniqueReactionsLoc[i]:
				loc = self.siteTypes[region]
				LMLogger.debug("RL[%d][%d] is being set cell %d in %s",rxnum, loc, rxnum * numSiteTypes + loc, region)
				RL[rxnum * numSiteTypes + loc]=1
			rxnum+=1

		LMLogger.debug("Transitions:")
		for t in self.transitionRates:
			s=t[0]
			via=t[1]
			to=t[2]
			rate=t[3]
			LMLogger.debug("\tD[%d][%d][%d] ([%d]) = %g", self.siteTypes[via],self.siteTypes[to],self.particleMap[s]-1, dix(self.siteTypes[via],self.siteTypes[to],self.particleMap[s]-1), rate)
			D[dix(self.siteTypes[via],self.siteTypes[to],self.particleMap[s]-1)]=rate

		LMLogger.debug("D: %s", D)
		LMLogger.debug("RL: %s", RL)
		for d in range(len(D)):
			dm.add_diffusion_matrix(D[d])

		for rx in range(len(RL)):
			dm.add_reaction_location_matrix(RL[rx])

		return dm
	
	def buildReactionModel(self):
		"""Return the Lattice Microbes ReactionModel object for fine-tuning

        Returns:
            The reaction model for this simulation
        """
		rm=lm.ReactionModel()

		numReactions=0
		for r in self.regions:
			numReactions += self.regions[r].getReactionCount()
		LMLogger.info("number of reactions = %d", numReactions)
		rm.set_number_reactions(numReactions)

		numSpecies=len(self.species_id)
		LMLogger.info("number of species   = %d", numSpecies)
		rm.set_number_species(numSpecies)

		for s in self.species_id:
			LMLogger.debug("\t set initial count of %s (id %d) to %d", s, self.particleMap[s], self.initial_counts[s])
			rm.add_initial_species_count(self.initial_counts[s])

		# Build reaction matricies
		rxtypes = [0] * numReactions;
		rxconst = [0] * numReactions;
		stoich  = [0] * numSpecies * numReactions;
		depmat  = [0] * numSpecies * numReactions;

		# Get unique reactions and the locations where they occur
		uniqueReactions = []
		for region in self.regions:
			for rx in self.regions[region].reactions:
				rxnHash = (rx[0],rx[1],rx[2])
				if rxnHash in uniqueReactions:
					idx = uniqueReactions.index(rxnHash)
				else:
					uniqueReactions.append(rxnHash)

		rnum=0
		for rx in uniqueReactions:
			# Define Reaction in stiochiometry, dependency and rate matricies
			reactant=rx[0]
			product=rx[1]
			rate=rx[2]
			LMLogger.debug("\t rx %d is %s -> %s at rate %g", rnum, reactant, product, rate)
			rxconst[rnum]=rate
			if(isinstance(reactant, tuple)):
				# Second Order
				if reactant[0] != reactant[1]:
					rxtypes[rnum]=2
					for sr in reactant:
						depmat[rnum + numReactions*(self.particleMap[sr]-1)]=1
						stoich[rnum + numReactions*(self.particleMap[sr]-1)] -= 1
				elif reactant[0] == reactant[1]:
				# Second Order Self Reaction
					rxtypes[rnum]=3
					depmat[rnum + numReactions*(self.particleMap[reactant[0]]-1)]=1
					stoich[rnum + numReactions*(self.particleMap[reactant[0]]-1)] -= 2
			else:
				# Zeroth order
				if reactant == '':
					rxtypes[rnum] = 0
					#stoich[rnum + numReactions*(self.particleMap[reactant]-1)] += 0
				else:
				# First order
					rxtypes[rnum]=1
					depmat[rnum + numReactions*(self.particleMap[reactant]-1)]=1
					stoich[rnum + numReactions*(self.particleMap[reactant]-1)] -= 1
				
			if(isinstance(product, tuple)):
				for sp in product:
					stoich[rnum + numReactions*(self.particleMap[sp]-1)] += 1
			else:
				if product != '':
					stoich[rnum + numReactions*(self.particleMap[product]-1)] += 1

			rnum += 1

		LMLogger.debug("rxtypes: %s", rxtypes)
		LMLogger.debug("rxconst: %s", rxconst)

		# Populate reaction object data
		for r in range(numReactions):
			rm.add_reaction()
			rm.mutable_reaction(r).set_type(rxtypes[r])
			rm.mutable_reaction(r).add_rate_constant(rxconst[r])

		LMLogger.debug("stoich: %s", stoich)
		for x in range(len(stoich)):
			rm.add_stoichiometric_matrix(stoich[x])

		LMLogger.debug("depmat: %s", depmat)
		for x in range(len(depmat)):
			rm.add_dependency_matrix(depmat[x])

		return rm
	
	def buildSpatialModel(self):
		"""Return the Lattice Microbes SpatialModel object for fine-tuning

        Returns:
            The spatial model (i.e. obstacles, sites, etc.) for this simulation
        """
		# Retreive the spatial model from the builder
		sm=lm.SpatialModel()
		self.lm_builder.getSpatialModel(sm)
		return sm

	def setWriteInterval(self, time):
		"""Set the time interval to write data at during simulation

        Args:
            time
                The length of time between data writes
        """
		self.parameters['writeInterval']=str(time)

	def setLatticeWriteInterval(self, time):
		"""Set the time interval to write latticedata at during simulation

        Args:
             time:
                The length of time between lattice data writes
        """
		self.parameters['latticeWriteInterval']=str(time)

	def setSimulationTime(self, time):
		"""Set the simulation time length

        Args:
            time:
                The length of simulation time 
        """
		self.parameters['maxTime']=str(time)

	def setTimestep(self, time):
		"""Set the simulation time step

        Args:
            time:
                The length of simulation timestep for RDME
        """
		self.parameters['timestep']=str(time)

	def setHookInterval(self, time):
		"""Set the hook interval

        Args:
            time:
                The interval to call hookSimulation
        """
		self.parameters['hookInterval']=str(time)

	def setRandomSeed(self, seed):
		"""Set a known random seed

        Args:
            seed:
                The seed value
        """
		self.parameters['seed']=str(seed)

	def setOverflowHandler(self, mode):
		self.parameters['rdme.mpd.overflowhandler']=str(mode)

	# ## Check that each simulation parameter is reasonable, (e.g. the diffusion rates are good, the particle numbers aren't too high, the grid isn't too fine or coarse, etc.)
	# # @param self
	# def check(self):
	# 	pass

	# ###########################
	# Lattice Modification Code #
	# ###########################
	def getLattice(self, force=False):
		"""Get a discretized version of the simulation domain.  Call this after building all spherical and capsule cells

        Args:
            force:
                Force a rediscretization of the lattice. This will throw away any additional changes made after the first call.

        Returns:
            A lattice object.  This function should usually only be called once.
        """
		if (not self.hasBeenDiscretized) or force:
			self.hasBeenDiscretized=True
			self.lattice=self.allocateLattice()
			self.lm_builder.discretizeTo(self.lattice, 0, 0.0)
		return self.lattice

	def allocateLattice(self): 
		lx,ly,lz = [int(round(cDim/self.latticeSpacing)) for cDim in self.continousDimensions]
		return lm.ByteLattice(lx, ly, lz, self.latticeSpacing, lm.getCompiledLatticeMaxOccupancy())

	def setLatticeSite(self, index, siteType):
		"""Set a particular lattice site type

        Args:
            index (x,y,z):
                A list of coordinates
            siteType:
                The type to set the lattice point to. This would be the name of a region that
                has previously been performed
        """
		# Bulletproofing
		if not self.hasBeenDiscretized:
			LMLogger.error('Trying to set a lattice site before calling getLattice()')
			return
		if len(index) != 3:
			LMLogger.error('Must pass a list of form (x,y,z) that has 3 elements, given: %s',index)
			return
	
		siteNum = self.siteTypes[siteType]
		self.lattice.setSiteType(index[0], index[1], index[2], siteNum)

	def getLatticeSite(self, index):
		"""Get a particular lattice site

        Args:
            index (x,y,z):
                a list of coordinates

        Returns:
            The type of the lattice site (string)
        """
		# Bulletproofing
		if not self.hasBeenDiscretized:
			LMLogger.error('Trying to set a lattice site before calling getLattice()')
			return
		if len(index) != 3:
			LMLogger.error('Must pass a list of form (x,y,z) that has 3 elements, given: %s',index)
			return
	
		siteType = self.lattice.getSiteType(index[0], index[1], index[2])
		for i in self.siteTypes:
			if siteType == i[1]:
				return i[0]

		LMLogger.warning('Unknown site type!')
		return "unknownSiteType"

	def addParticleAt(self, index, particleType):
		"""Add a particle/ to a particular site 

        Args:
            index (x,y,z):
                a list of spatial location
            specie:
                The specie type to add
        """
		# Bulletproofing
		if not self.hasBeenDiscretized:
			LMLogger.error('Trying to set a lattice site before calling getLattice()')
			return
		if len(index) != 3:
			LMLogger.error('Must pass a list of form (x,y,z) that has 3 elements, given: %s',index)
			return

		# Check bounds and convert to indices
		ix,iy,iz = [int(round(i/self.latticeSpacing)) for i in index]
		if ix < 0 or iy < 0 or iz < 0:
			LMLogger.error('All indices must be greater than 0')
			return
		lx,ly,lz = [int(round(cDim/self.latticeSpacing)) for cDim in self.continousDimensions]
		if ix > lx or iy > ly or iz > lz:
			LMLogger.error('All indices must be less than domain extents')
			return

		particleNum=self.particleMap[particleType]
		self.initial_counts[particleType]+=1
		# TODO: add try/catch
		self.lattice.addParticle(ix, iy, iz, particleNum)
		self.customAddedParticleList.append((particleType,1)) # This list is consumed in the "save" function in the "Custom" section


	## Create an HDF5 version of the simulation amenable for later running or stand-alone running
	# @param self
	# @param filename A file to write the simulation to
	def addParticleAtIdx(self, index, particleType):
		"""Add a particle/ to a particular site 

        Args:
            index (x,y,z):
                a list of lattice site indices
            specie:
                The specie type to add
        """
		# Bulletproofing
		if not self.hasBeenDiscretized:
			LMLogger.error('Trying to set a lattice site before calling getLattice()')
			return
		if len(index) != 3:
			LMLogger.error('Must pass a list of form (x,y,z) that has 3 elements, given: %s',index)
			return

		# Check bounds
		if index[0] < 0 or index[1] < 0 or index[2] < 0:
			LMLogger.error('All indices must be greater than 0')
			return
		lx,ly,lz = [int(round(cDim/self.latticeSpacing)) for cDim in self.continousDimensions]
		if index[0] > lx or index[1] > ly or index[2] > lz:
			LMLogger.error('All indices must be less than domain extents')
			return
		
		particleNum=self.particleMap[particleType]
		self.initial_counts[particleType]+=1
		# TODO: add try/catch
		self.lattice.addParticle(index[0], index[1], index[2], particleNum)
		self.customAddedParticleList.append((particleType,1)) # This list is consumed in the "save" function in the "Custom" section


	def save(self, filename):
		"""Create an HDF5 version of the simulation amenable for later running or stand-alone running

        Args:
            filename:
                A file to write the simulation to
        """
		lm.SimulationFile.create(filename)
		f=lm.SimulationFile(filename)

		rm=self.buildReactionModel()
		f.setReactionModel(rm)

		sm=self.buildSpatialModel()
		f.setSpatialModel(sm)

		dm=self.buildDiffusionModel()
		f.setDiffusionModel(dm)

		if not self.hasBeenDiscretized:
			llattice=self.allocateLattice()
			self.lm_builder.discretizeTo(llattice, 0, 0.0)
			f.setDiffusionModelLattice(dm, llattice)
			self.hasBeenDiscretized = True
		else:
			f.setDiffusionModelLattice(dm, self.lattice)

		for key in self.parameters:
			LMLogger.debug("Setting parameter %s = %s",key, self.parameters[key])
			f.setParameter(key, self.parameters[key])
		LMLogger.debug("speciesNames: %s", ",".join(self.species_id))
		f.setParameter("speciesNames", ",".join(self.species_id))
		typenames=",".join(sorted(self.siteTypes, key=self.siteTypes.get))
		LMLogger.debug("siteTypeNames: %s", typenames)
		f.setParameter("siteTypeNames", typenames)

		f.close()

		# Custom
		# This section performs any cleanup for custom setup the user has specified
		# i.e. adding particles at specific locations
		
		# Update custom added particle counts
		if len(self.customAddedParticleList) > 0:
			try:
				import h5py
				lmFile = h5py.File(filename)
				for i in self.customAddedParticleList:
					particleNumber = self.particleMap[i[0]] - 1
					lmFile['Model']['Reaction']['InitialSpeciesCounts'][particleNumber] += i[1]
				lmFile.close() # Clean up after ourselves
			except ImportError:
				LMLogger.error("h5py is required for addParticleAt, but is missing.")

	def run(self, filename, method, replicates=1, seed=None, cudaDevices=None, checkpointInterval=0):
		"""Run the simulation using the specified solver (e.g.  NextSubvolume, MultiParticleDiffusion, etc.) for the specified amount of time

        Args:
            filename:
                The HDF file to write to
            method:
                The class name for the solver to use (e.g., lm::rdme::MpdRdmeSolver")
            replicates:
                The number of replicates to serially run
            seed:
                A seed for the random number generator to use when running the simulation; None indicates default
        """

		if cudaDevices is None:
			cudaDevices = [0]

		if seed is not None:
			f = lm.SimulationFile(filename)
			f.setParameter("seed",str(seed))
			f.close()

		self.filename = filename
		for r in tqdm(range(1, replicates+1),ascii=True):
			lm.runSimulation(filename, r, method, cudaDevices, checkpointInterval)
			# update internal state with replicates that have been run
			self.replicates.append(r)
		# Update the filename
		self.filename = filename

	def runMPI(self, filename, method, replicates=1, driver="mpirun", ppe=1, seed=None):
		"""Run the simulation using a call to mpirun with the given options

        Args:
            filename:
                The HDF file to write to
            method:
                The class name for the solver to use (e.g., lm::cme::GillespieDSolver")
            replicates:
                The number of replicates to serially run
            driver:
                The program to execute the parallel run, e.g. "mpirun", "aprun", "ibrun", etc.
            ppe:
                The number of processing elements to use  (e.g. number of nodes)
            seed:
                A seed for the random number generator to use when running the simulation; None indicates default
        """
		if seed is not None:
			f = lm.SimulationFile(filename)
			f.setParameter("seed",str(seed))
			f.close()
		# Write out system command
		repStr = "1"
		if replicates != 1:
			repStr = "1-%d"%replicates
		cmdStr = '%s -np %d  lm -sp -r %s -sl "%s"  -f %s'%(driver, ppe, repStr, method, filename)
		LMLogger.debug("Running MPI LM job with command:\n\t%s"%(cmdStr))

		# Actually run command and check return code
		status = os.system(cmdStr)
		
		# Check that LM executed correctly
		if status != 0:
			LMLogger.error("Failed running MPI (status: %d) job running command:\n\t%s\n"%(status, cmdStr))

		# Update the internal state for replicates that have been run
		for r in range(1,replicates):
			self.replicates.append(r)
		# Update the filename
		self.filename = filename
	
	def runSolver(self, filename, solver, replicates=1, cudaDevices=None, checkpointInterval=0):
		"""Run a simulation with a given solver

        Args:
            filename:
                The HDF file to write to
            solver:
                An MESolver object
            replicates:
                The number of replicates to serially run
        """
		if cudaDevices is None:
			cudaDevices = [0]
		LMLogger.debug("Running Custom RDMESolver: %s"%(solver))
		for r in tqdm(range(1, replicates+1),ascii=True):
			LMLogger.debug("  Running replicate: %d"%(r))
			lm.runSolver(filename, r, solver, cudaDevices, checkpointInterval)
			# update internal state with replicates that have been run
			self.replicates.append(r)
		# Update the filename
		self.filename = filename
		

	def buildVector(self, arg1,arg2,arg3):
		v = lm.vector(arg1,arg2,arg3)
		v.thisown=0
		return v

	def buildPoint(self, arg1,arg2,arg3):
		p = lm.point(arg1,arg2,arg3)
		p.thisown=0
		return p

	def buildSphere(self, arg1,arg2,arg3):
		s = lm.Sphere(arg1,arg2,arg3)
		s.thisown=0
		return s

	def buildEllipse(self, arg1,arg2,arg3,arg4,arg5,arg6,arg7):
		e = lm.Ellipse(arg1,arg2,arg3,arg4,arg5,arg6,arg7)
		e.thisown=0
		return e

	def buildDifference(self, arg1,arg2,arg3):
		d = lm.Difference(arg1,arg2,arg3)
		d.thisown=0
		return d

	def _repr_html_(self):
		header = '<td style="text-align:%s"><b>%s</b></td>'
		row    = '<td style="text-align:%s">%s</td>'
		# Write out major parameters
		s  = ""
		s += "<h1> RDME Simulation: %s</h1>"%(self.name)
		s += "<h2> Simulation Parameters </h2>"
		s += "Lattice Spacing: %f<br/>"%float(self.latticeSpacing)
		s += "Dimensions: (%f,%f,%f)<br/>"%(self.continousDimensions[0], self.continousDimensions[1], self.continousDimensions[2])
		for k,v in self.parameters.items():
			s += "%s: %s<br/>"%(str(k), str(v))
		# Write completed simulations
		if len(self.replicates) > 0:
			s += "<h2> Simulations </h2>"
			s += "Filename: " + str(self.filename) + "<br/>"
			s += "Completed replicates: " + str(self.replicates) + "<br/>"
		# Write out regions
		if len(self.regions.keys()) > 0:
			s += "<h2> Regions </h2>"
			for k, v in self.regions.items():
				s += "%s<br/>"%(k)
		# Write out species
		if len(self.species_id) > 0:
			s += "<h2> Species </h2>"
			s += writeTable(["Specie","Particle ID","Initial Count"], [(specie,self.particleMap[specie],self.initial_counts[specie]) for specie in self.species_id])
		# Write out reactions
		if sum([len(reg.reactions) for reg in self.regions.values()]) > 0:
			s += "<h2> Reaction Model </h2>"
			rows = []
			for name, reg in self.regions.items():
				for rxn in reg.reactions:
					rxnstr, rate, units = getReactionString(rxn[0],rxn[1],rxn[2])
					rows.append( (rxnstr,name,rate,units) )
			s += writeTable(["Reaction","Region","Rate","Units"], rows)
		# Write out diffusion
		if len(self.regions) > 0:
			s += "<h2> Diffusion Model </h2>"
			s += "<h3> Default Rates </h3>" 
			s += writeTable(["Region","Diffusion Rate"], [(name, reg.defaultDiffusionRate) for name,reg in self.regions.items()])

			s += "<h3> Species Rates </h3>" 
			rows = []
			for name, reg in self.regions.items():
				for specie, diffCoeff in reg.diffusionRate.items():
					rows.append( (name, specie, diffCoeff) )
			s += writeTable(["Region","Species","Diffusion Rate"], rows)
			if len(self.transitionRates) > 0:
				s += "<h3> Transition Rates </h3>" 
				s += writeTable(["Specie","From Region", "To Region", "Diffusion Rate"], [(tr[0],tr[1],tr[2],tr[3]) for tr in self.transitionRates])
				s += '<table>'

		# HTML string for representing the model
		return s

class IntRDMESimulation(RDMESimulation):
	
	def __init__(self, dimensions, spacing, name="unnamed", defaultRegion="default"):
		RDMESimulation.__init__(self, dimensions, spacing, name, defaultRegion)
		self.bytesPerParticle=4;

	def allocateLattice(self): 
		lx,ly,lz = [int(round(cDim/self.latticeSpacing)) for cDim in self.continousDimensions]
		return lm.IntLattice(lx, ly, lz, self.latticeSpacing, lm.getCompiledLatticeMaxOccupancy())

