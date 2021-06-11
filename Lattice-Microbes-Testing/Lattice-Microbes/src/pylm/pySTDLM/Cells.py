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
import pyLM
import h5py
import math
import copy
import numpy
import matplotlib.pyplot as plt
from pyLM.LMLogger import *


class CellShape:
	"""Base class for a particular cell type for use in CellArranger (and later more)"""
	def __init__(self):
		self.name="Unknown"
		self.volume=-1.0
		self.mins=3*[0.]
		self.maxs=3*[0.]
		self.membraneThickness=0.0

		self.memtype="membrane"
		self.cyttype="cytoplasm"

	def getVolume(self):
		"""Volume occupied by cell

        When overriding this class, you must specify a "computeVolume()" function of no arguments

        Returns:
            Volume in intrinsic units
        """
		if self.volume <= 0.0:
			self.volume=self.computeVolume()
			if self.volume < 0.0:
				LMLogger.error("Volume computed was negative... call the developers.")
		return self.volume

	# TODO: Make it return an orientation as well
	def boundingBox(self):
		"""Return a bounding box for the cell"""
		return [self.mins, self.maxs]

	def addToSimulation(self, sim):
		"""Add the cell to the simulation

        Args:
            sim:
                An RDMESimulation object
        """
		self.internalAddToSim(sim)

	def translateCell(self, point):
		"""Shift a cell in space by the specified amount

        Args:
            point [x,y,z]:
                translation in space
        """
		self.internalTranslateCell(point)

	def createModelCell(self, attr):
		"""Create model cell

        Args:
            attrs:
                A dictionary of attributes
        """
		return self.internalModelCell(attr)

	def setRegions(self, membrane, cytoplasm):
		"""Set regions

        Args:
            membrane:
                The name of the region in which the membrane should be considered
            cytoplasm:
                The name of the region in which the cytoplasm should be considered
        """
		self.memtype = membrane
		self.cyttype = cytoplasm

class SphereCell(CellShape):
	'''A representation for a spherical cell. This cell requires an attribute dictionary that includes "radius" and "membraneThickness"'''
	def __init__(self, origin=3*[0.], radius=0.):
		CellShape.__init__(self)
		self.name="Sphere"
		self.origin=origin
		if len(origin) != 3:
			LMLogger.error("Must specify an origin for the cell in form [x,y,z].")
			return
		self.radius=radius
		if self.radius < 0.0:
			LMLogger.error("Radius must be positive.")

		self.mins = [x-self.radius for x in origin]
		self.maxs = [x+self.radius for x in origin]

	def internalModelCell(self, attr):
		if "radius" not in attr:
			LMLogger.error("When creating a model sphere cell, radius must be specified.")
			return
		if "membraneThickness" not in attr:
			LMLogger.error("When creating a model sphere cell, membraneThickness must be specified.")
			return

		# Return a model cell
		self.radius = float(attr['radius'])
		self.origin = [self.radius,self.radius,self.radius]
		self.membraneThickness = attr['membraneThickness']
		self.mins = [x-self.radius for x in self.origin]
		self.maxs = [x+self.radius for x in self.origin]
		return self


	def computeVolume(self):
		return self.radius*self.radius*self.radius*numpy.pi*4.0/3.0

	def internalAddToSim(self, sim):
		p1 = lm.point(*self.origin)
		LMLogger.info("Creating spheroid with diameter=%g, and thickness=%g", 2.0*self.radius, self.membraneThickness)
		
		membrane = lm.Sphere(p1, self.radius, sim.siteTypes[self.memtype])
		cytoplasm = lm.Sphere(p1, self.radius - self.membraneThickness, sim.siteTypes[self.cyttype])
		
		sim.lm_builder.addRegion(membrane)
		sim.lm_builder.addRegion(cytoplasm)

		# Prevent GC
		membrane.thisown=0
		cytoplasm.thisown=0

	def internalTranslateCell(self, point):
		self.origin = [self.origin[i] + point[i] for i in range(3)]

class BoxCell(CellShape):
	'''A representation for a spherical cell. This cell requires an attribute dictionary that includes "width", "height", "depth" and "membraneThickness"'''
	def __init__(self, p1=3*[0.], p2=3*[0.]):
		CellShape.__init__(self)
		self.name="Box"
		self.p1=p1
		self.p2=p2
		if len(p1) != 3 or len(p2):
			LMLogger.error("Must specify an points for the cell in form [x,y,z].")
			return

		self.mins = p1
		self.maxs = p2

		self.width=p2[0]-p1[0]
		self.height=p2[1]-p1[1]
		self.depth=p2[2]-p1[2]

	def internalModelCell(self, attr):
		if "width" not in attr:
			LMLogger.error("When creating a model box cell, width must be specified.")
			return
		if "height" not in attr:
			LMLogger.error("When creating a model box cell, heigh must be specified.")
			return
		if "depth" not in attr:
			LMLogger.error("When creating a model box cell, depth must be specified.")
			return
		if "membraneThickness" not in attr:
			LMLogger.error("When creating a model sphere cell, membraneThickness must be specified.")
			return

		# Return a model cell
		self.width = float(attr['width'])
		self.height = float(attr['height'])
		self.depth = float(attr['depth'])
		self.p1 = 3*[0.]
		self.p2 = [self.width,self.height,self.depth]
		self.membraneThickness = float(attr['membraneThickness'])
		self.mins = self.p1
		self.maxs = self.p2
		return self


	def computeVolume(self):
		return self.width*self.height*self.depth

	def internalAddToSim(self, sim):
		LMLogger.info("Creating box with p1=%s, p2=%s, and thickness=%g", self.p1, self.p2, self.membraneThickness)
		
		p1 = lm.point(*self.p1)
		p2 = lm.point(*self.p2)
		p3 = lm.point(*[x+self.membraneThickness for x in self.p1])
		p4 = lm.point(*[x-self.membraneThickness for x in self.p2])
		membrane = lm.Cuboid(p1, p2, sim.siteTypes[self.memtype])
		cytoplasm = lm.Cuboid(p3, p4, sim.siteTypes[self.cyttype])
		
		sim.lm_builder.addRegion(membrane)
		sim.lm_builder.addRegion(cytoplasm)

		# Prevent GC
		membrane.thisown=0
		cytoplasm.thisown=0

	def internalTranslateCell(self, point):
		self.p1 = [self.p1[i] + point[i] for i in range(3)]
		self.p2 = [self.p2[i] + point[i] for i in range(3)]

class CapsuleCell(CellShape):
	'''A representation for a spherical cell. This cell requires an attribute dictionary that includes "radius", "length", and "membraneThickness"'''
	def __init__(self, p1=3*[0.], p2=3*[0.], radius=0., length=0.):
		CellShape.__init__(self)
		self.name="Capsule"
		self.p1=p1
		self.p2=p2
		if len(p1) != 3 or len(p2) != 3:
			LMLogger.error("Must specify an points for the cell in form [x,y,z].")
			return
		self.radius=radius
		if self.radius < 0.0:
			LMLogger.error("Radius must be positive.")
		self.length=length
		if self.length < 0.0:
			LMLogger.error("Length must be positive.")

		min1 = [x-self.radius for x in p1]
		max1 = [x+self.radius for x in p1]
		min2 = [x-self.radius for x in p2]
		max2 = [x+self.radius for x in p2]
		self.mins = [min(min1[i],min2[i]) for i in range(3)]
		self.maxs = [max(max1[i],max2[i]) for i in range(3)]

	def internalModelCell(self, attr):
		if "radius" not in attr:
			LMLogger.error("When creating a model sphere cell, radius must be specified.")
			return
		if "length" not in attr:
			LMLogger.error("When creating a model sphere cell, length must be specified.")
			return
		if "membraneThickness" not in attr:
			LMLogger.error("When creating a model sphere cell, membraneThickness must be specified.")
			return

		# Return a model cell
		self.radius = float(attr['radius'])
		self.length = float(attr['length'])
		self.p1 = 3*[self.radius]
		self.p2 = [self.radius, self.radius, self.radius+(self.length-2.0*self.radius)]

		self.membraneThickness = float(attr['membraneThickness'])
		min1 = [x-self.radius for x in self.p1]
		max1 = [x+self.radius for x in self.p1]
		min2 = [x-self.radius for x in self.p2]
		max2 = [x+self.radius for x in self.p2]
		self.mins = [min(min1[i],min2[i]) for i in range(3)]
		self.maxs = [max(max1[i],max2[i]) for i in range(3)]
		return self


	def computeVolume(self):
		return numpy.pi*self.radius*self.radius*(self.length-2.0*self.radius) + 4.0/3.0*numpy.pi*self.radius*self.radius*self.radius

	def internalAddToSim(self, sim):
		LMLogger.info("Creating capsule with radius=%g, length=%g, and thickness=%g", self.radius, self.length, self.membraneThickness)

		p1 = lm.point(*self.p1)
		p2 = lm.point(*self.p2)
		membrane = lm.CapsuleShell(p1, p2, self.radius - self.membraneThickness, self.radius, sim.siteTypes[self.memtype])
		cytoplasm = lm.Capsule(p1, p2, self.radius - self.membraneThickness, sim.siteTypes[self.cyttype])
		
		sim.lm_builder.addRegion(membrane)
		sim.lm_builder.addRegion(cytoplasm)

		# Prevent GC
		membrane.thisown=0
		cytoplasm.thisown=0

	def internalTranslateCell(self, point):
		self.p1 = [self.p1[i] + point[i] for i in range(3)]
		self.p2 = [self.p2[i] + point[i] for i in range(3)]



class CapsuleShellCell(CapsuleCell):
	'''A representation for a spherical cell. This cell requires an attribute dictionary that includes "radius", "length" and "membraneThickness"'''
	def __init__(self, p1=3*[0.], p2=3*[0.], radius=0., length=0.):
		CapsuleCell.__init__(self, p1, p2, radius, length)
		self.name="CapsuleShell"

	def internalModelCell(self, attr):
		if "radius" not in attr:
			LMLogger.error("When creating a model sphere cell, radius must be specified.")
			return
		if "length" not in attr:
			LMLogger.error("When creating a model sphere cell, length must be specified.")
			return
		if "membraneThickness" not in attr:
			LMLogger.error("When creating a model sphere cell, membraneThickness must be specified.")
			return

		# Return a model cell
		self.radius = float(attr['radius'])
		self.length = float(attr['length'])
		self.p1 = 3*[self.radius]
		self.p2 = [self.radius, self.radius, self.radius+(self.length-2.0*self.radius)]

		self.membraneThickness = float(attr['membraneThickness'])
		min1 = [x-self.radius for x in self.p1]
		max1 = [x+self.radius for x in self.p1]
		min2 = [x-self.radius for x in self.p2]
		max2 = [x+self.radius for x in self.p2]
		self.mins = [min(min1[i],min2[i]) for i in range(3)]
		self.maxs = [max(max1[i],max2[i]) for i in range(3)]
		return self


	def computeVolume(self):
		rin = self.radius-self.membraneThickness
		return numpy.pi*self.radius*self.radius*(self.length-2.0*self.radius) + 4.0/3.0*numpy.pi*self.radius*self.radius*self.radius - numpy.pi*rin*rin*(self.length-2.0*rin) + 4.0/3.0*numpy.pi*rin*rin*rin

	def internalAddToSim(self, sim):
		LMLogger.info("Creating capsule with radius=%g, length=%g, and thickness=%g", self.radius, self.length, self.membraneThickness)

		p1 = lm.point(*self.p1)
		p2 = lm.point(*self.p2)
		membrane = lm.CapsuleShell(p1, p2, self.radius - self.membraneThickness, self.radius, sim.siteTypes[self.memtype])
		
		sim.lm_builder.addRegion(membrane)

		# Prevent GC
		membrane.thisown=0

