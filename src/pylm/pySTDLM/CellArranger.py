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
import random
import numpy
import matplotlib.pyplot as plt
from pyLM.LMLogger import *
from .Cells import *

	

class CellArranger:
	"""A base class for laying out cells in a 3D space.

    Args:
        cellType:
            The type of cell this packer represents, can be one of: ["Sphere", "Capsule", "Box", "CapsuleShell"]
        cellAttributes:
            A dictionary with the required attributes specified for the cell and the form
        arrangerAttributes:
            A dictionary with the required/optional requirements for the individiual arranger
    """
	def __init__(self, cellType=None, cellAttributes=None, arrangerAttributes=None, name='Unknown'):
		self.allowedCellTypes=["Sphere", "Capsule", "Box", "CapsuleShell"]
		self.allowedForms=["Random", "RandomSpherical", "TightGrid", "LooseGrid", "SkewedGrid"]		

		self.cellType=cellType
		self.cellAttributes=cellAttributes
		self.cells=[]

		self.packerName = name
		self.arrangerAttributes = arrangerAttributes

	# Internal function
	def cellFactory(self):
		if self.cellType == "Sphere":
			s=SphereCell()
			return s.createModelCell(self.cellAttributes)
		elif self.cellType =="Capsule":
			c=CapsuleCell()
			return c.createModelCell(self.cellAttributes)
		elif self.cellType =="CapsuleShell":
			c=CapsuleShellCell()
			return c.createModelCell(self.cellAttributes)
		elif self.cellType =="Box":
			b=BoxCell()
			return b.createModelCell(self.cellAttributes)
		else:
			LMLogger.error("Unknown cell type: %s"%self.cellType)
			return

	
	def packVolume(self, sim, volume=None):
		"""Pack a specified volume with cells

        Args:
            volume:
                A volume specified as [[xmin,ymin,zmin],[xmax,ymax,zmax]] (TODO: Add angles)
            form:
                The form of the packing, which can one of: ["Random", "RandomSpherical", "TightGrid", "LooseGrid", "SkewedGrid"]
        Returns:
            A tuple of the form [volumePacked, volumePercent, numberCells]
        """
		if not isinstance(sim, pyLM.RDME.RDMESimulation):
			LMLogger.error("CellArranger only works for RDMESimulations.")
			return
		if self.cellType not in self.allowedCellTypes:
			LMLogger.error("No cell type specified.")
			return
		if self.packerName not in self.allowedForms:
			LMLogger.warning("Name for CellArranger unknown: %s"%self.packerName)


		# Figure out extents
		mins=[0.0,0.0,0.0]
		maxs=sim.continousDimensions
		if volume is not None:
			mins = volume[0]
			maxs = volume[1]
		if len(mins) != 3 or len(maxs) != 3:
			LMLogger.error("Incorrect volume to pack: [%s,%s]"%(mins,maxs))
			return
		
		count = self.intPackVolume(sim, mins, maxs)

		# Compute final packing statistics
		packedVol=sum([c.getVolume() for c in self.cells])
		totVol=(maxs[0]-mins[0])*(maxs[1]-mins[1])*(maxs[2]-mins[2])

		return [packedVol, packedVol/totVol*100, count]

		
	## Add packed volume to RDME simulation
	# @param self
	# @param sim An RDMESimulation object
	def addToSimulation(self, sim):
		if not isinstance(sim, pyLM.RDME.RDMESimulation):
			LMLogger.error("CellArranger only works for RDMESimulations.")
			return
		
		# Add cell to the RDMESimulation
		for c in self.cells:
			# Add the regions if they have yet to be defined
			sim.addRegion(c.cyttype)
			sim.addRegion(c.memtype)
			c.addToSimulation(sim)

# ########################## #
# Subclasses of the Arranger #
# ########################## #
## @class RandomArranger
class RandomArranger(CellArranger):
	"""An arranger that randomly orients a number of cells in a regin until either a count or a volume fraction target is met.

    This class requires either "fraction" or "number" attributes specified in the arrangerAttributes.  It also has an optional
    parameter "padding" which maintains that much padding between the cells and their neighbors.
    """
	def __init__(self, cellType, cellAttributes, arrangerAttributes):
		CellArranger.__init__(self, cellType, cellAttributes, arrangerAttributes, 'Random')
		
	# internal functions
	def newBBox(self, bbox, padding):
		minsloc = bbox[0]
		maxsloc = bbox[1]
		
		newmins = [minsloc[i] - padding[i] for i in range(3)]
		newmaxs = [maxsloc[i] + padding[i] for i in range(3)]
	
		return [newmins, newmaxs]

	def doesIntersect(self, b1, b2):
		mid1 = [(b1[1][i] - b1[0][i])/2.0 + b1[0][i] for i in range(3)]
		mid2 = [(b2[1][i] - b1[0][i])/2.0 + b2[0][i] for i in range(3)]
		dim1 = [abs((b1[1][i] - b1[0][i])/2.0) for i in range(3)]	
		dim2 = [abs((b2[1][i] - b2[0][i])/2.0) for i in range(3)]	

		# distance between centers
		dist = [abs(mid2[i] - mid1[i]) for i in range(3)]

		if (dim1[0] + dim2[0] <= dist[0]) and (dim1[1] + dim2[1] <= dist[1]) and (dim1[2] + dim2[2] <= dist[2]):
			return True
		return False

	def outsideDomain(self, b1, mins, maxs):
		if b1[0][0] < mins[0] or b1[1][0] > maxs[0] or b1[0][1] < mins[1] or b1[1][1] > maxs[1] or b1[0][2] < mins[2] or b1[1][2] > maxs[2]:
			print("outside")
			return True	
		return False	
		

	def intPackVolume(self, sim, mins, maxs):

		# Pack volume
		count=0
		# Boudning box
		tempCell=self.cellFactory()
		bbox=tempCell.boundingBox()

		# Determine number of cells in each direction
		dimCell=[abs(bbox[0][i]-bbox[1][i]) for i in range(3)]
		dimDomain=[abs(maxs[i] - mins[i]) for i in range(3)]

		# bulletproofing
		if "fraction" not in self.arrangerAttributes and "number" not in self.arrangerAttributes:
			LMLogger.error("fraction or number must be specified for the RandomArranger")
			return
		fraction = 0.
		if "fraction" in self.arrangerAttributes:
			fraction = float(self.arrangerAttributes["fraction"])
			if fraction > 1.0 or fraction < 0.:
				LMLogger.error("fraction of volume to fill must be in range [0,1]")
				return
		number = 0
		if "number" in self.arrangerAttributes:
			number = int(self.arrangerAttributes["number"])
			if number < 0:
				LMLogger.error("number of cells to add must be positive")
				return
		padding=[0.,0.,0.]
		if "padding" in self.arrangerAttributes:
			padding = self.arrangerAttributes["padding"]
		if len(padding) != 3:
			LMLogger.error("padding must be a list of form [x,y,z]")
			return


		# begin Monte carlo
		notDone=True
		volAdded = 0.
		totVol = dimDomain[0] * dimDomain[1] * dimDomain[2]
		while notDone:
			# escape
			if number != 0:
				if count >= number:
					notDone = False
					break
			elif fraction != 0.0:
				if totVol*fraction <= volAdded:
					notDone = False
					break

			# get a random (x,y,z)
			pos = [random.uniform(mins[i],maxs[i]) for i in range(3)]
			#print pos		
	
			# set temporary cell to work in that volume
			trialCell = copy.deepcopy(tempCell)
			trialCell.translateCell([pos[0], pos[1], pos[2]])
			intersects = False
			for i in self.cells:
				if self.doesIntersect(self.newBBox(i.boundingBox(), padding), self.newBBox(trialCell.boundingBox(), padding)):
					intersects = True
					break

			if not intersects:
				if not self.outsideDomain(self.newBBox(trialCell.boundingBox(),padding), mins, maxs):
					self.cells.append(trialCell)
					count += 1
					volAdded += trialCell.getVolume()			
					LMLogger.info( "adding %g,%g,%g vol: %g"%(pos[0], pos[1], pos[2], trialCell.getVolume()))

		return count


class RandomSphericalArranger(CellArranger):
	"""An arranger that randomly orients a number of cells in a sphere inside a regin until either a count or a volume fraction target is met.

    This class requires either "fraction" or "number" attributes specified in the arrangerAttributes.  It also has an optional
    parameter "padding" which maintains that much padding between the cells and their neighbors.
    """
	def __init__(self, cellType, cellAttributes, arrangerAttributes):
		CellArranger.__init__(self, cellType, cellAttributes, arrangerAttributes, 'RandomSpherical')

	def intPackVolume(self, sim, mins, maxs):
		# Pack volume
		count=0
		# Boudning box
		tempCell=self.cellFactory()
		bbox=tempCell.boundingBox()

		# Determine number of cells in each direction
		dimCell=[abs(bbox[0][i]-bbox[1][i]) for i in range(3)]
		dimDomain=[abs(maxs[i] - mins[i]) for i in range(3)]

		return count


class TightGridArranger(CellArranger):
	"""An arranger that orients cells in square grid."""
	def __init__(self, cellType, cellAttributes, arrangerAttributes):
		CellArranger.__init__(self, cellType, cellAttributes, arrangerAttributes, 'TightGrid')

	def intPackVolume(self, sim, mins, maxs):
		# Pack volume
		count=0
		# Boudning box
		tempCell=self.cellFactory()
		bbox=tempCell.boundingBox()

		# Determine number of cells in each direction
		dimCell=[abs(bbox[0][i]-bbox[1][i]) for i in range(3)]
		dimDomain=[abs(maxs[i] - mins[i]) for i in range(3)]
		cellCounts=[int(math.floor(dimDomain[i]/dimCell[i])) for i in range(3)]
		
		for x in range(cellCounts[0]):
			for y in range(cellCounts[1]):
				for z in range(cellCounts[2]):
					# Add to internal representation
					addCell=copy.deepcopy(tempCell)
					addCell.translateCell([x*dimCell[0], y*dimCell[1], z*dimCell[2]])
					self.cells.append(addCell)
					count+=1

		return count


class LooseGridArranger(CellArranger):
	"""An arranger that orients cells in square grid with padding in between. 
    
    The attribute dictionary requires "offset" which is a list of the form [dx, dy, dz]."""
	def __init__(self, cellType, cellAttributes, arrangerAttributes):
		CellArranger.__init__(self, cellType, cellAttributes, arrangerAttributes, 'LooseGrid')

	def intPackVolume(self, sim, mins, maxs):
		# Pack volume
		count=0
		# Boudning box
		tempCell=self.cellFactory()
		bbox=tempCell.boundingBox()

		# Determine number of cells in each direction
		dimCell=[abs(bbox[0][i]-bbox[1][i]) for i in range(3)]
		dimDomain=[abs(maxs[i] - mins[i]) for i in range(3)]

		if "offset" not in self.cellAttributes:
			LMLogger.error("offset must be specified for LooseGrid arranger")
			return
		offset=self.cellAttributes['offset']
		if len(offset) != 3:
			LMLogger.error("offset must a vector [x,y,z]")
			return
		
		# compute number of cell in each direction
		cellCounts=[int(math.floor((dimDomain[i]+offset[i])/(dimCell[i]+offset[i]))) for i in range(3)]
		for x in range(cellCounts[0]):
			for y in range(cellCounts[1]):
				for z in range(cellCounts[2]):
					# Add to internal representation
					addCell=copy.deepcopy(tempCell)
					addCell.translateCell([x*(dimCell[0]+offset[0]), y*(dimCell[1]+offset[1]), z*(dimCell[2]+offset[2])])
					self.cells.append(addCell)
					count+=1

		return count

#
#class SkewedGridArranger(CellArranger):
#	def __init__(self, cellType, cellAttributes, arrangerAttributes):
#		CellArranger.__init__(self, cellType, cellAttributes, arrangerAttributes, 'SkewedGrid')
#
#	def intPackVolume(self, sim, mins, maxs):
#		# Pack volume
#		count=0
#		# Boudning box
#		tempCell=self.cellFactory()
#		bbox=tempCell.boundingBox()
#
#		# Determine number of cells in each direction
#		dimCell=[abs(bbox[0][i]-bbox[1][i]) for i in range(3)]
#		dimDomain=[abs(maxs[i] - mins[i]) for i in range(3)]
#
#		if "skew" not in self.cellAttributes.keys():
#			LMLogger.error("skew must be specified for a SkewedGrid arranger")
#			return
#		skew=self.cellAttributes['skew']
#		if len(skew) != 3:
#			LMLogger.error("skew must a vector [x,y,z]")
#			return
#
#		# compute the number of cells in each direction
#		cellCounts=[int(math.floor(dimDomain[i]/dimCell[i])) for i in range(3)]
#		
#		fracx=0.0
#		fracy=0.0
#		fracz=0.0
#		for z in range(cellCounts[2]):
#			for y in range(cellCounts[1]):
#				for x in range(cellCounts[0]):
#					# Add to internal representation
#					addCell=copy.deepcopy(tempCell)
#					addCell.translateCell([x*(dimCell[0]+fracx), y*(dimCell[1]+fracy), z*(dimCell[2])])
#					self.cells.append(addCell)
#					count+=1
#					print "adding cell %d"%count
#				fracx+=skew[0]
#				if fracx > 1.0:
#					fracx-=1.0
#			fracy+=skew[1]
#			if fracy > 1.0:
#				fracy-=1.0
#						
#		return count
#


