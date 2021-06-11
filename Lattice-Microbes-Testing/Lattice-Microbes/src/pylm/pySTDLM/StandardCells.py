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

import pyLM
from pyLM.RDME import *
from pyLM.units import *

# #########
# E. coli #
# #########
def buildEColiCell(sim, crowded=False, crowdedMembrane=False ):
	"""Build a standard E coli cell 0.5 microns across and 2 microns long.

    Args:
        sim:
            An RDMESimulation object
        crowded:
            Should the cytosol be crowded (default: false)
        crowdedMembrane:
            Should the membrane be crowded (default: false)
    """
	# Sanity check
	if not isinstance(sim, RDMESimulation):
		LMLogger.error('Could not build an E. coli cell in the sim which is not of type \'RDMESimulation\'')

	# Create the cell
	sim.buildCapsidCell(length=micron(2), diameter=micron(0.5), membraneThickness=nm(32))

	# Pack the cell
	if crowded is True:
		packFastGrowingEcoli(sim)

	# Pack the membrane
	if crowdedMembrane is True:
		pass
	
	# Give the simulation back
	return sim

def buildDividingEColiCell(sim, crowded=False, crowdedMembrane=False ):
	"""Build a standard E coli cell 0.5 microns across and 4 microns long.

    Args:
        sim:
            An RDMESimulation object
        crowded:
            Should the cytosol be crowded (default: false)
        crowdedMembrane:
            Should the membrane be crowded (default: false)
    """
	# Sanity check
	if not isinstance(sim, RDMESimulation):
		LMLogger.error('Could not build an E. coli cell in the sim which is not of type \'RDMESimulation\'')

	# Create the cell
	sim.buildCapsidCell(length=micron(4), diameter=micron(0.5), membraneThickness=nm(32))

	# Pack the cell
	if crowded is True:
		packFastGrowingEcoli(sim)

	# Pack the membrane
	if crowdedMembrane is True:
		pass
	
	# Give the simulation back
	return sim


def buildSphericalEColiCell(crowded=False, latticeSpacing=16, crowdedMembrane=False ):
	"""Build a Spherical E coli cell, one that lacks the gene for elongation.

    Args:
        sim:
            An RDMESimulation object
        crowded:
            Should the cytosol be crowded (default: false)
        crowdedMembrane:
            Should the membrane be crowded (default: false)
    """
	# Sanity check
	if not isinstance(sim, RDMESimulation):
		LMLogger.error('Could not build an E. coli cell in the sim which is not of type \'RDMESimulation\'')

	# Create the cell
	sim.buildSphericalCell(diameter=micron(0.5), membraneThickness=nm(32))

	# Pack the cell
	if crowded is True:
		packFastGrowingEcoli(sim)

	# Pack the membrane
	if crowdedMembrane is True:
		pass
	
	# Give the simulation back
	return sim


def buildFilamentousEColiCell(sim, length=8.0, crowded=False, crowdedMembrane=False ):
	"""Build a long filamentous E coli cell where the user can specify the length.

    Args:
        sim:
            An RDMESimulation object
        length:
            The length of the cell in microns (default: 8)
        crowded:
            Should the cytosol be crowded (default: false)
        crowdedMembrane:
            Should the membrane be crowded (default: false)
    """
	# Sanity check
	if not isinstance(sim, RDMESimulation):
		LMLogger.error('Could not build an E. coli cell in the sim which is not of type \'RDMESimulation\'')

	# Create the cell
	sim.buildCapsidCell(length=micron(length), diameter=micron(0.5), membraneThickness=nm(32))

	# Pack the cell
	if crowded is True:
		packFastGrowingEcoli(sim)

	# Pack the membrane
	if crowdedMembrane is True:
		pass
	
	# Give the simulation back
	return sim

def packFastGrowingEcoli(sim):
	"""Pack a cell with the protein distribution for fast growing E coli.

    Args:
        sim:
            A RDMESimulation object
    """
	if not isinstance(sim, RDMESimulation):
		return

	# get cytosol region
	region='cytoplasm'

	# pack the cell
	sim.packRegion(region, nm(10.4), 17.8, 1)
	sim.packRegion(region, nm(5.2),  18.6, 2)
	sim.packRegion(region, nm(4.3),  0.7, 3)
	sim.packRegion(region, nm(4.1),  0.3, 4)
	sim.packRegion(region, nm(4.0),  1.7, 5)
	sim.packRegion(region, nm(3.8),  1.2, 6)
	sim.packRegion(region, nm(3.5),  0.9, 7)
	sim.packRegion(region, nm(3.4),  2.5, 8)
	sim.packRegion(region, nm(3.0),  2.0, 9)
	sim.packRegion(region, nm(2.7),  2.0, 10)
	sim.packRegion(region, nm(2.3),  1.8, 11)
	sim.packRegion(region, nm(1.7),  0.4, 12)

	return sim

def packSlowGrowingEcoli(sim):
	"""Pack a cell with the protein distribution for fast growing E coli.

    Args:
        sim:
            A RDMESimulation object
    """
	if not isinstance(sim, RDMESimulation):
		return

	# get cytosol region
	region='cytoplasm'
	
	# pack the cell
	sim.packRegion(region, nm(10.4), 5.7, 1)
	sim.packRegion(region, nm(5.2),  25.5, 2)
	sim.packRegion(region, nm(4.3),  1.2, 3)
	sim.packRegion(region, nm(4.1),  0.4, 4)
	sim.packRegion(region, nm(4.0),  2.4, 5)
	sim.packRegion(region, nm(3.8),  1.8, 6)
	sim.packRegion(region, nm(3.5),  1.3, 7)
	sim.packRegion(region, nm(3.4),  3.6, 8)
	sim.packRegion(region, nm(3.0),  2.8, 9)
	sim.packRegion(region, nm(2.7),  2.7, 10)
	sim.packRegion(region, nm(2.3),  2.3, 11)
	sim.packRegion(region, nm(1.7),  0.6, 12)

	return sim

