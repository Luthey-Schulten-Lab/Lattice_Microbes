from pyLM import *
from pyLM.units import *

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()

import pySTDLM.StandardCells as sc

latticeSpacing = 16 #nm

# Create our simulation object
sim=RDME.RDMESimulation(dimensions=micron(1.024,1.024,4.096), spacing=nm(latticeSpacing))

# Build a capsid cell
#  This is a built in function that builds a spherocylinder (pill shaped) cell 
#  in the center of the domain with the specified length (total end to end) and
#  the diameter.  The membrane thickness is also specified.  This may need to 
#  be increased to about 2x the lattice spacing, otherwise holes can form in the
#  membrane.  The built in functions "buildCapsidCell" and "buildSphericalCell"
#  create 2 additional regions named "membrane" and "cytoplasm" which represent
#  what their name implies.  The extracellular space is the default region type.
sim.buildCapsidCell(length=micron(4), diameter=micron(1), membraneThickness=nm(32))

# This command packs obstacles into the E. coli cell based on 
#  a default packing used in the 2012 paper by Roberts et al.
sc.packFastGrowingEcoli(sim)

# define our chemical species
species = ['minDatp', 'minDadp', 'minE', 'minDm', 'minDEm']
sim.defineSpecies(species)

# Modify the cytoplasm to add diffusion rates and reactions
sim.modifyRegion('cytoplasm') \
	.setDefaultDiffusionRate(2.5e-12) \
	.addReaction(reactant='minDadp', product='minDatp', rate=0.5)

# Modify the membrane to add reactions
#  Here the diffusion rate of 'minDm' and 'minDEm' on the
#  membrane is set specifically, demonstrating how diffusion
#  coefficients can be customized for each region.  
sim.modifyRegion('membrane') \
	.setDefaultDiffusionRate(2.5e-12) \
	.setDiffusionRate(species='minDm',  rate=1e-14) \
	.setDiffusionRate(species='minDEm', rate=1e-14) \
	.addReaction(reactant='minDatp', product='minDm', rate=1e-12) \
	.addReaction(reactant=('minDatp','minDm'), product=('minDm','minDm'), rate=9e6) \
	.addReaction(reactant=('minDm','minE'), product='minDEm', rate=5.56e7) \
	.addReaction(reactant='minDEm', product=('minE','minDadp'), rate=0.7) 

# Set diffusive properties between regions
#  You can set the transition rate between two regions by calling "setTransitionRate"
#  on the simulation object, specifying the species and the source and destination regions.
#  This is a one way transition rate, and if the particle is to be able to move both ways
#  the reverse transition rate must be specified with a similar call with the "via" and "to"
#  arguments reversed.  Alternatively, a shortcut function "setTwoWayTransitionRate" will
#  set up diffusion both directions
sim.setTransitionRate(species='minDatp', via='cytoplasm', to='membrane', rate= 2.5e-12 )

# Populate the model with particles
#  All particles are placed in the cytoplasm
sim.addParticles(species='minDatp', region='cytoplasm', count=1758)
sim.addParticles(species='minDadp', region='cytoplasm', count=1758)
sim.addParticles(species='minE',    region='cytoplasm', count=914)


# Set simulation Parameters
sim.setTimestep(microsecond(50))
sim.setWriteInterval(ms(10))
sim.setLatticeWriteInterval(ms(10))
sim.setSimulationTime(ms(50))

sim.save(args.outputFile)

sim.run(args.outputFile, method='lm::rdme::MpdRdmeSolver', replicates=3)

