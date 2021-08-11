from pyLM import *
from pyLM.units import *

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-i', '--inputFile', default=None)
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()

import pySTDLM.StandardCells as sc

latticeSpacing = 16 #nm

# Create our simulation object
sim=RDME.RDMESimulation(dimensions=micron(1.024,1.024,4.096), spacing=nm(latticeSpacing))

# Build a capsid cell
sim.buildCapsidCell(length=micron(4), diameter=micron(1), membraneThickness=nm(32))

sc.packFastGrowingEcoli(sim)

# define our chemical species
species = ['minDatp', 'minDadp', 'minE', 'minDm', 'minDEm']
sim.defineSpecies(species)

# Modify the cytoplasm to add diffusion rates and reactions
sim.modifyRegion('cytoplasm') \
	.setDefaultDiffusionRate(2.5e-12) \
	.addReaction(reactant='minDadp', product='minDatp', rate=0.5)

# Modify the membrane to add reactions
sim.modifyRegion('membrane') \
	.setDefaultDiffusionRate(2.5e-12) \
	.setDiffusionRate(species='minDm',  rate=1e-14) \
	.setDiffusionRate(species='minDEm', rate=1e-14) \
	.addReaction(reactant='minDatp', product='minDm', rate=1e-12) \
	.addReaction(reactant=('minDatp','minDm'), product=('minDm','minDm'), rate=9e6) \
	.addReaction(reactant=('minDm','minE'), product='minDEm', rate=5.56e7) \
	.addReaction(reactant='minDEm', product=('minE','minDadp'), rate=0.7) 

# Set diffusive properties between regions
sim.setTransitionRate(species='minDatp', via='cytoplasm', to='membrane', rate= 2.5e-12 )


# Restart if necessary
#  If the user does specify an input file name, then the simulation state
#  from that last simulation will be copied to the beginning state of the
#  new simulation
if args.inputFile != None:
	# Open previous simulation
	oldSimFile = lm.SimulationFile(args.inputFile)

	# Read lattice from old simulation
	oldReplicate = 1
	latticeTimes = oldSimFile.getLatticeTimes(oldReplicate) # Replicate 1 is the default
	latticeDim = int(micron(1.024)/nm(16))
	latticeDimZ = int(micron(4.096)/nm(16))
	refLattice = lm.ByteLattice(latticeDim,latticeDim,latticeDimZ,nm(16),8)
	oldSimFile.getLattice(oldReplicate, len(latticeTimes)-1, refLattice)

	# Get old simulation time and find out how much additional runtime is needed
	#  for 2 seconds total time
	tend = (2.0-latticeTimes[-1])

	# Set simulation time
	sim.setSimulationTime(tend)
	# Copy lattice from final timestep of previous simulation
	newLattice = sim.getLattice()
	for z in range(latticeDimZ):
		for y in range(latticeDim):
			for x in range(latticeDim):
				# Copy the lattice site type
				#  Do this step only if you haven't set up the simulation as shown above
#				newLattice.setSiteType(x,y,z, refLattice.getSiteType(x,y,z))
				# Copy the particles fromt he old to new site
				for i in range(refLattice.getOccupancy(x,y,z)):
					newLattice.addParticle(x,y,z, refLattice.getParticle(x,y,z,i))
else:
	# Populate the model with particles
	sim.addParticles(species='minDatp', region='cytoplasm', count=1758)
	sim.addParticles(species='minDadp', region='cytoplasm', count=1758)
	sim.addParticles(species='minE',    region='cytoplasm', count=914)
	sim.setSimulationTime(2)



# Set simulation Parameters
sim.setTimestep(microsecond(50))
sim.setWriteInterval(ms(100))
sim.setLatticeWriteInterval(ms(100))

sim.save(args.outputFile)

sim.run(args.outputFile, method='lm::rdme::MpdRdmeSolver', replicates=1)

