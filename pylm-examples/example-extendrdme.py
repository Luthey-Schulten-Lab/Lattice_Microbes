# 
# RDME Solver Extension Example
#
# By creating a child class of lm.MpdRdmeSolver in python, we can override the
# virtual method 'hookSimulation' that is called at every frame write time.  We
# will be given a copy of the ByteLattice that we can modify.
#
# As a demonstration, below we create an RDME simulation that has a first order
# reaction of A->B.  We will then at every hook interval add one particle of C.
#

from pyLM import RDME
from pyLM.units import *

### Define our own solver class derived from MpdRdmeSolver
class MyOwnSolver(lm.MpdRdmeSolver):

	# The hookSimulation method defined here will be called at every frame write
	# time.  The return value is either 0 or 1, which will indicate if we
	# changed the state or not and need the lattice to be copied back to the GPU
	# before continuing.  If you do not return 1, your changes will not be
	# reflected.
	def hookSimulation(self, time, lattice):
		print("\n\nBegin hook of MyOwnSolver: Adding 1 particle of type C\n\n")
		lattice.addParticle(16,16,16,3);
		return 1

	def onBeginTrajectory(self):
		print(" *************************************")
		print("          H E L L O ! ! ! ! ! ! !")
		print(" *************************************")
		return 0


	def onEndTrajectory(self):
		print(" *************************************")
		print("    G O O D B Y E ! ! ! ! ! ! ! !")
		print(" *************************************")

		return 0

# End of MyOwnSolver class.  That's it!

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()

latticeSpacing = 32 #nm

# Create our simulation object
sim=RDME.RDMESimulation(dimensions=micron(1.024,1.024,1.024), spacing=nm(latticeSpacing))

# define our chemical species
species = ['A', 'B', 'C']
sim.defineSpecies(species)

# Modify the cytoplasm to add diffusion rates and reactions
sim.modifyRegion('default') \
	.setDefaultDiffusionRate(2.5e-12) \
	.addReaction(reactant='A', product='B', rate=0.5)

sim.addParticles(species='A', region='default', count=1000)
sim.addParticles(species='B', region='default', count=1000)

# Set simulation Parameters
sim.setTimestep(microsecond(50))
sim.setWriteInterval(ms(100))
sim.setLatticeWriteInterval(ms(100))
sim.setSimulationTime(1)
sim.setHookInterval(ms(10))

sim.save(args.outputFile)

# Create an instance of our local solver
solver=MyOwnSolver()
# Call the 'runSolver' method with the supplied solver
sim.runSolver(args.outputFile, solver=solver)

