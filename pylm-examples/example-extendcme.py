# 
# CME Extension Example
#
# By creating a child class of the solver in python, we can override the
# virtual method 'hookSimulation' that is called at every frame write time.  We
# will be given a copy of the ByteLattice that we can modify.
#
# As a demonstration, below we create a CME simulation that has a first order
# reaction of A->B.  We will then at every hook interval add one particle of C.
#

from pyLM import CME
from pyLM.units import *

### Define our own solver class derived from the LM solver:
class MyOwnSolver(lm.GillespieDSolver):

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


	# The hookSimulation method defined here will be called at every frame write
	# time.  The return value is either 0 or 1, which will indicate if we
	# changed the state or not and need the lattice to be copied back to the GPU
	# before continuing.  If you do not return 1, your changes will not be
	# reflected.
	def hookSimulation(self, time):
		print("In Hook Simulation:",time)
		speciesCounts = self.getSpeciesCountView()
		print("A: %d\tB: %d\tC: %d"%(speciesCounts[0],speciesCounts[1],speciesCounts[2]))
		if time > 0.5:
			speciesCounts[2] += 1
		return 1
# End of MyOwnSolver class.  That's it!

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()

# Create our simulation object
sim=CME.CMESimulation()

# define our chemical species
species = ['A', 'B', 'C']
sim.defineSpecies(species)

# Modify the cytoplasm to add diffusion rates and reactions
sim.addReaction(reactant='A', product='B', rate=0.5)

sim.addParticles(species='A', count=1000)
sim.addParticles(species='B', count=1000)

# Set simulation Parameters
sim.setWriteInterval(ms(100))
sim.setHookInterval(ms(250))
sim.setSimulationTime(1)

sim.save(args.outputFile)
lm.setVerbosityLevel(10)

# Create an instance of our local solver
solver=MyOwnSolver()
# Call the 'runSolver' method with the supplied solver
sim.runSolver(args.outputFile, solver=solver)



