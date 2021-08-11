from pyLM import *
from pyLM.units import *
from pySTDLM import *

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()


# Set up logging
import logging
LMLogger.setLMLoggerLevel(logging.INFO)

# Create our CME simulation object
sim=CME.CMESimulation()

# define our chemical species
species = ['DNA', 'RNA', 'Protein']
sim.defineSpecies(species)

# Add reactions to the simulation
sim.addReaction(reactant='DNA', product=('DNA','RNA'), rate=1.0/120.0)
sim.addReaction(reactant='RNA', product=('RNA','Protein'), rate=1.0/60.0)
sim.addReaction(reactant='RNA', product='', rate=1.0/120.0)
sim.addReaction(reactant='Protein', product='', rate=1.0/(15.0*60.0))

# Set our initial species counts
sim.addParticles(species='DNA', count=1)
sim.addParticles(species='RNA', count=2)
sim.addParticles(species='Protein', count=20)

# Define simulation parameters: run for 10 seconds, saving data every ms
sim.setWriteInterval(1.0)
sim.setSimulationTime(36000.0)
sim.save(args.outputFile)

# Run 100 replicates using the Gillespie solver
sim.run(filename=args.outputFile, method="lm::cme::GillespieDSolver", replicates=100)


# There is only one new function called in this example, "showAvgVarFromFile"
#  which is basically equivalent to "plotAvgVarFromFile" but that shows the
#  plot directly, instead of plotting it to a file.
PostProcessing.showAvgVarFromFile(args.outputFile, ['DNA','RNA','Protein']) 
PostProcessing.plotAvgVarFromFile(args.outputFile, ['DNA','RNA','Protein'], 'RNAProteinAvgVar.png')

