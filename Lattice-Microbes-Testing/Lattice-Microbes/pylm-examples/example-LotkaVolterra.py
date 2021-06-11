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

# Define our chemical species
species = ['P1', 'P2']
sim.defineSpecies(species)

# Add reactions to the simulation
sim.addReaction(reactant='P1',        product=('P1','P1'), rate=1)
sim.addReaction(reactant=('P1','P2'), product=('P2','P2'), rate=0.005)
sim.addReaction(reactant='P2',        product='',          rate=0.6)

# Set our initial species counts
sim.addParticles(species='P1', count=50)
sim.addParticles(species='P2', count=100)

# Define simulation parameters: run for 10 seconds, saving data every ms
sim.setWriteInterval(ms(1))
sim.setSimulationTime(100)
sim.save(args.outputFile)

# Run 50 replicates using the Gillespie solver
sim.run(filename=args.outputFile, method="lm::cme::GillespieDSolver", replicates=50)

# Plot average and variance
PostProcessing.plotAvgVarFromFile(args.outputFile, ['P1','P2'], 'LotkaVolterraDynamics.png')
PostProcessing.plotTraceFromFile(args.outputFile, ['P1','P2'], 1, 'LV.1.png')
PostProcessing.plotTraceFromFile(args.outputFile, ['P1','P2'], 2, 'LV.2.png')

