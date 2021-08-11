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
#  In order to specify concentrations of particles
#  instead of particle numbers, you must specify
#  the volume of the cell in terms Liters
sim=CME.CMESimulation(1e-15) # E. coli cell size

# Define our chemical species
species = ['A', 'B', 'C']
sim.defineSpecies(species)

# Add reactions to the simulation
sim.addReaction(reactant=('A','B'), product='C', rate=1.07e5)  # /M/s
sim.addReaction(reactant='C', product=('A','B'), rate=3.51e-1) # /s

# Set our initial species counts
#  Here we specify a concentration in units of Molar
sim.addConcentration(species='A', conc=1.66e-6)
sim.addConcentration(species='B', conc=1.66e-6)
sim.addParticles(species='C', count=0)

# Define simulation parameters: run for 10 seconds, saving data every ms
sim.setWriteInterval(ms(1))
sim.setSimulationTime(10)
sim.save(args.outputFile)

# Run 50 replicates using the Gillespie solver
sim.run(filename=args.outputFile, method="lm::cme::GillespieDSolver", replicates=50)

