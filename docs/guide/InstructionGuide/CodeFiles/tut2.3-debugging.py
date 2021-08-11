import numpy as np
import matplotlib.pyplot as plt
from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.PostProcessing import *

# Turn on logger
from pyLM import LMLogger
LMLogger.setLMLogConsole(logging.WARNING)

# Constants
V  = 1.0e-15       # L
NA = 6.022e23 	   # molecules/mole
kf = 1.07e5/(NA*V) # /M/s
kr = 0.351         # /s

# Create our CME simulation object
sim=CME.CMESimulation()

# define our chemical species
species = ['A', 'B', 'C']
sim.defineSpecies(species)

# Add reactions to the simulation
sim.addReaction(reactant=('A','B'), product='C', rate=kf)
sim.addReaction(reactant='C', product=('A','B'), rate=kr)

# Set our initial species counts
sim.addParticles(species='a', count=1000)
sim.addParticles(species='B', count=1000)
sim.addParticles(species='C', count=0)

# Define simulation parameters: run for 10 seconds, saving data every ms
sim.setWriteInterval(microsecond(30))
sim.setSimulationTime(30)
sim.save('T2.3-bimol.lm')

# Run 1 replicates using the Gillespie solver
sim.run(filename='T2.3-bimol.lm', method="lm::cme::GillespieDSolver", replicates=1)

# Plot the solution
plotTraceFromFile(filename='T2.3-bimol.lm', species=['A','C'], replicate=1, outfile='BimolecularStoch.png')
