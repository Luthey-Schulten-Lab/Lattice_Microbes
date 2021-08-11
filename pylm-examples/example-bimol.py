# Import Libraries #
# In order to use pyLM we need to import several libraries.  
#  The first is pyLM proper.  The second is a library with
#  a number of functions such as nm(), micron(), ms(),
#  microsecond(), etc. that allow cleaner definition of
#  units.  Finally, we import the pyLM standard library
#  of functionality "pySTDLM", which contains plotting
#  and post-processing commands.
from pyLM import *
from pyLM.units import *
from pySTDLM import *

# Setup Command Line #
# The next several lines set up the command line
#  parameters required to run this file.
#  Specifically, we require that an output file name
#  be specified.
import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()


# Set up logging #
#  Here we define the amount of logging information that
#  pyLM will print when it is running.  The amount of
#  information will increase with the levels:
#    DEBUG > INFO > WARNING > ERROR > CRITICAL
import logging
LMLogger.setLMLoggerLevel(logging.INFO)

# Specify Simulation #
# We begin by creating a CME Simulation that will 
#  include the whole problem definition.
sim=CME.CMESimulation()

# Here we register the chemical species with simulation;
#  first by specifing the names of the chemical species
#  followed by the actual command to register those species
#  with the simulation.  The defineSpecies command can
#  be called multiple times and will add any new names
#  to the list of species.
species = ['A', 'B', 'C']
sim.defineSpecies(species)

# Add reactions to the simulation adding a bimolecular 
#  association reaction and a unimolecular dissociation
#  reaction. When more than one reactant is involved,
#  the list of reactant names should be passed as a 
#  tuple.
sim.addReaction(reactant=('A','B'), product='C', rate=1.78e-4)
sim.addReaction(reactant='C', product=('A','B'), rate=3.51e-1)

# Set our initial species counts
sim.addParticles(species='A', count=1000)
sim.addParticles(species='B', count=1000)
sim.addParticles(species='C', count=0)

# Define simulation parameters: run for 10 seconds, 
#  saving data every ms. The simulation must be saved to 
#  disk using the save() command before the simulation can be run.
sim.setWriteInterval(ms(1))
sim.setSimulationTime(10)
sim.save(args.outputFile)

# Run the Simulation #
# Here we declare that the simulation will be run using the Gillespie 
#  algorithm and that we want 50 independent samplings of the
#  solution.
sim.run(filename=args.outputFile, method="lm::cme::GillespieDSolver", replicates=50)

