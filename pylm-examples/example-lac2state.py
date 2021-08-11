from pyLM import *
from pyLM.units import *
from pySTDLM import *

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()


# Create our simulation object
sim=CME.CMESimulation()

# Add the reactants
species = ['R2','O','R2O','IR2','IR2O','I2R2','I2R2O','mY','Y','I','Iex','YI']
sim.defineSpecies(species)

# Add the reactions
# Lac operon regulation
scalar = 2.076e-9
sim.addReaction(reactant=('R2','O'),   product='R2O',        rate=2.43e6*scalar)
sim.addReaction(reactant=('IR2','O'),  product='IR2O',       rate=1.21e6*scalar)
sim.addReaction(reactant=('I2R2','O'), product='I2R2O',      rate=2.43e4*scalar)
sim.addReaction(reactant='R2O',        product=('R2','O'),   rate=6.30e-4)
sim.addReaction(reactant='IR2O',       product=('IR2','O'),  rate=6.30e-4)
sim.addReaction(reactant='I2R2O',      product=('I2R2','O'), rate=3.15e-1)

# Transcription, translation, and degredation
sim.addReaction(reactant='O',          product=('O','mY'),   rate=1.26e-1)
sim.addReaction(reactant='mY',         product=('mY','Y'),   rate=4.44e-2)
sim.addReaction(reactant='mY',         product='',           rate=1.11e-2)
sim.addReaction(reactant='Y',          product='',           rate=2.10e-4)

# Inducer-repressor interactions
sim.addReaction(reactant=('I','R2'),   product='IR2',        rate=2.27e4*scalar)
sim.addReaction(reactant=('I','IR2'),  product='I2R2',       rate=1.14e4*scalar)
sim.addReaction(reactant=('I','R2O'),  product='IR2O',       rate=6.67e2*scalar)
sim.addReaction(reactant=('I','IR2O'), product='I2R2O',      rate=3.33e2*scalar)
sim.addReaction(reactant='IR2',        product=('I','R2'),   rate=2.00e-1)
sim.addReaction(reactant='I2R2',       product=('I','IR2'),  rate=4.00e-1)
sim.addReaction(reactant='IR2O',       product=('I','R2O'),  rate=1.00)
sim.addReaction(reactant='I2R2O',      product=('I','IR2O'), rate=2.00)

# Inducer transport
sim.addReaction(reactant='Iex',        product='I',          rate=2.33e-3)
sim.addReaction(reactant='I',          product='Iex',        rate=2.33e-3)
sim.addReaction(reactant=('Y','Iex'),  product=('YI','Iex'), rate=3.03e4*scalar)
sim.addReaction(reactant='YI',         product=('Y','Iex'),  rate=1.20e-1)
sim.addReaction(reactant='YI',         product=('Y','I'),    rate=1.20e+1)


# Populate the model with particles
sim.addParticles(species='R2',   count=9)
sim.addParticles(species='O',    count=1)
sim.addParticles(species='Y',    count=30)
sim.addParticles(species='I',    count=7224)
sim.addParticles(species='Iex',  count=7224)

# Set up the times
sim.setTimestep(ms(1))
sim.setWriteInterval(1)
sim.setSimulationTime(3600.0)


# Save the simulation state to a file
sim.save(args.outputFile)

# Run the simulation for 10 hours
#  This time we use the Next Reaction Solver instead of the Gillespie solver to speed
#  up the time to solution.
sim.run(filename=args.outputFile, method="lm::cme::NextReactionSolver", replicates=10)

