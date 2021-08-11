from pyLM import *
from pyLM.units import *
from pySTDLM.PostProcessing import *
from pySTDLM.NetworkVisualization import *
import math as m
import numpy as np
import matplotlib.pyplot as plt


# Create our simulation object
sim=RDME.RDMESimulation(dimensions=micron(1.024,1.024,4.096), spacing=nm(32))

# Build a capsid cell
sim.buildCapsidCell(length=micron(4), diameter=micron(1), \
                    membraneThickness=nm(64))

# define our chemical species
sim.defineSpecies(['minDadp', 'minDatp', 'minDm', 'minE', 'minDEm'])

# Modify the cytoplasm to add diffusion rates and reactions
sim.modifyRegion('cytoplasm') \
	.setDefaultDiffusionRate(2.5e-12) \
	.addReaction(reactant='minDadp', product='minDatp', rate=0.5)

# Modify the membrane to add reactions
V=3.14*4.0e-6*0.5e-6*0.5e-6*1000.0 # Liters
N_A=6.022e23                # molecules/mole
scalar=1.0/(N_A*V)

sim.modifyRegion('membrane') \
	.setDefaultDiffusionRate(2.5e-12) \
	.setDiffusionRate(species='minDm',  rate=1e-14) \
	.setDiffusionRate(species='minDEm', rate=1e-14) \
	.addReaction( 'minDatp', 'minDm',   rate=0.78) \
	.addReaction(('minDatp','minDm'), ('minDm','minDm'),rate=9e6*scalar) \
	.addReaction(('minDm','minE'),     'minDEm',     rate=5.56e7*scalar) \
	.addReaction( 'minDEm',           ('minE','minDadp'),rate=0.7) 

# Set diffusive properties between regions
sim.setTransitionRate(species='minDatp', via='cytoplasm', \
                      to='membrane', rate=2.5e-12)
sim.setTransitionRate(species='minDadp', via='cytoplasm', \
                      to='membrane', rate=2.5e-12)
sim.setTransitionRate(species='minE', via='cytoplasm', \
                      to='membrane', rate=2.5e-12)
sim.setTransitionRate(species='minDatp', to='cytoplasm', \
                      via='membrane', rate=2.5e-12)
sim.setTransitionRate(species='minDadp', to='cytoplasm', \
                      via='membrane', rate=2.5e-12)
sim.setTransitionRate(species='minE', to='cytoplasm', \
                      via='membrane', rate=2.5e-12)

# Populate the model with particles
sim.addParticles(species='minDatp', region='cytoplasm', count=1758)
sim.addParticles(species='minDadp', region='cytoplasm', count=1758)
sim.addParticles(species='minE',    region='cytoplasm', count=914)

# Set simulation Parameters
sim.setTimestep(microsecond(50))
sim.setWriteInterval(0.5)
sim.setLatticeWriteInterval(0.5)
sim.setSimulationTime(600)

# Run simulation
filename="T3.1-MinSystem.v3.lm"
sim.save(filename)
sim.run(filename, method='lm::rdme::MpdRdmeSolver', replicates=1)


###################
# Post-Processing #
###################
# Plot Kymogrpahs for MinE and MinD
plotOccupancyKymograph(filename, species='minDm', \
                                      filename='MinDmKymo.png')
plotOccupancyKymograph(filename, species='minE', \
                                      filename='MinEKymo.png')

# Calculate the Oscillation Frequency


