from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.PostProcessing import *
from pySTDLM.StandardReactionSystems import *
from pySTDLM.StandardCells import *

# Create our simulation object
latticeSpacing=nm(16)
sim=RDME.RDMESimulation(dimensions=micron(1.024,1.024,2.048), \
                        spacing=latticeSpacing)

# Set up cell geometry
sim.buildCapsidCell(length=micron(1.8), diameter=micron(0.5), \
                    membraneThickness=nm(32))

# Add the reactants and reactions
addLacTwoStateSystem(sim)

# Get handles to the different regions
mem='membrane'
cyt='cytoplasm'
ext='default'

# Crowd the cytoplasm
packFastGrowingEcoli(sim)

# Populate the model with particles
sim.addParticles(species='R2', region=cyt,  count=9)
sim.addParticles('O',   cyt, 1)
sim.addParticles('Y',   mem, 30)
sim.addParticles('I',   cyt, 7224*2)
sim.addParticles('Iex', ext, 7224*2)

# Set up the times
sim.setTimestep(microsecond(50))
sim.setWriteInterval(1)
sim.setLatticeWriteInterval(1)
sim.setSimulationTime(7200)

# Save the simulation state to a file
filename='T2.2-lac2state.lm'
sim.save(filename)
# Run the simulation for 10 hours
reps=10
sim.run(filename, "lm::rdme::MpdRdmeSolver", reps)

# Post-processing
for i in range(1,reps+1):
	print("plotting")
	plotTraceFromFile(filename, ['mY','Y'], i, "Trace.mY.Y.%d.png"%(i))

