from pyLM import RDME
from pyLM.units import *
from pySTDLM.CellArranger import *


# Create an RDME simulation of the correct size
latticeSpacing = 16 #nm
sim = RDME.RDMESimulation(dimensions=micron(1.024,4.096,4.096), spacing=nm(latticeSpacing))


# Define a CellArranger (Here we use a TightGrid. Other types are available)
tga = TightGridArranger(cellType="Sphere", cellAttributes={"radius":nm(512), "membraneThickness":nm(32.0)}, arrangerAttributes=None)

# Here we pack the cells into the volume (optionally, you can specify a subvolume into which to pack; this defaults to the whole volume)
tga.packVolume(sim)

# Here the cells are discretized to the lattice representing the RDMESimulation
tga.addToSimulation(sim)

# Here additional simulation parameters can be specified (particles, reactions, etc.)
species = ["A","B","C"]
sim.defineSpecies(species)

cyt = sim.modifyRegion("cytoplasm")
mem = sim.modifyRegion("membrane")

cyt.setDefaultDiffusionRate(2.5e-12)
mem.setDefaultDiffusionRate(2.5e-12)
sim.setTransitionRate(species='A', via='cytoplasm', to='membrane', rate= 2.5e-12 )
sim.setTransitionRate(species='B', via='cytoplasm', to='membrane', rate= 2.5e-12 )

cyt.addReaction(reactant=('A','B'), product='C', rate=0.5)
cyt.addReaction(reactant='C', product=('A','B'), rate=0.5)

sim.addParticles(species='A', region='cytoplasm', count=1000)
sim.addParticles(species='B', region='cytoplasm', count=1000)


# Save the simulation and simulate
sim.setTimestep(microsecond(50))
sim.setWriteInterval(ms(100))
sim.setLatticeWriteInterval(ms(100))
sim.setSimulationTime(1.0)

filename="TestArranger.lm"
sim.save(filename)

sim.run(filename, method='lm::rdme::MpdRdmeSolver', replicates=1)




