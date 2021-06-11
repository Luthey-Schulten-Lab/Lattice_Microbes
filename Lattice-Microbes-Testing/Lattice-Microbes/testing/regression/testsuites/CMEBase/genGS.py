from pyLM.units import *
import pyLM.CME as cme
import pySTDLM.PostProcessing as pp
import pySTDLM.StandardReactionSystems as srs

# Create simulation
sim = cme.CMESimulation()
srs.addLacTwoStateSystem(sim)

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

# Run a single replicate
fileN = "lacSeedTestGoldStandard.lm"
sim.save(fileN)
sim.run(fileN, method="lm::cme::GillespieDSolver", replicates=1, seed=123456789)

