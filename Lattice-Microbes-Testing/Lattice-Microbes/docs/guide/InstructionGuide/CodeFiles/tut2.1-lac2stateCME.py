from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.PostProcessing import *
from pySTDLM.NetworkVisualization import *
import math as m
import numpy as np
import matplotlib.pyplot as plt

# Create our simulation object
sim=CME.CMESimulation()

# Add the reactants
species = ['R2','O','R2O','IR2','IR2O','I2R2','I2R2O','mY','Y','I','Iex','YI']
sim.defineSpecies(species)

scalar = 2.076e-9 # Rate conversion from experiments to stochastic

# Add the reactions
# Lac operon regulation
sim.addReaction(reactant=('R2','O'), product='R2O',rate=2.43e6*scalar)
sim.addReaction(('IR2','O'),  'IR2O',       1.21e6*scalar)
sim.addReaction(('I2R2','O'), 'I2R2O',      2.43e4*scalar)
sim.addReaction('R2O',        ('R2','O'),   6.30e-4)
sim.addReaction('IR2O',       ('IR2','O'),  6.30e-4)
sim.addReaction('I2R2O',      ('I2R2','O'), 3.15e-1)

# Transcription, translation, and degredation
sim.addReaction('O',          ('O','mY'),   1.26e-1)
sim.addReaction('mY',         ('mY','Y'),   4.44e-2)
sim.addReaction('mY',         '',           1.11e-2)
sim.addReaction('Y',          '',           2.10e-4)

# Inducer-repressor interactions
sim.addReaction(('I','R2'),   'IR2',        2.27e4*scalar)
sim.addReaction(('I','IR2'),  'I2R2',       1.14e4*scalar)
sim.addReaction(('I','R2O'),  'IR2O',       6.67e2*scalar)
sim.addReaction(('I','IR2O'), 'I2R2O',      3.33e2*scalar)
sim.addReaction('IR2',        ('I','R2'),   2.00e-1)
sim.addReaction('I2R2',       ('I','IR2'),  4.00e-1)
sim.addReaction('IR2O',       ('I','R2O'),  1.00)
sim.addReaction('I2R2O',      ('I','IR2O'), 2.00)

# Inducer transport
sim.addReaction('Iex',        'I',          2.33e-3)
sim.addReaction('I',          'Iex',        2.33e-3)
sim.addReaction(('Y','Iex'),  ('YI','Iex'), 3.03e4*scalar)
sim.addReaction('YI',         ('Y','Iex'),  1.20e-1)
sim.addReaction('YI',         ('Y','I'),    1.20e+1)

# Populate the model with particles
sim.addParticles(species='R2',   count=9)
sim.addParticles('O',    1)
sim.addParticles('Y',    30)
sim.addParticles('I',    7224*4)
sim.addParticles('Iex',  7224*4)

# Set up the times
sim.setTimestep(ms(100))
sim.setWriteInterval(1)
sim.setSimulationTime(3600.0)

# Save the simulation state to a file
filename='T2.1-lac2state.lm'
sim.save(filename)

plotCMEReactionNetwork(sim, filename="lacGraph.gml")

# Run the simulation for 10 hours
reps=100
sim.run(filename, "lm::cme::GillespieDSolver", reps)

# Post-processing
for i in range(1,reps):
	plotTraceFromFile(filename, ['mY','Y'], i, "Trace.mY.Y.%d.png"%(i))

##########################
# Custom Post-Processing #
##########################
# Get out the data
fileHandle=openLMFile(filename)
times = getTimesteps(fileHandle)
simmY = reps*[]
simY  = reps*[]
for i in range(1,reps):
	simY.append(getSpecieTrace(fileHandle, 'Y', i))
	simmY.append(getSpecieTrace(fileHandle, 'mY', i))

# Find min and max of mY and Y
minY  = 1e9; maxY  = 0; minmY = 1e9; maxmY = 0
for i in range(reps-1):
	for y in simY[i]:
		minY = min(minY, y)
		maxY = max(maxY, y)
	for my in simmY[i]:
		minmY = min(minmY, my)
		maxmY = max(maxmY, my)

# Histogram the data (the hard way)
hist = np.zeros(shape=(maxmY,maxY))
for i in range(reps-1):
	for j in range(len(times)):
		hist[simmY[i][j]-1, simY[i][j]-1] += 1.0
# Compute logarithm of histogram value
histLog = np.zeros(shape=(maxmY,maxY))
for i in range(maxmY-1):
	for j in range(maxY-1):
		if hist[i,j] > 0:
			histLog[i,j] = m.log10(hist[i,j])

# Plot ourselves a histogram
plt.clf()
plt.imshow(histLog.transpose(), interpolation='nearest', origin='lower', aspect='auto')
plt.colorbar()
plt.xlabel('mRNA Count')
plt.ylabel('LacY Count')
plt.savefig('LacYmRNAHeatmap.png')

# Clean up after ourselves
closeLMFile(fileHandle)


