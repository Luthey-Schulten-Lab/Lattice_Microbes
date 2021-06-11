import os
import random
import time as tm
from pyLM import CME
from pyLM import RDME
from pyLM import LMLogger
from pyLM.units import *
import pySTDLM.PostProcessing as pp
import numpy as np

# Turn on Logging
import logging
LMLogger.setLMLogConsole(logging.INFO)
lmlog=LMLogger.LMLogger

### Define our own solver class derived from MpdHybridSolver
class MyOwnSolver(lm.MpdHybridSolver):
	# Timestep Variables
	tsnum=0
	curt=0.0
	delt=ms(250.0)

	# Species Variables
	numS=10000
	numE=0
	numES=0
	numP=0

	# Time Traces
	times=[]
	traceS=[]
	traceE=[]
	traceES=[]
	traceP=[]

	# The hookSimulation method defined here will be called at every frame 
	# write time.  The return value is either 0 or 1, which will indicate 
	# if we changed the state or not and need the lattice to be copied back 
	# to the GPU before continuing.  If you do not return 1, your changes 
	# will not be reflected.
	def hookSimulation(self, time, lattice):
		print("")
		print("")
		lmlog.info("Hook at time: %f sec"%time)
		lmlog.info("Creating CME simulation...")
		curtime=time
		# Create simulation
		csim=CME.CMESimulation()
		csim.defineSpecies(['E','S','ES','P'])

		# Count enzymes in cell
		parts=lattice.findParticles(4,4)
		Erdme=len(parts)

		# Add reactions and particles
		k1=0.001 # molecules/s
		k2=0.1 # /s
		k3=0.2 # /s
		csim.addReaction(('E','S'),'ES', k1)
		csim.addReaction('ES',('E','S'), k2)
		csim.addReaction('ES',('E','P'), k3)

		csim.addParticles('E',  int(Erdme-self.numES))
		csim.addParticles('ES', int(self.numES))
		csim.addParticles('P',  int(self.numP))
		csim.addParticles('S',  int(self.numS))

		# Set time data
		csim.setWriteInterval(ms(1))
		csim.setSimulationTime(self.delt)

		# Save and run simulation
		filename='cmeSim.%d.lm'%self.tsnum
		lmlog.info("Saving %s..."%filename)
		tm.sleep(1)
		csim.save(filename)
		lmlog.info("Running CME simulation...")
		os.system("lm -r 0-1 -ws -sl lm::cme::GillespieDSolver -f %s"%filename)
		tm.sleep(1)
		self.tsnum += 1

		# Read CME state
		lmlog.info("Postprocessing...")
		fHandle=pp.openLMFile(filename)
		S =pp.getSpecieTrace(fHandle, 'S')
		endidx=len(S)-1
		E =pp.getSpecieTrace(fHandle, 'E')
		ES=pp.getSpecieTrace(fHandle, 'ES')
		P =pp.getSpecieTrace(fHandle, 'P')

		ts = pp.getTimesteps(fHandle)
		tsShifted=[]
		for i in range(len(ts)):
			tsShifted.append(ts[i]+self.curt)
		self.curt = curtime
		self.times.extend(tsShifted)
		pp.closeLMFile(fHandle)

		# Add product to the RDME simulation
		for i in range(P[endidx]-self.numP):
			while True:
				x=random.randint(0,lattice.getXSize()-1)
				y=random.randint(0,lattice.getYSize()-1)
				z=random.randint(0,lattice.getZSize()-1)
			
				if lattice.getOccupancy(x,y,z) < 7:
					# Add a product particle
					lattice.addParticle(x,y,z,5)
					break

		# Update Solver internals
		lmlog.info("Updating internals...")
		self.numS=S[endidx]
		self.numE=E[endidx]
		self.numES=ES[endidx]
		self.numP=P[endidx]
		self.traceS.extend(S)
		self.traceE.extend(E)
		self.traceES.extend(ES)
		self.traceP.extend(P)

		lmlog.info("Resuming RDME simulation...")
		return 1

	def saveTraces(self):
		allTraces=[self.times, self.traceE, self.traceS, self.traceES, self.traceP]
		np.savetxt('T3.2-CMETraces.dat', np.transpose(allTraces))

	def plotTraces(self):
		plotStr="gnuplot plotter.gp"
		os.system(plotStr)
# End of MyOwnSolver class.  That's it!



# Create our simulation object
latticeSpacing = 32 #nm
sim=RDME.RDMESimulation(dimensions=micron(1.024,1.024,1.024), spacing=nm(latticeSpacing))

# define our chemical species
species = ['DNA', 'Polymerase', 'RNA', 'Enzyme', 'Product']
sim.defineSpecies(species)

# Modify the cytoplasm to add diffusion rates and reactions
reg=sim.modifyRegion('default')
reg.setDefaultDiffusionRate(2.5e-12) 
reg.setDiffusionRate('DNA', 0.0)
# 480 nt long with rate of 70 nt/sec
reg.addReaction(reactant=('DNA', 'Polymerase'), product=('DNA', 'Polymerase', 'RNA'), rate=0.5e-3)
# Protein length ~ 160 aa with rate of 40 aa/sec
reg.addReaction('RNA', ('RNA','Enzyme'), 0.25)
# RNA lifetime of 2 minutes
reg.addReaction('RNA','', 1.0/120.0)
# Enzyme is around for about 1/3 the cell cycle
reg.addReaction('Enzyme','', 1.0/1200.0)

# Add particles
sim.addParticles(species='DNA', region='default', count=1)
sim.addParticles('Polymerase', 'default', 460)

# Set simulation Parameters
sim.setTimestep(microsecond(50))
sim.setWriteInterval(ms(250))
sim.setLatticeWriteInterval(ms(250))
sim.setSimulationTime(60)

rdmeFilename="T3.2-mixedRDMECME.lm"
sim.save(rdmeFilename)

# Create an instance of our local solver
solver=MyOwnSolver()
# Call the 'runSolver' method with the supplied solver 
# then perform the post-processing
sim.runSolver(rdmeFilename, solver=solver)

# Post-process the data
solver.saveTraces()
rHandle=pp.openLMFile(rdmeFilename)
times=pp.getTimesteps(rHandle)
R=pp.getSpecieTrace(rHandle,'RNA')
E=pp.getSpecieTrace(rHandle,'Enzyme')
P=pp.getSpecieTrace(rHandle,'Product')
rdmeTraces=[times,R,E,P]
np.savetxt('T3.2-RDMETraces.dat',np.transpose(rdmeTraces))
solver.plotTraces()
pp.closeLMFile(rHandle)

