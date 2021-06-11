import sys, os, time
from pyLM.units import *
import pyLM.CME as cme
import pySTDLM.PostProcessing as pp
import pySTDLM.StandardReactionSystems as srs
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt

# Tests various reaction orders and compares
#  them to the exact solution for 1000 replicates.
#  Curently, it tests for 1% agreement.
#   0 -> bool, [error strings]
def testReactionOrders():
	didPass = True
	errors = []

	# Create simulation
	sim = cme.CMESimulation()
	sim.defineSpecies(["A","B","C","D","E"])

	# Zeroth order
	k0order=10.0
	sim.addReaction("","A",  k0order)
	# First order
	k1order=0.3
	sim.addReaction("B","C", k1order)
	# Second order
	k2order=0.05
	sim.addReaction(("D","D"),"E", k2order)

	# Add particles
	sim.addParticles("B", 1000)
	sim.addParticles("D", 500)

	# Run simulation
	sim.setWriteInterval(ms(1))
	sim.setSimulationTime(10)
	fileN = "testsuites/CMEBase/testReactionOrder.lm"
	sim.save(fileN)
	sim.run(fileN, method="lm::cme::GillespieDSolver", replicates=1000)

	# Open Data
	f = pp.openLMFile(fileN)

	# Fitting functions
	def zerothOrder(t,k):
		return k*t
	def firstOrder(t,k,C0):
		return C0*(1.0-np.exp(-k*t))
	def firstDegOrder(t,k,C0):
		return C0*np.exp(k*t)
	def secondOrderSelf(t,k,C0):
		return 1.0/(k*t + 1.0/C0)

	# Fit curves
	Aavg, Astd, times = pp.getAvgVarTrace(f, "A")
	popt,pcov = spo.curve_fit(zerothOrder, times,Aavg,1.0)
	if popt[0] > 1.01*k0order or popt[0] < 0.99*k0order:
		didPass = False
		errors.append("Zeroth order reaction was not within 1%%. Got %f, expected %f"%(popt[0], k0order))
		plt.semilogx(times,Aavg, label="Simulation")
		plt.semilogx(times,zerothOrder(np.array(times),k0order), label="Analytic")
		plt.legend()
		plt.savefig("testsuites/CMEBase/ZerothOrderDiagnostic.png")
		plt.clf()

	Bavg, Bstd, times = pp.getAvgVarTrace(f, "B")
	Cavg, Cstd, times = pp.getAvgVarTrace(f, "C")
	popt,pcov = spo.curve_fit(firstOrder, times,Cavg,[1.0,1000])
	print(popt)
	if popt[0] > 1.01*k1order or popt[0] < 0.99*k1order:
		didPass = False
		errors.append("First order production reaction was not within 1%%. Got %f, expected %f"%(popt[0], k1order))
		plt.semilogx(times,Cavg, label="Simulation")
		plt.semilogx(times,firstOrder(np.array(times),k1order,1000), label="Analytic")
		plt.legend()
		plt.savefig("testsuites/CMEBase/FirstOrderDiagnostic.png")
		plt.clf()
	popt,pcov = spo.curve_fit(firstDegOrder, times,Bavg,[-1.0,1000])
	print(popt)
	if popt[0] < -1.01*k1order or popt[0] > -0.99*k1order:
		didPass = False
		errors.append("First order degradation reaction was not within 1%%. Got %f, expected %f"%(popt[0], -k1order))
		plt.semilogx(times,Bavg, label="Simulation")
		plt.semilogx(times,firstDegOrder(np.array(times),-k1order,1000), label="Analytic")
		plt.legend()
		plt.savefig("testsuites/CMEBase/FirstOrderDegDiagnostic.png")
		plt.clf()
	
	Davg, Dstd, times = pp.getAvgVarTrace(f, "D")
	popt,pcov = spo.curve_fit(secondOrderSelf, times,Davg,[0.5,450])
	print(popt)
	if popt[0] > 1.01*k2order or popt[0] < 0.99*k2order:
		didPass = False
		errors.append("Second order degradation reaction was not within 1%%. Got %f, expected %f"%(popt[0], k2order))
		plt.semilogx(times,Davg, label="Simulation")
		plt.semilogx(times,secondOrderSelf(np.array(times),k2order,500), label="Analytic")
		plt.legend()
		plt.savefig("testsuites/CMEBase/SecondOrderSelfDiagnostic.png")
		plt.clf()

	# Clean up
	pp.closeLMFile(f)
	os.remove(fileN)

	return didPass, errors

# Tests for correct noise statistics 
#  in the prototypical example of 
#  gene expression from DNA.  Compares
#  the mean and variance with the Poisson
#  and Gamma distribution results.  Only
#  requires 10% accuracy due to the time
#  restraints.
#   0 -> bool, [error strings]
def testNoise():
	didPass = True
	errors = []

	# Create simulation
	sim = cme.CMESimulation()
	sim.defineSpecies(["D","M","P"])
	

	# Add constitutive expression
	k_trn = 0.02
	d_m   = 0.01
	k_tsl = 0.3
	d_p   = 0.0004
	sim.addReaction("D",("D","M"),  k_trn)
	sim.addReaction("M",("M","P"),  k_tsl)
	sim.addReaction("M","",  d_m)
	sim.addReaction("P","",  d_p)
	
	# Add Particles
	sim.addParticles("D",1)

	# Run simulation
	sim.setWriteInterval(1.0)
	sim.setSimulationTime(10000.0)
	fileN = "testsuites/CMEBase/testNoise.lm"
	sim.save(fileN)
	sim.run(fileN, method="lm::cme::GillespieDSolver", replicates=5000)

	# Analytic results
	avgM_ana = k_trn/d_m
	varM_ana = avgM_ana
	avgP_ana = k_tsl/d_p*avgM_ana
	varP_ana = k_trn/d_p*(k_tsl/d_m)**2
	
	# Read data and check results
	f = pp.openLMFile(fileN)
	avgM, varM, times = pp.getAvgVarTrace(f, "M")
	avgP, varP, times = pp.getAvgVarTrace(f, "P")
	avgM = np.average(avgM[len(avgM)/2.0:])
	varM = np.average(varM[len(varM)/2.0:])
	avgP = np.average(avgP[len(avgP)/2.0:])
	varP = np.average(varP[len(varP)/2.0:])

	if avgM > 1.1*avgM_ana or avgM < 0.9*avgM_ana:
		didPass = False
		errors.append("Average mRNA count wrong.  Got %f, expected %f"%(avgM, avgM_ana))
	if varM > 1.1*varM_ana or varM < 0.9*varM_ana:
		didPass = False
		errors.append("Variance in mRNA count wrong.  Got %f, expected %f"%(varM, varM_ana))
	if avgP > 1.1*avgP_ana or avgP < 0.9*avgP_ana:
		didPass = False
		errors.append("Average Protein count wrong.  Got %f, expected %f"%(avgP, avgP_ana))
	if varP > 1.1*varP_ana or varP < 0.9*varP_ana:
		didPass = False
		errors.append("Variance in Protein count wrong.  Got %f, expected %f"%(varP, varP_ana))

	# Clean up
	pp.closeLMFile(f)
	os.remove(fileN)

	return didPass, errors

# Tests that the exact same solution is generated by the new 
#  code compared to a gold standard file.
#   0 -> bool, [error strings]
def testSpeciesCounts():
	didPass = True
	errors = []

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
	fileN = "testsuites/CMEBase/lacSeedTest.lm"
	sim.save(fileN)
	sim.run(fileN, method="lm::cme::GillespieDSolver", replicates=1, seed=123456789)
	
	# Open old file
	f_old = pp.openLMFile("testsuites/CMEBase/lacSeedTestGoldStandard.lm")
	f_new = pp.openLMFile(fileN)

	# Compare time trace
	timesOld = pp.getTimesteps(f_old)
	timesNew = pp.getTimesteps(f_new)
	for i,t in enumerate(timesOld):
		if t != timesNew[i]:
			didPass = False
			errors.append("Failed timestamp comparisons at: old=%f, new=%f"%(t, timeNew[i]))
			break

	# Compare each species check for 1-to-1 comparison of time trace
	for s in sim.species_id:
		s_old = pp.getSpecieTrace(f_old, s)
		s_new = pp.getSpecieTrace(f_new, s)
		for i,t in enumerate(timesOld):
			if s_old[i] != s_new[i]:
				didPass = False
				errors.append("Failed specie comparisons for s='%s' at: old=%f, new=%f; given: old=%d, new=%d"%(s, t, timeNew[i], s_old[i], s_new[i]))
				break
			
	# Clean up
	pp.closeLMFile(f_old)
	pp.closeLMFile(f_new)
	if didPass:
		os.remove(fileN)

	return didPass, errors


# This function is imported
#  0 -> int : number passed, 
#       int : total tests,
#       []  : times for each test
#       {str : testname -> (str: Passed/Failed, [] : a list of error strings/exceptions)} : Test results 
def LMRegressionTester():
	# Return variables
	passedTests = 0
	testResults = {}
	timings = []
	
	# Tests array
	tests = []
	
	# Add all the tests in this file
	tests.append(("RxnOrder",testReactionOrders))
	tests.append(("Noise",testNoise))
	tests.append(("Counts",testSpeciesCounts))
	
	# Run tests
	for testName, test in tests:
		startT = time.time()
		try:
			res, reason = test()
			if res:
				testResults[testName] = ("Passed",[])
				passedTests += 1
			else:
				testResults[testName] = ("Failed",reason)
		except Exception as err:
			testResults[testName] = ("Failed",[err.message])
		endT = time.time()
		timings.append(endT - startT)

	return passedTests, len(tests), timings, testResults
	

