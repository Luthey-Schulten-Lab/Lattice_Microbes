# 
# University of Illinois Open Source License
# Copyright 2008-2018 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Michael J. Hallock and Joseph R. Peterson
# 
#

from pyLM import *
from pyLM.CME import *
from pyLM.RDME import *


# ######################
# E. coli MinDE System #
# ######################
def addMinDESystem(sim):
	"""Adds the MinDE system in Ecoli as a standard reaction network.

    Args:
        sim:
            A RDMESimulation object with the "cytoplasm" and "membrane" regions defined
    Returns:
        A pointer to the simulation object that was passed in
    """
	if not isinstance(sim, RDMESimulation):
		LMLogger.error('addMinDESystem Failed')
		raise Exception("ERROR: In addMinDESystem, argument is not RDMESimulation.")

	species = ['minDatp', 'minDadp', 'minE', 'minDm', 'minDEm']
	sim.defineSpecies(species)

	# Modify the cytoplasm to add diffusion rates and reactions
	sim.modifyRegion('cytoplasm') \
		.setDefaultDiffusionRate(2.5e-12) \
		.addReaction(reactant='minDadp', product='minDatp', rate=0.5)

	# Modify the membrane to add reactions
	sim.modifyRegion('membrane') \
		.setDefaultDiffusionRate(2.5e-12) \
		.setDiffusionRate(species='minDm',  rate= 0.0125e-6 / nm(sim.latticeSpacing)) \
		.setDiffusionRate(species='minDEm', rate=1e-14) \
		.setDiffusionRate(species='minE', rate=1e-14) \
		.addReaction(reactant='minDatp', product='minDm', rate=1e-12) \
		.addReaction(reactant=('minDatp','minDm'), product=('minDm','minDm'), rate=9e6) \
		.addReaction(reactant=('minDm','minE'), product='minDEm', rate=5.56e7) \
		.addReaction(reactant='minDEm', product=('minE','minDadp'), rate=0.7)

	# Set diffusive properties between regions
	sim.setTransitionRate(species='minDatp', via='cytoplasm', to='membrane', rate=2.5e-12 )

	# Populate the model with particles
	# TODO: Compute the volume of the cell and modify counts
	#sim.addParticles(species='minDatp', region='cytoplasm', count=1758)
	#sim.addParticles(species='minDadp', region='cytoplasm', count=1758)
	#sim.addParticles(species='minE',    region='cytoplasm', count=914)

	return sim




# #####################
# E. coli Lac Systems #
# #####################
def addLacTwoStateSystem(sim, inducerType = 'TMG'):
	"""Lac switch reaction with two states of the DNA.

    Reference: E. Roberts, A. Magis, J.O. Ortiz, W. Baumeister, Z. Luthey-Schulten.
      Noise Contributions in an Inducible Genetic Switch: A Whole-Cell
      Simulation Study. PLoS Comput. Bio. 7(3): 2011, e1002010.

    Args:
        sim:
            A CMESimulation or RDMESimulation object with the "cytoplasm" region defined
        inducerType:
            The type of inducer used for the switch, either TMG or IPTG

    Returns:
        A pointer to the simulation object that was passed in
    """
	# Bulletproofing
	if (inducerType != 'TMG') and (inducerType != 'IPTG'):
		raise Exception('ERROR: When attaching E. coli two state switch, second argument not valid string, either TMG or IPTG')

	# Define species
	species = ['R2','O','R2O','IR2','IR2O','I2R2','I2R2O','mY','Y','I','Iex','YI']
	sim.defineSpecies(species)
	
	# Get cytosol
	if isinstance(sim, RDMESimulation):
		cytregion=sim.modifyRegion('cytoplasm')
		extregion=sim.modifyRegion('default')
		memregion=sim.modifyRegion('membrane')
	elif isinstance(sim, CMESimulation):
		cytregion=sim
		extregion=sim
		memregion=sim

	# Add reactions #
	scalar = 2.076e-9 # Rate conversion from experiments to stochastic
	# Lac operon regulation
	cytregion.addReaction(reactant=('R2','O'),   product='R2O',        rate=2.43e6*scalar)
	cytregion.addReaction(reactant=('IR2','O'),  product='IR2O',       rate=1.21e6*scalar)
	cytregion.addReaction(reactant=('I2R2','O'), product='I2R2O',      rate=2.43e4*scalar)
	cytregion.addReaction(reactant='R2O',        product=('R2','O'),   rate=6.30e-4)
	cytregion.addReaction(reactant='IR2O',       product=('IR2','O'),  rate=6.30e-4)
	cytregion.addReaction(reactant='I2R2O',      product=('I2R2','O'), rate=3.15e-1)

	# Transcription, translation, and degredation
	cytregion.addReaction(reactant='O',          product=('O','mY'),   rate=1.26e-1)
	memregion.addReaction(reactant='mY',         product=('mY','Y'),   rate=4.44e-2)
	cytregion.addReaction(reactant='mY',         product='',           rate=1.11e-2)
	if isinstance(sim, RDMESimulation):
		memregion.addReaction(reactant='mY',         product='',           rate=1.11e-2)
	memregion.addReaction(reactant='Y',          product='',           rate=2.10e-4)

	# Inducer-repressor interactions
	if(inducerType == 'TMG'):
		cytregion.addReaction(reactant=('I','R2'),   product='IR2',        rate=2.27e4*scalar)
		cytregion.addReaction(reactant=('I','IR2'),  product='I2R2',       rate=1.14e4*scalar)
		cytregion.addReaction(reactant=('I','R2O'),  product='IR2O',       rate=6.67e2*scalar)
		cytregion.addReaction(reactant=('I','IR2O'), product='I2R2O',      rate=3.33e2*scalar)
	else:
		cytregion.addReaction(reactant=('I','R2'),   product='IR2',        rate=9.71e4*scalar)
		cytregion.addReaction(reactant=('I','IR2'),  product='I2R2',       rate=4.85e4*scalar)
		cytregion.addReaction(reactant=('I','R2O'),  product='IR2O',       rate=2.24e4*scalar)
		cytregion.addReaction(reactant=('I','IR2O'), product='I2R2O',      rate=1.12e4*scalar)
	cytregion.addReaction(reactant='IR2',        product=('I','R2'),   rate=2.00e-1)
	cytregion.addReaction(reactant='I2R2',       product=('I','IR2'),  rate=4.00e-1)
	cytregion.addReaction(reactant='IR2O',       product=('I','R2O'),  rate=1.00)
	cytregion.addReaction(reactant='I2R2O',      product=('I','IR2O'), rate=2.00)

	# Inducer transport
	memregion.addReaction(reactant='Iex',        product='I',          rate=2.33e-3)
	memregion.addReaction(reactant='I',          product='Iex',        rate=2.33e-3)
	memregion.addReaction(reactant=('Y','Iex'),  product='YI',         rate=3.03e4*scalar)
	memregion.addReaction(reactant='YI',         product=('Y','Iex'),  rate=1.20e-1)
	memregion.addReaction(reactant='YI',         product=('Y','I'),    rate=1.20e+1)

	# Set up diffusion rates
	if isinstance(sim, RDMESimulation):
		# Inducer Diffusion Processes
		cytregion.setDiffusionRate('I',1.28e-12)
		memregion.setDiffusionRate('I',1.28e-12)
		cytregion.setDiffusionRate('Iex',1.28e-12)
		extregion.setDiffusionRate('Iex',1.28e-12)
		sim.setTwoWayTransitionRate(species='I', one='cytoplasm', two='membrane', rate=1.28e-12 )
		sim.setTwoWayTransitionRate(species='Iex', one='membrane', two='default', rate=1.28e-12 )

		# Repressor Diffusion Processes 
		cytregion.setDiffusionRate('R2',1.0e-12)

		# mRNA Diffusion Processes
		memregion.setDiffusionRate('mY',0.1e-12)
		cytregion.setDiffusionRate('mY',0.1e-12)
		sim.setTwoWayTransitionRate(species='mY', one='cytoplasm', two='membrane', rate=0.1e-12 )

	return sim


def addPTSPathway(sim):
	"""Adds the Phosphoenolpyruvate-depenedent phosphotransferase reaction system.

    Reference: J.V. Rodriguez, J.A. Kaandorp, M. Dobrzynski, J.G. Blom
      Spatial stochastic modelling of the phosphoenolpyruvate-dependent
      phosphotransferase (PTS) pathway in Escherichia coli.
      Bioinform. 22:15 (2006), pp. 1895-1901.

    Args:
        sim:
            A RDMESimulation object with the "cytoplasm", "default" and "membrane" regions defined

    Returns:
        A pointer to the simulation object that was passed in
    """

	if not isinstance(sim, RDMESimulation):
		LMLogger.error('Tried to add the PTS pathway to a non-RDME simulation. Bailing out.')
		raise TypeError("Expected RDMESimulation")

	# Define species
	sim.defineSpecies(['EI','PEP','EIPEP','EIP','Pyr','HPr','HPrPEI','HPrP','IIA','IIAPHPr', 'IIAP','IICB','IICBPIIA','IICBP','Glc','GlcPIICB','GlcP'])

	# Get cytosol
	cyt=sim.modifyRegion('cytoplasm')
	mem=sim.modifyRegion('membrane')
	ext=sim.modifyRegion('default')

	# Set up diffusion rates
	um2minTom2s=1.0e-12/60.0
	diff1 = 197.8*um2minTom2s
	diff2 = 189.1*um2minTom2s
	diff3 = 378.0*um2minTom2s
	diff4 = 262.1*um2minTom2s
	diff5 = 300.0*um2minTom2s
	diff6 = 18000.0*um2minTom2s
	# Cytosol diffusion
	cyt.setDefaultDiffusionRate(0.0)
	cyt.setDiffusionRate('EI',     diff1)
	cyt.setDiffusionRate('EIP',    diff1)
	cyt.setDiffusionRate('EIPEP',  diff1)
	cyt.setDiffusionRate('EIPHPyr',diff2)
	cyt.setDiffusionRate('HPyr',   diff3)
	cyt.setDiffusionRate('HPyrP',  diff3)
	cyt.setDiffusionRate('HPrPIIA',diff4)
	cyt.setDiffusionRate('IIA',    diff5)
	cyt.setDiffusionRate('IIAP',   diff5)
	cyt.setDiffusionRate('IIAPHPr',diff5)
	cyt.setDiffusionRate('PEP',    diff6)
	cyt.setDiffusionRate('Pyr',    diff6)
	# Membrane diffusion
	mem.setDefaultDiffusionRate(0.0)
	mem.setDiffusionRate('IICB',     0.0)
	mem.setDiffusionRate('IICBPIIA', 0.0)
	mem.setDiffusionRate('IICBP',    0.0)
	mem.setDiffusionRate('GlcPIICB', 0.0)
	# Extracellular diffusion
	ext.setDefaultDiffusionRate(0.0)
	ext.setDiffusionRate('Glc', diff6)
	# EC-> Memebrane, Membrane-> Cytosol
	sim.setTransitionRate('Glc', 'default', 'membrane',    diff6)
	sim.setTransitionRate('Glc', 'membrane', 'diffusion',  diff6)
	sim.setTransitionRate('GlcP', 'membrane', 'cytoplasm', diff6)
	sim.setTransitionRate('GlcP', 'cytoplasm', 'membrane', diff6)
	sim.setTransitionRate('IIAP', 'membrane', 'cytoplasm', diff5)
	sim.setTransitionRate('IIAP', 'cytoplasm', 'membrane', diff5)
	sim.setTransitionRate('IIA',  'membrane', 'cytoplasm', diff5)
	sim.setTransitionRate('IIA',  'cytoplasm', 'membrane', diff5)

	# Set up reactions
	uMminTos=1.0/((sim.siteVolume()*6.022e23) * 1.0e-6)
	# Rxn 1
	cyt.addReaction(('EI','PEP'), 'EIPEP',  1960.0*uMminToS)
	cyt.addReaction('EIPEP', ('EI','PEP'),  480000.0/60.0)
	# Rxn 2
	cyt.addReaction('EIPEP', ('EIP','Pyr'), 108000.0/60.0)
	cyt.addReaction(('EIP','Pyr'), 'EIPEP', 294.0*uMminToS)
	# Rxn 3
	cyt.addReaction(('HPr','EIP'), 'HPrPEI', 14000.0*uMminToS)
	cyt.addReaction('HPrPEI', ('HPr','EIP'), 14000.0/60.0)
	# Rxn 4
	cyt.addReaction('HPrPEI', ('HPrP','EI'), 84000.0/60.0)
	cyt.addReaction(('HPrP','EI'), 'HPrPEI', 3360.0*uMminToS)
	# Rxn 5
	cyt.addReaction(('IIA','HPrP'), 'IIAPHPr', 21960.0*uMminToS)
	cyt.addReaction('IIAPHPr', ('IIA','HPrP'), 21960.0/60.0)
	# Rxn 6
	cyt.addReaction('IIAPHPr', ('IIAP','HPr'), 4392/60.0)
	cyt.addReaction(('IIAP','HPr'), 'IIAPHPr', 3384.0*uMminToS)
	# Rxn 7
	mem.addReaction(('IICB','IIAP'), 'IICBPIIA', 880.0*uMminToS)
	mem.addReaction('IICBPIIA', ('IIAP','IICIB'), 880.0/60.0)
	# Rxn 8
	mem.addReaction('IICBPIIA', ('IICBP','IIA'), 2640.0/60.0)
	mem.addReaction(('IICBP','IIA'), 'IICBPIIA', 960.0*uMminToS)
	# Rxn 9
	mem.addReaction(('IICBP','Glc'), 'GlcPIICB', 260.0*uMminToS)
	mem.addReaction('GlcPIICB', ('IICBP','Glc'), 389.0/60.0)
	# Rxn 10
	mem.addReaction('GlcPIICB', ('IICB','GlcP'), 4800.0/60.0)
	mem.addReaction(('IICB','GlcP'), 'GlcPIICB', 5.4e-3*uMminToS)



	# Set particle numbers
	sim.addParticles('EI',  cyt,1577)
	sim.addParticles('HPr', cyt,15766)
	sim.addParticles('IIA', cyt,12613)
	sim.addParticles('IICB',mem,3100)
	sim.addParticles('PEP', cyt,882890)
	sim.addParticles('Pyr', cyt,283780)
	sim.addParticles('Glc', ext,70767)
	sim.addParticles('GlcP',cyt,15766)

	return sim


