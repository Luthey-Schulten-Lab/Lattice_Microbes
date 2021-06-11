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
# Author(s): Joseph R. Peterson
#

import pyLM
import pyLM.CME
import pyLM.RDME
from pyLM.LMLogger import *
try:
    from libsbml import *
except ModuleNotFoundError:
    import sys
    if 'sphinx' not in sys.modules:
        raise


## Helper Functions ##
def importSBMLModel(filename):
	reader = SBMLReader()
	doc    = reader.readSBMLFromFile(filename)

	# Check errors
	if doc.getNumErrors() != 0:
		LMLogger.error("There were errors reading the SBML file: %s"%(filename))
		doc.printErrors()
		return None

	# Make sure it is correct version
	if doc.getLevel() != 3 or doc.getVersion() != 1:
		LMLogger.error("Document was incorrect version; must have Level 3, Version 1 document. Given: L%d, V%d" %(doc.getLevel(),doc.getVersion()))
		
	# Returns a pointer to the model
	model = doc.getModel()
	return model

def importSBMLModelL3V1(model, sim, region=None):
	# Get the units for the model
	modelSubstanceUnits = model.getSubstanceUnits()

	globalParameters = []
	globalParameterValues = {}

	# Get all parameter values
	for p in range(model.getNumParameters()-1):
		pID = model.getParameter(p).getId()
		globalParameters.append(pID)
		globalParameterValues[pID] = model.getParemeter(p).getValue()
	
	# Get all compartment sizes
	for c in range(model.getNumCompartments()-1):
		cID = model.getCompartment(c).getID()
		globalParameters.append(cID)
		globalParameterValues[cID] = model.getCompartment(c).getSize()

	# Not in the C implementation:
	speciesCounts = {}

	# Back to C implementation

	speciesIndices = {}
	numSpecies = model.getNumSpecies()
	for s in range(numSpecies-1):
		specie = model.getSpecies(s)
		speciesIndices[specie.getId()] = s

		# Get units for amount
		substanceUnits = modelSubstanceUnits;
		if specie.isSetSubstanceUnits():
			substanceUnits = specie.getSubstanceUnits()
		if substanceUnits != 'item':
			LMLogger.error('Unsupported species substance units: %s'%(substanceUnits))

		# Make sure we can process the specie
		if not specie.getHasOnlySubstranceUnits():
			LMLogger.error('Unsupported species property: "hasOnlySubstranceUnits" must be true')
		if specie.getBoundaryCondition():
			LMLogger.error('Unsupported species property: "boundaryCondition" must be false')
		if specie.getConstant():
			LMLogger.error('Unsupported species property: "constant" must be false')
		if specie.isSetConveersionFactor():
			LMLogger.error('Unsupported species property: "conversionFactor" must not be set')

		# Get initial count for the species
		if specie.isSetInitialAmount():
			speciesCounts[specie.getID()] = round(specie.getInitialAmount())
		elif specie.isSetInitialConcentration():
			LMLogger.error('Unsupported species property: "initialConcentration" must not be set')
		else:
			LMLogger.error('Unsupported species property: "initialAmount" or "initialConcentration" must be set')

		# Add the specie to the list of defined species
		sim.defineSpecies([specie.getName()])

	# Process reactions
	numReactions =  model.getNumReactions()
	for r in range(numReactions-1):
		rxn = model.getReaction(r)

		# Check that the reaction is supported
		if rxn.getReversible():
			LMLogger.error('Unsupported reaction property: "reversible" must be false')
		if not reaction.isSetKineticLaw():
			LMLogger.error('Unsupported reaction property: must have a kinetic law')

		# Get reactant list
		rcts = []
		for rct in range(rxn.getNumReactants()-1):
			lrct = rxn.getReactant(rct)
			
			if not lrct.isSetStoichiometry():
				LMLogger.error('Unsupported reaction property: "stoichiometry" for reactant must be set')
			if not lrct.getConstant():
				LMLogger.error('Unsupported reaction property: "stoichiometry" for reactant must be constant')

			for i in range(lrct.getStoichiometry()-1):
				rcts.append(lrct)


		# Get product list
		prds = []
		for prd in range(rxn.getNumProducts()-1):
			lprd = rxn.getProduct(prd)
			
			if not lprd.isSetStoichiometry():
				LMLogger.error('Unsupported reaction property: "stoichiometry" for product must be set')
			if not lprd.getConstant():
				LMLogger.error('Unsupported reaction property: "stoichiometry" for product must be constant')

			for i in range(lprd.getStoichiometry()-1):
				prds.append(lrct)

		rate = getRate(rxn, rcst, prds, sim, globalParameters, globalParameterValues)

		if region != None:
			region.addReaction(tuple(rcts), tuple(prds), rate)
		elif isinstance(sim, CMESimulation):
			sim.addReaction(tuple(rcts), tuple(prds), rate)
		elif isinstance(sim, RDMESimulation):
			pass
			# Determine compartment
			
		else:
			LMLogger.error('Unknown simulation type.')


def getRate(rxn, rcts, prds, sim, globalP, globalPV):
	kinetics = rxn.getKineticLaw()

	localParameters = []
	localParameterValues = {}

	# Copy global parameters to local parameters
	for i in globalParameters:
		localParameters.append(i)
		localParameterValues[i] = globalPV[i]

	# Get local parameters
	for i in range(kinetics.getNumLocalParameters()-1):
		lp = kinetics.getLocalParameter(i)

		if not lp.isSetValue():
			LMLogger.error('Unsupported reaction property: "value" for local parameter must be set')
		
		localParameters.append(lp.getId())
		localParameterValues[lp.getId()] = lp.getValue()

	# Determine reaction type
	m = kinetics.getMath()
	if isFirstOrder(m, localParameters):
		return importFirstOrder()
	elif isSecondOrder(m, localParameterS):
		return importSecondOrder()
	elif isSecondSelfOrder(m, localParameterS):
		return importSecondSelfOrder()
	else:
		LMLogger.error("Unsupported kinetic law: %s"%(kinetics.getFormula()))

	


def readSBMLtoCME(sim, filename):
	"""Read an SBML file for the reaction model for a CME simulation

    Args:
        sim:
            A CMESimulation object
        filename:
            A SBML filename/filepath

    Returns:
        The simulation object
    """
	# Bulletproofing
	if not isinstance(sim, CMESimulation):
		LMLogger.error("Must give a CMESimulation to 'readSBMLtoCME'.")

	model = getSBMLDoc(filename)

	importSBMLModelL3V1(model, sim)
	
	return sim

def readSBMLtoRDME(sim, filename):
	"""Read an SBML file for the reaction model for a RDME simulation

    Args:
        sim:
            A RDMESimulation object
        filename:
            A SBML filename/filepath

    Returns:
        The simulation object
    """
	# Bulletproofing
	if not isinstance(sim, RDMESimulation):
		LMLogger.error("Must give a RDMESimulation to 'readSBMLtoRDME'.")
	
	model = getSBMLDoc(filename)

	importSBMLModelL3V1(model, sim)

	return sim


def readSBMLtoRegion(sim, region, filename):
	"""Read an SBML file for the reaction model into a specific region of the RDME simulation

    Args:
        sim:
            A RDMESimulation object
        region:
            A region name that already exists in the RDMESimulation
        filename:
            A SBML filename/filepath

    Returns:
        The simulation object
    """
	# Bulletproofing
	if not isinstance(sim, RDMESimulation):
		LMLogger.error("Must give a RDMESimulation to 'readSBMLtoRegion'.")
	if not region in sim.regions:
		LMLogger.error("No region '%s' in RDMESimulation."%(region))

	reg = sim.modifyRegion(region)
	model = getSBMLDoc(filename)

	importSBMLModelL3V1(model, sim, region)
	
	return sim

