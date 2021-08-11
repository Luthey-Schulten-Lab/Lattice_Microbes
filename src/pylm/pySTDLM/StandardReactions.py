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

import pyLM
import pyLM.CME
import pyLM.RDME
from pyLM.LMLogger import LMLogger
try:
	import Bio
	from Bio import SeqIO
	hasBio = True
except ImportError:
	hasBio = False


def addMembraneTransporter(sim, transporter, number, name1, name2, dRate, kf, kr, region1='cytoplasm', region2='default', intoRegion='membrane'):
	"""Adds a membrane transport system to an RDME system

    Args:
        sim:
            An RDMESimulation to add the transporter system to
        transporter:
            Name of the transporter particle
        number:
            The number of transporter molecules
        name1:
            Name of the transported particle in region1
        name2:
            Name of the transported particle in region2
        dRate:
            The rate of diffusion from region1<->intoRegion and region2<->intoRegion
        kf:
            The forward reaction rate 
        kr:
            The reverse reaction rate
        region1:
            The region on the first side of the membrane
        region2:
            The region on the second side of the membrane
        intoRegion:
            The region representing the membrane
    Returns:
        The simulation object so this can be a chained call
    """
	if not isinstance(sim, pyLM.RDME.RDMESimulation):
		LMLogger.error("In 'addMembraneTransporter', must be an RDMESimulation")
		raise Exception("ERROR: In 'addMembraneTransporter', must be an RDMESimulation")

	# Add species to the system
	sim.defineSpecies([name1,name2,transporter]) # This will not add multiple types

	# Get handles to the region
	r1=sim.modifyRegion(region1)
	r2=sim.modifyRegion(region2)
	rm=sim.modifyRegion(intoRegion)

	# Add diffusion rates
	sim.setTransitionRate(name1, via=region1, to=intoRegion, rate=dRate)
	sim.setTransitionRate(name1, via=intoRegion, to=region1, rate=dRate)
	sim.setTransitionRate(name2, via=region2, to=intoRegion, rate=dRate)
	sim.setTransitionRate(name2, via=intoRegion, to=region2, rate=dRate)

	# Add reaction rates
	rm.addReaction(reactant=(transporter, name1), product=(transporter,name2), rate=af)
	rm.addReaction(reactant=(transporter, name2), product=(transporter,name1), rate=kr)

	# Add membrane particles to the system
	sim.addParticles(species=transporter, region=intoRegion, rate=number)

	return sim


def addPassiveTransport(sim, specie, dRate, region1='cytoplasm', region2='default', region3='membrane'):
	"""Adds a passive transport system to an RDME system

    Args:
        sim:
            An RDMESimulation to add the passive transport to
        specie:
            The species that can diffuse across the membrane
        dRate:
            The rate of diffusion across the membrane
        region1:
            The region on the first side of the membrane
        region2:
            The region on the second side of the membrane
        region3:
            The region representing the membrane

    Returns:
        The simulation object so this can be a chained call
    """
	if not isinstance(sim, pyLM.RDME.RDMESimulation):
		LMLogger.error("In 'addMembraneTransporter', must be an RDMESimulation")
		raise Exception("ERROR: In 'addMembraneTransporter', must be an RDMESimulation")

	# Add species to the system
	sim.defineSpecies(specie)

	# Get handles to the region
	r1=sim.modifyRegion(region1)
	r2=sim.modifyRegion(region2)
	rm=sim.modifyRegion(region3)

	# Add the diffusion rates
	sim.setTransitionRate(specie, via=r1,to=rm, rate=dRate)
	sim.setTransitionRate(specie, via=rm,to=r1, rate=dRate)
	sim.setTransitionRate(specie, via=r2,to=rm, rate=dRate)
	sim.setTransitionRate(specie, via=rm,to=r2, rate=dRate)

	return sim


def addMichaelisMenten(sim, reactant, enzyme, product, k1f, k1b, k2, region='cytoplasm'):
	"""Adds a Michaelis Menten Reaction

    Args:
        sim:
            The RDME or CME reaction
        reactant:
            The reactant that reacts with the enzyme
        enzyme:
            The enzyme catalyzing the reaction
        product:
            The product of the reaction
        k1f:
            The forward reaction rate
        k1b:
            The backward reaction rate
        k2:
            The second forward rate
        region:
            The region in which the reaction should occur (RDME only)

    Returns:
        The simulation object so this can be a chained call
    """
	# Specify region
	region=1
	if isinstance(sim, pyLM.RDME.RDMESimulation):
		region = sim.modifyRegion(region)
	elif isinstance(sim, pyLM.CME.CMESimulation):
		region = sim
	else:
		LMLogger.warning('Could not add Michaelis-Menten reaction to object that is not an RDMESimualtion or CMESimulation.')
		return sim

	# Add particles
	intermediate = 'intermediate%s%s' % (reactant,enzyme)
	sim.defineSpecies([reactant,enzyme,product,intermediate])

	# Add the reactions
	region.addReaction((reactant,enzyme), intermediate, k1f)
	region.addReaction(intermediate, (reactant,enzyme), k1b)
	region.addReaction(intermediate, (product,enzyme), k2)

	return sim

def addReversibleMichaelisMenten(sim, reactant, enzyme, product, k1f, k1b, k2f, k2b, k3f, k3b, region='cytoplasm'):
	"""Add a Reversible Michaelis Menten Reaction

    Args:
        sim:
            The RDME or CME reaction
        reactant:
            The reactant that reacts with the enzyme
        enzyme:
            The enzyme catalyzing the reaction
        product:
            The product of the reaction
        k1f:
            The forward reaction rate
        k1b:
            The backward reaction rate
        k2f:
            The second forward rate
        k2b:
            the second backward rate
        k3f:
            The third forward rate
        k3b:
            the third backward rate
        region:
            the region in which the reaction should occur (RDME only)

    Returns:
        The simulation object so this can be a chained call
    """
	# Specify region
	region=1
	if isinstance(sim, pyLM.RDME.RDMESimulation):
		region = sim.modifyRegion(region)
	elif isinstance(sim, pyLM.CME.CMESimulation):
		region = sim
	else:
		LMLogger.warning('Could not add reversible Michaelis-Menten reaction to object that is not an RDMESimualtion or CMESimulation.')
		return sim

	# Add particles
	intermediate1 = 'intermediate1%s%s' % (reactant,enzyme)
	intermediate2 = 'intermediate2%s%s' % (reactant,enzyme)
	sim.defineSpecies([reactant,enzyme,product,intermediate1,intermediate2])

	# Add the reactions
	region.addReaction(reactant=(reactant,enzyme),product=intermediate1, rate=k1f)
	region.addReaction(reactant=intermediate1,product=(reactant,enzyme), rate=k1b)
	region.addReaction(reactant=intermediate1, product=intermediate2, rate=k2f)
	region.addReaction(reactant=intermediate2, product=intermediate1, rate=k2b)
	region.addReaction(reactant=intermediate2, product=(product,enzyme), rate=k3f)
	region.addReaction(reactant=(product,enzyme), product=intermediate2, rate=k3b)

	return sim

def createExpressionModel(sim, gb, kt, kd, kr, kdil = None, regions=None):
	"""Create a set of gene/mRNA/protein reactions based on a genebank file and a set of rates assuming constitutive expression

    Args:
        sim:
            The RDME or CME simulation
        gb:
            The genbank filename. File should be readable by BioPython.
        kt:
            mRNA transcription rate dictionary { locusTag -> rate }
        kd:
            mRNA degradation rate dictionary { locusTag -> rate }
        kr:
            Protein transcription rate dictionary {locusTag -> rate }
        kdil:
            Protein dillution/degradation rate dictionary {locusTag -> rate} (Optional; Default "None", meaning no dillution reaction will be specified)
        regions:
            Regions for the reactions to occur. Degradation is allowed in both regions. {locusTag -> (transcriptionRegion, translationRegion)} (required for RDMESimulations)
    
    Returns:
        SeqIO representation of the Genbank file with qualifiers added: qualifiers["dna_id"->(str,int),"rna_id"->(str,int), "protein_id"->(str,int)]
    """
	# Bulletproofing
	if not hasBio:
		raise Exception("Cannot read genbank file: BioPython not installed.")
	if not isinstance(sim, pyLM.RDME.RDMESimulation) and not isinstance(sim, pyLM.CME.CMESimulation):
		raise Exception("A CME or RDME simulation object must be specified")
	if len(set(kt.keys()).difference(set(kd.keys())).difference(set(kr.keys()))) != 0:
		raise Exception("kt, kd, and kr values must exist for the same set of locus tags")
	if kdil != None:
		if len(set(kt.keys()).difference(set(kdil.keys()))) != 0:
			raise Exception("kt, kd, kr, and kdil values must exist for the same set of locus tags")
	if regions != None:
		if len(set(kt.keys()).differene(set(regions.keys()))) != 0:
			raise Exception("rates and regions of rections must exist for the same set of locus tags")
	if isinstance(sim, pyLM.RDME.RDMESimulation) and regions == None:
		raise Exception("Regions must be specified when creating an expression model for an RDMESimulation")

	# Read simulation file
	loci = [k for k in kt.keys()]
	readLoci = []
	genome = SeqIO.read(gb,"genbank")
	for feature in genome.features:
		if feature.type in ["gene","rRNA"]:
			locus_tag = feature.qualifiers["locus_tag"][0]
			if locus_tag in loci and locus_tag not in readLoci:
				readLoci.append(locus_tag)
				# Add reactants
				dna = "D_%s"%(locus_tag)
				rna = "R_%s"%(locus_tag)
				protein = "P_%s"%(locus_tag)
				sim.defineSpecies([dna,rna,protein])
				feature.qualifiers["dna_id"]     = (dna,     sim.particleMap[dna])
				feature.qualifiers["rna_id"]     = (rna,     sim.particleMap[rna])
				feature.qualifiers["protein_id"] = (protein, sim.particleMap[protein])

				# Add reactions
				if isinstance(sim, pyLM.RDME.RDMESimulation):
					transcriptionRegion = sim.modifyRegion(regions[locus_tag][0])
					translationRegion   = sim.modifyRegion(regions[locus_tag][1])

					transcriptionRegion.addReaction(dna, (dna,rna), kt[locus_tag]) # Transcription
					transcriptionRegion.addReaction(rna, '',        kd[locus_tag]) # Degradation
					translationRegion.addReaction(rna,   '',        kd[locus_tag]) # Degradation
					translationRegion.addReaction(rna, (rna,protein), kr[locus_tag]) # Translation
					if kdil != None:
						transcriptionRegion.addReaction(protein, '',  kdil[locus_tag]) # Dillution
						translationRegion.addReaction(protein, '',    kdil[locus_tag]) # Dillution				
				else:
					sim.addReaction(dna, (dna,rna), kt[locus_tag]) # Transcription
					sim.addReaction(rna, '',        kd[locus_tag]) # Degradation
					sim.addReaction(rna, (rna,protein), kr[locus_tag]) # Translation
					if kdil != None:
						sim.addReaction(protein, '',    kdil[locus_tag]) # Dillution				

	unassigned = set(loci).difference(set(readLoci))
	if len(unassigned) != 0:
		LMLogger.warning("Not all loci found when creating expression model: %s"%(unassigned))

	return genome
