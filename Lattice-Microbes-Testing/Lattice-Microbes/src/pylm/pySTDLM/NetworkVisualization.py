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
from pyLM.LMLogger import *
try:
    from igraph import *
except ModuleNotFoundError:
    import sys
    if 'sphinx' not in sys.modules:
        raise

# For HDF5 file support:
import h5py
from lxml import etree
# For dynamic graphs:
try:
	import gexf
	hasGexf = True 
except ImportError:
	hasGexf = False

# #################################
# Reaction Network Visualizations #
# #################################
def plotReactionNetwork(sim, filename, collapseReversible=False):
	"""Plot the reaction scheme as a network 

    Args:
        sim:
            An RDMESimulation or CMESimulation object
        filename:
            A file to which to output
        collapseReversible:
            Collapsee reversible reactions into one node
    """
	# Bulletproofing
	if not (isinstance(sim, RDMESimulation) or isinstance(sim, CMESimulation)):
		LMLogger.warning("When trying to plot reaction network, was given object that was not RDMESimulation or CMESimulation")
		return

	# Call helper functions
	if isinstance(sim, RDMESimulation):
		plotRDMEReactionNetwork(sim,filename, collapseReversible)
	elif isinstance(sim, CMESimulation):
		plotCMEReactionNetwork(sim,filename, collapseReversible)
		

def plotCMEReactionNetwork(sim, filename, withRxnNodes=False, collapseReversible=False):
	"""Plot the reaction scheme as a network 

    Args:
        sim:
            An CMESimulation object
        filename:
            A file to which to output
        withRxnNodes:
            Plot the graph with reactions as nodes (default: false)
        collapseReversible:
            Collapsee reversible reactions into one node (default: false)
    """
	# Get handles to CME internals
	rs=sim.reactions
	nr=0
	ss=sim.species_id
	ns=len(ss)
	sc=sim.initial_counts

	# Create Graph
	g = Graph(directed=True)
	totNodes=ns
	if withRxnNodes:
		nr=len(rs)
		totNodes += nr
	g.add_vertices(totNodes)

	# Add reactions and species to graph
	count=0
	for v in g.vs:
		if count < nr:
			g.vs[v.index]["type"]="reaction"
			g.vs[v.index]["label"]=str(rs[count][0])+str("->")+str(rs[count][1])
			g.vs[v.index]["rate"]=rs[count][2]
		else:
			g.vs[v.index]["type"]="specie"
			g.vs[v.index]["label"]=ss[count-nr]
			g.vs[v.index]["icount"]=sc[ss[count-nr]]
			
		count+=1
	
	# Add edge dependencies
	count=0
	if withRxnNodes:
		for rxn in rs:
			# Add reactants
			if isinstance(rxn[0], tuple):
				rctnum=ss.index(r)
				g.add_edges([(rctnum+nr,count)])
			else:
				if rxn[0] != '': # FIXME
					rctnum=ss.index(rxn[0])
					g.add_edges([(rctnum+nr,count)])
			# Add products
			if isinstance(rxn[1], tuple):
				for p in rxn[1]:
					pdtnum=ss.index(p)
					g.add_edges([(count, pdtnum+nr)])
			else:
				if rxn[1] != '': # FIXME
					pdtnum=ss.index(rxn[1])
					g.add_edges([(count, pdtnum+nr)])

			count+=1	
	else:
		for rxn in rs:
			if isinstance(rxn[0],tuple):
				for r in rxn[0]:
					rctnum=ss.index(r)
					if isinstance(rxn[1],tuple):
						for p in rxn[1]:
							prdnum=ss.index(p)
							g.add_edges([(rctnum,prdnum)])
					else:
						if rxn[1] != '': # FIXME
							prdnum=ss.index(rxn[1])
							g.add_edges([(rctnum,prdnum)])
			else:
				rctnum=ss.index(rxn[0])
				if isinstance(rxn[1],tuple):
					for p in rxn[1]:
						prdnum=ss.index(p)
						g.add_edges([(rctnum,prdnum)])
				else:
					if rxn[1] != '': # FIXME
						prdnum=ss.index(rxn[1])
						g.add_edges([(rctnum,prdnum)])
				

	# Save the graph...
	g.save(filename,format="gml")

def plotRDMEReactionNetwork(sim, filename, collapseReversible=False):
	"""Plot the reaction scheme as a network 

    Args:
        sim:
            An RDMESimulation object
        filename:
            A file to which to output
        collapseReversible:
            Collapsee reversible reactions into one node
    """
	LMLogger.error("plotRDMEReactionNetwork unimplemented")
	return
	# Get handles to RDME internals
	regs=sim.regions
	sc=sim.initial_counts

	nr=0
	ns=len(sim.species_id)
	for r in regs:
		nr += len(r.reactions)

	# Create Graph
	g = Graph(directed=True)
	g.add_vertices(nr+ns)

	# Add reactions and species to graph
	count=0
	for v in g.vs:
		if count < nr:
			g.vs[v.index]["type"]="reaction"
			g.vs[v.index]["label"]=str(rs[count][0])+str("->")+str(rs[count][1])
			g.vs[v.index]["rate"]=rs[count][2]
		else:
			g.vs[v.index]["type"]="specie"
			g.vs[v.index]["label"]=ss[count-nr]
			g.vs[v.index]["icount"]=sc[ss[count-nr]]
			
		count+=1
	
	# Add edge dependencies
	count=0
	for rxn in rs:
		print(rxn)
		# Add reactants
		if isinstance(rxn[0], tuple):
			for r in rxn[0]:
				print(r)
				rctnum=ss.index(r)
				g.add_edges((rctnum+nr,count))
		else:
			rctnum=ss.index(rxn[0])
			g.add_edges((rctnum+nr,count))
		# Add products
		if isinstance(rxn[1], tuple):
			for p in rxn[1]:
				print(p)
				pdtnum=ss.index(p)
				g.add_edges((count, pdtnum+nr))
		else:
			pdtnum=ss.index(rxn[1])
			g.add_edges((count, pdtnum+nr))

		count+=1	

	# Save the graph...
	g.save(filename,format="gml")

def plotCMEDynamicReactionNetwork(sim, filename, outfile, stride = 1, showMax=False, showMin=False, threshold=-1, replicate=1):
	"""Plot the dynamics of the species on the network into a dynamic graph file (extension: .gexf)

    Args:
        sim:
            A CMESimulation object
        filename:
            An file from which to read
        outfile:
            A filename for which to output without the extension (.gexf)
        stride:
            Stride through the times (Default=1)
        showMax:
            Show the maximum value achieved over the timecourse on the node (Default=False)
        showMin:
            Show the maximum value achieved over the timecourse on the node (Default=False)
        threshold:
            Only show nodes that attain this value for at least one timepoint of the simulation. 
            Set -1 to show all nodes (Default=-1) replicate The replicate to show (Default=1)
    """
	if not hasGexf:
		raise Exception("pygexf not installed, cannot use plotCMEDynamicReactionNetwork.")

	# Get handles to the internals
	species=sim.species_id
	speciesNumbers=sim.particleMap
	reactions=sim.reactions

	# Open HDF5 file
	f=h5py.File(filename)
	repStr = sorted(f['Simulations'].keys())[replicate-1]
	rep = f['Simulations'][repStr]
	
	# Read time related variables
	times=rep['SpeciesCountTimes']
	       
	# Create a graph object
	gexfGraph=gexf.Graph(type="directed", mode="dynamic", label="AutomaticGeneratedTimeCounts", time_format="integer", start=str(0), end=str(times[len(times)-1]))
	# Add Node attributes
	countAttrID=gexfGraph.addNodeAttribute(title="count",type="integer",mode="dynamic", force_id="count")
	if showMin:
		minAttrID=gexfGraph.addNodeAttribute(title="minCount",type="integer",mode="static", force_id="minCount")
		minTAttrID=gexfGraph.addNodeAttribute(title="minCountTime",type="integer",mode="static", force_id="minCountTime")
	if showMax:
		maxAttrID=gexfGraph.addNodeAttribute(title="maxCount",type="integer",mode="static", force_id="maxCount")
		maxTAttrID=gexfGraph.addNodeAttribute(title="maxCountTime",type="integer",mode="static", force_id="maxCountTime")
	# Add Edge attributes
#	gexfGraph.addEdgeAttribute(title="rate",defaultValue="0",type="float",mode="static",force_id=99)

	# Add nodes to the graph
	for i in range(len(speciesNumbers)):
		gexfGraph.addNode(id=str(i),label=species[i])
	
	# Add edges to the graph
	countID=0
	for i in range(len(reactions)):
		rct=reactions[i][0]
		prd=reactions[i][1]
		if isinstance(rct, tuple):
			if isinstance(prd, tuple):
				for r in rct:
					for p in prd:
						gexfGraph.addEdge(id=str(countID), source=str(speciesNumbers[r]-1),target=str(speciesNumbers[p]-1))
						countID+=1	
			else:
				for r in rct:
					gexfGraph.addEdge(id=str(countID), source=str(speciesNumbers[r]-1),target=str(speciesNumbers[prd]-1))
					countID+=1	
					
		else:
			if isinstance(prd, tuple):
				for p in prd:
					gexfGraph.addEdge(id=str(countID), source=str(speciesNumbers[rct]-1),target=str(speciesNumbers[p]-1))
					countID+=1	
			else:
				gexfGraph.addEdge(id=str(countID), source=str(speciesNumbers[rct]-1),target=str(speciesNumbers[prd]-1))
				countID+=1	
					
		

	# Declare the time based attribute on the nodes and fill in values
	specCounts=rep['SpeciesCounts']
	tt=0
	mins=len(specCounts)*[666666666]
	minT=len(specCounts)*[0]
	maxs=len(specCounts)*[0]
	maxT=len(specCounts)*[0]
	# Account for thresholding
	exceedsThreshold=len(specCounts)*[False]
	if threshold > 0:
		for t in range(0,len(times),stride):
			for s in range(len(specCounts[t])):
				# Get mins and maxs
				sc=specCounts[t][s]
				if sc >= threshold:
					exceedsThreshold[s]=True
				if sc < mins[s]:
					mins[s] = sc
					minT[s] = t
				if sc > maxs[s]:
					maxs[s] = sc
					maxT[s] = t
				tt+=1
	else: # Don't account for thresholding
		for t in range(0,len(times),stride):
			for s in range(len(specCounts[t])):
				# Get mins and maxs
				sc=specCounts[t][s]
				if sc < mins[s]:
					mins[s] = sc
					minT[s] = t
				if sc > maxs[s]:
					maxs[s] = sc
					maxT[s] = t
				tt+=1

	# Add minimum and maximum values
	if showMax:
		for s in range(len(specCounts[0])):
			localNode=gexfGraph.nodes[str(s)]
			localNode.addAttribute(id="maxCount", value=str(maxs[s]))
			localNode.addAttribute(id="maxCountTime", value=str(maxT[s]))
	if showMin:
		for s in range(len(specCounts[0])):
			localNode=gexfGraph.nodes[str(s)]
			localNode.addAttribute(id="minCount", value=str(mins[s]))
			localNode.addAttribute(id="minCountTime", value=str(minT[s]))

	# Apply thresholding
	if threshold != -1:
		for s in range(len(specCounts[t])):
			localNode=gexfGraph.nodes[str(s)]
			if exceedsThreshold[s]:
				localNode.start=str(0)
				localNode.end=str(times[len(times)-1])
			else:
				localNode.start=str(0)
				localNode.end=str(0)
			


	# Save the graph
	with open(outfile+".gexf","w") as text_file:
		text_file.write(etree.tostring(gexfGraph.getXML(),pretty_print=True,encoding='utf-8',xml_declaration=True))

