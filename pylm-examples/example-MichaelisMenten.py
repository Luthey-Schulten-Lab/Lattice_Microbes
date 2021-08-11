from pyLM import *
from pyLM.units import *

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()

import logging

LMLogger.setLMLogConsole(logging.WARNING)

# Create our RDME Simulation object
#  Here we define an RDME simulation including the inital dimensions and
#  the lattice spacing. The dimensions argument should be a tuple with 
#  three components specifying the x, y and z dimensions of the lattice.
#  The "spacing" argument defines the sidelength of the cubic lattice
#  sites.  Here we wrap the numbers in decorators that scale the numbers
#  appropriately (e.g. micron multiplies by 1e-6).  We find that these
#  generally improve the readability of the code.
sim=RDME.RDMESimulation(dimensions=micron(4.5,4.5,4.5), spacing=nm(4500/64.0))


# Define our chemical species
#  Species are defined exactly the same as in CME simulations
species = ['E', 'S', 'ES', 'P']
sim.defineSpecies(species)

# Add reactions to the simulation
#  In RDME, there is a notion of regions.  These are not spatial regions
#  so much as they are a type of region.  For example, these might be 
#  the cytoplasm of the cells in the simulation, or the extracellular space.
#  Reactions are defined per region.  To do this, first get a handle to 
#  the region by calling "modifyRegion()" on the simulation passing the
#  name of the region.  All RDME simulations are initialized with a region
#  called "default".  Once you have a handle to the region you do things
#  like specify a default diffusion rate for particles, specify specific
#  diffusion rates for particles, and define reactions that occur in those
#  regions.  These behaviors are shown below.
region=sim.modifyRegion('default')
region.setDefaultDiffusionRate(1e-12)
region.addReaction(reactant=('E','S'), product='ES', rate=1.07e-4)
region.addReaction(reactant='ES', product=('E','S'), rate=1e0)
region.addReaction(reactant='ES', product=('E','P'), rate=1e0)

# Set our initial species counts
#  Particles can be added to the simulation, with an optional
#  region to put them in using the following commands.  These will
#  randomly place the count of particles in any part of the domain
#  that is specified as that region type.  Other functionality in 
#  pySTDLM allows for placing particles in only certain spatial 
#  portions of the simulation domain.  In the case that no
#  region is specified, the default region is assumed.
sim.addParticles(species='E', count=909)
sim.addParticles(species='S', count=9091)
sim.addParticles(species='ES', count=0)
sim.addParticles(species='P', count=0)

# Define simulation parameters: run for 10 seconds, saving data every ms
sim.setTimestep(ms(1))
# The following command will tell the simulation how often to save
#  the lattice.  This must be a multiple of the intervale specified
#  in "setWriteInterval".  This can be less frequently as well, to 
#  reduce the size of the output file.
sim.setLatticeWriteInterval(0.01)
sim.setWriteInterval(0.01)
sim.setSimulationTime(10.0)
sim.save(args.outputFile)

# Here we run a single replicate.  The GPU accelerated MPDRDME Solver is
#  used.  
sim.run(filename=args.outputFile, method="lm::rdme::MpdRdmeSolver", replicates=1)


