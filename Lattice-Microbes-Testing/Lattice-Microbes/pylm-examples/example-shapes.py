from pyLM import *
from pyLM.units import *

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-o', '--outputFile', required=True)
args = ap.parse_args()

import pySTDLM.StandardCells as sc

latticeSpacing = 16 #nm

# Create our simulation object
sim=RDME.RDMESimulation(dimensions=micron(1.024,1.024,8.192), spacing=nm(latticeSpacing))

# Build a capsid cell
species = ['A','B','C']
sim.defineSpecies(species)

# These commands show how to define new regions in the simulation.
#  Specifically, we define a cytoplasm and a membrane region type
#  which can be used to delimit cells.
sim.addRegion('cytoplasm')
sim.addRegion('membrane')

##################
# Set Operations #
##################
# Test the Intersection Operation
#  Here we define two sphere shaed objects and then define an an object based on
#  the intersection of the two spheres.  Spheres require a center point to be defined,
#  a radius, and the number of the region.  The number of the region can be extracted
#  from the simulation by passing the region name to the "siteTypes" dictionary
sphere1 = lm.Sphere(lm.point(*micron(0.512, 0.512, 1.024)), micron(0.3), sim.siteTypes['cytoplasm']) # si_dist_t radius, site_t type)
sphere2 = lm.Sphere(lm.point(*micron(0.512, 0.512, 1.224)), micron(0.3), sim.siteTypes['cytoplasm']) # si_dist_t radius, site_t type)
intersection1 = lm.Intersection(sphere1,sphere2, sim.siteTypes['cytoplasm'])

# Test the Union Operation
#  Here we define a union between the previous intersection and a new sphere
sphere3 = lm.Sphere(lm.point(*micron(0.512, 0.7, 1.124)), micron(0.3), sim.siteTypes['cytoplasm'])
union1 = lm.Union(intersection1,sphere3,sim.siteTypes['cytoplasm'])
# Once a shape is defined, it must be added to the LM RDME simulation builder.
#  this is accomplished by calling "addRegion" on the "lm_builder" object that
#  is part of the RDMESimulation object
sim.lm_builder.addRegion(union1)

# Test Difference Operation
#  A cube is defined (lower and upper corners, followed by the region number and the up vector).
#  Also a capsule is defined based on the points of origins of the spheres at each end of the capsule and
#  the radius of the capsule.
box1     = lm.Cuboid(lm.point(*micron(0.1, 0.1, 2.048)),  lm.point(*micron(0.9, 0.9, 3.0)),  sim.siteTypes['cytoplasm'])
capsule1 = lm.Capsule(lm.point(*micron(0.512, 0.512, 2.048)),  lm.point(*micron(0.512, 0.512, 3.0)), micron(0.256),  sim.siteTypes['cytoplasm'])
# This difference is an assymetric difference, meaning it is the box with the cylinder cut out of it (b1-c1)
difference1 = lm.Difference(box1,capsule1, sim.siteTypes['membrane'])
sim.lm_builder.addRegion(difference1)

# Test simple difference
box2    = lm.Cuboid(lm.point(*micron(0.0,  0.0,  3.75)),  lm.point(*micron(0.25, 0.25, 4.0)),  sim.siteTypes['cytoplasm'])
sphere4 = lm.Sphere(lm.point(*micron(0.25, 0.25, 4.0)), micron(0.125), sim.siteTypes['cytoplasm'])
difference2 = lm.Difference(box2,sphere4, sim.siteTypes['cytoplasm'])
sim.lm_builder.addRegion(difference2)

#################
# Simple Shapes #
#################
# Draw a small Ellipsoid
#  Create an ellipse object by defining the center, followed by the radii of the three axes.  The ellipse also
#  takes an up vector and a pointing vector
ellipsoid1 = lm.Ellipse(lm.point(*micron(0.512, 0.512, 2.75)), micron(0.1), micron(0.05), micron(0.2), sim.siteTypes['membrane'], lm.vector(0.0,0.0,1.0), lm.vector(1.0,0.0,0.0))
sim.lm_builder.addRegion(ellipsoid1)

# Draw a Torus
#  Define a torus based on the center, the radius of the ring and the radius of the torus solid
torus1 = lm.Torus(lm.point(*micron(0.512, 0.512, 3.5)), micron(0.3), micron(0.1), sim.siteTypes['membrane'],lm.vector(0.0,0.0,1.0))
sim.lm_builder.addRegion(torus1)

# Draw some Cylinders
cylinder1 = lm.Cylinder(lm.point(micron(0.512),micron(0.512),micron(5.0)), lm.point(micron(0.512),micron(0.512),micron(5.5)), micron(0.3), sim.siteTypes['cytoplasm'] )
cylinder2 = lm.Cylinder(lm.point(micron(0.512),micron(0.512),micron(5.5)), lm.point(micron(0.512),micron(0.712),micron(6.0)), micron(0.3), sim.siteTypes['cytoplasm'] )
cylinder3 = lm.Cylinder(lm.point(micron(0.712),micron(0.512),micron(6.5)), lm.point(micron(0.512),micron(0.512),micron(6.5)), micron(0.3), sim.siteTypes['membrane'] )
sim.lm_builder.addRegion(cylinder1)
sim.lm_builder.addRegion(cylinder2)
sim.lm_builder.addRegion(cylinder3)

# Draw some Cones
cone1 = lm.Cone(lm.point(*micron(0.712,0.712,4.4)), micron(0.2), micron(0.5), sim.siteTypes['membrane'], lm.vector(0,0,1)) # z
cone2 = lm.Cone(lm.point(*micron(0.712,0.212,4.4)), micron(0.2), micron(0.4), sim.siteTypes['membrane'], lm.vector(0,1,0)) # y
cone3 = lm.Cone(lm.point(*micron(0.212,0.712,4.4)), micron(0.2), micron(0.3), sim.siteTypes['membrane'], lm.vector(1,0,0)) # x 
cone4 = lm.Cone(lm.point(*micron(0.212,0.212,4.4)), micron(0.2), micron(0.2), sim.siteTypes['membrane'], lm.vector(1,1,1)) # 1,1,1
cone5 = lm.Cone(lm.point(*micron(0.212,0.212,4.4)), micron(0.2), micron(0.2), sim.siteTypes['membrane'], lm.vector(-1,-1,-1)) # 1,1,1
sim.lm_builder.addRegion(cone1)
sim.lm_builder.addRegion(cone2)
sim.lm_builder.addRegion(cone3)
sim.lm_builder.addRegion(cone4)
sim.lm_builder.addRegion(cone5)

# Draw some Capsules
capsule1 = lm.Capsule(lm.point(micron(0.212),micron(0.212),micron(6.0)), lm.point(micron(0.212),micron(0.212),micron(6.5)), micron(0.1), sim.siteTypes['cytoplasm'] )
capsule2 = lm.Capsule(lm.point(micron(0.512),micron(0.512),micron(6.5)), lm.point(micron(0.212),micron(0.712),micron(7.0)), micron(0.2), sim.siteTypes['cytoplasm'] )
capsule3 = lm.Capsule(lm.point(micron(0.712),micron(0.512),micron(7.5)), lm.point(micron(0.512),micron(0.512),micron(7.5)), micron(0.1), sim.siteTypes['membrane'] )
sim.lm_builder.addRegion(capsule1)
sim.lm_builder.addRegion(capsule2)
sim.lm_builder.addRegion(capsule3)

# Draw some Capsule Shells
capsuleShell1 = lm.CapsuleShell(lm.point(micron(0.212),micron(0.212),micron(7.25)), lm.point(micron(0.212),micron(0.712),micron(7.25)), micron(0.1), micron(0.2), sim.siteTypes['membrane'] )
#capsuleShell2 = lm.Capsule(lm.point(micron(0.512),micron(0.512),micron()), lm.point(micron(0.212),micron(0.712),micron(7.0)), micron(0.4), micron(0.1), sim.siteTypes['cytoplasm'] )
sim.lm_builder.addRegion(capsuleShell1)
#sim.lm_builder.addRegion(capsuleShell2)

# Draw a Pyramid

# Draw a Prism




########################################
# Check objects that breach boundaries #
########################################
# X- 
sphere5 = lm.Sphere(lm.point(*micron(0.0, 0.512, 4.096)), micron(0.25), sim.siteTypes['cytoplasm']) # si_dist_t radius, site_t type)
sim.lm_builder.addRegion(sphere5)
# X+ 
sphere6 = lm.Sphere(lm.point(*micron(1.024, 0.512, 5.120)), micron(0.25), sim.siteTypes['cytoplasm']) # si_dist_t radius, site_t type)
sim.lm_builder.addRegion(sphere6)
# X+/Y+
sphere7 = lm.Sphere(lm.point(*micron(1.024, 1.024, 6.144)), micron(0.25), sim.siteTypes['cytoplasm']) # si_dist_t radius, site_t type)
sim.lm_builder.addRegion(sphere7)
# X+/Y+
sphere8 = lm.Sphere(lm.point(*micron(1.024, 1.024, 8.192)), micron(0.25), sim.siteTypes['cytoplasm']) # si_dist_t radius, site_t type)
sim.lm_builder.addRegion(sphere8)


# Modify the cytoplasm to add diffusion rates and reactions
sim.modifyRegion('cytoplasm') \
	.setDefaultDiffusionRate(2.5e-12) \
	.addReaction(reactant=('A','B'), product='C', rate=0.5) \
	.addReaction(reactant='C', product=('A','B'), rate=0.5) \


# Modify the membrane to add reactions
sim.modifyRegion('membrane') \
	.setDefaultDiffusionRate(2.5e-12) 

# Set diffusive properties between regions
sim.setTransitionRate(species='A', via='cytoplasm', to='membrane', rate= 2.5e-12)
sim.setTransitionRate(species='B', via='cytoplasm', to='membrane', rate= 2.5e-12)

# Populate the model with particles
sim.addParticles(species='A', region='cytoplasm', count=100)
sim.addParticles(species='B', region='membrane', count=100)

# Set simulation Parameters
sim.setTimestep(microsecond(50))
sim.setWriteInterval(ms(100))
sim.setLatticeWriteInterval(ms(100))
sim.setSimulationTime(ms(100))

# Save simulation
sim.save(args.outputFile)

