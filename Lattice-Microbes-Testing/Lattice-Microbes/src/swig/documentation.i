
// File: index.xml

// File: structlm_1_1builder_1_1bounding__box.xml


%feature("docstring") lm::builder::bounding_box "

The bounds for a shape represented as a rectangular box.  

C++ includes: Shape.h
";

%feature("docstring") lm::builder::bounding_box::bounding_box "
`bounding_box(x1=0.0, y1=0.0, z1=0.0, x2=0.0, y2=0.0, z2=0.0)`  

Create a bounding box for coordinates.  

Parameters
----------
* `x1` :  
    Bottom x  
* `y1` :  
    Bottom y  
* `z1` :  
    Bottom z  
* `x2` :  
    Top x  
* `y2` :  
    Top y  
* `z2` :  
    Top z  
";

%feature("docstring") lm::builder::bounding_box::bounding_box "
`bounding_box(min, max)`  

Create a bounding box from points.  

Parameters
----------
* `min` :  
    Minimum coordinate  
* `max` :  
    Maximum coordinate  
";

%feature("docstring") lm::builder::bounding_box::joinWith "
`joinWith(j) -> bounding_box`  

Return a bounding box spanning the two bounding boxes.  

Parameters
----------
* `j` :  
    Other bounding box  

Returns
-------
new bounding box  
";

%feature("docstring") lm::builder::bounding_box::volume "
`volume() -> double`  

Returns the bounding box volume.  

Returns
-------
volume  
";

%feature("docstring") lm::builder::bounding_box::randomPointInside "
`randomPointInside() -> point`  

Retun a random point inside the bounding box.  

Returns
-------
point  
";

// File: classlm_1_1rdme_1_1ByteLattice.xml


%feature("docstring") lm::rdme::ByteLattice "

A Lattice that is based on packed bytes of memory, i.e. one byte per lattice site to hold particles.  

C++ includes: ByteLattice.h
";

%feature("docstring") lm::rdme::ByteLattice::nativeSerialize "
`nativeSerialize(destBuffer, lattice, latticeSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::nativeSerializeSites "
`nativeSerializeSites(destBuffer, lattice, latticeSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::copyNativeToRowMajorByte "
`copyNativeToRowMajorByte(destBuffer, sourceBuffer, xSize, ySize, zSize, particlesPerSite, bufferSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::copyRowMajorByteToNative "
`copyRowMajorByteToNative(destBuffer, sourceBuffer, xSize, ySize, zSize, particlesPerSite, bufferSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::copySitesRowMajorByteToNative "
`copySitesRowMajorByteToNative(destBuffer, sourceBuffer, xSize, ySize, zSize, bufferSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::copySitesNativeToRowMajorByte "
`copySitesNativeToRowMajorByte(destBuffer, sourceBuffer, xSize, ySize, zSize, bufferSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::getMaxSiteType "
`getMaxSiteType() -> site_t`  

Get the maximum number of site types possible in the lattice.  
";

%feature("docstring") lm::rdme::ByteLattice::getMaxParticle "
`getMaxParticle() -> particle_t`  

Get the maximum number of particle types possible in the lattice.  
";

%feature("docstring") lm::rdme::ByteLattice::getMaxOccupancy "
`getMaxOccupancy() -> site_size_t`  

Get the maximum number of particles that can live in a site.  
";

%feature("docstring") lm::rdme::ByteLattice::ByteLattice "
`ByteLattice(size, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::ByteLattice::ByteLattice "
`ByteLattice(xSize, ySize, zSize, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::ByteLattice::~ByteLattice "
`~ByteLattice()`  
";

%feature("docstring") lm::rdme::ByteLattice::getNeighboringSites "
`getNeighboringSites(index, neighboringIndices)`  

Get the sites that are neighbor to the indicated site.  

Parameters
----------
* `index` :  
    Index of the site for which to get neighbors  
* `neighboringIndicies` :  
    An array to hold the indicies of the neighbor sites  
";

%feature("docstring") lm::rdme::ByteLattice::getSiteType "
`getSiteType(x, y, z) -> site_t`  

Get the site type at the specified location.  
";

%feature("docstring") lm::rdme::ByteLattice::getSiteType "
`getSiteType(index) -> site_t`  

Get the site type at the specified location.  
";

%feature("docstring") lm::rdme::ByteLattice::setSiteType "
`setSiteType(x, y, z, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::ByteLattice::setSiteType "
`setSiteType(index, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::ByteLattice::copySites "
`copySites(destBuffer, latticeSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::getOccupancy "
`getOccupancy(x, y, z) -> site_size_t`  

Get the number of particles in the specified lattice site.  
";

%feature("docstring") lm::rdme::ByteLattice::getOccupancy "
`getOccupancy(index) -> site_size_t`  

Get the number of particles in the specified lattice site.  
";

%feature("docstring") lm::rdme::ByteLattice::getParticle "
`getParticle(x, y, z, particleIndex) -> particle_t`  

Gets the current state of a give site in the lattice.  

Parameters
----------
* `x` :  
    The zero based x index of the site to retrieve.  
* `y` :  
    The zero based y index of the site to retrieve.  
* `z` :  
    The zero based z index of the site to retrieve.  

Returns
-------
The value in the lattice at the specified site.  
";

%feature("docstring") lm::rdme::ByteLattice::getParticle "
`getParticle(index, particleIndex) -> particle_t`  

Get the particle at the specified site with at the specified number in the particle list.  
";

%feature("docstring") lm::rdme::ByteLattice::addParticle "
`addParticle(x, y, z, particle)`  

Sets the current state of a give site in the lattice.  

Parameters
----------
* `x` :  
    The zero based x index of the site to set.  
* `y` :  
    The zero based y index of the site to set.  
* `z` :  
    The zero based z index of the site to set.  
*  :  
";

%feature("docstring") lm::rdme::ByteLattice::addParticle "
`addParticle(index, particle)`  

Add a particle to the specified site.  
";

%feature("docstring") lm::rdme::ByteLattice::removeParticles "
`removeParticles(x, y, z)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::ByteLattice::removeParticles "
`removeParticles(index)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::ByteLattice::removeAllParticles "
`removeAllParticles()`  

Empty all particles from the specified site.  
";

%feature("docstring") lm::rdme::ByteLattice::getParticleCounts "
`getParticleCounts() -> std::map< particle_t, uint >`  

Get the number of each particle type in the lattice.  

Particle searching methods.  
";

%feature("docstring") lm::rdme::ByteLattice::findParticles "
`findParticles(minParticleType, maxParticleType) -> std::vector< particle_loc_t >`  

Get the number of the specified particles types in the lattice.  
";

%feature("docstring") lm::rdme::ByteLattice::setFromRowMajorByteData "
`setFromRowMajorByteData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::setSitesFromRowMajorByteData "
`setSitesFromRowMajorByteData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::ByteLattice::getParticleLatticeView "
`getParticleLatticeView(particleLattice, Nw, Nz, Ny, Nx, Np)`  
";

%feature("docstring") lm::rdme::ByteLattice::getSiteLatticeView "
`getSiteLatticeView(siteLattice, Nz, Ny, Nx)`  
";

%feature("docstring") lm::rdme::ByteLattice::getLatticeMemorySize "
`getLatticeMemorySize() -> size_t`  
";

// File: classlm_1_1builder_1_1Capsule.xml


%feature("docstring") lm::builder::Capsule "

A Shape representing a cylinder with hemispherical ends.  

C++ includes: Capsule.h
";

%feature("docstring") lm::builder::Capsule::Capsule "
`Capsule(p1, p2, radius, type)`  

Create a Capsule.  

Parameters
----------
* `p1` :  
    Point of first end of cylinder  
* `p2` :  
    Point of second end of cylinder  
* `radius` :  
    Radius of the capsule cylinder/hemispheres ends  
* `type` :  
    Site type of the capsule  
";

%feature("docstring") lm::builder::Capsule::~Capsule "
`~Capsule()`  

Destroy the Capsule.  
";

%feature("docstring") lm::builder::Capsule::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Capsule::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Capsule::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Capsule::getP1 "
`getP1() -> point`  

Get point of first end of cylinder.  
";

%feature("docstring") lm::builder::Capsule::getP2 "
`getP2() -> point`  

Get point of second end of cylinder.  
";

%feature("docstring") lm::builder::Capsule::getRadius "
`getRadius() -> si_dist_t`  

Get the radius of the cylinder/hemisphere ends.  
";

%feature("docstring") lm::builder::Capsule::getVolume "
`getVolume() -> double`  

Get the total internal volume of the Capsule.  
";

// File: classlm_1_1builder_1_1CapsuleShell.xml


%feature("docstring") lm::builder::CapsuleShell "

A capsule shell Shape with a hollow inside.  

C++ includes: CapsuleShell.h
";

%feature("docstring") lm::builder::CapsuleShell::CapsuleShell "
`CapsuleShell(p1, p2, innerRadius, outerRadius, type)`  

Create a Capsule.  

Parameters
----------
* `p1` :  
    Point of first end of cylinder  
* `p2` :  
    Point of second end of cylinder  
* `innerRadius` :  
    Inner radius of the capsule shell  
* `outerRadius` :  
    Outer radius of the capsule shell  
* `type` :  
    Type of the sites in the capsule shell  
";

%feature("docstring") lm::builder::CapsuleShell::~CapsuleShell "
`~CapsuleShell()`  
";

%feature("docstring") lm::builder::CapsuleShell::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::CapsuleShell::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::CapsuleShell::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::CapsuleShell::getP1 "
`getP1() -> point`  

Get point of first end of cylinder.  
";

%feature("docstring") lm::builder::CapsuleShell::getP2 "
`getP2() -> point`  

Get point of second end of cylinder.  
";

%feature("docstring") lm::builder::CapsuleShell::getInnerRadius "
`getInnerRadius() -> si_dist_t`  

Get the inner radius of the cylinder.  
";

%feature("docstring") lm::builder::CapsuleShell::getOuterRadius "
`getOuterRadius() -> si_dist_t`  

Get the outer radius of the cylinder.  
";

%feature("docstring") lm::builder::CapsuleShell::getVolume "
`getVolume() -> double`  

Get the total enclosed voluem of the CapsuleShell.  
";

// File: classlm_1_1main_1_1CheckpointSignaler.xml


%feature("docstring") lm::main::CheckpointSignaler "

A type of worker thread that checkpoints at a specified interval.  

C++ includes: CheckpointSignaler.h
";

%feature("docstring") lm::main::CheckpointSignaler::CheckpointSignaler "
`CheckpointSignaler()`  

Creates a new CheckpointSignaler Worker thread.  
";

%feature("docstring") lm::main::CheckpointSignaler::~CheckpointSignaler "
`~CheckpointSignaler()`  
";

%feature("docstring") lm::main::CheckpointSignaler::startCheckpointing "
`startCheckpointing(checkpointInterval)`  

Tells the thread to start checkpointing every checkpointInterval seconds.  

Parameters
----------
* `checkpointInterval` :  
    How often in seconds to checkpoint  
";

%feature("docstring") lm::main::CheckpointSignaler::stopCheckpointing "
`stopCheckpointing()`  

Tell the thread to stop checkpointing.  
";

%feature("docstring") lm::main::CheckpointSignaler::wake "
`wake()`  

Wake the thread if inactive.  
";

// File: classlm_1_1cme_1_1CMESolver.xml


%feature("docstring") lm::cme::CMESolver "

C++ includes: CMESolver.h
";

%feature("docstring") lm::cme::CMESolver::CMESolver "
`CMESolver(neededDists)`  
";

%feature("docstring") lm::cme::CMESolver::~CMESolver "
`~CMESolver()`  
";

%feature("docstring") lm::cme::CMESolver::initialize "
`initialize(replicate, parameters, resources)`  

Initialize the simulation.  

Parameters
----------
* `replicate` :  
    Replicate number out of total replicates  
* `parameters` :  
    A map of all the parameters for the simulation  
* `A` :  
    list of resources assigned to the simulation  
";

%feature("docstring") lm::cme::CMESolver::setReactionModel "
`setReactionModel(reactionModel)`  
";

%feature("docstring") lm::cme::CMESolver::buildModel "
`buildModel(numberSpecies, numberReactions, initialSpeciesCounts, reactionTypesA, k, S, D, kCols=1)`  
";

%feature("docstring") lm::cme::CMESolver::setModelPropensityFunction "
`setModelPropensityFunction(reaction, propensityFunction, propensityFunctionArg)`  
";

%feature("docstring") lm::cme::CMESolver::setSpeciesUpperLimit "
`setSpeciesUpperLimit(species, limit)`  
";

%feature("docstring") lm::cme::CMESolver::setSpeciesLowerLimit "
`setSpeciesLowerLimit(species, limit)`  
";

%feature("docstring") lm::cme::CMESolver::setFptTrackingList "
`setFptTrackingList(speciesList)`  
";

%feature("docstring") lm::cme::CMESolver::addToParameterTrackingList "
`addToParameterTrackingList(parameter)`  
";

%feature("docstring") lm::cme::CMESolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

%feature("docstring") lm::cme::CMESolver::getSpeciesCountView "
`getSpeciesCountView(counts, number)`  
";

%feature("docstring") lm::cme::CMESolver::getReactionRateConstantsView "
`getReactionRateConstantsView(reactionNumber, rates, rateConstantCount)`  
";

// File: classlm_1_1CommandLineArgumentException.xml


%feature("docstring") lm::CommandLineArgumentException "

Exception from a command line argument.  

C++ includes: Exceptions.h
";

%feature("docstring") lm::CommandLineArgumentException::CommandLineArgumentException "
`CommandLineArgumentException(message)`  
";

%feature("docstring") lm::CommandLineArgumentException::CommandLineArgumentException "
`CommandLineArgumentException(message, arg1)`  
";

// File: structlm_1_1cme_1_1CMESolver_1_1CompetitiveMMPropensityArgs.xml

// File: classlm_1_1main_1_1ResourceAllocator_1_1ComputeResources.xml


%feature("docstring") lm::main::ResourceAllocator::ComputeResources "

A representation for the resources for a given node.  

C++ includes: ResourceAllocator.h
";

%feature("docstring") lm::main::ResourceAllocator::ComputeResources::toString "
`toString() -> string`  
";

// File: classlm_1_1builder_1_1Cone.xml


%feature("docstring") lm::builder::Cone "

A Shape that represents a Cone.  

C++ includes: Cone.h
";

%feature("docstring") lm::builder::Cone::Cone "
`Cone(center, radius, height, type, normal=vector(1.0, 0.0, 0.0))`  

Create a Cone.  

Parameters
----------
* `center` :  
    Point center of the circle of the base  
* `radius` :  
    Radius of the base  
* `height` :  
    Height of the cone  
* `type` :  
    The type of the sites within the cone  
* `normal` :  
    Normal to the center of the cone base  
";

%feature("docstring") lm::builder::Cone::~Cone "
`~Cone()`  

Destroy the Sphere.  
";

%feature("docstring") lm::builder::Cone::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cone::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cone::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cone::getCenter "
`getCenter() -> point`  

Get the center of the cone.  
";

%feature("docstring") lm::builder::Cone::getRadius "
`getRadius() -> si_dist_t`  

Get the radius of the cone.  
";

%feature("docstring") lm::builder::Cone::getHeight "
`getHeight() -> si_dist_t`  

Get the height of the cone.  
";

%feature("docstring") lm::builder::Cone::getVolume "
`getVolume() -> double`  

Get the volume bounded by the cone.  
";

// File: classlm_1_1builder_1_1Cuboid.xml


%feature("docstring") lm::builder::Cuboid "

A cube-like Shape.  

C++ includes: Cuboid.h
";

%feature("docstring") lm::builder::Cuboid::Cuboid "
`Cuboid(p1, p2, type)`  

Create a Cuboid that is axis aligned.  

Parameters
----------
* `p1` :  
    Point of the lower corner  
* `p2` :  
    Point of the upper corner  
* `type` :  
    The type of the sites within the cuboid  
";

%feature("docstring") lm::builder::Cuboid::Cuboid "
`Cuboid(center, w, h, d, type, at=vector(1.0, 0.0, 0.0), up=vector(0.0, 1.0, 0.0))`  

Create a Cuboid based on two orientations, and a width, height, depth.  

center The definition of the center of the cuboid  

Parameters
----------
* `w` :  
    Width (along x axis when unrotated)  
* `h` :  
    Height (along y axis when unrotated)  
* `d` :  
    Depth (along z axis when unrotated)  
* `at` :  
    A vector representing where the cuboid is pointing (default along x axis)  
* `at` :  
    A vector representing where the cuboid up is (default along y axis)  
";

%feature("docstring") lm::builder::Cuboid::~Cuboid "
`~Cuboid()`  

Destroy the Cuboid.  
";

%feature("docstring") lm::builder::Cuboid::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cuboid::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cuboid::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cuboid::getP1 "
`getP1() -> point`  

Get the lower point.  
";

%feature("docstring") lm::builder::Cuboid::getP2 "
`getP2() -> point`  

Get the upper point.  
";

%feature("docstring") lm::builder::Cuboid::getVolume "
`getVolume() -> double`  

Get the volume bounded by the cuboid shape.  
";

// File: classlm_1_1CUDA.xml


%feature("docstring") lm::CUDA "

Class for accessing CUDA functions.  

C++ includes: lm_cuda.h
";

%feature("docstring") lm::CUDA::getNumberDevices "
`getNumberDevices() -> int`  
";

%feature("docstring") lm::CUDA::setCurrentDevice "
`setCurrentDevice(device)`  
";

%feature("docstring") lm::CUDA::getCurrentDevice "
`getCurrentDevice() -> int`  
";

%feature("docstring") lm::CUDA::getFreeMemory "
`getFreeMemory(device) -> size_t`  
";

%feature("docstring") lm::CUDA::getComputeCapabilities "
`getComputeCapabilities(device, major, minor)`  
";

%feature("docstring") lm::CUDA::printCapabilities "
`printCapabilities(device)`  
";

%feature("docstring") lm::CUDA::getCapabilitiesString "
`getCapabilitiesString(device) -> std::string`  
";

// File: classlm_1_1rdme_1_1CudaByteLattice.xml


%feature("docstring") lm::rdme::CudaByteLattice "

C++ includes: CudaByteLattice.h
";

%feature("docstring") lm::rdme::CudaByteLattice::CudaByteLattice "
`CudaByteLattice(size, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::CudaByteLattice::CudaByteLattice "
`CudaByteLattice(xSize, ySize, zSize, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::CudaByteLattice::~CudaByteLattice "
`~CudaByteLattice()`  
";

%feature("docstring") lm::rdme::CudaByteLattice::copyToGPU "
`copyToGPU()`  
";

%feature("docstring") lm::rdme::CudaByteLattice::copyFromGPU "
`copyFromGPU()`  
";

%feature("docstring") lm::rdme::CudaByteLattice::getGPUMemorySrc "
`getGPUMemorySrc() -> void *`  
";

%feature("docstring") lm::rdme::CudaByteLattice::getGPUMemoryDest "
`getGPUMemoryDest() -> void *`  
";

%feature("docstring") lm::rdme::CudaByteLattice::swapSrcDest "
`swapSrcDest()`  
";

%feature("docstring") lm::rdme::CudaByteLattice::getGPUMemorySiteTypes "
`getGPUMemorySiteTypes() -> void *`  
";

%feature("docstring") lm::rdme::CudaByteLattice::getParticleMemorySize "
`getParticleMemorySize() -> size_t`  
";

%feature("docstring") lm::rdme::CudaByteLattice::setSiteType "
`setSiteType(x, y, z, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::CudaByteLattice::setSiteType "
`setSiteType(index, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::CudaByteLattice::addParticle "
`addParticle(x, y, z, particle)`  

Sets the current state of a give site in the lattice.  

Parameters
----------
* `x` :  
    The zero based x index of the site to set.  
* `y` :  
    The zero based y index of the site to set.  
* `z` :  
    The zero based z index of the site to set.  
*  :  
";

%feature("docstring") lm::rdme::CudaByteLattice::addParticle "
`addParticle(index, particle)`  

Add a particle to the specified site.  
";

%feature("docstring") lm::rdme::CudaByteLattice::removeParticles "
`removeParticles(x, y, z)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::CudaByteLattice::removeParticles "
`removeParticles(index)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::CudaByteLattice::removeAllParticles "
`removeAllParticles()`  

Empty all particles from the specified site.  
";

%feature("docstring") lm::rdme::CudaByteLattice::setFromRowMajorByteData "
`setFromRowMajorByteData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::CudaByteLattice::getSiteLatticeView "
`getSiteLatticeView(siteLattice, Nz, Ny, Nx)`  
";

%feature("docstring") lm::rdme::CudaByteLattice::getParticleLatticeView "
`getParticleLatticeView(particleLattice, Nw, Nz, Ny, Nx, Np)`  
";

// File: classlm_1_1CUDADriverException.xml


%feature("docstring") lm::CUDADriverException "

CUDA driver exception.  

C++ includes: lm_cuda.h
";

%feature("docstring") lm::CUDADriverException::CUDADriverException "
`CUDADriverException(error, file, line)`  
";

// File: classlm_1_1CUDAException.xml


%feature("docstring") lm::CUDAException "

CUDA runtime exception.  

C++ includes: lm_cuda.h
";

%feature("docstring") lm::CUDAException::CUDAException "
`CUDAException(error, file, line)`  
";

// File: classlm_1_1rdme_1_1CudaIntLattice.xml


%feature("docstring") lm::rdme::CudaIntLattice "

C++ includes: CudaIntLattice.h
";

%feature("docstring") lm::rdme::CudaIntLattice::CudaIntLattice "
`CudaIntLattice(size, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::CudaIntLattice::CudaIntLattice "
`CudaIntLattice(xSize, ySize, zSize, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::CudaIntLattice::~CudaIntLattice "
`~CudaIntLattice()`  
";

%feature("docstring") lm::rdme::CudaIntLattice::copyToGPU "
`copyToGPU()`  
";

%feature("docstring") lm::rdme::CudaIntLattice::copyFromGPU "
`copyFromGPU()`  
";

%feature("docstring") lm::rdme::CudaIntLattice::getGPUMemorySrc "
`getGPUMemorySrc() -> void *`  
";

%feature("docstring") lm::rdme::CudaIntLattice::getGPUMemoryDest "
`getGPUMemoryDest() -> void *`  
";

%feature("docstring") lm::rdme::CudaIntLattice::swapSrcDest "
`swapSrcDest()`  
";

%feature("docstring") lm::rdme::CudaIntLattice::getGPUMemorySiteTypes "
`getGPUMemorySiteTypes() -> void *`  
";

%feature("docstring") lm::rdme::CudaIntLattice::setSiteType "
`setSiteType(x, y, z, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::CudaIntLattice::setSiteType "
`setSiteType(index, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::CudaIntLattice::addParticle "
`addParticle(x, y, z, particle)`  

Sets the current state of a give site in the lattice.  

Parameters
----------
* `x` :  
    The zero based x index of the site to set.  
* `y` :  
    The zero based y index of the site to set.  
* `z` :  
    The zero based z index of the site to set.  
*  :  
";

%feature("docstring") lm::rdme::CudaIntLattice::addParticle "
`addParticle(index, particle)`  

Add a particle to the specified site.  
";

%feature("docstring") lm::rdme::CudaIntLattice::removeParticles "
`removeParticles(x, y, z)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::CudaIntLattice::removeParticles "
`removeParticles(index)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::CudaIntLattice::removeAllParticles "
`removeAllParticles()`  

Empty all particles from the specified site.  
";

%feature("docstring") lm::rdme::CudaIntLattice::setFromRowMajorByteData "
`setFromRowMajorByteData(buffer, bufferSize)`  
";

// File: classlm_1_1builder_1_1Cylinder.xml


%feature("docstring") lm::builder::Cylinder "

A Shape representing a cylinder.  

C++ includes: Cylinder.h
";

%feature("docstring") lm::builder::Cylinder::Cylinder "
`Cylinder(p1, p2, radius, type)`  

Create a Cylinder.  

Parameters
----------
* `p1` :  
    Point of first end of cylinder  
* `p2` :  
    Point of second end of cylinder  
* `radius` :  
    Radius of the capsule cylinder/hemispheres ends  
* `type` :  
    Site type of the capsule  
";

%feature("docstring") lm::builder::Cylinder::~Cylinder "
`~Cylinder()`  

Destroy the Capsule.  
";

%feature("docstring") lm::builder::Cylinder::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cylinder::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cylinder::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Cylinder::getP1 "
`getP1() -> point`  

Get point of first end of cylinder.  
";

%feature("docstring") lm::builder::Cylinder::getP2 "
`getP2() -> point`  

Get point of second end of cylinder.  
";

%feature("docstring") lm::builder::Cylinder::getRadius "
`getRadius() -> si_dist_t`  

Get the radius of the cylinder/hemisphere ends.  
";

%feature("docstring") lm::builder::Cylinder::getVolume "
`getVolume() -> double`  

Get the total internal volume of the Capsule.  
";

// File: classlm_1_1main_1_1DataOutputQueue.xml


%feature("docstring") lm::main::DataOutputQueue "

A queue class that writes out data in the order it is recieved.  

C++ includes: DataOutputQueue.h
";

%feature("docstring") lm::main::DataOutputQueue::setInstance "
`setInstance(instance)`  

Sets the current DataOutputQueue that is active.  

Parameters
----------
* `instance` :  
    The new \"current\" queue  
";

%feature("docstring") lm::main::DataOutputQueue::getInstance "
`getInstance() -> DataOutputQueue *`  

Get the current DataOuputQueue tht is active.  

Returns
-------
Active output queue  
";

%feature("docstring") lm::main::DataOutputQueue::DataOutputQueue "
`DataOutputQueue()`  

Create a DataOutputQueue.  
";

%feature("docstring") lm::main::DataOutputQueue::~DataOutputQueue "
`~DataOutputQueue()`  
";

%feature("docstring") lm::main::DataOutputQueue::pushDataSet "
`pushDataSet(type, replicate, message, payload=NULL, payloadSize=0, payloadSerializer=NULL)`  

Put some data on the queue to be output.  

Parameters
----------
* `type` :  
    Type of data that is written to the data as a metaheader  
* `message` :  
    A message to put to protocol buffers  
* `payload` :  
    The data to be output  
* `payloadSize` :  
    The size of the payload in bytes  
* `payloadSerializer` :  
    A function that converts the object/data to a stream of bytes (If this is NULL the object is written out as pure bytes)  
";

%feature("docstring") lm::main::DataOutputQueue::pushDataSet "
`pushDataSet(data, dataSize)`  

Put some data on the queue.  

Parameters
----------
* `data` :  
    The data to be output  
* `dataSize` :  
    Size of the data in bytes  
";

%feature("docstring") lm::main::DataOutputQueue::pushDataSet "
`pushDataSet(dataSet)`  

Put a DataSet on the queue.  

Parameters
----------
* `dataSet` :  
    The data set to be put on the queue for output  
";

%feature("docstring") lm::main::DataOutputQueue::pushDataSet "
`pushDataSet(md, replicate, payload)`  
";

%feature("docstring") lm::main::DataOutputQueue::queryH5 "
`queryH5(mode, hdr, payload, replicate, path, attr=\"\")`  
";

%feature("docstring") lm::main::DataOutputQueue::queryH5attr "
`queryH5attr(replicate, path, attr) -> T`  
";

// File: classlm_1_1main_1_1DataOutputQueue_1_1DataSet.xml

// File: classlm_1_1builder_1_1Difference.xml


%feature("docstring") lm::builder::Difference "

A Shape that represents a Difference between the first and second object.  

C++ includes: Difference.h
";

%feature("docstring") lm::builder::Difference::Difference "
`Difference(s1, s2, type, symmetric=false)`  

Create a Difference.  

Parameters
----------
* `s1` :  
    The first shape to Difference  
* `s2` :  
    The second shape to Difference  
* `type` :  
    The type of the sites within the difference  
* `symmetric` :  
    Determine if the difference is symmetric. If false, difference is the 1st shape minus the second shape  
";

%feature("docstring") lm::builder::Difference::~Difference "
`~Difference()`  

Destroy the Sphere.  
";

%feature("docstring") lm::builder::Difference::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Difference::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Difference::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Difference::getVolume "
`getVolume(reintegrate=false) -> double`  

Get the volume bounded by the sphere.  
";

%feature("docstring") lm::builder::Difference::getVolume "
`getVolume() -> double`  

Get the total internal volume of the shape.  
";

// File: classlm_1_1builder_1_1Ellipse.xml


%feature("docstring") lm::builder::Ellipse "

C++ includes: Ellipse.h
";

%feature("docstring") lm::builder::Ellipse::Ellipse "
`Ellipse(center, r1, r2, r3, type, orientation1=vector(0.0, 0.0, 1.0), orientation2=vector(0.0, 1.0, 0.0))`  

Create a Ellipse.  

Parameters
----------
* `center` :  
    Point center of the circle of the slice plane through the sphere  
* `r1` :  
    First axis radius of the Ellipse  
* `r2` :  
    Second axis radius of the Ellipse  
* `r3` :  
    Third axis radius of the Ellipse  
* `type` :  
    The type of the sites within the sphere  
* `orientation1` :  
    The direction along which to lay the first axis, default: along the z-axis  
* `orientation2` :  
    The direction along which to lay the second axis, default: along the y-axis  
";

%feature("docstring") lm::builder::Ellipse::~Ellipse "
`~Ellipse()`  

Destroy the Torus.  
";

%feature("docstring") lm::builder::Ellipse::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Ellipse::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Ellipse::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Ellipse::setCenter "
`setCenter(center)`  

Set the center of the sphere.  

Parameters
----------
* `center` :  
    Point of the center  
";

%feature("docstring") lm::builder::Ellipse::getCenter "
`getCenter() -> point`  

Get the center of the sphere.  
";

%feature("docstring") lm::builder::Ellipse::getRadius1 "
`getRadius1() -> si_dist_t`  

Get the first radius of the sphere.  
";

%feature("docstring") lm::builder::Ellipse::getRadius2 "
`getRadius2() -> si_dist_t`  

Get the second radius of the sphere.  
";

%feature("docstring") lm::builder::Ellipse::getRadius3 "
`getRadius3() -> si_dist_t`  

Get the third radius of the sphere.  
";

%feature("docstring") lm::builder::Ellipse::getVolume "
`getVolume() -> double`  

Get the volume bounded by the sphere.  
";

// File: classlm_1_1Exception.xml


%feature("docstring") lm::Exception "

Base class for exceptions.  

A class for writing exceptions to a buffer that can then be written to either a console or a stream.  

C++ includes: Exceptions.h
";

%feature("docstring") lm::Exception::Exception "
`Exception(message=\"\")`  

Create an Exception.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg)`  

Create and Exception with one integer error code.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg1, arg2)`  

Create and Exception with two integer error codes.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg1, arg2, arg3)`  

Create and Exception with three integer error codes.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg)`  

Create and Exception with one error string.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg1, arg2)`  

Create and Exception with two error strings.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg1, arg2, arg3)`  

Create and Exception with three error strings.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg1, arg2)`  

Create and Exception with one integer error code and one error string.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg1, arg2, arg3)`  

Create and Exception with two integer error codes and one error string.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg, file, line)`  

Create and Exception with one integer error code, a file and a line number.  
";

%feature("docstring") lm::Exception::Exception "
`Exception(message, arg, file, line)`  

Create and Exception with one error string, a file and a line number.  
";

%feature("docstring") lm::Exception::~Exception "
`~Exception()`  

Destroy the Exception.  
";

%feature("docstring") lm::Exception::what "
`what() -> const char *`  

Get the error string.  

Returns
-------
messageBuffer A pointer to the error message  
";

// File: structlm_1_1cme_1_1CMESolver_1_1FirstOrderPropensityArgs.xml

// File: classlm_1_1cme_1_1FluctuatingNRSolver.xml


%feature("docstring") lm::cme::FluctuatingNRSolver "

C++ includes: FluctuatingNRSolver.h
";

%feature("docstring") lm::cme::FluctuatingNRSolver::FluctuatingNRSolver "
`FluctuatingNRSolver()`  
";

%feature("docstring") lm::cme::FluctuatingNRSolver::~FluctuatingNRSolver "
`~FluctuatingNRSolver()`  
";

%feature("docstring") lm::cme::FluctuatingNRSolver::setReactionModel "
`setReactionModel(rm)`  
";

%feature("docstring") lm::cme::FluctuatingNRSolver::buildModel "
`buildModel(numberSpecies, numberReactions, initialSpeciesCounts, reactionType, k, S, D, kCols=1)`  
";

%feature("docstring") lm::cme::FluctuatingNRSolver::buildModel "
`buildModel(numberSpecies, numberReactions, initialSpeciesCounts, reactionType, K, S, D, nvar, ntau, noiseRecalcFraction, kCols=1)`  
";

// File: structlm_1_1cme_1_1CMESolver_1_1FPTTracking.xml

// File: classlm_1_1cme_1_1GillespieDSolver.xml


%feature("docstring") lm::cme::GillespieDSolver "

C++ includes: GillespieDSolver.h
";

%feature("docstring") lm::cme::GillespieDSolver::GillespieDSolver "
`GillespieDSolver()`  
";

%feature("docstring") lm::cme::GillespieDSolver::~GillespieDSolver "
`~GillespieDSolver()`  
";

%feature("docstring") lm::cme::GillespieDSolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::GillespieDSolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::GillespieDSolver::buildModel "
`buildModel(numberSpecies, numberReactions, initialSpeciesCounts, reactionType, k, S, D, kCols=1)`  
";

%feature("docstring") lm::cme::GillespieDSolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: structgpu__info.xml


%feature("docstring") gpu_info "

C++ includes: ZDivMultiGPUMapper.h
";

// File: structlm_1_1rdme_1_1gpu__worker__thread__params.xml


%feature("docstring") lm::rdme::gpu_worker_thread_params "

C++ includes: MGPUMpdRdmeSolver.h
";

// File: structH5Lookup.xml


%feature("docstring") H5Lookup "

C++ includes: ArbitraryH5.h
";

// File: structH5MetaData.xml


%feature("docstring") H5MetaData "

C++ includes: ArbitraryH5.h
";

// File: classlm_1_1io_1_1hdf5_1_1HDF5Exception.xml


%feature("docstring") lm::io::hdf5::HDF5Exception "

Exception due to a failure in writing the HDF5 file.  

C++ includes: lm_hdf5.h
";

%feature("docstring") lm::io::hdf5::HDF5Exception::HDF5Exception "
`HDF5Exception(errorCode, file, line)`  

Create the Exception.  

Parameters
----------
* `errorCode` :  
    Code for the HDF5 error  
* `file` :  
    Source file in which the error occurred  
* `line` :  
    Line in the source file at which the error occurred  
";

%feature("docstring") lm::io::hdf5::HDF5Exception::printStackTrace "
`printStackTrace()`  
";

// File: classlm_1_1builder_1_1Hemisphere.xml


%feature("docstring") lm::builder::Hemisphere "

A hemisphere Shape.  

C++ includes: Hemisphere.h
";

%feature("docstring") lm::builder::Hemisphere::Hemisphere "
`Hemisphere(center, radius, orientation, type)`  

Create a Hemisphere.  

Parameters
----------
* `center` :  
    Point center of the circle of the slice plane through the Hemisphere  
* `radius` :  
    Radius of the Hemisphere  
* `orientation` :  
    Orientation normal to the center point of the center point in the direction of the curved part of the hemisphere  
* `type` :  
    The type of the sites within the hemisphere  
";

%feature("docstring") lm::builder::Hemisphere::~Hemisphere "
`~Hemisphere()`  

Destroy the Hemisphere.  
";

%feature("docstring") lm::builder::Hemisphere::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Hemisphere::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Hemisphere::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Hemisphere::getCenter "
`getCenter() -> point`  

Get the center of the Hemisphere.  
";

%feature("docstring") lm::builder::Hemisphere::getRadius "
`getRadius() -> si_dist_t`  

Get the radius of the Hemisphere.  
";

%feature("docstring") lm::builder::Hemisphere::getOrientation "
`getOrientation() -> vector`  

Get the orientation vector of the Hemisphere.  
";

%feature("docstring") lm::builder::Hemisphere::getVolume "
`getVolume() -> double`  

Get the volume contained by the Hemisphere.  
";

// File: classlm_1_1cme_1_1HillSwitch.xml


%feature("docstring") lm::cme::HillSwitch "

C++ includes: HillSwitch.h
";

%feature("docstring") lm::cme::HillSwitch::HillSwitch "
`HillSwitch()`  
";

%feature("docstring") lm::cme::HillSwitch::~HillSwitch "
`~HillSwitch()`  
";

%feature("docstring") lm::cme::HillSwitch::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::HillSwitch::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::HillSwitch::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1builder_1_1Intersection.xml


%feature("docstring") lm::builder::Intersection "

A Shape that represents a Intersection of two shapes.  

C++ includes: Intersection.h
";

%feature("docstring") lm::builder::Intersection::Intersection "
`Intersection(s1, s2, type)`  

Create a Intersection.  

Parameters
----------
* `s1` :  
    The first shape to Intersection  
* `s2` :  
    The second shape to Intersection  
* `type` :  
    The type of the sites within the intersection  
";

%feature("docstring") lm::builder::Intersection::~Intersection "
`~Intersection()`  

Destroy the Sphere.  
";

%feature("docstring") lm::builder::Intersection::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Intersection::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Intersection::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Intersection::getVolume "
`getVolume(reintegrate=false) -> double`  

Get the volume bounded by the sphere.  
";

%feature("docstring") lm::builder::Intersection::getVolume "
`getVolume() -> double`  

Get the total internal volume of the shape.  
";

// File: classlm_1_1rdme_1_1IntLattice.xml


%feature("docstring") lm::rdme::IntLattice "

A Lattice that is based on one particle per word, with sites strided per particle.  

C++ includes: IntLattice.h
";

%feature("docstring") lm::rdme::IntLattice::nativeSerialize "
`nativeSerialize(destBuffer, lattice, latticeSize)`  
";

%feature("docstring") lm::rdme::IntLattice::copyNativeToRowMajor "
`copyNativeToRowMajor(destBuffer, sourceBuffer, xSize, ySize, zSize, particlesPerSite, bufferSize)`  
";

%feature("docstring") lm::rdme::IntLattice::copyRowMajorToNative "
`copyRowMajorToNative(destBuffer, sourceBuffer, xSize, ySize, zSize, particlesPerSite, bufferSize)`  
";

%feature("docstring") lm::rdme::IntLattice::copySitesRowMajorByteToNative "
`copySitesRowMajorByteToNative(destBuffer, sourceBuffer, xSize, ySize, zSize, bufferSize)`  
";

%feature("docstring") lm::rdme::IntLattice::getMaxSiteType "
`getMaxSiteType() -> site_t`  

Get the maximum number of site types possible in the lattice.  
";

%feature("docstring") lm::rdme::IntLattice::getMaxParticle "
`getMaxParticle() -> particle_t`  

Get the maximum number of particle types possible in the lattice.  
";

%feature("docstring") lm::rdme::IntLattice::getMaxOccupancy "
`getMaxOccupancy() -> site_size_t`  

Get the maximum number of particles that can live in a site.  
";

%feature("docstring") lm::rdme::IntLattice::IntLattice "
`IntLattice(size, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::IntLattice::IntLattice "
`IntLattice(xSize, ySize, zSize, spacing, particlesPerSite)`  
";

%feature("docstring") lm::rdme::IntLattice::~IntLattice "
`~IntLattice()`  
";

%feature("docstring") lm::rdme::IntLattice::getNeighboringSites "
`getNeighboringSites(index, neighboringIndices)`  

Get the sites that are neighbor to the indicated site.  

Parameters
----------
* `index` :  
    Index of the site for which to get neighbors  
* `neighboringIndicies` :  
    An array to hold the indicies of the neighbor sites  
";

%feature("docstring") lm::rdme::IntLattice::getSiteType "
`getSiteType(x, y, z) -> site_t`  

Get the site type at the specified location.  
";

%feature("docstring") lm::rdme::IntLattice::getSiteType "
`getSiteType(index) -> site_t`  

Get the site type at the specified location.  
";

%feature("docstring") lm::rdme::IntLattice::setSiteType "
`setSiteType(x, y, z, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::IntLattice::setSiteType "
`setSiteType(index, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::IntLattice::getOccupancy "
`getOccupancy(x, y, z) -> site_size_t`  

Get the number of particles in the specified lattice site.  
";

%feature("docstring") lm::rdme::IntLattice::getOccupancy "
`getOccupancy(index) -> site_size_t`  

Get the number of particles in the specified lattice site.  
";

%feature("docstring") lm::rdme::IntLattice::getParticle "
`getParticle(x, y, z, particleIndex) -> particle_t`  

Gets the current state of a give site in the lattice.  

Parameters
----------
* `x` :  
    The zero based x index of the site to retrieve.  
* `y` :  
    The zero based y index of the site to retrieve.  
* `z` :  
    The zero based z index of the site to retrieve.  

Returns
-------
The value in the lattice at the specified site.  
";

%feature("docstring") lm::rdme::IntLattice::getParticle "
`getParticle(index, particleIndex) -> particle_t`  

Get the particle at the specified site with at the specified number in the particle list.  
";

%feature("docstring") lm::rdme::IntLattice::addParticle "
`addParticle(x, y, z, particle)`  

Sets the current state of a give site in the lattice.  

Parameters
----------
* `x` :  
    The zero based x index of the site to set.  
* `y` :  
    The zero based y index of the site to set.  
* `z` :  
    The zero based z index of the site to set.  
*  :  
";

%feature("docstring") lm::rdme::IntLattice::addParticle "
`addParticle(index, particle)`  

Add a particle to the specified site.  
";

%feature("docstring") lm::rdme::IntLattice::removeParticles "
`removeParticles(x, y, z)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::IntLattice::removeParticles "
`removeParticles(index)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::IntLattice::removeAllParticles "
`removeAllParticles()`  

Empty all particles from the specified site.  
";

%feature("docstring") lm::rdme::IntLattice::getParticleCounts "
`getParticleCounts() -> std::map< particle_t, uint >`  

Get the number of each particle type in the lattice.  

Particle searching methods.  
";

%feature("docstring") lm::rdme::IntLattice::findParticles "
`findParticles(minParticleType, maxParticleType) -> std::vector< particle_loc_t >`  

Get the number of the specified particles types in the lattice.  
";

%feature("docstring") lm::rdme::IntLattice::setFromRowMajorByteData "
`setFromRowMajorByteData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::IntLattice::setFromRowMajorData "
`setFromRowMajorData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::IntLattice::setSitesFromRowMajorByteData "
`setSitesFromRowMajorByteData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::IntLattice::getParticleLatticeView "
`getParticleLatticeView(particleLattice, Nw, Nz, Ny, Nx, Np)`  
";

%feature("docstring") lm::rdme::IntLattice::getSiteLatticeView "
`getSiteLatticeView(siteLattice, Nz, Ny, Nx)`  
";

%feature("docstring") lm::rdme::IntLattice::getLatticeMemorySize "
`getLatticeMemorySize() -> size_t`  
";

// File: classlm_1_1rdme_1_1IntMpdRdmeSolver.xml


%feature("docstring") lm::rdme::IntMpdRdmeSolver "

C++ includes: IntMpdRdmeSolver.h
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::IntMpdRdmeSolver "
`IntMpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::~IntMpdRdmeSolver "
`~IntMpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::initialize "
`initialize(replicate, parameters, resources)`  

Initialize the simulation.  

Parameters
----------
* `replicate` :  
    Replicate number out of total replicates  
* `parameters` :  
    A map of all the parameters for the simulation  
* `A` :  
    list of resources assigned to the simulation  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::buildModel "
`buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypeA, kA, SA, DA, kCols=1)`  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::buildDiffusionModel "
`buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData=true)`  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1InvalidArgException.xml


%feature("docstring") lm::InvalidArgException "

Exception for an argument to a function.  

C++ includes: Exceptions.h
";

%feature("docstring") lm::InvalidArgException::InvalidArgException "
`InvalidArgException(argMessage)`  
";

%feature("docstring") lm::InvalidArgException::InvalidArgException "
`InvalidArgException(arg, argMessage)`  
";

%feature("docstring") lm::InvalidArgException::InvalidArgException "
`InvalidArgException(arg, argMessage, argMessageParameter)`  
";

%feature("docstring") lm::InvalidArgException::InvalidArgException "
`InvalidArgException(arg, argMessage, argMessageParameter)`  
";

%feature("docstring") lm::InvalidArgException::InvalidArgException "
`InvalidArgException(arg, argMessage, argMessageParameter1, argMessageParameter2)`  
";

// File: classlm_1_1rdme_1_1InvalidParticleException.xml


%feature("docstring") lm::rdme::InvalidParticleException "

InvalidArgException for when a particle is.  

C++ includes: Lattice.h
";

%feature("docstring") lm::rdme::InvalidParticleException::InvalidParticleException "
`InvalidParticleException(particleIndex)`  

Create an eception based on the particle index.  

Parameters
----------
* `particleIndex` :  
    Index of the particle  
";

// File: classlm_1_1rdme_1_1InvalidSiteException.xml


%feature("docstring") lm::rdme::InvalidSiteException "

InvalidArgException for when an index into the lattice is out of bound or does not exist.  

C++ includes: Lattice.h
";

%feature("docstring") lm::rdme::InvalidSiteException::InvalidSiteException "
`InvalidSiteException(x, y, z)`  

Create an exception based on the lattice location.  

Parameters
----------
* `x` :  
    Lattice x point  
* `y` :  
    Lattice y point  
* `z` :  
    Lattice z Point  
";

%feature("docstring") lm::rdme::InvalidSiteException::InvalidSiteException "
`InvalidSiteException(index)`  

Create an exception based on the lattice index.  

Parameters
----------
* `index` :  
    Lattice index  
";

// File: classlm_1_1IOException.xml


%feature("docstring") lm::IOException "

Exception when input or output fails.  

C++ includes: Exceptions.h
";

%feature("docstring") lm::IOException::IOException "
`IOException(message, arg)`  
";

%feature("docstring") lm::IOException::IOException "
`IOException(message, arg)`  
";

// File: structlm_1_1cme_1_1TwoStateHillSwitch_1_1KHillPropensityArgs.xml

// File: structlm_1_1cme_1_1CMESolver_1_1KHillPropensityArgs.xml

// File: structlm_1_1cme_1_1CMESolver_1_1KHillTransportPropensityArgs.xml

// File: classlm_1_1cme_1_1LacHillSwitch.xml


%feature("docstring") lm::cme::LacHillSwitch "

C++ includes: LacHillSwitch.h
";

%feature("docstring") lm::cme::LacHillSwitch::LacHillSwitch "
`LacHillSwitch()`  
";

%feature("docstring") lm::cme::LacHillSwitch::~LacHillSwitch "
`~LacHillSwitch()`  
";

%feature("docstring") lm::cme::LacHillSwitch::setReactionModel "
`setReactionModel(rm)`  
";

// File: classlm_1_1rdme_1_1Lattice.xml


%feature("docstring") lm::rdme::Lattice "

Base class for lattice type objects.  

C++ includes: Lattice.h
";

%feature("docstring") lm::rdme::Lattice::rowMajorByteSerialize "
`rowMajorByteSerialize(destBuffer, lattice, bufferSize)`  
";

%feature("docstring") lm::rdme::Lattice::rowMajorIntSerialize "
`rowMajorIntSerialize(destBuffer, lattice, bufferSize)`  
";

%feature("docstring") lm::rdme::Lattice::rowMajorByteSerializeSites "
`rowMajorByteSerializeSites(destBuffer, lattice, bufferSize)`  
";

%feature("docstring") lm::rdme::Lattice::getMaxSiteType "
`getMaxSiteType() -> site_t`  

Get the maximum number of site types possible in the lattice.  
";

%feature("docstring") lm::rdme::Lattice::getMaxParticle "
`getMaxParticle() -> particle_t`  

Get the maximum number of particle types possible in the lattice.  
";

%feature("docstring") lm::rdme::Lattice::getMaxOccupancy "
`getMaxOccupancy() -> site_size_t`  

Get the maximum number of particles that can live in a site.  
";

%feature("docstring") lm::rdme::Lattice::Lattice "
`Lattice(size, spacing)`  

Create a Lattice object.  

Parameters
----------
* `size` :  
    Size of the lattice as (x,y,z)  
* `spacing` :  
    Physical spacing between lattice sites  
";

%feature("docstring") lm::rdme::Lattice::Lattice "
`Lattice(xSize, ySize, zSize, spacing)`  

Create a Lattice object.  

Parameters
----------
* `xSize` :  
    Number of sites in x  
* `ySize` :  
    Number of sites in y  
* `zSize` :  
    Number of sites in z  
* `spacing` :  
    Physical spacing between lattice sites  
";

%feature("docstring") lm::rdme::Lattice::~Lattice "
`~Lattice()`  

Destroy the Lattice object.  
";

%feature("docstring") lm::rdme::Lattice::getSize "
`getSize() -> lattice_coord_t`  

Get size of the Lattice.  
";

%feature("docstring") lm::rdme::Lattice::getXSize "
`getXSize() -> lattice_size_t`  

Get x dimension of the Lattice.  
";

%feature("docstring") lm::rdme::Lattice::getYSize "
`getYSize() -> lattice_size_t`  

Get y dimension of the Lattice.  
";

%feature("docstring") lm::rdme::Lattice::getZSize "
`getZSize() -> lattice_size_t`  

Get z dimension of the Lattice.  
";

%feature("docstring") lm::rdme::Lattice::getNumberSites "
`getNumberSites() -> lattice_size_t`  

Get total number of sites in the Lattice.  
";

%feature("docstring") lm::rdme::Lattice::getSpacing "
`getSpacing() -> si_dist_t`  

Get spacing between lattice sites.  
";

%feature("docstring") lm::rdme::Lattice::getNeighboringSites "
`getNeighboringSites(index, neighboringIndices)`  

Get the sites that are neighbor to the indicated site.  

Parameters
----------
* `index` :  
    Index of the site for which to get neighbors  
* `neighboringIndicies` :  
    An array to hold the indicies of the neighbor sites  
";

%feature("docstring") lm::rdme::Lattice::getSiteType "
`getSiteType(x, y, z) -> site_t`  

Get the site type at the specified location.  
";

%feature("docstring") lm::rdme::Lattice::getSiteType "
`getSiteType(index) -> site_t`  

Get the site type at the specified location.  
";

%feature("docstring") lm::rdme::Lattice::setSiteType "
`setSiteType(x, y, z, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::Lattice::setSiteType "
`setSiteType(index, site)`  

Set the site type at the specified location.  
";

%feature("docstring") lm::rdme::Lattice::getNearbySites "
`getNearbySites(xc, yc, zc, minDistance, maxDistance) -> std::vector< lattice_coord_t >`  

Get a list of sites near the specified site within a certain distance.  

Parameters
----------
* `xc` :  
    X location of the center site  
* `yc` :  
    Y location of the center site  
* `zc` :  
    Z location of the center site  
* `minDistance` :  
    Minimum distance for halo of sites  
* `maxDistance` :  
    Maximum distance for halo of sites  
";

%feature("docstring") lm::rdme::Lattice::getOccupancy "
`getOccupancy(x, y, z) -> site_size_t`  

Get the number of particles in the specified lattice site.  
";

%feature("docstring") lm::rdme::Lattice::getOccupancy "
`getOccupancy(index) -> site_size_t`  

Get the number of particles in the specified lattice site.  
";

%feature("docstring") lm::rdme::Lattice::getParticle "
`getParticle(x, y, z, particleIndex) -> particle_t`  

Get the particle at the specified site with at the specified number in the particle list.  
";

%feature("docstring") lm::rdme::Lattice::getParticle "
`getParticle(index, particleIndex) -> particle_t`  

Get the particle at the specified site with at the specified number in the particle list.  
";

%feature("docstring") lm::rdme::Lattice::addParticle "
`addParticle(x, y, z, particle)`  

Add a particle to the specified site.  
";

%feature("docstring") lm::rdme::Lattice::addParticle "
`addParticle(index, particle)`  

Add a particle to the specified site.  
";

%feature("docstring") lm::rdme::Lattice::removeParticles "
`removeParticles(x, y, z)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::Lattice::removeParticles "
`removeParticles(index)`  

Remove a particle to the specified site.  
";

%feature("docstring") lm::rdme::Lattice::removeAllParticles "
`removeAllParticles()`  

Empty all particles from the specified site.  
";

%feature("docstring") lm::rdme::Lattice::getParticleCounts "
`getParticleCounts() -> std::map< particle_t, uint >`  

Get the number of each particle type in the lattice.  

Particle searching methods.  
";

%feature("docstring") lm::rdme::Lattice::findParticles "
`findParticles(minParticleType, maxParticleType) -> std::vector< particle_loc_t >`  

Get the number of the specified particles types in the lattice.  
";

%feature("docstring") lm::rdme::Lattice::print "
`print()`  

Print the lattice to the console.  
";

%feature("docstring") lm::rdme::Lattice::setFromRowMajorByteData "
`setFromRowMajorByteData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::Lattice::setSitesFromRowMajorByteData "
`setSitesFromRowMajorByteData(buffer, bufferSize)`  
";

%feature("docstring") lm::rdme::Lattice::getLatticeMemorySize "
`getLatticeMemorySize() -> size_t`  
";

// File: structlattice__coord__t.xml


%feature("docstring") lattice_coord_t "

Type to store a lattice coordinate.  

C++ includes: Lattice.h
";

%feature("docstring") lattice_coord_t::lattice_coord_t "
`lattice_coord_t(x=0, y=0, z=0)`  

Create a lattice coordinate.  

Parameters
----------
* `x` :  
    Lattice x point  
* `y` :  
    Lattice y point  
* `z` :  
    Lattice z point  
";

// File: classlm_1_1builder_1_1LatticeBuilder.xml


%feature("docstring") lm::builder::LatticeBuilder "

A class that defines regions of a lattice based on a set of geometries defined by shapes. It also allows packing different types of particles into different regions.  

C++ includes: LatticeBuilder.h
";

%feature("docstring") lm::builder::LatticeBuilder::LatticeBuilder "
`LatticeBuilder(xLen, yLen, zLen, collisionGridSpacing, seedTop, seedBottom=0)`  

Create a Lattice Builder.  

Parameters
----------
* `xLen` :  
    Length of the domain along x-axis  
* `yLen` :  
    Length of the domain along y-axis  
* `zLen` :  
    Length of the domain along z-axis  
* `collisionGridSpacing` :  
    The spacing for collision objects  
* `seedTop` :  
    High 32 bits of the seed (allows a constant seed for debugging)  
* `seedBottom` :  
    Low 32 bits of the seed  
";

%feature("docstring") lm::builder::LatticeBuilder::~LatticeBuilder "
`~LatticeBuilder()`  
";

%feature("docstring") lm::builder::LatticeBuilder::addRegion "
`addRegion(shape)`  

Add a region to the lattice.  

Parameters
----------
* `shape` :  
    A Shape object to add as a region  
";

%feature("docstring") lm::builder::LatticeBuilder::placeObject "
`placeObject(shape) -> bool`  

Add an shape to the lattice.  

Parameters
----------
* `shape` :  
    A Shape object to add as a region  

Returns
-------
true if the object to place does not intersect another object  
";

%feature("docstring") lm::builder::LatticeBuilder::removeObject "
`removeObject(s)`  

Remove the shape from the lattice.  

s Shape to remove  
";

%feature("docstring") lm::builder::LatticeBuilder::placeSphere "
`placeSphere(center, radius, type) -> bool`  

Place a sphere in the lattice (for obstacles)  

Parameters
----------
* `center` :  
    The center point of sphere obstacle  
* `radius` :  
    Radius of the sphere obstacle  
* `type` :  
    The type of site in which to place sphere  

Returns
-------
true if the sphere did not intersect  
";

%feature("docstring") lm::builder::LatticeBuilder::removeSphere "
`removeSphere(center, radius, type)`  

Remove a sphere in the lattice (for obstacles)  

Parameters
----------
* `center` :  
    The center point of sphere obstacle  
* `radius` :  
    Radius of the sphere obstacle  
* `type` :  
    The type of site in which to place sphere  
";

%feature("docstring") lm::builder::LatticeBuilder::placeRandomSphere "
`placeRandomSphere(radius, type, region) -> uint`  

Place a sphere randomly in the lattice (for obstacles)  

Parameters
----------
* `radius` :  
    Radius of the sphere obstacle  
* `type` :  
    The type of site in which to place sphere  
* `region` :  
    The region in which to place obstacle randomly  

Returns
-------
number of times a placement operation occured  
";

%feature("docstring") lm::builder::LatticeBuilder::placeRandomSpheres "
`placeRandomSpheres(count, radius, type, region)`  

Place many spheres randomly in the lattice (for obstacles)  

Parameters
----------
* `radius` :  
    Radius of the sphere obstacle  
* `type` :  
    The type of site in which to place sphere  
* `region` :  
    The region in which to place obstacle randomly  

Returns
-------
number of times a placement operation occured  
";

%feature("docstring") lm::builder::LatticeBuilder::fillWithRandomSpheres "
`fillWithRandomSpheres(volumeFraction, radius, type, region)`  

Fill a region with random spheres to a specified volume fraction.  

Parameters
----------
* `volumeFraction` :  
    Total fraction of volume that should be filled with spheres  
* `radius` :  
    Radius of spheres  
* `type` :  
    The type of site to fill (i.e. the type of site to exclude other objects from)  
* `region` :  
    The region of the lattice to place spheres into  
";

%feature("docstring") lm::builder::LatticeBuilder::getSpatialModel "
`getSpatialModel(spatialModel)`  

Gets a spatial model of the lattice for interface with python. NOTE: this operation clears the object passed in from python.  

Parameters
----------
* `spatialModel` :  
    An object of a spatial model for interaction in python or HDF5. The model will be filled with the current lattice  
";

%feature("docstring") lm::builder::LatticeBuilder::addParticles "
`addParticles(particleType, siteType, count)`  

Add particles of a given type.  

Parameters
----------
* `particleType` :  
    The type of particles to randomly place in the lattice  
* `siteType` :  
    Type of lattice site into which to place  
* `count` :  
    Number of particles to place  
";

%feature("docstring") lm::builder::LatticeBuilder::discretizeTo "
`discretizeTo(lattice, obstacleSiteType, fractionObstacleSitesOccupied)`  

Discretizes the regions to a square lattice.  

Parameters
----------
* `lattice` :  
    Lattice object into which to place  
* `obstacleSiteType` :  
    An identifier for obstacle sites in the lattice  
* `fractionObstacleSitesOccupied` :  
    Percentage of obstacle sites to be filled  
";

// File: classlm_1_1main_1_1LocalDataOutputWorker.xml


%feature("docstring") lm::main::LocalDataOutputWorker "

A class that handles output on a particular core.  

C++ includes: LocalDataOutputWorker.h
";

%feature("docstring") lm::main::LocalDataOutputWorker::LocalDataOutputWorker "
`LocalDataOutputWorker(file)`  

Create a new LocalDataOuputWorker.  

Parameters
----------
* `file` :  
    The file in which to save the current data stream  
";

%feature("docstring") lm::main::LocalDataOutputWorker::~LocalDataOutputWorker "
`~LocalDataOutputWorker()`  
";

%feature("docstring") lm::main::LocalDataOutputWorker::pushDataSet "
`pushDataSet(dataSet)`  

Push a DataSet for writing.  

Parameters
----------
* `dataSet` :  
    DataSet to add to the queue for writing  
";

%feature("docstring") lm::main::LocalDataOutputWorker::wake "
`wake()`  

Wake the thread from sleeping state.  
";

%feature("docstring") lm::main::LocalDataOutputWorker::abort "
`abort()`  

Abort the thread.  
";

%feature("docstring") lm::main::LocalDataOutputWorker::checkpoint "
`checkpoint()`  

Tell the thread to checkpoint at the next available times.  
";

// File: structlm_1_1builder_1_1matrix.xml


%feature("docstring") lm::builder::matrix "

A matrix used for rotations.  

C++ includes: Shape.h
";

%feature("docstring") lm::builder::matrix::matrix "
`matrix(m11=0.0, m12=0.0, m13=0.0, m21=0.0, m22=0.0, m23=0.0, m31=0.0, m32=0.0, m33=0.0)`  

Create a matrix with the speficied elements.  
";

%feature("docstring") lm::builder::matrix::transpose "
`transpose() -> matrix`  
";

%feature("docstring") lm::builder::matrix::determinant "
`determinant() -> si_dist_t`  
";

%feature("docstring") lm::builder::matrix::trace "
`trace() -> si_dist_t`  
";

%feature("docstring") lm::builder::matrix::mult "
`mult(r) -> vector`  

Multiply by a vector.  
";

%feature("docstring") lm::builder::matrix::mult "
`mult(r) -> matrix`  

Multiply by a matrix.  
";

%feature("docstring") lm::builder::matrix::eulerMatrixFromAngles "
`eulerMatrixFromAngles(phi, theta, psi) -> matrix`  

Get a forward rotation matrix from angles.  

Parameters
----------
* `phi` :  
    The angle around x (in radians)  
* `theta` :  
    The angle around y (in radians)  
* `psi` :  
    The angle around z (in radians)  
";

%feature("docstring") lm::builder::matrix::Identity "
`Identity() -> matrix`  

Get an identity matrix.  
";

// File: classlm_1_1builder_1_1Mesh.xml


%feature("docstring") lm::builder::Mesh "

A Shape that represents a mesh.  

C++ includes: Mesh.h
";

%feature("docstring") lm::builder::Mesh::Mesh "
`Mesh(center, radius, type)`  

Create a Mesh.  

Parameters
----------
* `center` :  
    Point center of the circle of the slice plane through the sphere  
* `radius` :  
    Radius of the sphere  
* `type` :  
    The type of the sites within the sphere  
";

%feature("docstring") lm::builder::Mesh::~Sphere "
`~Sphere()`  

Destroy the Sphere.  
";

%feature("docstring") lm::builder::Mesh::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Mesh::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Mesh::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Mesh::setCenter "
`setCenter(center)`  

Set the center of the sphere.  

Parameters
----------
* `center` :  
    Point of the center  
";

%feature("docstring") lm::builder::Mesh::getCenter "
`getCenter() -> point`  

Get the center of the sphere.  
";

%feature("docstring") lm::builder::Mesh::getRadius "
`getRadius() -> si_dist_t`  

Get the radius of the sphere.  
";

%feature("docstring") lm::builder::Mesh::getVolume "
`getVolume() -> double`  

Get the volume bounded by the sphere.  
";

// File: classlm_1_1me_1_1MESolver.xml


%feature("docstring") lm::me::MESolver "

An abstract base class for all Master Equation solvers, this is essentially a representation of \"the simulation instance\".  

C++ includes: MESolver.h
";

%feature("docstring") lm::me::MESolver::MESolver "
`MESolver()`  

Create the MESolver.  
";

%feature("docstring") lm::me::MESolver::~MESolver "
`~MESolver()`  
";

%feature("docstring") lm::me::MESolver::initialize "
`initialize(replicate, parameters, resources)`  

Initialize the simulation.  

Parameters
----------
* `replicate` :  
    Replicate number out of total replicates  
* `parameters` :  
    A map of all the parameters for the simulation  
* `A` :  
    list of resources assigned to the simulation  
";

%feature("docstring") lm::me::MESolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::me::MESolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::me::MESolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1me_1_1MESolverFactory.xml


%feature("docstring") lm::me::MESolverFactory "

A factory object used to create Master Equation simulation instances.  

C++ includes: MESolverFactory.h
";

%feature("docstring") lm::me::MESolverFactory::MESolverFactory "
`MESolverFactory()`  

Create the factory.  
";

%feature("docstring") lm::me::MESolverFactory::setSolver "
`setSolver(solver)`  

Set the type of solver the factory makes.  

Parameters
----------
* `solver` :  
    The name of the solver from the command line  
";

%feature("docstring") lm::me::MESolverFactory::needsReactionModel "
`needsReactionModel() -> bool`  

Tell if solver needs a reaction model.  
";

%feature("docstring") lm::me::MESolverFactory::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tell if solver needs a diffusion model.  
";

%feature("docstring") lm::me::MESolverFactory::instantiate "
`instantiate() -> MESolver *`  

Get an instance of a MESolver object for running a replicate.  

Returns
-------
The new solver instance  
";

// File: classlm_1_1rdme_1_1MGPUMpdRdmeSolver.xml


%feature("docstring") lm::rdme::MGPUMpdRdmeSolver "

C++ includes: MGPUMpdRdmeSolver.h
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::MGPUMpdRdmeSolver "
`MGPUMpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::~MGPUMpdRdmeSolver "
`~MGPUMpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::initialize "
`initialize(replicate, parameters, resources)`  

Initialize the simulation.  

Parameters
----------
* `replicate` :  
    Replicate number out of total replicates  
* `parameters` :  
    A map of all the parameters for the simulation  
* `A` :  
    list of resources assigned to the simulation  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::buildModel "
`buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypeA, kA, SA, DA, kCols=1)`  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::buildDiffusionModel "
`buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData=true)`  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

%feature("docstring") lm::rdme::MGPUMpdRdmeSolver::setReactionRate "
`setReactionRate(rxid, rate)`  
";

// File: structlm_1_1cme_1_1CMESolver_1_1MichaelisMentenPropensityArgs.xml

// File: classlm_1_1rdme_1_1MpdRdmeSolver.xml


%feature("docstring") lm::rdme::MpdRdmeSolver "

C++ includes: MpdRdmeSolver.h
";

%feature("docstring") lm::rdme::MpdRdmeSolver::MpdRdmeSolver "
`MpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver::~MpdRdmeSolver "
`~MpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver::initialize "
`initialize(replicate, parameters, resources)`  

Initialize the simulation.  

Parameters
----------
* `replicate` :  
    Replicate number out of total replicates  
* `parameters` :  
    A map of all the parameters for the simulation  
* `A` :  
    list of resources assigned to the simulation  
";

%feature("docstring") lm::rdme::MpdRdmeSolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::MpdRdmeSolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::MpdRdmeSolver::buildModel "
`buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypeA, kA, SA, DA, kCols=1)`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver::buildDiffusionModel "
`buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData=true)`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1rdme_1_1MpdTestHarness.xml


%feature("docstring") lm::rdme::MpdTestHarness "

C++ includes: MpdTestHarness.h
";

%feature("docstring") lm::rdme::MpdTestHarness::MpdTestHarness "
`MpdTestHarness()`  
";

%feature("docstring") lm::rdme::MpdTestHarness::~MpdTestHarness "
`~MpdTestHarness()`  
";

%feature("docstring") lm::rdme::MpdTestHarness::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1MPI.xml


%feature("docstring") lm::MPI "

Handles the MPI capabilities and properties of the simulation.  

C++ includes: lm_mpi.h
";

%feature("docstring") lm::MPI::init "
`init(argc, argv)`  

Initialize the MPI runtime.  

Parameters
----------
* `argc` :  
    Command line number of arguments  
* `argv` :  
    Command line array of arguments  
";

%feature("docstring") lm::MPI::printCapabilities "
`printCapabilities()`  

Print the capabilities of the current processing device.  
";

%feature("docstring") lm::MPI::finalize "
`finalize()`  

Close and cleanup the MPI runtime.  
";

// File: classlm_1_1MPIException.xml


%feature("docstring") lm::MPIException "

An Exception thrown when an MPI call fails.  

C++ includes: lm_mpi.h
";

%feature("docstring") lm::MPIException::MPIException "
`MPIException(error)`  

Create the MPIException.  

Parameters
----------
* `error` :  
    MPI error code  
";

// File: classlm_1_1rdme_1_1MPIMpdRdmeSolver.xml


%feature("docstring") lm::rdme::MPIMpdRdmeSolver "

C++ includes: MPIMpdRdmeSolver.h
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::MPIMpdRdmeSolver "
`MPIMpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::~MPIMpdRdmeSolver "
`~MPIMpdRdmeSolver()`  
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::initialize "
`initialize(replicate, parameters, resources)`  

Initialize the simulation.  

Parameters
----------
* `replicate` :  
    Replicate number out of total replicates  
* `parameters` :  
    A map of all the parameters for the simulation  
* `A` :  
    list of resources assigned to the simulation  
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::buildModel "
`buildModel(numberSpeciesA, numberReactionsA, initialSpeciesCountsA, reactionTypeA, kA, SA, DA, kCols=1)`  
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::buildDiffusionModel "
`buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData=true)`  
";

%feature("docstring") lm::rdme::MPIMpdRdmeSolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1main_1_1MPIRemoteDataOutputQueue.xml


%feature("docstring") lm::main::MPIRemoteDataOutputQueue "

A DataOutputQueue that works for MPI jobs across various different nodes that outputs only on the \"master\" thread.  

C++ includes: MPIRemoteDataOutputQueue.h
";

%feature("docstring") lm::main::MPIRemoteDataOutputQueue::MPIRemoteDataOutputQueue "
`MPIRemoteDataOutputQueue()`  

Create a new MPIRemoteDataOutputQueue.  
";

%feature("docstring") lm::main::MPIRemoteDataOutputQueue::~MPIRemoteDataOutputQueue "
`~MPIRemoteDataOutputQueue()`  
";

%feature("docstring") lm::main::MPIRemoteDataOutputQueue::isEmpty "
`isEmpty() -> bool`  

Tell if the queue is empty.  

Returns
-------
true if the queue is empty  
";

%feature("docstring") lm::main::MPIRemoteDataOutputQueue::popDataSetIntoBuffer "
`popDataSetIntoBuffer(buffer, bufferSize) -> size_t`  

Pop data from queue and write to the buffer.  

Parameters
----------
* `buffer` :  
    The data buffer into which the next DataSet object should be written  
* `bufferSize` :  
    The size of the buffer  

Returns
-------
Number of bytes written to the buffer  
";

// File: classMultiGPUMapper.xml


%feature("docstring") MultiGPUMapper "

C++ includes: MultiGPUMapper.h
";

%feature("docstring") MultiGPUMapper::MultiGPUMapper "
`MultiGPUMapper(ldim, cellsize, apron, overlap, num_gpus, devices, pages)`  
";

%feature("docstring") MultiGPUMapper::~MultiGPUMapper "
`~MultiGPUMapper()`  
";

%feature("docstring") MultiGPUMapper::get_num_gpus "
`get_num_gpus() -> int`  
";

%feature("docstring") MultiGPUMapper::use "
`use(gpu) -> bool`  
";

%feature("docstring") MultiGPUMapper::get_overlap "
`get_overlap() -> int`  
";

%feature("docstring") MultiGPUMapper::get_apron "
`get_apron() -> int`  
";

%feature("docstring") MultiGPUMapper::set_affinity "
`set_affinity(arg1)`  
";

%feature("docstring") MultiGPUMapper::get_affinity "
`get_affinity() -> int`  
";

%feature("docstring") MultiGPUMapper::get_lattice_dim "
`get_lattice_dim() -> dim3`  
";

%feature("docstring") MultiGPUMapper::getSegmentDescriptor "
`getSegmentDescriptor(gpu) -> SegmentDescriptor_s *`  
";

%feature("docstring") MultiGPUMapper::get_global_size "
`get_global_size() -> size_t`  
";

%feature("docstring") MultiGPUMapper::record_execution_cost "
`record_execution_cost(arg1, arg2)`  
";

%feature("docstring") MultiGPUMapper::rebalance "
`rebalance() -> bool`  
";

%feature("docstring") MultiGPUMapper::numa_bind_thread "
`numa_bind_thread(arg1) -> bool`  
";

%feature("docstring") MultiGPUMapper::get_global_dim "
`get_global_dim(gpu) -> dim3`  
";

%feature("docstring") MultiGPUMapper::get_local_dim "
`get_local_dim(gpu) -> dim3`  
";

%feature("docstring") MultiGPUMapper::get_global_offset "
`get_global_offset(gpu) -> int3`  
";

%feature("docstring") MultiGPUMapper::get_local_size "
`get_local_size(gpu) -> size_t`  
";

%feature("docstring") MultiGPUMapper::get_authority_size "
`get_authority_size(gpu) -> size_t`  
";

%feature("docstring") MultiGPUMapper::get_global_input_offset "
`get_global_input_offset(gpu) -> ssize_t`  
";

%feature("docstring") MultiGPUMapper::get_global_output_offset "
`get_global_output_offset(gpu) -> size_t`  
";

%feature("docstring") MultiGPUMapper::get_authority_offset "
`get_authority_offset(gpu) -> size_t`  
";

%feature("docstring") MultiGPUMapper::stage_in "
`stage_in(gpu, dptr, hptr)`  
";

%feature("docstring") MultiGPUMapper::stage_in_sites "
`stage_in_sites(gpu, dptr, hptr)`  
";

%feature("docstring") MultiGPUMapper::stage_out "
`stage_out(gpu, hptr, dptr)`  
";

%feature("docstring") MultiGPUMapper::publish_state "
`publish_state(gpu, dptr, timestamp)`  
";

%feature("docstring") MultiGPUMapper::refresh "
`refresh(gpu, dptr, timestamp)`  
";

%feature("docstring") MultiGPUMapper::schedule_send "
`schedule_send(gpu, dptr, timestamp, neighbor, stream)`  
";

%feature("docstring") MultiGPUMapper::schedule_recv "
`schedule_recv(gpu, dptr, timestamp, neighbor, stream)`  
";

%feature("docstring") MultiGPUMapper::map_index_to_gpu "
`map_index_to_gpu(index) -> int`  
";

%feature("docstring") MultiGPUMapper::initialize_gpu "
`initialize_gpu(gpu)`  
";

%feature("docstring") MultiGPUMapper::determine_load_balance "
`determine_load_balance() -> bool`  
";

// File: structneighbor__buffer.xml


%feature("docstring") neighbor_buffer "

C++ includes: ZDivMPIGPUMapper.h
";

// File: classlm_1_1cme_1_1NextReactionSolver.xml


%feature("docstring") lm::cme::NextReactionSolver "

C++ includes: NextReactionSolver.h
";

%feature("docstring") lm::cme::NextReactionSolver::NextReactionSolver "
`NextReactionSolver()`  
";

%feature("docstring") lm::cme::NextReactionSolver::NextReactionSolver "
`NextReactionSolver(neededDists)`  
";

%feature("docstring") lm::cme::NextReactionSolver::~NextReactionSolver "
`~NextReactionSolver()`  
";

%feature("docstring") lm::cme::NextReactionSolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::NextReactionSolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::NextReactionSolver::buildModel "
`buildModel(numberSpecies, numberReactions, initialSpeciesCounts, reactionType, k, S, D, kCols=1)`  
";

%feature("docstring") lm::cme::NextReactionSolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1rdme_1_1NextSubvolumeSolver.xml


%feature("docstring") lm::rdme::NextSubvolumeSolver "

C++ includes: NextSubvolumeSolver.h
";

%feature("docstring") lm::rdme::NextSubvolumeSolver::NextSubvolumeSolver "
`NextSubvolumeSolver()`  
";

%feature("docstring") lm::rdme::NextSubvolumeSolver::NextSubvolumeSolver "
`NextSubvolumeSolver(neededDists)`  
";

%feature("docstring") lm::rdme::NextSubvolumeSolver::~NextSubvolumeSolver "
`~NextSubvolumeSolver()`  
";

%feature("docstring") lm::rdme::NextSubvolumeSolver::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::NextSubvolumeSolver::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::rdme::NextSubvolumeSolver::buildDiffusionModel "
`buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData=true)`  
";

%feature("docstring") lm::rdme::NextSubvolumeSolver::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: structlm_1_1cme_1_1CMESolver_1_1NoncompetitiveMMPropensityArgs.xml

// File: structlm_1_1cme_1_1SelfRegulatingGeneSwitch_1_1OUKHillPropensityArgs.xml

// File: structlm_1_1cme_1_1FluctuatingNRSolver_1_1OUPropensityArgs.xml

// File: structparticle__loc__t.xml


%feature("docstring") particle_loc_t "

Type to store a particle and position.  

C++ includes: Lattice.h
";

%feature("docstring") particle_loc_t::particle_loc_t "
`particle_loc_t(p=0, x=0, y=0, z=0, index=0)`  

Create a particle location.  

Parameters
----------
* `p` :  
    Particle type at the location  
* `x` :  
    Lattice x point  
* `y` :  
    Lattice y point  
* `z` :  
    Lattice z point  
* `index` :  
    Lattice index  
";

// File: structlm_1_1builder_1_1LatticeBuilder_1_1ParticlePlacement.xml

// File: structlm_1_1builder_1_1point.xml


%feature("docstring") lm::builder::point "

Type to store a position in space.  

C++ includes: Shape.h
";

%feature("docstring") lm::builder::point::point "
`point(x=0.0, y=0.0, z=0.0)`  

Create a point.  

Parameters
----------
* `x` :  
    X location  
* `y` :  
    Y location  
* `z` :  
    Z location  
";

%feature("docstring") lm::builder::point::point "
`point(p2)`  

Create a point from another point.  

Parameters
----------
* `p2` :  
    The point to copy  
";

%feature("docstring") lm::builder::point::distanceSquared "
`distanceSquared(p2) -> si_dist_t`  

Determine the distance squared between two points.  

Parameters
----------
* `p2` :  
    Second point  

Returns
-------
scalar distance squared  
";

%feature("docstring") lm::builder::point::minus "
`minus(r) -> point`  

Subtract the vector from the point.  

Parameters
----------
* `v` :  
    The vector to subtract  

Returns
-------
A new point at that location  
";

%feature("docstring") lm::builder::point::plus "
`plus(r) -> point`  

Add the vector from the point.  

Parameters
----------
* `v` :  
    The vector to add  

Returns
-------
A new point at that location  
";

%feature("docstring") lm::builder::point::distance "
`distance(p2) -> si_dist_t`  

Determine the distance to another point.  

Parameters
----------
* `p2` :  
    Second point  

Returns
-------
scalar distance  
";

// File: classlm_1_1Print.xml


%feature("docstring") lm::Print "

Print messages to the console at varying levels of verbosity.  

C++ includes: Print.h
";

%feature("docstring") lm::Print::printf "
`printf(level, fmt, arg3)`  

Prints to the console at varying levels of verbosity or \"Severity\".  

Parameters
----------
* `level` :  
    The level of \"Severity\" as indicated above  
* `fmt` :  
    A printf style format string  
* `Arguments` :  
    to the printf format string  
";

%feature("docstring") lm::Print::printDateTimeString "
`printDateTimeString()`  

Print the date and time to the console.  
";

%feature("docstring") lm::Print::getDateTimeString "
`getDateTimeString() -> std::string`  

Gets the date and time as a string.  

Returns
-------
Date and time as a string  
";

%feature("docstring") lm::Print::verbosityLevel "
`verbosityLevel() -> int`  

Gets the verbosity level. Lower numbers are higher priority.  

Returns
-------
Verbosity level  
";

%feature("docstring") lm::Print::verbosityLevel "
`verbosityLevel(x)`  

Sets the verbosity level. Lower numbers are higher priority.  

Parameters
----------
* `x` :  
    Verbosity level  
";

%feature("docstring") lm::Print::ifPrinted "
`ifPrinted(x) -> bool`  

True if a call to `Print::printf(x, ...)` would print to stdout.  

Parameters
----------
* `x` :  
    Verbosity level  

Returns
-------
Predicate  
";

// File: structlm_1_1cme_1_1CMESolver_1_1PropensityArgs.xml

// File: classlm_1_1thread_1_1PthreadException.xml


%feature("docstring") lm::thread::PthreadException "

An Exception class for handling pthread exceptions.  

C++ includes: Thread.h
";

%feature("docstring") lm::thread::PthreadException::PthreadException "
`PthreadException(error, file, line)`  
";

// File: classlm_1_1rng_1_1RandomGenerator.xml


%feature("docstring") lm::rng::RandomGenerator "

Base class for random number generators in Lattice Microbes.  

C++ includes: RandomGenerator.h
";

%feature("docstring") lm::rng::RandomGenerator::RandomGenerator "
`RandomGenerator(seedTop, seedBottom, availableDists=(Distributions)(ALL))`  
";

%feature("docstring") lm::rng::RandomGenerator::~RandomGenerator "
`~RandomGenerator()`  
";

%feature("docstring") lm::rng::RandomGenerator::getSeed "
`getSeed() -> uint64_t`  

Get the current seed.  

Returns
-------
seed  
";

%feature("docstring") lm::rng::RandomGenerator::getRandom "
`getRandom() -> uint32_t`  

Get a random integer.  

Returns
-------
random A random integer  
";

%feature("docstring") lm::rng::RandomGenerator::getRandomDouble "
`getRandomDouble() -> double`  

Get a random double.  

Returns
-------
random A random double  
";

%feature("docstring") lm::rng::RandomGenerator::getExpRandomDouble "
`getExpRandomDouble() -> double`  

Get a random exponentially distributed double.  

Returns
-------
random A random double  
";

%feature("docstring") lm::rng::RandomGenerator::getNormRandomDouble "
`getNormRandomDouble() -> double`  

Get a random normally distributed double.  

Returns
-------
random A random double  
";

%feature("docstring") lm::rng::RandomGenerator::getRandomDoubles "
`getRandomDoubles(rngs, numberRNGs)`  

Get a number of random doubles.  

Parameters
----------
* `rngs` :  
    A preallocated array for random numbers  
* `numberRNGs` :  
    Number or randoms to put in array  
";

%feature("docstring") lm::rng::RandomGenerator::getExpRandomDoubles "
`getExpRandomDoubles(rngs, numberRNGs)`  

Get a number of random exponentially distiributed doubles.  

Parameters
----------
* `rngs` :  
    A preallocated array for random numbers  
* `numberRNGs` :  
    Number or randoms to put in array  
";

%feature("docstring") lm::rng::RandomGenerator::getNormRandomDoubles "
`getNormRandomDoubles(rngs, numberRNGs)`  

Get a number of random normally distributed doubles.  

Parameters
----------
* `rngs` :  
    A preallocated array for random numbers  
* `numberRNGs` :  
    Number or randoms to put in array  
";

// File: classlm_1_1rdme_1_1RDMESolver.xml


%feature("docstring") lm::rdme::RDMESolver "

C++ includes: RDMESolver.h
";

%feature("docstring") lm::rdme::RDMESolver::RDMESolver "
`RDMESolver(neededDists)`  
";

%feature("docstring") lm::rdme::RDMESolver::~RDMESolver "
`~RDMESolver()`  
";

%feature("docstring") lm::rdme::RDMESolver::setDiffusionModel "
`setDiffusionModel(dm, lattice, latticeSize, latticeSites, latticeSitesSize)`  
";

%feature("docstring") lm::rdme::RDMESolver::buildDiffusionModel "
`buildDiffusionModel(numberSiteTypesA, DFA, RLA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing, latticeData, latticeSitesData, rowMajorData=true)`  
";

// File: structlm_1_1reaction_1_1ReactionQueue_1_1ReactionEvent.xml


%feature("docstring") lm::reaction::ReactionQueue::ReactionEvent "

Definition of a reaction event with the time to reaction and the propensity for the reaction to occur.  

C++ includes: ReactionQueue.h
";

// File: classlm_1_1reaction_1_1ReactionQueue.xml


%feature("docstring") lm::reaction::ReactionQueue "

A queue that contains information on reaction events.  

C++ includes: ReactionQueue.h
";

%feature("docstring") lm::reaction::ReactionQueue::ReactionQueue "
`ReactionQueue(numberReactions)`  

Initialize the reaction queue.  

Parameters
----------
* `numberReactions` :  
    The number of reactions that the queue should be able to hold  
";

%feature("docstring") lm::reaction::ReactionQueue::~ReactionQueue "
`~ReactionQueue()`  

Destroy the reaction queue.  
";

%feature("docstring") lm::reaction::ReactionQueue::getNextReaction "
`getNextReaction() -> uint`  

get the next reaction  

Returns
-------
reaction Number of the next reaction  
";

%feature("docstring") lm::reaction::ReactionQueue::getReactionEvent "
`getReactionEvent(reactionIndex) -> ReactionEvent`  

Get the next reaction event.  

Parameters
----------
* `reactionIndex` :  
    The index of the reaction event reactionEvent The specified reaction event  
";

%feature("docstring") lm::reaction::ReactionQueue::updateReactionEvent "
`updateReactionEvent(reactionIndex, newTime, newPropensity)`  

Update the reaction event queue.  

Parameters
----------
* `reactionIndex` :  
    Index of the reaction to update  
* `newTime` :  
    The new time at which the reaction will occur  
* `newPropensity` :  
    The new propensity for the reaction to occur  
";

// File: structlm_1_1io_1_1hdf5_1_1SimulationFile_1_1ReplicateHandles.xml


%feature("docstring") lm::io::hdf5::SimulationFile::ReplicateHandles "

A handle for the different replicates that may be stored in the HDF5 file.  

C++ includes: SimulationFile.h
";

%feature("docstring") lm::io::hdf5::SimulationFile::ReplicateHandles::ReplicateHandles "
`ReplicateHandles()`  

Set up the ReplicateHandles.  
";

// File: classlm_1_1main_1_1ReplicateRunner.xml


%feature("docstring") lm::main::ReplicateRunner "

A thread that launches all the various replicates requested.  

C++ includes: ReplicateRunner.h
";

%feature("docstring") lm::main::ReplicateRunner::ReplicateRunner "
`ReplicateRunner(replicate, solverFactory, parameters, reactionModel, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resources)`  

Create a new replicate runner.  

Parameters
----------
* `replicate` :  
    The number of the first replicate  
* `solverFactory` :  
    The factory used to create a solver object for each replicate  
* `parameters` :  
    A map of parameters for the simulation  
* `reactionModel` :  
    An object to take care of reactions  
* `diffusionModel` :  
    An object to take care of diffusion  
* `lattice` :  
    The actual data representing the lattice  
* `latticeSize` :  
    The size of the lattice in bytes  
* `latticeSitesSize` :  
    The size of the lattice sites in bytes  
* `resources` :  
    A manager for the resources that handles GPUs and CPUs  
";

%feature("docstring") lm::main::ReplicateRunner::~ReplicateRunner "
`~ReplicateRunner()`  
";

%feature("docstring") lm::main::ReplicateRunner::wake "
`wake()`  

Wake the thread from sleep.  
";

%feature("docstring") lm::main::ReplicateRunner::run "
`run() -> int`  

Run the thread.  

Returns
-------
0 on success, -1 on failure  
";

%feature("docstring") lm::main::ReplicateRunner::getReplicate "
`getReplicate() -> int`  

Get the current replicate number.  
";

%feature("docstring") lm::main::ReplicateRunner::hasReplicateFinished "
`hasReplicateFinished() -> bool`  

Tell whether the previously running replicate has finished running.  
";

%feature("docstring") lm::main::ReplicateRunner::getReplicateExitCode "
`getReplicateExitCode() -> int`  

Get the exit code for the last replicate to finish.  
";

// File: classlm_1_1main_1_1ResourceAllocator.xml


%feature("docstring") lm::main::ResourceAllocator "

An object that tracks the available resources for the main simulation runner.  

C++ includes: ResourceAllocator.h
";

%feature("docstring") lm::main::ResourceAllocator::ResourceAllocator "
`ResourceAllocator(processNumber, numberCpuCores, cpuCoresPerReplicate)`  

Create a ResourceAllocator.  

Parameters
----------
* `processNumber` :  
    The process identifier  
* `numberCpuCores` :  
    The number of cores on teh resource  
* `cpuCoresPerPrelicate` :  
    The number of cores to be used by each simulation replicate  
";

%feature("docstring") lm::main::ResourceAllocator::ResourceAllocator "
`ResourceAllocator(processNumber, numberCpuCores, cpuCoresPerReplicate, cudaDevices, cudaDevicesPerReplicate)`  

Create a ResourceAllocator.  

Parameters
----------
* `processNumber` :  
    The process identifier  
* `numberCpuCores` :  
    The number of cores on teh resource  
* `cpuCoresPerPrelicate` :  
    The number of cores to be used by each simulation replicate  
* `cudaDevices` :  
    A set of the identifiers for the available CUDA devices  
* `cudaDevicesPerReplicate` :  
    The number of GPUs to be used by each simulation replicate  
";

%feature("docstring") lm::main::ResourceAllocator::~ResourceAllocator "
`~ResourceAllocator()`  
";

%feature("docstring") lm::main::ResourceAllocator::getMaxSimultaneousReplicates "
`getMaxSimultaneousReplicates() -> int`  

Get maximum number of replicates that can run at a time on the available resources based on cudaDevices and numberCpuCores.  
";

%feature("docstring") lm::main::ResourceAllocator::assignReplicate "
`assignReplicate(replicate) -> ComputeResources`  

Assign a replicate to free resources.  

Parameters
----------
* `replicate` :  
    The replicate identifier  

Returns
-------
A class representing the compute resources given to the replicate  
";

%feature("docstring") lm::main::ResourceAllocator::removeReplicate "
`removeReplicate(replicate)`  

Remove a replicate from those that are running.  

Parameters
----------
* `replicate` :  
    The replicate identifier  
";

%feature("docstring") lm::main::ResourceAllocator::reserveCpuCore "
`reserveCpuCore() -> int`  

Reserve a particular core.  

Returns
-------
The identifer fo the core that has been reserved  
";

// File: structlm_1_1cme_1_1CMESolver_1_1SecondOrderPropensityArgs.xml

// File: structlm_1_1cme_1_1CMESolver_1_1SecondOrderSelfPropensityArgs.xml

// File: structsegmentDescriptor.xml


%feature("docstring") segmentDescriptor "

C++ includes: SegmentDescriptor.h
";

// File: classlm_1_1cme_1_1SelfRegulatingGeneSwitch.xml


%feature("docstring") lm::cme::SelfRegulatingGeneSwitch "

C++ includes: SelfRegulatingGeneSwitch.h
";

%feature("docstring") lm::cme::SelfRegulatingGeneSwitch::SelfRegulatingGeneSwitch "
`SelfRegulatingGeneSwitch()`  
";

%feature("docstring") lm::cme::SelfRegulatingGeneSwitch::~SelfRegulatingGeneSwitch "
`~SelfRegulatingGeneSwitch()`  
";

%feature("docstring") lm::cme::SelfRegulatingGeneSwitch::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::SelfRegulatingGeneSwitch::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::SelfRegulatingGeneSwitch::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1builder_1_1Shape.xml


%feature("docstring") lm::builder::Shape "

Abstract base class for all the shapes in Lattice Microbes simulation builder.  

C++ includes: Shape.h
";

%feature("docstring") lm::builder::Shape::Shape "
`Shape(shapeType, boundingBox, type, at=vector(0.0, 0.0, 1.0), up=vector(0.0, 1.0, 0.0))`  

Create a Shape.  

Parameters
----------
* `shapeType` :  
    The type of shape should be of type ShapeType  
* `boundingBox` :  
    The extents of the object used for fast collision/contains checking  
* `type` :  
    The type of site that the object represents  
* `at` :  
    A vector describing what direction the object is pointing \"at\"; default: (0,0,1)  
* `up` :  
    A vector describing what direction is \"up\" relative to \"at\"; default: (0,1,0)  
";

%feature("docstring") lm::builder::Shape::~Shape "
`~Shape()`  

Destroy the Shape.  
";

%feature("docstring") lm::builder::Shape::boundingBoxesIntersect "
`boundingBoxesIntersect(query) -> bool`  

Checks if another shape's bounding box interstects with this shape's bounding box.  

Parameters
----------
* `query` :  
    The other shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Shape::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Shape::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Shape::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Shape::getBoundingBox "
`getBoundingBox() -> bounding_box`  

Get the bounding box.  
";

%feature("docstring") lm::builder::Shape::getType "
`getType() -> site_t`  

Get the site type associated with the shape.  
";

%feature("docstring") lm::builder::Shape::getShapeType "
`getShapeType() -> ShapeType`  

Get the shape type.  
";

%feature("docstring") lm::builder::Shape::getVolume "
`getVolume() -> double`  

Get the total internal volume of the shape.  
";

%feature("docstring") lm::builder::Shape::discretizeTo "
`discretizeTo(lattice)`  

Discretize the object to the specified lattice.  

Parameters
----------
* `lattice` :  
    Lattice in which to discretize the shape  
";

// File: classlm_1_1main_1_1SignalHandler.xml


%feature("docstring") lm::main::SignalHandler "

A thread safe Worker thread type that captures SIG* type signals over the various threads.  

C++ includes: SignalHandler.h
";

%feature("docstring") lm::main::SignalHandler::SignalHandler "
`SignalHandler()`  

Create a SignalHandler.  
";

%feature("docstring") lm::main::SignalHandler::~SignalHandler "
`~SignalHandler()`  
";

%feature("docstring") lm::main::SignalHandler::wake "
`wake()`  

Wake a sleeping Signalhandler.  
";

// File: classlm_1_1io_1_1hdf5_1_1SimulationFile.xml


%feature("docstring") lm::io::hdf5::SimulationFile "

A representation of the simulation that is used to input or output from an HDF5 file.  

C++ includes: SimulationFile.h
";

%feature("docstring") lm::io::hdf5::SimulationFile::isValidFile "
`isValidFile(filename) -> bool`  

Check that the HDF5 file is valid.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::isValidFile "
`isValidFile(filename) -> bool`  

Check that the HDF5 file is valid.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::create "
`create(filename)`  

Create an HDF5 file.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::create "
`create(filename)`  

Create an HDF5 file.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::create "
`create(filename, numberSpecies)`  

Create an HDF5 file with the specified number of species.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::create "
`create(filename, numberSpecies)`  

Create an HDF5 file with the specified number of species.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::create "
`create(filename, initializeModel, numberSpecies=0)`  

Create an HDF5 file and initialize the model with an optional number of species.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::SimulationFile "
`SimulationFile(filename)`  

Create a SimulationFile with the specified filename.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::SimulationFile "
`SimulationFile(filename)`  

Create a SimulationFile with the specified filename.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::~SimulationFile "
`~SimulationFile()`  

Destroy the SimulationFile.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::close "
`close()`  

Close the file.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::flush "
`flush()`  

Flush the current buffered data to the file.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::checkpoint "
`checkpoint() -> string`  

Create the current checkpoint to the HDF5 file.  

Parameters
----------
* `checkpointFilename` :  
    Name of the new checkpoint file  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getParameters "
`getParameters() -> map< string, string >`  

Get all the parameters from the current replicate.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getParameter "
`getParameter(key, defaultValue=\"\") -> string`  

Get the specified parameter from the current replicate.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setParameter "
`setParameter(key, value)`  

Set the specified parameter to the value.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getReactionModel "
`getReactionModel(reactionModel)`  

Pops the reaction model from the HDF5 file.  

Parameters
----------
* `reactionModel` :  
    Memory at which to place the reactionModel  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setReactionModel "
`setReactionModel(reactionModel)`  

Pushes the reaction model into the HDF5 file.  

Parameters
----------
* `reactionModel` :  
    Memory from which to place the reactionModel  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getDiffusionModel "
`getDiffusionModel(diffusionModel)`  

Pops the diffusion model from the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memory at which to place the diffusionModel  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getDiffusionModelLattice "
`getDiffusionModelLattice(diffusionModel, lattice, latticeMaxSize, latticeSites, latticeSitesMaxSize)`  

Pops the diffusion lattice from the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memoryfrom which to get the lattice  
* `lattice` :  
    Memory in which to place lattice  
* `latticeMaxSize` :  
    Total size of the lattice (i.e. number of bytes)  
* `latticeSites` :  
    Actual lattice site data in byte format  
* `latticeSitesMaxSize` :  
    Max size of a lattice site (i.e. number of particles it can contain)  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getDiffusionModelLattice "
`getDiffusionModelLattice(diffusionModel, lattice)`  

Pops the diffusion lattice from the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memoryfrom which to get the lattice  
* `lattice` :  
    Memory at which to place the lattice  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setDiffusionModel "
`setDiffusionModel(diffusionModel)`  

Pushes the diffusion model into the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memory from which to place the diffusionModel  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setDiffusionModelLattice "
`setDiffusionModelLattice(m, lattice, latticeSites)`  

Pushes the diffusion lattice into the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memoryfrom which to get the lattice  
* `lattice` :  
    Memory from which to place the lattice - uint8_t particle lattice  
* `latticeSites` :  
    Memory from which to place lattice contents  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setDiffusionModelLattice "
`setDiffusionModelLattice(m, lattice, latticeSites)`  

Pushes the diffusion lattice into the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memoryfrom which to get the lattice  
* `lattice` :  
    Memory from which to place the lattice - uint32_t particle lattice  
* `latticeSites` :  
    Memory from which to place lattice contents  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setDiffusionModelLattice "
`setDiffusionModelLattice(m, lattice)`  

Pushes the diffusion lattice into the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memory from which to get the lattice  
* `lattice` :  
    Memory from which to place the lattice  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setDiffusionModelLattice "
`setDiffusionModelLattice(m, lattice)`  

Pushes the diffusion lattice into the HDF5 file.  

Parameters
----------
* `diffusionModel` :  
    Memory from which to get the lattice  
* `lattice` :  
    Memory from which to place the lattice  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setSpatialModel "
`setSpatialModel(model)`  

Pushes the spacial model (i.e. obstacles) into the HDF5 file.  

Parameters
----------
* `model` :  
    The model object  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getSpatialModel "
`getSpatialModel(model)`  

Pops the spacial model (i.e. obstacles) from the HDF5 file.  

Parameters
----------
* `model` :  
    The model in which to store the object  
";

%feature("docstring") lm::io::hdf5::SimulationFile::replicateExists "
`replicateExists(replicate) -> bool`  

Checks if the specified replicate exists.  

Parameters
----------
* `replicate` :  
    Replicate number  

Returns
-------
true/false  
";

%feature("docstring") lm::io::hdf5::SimulationFile::openReplicate "
`openReplicate(replicate)`  

Opens the specified replicate for reading.  

Parameters
----------
* `replicate` :  
    Replicate number  
";

%feature("docstring") lm::io::hdf5::SimulationFile::appendSpeciesCounts "
`appendSpeciesCounts(replicate, speciesCounts)`  

Appends the species counts in the various sites into the replicate.  

Parameters
----------
* `replicate` :  
    Replicate number  
* `speciesCounts` :  
    Protocol buffers object with species counts in lattice  
";

%feature("docstring") lm::io::hdf5::SimulationFile::appendLattice "
`appendLattice(replicate, lattice, latticeData, latticeDataSize)`  

Appends the lattice to the current replicate.  

Parameters
----------
* `replicate` :  
    Replicate number  
* `lattice` :  
    Lattice to add to replicate  
* `latticeData` :  
    The actual data of the lattice  
* `latticeDataSize` :  
    The size of the lattice  
";

%feature("docstring") lm::io::hdf5::SimulationFile::appendLattice_U32LE "
`appendLattice_U32LE(replicate, lattice, latticeData, latticeDataSize)`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::arbitraryH5 "
`arbitraryH5(ptr)`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::arbitraryH5Lookup "
`arbitraryH5Lookup(in, out, out_sz)`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::appendSites "
`appendSites(replicate, lattice, siteData, siteDataSize)`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::appendParameterValues "
`appendParameterValues(replicate, parameterValues)`  

Adds all the parameter values to the replicate.  

Parameters
----------
* `replicate` :  
    Replicate number  
* `parameterValues` :  
    Actual values to store  
";

%feature("docstring") lm::io::hdf5::SimulationFile::setFirstPassageTimes "
`setFirstPassageTimes(replicate, speciesCounts)`  

Adds all the first passage times to the replicate.  

Parameters
----------
* `replicate` :  
    Replicate number  
* `speciesCounts` :  
    Actual passage times to store  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getLatticeTimes "
`getLatticeTimes(replicate) -> vector< double >`  

Get the timestep times for the replicate.  

Parameters
----------
* `replicate` :  
    Replicate number  

Returns
-------
A set of timestep times  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getLattice "
`getLattice(replicate, latticeIndex, lattice)`  

Get the lattice for the replicate.  

Parameters
----------
* `replicate` :  
    Replicate number  
* `latticeIndex` :  
    Seems to be unused...  
* `lattice` :  
    The lattice object in which to store the data  
";

%feature("docstring") lm::io::hdf5::SimulationFile::closeReplicate "
`closeReplicate(replicate)`  

Close the specified replicate.  

Parameters
----------
* `replicate` :  
    Replicate number  
";

%feature("docstring") lm::io::hdf5::SimulationFile::closeAllReplicates "
`closeAllReplicates()`  

Close all the replicates.  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getParticleCounts "
`getParticleCounts(replicate, latticeIndex) -> std::map< uint32_t, uint >`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getSpeciesCounts "
`getSpeciesCounts(replicate) -> std::map< double, vector< int > >`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getSpeciesCountTimes "
`getSpeciesCountTimes(replicate) -> vector< double >`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getSpeciesNames "
`getSpeciesNames() -> std::map< uint, string >`  
";

%feature("docstring") lm::io::hdf5::SimulationFile::getSiteTypeNames "
`getSiteTypeNames() -> std::map< uint, string >`  
";

// File: classlm_1_1io_1_1SimulationParameters.xml


%feature("docstring") lm::io::SimulationParameters "

C++ includes: SimulationParameters.h
";

%feature("docstring") lm::io::SimulationParameters::fromMessage "
`fromMessage(message) -> map< string, string >`  
";

%feature("docstring") lm::io::SimulationParameters::intoMessage "
`intoMessage(message, parameters)`  
";

// File: structlm_1_1cme_1_1CMESolver_1_1SpeciesLimit.xml

// File: classlm_1_1builder_1_1Sphere.xml


%feature("docstring") lm::builder::Sphere "

A Shape that represents a Sphere.  

C++ includes: Sphere.h
";

%feature("docstring") lm::builder::Sphere::Sphere "
`Sphere(center, radius, type)`  

Create a Sphere.  

Parameters
----------
* `center` :  
    Point center of the circle of the slice plane through the sphere  
* `radius` :  
    Radius of the sphere  
* `type` :  
    The type of the sites within the sphere  
";

%feature("docstring") lm::builder::Sphere::~Sphere "
`~Sphere()`  

Destroy the Sphere.  
";

%feature("docstring") lm::builder::Sphere::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Sphere::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Sphere::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Sphere::setCenter "
`setCenter(center)`  

Set the center of the sphere.  

Parameters
----------
* `center` :  
    Point of the center  
";

%feature("docstring") lm::builder::Sphere::getCenter "
`getCenter() -> point`  

Get the center of the sphere.  
";

%feature("docstring") lm::builder::Sphere::getRadius "
`getRadius() -> si_dist_t`  

Get the radius of the sphere.  
";

%feature("docstring") lm::builder::Sphere::getVolume "
`getVolume() -> double`  

Get the volume bounded by the sphere.  
";

// File: classlm_1_1thread_1_1Thread.xml


%feature("docstring") lm::thread::Thread "

A base class wrapping pthreads.  

C++ includes: Thread.h
";

%feature("docstring") lm::thread::Thread::Thread "
`Thread()`  

Creates a pthread locking mechanism and initializes \"Thread\".  
";

%feature("docstring") lm::thread::Thread::~Thread "
`~Thread()`  

Destory the Thread.  
";

%feature("docstring") lm::thread::Thread::start "
`start()`  

If no thread exists, creates a new thread and begins execution.  
";

%feature("docstring") lm::thread::Thread::stop "
`stop()`  

Joins the thread with the parent waiting if necessary.  
";

%feature("docstring") lm::thread::Thread::wake "
`wake()`  

Wakes a sleeping thead.  
";

%feature("docstring") lm::thread::Thread::getId "
`getId() -> pthread_t`  

Returns the pthread based id for the Thread.  

Returns
-------
pthread assigned id  
";

%feature("docstring") lm::thread::Thread::setAffinity "
`setAffinity(cpuNumber)`  

Binds the thread to a CPU core.  

Parameters
----------
* `The` :  
    cpu core for which to bind the thread  
";

// File: classlm_1_1builder_1_1Torus.xml


%feature("docstring") lm::builder::Torus "

A Shape that represents a Torus.  

C++ includes: Torus.h
";

%feature("docstring") lm::builder::Torus::Torus "
`Torus(center, r1, r2, type, orientation=vector(0.0, 0.0, 1.0), up=vector(0.0, 1.0, 0.0))`  

Create a Torus.  

Parameters
----------
* `center` :  
    Point center of the circle of the slice plane through the sphere  
* `r1` :  
    Large radius of the Torus  
* `r2` :  
    Small radius of the Torus  
* `type` :  
    The type of the sites within the sphere  
* `orientation` :  
    The direction vector around which the ring is made, default: (0,0,1)  
* `up` :  
    A vector to define the y axis of the object's coordinate system , default: (0,1,0)  
";

%feature("docstring") lm::builder::Torus::~Torus "
`~Torus()`  

Destroy the Torus.  
";

%feature("docstring") lm::builder::Torus::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Torus::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Torus::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Torus::setCenter "
`setCenter(center)`  

Set the center of the sphere.  

Parameters
----------
* `center` :  
    Point of the center  
";

%feature("docstring") lm::builder::Torus::getCenter "
`getCenter() -> point`  

Get the center of the sphere.  
";

%feature("docstring") lm::builder::Torus::getRadius1 "
`getRadius1() -> si_dist_t`  

Get the large radius of the sphere.  
";

%feature("docstring") lm::builder::Torus::getRadius2 "
`getRadius2() -> si_dist_t`  

Get the small radius of the sphere.  
";

%feature("docstring") lm::builder::Torus::getVolume "
`getVolume() -> double`  

Get the volume bounded by the sphere.  
";

// File: structlm_1_1cme_1_1CMESolver_1_1TrackedParameter.xml

// File: classlm_1_1cme_1_1TwoStateExpression.xml


%feature("docstring") lm::cme::TwoStateExpression "

C++ includes: TwoStateExpression.h
";

%feature("docstring") lm::cme::TwoStateExpression::TwoStateExpression "
`TwoStateExpression()`  
";

%feature("docstring") lm::cme::TwoStateExpression::~TwoStateExpression "
`~TwoStateExpression()`  
";

%feature("docstring") lm::cme::TwoStateExpression::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::TwoStateExpression::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::TwoStateExpression::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1cme_1_1TwoStateHillLoopSwitch.xml


%feature("docstring") lm::cme::TwoStateHillLoopSwitch "

C++ includes: TwoStateHillLoopSwitch.h
";

%feature("docstring") lm::cme::TwoStateHillLoopSwitch::TwoStateHillLoopSwitch "
`TwoStateHillLoopSwitch()`  
";

%feature("docstring") lm::cme::TwoStateHillLoopSwitch::~TwoStateHillLoopSwitch "
`~TwoStateHillLoopSwitch()`  
";

%feature("docstring") lm::cme::TwoStateHillLoopSwitch::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::TwoStateHillLoopSwitch::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::TwoStateHillLoopSwitch::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

// File: classlm_1_1cme_1_1TwoStateHillSwitch.xml


%feature("docstring") lm::cme::TwoStateHillSwitch "

C++ includes: TwoStateHillSwitch.h
";

%feature("docstring") lm::cme::TwoStateHillSwitch::TwoStateHillSwitch "
`TwoStateHillSwitch()`  
";

%feature("docstring") lm::cme::TwoStateHillSwitch::~TwoStateHillSwitch "
`~TwoStateHillSwitch()`  
";

%feature("docstring") lm::cme::TwoStateHillSwitch::needsReactionModel "
`needsReactionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::TwoStateHillSwitch::needsDiffusionModel "
`needsDiffusionModel() -> bool`  

Tells whether the solver needs a reaction model.  
";

%feature("docstring") lm::cme::TwoStateHillSwitch::generateTrajectory "
`generateTrajectory()`  

Actually run the simulation.  
";

%feature("docstring") lm::cme::TwoStateHillSwitch::runConstantProtein "
`runConstantProtein()`  
";

%feature("docstring") lm::cme::TwoStateHillSwitch::runGeometricProtein "
`runGeometricProtein()`  
";

// File: structlm_1_1cme_1_1CMESolver_1_1UncompetitiveMMPropensityArgs.xml

// File: classlm_1_1builder_1_1Union.xml


%feature("docstring") lm::builder::Union "

A Shape that represents a Union between two shapes.  

C++ includes: Union.h
";

%feature("docstring") lm::builder::Union::Union "
`Union(s1, s2, type)`  

Create a Union.  

Parameters
----------
* `s1` :  
    The first shape to Union  
* `s2` :  
    The second shape to Union  
* `type` :  
    The type of the sites within the union  
";

%feature("docstring") lm::builder::Union::~Union "
`~Union()`  

Destroy the Sphere.  
";

%feature("docstring") lm::builder::Union::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Union::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Union::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::Union::getVolume "
`getVolume(reintegrate=false) -> double`  

Get the volume bounded by the sphere.  
";

%feature("docstring") lm::builder::Union::getVolume "
`getVolume() -> double`  

Get the total internal volume of the shape.  
";

// File: classlm_1_1builder_1_1UnionSet.xml


%feature("docstring") lm::builder::UnionSet "

A Shape that represents a Union of many shapes.  

C++ includes: UnionSet.h
";

%feature("docstring") lm::builder::UnionSet::UnionSet "
`UnionSet(type)`  

Create a UnionSet.  

Parameters
----------
* `type` :  
    The type of the sites within the union  
";

%feature("docstring") lm::builder::UnionSet::~UnionSet "
`~UnionSet()`  

Destroy the Sphere.  
";

%feature("docstring") lm::builder::UnionSet::addShape "
`addShape(s)`  

Add a shape to the union.  

Parameters
----------
* `shape` :  
    A Shape object  
";

%feature("docstring") lm::builder::UnionSet::intersects "
`intersects(query) -> bool`  

Check if two shapes intersect.  

Parameters
----------
* `query` :  
    The other shape to check  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::UnionSet::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified point.  

Parameters
----------
* `query` :  
    Point to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::UnionSet::contains "
`contains(query) -> bool`  

Determine if the shape contains the specified shape.  

Parameters
----------
* `query` :  
    Shape to test  

Returns
-------
true/false  
";

%feature("docstring") lm::builder::UnionSet::getVolume "
`getVolume(reintegrate=false) -> double`  

Get the volume bounded by the sphere.  
";

%feature("docstring") lm::builder::UnionSet::getVolume "
`getVolume() -> double`  

Get the total internal volume of the shape.  
";

// File: structlm_1_1builder_1_1vector.xml


%feature("docstring") lm::builder::vector "

A vector which points in a direction.  

C++ includes: Shape.h
";

%feature("docstring") lm::builder::vector::vector "
`vector(x=0.0, y=0.0, z=0.0)`  

Create a vector.  

Parameters
----------
* `x` :  
    Distance in x  
* `y` :  
    Distance in y  
* `z` :  
    Distance in z  
";

%feature("docstring") lm::builder::vector::vector "
`vector(from, to)`  

Construct a vector pointing from one point to another.  

Parameters
----------
* `from` :  
    The point at the root of the vector  
* `to` :  
    The point at the end of the vector  
";

%feature("docstring") lm::builder::vector::vector "
`vector(v)`  

Construct a vector from another vector.  

Parameters
----------
* `v` :  
    The vector to copy  
";

%feature("docstring") lm::builder::vector::vector "
`vector(p)`  

Construct a vector from a point.  

Parameters
----------
* `p` :  
    The point to cast  
";

%feature("docstring") lm::builder::vector::toPoint "
`toPoint() -> point`  

convert the vector to a point  
";

%feature("docstring") lm::builder::vector::length "
`length() -> si_dist_t`  

Get the vector length.  
";

%feature("docstring") lm::builder::vector::unitize "
`unitize() -> vector`  

Get a unit vector pointing in the same direction as this vector.  
";

%feature("docstring") lm::builder::vector::dot "
`dot(r) -> si_dist_t`  

Compute the dot product between vectors.  

Parameters
----------
* `r` :  
    The right hand vector  
";

%feature("docstring") lm::builder::vector::cross "
`cross(r) -> vector`  

Compute the cross product between vectors.  

Parameters
----------
* `r` :  
    The right hand vector  
";

%feature("docstring") lm::builder::vector::scale "
`scale(r) -> vector`  

Scale the vector by a constant.  

Parameters
----------
* `r` :  
    The scalar  
";

%feature("docstring") lm::builder::vector::mult "
`mult(r) -> matrix`  

Multiply with another vector (taking this vector to be a row vector) to form a matrix.  

Parameters
----------
* `r` :  
    The column vector  
";

%feature("docstring") lm::builder::vector::mult "
`mult(r) -> vector`  

Multiply with a matrix to form a new vector (assuming this vector is a row vector)  

Parameters
----------
* `r` :  
    The matrix  
";

%feature("docstring") lm::builder::vector::findPerpendicularVector "
`findPerpendicularVector(v) -> vector`  
";

// File: classlm_1_1thread_1_1Worker.xml


%feature("docstring") lm::thread::Worker "

An actual worker Thread that runs a simulation replicate.  

C++ includes: Worker.h
";

%feature("docstring") lm::thread::Worker::Worker "
`Worker()`  

Creates thread and attaches it to the manager.  
";

%feature("docstring") lm::thread::Worker::~Worker "
`~Worker()`  

Removes thread from manager and deletes thread.  
";

%feature("docstring") lm::thread::Worker::abort "
`abort()`  

Set the aborted status for a thread and wakes thread.  
";

%feature("docstring") lm::thread::Worker::checkpoint "
`checkpoint()`  

Checkpoint (currently unimplemented)  
";

// File: classlm_1_1thread_1_1WorkerManager.xml


%feature("docstring") lm::thread::WorkerManager "

A singleton manager that creates and manages workers that run on a processing node.  

C++ includes: WorkerManager.h
";

%feature("docstring") lm::thread::WorkerManager::getInstance "
`getInstance() -> WorkerManager *`  

Get the global worker manager.  

Returns
-------
globalWorkerManager Global thread manager instance handle  
";

%feature("docstring") lm::thread::WorkerManager::WorkerManager "
`WorkerManager()`  

Create the WorkerManager.  
";

%feature("docstring") lm::thread::WorkerManager::~WorkerManager "
`~WorkerManager()`  

Destroy the WorkerManager.  
";

%feature("docstring") lm::thread::WorkerManager::addWorker "
`addWorker(worker)`  

Adds a worker to the manager.  

Parameters
----------
* `worker` :  
    Worker to add  
";

%feature("docstring") lm::thread::WorkerManager::removeWorker "
`removeWorker(worker)`  

Removes a worker from the manager.  

Parameters
----------
* `worker` :  
    Worker handle for which to remove  
";

%feature("docstring") lm::thread::WorkerManager::abortWorkers "
`abortWorkers()`  

Aborts all the worker threads.  
";

%feature("docstring") lm::thread::WorkerManager::checkpointWorkers "
`checkpointWorkers()`  

Causes all the worker threads to checkpoint (checkpointing is currently unimplemented)  
";

%feature("docstring") lm::thread::WorkerManager::stopWorkers "
`stopWorkers()`  

Stops all the worker threads (e.g. merges them with the master)  
";

// File: classlm_1_1rng_1_1XORShift.xml


%feature("docstring") lm::rng::XORShift "

A random number generator that works for a grid based on the index.  

C++ includes: XORShift.h
";

%feature("docstring") lm::rng::XORShift::XORShift "
`XORShift(seedTop, seedBottom)`  

Create the XORShift random number generator.  

Parameters
----------
* `seedTop` :  
    Top 32 bits of the random seed  
* `seedBottom` :  
    Bottom 32 bits of the random seed  
";

%feature("docstring") lm::rng::XORShift::~XORShift "
`~XORShift()`  

Destory the XORShift object.  
";

%feature("docstring") lm::rng::XORShift::getRandom "
`getRandom() -> uint32_t`  

Get a random integer.  

Returns
-------
random A random integer  
";

%feature("docstring") lm::rng::XORShift::getRandomDouble "
`getRandomDouble() -> double`  

Get a random double.  

Returns
-------
random A random double  

Returns a float value in the range [0.0 1.0).  
";

%feature("docstring") lm::rng::XORShift::getExpRandomDouble "
`getExpRandomDouble() -> double`  

Get a random exponentially distributed double.  

Returns
-------
random A random double  

Returns an exponentially distributed value.  
";

%feature("docstring") lm::rng::XORShift::getNormRandomDouble "
`getNormRandomDouble() -> double`  

Get a random normally distributed double.  

Returns
-------
random A random double  

Returns an normally distributed value.  
";

%feature("docstring") lm::rng::XORShift::getRandomDoubles "
`getRandomDoubles(rngs, numberRNGs)`  

Get a set of random doubles.  

Parameters
----------
* `rngs` :  
    Pointer to memory for random numbers  
* `numberRNGs` :  
    The size of the rngs variable  
";

%feature("docstring") lm::rng::XORShift::getExpRandomDoubles "
`getExpRandomDoubles(rngs, numberRNGs)`  

Get a set of random exponentially distributed doubles.  

Parameters
----------
* `rngs` :  
    Pointer to memory for random numbers  
* `numberRNGs` :  
    The size of the rngs variable  
";

%feature("docstring") lm::rng::XORShift::getNormRandomDoubles "
`getNormRandomDoubles(rngs, numberRNGs)`  

Get a set of random normally distributed doubles.  

Parameters
----------
* `rngs` :  
    Pointer to memory for random numbers  
* `numberRNGs` :  
    The size of the rngs variable  
";

// File: classlm_1_1rng_1_1XORWow.xml


%feature("docstring") lm::rng::XORWow "

C++ includes: XORWow.h
";

%feature("docstring") lm::rng::XORWow::XORWow "
`XORWow(cudaDevice, seedTop, seedBottom, availableDists)`  
";

%feature("docstring") lm::rng::XORWow::~XORWow "
`~XORWow()`  
";

%feature("docstring") lm::rng::XORWow::getRandom "
`getRandom() -> uint32_t`  

Get a random integer.  

Returns
-------
random A random integer  
";

%feature("docstring") lm::rng::XORWow::getRandomDouble "
`getRandomDouble() -> double`  

Get a random double.  

Returns
-------
random A random double  
";

%feature("docstring") lm::rng::XORWow::getExpRandomDouble "
`getExpRandomDouble() -> double`  

Get a random exponentially distributed double.  

Returns
-------
random A random double  
";

%feature("docstring") lm::rng::XORWow::getNormRandomDouble "
`getNormRandomDouble() -> double`  

Get a random normally distributed double.  

Returns
-------
random A random double  
";

%feature("docstring") lm::rng::XORWow::getRandomDoubles "
`getRandomDoubles(rngs, numberRNGs)`  

Get a number of random doubles.  

Parameters
----------
* `rngs` :  
    A preallocated array for random numbers  
* `numberRNGs` :  
    Number or randoms to put in array  
";

%feature("docstring") lm::rng::XORWow::getExpRandomDoubles "
`getExpRandomDoubles(rngs, numberRNGs)`  

Get a number of random exponentially distiributed doubles.  

Parameters
----------
* `rngs` :  
    A preallocated array for random numbers  
* `numberRNGs` :  
    Number or randoms to put in array  
";

%feature("docstring") lm::rng::XORWow::getNormRandomDoubles "
`getNormRandomDoubles(rngs, numberRNGs)`  

Get a number of random normally distributed doubles.  

Parameters
----------
* `rngs` :  
    A preallocated array for random numbers  
* `numberRNGs` :  
    Number or randoms to put in array  
";

// File: classZDivMPIGPUMapper.xml


%feature("docstring") ZDivMPIGPUMapper "

C++ includes: ZDivMPIGPUMapper.h
";

%feature("docstring") ZDivMPIGPUMapper::ZDivMPIGPUMapper "
`ZDivMPIGPUMapper(x, y, z, arg4, arg5, arg6, arg7, pz=false, pages=1)`  
";

%feature("docstring") ZDivMPIGPUMapper::~ZDivMPIGPUMapper "
`~ZDivMPIGPUMapper()`  
";

%feature("docstring") ZDivMPIGPUMapper::initialize "
`initialize() -> SegmentDescriptor_s *`  
";

%feature("docstring") ZDivMPIGPUMapper::initialize_gpu "
`initialize_gpu()`  
";

%feature("docstring") ZDivMPIGPUMapper::get_global_dim "
`get_global_dim() -> dim3`  
";

%feature("docstring") ZDivMPIGPUMapper::get_local_dim "
`get_local_dim() -> dim3`  
";

%feature("docstring") ZDivMPIGPUMapper::get_global_offset "
`get_global_offset() -> int3`  
";

%feature("docstring") ZDivMPIGPUMapper::get_local_size "
`get_local_size() -> size_t`  
";

%feature("docstring") ZDivMPIGPUMapper::get_authority_size "
`get_authority_size() -> size_t`  
";

%feature("docstring") ZDivMPIGPUMapper::get_global_input_offset "
`get_global_input_offset() -> ssize_t`  
";

%feature("docstring") ZDivMPIGPUMapper::get_global_output_offset "
`get_global_output_offset() -> size_t`  
";

%feature("docstring") ZDivMPIGPUMapper::get_authority_offset "
`get_authority_offset() -> size_t`  
";

%feature("docstring") ZDivMPIGPUMapper::stage_in "
`stage_in(dptr, hptr)`  
";

%feature("docstring") ZDivMPIGPUMapper::stage_in_sites "
`stage_in_sites(dptr, hptr)`  
";

%feature("docstring") ZDivMPIGPUMapper::stage_out "
`stage_out(hptr, dptr)`  
";

%feature("docstring") ZDivMPIGPUMapper::communicate_edges "
`communicate_edges(dptr, key)`  
";

%feature("docstring") ZDivMPIGPUMapper::recv_edge "
`recv_edge(key, neigh)`  
";

%feature("docstring") ZDivMPIGPUMapper::send_edge "
`send_edge(key, neigh)`  
";

%feature("docstring") ZDivMPIGPUMapper::copy_edge_to_device "
`copy_edge_to_device(dptr, key, neighbor, s)`  
";

%feature("docstring") ZDivMPIGPUMapper::copy_edge_to_host "
`copy_edge_to_host(dptr, key, neighbor, s)`  
";

%feature("docstring") ZDivMPIGPUMapper::stage_in_from_slave "
`stage_in_from_slave(rank, hptr)`  
";

%feature("docstring") ZDivMPIGPUMapper::stage_out_to_master "
`stage_out_to_master(hptr)`  
";

// File: classZDivMultiGPUMapper.xml


%feature("docstring") ZDivMultiGPUMapper "

C++ includes: ZDivMultiGPUMapper.h
";

%feature("docstring") ZDivMultiGPUMapper::ZDivMultiGPUMapper "
`ZDivMultiGPUMapper(arg1, arg2, arg3, arg4, arg5, gl=NULL, pz=false, pages=1)`  
";

%feature("docstring") ZDivMultiGPUMapper::ZDivMultiGPUMapper "
`ZDivMultiGPUMapper(x, y, z, arg4, arg5, arg6, arg7, gl=NULL, pz=false, pages=1)`  
";

%feature("docstring") ZDivMultiGPUMapper::~ZDivMultiGPUMapper "
`~ZDivMultiGPUMapper()`  
";

%feature("docstring") ZDivMultiGPUMapper::initialize_gpu "
`initialize_gpu(gpu)`  
";

%feature("docstring") ZDivMultiGPUMapper::get_global_dim "
`get_global_dim(gpu) -> dim3`  
";

%feature("docstring") ZDivMultiGPUMapper::get_local_dim "
`get_local_dim(gpu) -> dim3`  
";

%feature("docstring") ZDivMultiGPUMapper::get_global_offset "
`get_global_offset(gpu) -> int3`  
";

%feature("docstring") ZDivMultiGPUMapper::get_local_size "
`get_local_size(gpu) -> size_t`  
";

%feature("docstring") ZDivMultiGPUMapper::get_authority_size "
`get_authority_size(gpu) -> size_t`  
";

%feature("docstring") ZDivMultiGPUMapper::get_global_input_offset "
`get_global_input_offset(gpu) -> ssize_t`  
";

%feature("docstring") ZDivMultiGPUMapper::get_global_output_offset "
`get_global_output_offset(gpu) -> size_t`  
";

%feature("docstring") ZDivMultiGPUMapper::get_authority_offset "
`get_authority_offset(gpu) -> size_t`  
";

%feature("docstring") ZDivMultiGPUMapper::stage_in "
`stage_in(gpu, dptr, hptr)`  
";

%feature("docstring") ZDivMultiGPUMapper::stage_in_sites "
`stage_in_sites(gpu, dptr, hptr)`  
";

%feature("docstring") ZDivMultiGPUMapper::stage_out "
`stage_out(gpu, hptr, dptr)`  
";

%feature("docstring") ZDivMultiGPUMapper::publish_state "
`publish_state(gpu, dptr, key)`  
";

%feature("docstring") ZDivMultiGPUMapper::refresh "
`refresh(gpu, dptr, key)`  
";

%feature("docstring") ZDivMultiGPUMapper::schedule_send "
`schedule_send(gpu, dptr, timestamp, neighbor, stream)`  
";

%feature("docstring") ZDivMultiGPUMapper::schedule_recv "
`schedule_recv(gpu, dptr, timestamp, neighbor, stream)`  
";

%feature("docstring") ZDivMultiGPUMapper::map_index_to_gpu "
`map_index_to_gpu(index) -> int`  
";

%feature("docstring") ZDivMultiGPUMapper::getinfo "
`getinfo(gpu) -> gpu_info *`  
";

%feature("docstring") ZDivMultiGPUMapper::gettbuf "
`gettbuf(gpu, key, neighbor) -> unsigned int *`  
";

%feature("docstring") ZDivMultiGPUMapper::getrbuf "
`getrbuf(gpu, key, neighbor) -> unsigned int *`  
";

%feature("docstring") ZDivMultiGPUMapper::send_tbuf "
`send_tbuf(gpu, key, neighbor, s)`  
";

// File: structlm_1_1cme_1_1CMESolver_1_1ZerothOrderHeavisidePropensityArgs.xml

// File: structlm_1_1cme_1_1CMESolver_1_1ZerothOrderPropensityArgs.xml

// File: namespacelm.xml

// File: namespacelm_1_1builder.xml

// File: namespacelm_1_1cme.xml

// File: namespacelm_1_1io.xml

// File: namespacelm_1_1io_1_1hdf5.xml

// File: namespacelm_1_1main.xml

// File: namespacelm_1_1me.xml

// File: namespacelm_1_1message.xml

// File: namespacelm_1_1rdme.xml

%feature("docstring") lm::rdme::IntMpdRdmeSolver_x_kernel "
`IntMpdRdmeSolver_x_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver_y_kernel "
`IntMpdRdmeSolver_y_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver_z_kernel "
`IntMpdRdmeSolver_z_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::IntMpdRdmeSolver_reaction_kernel "
`IntMpdRdmeSolver_reaction_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::getCompiledLatticeMaxOccupancy "
`getCompiledLatticeMaxOccupancy() -> unsigned int`  
";

%feature("docstring") lm::rdme::MGPU_x_kernel "
`MGPU_x_kernel(inLattice, inSites, outLattice, z_start, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MGPU_y_kernel "
`MGPU_y_kernel(inLattice, inSites, outLattice, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MGPU_z_kernel "
`MGPU_z_kernel(inLattice, inSites, outLattice, timestepHash, siteOverflowList, z_start) -> __global__ void`  
";

%feature("docstring") lm::rdme::MGPU_reaction_kernel "
`MGPU_reaction_kernel(inLattice, inSites, outLattice, timestepHash, siteOverflowList, z_start) -> __global__ void`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver_x_kernel "
`MpdRdmeSolver_x_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver_y_kernel "
`MpdRdmeSolver_y_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver_z_kernel "
`MpdRdmeSolver_z_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MpdRdmeSolver_reaction_kernel "
`MpdRdmeSolver_reaction_kernel(inLattice, inSites, outLattice, gridXSize, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MPI_x_kernel "
`MPI_x_kernel(inLattice, inSites, outLattice, z_start, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MPI_y_kernel "
`MPI_y_kernel(inLattice, inSites, outLattice, timestepHash, siteOverflowList) -> __global__ void`  
";

%feature("docstring") lm::rdme::MPI_z_kernel "
`MPI_z_kernel(inLattice, inSites, outLattice, timestepHash, siteOverflowList, z_start) -> __global__ void`  
";

%feature("docstring") lm::rdme::MPI_reaction_kernel "
`MPI_reaction_kernel(inLattice, inSites, outLattice, timestepHash, siteOverflowList, z_start) -> __global__ void`  
";

// File: namespacelm_1_1reaction.xml

// File: namespacelm_1_1rng.xml

%feature("docstring") lm::rng::xorwow_init_kernel "
`xorwow_init_kernel(seed, rngState) -> __global__ void`  
";

%feature("docstring") lm::rng::xorwow_generate_kernel "
`xorwow_generate_kernel(state, randomValues, expRandomValues, normRandomValues, iterations) -> __global__ void`  
";

// File: namespacelm_1_1thread.xml

// File: Capsule_8cpp.xml

// File: Capsule_8h.xml

// File: CapsuleShell_8cpp.xml

// File: CapsuleShell_8h.xml

// File: Cone_8cpp.xml

// File: Cone_8h.xml

// File: Cuboid_8cpp.xml

// File: Cuboid_8h.xml

// File: Cylinder_8cpp.xml

// File: Cylinder_8h.xml

// File: Difference_8cpp.xml

// File: Difference_8h.xml

// File: Ellipse_8cpp.xml

// File: Ellipse_8h.xml

// File: Hemisphere_8cpp.xml

// File: Hemisphere_8h.xml

// File: Intersection_8cpp.xml

// File: Intersection_8h.xml

// File: LatticeBuilder_8cpp.xml

// File: LatticeBuilder_8h.xml

// File: Mesh_8cpp.xml

// File: Mesh_8h.xml

// File: Shape_8cpp.xml

// File: Shape_8h.xml

// File: Sphere_8cpp.xml

// File: Sphere_8h.xml

// File: Torus_8cpp.xml

// File: Torus_8h.xml

// File: Union_8cpp.xml

// File: Union_8h.xml

// File: UnionSet_8cpp.xml

// File: UnionSet_8h.xml

// File: common_8cpp.xml

%feature("docstring") getPhysicalCpuCores "
`getPhysicalCpuCores() -> int`  

Prints the copyright notice. Gets the number of physical cpu cores on the system.  
";

%feature("docstring") parseArguments "
`parseArguments(argc, argv, defaultSolver)`  

Parses the command line arguments.  
";

%feature("docstring") parseIntListArg "
`parseIntListArg(list, arg)`  
";

%feature("docstring") parseTimeArg "
`parseTimeArg(arg) -> time_t`  
";

%feature("docstring") parseIntReciprocalArg "
`parseIntReciprocalArg(arg) -> float`  
";

%feature("docstring") printUsage "
`printUsage(argc, argv)`  

Prints the usage for the program.  
";

// File: common_8h.xml

%feature("docstring") printCopyright "
`printCopyright(argc, argv)`  

This function prints the copyright notice.  
";

%feature("docstring") getPhysicalCpuCores "
`getPhysicalCpuCores() -> int`  

Prints the copyright notice. Gets the number of physical cpu cores on the system.  
";

%feature("docstring") parseArguments "
`parseArguments(argc, argv, defaultSolver)`  

Parses the command line arguments.  
";

%feature("docstring") parseIntListArg "
`parseIntListArg(list, option)`  
";

%feature("docstring") parseTimeArg "
`parseTimeArg(option) -> time_t`  
";

%feature("docstring") parseIntReciprocalArg "
`parseIntReciprocalArg(option) -> float`  
";

%feature("docstring") printUsage "
`printUsage(argc, argv)`  

Prints the usage for the program.  
";

%feature("docstring") discoverEnvironment "
`discoverEnvironment()`  

Figures out the environment for the program (directories, files, etc).  
";

%feature("docstring") initPython "
`initPython()`  

Starts the python interpreter.  
";

%feature("docstring") finalizePython "
`finalizePython()`  
";

%feature("docstring") startInterpreter "
`startInterpreter()`  
";

%feature("docstring") executeScript "
`executeScript(filename, arguments, replicate=0)`  
";

// File: lm__python_8cpp.xml

%feature("docstring") mkwStr "
`mkwStr(str) -> wstring`  
";

%feature("docstring") parseArguments "
`parseArguments(argc, argv)`  

Parses the command line arguments.  
";

%feature("docstring") printUsage "
`printUsage(argc, argv)`  

Prints the usage for the program.  
";

%feature("docstring") discoverEnvironment "
`discoverEnvironment()`  

Figures out the environment for the program (directories, files, etc).  
";

%feature("docstring") initPython "
`initPython()`  

Starts the python interpreter.  
";

%feature("docstring") finalizePython "
`finalizePython()`  
";

%feature("docstring") startInterpreter "
`startInterpreter()`  
";

%feature("docstring") executeScript "
`executeScript(filename, arguments)`  
";

%feature("docstring") executeScript "
`executeScript(filename, arguments)`  
";

%feature("docstring") PyInit__lm "
`PyInit__lm(arg1) -> PyObject *`  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

// File: lm__sbml__import_8cpp.xml

%feature("docstring") parseArguments "
`parseArguments(argc, argv)`  

Parses the command line arguments.  
";

%feature("docstring") printUsage "
`printUsage(argc, argv)`  

Prints the usage for the program.  
";

%feature("docstring") importSBMLModel "
`importSBMLModel(lmFile, sbmlFilename)`  
";

%feature("docstring") importSBMLModelL3V1 "
`importSBMLModelL3V1(lmModel, sbmlModel)`  
";

%feature("docstring") importSBMLModelL3V1Kinetics "
`importSBMLModelL3V1Kinetics(reactionIndex, kinetics, lmModel, D, speciesIndices, numberReactions, globalParameters, globalParameterValues)`  
";

%feature("docstring") isFirstOrderReaction "
`isFirstOrderReaction(root, parameters) -> bool`  
";

%feature("docstring") importFirstOrderReaction "
`importFirstOrderReaction(root, parameters, parameterValues, reactionIndex, numberReactions, lmModel, D, speciesIndices)`  
";

%feature("docstring") isSecondOrderReaction "
`isSecondOrderReaction(root, parameters) -> bool`  
";

%feature("docstring") importSecondOrderReaction "
`importSecondOrderReaction(root, parameters, parameterValues, reactionIndex, numberReactions, lmModel, D, speciesIndices)`  
";

%feature("docstring") isSecondOrderSelfReaction "
`isSecondOrderSelfReaction(root, parameters) -> bool`  
";

%feature("docstring") importSecondOrderSelfReaction "
`importSecondOrderSelfReaction(root, parameters, parameterValues, reactionIndex, numberReactions, lmModel, D, speciesIndices)`  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

%feature("docstring") getSpeciesUsedInExpression "
`getSpeciesUsedInExpression(speciesUsed, node, parameters)`  
";

%feature("docstring") getOperatorsUsedInExpression "
`getOperatorsUsedInExpression(speciesUsed, node)`  
";

%feature("docstring") getFirstExpressionOfType "
`getFirstExpressionOfType(node, type) -> const ASTNode *`  
";

%feature("docstring") calculateMultiplierInExpression "
`calculateMultiplierInExpression(node, parameterValues, ignoreSpeciesMinusOne=false) -> double`  
";

// File: lm__setdm_8cpp.xml

%feature("docstring") parseArguments "
`parseArguments(argc, argv)`  

Parses the command line arguments.  
";

%feature("docstring") printUsage "
`printUsage(argc, argv)`  

Prints the usage for the program.  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

// File: lm__setp_8cpp.xml

%feature("docstring") parseArguments "
`parseArguments(argc, argv)`  

Parses the command line arguments.  
";

%feature("docstring") printUsage "
`printUsage(argc, argv)`  

Prints the usage for the program.  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

// File: lm__setrm_8cpp.xml

%feature("docstring") parseArguments "
`parseArguments(argc, argv)`  

Parses the command line arguments.  
";

%feature("docstring") printUsage "
`printUsage(argc, argv)`  

Prints the usage for the program.  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

// File: Main-MPMPD_8cpp.xml

%feature("docstring") listDevicesMPI "
`listDevicesMPI()`  
";

%feature("docstring") executeSimulationMPI "
`executeSimulationMPI()`  
";

%feature("docstring") executeSimulationMPISingleMaster "
`executeSimulationMPISingleMaster()`  
";

%feature("docstring") executeSimulationMPISingleSlave "
`executeSimulationMPISingleSlave()`  
";

%feature("docstring") broadcastSimulationParameters "
`broadcastSimulationParameters(staticDataBuffer, simulationParameters)`  
";

%feature("docstring") broadcastReactionModel "
`broadcastReactionModel(staticDataBuffer, reactionModel)`  
";

%feature("docstring") broadcastDiffusionModel "
`broadcastDiffusionModel(staticDataBuffer, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize)`  
";

%feature("docstring") receiveSimulationParameters "
`receiveSimulationParameters(staticDataBuffer) -> map< string, string >`  
";

%feature("docstring") receiveReactionModel "
`receiveReactionModel(staticDataBuffer, reactionModel)`  
";

%feature("docstring") receiveDiffusionModel "
`receiveDiffusionModel(staticDataBuffer, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize)`  
";

%feature("docstring") runReplicateNOW "
`runReplicateNOW(replicate, solverFactory, simulationParameters, reactionModel, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resourceAllocator) -> int`  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

// File: MainMPI_8cpp.xml

%feature("docstring") listDevicesMPI "
`listDevicesMPI()`  
";

%feature("docstring") executeSimulationMPI "
`executeSimulationMPI()`  
";

%feature("docstring") executeSimulationMPISingleMaster "
`executeSimulationMPISingleMaster()`  
";

%feature("docstring") executeSimulationMPISingleSlave "
`executeSimulationMPISingleSlave()`  
";

%feature("docstring") broadcastSimulationParameters "
`broadcastSimulationParameters(staticDataBuffer, simulationParameters)`  
";

%feature("docstring") broadcastReactionModel "
`broadcastReactionModel(staticDataBuffer, reactionModel)`  
";

%feature("docstring") broadcastDiffusionModel "
`broadcastDiffusionModel(staticDataBuffer, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize)`  
";

%feature("docstring") receiveSimulationParameters "
`receiveSimulationParameters(staticDataBuffer) -> map< string, string >`  
";

%feature("docstring") receiveReactionModel "
`receiveReactionModel(staticDataBuffer, reactionModel)`  
";

%feature("docstring") receiveDiffusionModel "
`receiveDiffusionModel(staticDataBuffer, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize)`  
";

%feature("docstring") startReplicate "
`startReplicate(replicate, solverFactory, simulationParameters, reactionModel, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resourceAllocator) -> ReplicateRunner *`  
";

%feature("docstring") popNextFinishedReplicate "
`popNextFinishedReplicate(runningReplicates, resourceAllocator) -> ReplicateRunner *`  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

// File: MainSA_8cpp.xml

%feature("docstring") listDevices "
`listDevices()`  
";

%feature("docstring") executeSimulation "
`executeSimulation()`  
";

%feature("docstring") startReplicate "
`startReplicate(replicate, solverFactory, simulationParameters, reactionModel, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resourceAllocator) -> ReplicateRunner *`  
";

%feature("docstring") popNextFinishedReplicate "
`popNextFinishedReplicate(runningReplicates, resourceAllocator) -> ReplicateRunner *`  
";

%feature("docstring") main "
`main(argc, argv) -> int`  
";

// File: CMESolver_8cpp.xml

// File: CMESolver_8h.xml

// File: FluctuatingNRSolver_8cpp.xml

// File: FluctuatingNRSolver_8h.xml

// File: GillespieDSolver_8cpp.xml

// File: GillespieDSolver_8h.xml

// File: HillSwitch_8cpp.xml

// File: HillSwitch_8h.xml

// File: LacHillSwitch_8cpp.xml

// File: LacHillSwitch_8h.xml

// File: NextReactionSolver_8cpp.xml

// File: NextReactionSolver_8h.xml

// File: SelfRegulatingGeneSwitch_8cpp.xml

// File: SelfRegulatingGeneSwitch_8h.xml

// File: TwoStateExpression_8cpp.xml

// File: TwoStateExpression_8h.xml

// File: TwoStateHillLoopSwitch_8cpp.xml

// File: TwoStateHillLoopSwitch_8h.xml

// File: TwoStateHillSwitch_8cpp.xml

// File: TwoStateHillSwitch_8h.xml

// File: CheckpointSignaler_8cpp.xml

// File: CheckpointSignaler_8h.xml

// File: DataOutputQueue_8cpp.xml

// File: DataOutputQueue_8h.xml

// File: Exceptions_8h.xml

// File: Globals_8cpp.xml

// File: Globals_8h.xml

// File: LocalDataOutputWorker_8cpp.xml

// File: LocalDataOutputWorker_8h.xml

// File: Math_8h.xml

%feature("docstring") isPower2 "
`isPower2(x) -> bool`  
";

%feature("docstring") isPower2 "
`isPower2(x) -> bool`  
";

%feature("docstring") isPower2 "
`isPower2(x) -> bool`  
";

%feature("docstring") log2 "
`log2(x) -> unsigned int`  
";

%feature("docstring") log2 "
`log2(x) -> unsigned int`  
";

%feature("docstring") log2 "
`log2(x) -> unsigned int`  
";

// File: MPIRemoteDataOutputQueue_8cpp.xml

// File: MPIRemoteDataOutputQueue_8h.xml

// File: Print_8cpp.xml

// File: Print_8h.xml

// File: ReplicateRunner_8cpp.xml

// File: ReplicateRunner_8h.xml

// File: ResourceAllocator_8cpp.xml

// File: ResourceAllocator_8h.xml

// File: runSimulation_8cpp.xml

%feature("docstring") abortHandler "
`abortHandler(x)`  
";

%feature("docstring") runSimulation "
`runSimulation(simulationFilename, replicate, solverClass, cudaDevices, checkpointInterval)`  
";

%feature("docstring") runSolver "
`runSolver(simulationFilename, replicate, solver, cudaDevices, checkpointInterval)`  
";

// File: runSimulation_8h.xml

%feature("docstring") runSimulation "
`runSimulation(simulationFilename, replicate, solverClass, cudaDevices, checkpointInterval)`  
";

%feature("docstring") runSolver "
`runSolver(simulationFilename, replicate, solver, cudaDevices, checkpointInterval)`  
";

// File: SignalHandler_8cpp.xml

// File: SignalHandler_8h.xml

// File: Timer_8h.xml

// File: Types_8h.xml

// File: util_8cpp.xml

%feature("docstring") parseIndices "
`parseIndices(s, matrixDims) -> vector< int >`  
";

%feature("docstring") parseRange "
`parseRange(rangeStr, maxDim) -> pair< int, int >`  
";

%feature("docstring") parseValues "
`parseValues(s) -> vector< double >`  
";

%feature("docstring") printCopyright "
`printCopyright(argc, argv)`  

This function prints the copyright notice.  
";

%feature("docstring") printBuildConfig "
`printBuildConfig()`  
";

// File: util_8h.xml

%feature("docstring") parseIndices "
`parseIndices(s, matrixDims) -> vector< int >`  
";

%feature("docstring") parseRange "
`parseRange(arg1, maxDim) -> pair< int, int >`  
";

%feature("docstring") parseValues "
`parseValues(s) -> vector< double >`  
";

%feature("docstring") printCopyright "
`printCopyright(argc, argv)`  

This function prints the copyright notice.  
";

%feature("docstring") printBuildConfig "
`printBuildConfig()`  
";

// File: constant_8cu.xml

// File: constant_8cuh.xml

// File: ldg_8h.xml

// File: lm__cuda_8cu.xml

// File: lm__cuda_8h.xml

// File: ArbitraryH5_8h.xml

%feature("docstring") get_h5_type_id "
`get_h5_type_id() -> hid_t`  
";

%feature("docstring") get_h5_type_id< char > "
`get_h5_type_id< char >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< unsigned char > "
`get_h5_type_id< unsigned char >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< short > "
`get_h5_type_id< short >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< unsigned short > "
`get_h5_type_id< unsigned short >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< int > "
`get_h5_type_id< int >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< unsigned int > "
`get_h5_type_id< unsigned int >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< long > "
`get_h5_type_id< long >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< unsigned long > "
`get_h5_type_id< unsigned long >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< long long > "
`get_h5_type_id< long long >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< unsigned long long > "
`get_h5_type_id< unsigned long long >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< float > "
`get_h5_type_id< float >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< double > "
`get_h5_type_id< double >() -> hid_t`  
";

%feature("docstring") get_h5_type_id< long double > "
`get_h5_type_id< long double >() -> hid_t`  
";

%feature("docstring") make_H5_meta "
`make_H5_meta(mode, shape, name) -> H5MetaData`  
";

%feature("docstring") make_H5_lookup "
`make_H5_lookup(mode, path, attr=\"\") -> H5Lookup`  
";

// File: lm__hdf5_8h.xml

// File: SimulationFile_8cpp.xml

// File: SimulationFile_8h.xml

// File: SimulationFile__create_8cpp.xml

// File: SimulationParameters_8cpp.xml

// File: SimulationParameters_8h.xml

// File: Profile_8h.xml

// File: Profile__NVTX_8h.xml

// File: ProfileCodes_8h.xml

// File: MESolver_8cpp.xml

// File: MESolver_8h.xml

// File: MESolverFactory_8cpp.xml

// File: MESolverFactory_8h.xml

// File: lm__mpi_8cpp.xml

// File: lm__mpi_8h.xml

%feature("docstring") lm::MPIErrorHandler "
`MPIErrorHandler(arg1, rc, arg3)`  
";

// File: ByteLattice_8cpp.xml

// File: ByteLattice_8h.xml

// File: CudaByteLattice_8cu.xml

// File: CudaByteLattice_8h.xml

// File: CudaIntLattice_8cu.xml

// File: CudaIntLattice_8h.xml

// File: bit__packed__diffusion__1d__dev_8cu.xml

// File: byte__diffusion__1d__dev_8cu.xml

// File: byte__diffusion__1d__mgpu__dev_8cu.xml

// File: byte__diffusion__2d__dev_8cu.xml

// File: byte__diffusion__3dw__dev_8cu.xml

// File: byte__reaction__dev_8cu.xml

// File: lattice__sim__1d__dev_8cu.xml

// File: lattice__sim__2d__dev_8cu.xml

// File: lattice__sim__3d__dev_8cu.xml

// File: lattice__sim__3dw__dev_8cu.xml

// File: word__diffusion__1d__dev_8cu.xml

// File: word__reaction__dev_8cu.xml

// File: xor__random__dev_8cu.xml

// File: MultiGPUMapper_8cu.xml

// File: MultiGPUMapper_8h.xml

// File: osx__barrier_8cu.xml

// File: osx__barrier_8h.xml

// File: SegmentDescriptor_8h.xml

// File: ZDivMPIGPUMapper_8cu.xml

// File: ZDivMPIGPUMapper_8h.xml

// File: ZDivMultiGPUMapper_8cu.xml

// File: ZDivMultiGPUMapper_8h.xml

// File: IntLattice_8cpp.xml

// File: IntLattice_8h.xml

// File: IntMpdRdmeSolver_8cu.xml

// File: IntMpdRdmeSolver_8h.xml

// File: Lattice_8cpp.xml

// File: Lattice_8h.xml

// File: MGPUMpdRdmeSolver_8cu.xml

// File: MGPUMpdRdmeSolver_8h.xml

// File: MpdRdmeSolver_8cu.xml

// File: MpdRdmeSolver_8h.xml

// File: MpdTestHarness_8cu.xml

// File: MpdTestHarness_8h.xml

// File: MPIMpdRdmeSolver_8cu.xml

// File: MPIMpdRdmeSolver_8h.xml

// File: NextSubvolumeSolver_8cpp.xml

// File: NextSubvolumeSolver_8h.xml

// File: RDMESolver_8cpp.xml

// File: RDMESolver_8h.xml

// File: ReactionQueue_8h.xml

// File: RandomGenerator_8cpp.xml

// File: RandomGenerator_8h.xml

// File: XORShift_8cpp.xml

// File: XORShift_8h.xml

// File: XORWow_8cu.xml

// File: XORWow_8h.xml

// File: Thread_8cpp.xml

// File: Thread_8h.xml

// File: Worker_8cpp.xml

// File: Worker_8h.xml

// File: WorkerManager_8cpp.xml

// File: WorkerManager_8h.xml

// File: dir_98210ece4f8a652fe0a2321bd18e8fca.xml

// File: dir_50a28221d0f5027093306cfc52a57d36.xml

// File: dir_ebb70c7f998c50b7923b52bdf51b4230.xml

// File: dir_4270bfced15e0e73154b13468c7c9ad9.xml

// File: dir_7f72761faf76ad5b97d7f36c114afa29.xml

// File: dir_2a518cd58ab2be4f554f5d42f7b4d2e1.xml

// File: dir_b105aaf15c5b45ddf8775d85d576d4a8.xml

// File: dir_bc161955dc3a3d2485839eba21420d01.xml

// File: dir_1ea1ebb12869271316950dbce045adaf.xml

// File: dir_b4974efccefa95895ac72668784e5964.xml

// File: dir_13620b7f1330796542fcb7e43593aaf8.xml

// File: dir_b34cadae8af9b4163446ad9749a9a622.xml

// File: dir_d0cc35b808c27df60b208baf27098539.xml

// File: dir_a2f2002add987e269fffac77984b5ad9.xml

// File: dir_01284e59d658032137ac90170bc51d5c.xml

