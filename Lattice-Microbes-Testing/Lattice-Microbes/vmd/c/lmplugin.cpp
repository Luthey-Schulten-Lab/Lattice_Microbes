/*
 * University of Illinois Open Source License
 * Copyright 2008-2011 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 * of the Software, and to permit persons to whom the Software is furnished to 
 * do so, subject to the following conditions:
 * 
 * - Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimers.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimers in the documentation 
 * and/or other materials provided with the distribution.
 * 
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts
 */

#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <vector>
#include <stdexcept>
#include "molfile_plugin.h"
#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/Version.h"
//#include "lm/Print.h"
#include "lm/builder/Shape.h"
#include "lm/io/DiffusionModel.pb.h"
#include "lm/io/SpatialModel.pb.h"
#include "lm/io/hdf5/HDF5.h"
#include "lm/io/hdf5/SimulationFile.h"
#include "lm/rdme/Lattice.h"
#include "lm/rdme/ByteLattice.h"
#include <sys/time.h>

using std::map;
using std::vector;
using lm::Exception;
using lm::io::DiffusionModel;
using lm::io::SpatialModel;
using lm::io::hdf5::HDF5Exception;
using lm::io::hdf5::SimulationFile;
using lm::rdme::ByteLattice;

#define USDIFF(a,b) ((b.tv_sec*1000000 + b.tv_usec) - (a.tv_sec*1000000 + a.tv_usec))
struct timer
{
	float acc;
	struct timeval start;
	struct timeval end;
};

void init_timer(struct timer *t)
{
	t->acc=0.0f;
}
void start_timer(struct timer *t)
{
	gettimeofday(&t->start,NULL);
}
void stop_timer(struct timer *t)
{
	gettimeofday(&t->end,NULL);
	t->acc += USDIFF(t->start, t->end);
}
void report_timer(struct timer *t, const char *str)
{
	if(t->acc < 1e6)
		printf(" *** Time for %s: %g us\n", str, t->acc);	
	else
		printf(" *** Time for %s %g S\n", str, t->acc / 1000000.0f);
}

struct timer timer_tsload;
struct timer timer_getlattice;
struct timer timer_movevalid;
struct timer timer_moveinvalid;

struct lmobjects
{
    lmobjects() :
      file(NULL)
     ,replicate(1)
     ,numberFrames(0)
     ,nextFrame(0)
     ,usingObstacleAtoms(false)
     ,usingSiteAtoms(false)
     ,usingFractionalSites(true)
     ,lattice(NULL)
     ,particleStartAtom(NULL)
     ,siteStartAtom(0)
     ,spatialModelGeometryCount(0)
     ,spatialModelGeometry(NULL)
	{
		bond_from[0] = 0;
		bond_to[0] = 0;
	}

    ~lmobjects()
    {
        if (file != NULL) {delete file; file=NULL;}
        if (lattice != NULL) {delete lattice; lattice=NULL;}
        if (spatialModelGeometry != NULL) {delete spatialModelGeometry; spatialModelGeometry=NULL;}
        if (volumetricMetadata != NULL) {delete volumetricMetadata;}
		if(particleStartAtom!=NULL) free(particleStartAtom);
    }

    int bond_from[1];
    int bond_to[1];
    SimulationFile * file;
    uint replicate;
    uint numberFrames;
    vector<si_time_t> frameTimes;
    uint nextFrame;
    bool usingObstacleAtoms;
    bool usingSiteAtoms;
    bool usingFractionalSites;
    map<particle_t,uint> maxParticleCounts;
    DiffusionModel diffusionModel;
    SpatialModel spatialModel;
    ByteLattice * lattice;
	uint* particleStartAtom;
	uint obstacleStartAtom;
	uint siteStartAtom;
    int spatialModelGeometryCount;
	molfile_graphics_t * spatialModelGeometry;
	vector<site_t> volumetricSiteTypes;
	vector<std::string> volumetricNames;
	molfile_volumetric_t* volumetricMetadata;

};

static map<particle_t,uint> lm_get_max_particle_counts(lmobjects* objects)
{
	struct timer t; init_timer(&t);
	start_timer(&t);
	map<particle_t,uint> maxCounts;

	if (objects->file->replicateExists(objects->replicate))
	{
		vector<double> latticeTimes=objects->file->getLatticeTimes(objects->replicate);
		std::map<double, vector<int> >speciesCounts=objects->file->getSpeciesCounts(objects->replicate);
		
		for (uint i=0; i<objects->numberFrames; i++)
		{
			//printf("Scanning frame %d (time %g) for max particles\n", i, latticeTimes[i]);

			try
			{
				vector<int> counts=speciesCounts.at(latticeTimes[i]);
				for(uint i=0; i<counts.size(); i++)
				{
					if((uint)counts[i] > maxCounts[i])
						maxCounts[i] = counts[i];
				}
					
			}
			catch(const std::out_of_range &e)
			{
				printf("WARNING: No speciesCount value at time %g.  Reverting to slow count\n", latticeTimes[i]);
				map<particle_t,uint> counts = objects->file->getParticleCounts(objects->replicate, i);
				for (map<particle_t,uint>::iterator it=counts.begin(); it != counts.end(); it++)
				{
					particle_t key = it->first;
					if (maxCounts.count(key) == 0 || counts[key] > maxCounts[key])
					{
						maxCounts[key] = counts[key];
					}
				}
			}
		}
	}
	else
	{
		// Just load the initial lattice from the model.
		objects->file->getDiffusionModelLattice(&objects->diffusionModel, objects->lattice);
		maxCounts = objects->lattice->getParticleCounts();
	}

	stop_timer(&t);
	report_timer(&t, "lm_get_max_particle_counts");
	return maxCounts;
}

static void *lm_open_read(const char *filename, const char *filetype, int *natoms)
{	
	//if (getenv("LM_VERBOSITY") != NULL && atoi(getenv("LM_VERBOSITY")) > 0) 
	//	lm::Print::verbosityLevel(atoi(getenv("LM_VERBOSITY")));

	struct timer t; init_timer(&t);
	start_timer(&t);
	
	// Allocate an structure to store our data between callbacks.
	lmobjects * objects = new lmobjects();
	
	// init long-range timers
	init_timer(&timer_tsload);
	init_timer(&timer_getlattice);
	init_timer(&timer_movevalid);
	init_timer(&timer_moveinvalid);

	try
	{

		// Get the replicate to open.
		if (getenv("LM_REPLICATE") != NULL) objects->replicate = atoi(getenv("LM_REPLICATE"));

		// See if we should use obstacle atoms.
		if (getenv("LM_CREATE_OBSTACLE_ATOMS") != NULL && atoi(getenv("LM_CREATE_OBSTACLE_ATOMS")) > 0) objects->usingObstacleAtoms = true;

		// See if we should use site atoms.
		if (getenv("LM_CREATE_SITE_ATOMS") != NULL && atoi(getenv("LM_CREATE_SITE_ATOMS")) > 0) objects->usingSiteAtoms = true;

        // Open the file.
        objects->file = new SimulationFile(filename);
		printf("LMplugin Info) Opened file: %s, replicate %d, using options LM_CREATE_OBSTACLE_ATOMS=%d, LM_CREATE_SITE_ATOMS=%d.\n", filename, objects->replicate, objects->usingObstacleAtoms?1:0, objects->usingSiteAtoms?1:0);

        // Read the diffusion model.
        objects->file->getDiffusionModel(&objects->diffusionModel);
		printf("LMplugin Info) Lattice size: %u %u %u, particles per site %u, spacing: %e m.\n", objects->diffusionModel.lattice_x_size(), objects->diffusionModel.lattice_y_size(), objects->diffusionModel.lattice_z_size(), objects->diffusionModel.particles_per_site(), objects->diffusionModel.lattice_spacing());

        // Read the spatial model.
        objects->file->getSpatialModel(&objects->spatialModel);

		// Allocate a new lattice.
		objects->lattice = new ByteLattice(objects->diffusionModel.lattice_x_size(), objects->diffusionModel.lattice_y_size(), objects->diffusionModel.lattice_z_size(), objects->diffusionModel.lattice_spacing(), objects->diffusionModel.particles_per_site());

        // Reset the atom count to zero.
        *natoms = 0;

		// Allocate atom count structures
		objects->particleStartAtom=(uint*)calloc(objects->lattice->getMaxParticle()+1, sizeof(uint));

		// See if this is a valid replicate.
		if (objects->file->replicateExists(objects->replicate))
		{
			// Figure out how many frames there are.
			objects->frameTimes = objects->file->getLatticeTimes(objects->replicate);
			objects->numberFrames = objects->frameTimes.size();

		}
		else
		{
			// We need one frame to display the model.
			objects->frameTimes.push_back(0.0);
			objects->numberFrames = objects->frameTimes.size();
		}

		// Figure out how the maximum number of particles we have in any one frame.
		uint nParticleAtoms=0;
		objects->maxParticleCounts = lm_get_max_particle_counts(objects);
		for (map<particle_t,uint>::const_iterator it=objects->maxParticleCounts.begin(); it != objects->maxParticleCounts.end(); it++)
			nParticleAtoms += it->second;
		printf("LMplugin Info) Using %d particle atoms.\n", nParticleAtoms);

		// Figure out how many atoms we need for obstacles.
		uint nObstacleAtoms=0;
		if (objects->usingObstacleAtoms)
		{
			for (uint i=0; i<(uint)objects->spatialModel.obstacle_size(); i++)
			{
				if (objects->spatialModel.obstacle(i).shape() == lm::builder::Shape::SPHERE)
				{
					nObstacleAtoms++;
				}
			}
		}
	    printf("LMplugin Info) Using %d obstacle atoms.\n", nObstacleAtoms);

		// Figure out how the maximum number of sites we have in any one frame.
		uint nSiteAtoms=0;
		if (objects->usingSiteAtoms)
		{
			nSiteAtoms = objects->diffusionModel.lattice_x_size()*objects->diffusionModel.lattice_y_size()*objects->diffusionModel.lattice_z_size();
		    printf("LMplugin Info) Using %d site atoms.\n", nSiteAtoms);
		}

		// Set the total number of atoms we need.
		*natoms = nParticleAtoms + nObstacleAtoms + nSiteAtoms;

		std::map<uint, std::string> names=objects->file->getSiteTypeNames();
		// Figure out the which site types should have volumetric data laoded.
		for (size_t i=0; i < objects->diffusionModel.number_site_types(); i++)
		{
				objects->volumetricSiteTypes.push_back(i);
				std::string name;

			try
			{
				name=names.at(i);
			}
			catch(const std::out_of_range &e)
			{
				name="lattice-site-";
				name += (char)('0'+i);
			}

			objects->volumetricNames.push_back(name);
		}

		// Create the volumetric metadata objects.
		objects->volumetricMetadata = (molfile_volumetric_t*)malloc(sizeof(molfile_volumetric_t)*objects->volumetricSiteTypes.size());
		if (objects->volumetricMetadata == NULL) throw Exception("Unable to allocate lkobjects structure.");
		for (uint i=0; i<objects->volumetricSiteTypes.size(); i++)
		{

			strncpy(objects->volumetricMetadata[i].dataname, objects->volumetricNames[i].c_str(),255);
			objects->volumetricMetadata[i].dataname[255]='\0';
			int xd=objects->diffusionModel.lattice_x_size()+2;
			int yd=objects->diffusionModel.lattice_y_size()+2;
			int zd=objects->diffusionModel.lattice_z_size()+2;
			float lambdae10=objects->lattice->getSpacing() * 1e10;
			objects->volumetricMetadata[i].origin[0]=-1.5*lambdae10;
			objects->volumetricMetadata[i].origin[1]=-1.5*lambdae10;
			objects->volumetricMetadata[i].origin[2]=-1.5*lambdae10;
			objects->volumetricMetadata[i].xaxis[0]=xd*lambdae10;
			objects->volumetricMetadata[i].xaxis[1]=0.0;
			objects->volumetricMetadata[i].xaxis[2]=0.0;
			objects->volumetricMetadata[i].yaxis[0]=0.0;
			objects->volumetricMetadata[i].yaxis[1]=yd*lambdae10;
			objects->volumetricMetadata[i].yaxis[2]=0.0;
			objects->volumetricMetadata[i].zaxis[0]=0.0;
			objects->volumetricMetadata[i].zaxis[1]=0.0;
			objects->volumetricMetadata[i].zaxis[2]=zd*lambdae10;
			objects->volumetricMetadata[i].xsize=xd;
			objects->volumetricMetadata[i].ysize=yd;
			objects->volumetricMetadata[i].zsize=zd;
			objects->volumetricMetadata[i].has_color=0;

		}
	}
    catch (HDF5Exception & e)
    {
        printf("LMplugin Error) Caught HDF5 exception during read: %s\n", e.what());
        e.printStackTrace();
		delete objects; objects=NULL;
    }
    catch (Exception & e)
    {
		printf("LMplugin Error) Caught exception during read: %s\n", e.what());
		delete objects; objects=NULL;
    }
	catch(std::exception& e)
	{
		printf("LMplugin Error) Caught exception during read: %s\n", e.what());
		delete objects; objects=NULL;
	}
	catch(...)
	{
		printf("LMplugin Error) Caught exception during read.\n");
		delete objects; objects=NULL;
	}

	stop_timer(&t);
	report_timer(&t, "lm_open_read");
	return objects;;
}

static int lm_read_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
	struct timer t; init_timer(&t);
	start_timer(&t);
	try
	{
		lmobjects* objects = (lmobjects*)mydata;
	
		// Get the species Names
		std::map<uint, std::string> names=objects->file->getSpeciesNames();

		// Index to track the next atom to create.
		int atomIndex=0;

		// Create the atoms for the particles.
		for (map<particle_t,uint>::const_iterator it=objects->maxParticleCounts.begin(); it != objects->maxParticleCounts.end(); it++)
		{
			particle_t type = it->first;
			uint count = it->second;

			std::string name;
			char cname[16];
			try
			{
				name=names.at(type);
			}
			catch(const std::out_of_range &e)
			{
				name="particle";
			}

			if(name.length() > 15)
			{
				snprintf(cname, 16, "%s~%d", name.substr(0,10).c_str(), type);
				printf("WARNING: particle type %d name is truncated from %s to %s\n", type, name.c_str(), cname);
			}
			else
				strncpy(cname, name.c_str(),16);


			// Record the starting atom index for particles of this type.
			objects->particleStartAtom[type] = atomIndex;

			// Create one atom for each possible particle of this type.
			for (uint i=0; i<count; i++)
			{
				strncpy(atoms[atomIndex].name, cname, 16);
				snprintf(atoms[atomIndex].type, 16, "%d", type);
				snprintf(atoms[atomIndex].resname, 8, "%d", type);
				atoms[atomIndex].resid = type;
				strcpy(atoms[atomIndex].segid, "");
				strcpy(atoms[atomIndex].chain, "");
				if (objects->usingFractionalSites)
					atoms[atomIndex].radius = ((float)objects->lattice->getSpacing())*1e10/4.0;
				else
					atoms[atomIndex].radius = ((float)objects->lattice->getSpacing())*1e10/2.0;
				atomIndex++;
			}
		}

		// Create atoms for the obstacles.
		if (objects->usingObstacleAtoms)
		{
			objects->obstacleStartAtom = atomIndex;
			for (uint i=0; i<(uint)objects->spatialModel.obstacle_size(); i++)
			{
				if (objects->spatialModel.obstacle(i).shape() == lm::builder::Shape::SPHERE)
				{
					strcpy(atoms[atomIndex].name, "obstacle");
					snprintf(atoms[atomIndex].type, 16, "%d", objects->spatialModel.obstacle(i).site_type());
					snprintf(atoms[atomIndex].resname, 8, "%d", objects->spatialModel.obstacle(i).site_type());
					atoms[atomIndex].resid = i;
					strcpy(atoms[atomIndex].segid, "");
					strcpy(atoms[atomIndex].chain, "");
					atoms[atomIndex].radius = (float)(objects->spatialModel.obstacle(i).shape_parameter(3)*1e10);
					atomIndex++;
				}
			}
		}

		// Load the diffusion model
		// this gets the site types populated in objects->lattice
		objects->file->getDiffusionModelLattice(&objects->diffusionModel, objects->lattice);

		// Create the atoms for the sites.
		if (objects->usingSiteAtoms)
		{
			objects->siteStartAtom = atomIndex;

			// Get the lattice configuration.
			// TODO: Use objexts->lattice instead of loading again
			size_t sitesSize = objects->diffusionModel.lattice_x_size()*objects->diffusionModel.lattice_y_size()*objects->diffusionModel.lattice_z_size();
			byte * siteTypes = new byte[sitesSize];
			objects->file->getDiffusionModelLattice(&objects->diffusionModel, NULL, 0, siteTypes, sitesSize);

			// Go through the sites and create the atoms.
			for (lattice_size_t x=0, index=0; x<objects->diffusionModel.lattice_x_size(); x++)
				for (lattice_size_t y=0; y<objects->diffusionModel.lattice_y_size(); y++)
					for (lattice_size_t z=0; z<objects->diffusionModel.lattice_z_size(); z++, index++)
					{
						strcpy(atoms[atomIndex].name, "site");
						snprintf(atoms[atomIndex].type, 16, "%d", siteTypes[index]);
						snprintf(atoms[atomIndex].resname, 8, "%d", siteTypes[index]);
						atoms[atomIndex].resid = siteTypes[index];
						strcpy(atoms[atomIndex].segid, "");
						strcpy(atoms[atomIndex].chain, "");
						atoms[atomIndex].radius = ((float)objects->lattice->getSpacing())*1e10/2.0;
						atomIndex++;
					}

			delete [] siteTypes;
		}
		*optflags = MOLFILE_RADIUS;
		
		stop_timer(&t);
		report_timer(&t, "lm_read_structure");
		return MOLFILE_SUCCESS;
	}
	catch (lm::Exception & e)
	{
		printf("LMplugin Error) Caught exception during structure create: %s\n", e.what());
	}
	catch(std::exception & e)
	{
		printf("LMplugin Error) Caught exception during structure create: %s\n", e.what());
	}
	catch(...)
	{
		printf("LMplugin Error) Caught exception during structure create.\n");
	}
	return MOLFILE_ERROR;
}

static int lm_read_volumetric_metadata(void *mydata, int *nsets, molfile_volumetric_t **metadata)
{
	lmobjects* objects = (lmobjects*)mydata;

	*nsets=objects->volumetricSiteTypes.size();
	*metadata=objects->volumetricMetadata;
	return MOLFILE_SUCCESS;
}

static int lm_read_volumetric_data(void *mydata, int set, float *datablock, float *colorblock)
{
	lmobjects* objects = (lmobjects*)mydata;
	site_t siteType = objects->volumetricSiteTypes[set];

	printf("LKplugin Info) Reading volumetric data set %d for site type %d.\n", set, (int)siteType);
	int datablockIndex=0;
	int count_type=0;
	int count_not=0;
	int count_border=0;

	for (lattice_size_t z=0; z<objects->lattice->getZSize()+2; z++)
    	for (lattice_size_t y=0; y<objects->lattice->getYSize()+2; y++)
        	for (lattice_size_t x=0; x<objects->lattice->getXSize()+2; x++,datablockIndex++)
        	{
        		if (x >= 1 && x <=objects->lattice->getXSize() && y >= 1 && y <=objects->lattice->getYSize() && z >= 1 && z <=objects->lattice->getZSize())
        		{
					if (objects->lattice->getSiteType(x-1, y-1, z-1) == siteType)
					{
						datablock[datablockIndex] = (float)1.0;
						count_type++;
					}
					else
					{
						datablock[datablockIndex] = (float)0.0;
						count_not++;
					}
				}
        		else
        		{
					datablock[datablockIndex] = (float)0.0;
					count_border++;
        		}
        	}
	//lm::Print::printf(lm::Print::DEBUG, "VolumetriSite type %d: voxels inside: %d, outside: %d, border: %d.\n", siteType, count_type, count_not, count_border);
	return MOLFILE_SUCCESS;
}

static int lm_read_bonds(void *mydata, int *nbonds, int **from, int **to, float **bondorder, int **bondtype, int *nbondtypes, char ***bondtypename)
{
	lmobjects* objects = (lmobjects*)mydata;

	*nbonds = 0;
	*from = objects->bond_from;
	*to = objects->bond_to;
	*bondorder = NULL;
	*bondtype = NULL;
	*nbondtypes = 0;
	*bondtypename = NULL;
	return MOLFILE_SUCCESS;
}

static int lm_read_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
{
	start_timer(&timer_tsload);
	try
	{
		lmobjects* objects = (lmobjects*)mydata;
		
		float lambdae10=objects->lattice->getSpacing() * 1e10;
		// TODO dont hardcode known number of particles per site
		float xfraction[8], yfraction[8], zfraction[8];
		for(int i=0; i<8; i++)
		{
			if (objects->usingFractionalSites)
			{
				xfraction[i]= ((i%2)==0)?0.35f:0.65f;
				yfraction[i]= ((i%4)<2)?0.35f:0.65f;
				zfraction[i]= (i<4)?0.35f:0.65f;
			}
			else
			{
				xfraction[i]=yfraction[i]=zfraction[i]=0.5f;
			}
		}

		if (objects->nextFrame < objects->numberFrames)
		{
			if (ts != NULL)
			{
				// See if the replicate exists.
				if (objects->file->replicateExists(objects->replicate))
				{
					// Load the frame.
					start_timer(&timer_getlattice);
					objects->file->getLattice(objects->replicate, objects->nextFrame, objects->lattice);
					stop_timer(&timer_getlattice);
				}
				else
				{
					// Load the model.
					objects->file->getDiffusionModelLattice(&objects->diffusionModel, objects->lattice);
				}

				// Map for tracking how many atoms of each particle type have been used.
				// TODO: remove hardcoded species count
				uint *particleAtomsUsed=(uint*)calloc(256, sizeof(uint)); // max species plus one (since 1-based in lattice)

				start_timer(&timer_movevalid);
				// Move the particle atoms to their proper position.
				lattice_size_t size_x=objects->lattice->getXSize();
				lattice_size_t size_y=objects->lattice->getYSize();
				lattice_size_t size_z=objects->lattice->getZSize();

				for (lattice_size_t z=0; z<size_z; z++)
				{
					for (lattice_size_t y=0; y<size_y; y++)
					{
						for (lattice_size_t x=0; x<size_x; x++)
						{
							uint particles = objects->lattice->getOccupancy(x,y,z);
							for (uint particleIndex=0; particleIndex < particles; particleIndex++)
							{
								// Get the type of particle.
								//particle_t particleType = objects->lattice->getParticleUnsafe(x, y, z, particleIndex);
								particle_t particleType = objects->lattice->getParticle(x, y, z, particleIndex);

								// Make sure we have not already moved too many atoms for this type.
#ifdef EXTRA_SANITY_CHECKS
								if (particleAtomsUsed[particleType] >= objects->maxParticleCounts[particleType])
									throw Exception("Too many particles of the given type.");
#endif

								// Get the next atom index for this type of particle.
								uint atomIndex = objects->particleStartAtom[particleType]+particleAtomsUsed[particleType];
								particleAtomsUsed[particleType]++;

								// Move the atom.
								ts->coords[atomIndex*3+0] = (((float)x)+xfraction[particleIndex])*lambdae10;
								ts->coords[atomIndex*3+1] = (((float)y)+yfraction[particleIndex])*lambdae10;
								ts->coords[atomIndex*3+2] = (((float)z)+zfraction[particleIndex])*lambdae10;
							}
						}
					}
				}

				stop_timer(&timer_movevalid);
				start_timer(&timer_moveinvalid);
				// Move any unused particle to -1, -1, -1.
				for(uint i=0; i <= objects->lattice->getMaxParticle(); i++)
				{
					particle_t particleType = i;
					uint startAtom = objects->particleStartAtom[i];
					uint usedAtoms = particleAtomsUsed[particleType];
					uint maxAtoms = objects->maxParticleCounts[particleType];
					if (usedAtoms < maxAtoms)
					{
						//printf("LMplugin Debug) Unused particle atoms of type %u: %d\n", particleType, maxAtoms-usedAtoms);
						for (uint i=usedAtoms; i<maxAtoms; i++)
						{
							uint atomIndex = startAtom+i;
							ts->coords[atomIndex*3+0] = -1.0f*lambdae10;
							ts->coords[atomIndex*3+1] = -1.0f*lambdae10;
							ts->coords[atomIndex*3+2] = -1.0f*lambdae10;
						}

					}
				}
				free(particleAtomsUsed);
				stop_timer(&timer_moveinvalid);
				
				// Move the obstacles to their proper positions.
				if (objects->usingObstacleAtoms)
				{
					uint atomIndex = objects->obstacleStartAtom;
					for (uint i=0; i<(uint)objects->spatialModel.obstacle_size(); i++)
					{
						if (objects->spatialModel.obstacle(i).shape() == lm::builder::Shape::SPHERE)
						{
							ts->coords[atomIndex*3+0] = (float)(objects->spatialModel.obstacle(i).shape_parameter(0)*1e10);
							ts->coords[atomIndex*3+1] = (float)(objects->spatialModel.obstacle(i).shape_parameter(1)*1e10);
							ts->coords[atomIndex*3+2] = (float)(objects->spatialModel.obstacle(i).shape_parameter(2)*1e10);
							atomIndex++;
						}
					}
				}

	            // Move the site atoms to their proper positions.
				if (objects->usingSiteAtoms)
				{
					uint atomIndex = objects->siteStartAtom;
					for (lattice_size_t x=0; x<objects->diffusionModel.lattice_x_size(); x++)
						for (lattice_size_t y=0; y<objects->diffusionModel.lattice_y_size(); y++)
							for (lattice_size_t z=0; z<objects->diffusionModel.lattice_z_size(); z++)
							{
								ts->coords[atomIndex*3+0] = (((float)x)+0.5f)*((float)objects->lattice->getSpacing()*1e10);
								ts->coords[atomIndex*3+1] = (((float)y)+0.5f)*((float)objects->lattice->getSpacing()*1e10);
								ts->coords[atomIndex*3+2] = (((float)z)+0.5f)*((float)objects->lattice->getSpacing()*1e10);
								atomIndex++;
							}
				}
				// printf("LMplugin Info) Loaded time step: %e s\n", objects->frameTimes[objects->nextFrame]);
			}
			objects->nextFrame++;
			stop_timer(&timer_tsload);
			return MOLFILE_SUCCESS;		
		}
		else
		{
			return MOLFILE_EOF;
		}
	}
	catch (lm::Exception & e)
	{
		printf("LMplugin Error) Caught exception during timestep read: %s\n", e.what());
	}
	catch(std::exception & e)
	{
		printf("LMplugin Error) Caught exception during timestep read: %s\n", e.what());
	}
	catch(...)
	{
		printf("LMplugin Error) Caught exception during timestep read.\n");
	}
    return MOLFILE_ERROR;
}

static int lm_read_rawgraphics(void *mydata, int *nelem, const molfile_graphics_t **data)
{
	lmobjects* objects = (lmobjects *)mydata;

    // Figure out how many geometry objects we need.
    objects->spatialModelGeometryCount = 0;
    for (uint i=0; i<(uint)objects->spatialModel.region_size(); i++)
    {
        if (objects->spatialModel.region(i).shape() == lm::builder::Shape::SPHERE)
        {
            objects->spatialModelGeometryCount += 1;
        }
        else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CYLINDER)
        {
            objects->spatialModelGeometryCount += 1;
        }
        else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CAPSULE)
        {
            objects->spatialModelGeometryCount += 3;
        }
        else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CAPSULE_SHELL)
        {
            objects->spatialModelGeometryCount += 3;
        }
        else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CUBOID)
        {
        }
        else
        {
            printf("LMplugin Error) Unknown spatial model shape type: %d\n", objects->spatialModel.region(i).shape());
        }
    }

    // Allocate memory for the geometry objects.
    objects->spatialModelGeometry = (molfile_graphics_t*)malloc(sizeof(molfile_graphics_t)*objects->spatialModelGeometryCount);

    // Fill in the geometry object table.
    int nextGeometryObject = 0;
    for (uint i=0; i<(uint)objects->spatialModel.region_size() && nextGeometryObject<objects->spatialModelGeometryCount; i++)
    {
        if (objects->spatialModel.region(i).shape() == lm::builder::Shape::SPHERE)
        {
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_SPHERE;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(0)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(1)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(2)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(3)*1e10;
            nextGeometryObject += 1;
        }
/*        else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CYLINDER)
        {
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_CYLINDER;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(0)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(1)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(2)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[3] = objects->spatialModel.region(i).shape_parameter(3)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[4] = objects->spatialModel.region(i).shape_parameter(4)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[5] = objects->spatialModel.region(i).shape_parameter(5)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(6)*1e10;
            nextGeometryObject += 1;
        }
*/      else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CAPSULE)
        {
        	printf("Creating cylinder with radius %e\n",objects->spatialModel.region(i).shape_parameter(6)*1e10);
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_CYLINDER;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(0)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(1)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(2)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[3] = objects->spatialModel.region(i).shape_parameter(3)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[4] = objects->spatialModel.region(i).shape_parameter(4)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[5] = objects->spatialModel.region(i).shape_parameter(5)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(6)*1e10;
            nextGeometryObject += 1;
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_SPHERE;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(0)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(1)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(2)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(6)*1e10;
            nextGeometryObject += 1;
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_SPHERE;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(3)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(4)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(5)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(6)*1e10;
            nextGeometryObject += 1;
        }
        else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CAPSULE_SHELL)
        {
        	printf("Creating cylinder shell with radius %e\n",objects->spatialModel.region(i).shape_parameter(6)*1e10);
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_CYLINDER;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(0)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(1)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(2)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[3] = objects->spatialModel.region(i).shape_parameter(3)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[4] = objects->spatialModel.region(i).shape_parameter(4)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[5] = objects->spatialModel.region(i).shape_parameter(5)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(6)*1e10;
            nextGeometryObject += 1;
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_SPHERE;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(0)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(1)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(2)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(6)*1e10;
            nextGeometryObject += 1;
            objects->spatialModelGeometry[nextGeometryObject].type = MOLFILE_SPHERE;
            objects->spatialModelGeometry[nextGeometryObject].style = 100;
            objects->spatialModelGeometry[nextGeometryObject].data[0] = objects->spatialModel.region(i).shape_parameter(3)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[1] = objects->spatialModel.region(i).shape_parameter(4)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].data[2] = objects->spatialModel.region(i).shape_parameter(5)*1e10;
            objects->spatialModelGeometry[nextGeometryObject].size = objects->spatialModel.region(i).shape_parameter(6)*1e10;
            nextGeometryObject += 1;
        }
        else if (objects->spatialModel.region(i).shape() == lm::builder::Shape::CUBOID)
        {
        }
		else
		{
			printf("LMplugin ERROR) Shape %d at index %d unsupported\n", objects->spatialModel.region(i).shape(), i);
		}
			
    }


    printf("LMplugin Info) Created %d geometry objects for the spatial model.\n", objects->spatialModelGeometryCount);

	// Return the data to the caller.
	*nelem = objects->spatialModelGeometryCount;
	*data = objects->spatialModelGeometry;
	return MOLFILE_SUCCESS;
}

static void lm_close_read(void *mydata) {
	report_timer(&timer_tsload, "total timestep load");
	report_timer(&timer_getlattice, "total lattice load");
	report_timer(&timer_movevalid, "atom move");
	report_timer(&timer_moveinvalid, "atom move inactive");
	lmobjects* objects = (lmobjects *)mydata;
	if (objects != NULL) delete objects;
}


/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
	printf("LMplugin Info) version %d.%d build %s\n", VERSION_NUM_MAJOR, VERSION_NUM_MINOR, BUILD_INFO);
	memset(&plugin, 0, sizeof(molfile_plugin_t));
	plugin.abiversion = vmdplugin_ABIVERSION;
	plugin.type = MOLFILE_PLUGIN_TYPE;
	plugin.name = "lm";
	plugin.prettyname = "Lattice Microbes";
	plugin.author = "Elijah Roberts";
	plugin.majorv = VERSION_NUM_MAJOR;
	plugin.minorv = VERSION_NUM_MINOR;
	plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
	plugin.filename_extension = "lm,lm5";
	plugin.open_file_read = lm_open_read;
	plugin.read_structure = lm_read_structure;
	plugin.read_bonds = lm_read_bonds;
	plugin.read_next_timestep = lm_read_timestep;
	plugin.close_file_read = lm_close_read;
	plugin.read_volumetric_metadata = lm_read_volumetric_metadata;
	plugin.read_volumetric_data = lm_read_volumetric_data;
	plugin.read_rawgraphics = lm_read_rawgraphics;
	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
	(*cb)(v, (vmdplugin_t *)&plugin);
	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
	return VMDPLUGIN_SUCCESS;
}
