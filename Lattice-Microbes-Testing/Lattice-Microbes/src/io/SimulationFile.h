/*
 * University of Illinois Open Source License
 * Copyright 2010-2018 Luthey-Schulten Group,
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

#ifndef LM_IO_HDF5_SIMULATIONFILE_H_
#define LM_IO_HDF5_SIMULATIONFILE_H_

#include <string>
#include <map>
#include <vector>
#include "core/Exceptions.h"
#include "core/Types.h"
#include "io/lm_hdf5.h"
#include "io/ArbitraryH5.h"

namespace lm {

namespace rdme{
class Lattice;
class ByteLattice;
class IntLattice;
}

namespace io {

class DiffusionModel;
class ReactionModel;
class Lattice;
class ParameterValues;
class SpeciesCounts;
class SpatialModel;
class FirstPassageTimes;

namespace hdf5 {

using std::string;
using std::map;
using std::vector;
using lm::IOException;

/// @class SimulationFile
/// @brief A representation of the simulation that is used to input or output from an HDF5 file
class SimulationFile
{
// Static variables and functions
public:
    static const uint MIN_VERSION;
    static const uint CURRENT_VERSION;
    static const uint MAX_REACTION_RATE_CONSTANTS;
    static const uint MAX_SHAPE_PARAMETERS;

public:
    /// @brief Check that the HDF5 file is valid
    static bool isValidFile(const string filename);
    /// @brief Check that the HDF5 file is valid
    static bool isValidFile(const char * filename);
    /// @brief Create an HDF5 file
    static void create(const string filename);
    /// @brief Create an HDF5 file
    static void create(const char *  filename);
    /// @brief Create an HDF5 file with the specified number of species
    static void create(const string filename, unsigned int numberSpecies);
    /// @brief Create an HDF5 file with the specified number of species
    static void create(const char *  filename, unsigned int numberSpecies);
    /// @brief Create an HDF5 file and initialize the model with an optional number of species
    static void create(const char * filename, bool initializeModel, unsigned int numberSpecies=0);

protected:
    static herr_t parseParameter(hid_t location_id, const char *attr_name, const H5A_info_t *ainfo, void *op_data);
	void writeParticleLattice_U8LE(const char* path,  uint8_t *data, unsigned int x, unsigned int y, unsigned int z, unsigned int w);
	void writeParticleLattice_U32LE(hid_t handle, const char* path,  unsigned int *data, unsigned int x, unsigned int y, unsigned int z, unsigned int w);
	void write3DLattice_U8LE(const char* path,  uint8_t *data, unsigned int x, unsigned int y, unsigned int z);
    void appendToDataset(H5MetaData *header, byte *data);
    void makeNewDataset(H5MetaData *header, byte *data);
    void appendDsToGroup(H5MetaData *header, byte *data);
    void readDataset(H5Lookup *req, byte *&data, size_t &sz);
    void readAttr(H5Lookup *req, byte *&data, size_t &sz);
    
    
// Non static functions and variables
public:
    /// @brief Create a SimulationFile with the specified filename
    SimulationFile(const string filename);
    /// @brief Create a SimulationFile with the specified filename
    SimulationFile(const char* filename);
    /// @brief Destroy the SimulationFile
	virtual ~SimulationFile();
    
    /// @brief Close the file
    virtual void close();
    /// @brief Flush the current buffered data to the file
    virtual void flush();
    /// @brief Create the current checkpoint to the HDF5 file
    /// @param checkpointFilename Name of the new checkpoint file
    virtual string checkpoint();

    // Methods for working with parameters.
    /// @brief Get all the parameters from the current replicate
    virtual map<string,string> getParameters();
    /// @brief Get the specified parameter from the current replicate
    virtual string getParameter(string key, string defaultValue="");
    /// @brief Set the specified parameter to the value
    virtual void setParameter(string key, string value);

    // Methods for working with the model.
    /// @brief Pops the reaction model from the HDF5 file
    /// @param reactionModel Memory at which to place the reactionModel
    virtual void getReactionModel(lm::io::ReactionModel * reactionModel);
    /// @brief Pushes the reaction model into the HDF5 file
    /// @param reactionModel Memory from which to place the reactionModel
    virtual void setReactionModel(lm::io::ReactionModel * reactionModel);
    /// @brief Pops the diffusion model from the HDF5 file
    /// @param diffusionModel Memory at which to place the diffusionModel
    virtual void getDiffusionModel(lm::io::DiffusionModel * diffusionModel);
    /// @brief Pops the diffusion lattice from the HDF5 file
    /// @param diffusionModel Memoryfrom which to get the lattice
    /// @param lattice Memory in which to place lattice
    /// @param latticeMaxSize Total size of the lattice (i.e. number of bytes)
    /// @param latticeSites Actual lattice site data in byte format
    /// @param latticeSitesMaxSize Max size of a lattice site (i.e. number of particles it can contain)
    virtual void getDiffusionModelLattice(lm::io::DiffusionModel * diffusionModel, byte * lattice, size_t latticeMaxSize, byte * latticeSites, size_t latticeSitesMaxSize);
    /// @brief Pops the diffusion lattice from the HDF5 file
    /// @param diffusionModel Memoryfrom which to get the lattice
    /// @param lattice Memory at which to place the lattice
    virtual void getDiffusionModelLattice(lm::io::DiffusionModel * diffusionModel, lm::rdme::Lattice * lattice);
    /// @brief Pushes the diffusion model into the HDF5 file
    /// @param diffusionModel Memory from which to place the diffusionModel
    virtual void setDiffusionModel(lm::io::DiffusionModel * diffusionModel);

    /// @brief Pushes the diffusion lattice into the HDF5 file
    /// @param diffusionModel Memoryfrom which to get the lattice
    /// @param lattice Memory from which to place the lattice - uint8_t particle lattice
    /// @param latticeSites Memory from which to place lattice contents
    virtual void setDiffusionModelLattice(lm::io::DiffusionModel * m, uint8_t * lattice, uint8_t * latticeSites);

    /// @brief Pushes the diffusion lattice into the HDF5 file
    /// @param diffusionModel Memoryfrom which to get the lattice
    /// @param lattice Memory from which to place the lattice - uint32_t particle lattice
    /// @param latticeSites Memory from which to place lattice contents
    virtual void setDiffusionModelLattice(lm::io::DiffusionModel * m, uint32_t * lattice, uint8_t * latticeSites);

    /// @brief Pushes the diffusion lattice into the HDF5 file
    /// @param diffusionModel Memory from which to get the lattice
    /// @param lattice Memory from which to place the lattice
    virtual void setDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::ByteLattice * lattice);

    /// @brief Pushes the diffusion lattice into the HDF5 file
    /// @param diffusionModel Memory from which to get the lattice
    /// @param lattice Memory from which to place the lattice
    virtual void setDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::IntLattice * lattice);

    /// @brief Pushes the spacial model (i.e. obstacles) into the HDF5 file
    /// @param model The model object
    virtual void setSpatialModel(lm::io::SpatialModel * model);
    /// @brief Pops the spacial model (i.e. obstacles) from the HDF5 file
    /// @param model The model in which to store the object
    virtual void getSpatialModel(lm::io::SpatialModel * model);


    // Methods for working with a replicate.
    /// @brief Checks if the specified replicate exists
    /// @param replicate Replicate number
    /// @return true/false
    virtual bool replicateExists(unsigned int replicate);
    /// @brief Opens the specified replicate for reading
    /// @param replicate Replicate number
    virtual void openReplicate(unsigned int replicate);
    /// @brief Appends the species counts in the various sites into the replicate
    /// @param replicate Replicate number
    /// @param speciesCounts Protocol buffers object with species counts in lattice
    virtual void appendSpeciesCounts(unsigned int replicate, lm::io::SpeciesCounts * speciesCounts);
    /// @brief Appends the lattice to the current replicate
    /// @param replicate Replicate number
    /// @param lattice Lattice to add to replicate
    /// @param latticeData The actual data of the lattice
    /// @param latticeDataSize The size of the lattice
    virtual void appendLattice(unsigned int replicate, lm::io::Lattice * lattice, byte * latticeData, size_t latticeDataSize);
    virtual void appendLattice_U32LE(unsigned int replicate, lm::io::Lattice * lattice, uint32_t * latticeData, size_t latticeDataSize);

    virtual void arbitraryH5(byte *ptr);
    virtual void arbitraryH5Lookup(byte *in, byte *&out, size_t &out_sz);
	virtual void appendSites(unsigned int replicate, lm::io::Lattice * lattice, uint8_t * siteData, size_t siteDataSize);
    /// @brief Adds all the parameter values to the replicate
    /// @param replicate Replicate number
    /// @param parameterValues Actual values to store
    virtual void appendParameterValues(unsigned int replicate, lm::io::ParameterValues * parameterValues);
    /// @brief Adds all the first passage times to the replicate
    /// @param replicate Replicate number
    /// @param speciesCounts Actual passage times to store
    virtual void setFirstPassageTimes(unsigned int replicate, lm::io::FirstPassageTimes * speciesCounts);
    /// @brief Get the timestep times for the replicate
    /// @param replicate Replicate number
    /// @return A set of timestep times
    virtual vector<double> getLatticeTimes(unsigned int replicate);
    /// @brief Get the lattice for the replicate
    /// @param replicate Replicate number
    /// @param latticeIndex Seems to be unused...
    /// @param lattice The lattice object in which to store the data
    virtual void getLattice(unsigned int replicate, unsigned int latticeIndex, lm::rdme::Lattice * lattice);
    /// @brief Close the specified replicate
    /// @param replicate Replicate number
    virtual void closeReplicate(unsigned int replicate);
    /// @brief Close all the replicates
    virtual void closeAllReplicates();

    //virtual void appendSpatialModelObjects(unsigned int replicate, lm::io::SpatialModel * model);
    //virtual void getSpatialModelObjects(unsigned int replicate, lm::io::SpatialModel * model);

	/*
    virtual lattice_coord_t getLatticeSize() const;
	virtual nmdist_t getLatticeSpacing() const;
	virtual uint getMaxParticlesPerSite() const;
	virtual lattice_particle_t getMaxParticleType() const;
	virtual lattice_site_t getMaxSiteType() const;
	virtual const std::map<uint64,uint64> getMaxParticleCounts() const;
	virtual const std::map<uint64,uint64> getMaxSiteCounts() const;
	
	virtual uint64 getNumberFrames() const;	
	virtual const std::vector<nstime_t> getFrameTimes() const;
	virtual void loadFrame(uint64 frameIndex, Lattice* lattice, nstime_t* time=NULL) const;
	
	virtual uint64 getNumberLatticeConfigurations() const;
	virtual const std::vector<nstime_t> getLatticeConfigurationTimes() const;
	virtual void loadLatticeConfiguration(uint64 latticeIndex, Lattice* lattice, nstime_t* time=NULL) const;
    */

	std::map<uint32_t, uint> getParticleCounts(unsigned int replicate, unsigned int latticeIndex);
	std::map<double, vector<int> > getSpeciesCounts(unsigned int replicate);
	vector<double> getSpeciesCountTimes(unsigned int replicate);
	std::map<uint, string> getSpeciesNames();
	std::map<uint, string> getSiteTypeNames();

	
public:

    /// @struct ReplicateHandles
    /// @brief A handle for the different replicates that may be stored in the HDF5 file
    struct ReplicateHandles
    {
        hid_t group;                    // ID of the group of replicates
        hid_t speciesCountsDataset;     // ?
        hid_t speciesCountTimesDataset; // ?
        
        /// @brief Set up the ReplicateHandles
        ReplicateHandles():group(H5I_INVALID_HID),speciesCountsDataset(H5I_INVALID_HID),speciesCountTimesDataset(H5I_INVALID_HID) {}
    };

	
protected:
    // Protected functions
    virtual void open();
    virtual void openGroups();
    virtual void loadParameters();
    virtual void loadModel();
    virtual ReplicateHandles * openReplicateHandles(unsigned int replicate);
    virtual ReplicateHandles * createReplicateHandles(string replicateString);
    virtual void closeReplicateHandles(ReplicateHandles * handles);
	
protected:
    string          filename;       // HDF5 file name
    hid_t           file;           // file handle
    unsigned int    version;        // HDF5 file version

    // Main group handles.
    hid_t           parametersGroup, modelGroup, simulationsGroup;

    // The parameters.
    map<string,string> parameterMap;

    // The model.
    bool            modelLoaded;
    unsigned int    numberSpecies;

    // Handles for each replicate that is open.
    map<unsigned int,ReplicateHandles *> openReplicates;

};

}
}
}

#endif
