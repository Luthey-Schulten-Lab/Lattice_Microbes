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

#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <sys/stat.h>
#include "config.h"
#include "core/Exceptions.h"
#include "core/Math.h"
#include "core/Print.h"
#include "io/ArbitraryH5.h"
#include "DiffusionModel.pb.h"
#include "FirstPassageTimes.pb.h"
#include "Lattice.pb.h"
#include "ParameterValues.pb.h"
#include "ReactionModel.pb.h"
#include "SpatialModel.pb.h"
#include "SpeciesCounts.pb.h"
#include "io/lm_hdf5.h"
#include "io/SimulationFile.h"
#include "rdme/Lattice.h"
#include "rdme/ByteLattice.h"
#include "rdme/IntLattice.h"

using std::list;
using std::map;
using std::string;
using std::vector;
using lm::IOException;

namespace lm {
namespace io {
namespace hdf5 {

const uint SimulationFile::MIN_VERSION                   = 2;
const uint SimulationFile::CURRENT_VERSION               = 4;
const uint SimulationFile::MAX_REACTION_RATE_CONSTANTS   = 10;
const uint SimulationFile::MAX_SHAPE_PARAMETERS		     = 10;


SimulationFile::SimulationFile(const string filename)
:filename(filename),file(H5I_INVALID_HID),version(0),parametersGroup(H5I_INVALID_HID),modelGroup(H5I_INVALID_HID),simulationsGroup(H5I_INVALID_HID),modelLoaded(false),numberSpecies(0)
{
	open();
}

SimulationFile::SimulationFile(const char* filename)
:filename(filename),file(H5I_INVALID_HID),version(0),parametersGroup(H5I_INVALID_HID),modelGroup(H5I_INVALID_HID),simulationsGroup(H5I_INVALID_HID),modelLoaded(false),numberSpecies(0)
{
	open();
}

void SimulationFile::open()
{
    // Make sure gzip is supported.
    unsigned int filter_info;
    if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) throw lm::Exception("The HDF5 library does not support gzip compression.");
    HDF5_EXCEPTION_CHECK(H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info));
    if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) || !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED)) throw lm::Exception("The HDF5 library does not support gzip filtering for both encoding and decoding.");

    // Turn of error printing.
    HDF5_EXCEPTION_CHECK(H5Eset_auto2(H5E_DEFAULT, NULL, NULL));

    // Make sure the file exists.
    struct stat fileStats;
    if (stat(filename.c_str(), &fileStats) != 0) throw lm::IOException("the specified file did not exist", filename.c_str());

    // Make sure it is a regular file.
    if (!S_ISREG(fileStats.st_mode)) throw lm::IOException("the specified file was not a regular file", filename.c_str());

    // Make sure the file has the correct magic and hdf headers.
    if (!isValidFile(filename)) throw lm::IOException("the specified file was not of the correct format", filename.c_str());

    // Open the file.
    HDF5_EXCEPTION_CALL(file,H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));

    // Get the size of the user block.
    hid_t creationProperties;
    hsize_t userblockSize;
    HDF5_EXCEPTION_CALL(creationProperties,H5Fget_create_plist(file));
    HDF5_EXCEPTION_CHECK(H5Pget_userblock(creationProperties, &userblockSize));
    if (userblockSize < 4) throw lm::IOException("the specified file did not have the correct user block size", userblockSize);

    // Make sure the version is supported.
    HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/", "formatVersion", &version));
    if (version < MIN_VERSION || version > CURRENT_VERSION) throw lm::IOException("the specified file format version is not supported", version);

    // Open the groups.
    openGroups();

    // Loaded the parameters.
    loadParameters();
}

void SimulationFile::openGroups()
{
    HDF5_EXCEPTION_CALL(parametersGroup,H5Gopen2(file, "/Parameters", H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(modelGroup,H5Gopen2(file, "/Model", H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(simulationsGroup,H5Gopen2(file, "/Simulations", H5P_DEFAULT));
}

SimulationFile::~SimulationFile()
{
	//Close the file, if it is still open.
	close();
}

void SimulationFile::flush()
{
    if (file != H5I_INVALID_HID)
    {
        lm::Print::printf(Print::DEBUG, "Flushing file %s.", filename.c_str());
        HDF5_EXCEPTION_CHECK(H5Fflush(file, H5F_SCOPE_GLOBAL));
    }
}

string SimulationFile::checkpoint()
{
    // Close the file.
    close();

    // Copy the file to a checkpoint file.
    string checkpointFilename(filename+".chk");
    FILE * out=NULL, * in=NULL;
    in=fopen(filename.c_str(), "rb");
    if (in == NULL) throw lm::IOException("the checkpoint input file could not be opened", filename.c_str());
    out=fopen(checkpointFilename.c_str(), "wb");
    if (out == NULL) throw lm::IOException("the checkpoint output file could not be opened", checkpointFilename.c_str());
    size_t bufferSize = 1024*1024;
    char * buffer = new char[bufferSize];
    while (!feof(in))
    {
        size_t bytesRead = fread(buffer, sizeof(char), bufferSize, in);
        if (bytesRead != bufferSize && ferror(in)) throw lm::IOException("the checkpoint input file could not be read", filename.c_str());
        size_t bytesWritten = fwrite(buffer, sizeof(char), bytesRead, out);
        if (bytesWritten != bytesRead) throw lm::IOException("the checkpoint input file could not be read", checkpointFilename.c_str());
    }
    if (fclose(in) != 0)  throw lm::IOException("the checkpoint input file could not be closed", filename.c_str());
    if (fclose(out) != 0)  throw lm::IOException("the checkpoint output file could not be closed", checkpointFilename.c_str());
    delete[] buffer;

    // Open the file again.
    open();

    return checkpointFilename;
}

void SimulationFile::close()
{
    // Close any open replicate groups.
    closeAllReplicates();

    // Close any open groups.
    if (parametersGroup != H5I_INVALID_HID)
    {
        HDF5_EXCEPTION_CHECK(H5Gclose(parametersGroup));
        parametersGroup = H5I_INVALID_HID;
    }
    if (modelGroup != H5I_INVALID_HID)
    {
        HDF5_EXCEPTION_CHECK(H5Gclose(modelGroup));
        modelGroup = H5I_INVALID_HID;
    }
    if (simulationsGroup != H5I_INVALID_HID)
    {
        HDF5_EXCEPTION_CHECK(H5Gclose(simulationsGroup));
        simulationsGroup = H5I_INVALID_HID;
    }

    // Close the file.
    if (file != H5I_INVALID_HID)
    {
        lm::Print::printf(Print::DEBUG, "Closing file %s, %d open objects remaining.", filename.c_str(), H5Fget_obj_count(file, H5F_OBJ_ALL)-1);
        HDF5_EXCEPTION_CHECK(H5Fclose(file));
        file = H5I_INVALID_HID;
    }
}

bool SimulationFile::isValidFile(const char * filename)
{
    return isValidFile(string(filename));
}

bool SimulationFile::isValidFile(const string filename)
{
    // Make sure the file has the right magic.
    FILE * fp = NULL;
    fp=fopen(filename.c_str(), "r");
    if (fp == NULL) throw lm::IOException("the specified file could not be opened", filename.c_str());
    char magic[5];
    magic[4] = '\0';
    for (int i=0; i<4; i++)
    {
        magic[i] = fgetc(fp);
        if (magic[i] == EOF)
        {
            magic[i] = '\0';
            break;
        }
    }
    if (fclose(fp) != 0)  throw lm::IOException("the specified file could not be closed", filename.c_str());
    if (strncmp(magic, "LMH5", 4) != 0) return false;

    // Make sure it is an hdf5 file.
    int isHdf5 = false;
    HDF5_EXCEPTION_CALL(isHdf5,H5Fis_hdf5(filename.c_str()));
    if (!isHdf5) return false;

    return true;

}

void SimulationFile::loadParameters()
{
    hsize_t n=0;
    HDF5_EXCEPTION_CHECK(H5Aiterate2(parametersGroup, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, &n, &SimulationFile::parseParameter, this));
}

herr_t SimulationFile::parseParameter(hid_t location_id, const char *attr_name, const H5A_info_t *ainfo, void *op_data)
{
    SimulationFile * file = reinterpret_cast<SimulationFile *>(op_data);

    // Open the attribute.
    hid_t attr, type;
    H5T_class_t typeClass;
    hsize_t size;
    HDF5_EXCEPTION_CALL(attr,H5Aopen(location_id, attr_name, H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(type,H5Aget_type(attr));
    typeClass=H5Tget_class(type);
    HDF5_EXCEPTION_CALL(size,H5Aget_storage_size(attr));
    if (file->version == 2 && typeClass==H5T_FLOAT && size == sizeof(double))
    {
        double value;
        char buffer[33];
        HDF5_EXCEPTION_CHECK(H5Aread(attr, H5T_NATIVE_DOUBLE, &value));
        snprintf(buffer, sizeof(buffer), "%g", value);
        file->parameterMap[attr_name] = buffer;
    }
    else if (file->version >= 3 && typeClass==H5T_STRING)
    {
        // Get the dataspace.
        hid_t space;
        HDF5_EXCEPTION_CALL(space,H5Aget_space(attr));

        // Create the memory datatype.
        hid_t memtype;
        HDF5_EXCEPTION_CALL(memtype,H5Tcopy(H5T_C_S1));
        HDF5_EXCEPTION_CHECK(H5Tset_size(memtype, size));

        // Read the data.
        char * value = new char[size];
        HDF5_EXCEPTION_CHECK(H5Aread(attr, memtype, value));

        // Add the parameter to the map.
        file->parameterMap[attr_name] = value;

        // Reclaim the memory.
        delete [] value;
        HDF5_EXCEPTION_CHECK(H5Sclose(space));
        HDF5_EXCEPTION_CHECK(H5Tclose(memtype));
    }
    HDF5_EXCEPTION_CHECK(H5Tclose(type));
    HDF5_EXCEPTION_CHECK(H5Aclose(attr));
    return 0;
}

map<string,string> SimulationFile::getParameters()
{
    return parameterMap;
}

string SimulationFile::getParameter(string key, string defaultValue)
{
    // If we didn't find the key, return the default value.
    map<string,string>::iterator it = parameterMap.find(key);
    if (it == parameterMap.end()) return defaultValue;

    return it->second;
}

void SimulationFile::setParameter(string key, string value)
{
    // Set the parameter in the map.
    parameterMap[key] = value;

    // Save the parameter in the file.
    if (version <= 2)
    {
        double dvalue = atof(value.c_str());
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_double(file, "/Parameters", key.c_str(), &dvalue, 1));
    }
    else
    {
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_string(file, "/Parameters", key.c_str(), value.c_str()));
    }
}

void SimulationFile::loadModel()
{
    if (!modelLoaded)
    {
        if (version <= 3)
        {
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model", "numberSpecies", &numberSpecies));
            modelLoaded = true;
        }
        else
        {
            if (H5Lexists(file, "/Model/Reaction", H5P_DEFAULT) > 0)
            {
                HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Reaction", "numberSpecies", &numberSpecies));
                modelLoaded = true;
            }
            else
            {
                throw Exception("No model has been defined.");
            }
        }
    }
}

void SimulationFile::getReactionModel(lm::io::ReactionModel * reactionModel)
{
    // Make sure the model is not null and then clear it.
    if (reactionModel == NULL) throw InvalidArgException("reactionModel", "cannot be null");
    reactionModel->Clear();

    if (H5Lexists(file, "/Model/Reaction", H5P_DEFAULT))
    {
        // Read at least the numbers of species.
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Reaction", "numberSpecies", &numberSpecies));
        reactionModel->set_number_species(numberSpecies);
        reactionModel->set_number_reactions(0);

        // If we have the number of reactions, we must have a full model so read it.
        if (H5Aexists_by_name(file, "/Model/Reaction", "numberReactions", H5P_DEFAULT) > 0)
        {
            hsize_t dims[2];
            H5T_class_t type;
            size_t size;

            // Read the initial species counts.
            H5LTget_dataset_info(file, "/Model/Reaction/InitialSpeciesCounts", dims, &type, &size);
            if (dims[0] != numberSpecies || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/InitialSpeciesCounts");

	    // Read the initial species counts.
	    int * sbuf = new int[numberSpecies];
	    H5LTread_dataset_int(file, "/Model/Reaction/InitialSpeciesCounts", sbuf);
	    for (uint i=0; i<numberSpecies; i++) reactionModel->add_initial_species_count((uint)sbuf[i]);
	    delete [] sbuf;

            // Read the number of reactions.
            uint numberReactions;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Reaction", "numberReactions", &numberReactions));
            reactionModel->set_number_reactions(numberReactions);

            // Read the reaction tables.
            if (numberReactions > 0)
            {
                // Make sure all of the data sets are the correct size.
                H5LTget_dataset_info(file, "/Model/Reaction/ReactionTypes", dims, &type, &size);
                if (dims[0] != numberReactions || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/ReactionTypes");
                H5LTget_dataset_info(file, "/Model/Reaction/ReactionRateConstants", dims, &type, &size);
                if (dims[0] != numberReactions || dims[1] != MAX_REACTION_RATE_CONSTANTS || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/ReactionRateConstants");
                H5LTget_dataset_info(file, "/Model/Reaction/StoichiometricMatrix", dims, &type, &size);
                if (dims[0] != numberSpecies || dims[1] != numberReactions || size != sizeof(int)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/StoichiometricMatrix");
                H5LTget_dataset_info(file, "/Model/Reaction/DependencyMatrix", dims, &type, &size);
                if (dims[0] != numberSpecies || dims[1] != numberReactions || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/DependencyMatrix");

                // If we have rate noise terms, make sure they are the correct size.
                const uint NUMBER_NOISE_COLS = 2;
                bool hasNoiseTable = false;
                if (H5Lexists(file, "/Model/Reaction/ReactionRateNoise", H5P_DEFAULT))
                {
                    hasNoiseTable = true;
                    H5LTget_dataset_info(file, "/Model/Reaction/ReactionRateNoise", dims, &type, &size);
                    if (dims[0] != numberReactions || dims[1] != NUMBER_NOISE_COLS || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/ReactionRateNoise");
                }

                // Allocate some buffers for reading the data.
                int * intBuffer = new int[numberSpecies*numberReactions];
                double * doubleBuffer = new double[numberReactions*MAX_REACTION_RATE_CONSTANTS];
                double * noiseBuffer = new double[numberReactions*NUMBER_NOISE_COLS];

                // Read the reaction info.
                H5LTread_dataset_int(file, "/Model/Reaction/ReactionTypes", intBuffer);
                H5LTread_dataset_double(file, "/Model/Reaction/ReactionRateConstants", doubleBuffer);
                if (hasNoiseTable) H5LTread_dataset_double(file, "/Model/Reaction/ReactionRateNoise", noiseBuffer);
                for (uint i=0; i<numberReactions; i++)
                {
                    reactionModel->add_reaction();
                    reactionModel->mutable_reaction(i)->set_type((uint)intBuffer[i]);
                    for (uint j=0; j<MAX_REACTION_RATE_CONSTANTS; j++)
                    {
                        double k = doubleBuffer[i*MAX_REACTION_RATE_CONSTANTS+j];
                        if (!std::isnan(k))
                            reactionModel->mutable_reaction(i)->add_rate_constant(k);
                        else
                            break;
                    }

                    // If we have noise terms, set them.
                    if (hasNoiseTable)
                    {
                        double nvar = noiseBuffer[i*NUMBER_NOISE_COLS];
                        double ntau = noiseBuffer[i*NUMBER_NOISE_COLS+1];
                        if (nvar > 0.0 && ntau > 0.0 && !std::isnan(nvar) && !std::isnan(ntau))
                        {
                            reactionModel->mutable_reaction(i)->set_rate_has_noise(true);
                            reactionModel->mutable_reaction(i)->set_rate_noise_variance(nvar);
                            reactionModel->mutable_reaction(i)->set_rate_noise_tau(ntau);
                        }
                    }
                }

                // Read the matrices.
                H5LTread_dataset_int(file, "/Model/Reaction/StoichiometricMatrix", intBuffer);
                for (uint i=0; i<numberSpecies*numberReactions; i++) reactionModel->add_stoichiometric_matrix(intBuffer[i]);
                H5LTread_dataset_int(file, "/Model/Reaction/DependencyMatrix", intBuffer);
                for (uint i=0; i<numberSpecies*numberReactions; i++) reactionModel->add_dependency_matrix((uint)intBuffer[i]);

                // Free the buffers.
                delete [] noiseBuffer;
                delete [] doubleBuffer;
                delete [] intBuffer;
            }
        }
    }
}

void SimulationFile::setReactionModel(lm::io::ReactionModel * reactionModel)
{
    // Validate that the model is consistent.
    if (reactionModel == NULL) throw InvalidArgException("reactionModel", "cannot be NULL");
    if (reactionModel->number_species() == 0) throw InvalidArgException("reactionModel.number_species", "cannot be zero");
    if (reactionModel->initial_species_count_size() != (int)reactionModel->number_species()) throw InvalidArgException("reactionModel.initial_species_count", "inconsistent size");
    if (reactionModel->reaction_size() != (int)reactionModel->number_reactions()) throw InvalidArgException("reactionModel.reaction", "inconsistent size");
    if (reactionModel->stoichiometric_matrix_size() != (int)(reactionModel->number_species()*reactionModel->number_reactions())) throw InvalidArgException("reactionModel.stoichiometric_matrix", "inconsistent size");
    if (reactionModel->dependency_matrix_size() != (int)(reactionModel->number_species()*reactionModel->number_reactions())) throw InvalidArgException("reactionModel.dependency_matrix", "inconsistent size");

    // If a reaction model already exists, delete it.
    if (H5Lexists(file, "/Model/Reaction", H5P_DEFAULT))
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/Model/Reaction", H5P_DEFAULT));
    }

    // Create the group for the reaction model.
    hid_t group;
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Reaction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));

    // Write the numbers of species and reactions.
    numberSpecies = reactionModel->number_species();
    uint numberReactions = reactionModel->number_reactions();
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Reaction", "numberSpecies", &numberSpecies, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Reaction", "numberReactions", &numberReactions, 1));

    hsize_t dims[2];

    // Write the initial species counts.
    dims[0] = numberSpecies;
    HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/InitialSpeciesCounts", 1, dims, H5T_STD_U32LE, reactionModel->initial_species_count().data()));

    // If we have any reactions, write out the reaction tables.
    if (reactionModel->number_reactions())
    {
        // Write the reaction tables.
        uint * types = new uint[numberReactions];
        double * constants = new double[numberReactions*MAX_REACTION_RATE_CONSTANTS];
        bool hasNoiseTable = false;
        const uint NUMBER_NOISE_COLS = 2;
        double * noiseTerms = new double[numberReactions*NUMBER_NOISE_COLS];
        for (uint i=0; i<numberReactions*NUMBER_NOISE_COLS; i++) noiseTerms[i] = 0.0;
        for (uint i=0; i<numberReactions; i++)
        {
            types[i] = reactionModel->reaction(i).type();
            uint j=0;
            for (; j<(uint)(reactionModel->reaction(i).rate_constant_size()) && j<MAX_REACTION_RATE_CONSTANTS; j++)
                constants[i*MAX_REACTION_RATE_CONSTANTS+j] = reactionModel->reaction(i).rate_constant(j);
            for (; j<MAX_REACTION_RATE_CONSTANTS; j++)
                constants[i*MAX_REACTION_RATE_CONSTANTS+j] = NAN;

            // If we have noise terms, fill them in and mark that we need to create the noise table.
            if (reactionModel->reaction(i).rate_has_noise())
            {
                hasNoiseTable = true;
                noiseTerms[i*NUMBER_NOISE_COLS] = reactionModel->reaction(i).rate_noise_variance();
                noiseTerms[i*NUMBER_NOISE_COLS+1] = reactionModel->reaction(i).rate_noise_tau();
            }
        }
        dims[0] = numberReactions;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/ReactionTypes", 1, dims, H5T_STD_U32LE, types));
        dims[0] = numberReactions;
        dims[1] = MAX_REACTION_RATE_CONSTANTS;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/ReactionRateConstants", 2, dims, H5T_IEEE_F64LE, constants));

        // If necessary, create the noise table.
        if (hasNoiseTable)
        {
            dims[0] = numberReactions;
            dims[1] = NUMBER_NOISE_COLS;
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/ReactionRateNoise", 2, dims, H5T_IEEE_F64LE, noiseTerms));
        }
        delete [] noiseTerms;
        delete [] constants;
        delete [] types;

        // Write the matrices.
        dims[0] = numberSpecies;
        dims[1] = numberReactions;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/StoichiometricMatrix", 2, dims, H5T_STD_I32LE, reactionModel->stoichiometric_matrix().data()));
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/DependencyMatrix", 2, dims, H5T_STD_U32LE, reactionModel->dependency_matrix().data()));
    }
}

void SimulationFile::getDiffusionModel(lm::io::DiffusionModel * diffusionModel)
{
    // Make sure the model is not null and then clear it.
    if (diffusionModel == NULL) throw InvalidArgException("diffusionModel", "cannot be null");
    diffusionModel->Clear();

    if (H5Lexists(file, "/Model/Diffusion", H5P_DEFAULT))
    {
        // Read the diffusion model attributes.
        uint numberSpecies, numberReactions, numberSiteTypes, latticeXSize, latticeYSize, latticeZSize, particlesPerSite;
        double latticeSpacing;
		uint bytes_per_particle;
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "numberSpecies", &numberSpecies));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "numberReactions", &numberReactions));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "numberSiteTypes", &numberSiteTypes));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_double(file, "/Model/Diffusion", "latticeSpacing", &latticeSpacing));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "latticeXSize", &latticeXSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "latticeYSize", &latticeYSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "latticeZSize", &latticeZSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "particlesPerSite", &particlesPerSite));
		if(H5LTget_attribute_uint(file, "/Model/Diffusion", "bytes_per_particle", &bytes_per_particle) < 0)
		{
			// Must not have existed, assume byte lattice
			bytes_per_particle=1;
		}

        // Fill in the model.
        diffusionModel->set_number_species(numberSpecies);
        diffusionModel->set_number_reactions(numberReactions);
        diffusionModel->set_number_site_types(numberSiteTypes);
        diffusionModel->set_lattice_spacing(latticeSpacing);
        diffusionModel->set_lattice_x_size(latticeXSize);
        diffusionModel->set_lattice_y_size(latticeYSize);
        diffusionModel->set_lattice_z_size(latticeZSize);
        diffusionModel->set_particles_per_site(particlesPerSite);
        diffusionModel->set_bytes_per_particle(bytes_per_particle);

        hsize_t dims[3];
        H5T_class_t type;
        size_t size;

        // Read the diffusion matrix.
        H5LTget_dataset_info(file, "/Model/Diffusion/DiffusionMatrix", dims, &type, &size);
        if (dims[0] != numberSiteTypes || dims[1] != numberSiteTypes || dims[2] != numberSpecies || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/DiffusionMatrix");
        double * D = new double[numberSiteTypes*numberSiteTypes*numberSpecies];
        H5LTread_dataset_double(file, "/Model/Diffusion/DiffusionMatrix", D);
        for (uint i=0; i<numberSiteTypes*numberSiteTypes*numberSpecies; i++) diffusionModel->add_diffusion_matrix(D[i]);
        delete [] D;

        // If we have reactions, read the reaction location matrix.
        if (numberReactions > 0)
        {
            H5LTget_dataset_info(file, "/Model/Diffusion/ReactionLocationMatrix", dims, &type, &size);
            if (dims[0] != numberReactions || dims[1] != numberSiteTypes || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/ReactionLocationMatrix");
            uint * RL = new uint[numberReactions*numberSiteTypes];
			H5LTread_dataset(file, "/Model/Diffusion/ReactionLocationMatrix", H5T_STD_U32LE, RL);
			for (uint i=0; i<numberReactions*numberSiteTypes; i++) diffusionModel->add_reaction_location_matrix(RL[i]);
	        delete [] RL;
        }
    }
}

void SimulationFile::getDiffusionModelLattice(lm::io::DiffusionModel * m, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize)
{
    int ndims;
    hsize_t dims[4];
    H5T_class_t type;
    size_t size;

    // Read the lattice.
    if (lattice != NULL)
    {
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Model/Diffusion/Lattice", &ndims));
		if (ndims != 4) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Model/Diffusion/Lattice",dims, &type, &size));
		switch(m->bytes_per_particle())
		{
			case 1:
				if (m->lattice_x_size() != dims[0] || m->lattice_y_size() != dims[1] || m->lattice_z_size() != dims[2] || m->particles_per_site() != dims[3] || size != sizeof(uint8_t)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
				if (dims[0]*dims[1]*dims[2]*dims[3] != latticeSize) throw InvalidArgException("lattice", "incorrect lattice size");
				HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/Lattice", H5T_NATIVE_UINT8, lattice));
				break;
			case 4:
				if (m->lattice_x_size() != dims[0] || m->lattice_y_size() != dims[1] || m->lattice_z_size() != dims[2] || m->particles_per_site() != dims[3] || size != sizeof(uint32_t)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
				if (dims[0]*dims[1]*dims[2]*dims[3]*size != latticeSize) throw InvalidArgException("lattice", "incorrect lattice size");
				HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/Lattice", H5T_NATIVE_UINT32, lattice));
				break;
			default:
				throw Exception("Unknown bytes per particle", m->bytes_per_particle());
		}
    }

    // Read the lattice sites.
    if (latticeSites != NULL)
    {
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Model/Diffusion/LatticeSites", &ndims));
		if (ndims != 3) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/LatticeSites");
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Model/Diffusion/LatticeSites",dims, &type, &size));
		if (m->lattice_x_size() != dims[0] || m->lattice_y_size() != dims[1] || m->lattice_z_size() != dims[2] || size != sizeof(uint8_t)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/LatticeSites");
		if (dims[0]*dims[1]*dims[2] != latticeSitesSize) throw InvalidArgException("latticeSites", "insufficient size");
		HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/LatticeSites", H5T_NATIVE_UINT8, latticeSites));
    }
}

void SimulationFile::getDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::Lattice * lattice)
{
	if (lattice->getXSize() != m->lattice_x_size() || lattice->getYSize() != m->lattice_y_size() || lattice->getZSize() != m->lattice_z_size() || lattice->getMaxOccupancy() != m->particles_per_site()) throw InvalidArgException("lattice", "lattice size not consistent with the diffusion model");

	int ndims;
    hsize_t dims[4];
    H5T_class_t type;
    size_t size;

    if (lattice != NULL)
    {
        // Read the lattice.
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Model/Diffusion/Lattice", &ndims));
		if (ndims != 4) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Model/Diffusion/Lattice",dims, &type, &size));
		if (m->lattice_x_size() != dims[0] || m->lattice_y_size() != dims[1] || m->lattice_z_size() != dims[2] || m->particles_per_site() != dims[3] || size != sizeof(uint8_t)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
	    size_t latticeDataSize = dims[0]*dims[1]*dims[2]*dims[3];
	    uint8_t * latticeData = new uint8_t[latticeDataSize];
		HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/Lattice", H5T_NATIVE_UINT8, latticeData));
		lattice->setFromRowMajorByteData(latticeData, latticeDataSize);
		delete [] latticeData;

		// Read the lattice sites.
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Model/Diffusion/LatticeSites", &ndims));
		if (ndims != 3) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/LatticeSites");
		HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Model/Diffusion/LatticeSites",dims, &type, &size));
		if (m->lattice_x_size() != dims[0] || m->lattice_y_size() != dims[1] || m->lattice_z_size() != dims[2] || size != sizeof(uint8_t)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/LatticeSites");
	    size_t latticeSitesSize = dims[0]*dims[1]*dims[2];
	    uint8_t * latticeSites = new uint8_t[latticeSitesSize];
		HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/LatticeSites", H5T_NATIVE_UINT8, latticeSites));
		lattice->setSitesFromRowMajorByteData(latticeSites, latticeSitesSize);
		delete [] latticeSites;
    }
}

void SimulationFile::setDiffusionModel(lm::io::DiffusionModel * diffusionModel)
{
    // Validate that the model is consistent.
    if (diffusionModel == NULL) throw InvalidArgException("diffusionModel", "cannot be NULL");
    if (diffusionModel->number_species() == 0) throw InvalidArgException("diffusionModel.number_species", "cannot be zero");
    if (diffusionModel->number_site_types() == 0) throw InvalidArgException("diffusionModel.number_site_types", "cannot be zero");
    if (diffusionModel->diffusion_matrix_size() != (int)(diffusionModel->number_site_types()*diffusionModel->number_site_types()*diffusionModel->number_species())) throw InvalidArgException("diffusion.diffusion_matrix", "inconsistent size");
    if (diffusionModel->reaction_location_matrix_size() != (int)(diffusionModel->number_reactions()*diffusionModel->number_site_types())) throw InvalidArgException("diffusion.reaction_location_matrix", "inconsistent size");
    if (diffusionModel->lattice_spacing() <= 0.0) throw InvalidArgException("diffusionModel.lattice_spacing", "must be greater than zero");
    if (diffusionModel->lattice_x_size() == 0) throw InvalidArgException("diffusionModel.lattice_x_size", "cannot be zero");
    if (diffusionModel->lattice_y_size() == 0) throw InvalidArgException("diffusionModel.lattice_y_size", "cannot be zero");
    if (diffusionModel->lattice_z_size() == 0) throw InvalidArgException("diffusionModel.lattice_z_size", "cannot be zero");
    if (diffusionModel->particles_per_site() == 0) throw InvalidArgException("diffusionModel.particles_per_site", "cannot be zero");
    if (diffusionModel->bytes_per_particle() == 0) throw InvalidArgException("diffusionModel.bytes_per_particle", "cannot be zero");

    // If a diffusion model already exists, delete it.
    if (H5Lexists(file, "/Model/Diffusion", H5P_DEFAULT))
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/Model/Diffusion", H5P_DEFAULT));
    }

    // Create the group for the reaction model.
    hid_t group;
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Diffusion", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));

    // Write the attributes.
    numberSpecies = diffusionModel->number_species();
    uint numberReactions = diffusionModel->number_reactions();
    uint numberSiteTypes = diffusionModel->number_site_types();
    double latticeSpacing = diffusionModel->lattice_spacing();
    uint latticeXSize = diffusionModel->lattice_x_size();
    uint latticeYSize = diffusionModel->lattice_y_size();
    uint latticeZSize = diffusionModel->lattice_z_size();
    uint particlesPerSite = diffusionModel->particles_per_site();
    uint bytes_per_particle = diffusionModel->bytes_per_particle();
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "numberSpecies", &numberSpecies, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "numberReactions", &numberReactions, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "numberSiteTypes", &numberSiteTypes, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_double(file, "/Model/Diffusion", "latticeSpacing", &latticeSpacing, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "latticeXSize", &latticeXSize, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "latticeYSize", &latticeYSize, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "latticeZSize", &latticeZSize, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "particlesPerSite", &particlesPerSite, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "bytes_per_particle", &bytes_per_particle, 1));

    // Write the diffusion matrix.
    {
		const unsigned int RANK=3;
		hsize_t dims[RANK];
		dims[0] = numberSiteTypes;
		dims[1] = numberSiteTypes;
		dims[2] = numberSpecies;
		HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Diffusion/DiffusionMatrix", RANK, dims, H5T_IEEE_F64LE, diffusionModel->diffusion_matrix().data()));
	}

    // Write the reaction location matrix.
    if (numberReactions > 0)
    {
		const unsigned int RANK=2;
		hsize_t dims[RANK];
		dims[0] = numberReactions;
		dims[1] = numberSiteTypes;
		HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Diffusion/ReactionLocationMatrix", RANK, dims, H5T_STD_U32LE, diffusionModel->reaction_location_matrix().data()));
    }
}

void SimulationFile::writeParticleLattice_U8LE(const char* path,  uint8_t *data, unsigned int x, unsigned int y, unsigned int z, unsigned int w)
{
	unsigned int RANK=4;
	hsize_t dims[RANK], chunk[RANK];
	dims[0] = x;
	dims[1] = y;
	dims[2] = z;
	dims[3] = w;
	chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,x);
	chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,y);
	chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,z);
	chunk[3] = w;
	hid_t dataspaceHandle, dcplHandle, datasetHandle;
	HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
	HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
	HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
	HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
	HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(file, path, H5T_STD_U8LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
	HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
	HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
	HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
	HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
}

void SimulationFile::writeParticleLattice_U32LE(hid_t handle, const char* path,  unsigned int *data, unsigned int x, unsigned int y, unsigned int z, unsigned int w)
{
	unsigned int RANK=4;
	hsize_t dims[RANK], chunk[RANK];
	dims[0] = x;
	dims[1] = y;
	dims[2] = z;
	dims[3] = w;
	chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,x);
	chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,y);
	chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,z);
	chunk[3] = w;
	hid_t dataspaceHandle, dcplHandle, datasetHandle;
	HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
	HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
	HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
	HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
	HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(handle, path, H5T_STD_U32LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
	HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
	HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
	HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
	HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
}

void SimulationFile::write3DLattice_U8LE(const char* path,  uint8_t *data, unsigned int x, unsigned int y, unsigned int z)
{
	unsigned int RANK=3;
	hsize_t dims[RANK], chunk[RANK];
	dims[0] = x;
	dims[1] = y;
	dims[2] = z;
	chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,x);
	chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,y);
	chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,z);
	hid_t dataspaceHandle, dcplHandle, datasetHandle;
	HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
	HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
	HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
	HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
	HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(file, path, H5T_STD_U8LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
	HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
	HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
	HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
	HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
}

void SimulationFile::setDiffusionModelLattice(lm::io::DiffusionModel * m, uint8_t * lattice, uint8_t * latticeSites)
{
	writeParticleLattice_U8LE("/Model/Diffusion/Lattice", lattice, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size(), m->particles_per_site());
	write3DLattice_U8LE("/Model/Diffusion/LatticeSites", latticeSites, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size());
}

void SimulationFile::setDiffusionModelLattice(lm::io::DiffusionModel * m, uint32_t * lattice, uint8_t * latticeSites)
{
	writeParticleLattice_U32LE(file, "/Model/Diffusion/Lattice", lattice, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size(), m->particles_per_site());
	write3DLattice_U8LE("/Model/Diffusion/LatticeSites", latticeSites, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size());
}

void SimulationFile::setDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::ByteLattice * lattice)
{
	size_t latticeDataSize = m->lattice_x_size()*m->lattice_y_size()*m->lattice_z_size()*m->particles_per_site();
	byte * latticeData =  new byte[latticeDataSize];
	lm::rdme::Lattice::rowMajorByteSerialize(latticeData, lattice, latticeDataSize);
	writeParticleLattice_U8LE("/Model/Diffusion/Lattice", latticeData, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size(), m->particles_per_site());
	delete [] latticeData;

	size_t latticeSitesDataSize = m->lattice_x_size()*m->lattice_y_size()*m->lattice_z_size();
	byte * latticeSitesData =  new byte[latticeSitesDataSize];
	lm::rdme::Lattice::rowMajorByteSerializeSites(latticeSitesData, lattice, latticeSitesDataSize);
	write3DLattice_U8LE("/Model/Diffusion/LatticeSites", latticeSitesData, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size());
	delete [] latticeSitesData;
}

void SimulationFile::setDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::IntLattice * lattice)
{
	size_t latticeDataSize = m->lattice_x_size()*m->lattice_y_size()*m->lattice_z_size()*m->particles_per_site();
	uint32_t * latticeData =  new uint32_t[latticeDataSize];
	lm::rdme::Lattice::rowMajorIntSerialize(latticeData, lattice, latticeDataSize);
	writeParticleLattice_U32LE(file, "/Model/Diffusion/Lattice", latticeData, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size(), m->particles_per_site());
	delete [] latticeData;
			
	size_t latticeSitesDataSize = m->lattice_x_size()*m->lattice_y_size()*m->lattice_z_size();
	byte * latticeSitesData =  new byte[latticeSitesDataSize];
	lm::rdme::Lattice::rowMajorByteSerializeSites(latticeSitesData, lattice, latticeSitesDataSize);
	write3DLattice_U8LE("/Model/Diffusion/LatticeSites", latticeSitesData, m->lattice_x_size(), m->lattice_y_size(), m->lattice_z_size());
	delete [] latticeSitesData;
}

void SimulationFile::getSpatialModel(lm::io::SpatialModel * spatialModel)
{
    // Make sure the model is not null and then clear it.
    if (spatialModel == NULL) throw InvalidArgException("spatialModel", "cannot be null");
    spatialModel->Clear();

    if (H5Lexists(file, "/Model/Spatial", H5P_DEFAULT))
    {
        // Read the diffusion model attributes.
        uint numberRegions, numberObstacles;
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Spatial/Regions", "numberRegions", &numberRegions));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Spatial/Obstacles", "numberObstacles", &numberObstacles));

        hsize_t dims[2];
        H5T_class_t type;
        size_t size;

        if (numberRegions > 0)
        {
            // Read the region types table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Regions/Types", dims, &type, &size);
                if (dims[0] != numberRegions || dims[1] != 2 || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Regions/Types");
                uint * types = new uint[dims[0]*dims[1]];
                H5LTread_dataset(file, "/Model/Spatial/Regions/Types", H5T_STD_U32LE, types);
                for (uint i=0; i<dims[0]; i++)
                {
                    spatialModel->add_region();
                    spatialModel->mutable_region(i)->set_shape(types[i*dims[1]]);
                    spatialModel->mutable_region(i)->set_site_type(types[i*dims[1]+1]);
                }
                delete [] types;
            }

            // Read the shape parameters table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Regions/ShapeParameters", dims, &type, &size);
                if (dims[0] != numberRegions || dims[1] != MAX_SHAPE_PARAMETERS || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Regions/ShapeParameters");
                double * shapeParams = new double[dims[0]*dims[1]];
                H5LTread_dataset_double(file, "/Model/Spatial/Regions/ShapeParameters", shapeParams);
                for (uint i=0; i<dims[0]; i++)
                {
                    for (uint j=0; j<dims[1]; j++)
                    {
                        double shapeParam = shapeParams[i*dims[1]+j];
                        if (!std::isnan(shapeParam))
                            spatialModel->mutable_region(i)->add_shape_parameter(shapeParam);
                        else
                            break;
                    }
                }
                delete [] shapeParams;
            }
        }

        if (numberObstacles > 0)
        {
            // Read the obstacle types table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Obstacles/Types", dims, &type, &size);
                if (dims[0] != numberObstacles || dims[1] != 2 || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Obstacles/Types");
                uint * types = new uint[dims[0]*dims[1]];
                H5LTread_dataset(file, "/Model/Spatial/Obstacles/Types", H5T_STD_U32LE, types);
                for (uint i=0; i<dims[0]; i++)
                {
                    spatialModel->add_obstacle();
                    spatialModel->mutable_obstacle(i)->set_shape(types[i*dims[1]]);
                    spatialModel->mutable_obstacle(i)->set_site_type(types[i*dims[1]+1]);
                }
                delete [] types;
            }

            // Read the obstacle parameters table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Obstacles/ShapeParameters", dims, &type, &size);
                if (dims[0] != numberObstacles || dims[1] != MAX_SHAPE_PARAMETERS || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Obstacles/ShapeParameters");
                double * shapeParams = new double[dims[0]*dims[1]];
                H5LTread_dataset_double(file, "/Model/Spatial/Obstacles/ShapeParameters", shapeParams);
                for (uint i=0; i<dims[0]; i++)
                {
                    for (uint j=0; j<dims[1]; j++)
                    {
                        double shapeParam = shapeParams[i*dims[1]+j];
                        if (!std::isnan(shapeParam))
                            spatialModel->mutable_obstacle(i)->add_shape_parameter(shapeParam);
                        else
                            break;
                    }
                }
                delete [] shapeParams;
            }
        }
    }
}

void SimulationFile::setSpatialModel(lm::io::SpatialModel * spatialModel)
{
    // Validate that the model is consistent.
    if (spatialModel == NULL) throw InvalidArgException("spatialModel", "cannot be NULL");

    // If a diffusion model already exists, delete it.
    if (H5Lexists(file, "/Model/Spatial", H5P_DEFAULT))
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/Model/Spatial", H5P_DEFAULT));
    }

    // Create the group for the spatial model.
    hid_t group;
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Spatial", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Spatial/Regions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Spatial/Obstacles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));

    // Write the attributes.
    uint numberRegions = spatialModel->region_size();
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Spatial/Regions", "numberRegions", &numberRegions, 1));
    uint numberObstacles = spatialModel->obstacle_size();
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Spatial/Obstacles", "numberObstacles", &numberObstacles, 1));

    if (numberRegions > 0)
    {
        // Write the region types table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberRegions;
            dims[1] = 2;
            uint * types = new uint[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                types[i*dims[1]] = spatialModel->region(i).shape();
                types[i*dims[1]+1] = spatialModel->region(i).site_type();
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Regions/Types", RANK, dims, H5T_STD_U32LE, types));
            delete [] types;
        }

        // Write the shape parameters table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberRegions;
            dims[1] = MAX_SHAPE_PARAMETERS;
            double * shapeParams = new double[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                for (uint j=0; j<dims[1]; j++)
                {
                    if (j < (uint)spatialModel->region(i).shape_parameter_size())
                        shapeParams[i*dims[1]+j] = spatialModel->region(i).shape_parameter(j);
                    else
                        shapeParams[i*dims[1]+j] = NAN;
                }
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Regions/ShapeParameters", RANK, dims, H5T_IEEE_F64LE, shapeParams));
            delete [] shapeParams;
        }
    }

    if (numberObstacles > 0)
    {
        // Write the obstacle types table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberObstacles;
            dims[1] = 2;
            uint * types = new uint[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                types[i*dims[1]] = spatialModel->obstacle(i).shape();
                types[i*dims[1]+1] = spatialModel->obstacle(i).site_type();
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Obstacles/Types", RANK, dims, H5T_STD_U32LE, types));
            delete [] types;
        }

        // Write the shape parameters table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberObstacles;
            dims[1] = MAX_SHAPE_PARAMETERS;
            double * shapeParams = new double[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                for (uint j=0; j<dims[1]; j++)
                {
                    if (j < (uint)spatialModel->obstacle(i).shape_parameter_size())
                        shapeParams[i*dims[1]+j] = spatialModel->obstacle(i).shape_parameter(j);
                    else
                        shapeParams[i*dims[1]+j] = NAN;
                }
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Obstacles/ShapeParameters", RANK, dims, H5T_IEEE_F64LE, shapeParams));
            delete [] shapeParams;
        }
    }
}

bool SimulationFile::replicateExists(unsigned int replicate)
{
    char replicateName[8];
    snprintf(replicateName, sizeof(replicateName), "%07d", replicate);
    if (H5Lexists(simulationsGroup, replicateName, H5P_DEFAULT) > 0) return true;
    return false;
}

void SimulationFile::openReplicate(unsigned int replicate)
{
    openReplicateHandles(replicate);
}

void SimulationFile::appendSpeciesCounts(unsigned int replicate, lm::io::SpeciesCounts * speciesCounts)
{
    ReplicateHandles * handles = openReplicateHandles(replicate);

    // Update the species counts dataset.
    {
        // Get the current size of the dataset.
        unsigned int RANK=2;
        hsize_t dims[RANK];
        hid_t dataspace_id;
        int result;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountsDataset));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

        // Extend the dataset by the number of rows in the data set.
        dims[0] += speciesCounts->number_entries();
        HDF5_EXCEPTION_CHECK(H5Dset_extent(handles->speciesCountsDataset, dims));

        // Create the memory dataset.
        hid_t memspace_id;
        hsize_t memDims[RANK];
        memDims[0] = speciesCounts->number_entries();
        memDims[1] = speciesCounts->number_species();
        HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountsDataset));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-speciesCounts->number_entries();
        start[1] = 0;
        count[0] = memDims[0];
        count[1] = memDims[1];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(handles->speciesCountsDataset, H5T_NATIVE_INT32, memspace_id, dataspace_id, H5P_DEFAULT, speciesCounts->species_count().data()));

        // Cleanup some resources.
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    }

    // Update the species count times dataset.
    {
        // Get the current size of the dataset.
        unsigned int RANK=1;
        hsize_t dims[RANK];
        hid_t dataspace_id;
        int result;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountTimesDataset));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

        // Extend the dataset by the number of rows in the data set.
        dims[0] += speciesCounts->number_entries();
        HDF5_EXCEPTION_CHECK(H5Dset_extent(handles->speciesCountTimesDataset, dims));

        // Create the memory dataset.
        hid_t memspace_id;
        hsize_t memDims[RANK];
        memDims[0] = speciesCounts->number_entries();
        HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountTimesDataset));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-speciesCounts->number_entries();
        count[0] = memDims[0];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(handles->speciesCountTimesDataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, speciesCounts->time().data()));

        // Cleanup some resources.
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    }
}

void SimulationFile::appendLattice(unsigned int replicate, lm::io::Lattice * lattice, uint8_t * latticeData, size_t latticeDataSize)
{
    if (size_t(lattice->lattice_x_size())*size_t(lattice->lattice_y_size())*size_t(lattice->lattice_z_size())*size_t(lattice->particles_per_site())*sizeof(uint8_t) != latticeDataSize) throw InvalidArgException("lattice", "incorrect lattice size");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the lattice group.
    hid_t latticeGroupHandle, latticeTimesDatasetHandle;
    if (H5Lexists(replicateHandles->group, "Lattice", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CALL(latticeGroupHandle,H5Gopen2(replicateHandles->group, "Lattice", H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dopen2(replicateHandles->group, "LatticeTimes", H5P_DEFAULT));
    }
    else
    {
        HDF5_EXCEPTION_CALL(latticeGroupHandle,H5Gcreate2(replicateHandles->group, "Lattice", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Create the lattice times dataset.
        const unsigned int RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspaceHandle, propsHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(propsHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(propsHandle, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dcreate2(replicateHandles->group, "LatticeTimes", H5T_IEEE_F64LE, dataspaceHandle, H5P_DEFAULT, propsHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(propsHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    }

    // Append the time to the times data set.
    uint latticeIndex;
    {
        // Get the current size of the dataset.
        unsigned int RANK=1;
        hsize_t dims[RANK];
        hid_t dataspaceHandle;
        int result;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Dget_space(latticeTimesDatasetHandle));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspaceHandle, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));

        // Extend the dataset by the number of rows in the data set.
        dims[0] += 1;
        HDF5_EXCEPTION_CHECK(H5Dset_extent(latticeTimesDatasetHandle, dims));

        // Create the memory dataset.
        hid_t memspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = 1;
        HDF5_EXCEPTION_CALL(memspaceHandle,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Dget_space(latticeTimesDatasetHandle));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-1;
        latticeIndex=start[0];
        count[0] = memDims[0];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, start, NULL, count, NULL));
        double time = lattice->time();
        HDF5_EXCEPTION_CHECK(H5Dwrite(latticeTimesDatasetHandle, H5T_NATIVE_DOUBLE, memspaceHandle, dataspaceHandle, H5P_DEFAULT, &time));

        // Cleanup some resources.
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspaceHandle));
    }

    // Create the lattice data set.
    {
        char latticeDatasetName[11];
        snprintf(latticeDatasetName, sizeof(latticeDatasetName), "%010d", latticeIndex);
        const unsigned int RANK=4;
        hsize_t dims[RANK], chunk[RANK];
        dims[0] = lattice->lattice_x_size();
        dims[1] = lattice->lattice_y_size();
        dims[2] = lattice->lattice_z_size();
        dims[3] = lattice->particles_per_site();
        chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,lattice->lattice_x_size());
        chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,lattice->lattice_y_size());
        chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,lattice->lattice_z_size());
        chunk[3] = lattice->particles_per_site();
        hid_t dataspaceHandle, dcplHandle, datasetHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
        HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
        HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(latticeGroupHandle, latticeDatasetName, H5T_STD_U8LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, latticeData));
        HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    }

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Dclose(latticeTimesDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(latticeGroupHandle));

}

void SimulationFile::appendLattice_U32LE(unsigned int replicate, lm::io::Lattice * lattice, uint32_t * latticeData, size_t latticeDataSize)
{
    if (lattice->lattice_x_size()*lattice->lattice_y_size()*lattice->lattice_z_size()*lattice->particles_per_site()*sizeof(uint32_t) != latticeDataSize) throw InvalidArgException("lattice", "incorrect lattice size");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the lattice group.
    hid_t latticeGroupHandle, latticeTimesDatasetHandle;
    if (H5Lexists(replicateHandles->group, "Lattice", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CALL(latticeGroupHandle,H5Gopen2(replicateHandles->group, "Lattice", H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dopen2(replicateHandles->group, "LatticeTimes", H5P_DEFAULT));
    }
    else
    {
        HDF5_EXCEPTION_CALL(latticeGroupHandle,H5Gcreate2(replicateHandles->group, "Lattice", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Create the lattice times dataset.
        const unsigned int RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspaceHandle, propsHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(propsHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(propsHandle, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dcreate2(replicateHandles->group, "LatticeTimes", H5T_IEEE_F64LE, dataspaceHandle, H5P_DEFAULT, propsHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(propsHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    }

    // Append the time to the times data set.
    uint latticeIndex;
    {
        // Get the current size of the dataset.
        unsigned int RANK=1;
        hsize_t dims[RANK];
        hid_t dataspaceHandle;
        int result;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Dget_space(latticeTimesDatasetHandle));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspaceHandle, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));

        // Extend the dataset by the number of rows in the data set.
        dims[0] += 1;
        HDF5_EXCEPTION_CHECK(H5Dset_extent(latticeTimesDatasetHandle, dims));

        // Create the memory dataset.
        hid_t memspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = 1;
        HDF5_EXCEPTION_CALL(memspaceHandle,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Dget_space(latticeTimesDatasetHandle));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-1;
        latticeIndex=start[0];
        count[0] = memDims[0];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, start, NULL, count, NULL));
        double time = lattice->time();
        HDF5_EXCEPTION_CHECK(H5Dwrite(latticeTimesDatasetHandle, H5T_NATIVE_DOUBLE, memspaceHandle, dataspaceHandle, H5P_DEFAULT, &time));

        // Cleanup some resources.
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspaceHandle));
    }

	char latticeDatasetName[11];
	snprintf(latticeDatasetName, sizeof(latticeDatasetName), "%010d", latticeIndex);
	writeParticleLattice_U32LE(latticeGroupHandle, latticeDatasetName, latticeData, lattice->lattice_x_size(),lattice->lattice_y_size(),lattice->lattice_z_size(),lattice->particles_per_site());
	

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Dclose(latticeTimesDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(latticeGroupHandle));

}

void SimulationFile::appendSites(unsigned int replicate, lm::io::Lattice * lattice, uint8_t * siteData, size_t siteDataSize)
{
    if (lattice->lattice_x_size()*lattice->lattice_y_size()*lattice->lattice_z_size()*sizeof(uint8_t) != siteDataSize) throw InvalidArgException("siteData", "incorrect lattice size");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the sites group
    hid_t siteGroupHandle, latticeTimesDatasetHandle;
    if (H5Lexists(replicateHandles->group, "Sites", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CALL(siteGroupHandle,H5Gopen2(replicateHandles->group, "Sites", H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dopen2(replicateHandles->group, "SiteTimes", H5P_DEFAULT));
    }
    else
    {
        HDF5_EXCEPTION_CALL(siteGroupHandle,H5Gcreate2(replicateHandles->group, "Sites", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Create the SiteTimes dataset.
        const unsigned int RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspaceHandle, propsHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(propsHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(propsHandle, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dcreate2(replicateHandles->group, "SiteTimes", H5T_IEEE_F64LE, dataspaceHandle, H5P_DEFAULT, propsHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(propsHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    }

    // Append the time to the times data set.
    uint latticeIndex;
    {
        // Get the current size of the dataset.
        unsigned int RANK=1;
        hsize_t dims[RANK];
        hid_t dataspaceHandle;
        int result;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Dget_space(latticeTimesDatasetHandle));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspaceHandle, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));

        // Extend the dataset by the number of rows in the data set.
        dims[0] += 1;
        HDF5_EXCEPTION_CHECK(H5Dset_extent(latticeTimesDatasetHandle, dims));

        // Create the memory dataset.
        hid_t memspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = 1;
        HDF5_EXCEPTION_CALL(memspaceHandle,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Dget_space(latticeTimesDatasetHandle));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-1;
        latticeIndex=start[0];
        count[0] = memDims[0];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, start, NULL, count, NULL));
        double time = lattice->time();
        HDF5_EXCEPTION_CHECK(H5Dwrite(latticeTimesDatasetHandle, H5T_NATIVE_DOUBLE, memspaceHandle, dataspaceHandle, H5P_DEFAULT, &time));

        // Cleanup some resources.
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspaceHandle));
    }

    // Create the lattice data set.
    {
        char latticeDatasetName[11];
        snprintf(latticeDatasetName, sizeof(latticeDatasetName), "%010d", latticeIndex);
        const unsigned int RANK=3;
        hsize_t dims[RANK], chunk[RANK];
        dims[0] = lattice->lattice_x_size();
        dims[1] = lattice->lattice_y_size();
        dims[2] = lattice->lattice_z_size();
        chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,lattice->lattice_x_size());
        chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,lattice->lattice_y_size());
        chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,lattice->lattice_z_size());
        hid_t dataspaceHandle, dcplHandle, datasetHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
        HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
        HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(siteGroupHandle, latticeDatasetName, H5T_STD_U8LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, siteData));
        HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    }

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Dclose(latticeTimesDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(siteGroupHandle));

}

void SimulationFile::appendParameterValues(unsigned int replicate, lm::io::ParameterValues * parameterValues)
{
    if (parameterValues->value_size() != parameterValues->time_size()) throw InvalidArgException("parameterValues", "inconsistent number of entries");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the parameter values group.
    hid_t pvGroupHandle;
    if ((pvGroupHandle=H5Gopen2(replicateHandles->group, "ParameterValues", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(pvGroupHandle,H5Gcreate2(replicateHandles->group, "ParameterValues", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    // Open or create the parameter values dataset.
    hid_t pvDatasetHandle;
    if ((pvDatasetHandle=H5Dopen2(pvGroupHandle, parameterValues->parameter().c_str(), H5P_DEFAULT)) < 0)
    {
        // Create the datasets.
        uint RANK=2;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        dims[1] = 2;
        maxDims[0] = H5S_UNLIMITED;
        maxDims[1] = 2;
        chunkDims[0] = 100;
        chunkDims[1] = 2;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(pvDatasetHandle,H5Dcreate2(pvGroupHandle, parameterValues->parameter().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Get the current size of the dataset.
    unsigned int RANK=2;
    hsize_t dims[RANK];
    hid_t dataspace_id;
    int result;
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(pvDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

    // Extend the dataset by the number of rows in the data set.
    dims[0] += parameterValues->value_size();
    HDF5_EXCEPTION_CHECK(H5Dset_extent(pvDatasetHandle, dims));

    // Create the memory dataset.
    hid_t memspace_id;
    hsize_t memDims[RANK];
    memDims[0] = parameterValues->value_size();
    memDims[1] = 1;
    HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

    // Write the new times.
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(pvDatasetHandle));
    hsize_t start[RANK], count[RANK];
    start[0] = dims[0]-parameterValues->value_size();
    count[0] = parameterValues->value_size();
    start[1] = 0;
    count[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(pvDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, parameterValues->time().data()));
    start[1] = 1;
    count[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(pvDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, parameterValues->value().data()));

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    HDF5_EXCEPTION_CHECK(H5Dclose(pvDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(pvGroupHandle));
}


void SimulationFile::setFirstPassageTimes(unsigned int replicate, lm::io::FirstPassageTimes * firstPassageTimes)
{
    // Make sure the data is consistent.
    if (firstPassageTimes->species_count_size() == 0) throw InvalidArgException("firstPassageTimes", "no entries to save");
    if (firstPassageTimes->species_count_size() != firstPassageTimes->first_passage_time_size()) throw InvalidArgException("firstPassageTimes", "inconsistent number of first passage time entries");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open the first passage time group.
    hid_t fptGroupHandle;
    if ((fptGroupHandle=H5Gopen2(replicateHandles->group, "FirstPassageTimes", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(fptGroupHandle,H5Gcreate2(replicateHandles->group, "FirstPassageTimes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    // Construct a string representation of the species name.
    std::stringstream ss;
    ss.fill('0');
    ss.width(2);
    ss << firstPassageTimes->species();
    string speciesString = ss.str();

    // Open the species group.
    hid_t speciesGroupHandle, countsDatasetHandle, timesDatasetHandle;
    if ((speciesGroupHandle=H5Gopen2(fptGroupHandle, speciesString.c_str(), H5P_DEFAULT)) >= 0)
    {
        // Open the data sets.
        HDF5_EXCEPTION_CALL(countsDatasetHandle,H5Dopen2(speciesGroupHandle, "Counts", H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(timesDatasetHandle,H5Dopen2(speciesGroupHandle, "Times", H5P_DEFAULT));
    }
    else
    {
        // Create the group.
        HDF5_EXCEPTION_CALL(speciesGroupHandle,H5Gcreate2(fptGroupHandle, speciesString.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Create the datasets.
        uint RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(countsDatasetHandle,H5Dcreate2(speciesGroupHandle, "Counts", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(timesDatasetHandle,H5Dcreate2(speciesGroupHandle, "Times", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Read the lowest and highest first passage times in the file.
    uint RANK=1;
    int result;
    hid_t countsDataspace;
    hsize_t countsDims[RANK];
    HDF5_EXCEPTION_CALL(countsDataspace,H5Dget_space(countsDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(countsDataspace, countsDims, NULL));
    uint minCount=0, maxCount=0;
    if (countsDims[0] == 1)
    {
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &minCount));
        maxCount = minCount;
    }
    else if (countsDims[0] > 1)
    {
        // Create the memory dataspace.
        hid_t memDataspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = 1;
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));

        hsize_t start[RANK], count[RANK];
        count[0] = 1;
        start[0] = 0;
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(countsDataspace, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_UINT32, memDataspaceHandle, countsDataspace, H5P_DEFAULT, &minCount));
        start[0] = countsDims[0]-1;
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(countsDataspace, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_UINT32, memDataspaceHandle, countsDataspace, H5P_DEFAULT, &maxCount));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));
    }
    HDF5_EXCEPTION_CHECK(H5Sclose(countsDataspace));

    // Figure out the min and the max from the new counts.
    uint32_t minNewCount=firstPassageTimes->species_count(0), maxNewCount=firstPassageTimes->species_count(0);
    for (int i=1; i<firstPassageTimes->species_count_size(); i++)
    {
        if (firstPassageTimes->species_count(i) < minNewCount)
            minNewCount = firstPassageTimes->species_count(i);
        else if (firstPassageTimes->species_count(i) > maxNewCount)
            maxNewCount = firstPassageTimes->species_count(i);
    }

    // If there are no existing records, just insert one large block.
    if (countsDims[0] == 0)
    {
        // Allocate the block.
        countsDims[0] = maxNewCount-minNewCount+1;
        uint32_t * counts = new uint32_t[countsDims[0]];
        double * times = new double[countsDims[0]];
        for (uint i=0; i<countsDims[0]; i++)
        {
            counts[i] = 0;
            times[i] = -1.0;
        }

        // Fill in the block.
        for (int i=0; i<firstPassageTimes->species_count_size(); i++)
        {
            uint32_t count = firstPassageTimes->species_count(i);
            counts[count-minNewCount] = count;
            times[count-minNewCount] = firstPassageTimes->first_passage_time(i);
        }

        // Extend the datasets and write the block.
        HDF5_EXCEPTION_CHECK(H5Dset_extent(countsDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(countsDatasetHandle, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, counts));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(timesDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(timesDatasetHandle, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, times));

        // Free the buffers.
        delete [] counts;
        delete [] times;
    }

    // Otherwise if there are only new records to append to the end, just add them
    else if (maxNewCount > maxCount && minNewCount > maxCount)
    {
        // Allocate the new buffer.
        hsize_t endingStart[RANK], endingCount[RANK];
        endingCount[0] = maxNewCount-maxCount;
        uint32_t * newEndingCounts = new uint32_t[endingCount[0]];
        double * newEndingTimes = new double[endingCount[0]];
        for (uint i=0; i<endingCount[0]; i++)
        {
            newEndingCounts[i] = 0;
            newEndingTimes[i] = -1.0;
        }

        // Fill in the buffer.
        for (int i=0; i<firstPassageTimes->species_count_size(); i++)
        {
            uint32_t count = firstPassageTimes->species_count(i);
            newEndingCounts[count-maxCount-1] = count;
            newEndingTimes[count-maxCount-1] = firstPassageTimes->first_passage_time(i);
        }

        // Create the memory dataspace.
        hid_t memDataspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = endingCount[0];
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));

        // Append the records.
        hid_t countsDataspaceHandle, timesDataspaceHandle;
        endingStart[0] = countsDims[0];
        countsDims[0] += endingCount[0];
        HDF5_EXCEPTION_CHECK(H5Dset_extent(countsDatasetHandle, countsDims));
        HDF5_EXCEPTION_CALL(countsDataspaceHandle,H5Dget_space(countsDatasetHandle));
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(countsDataspaceHandle, H5S_SELECT_SET, endingStart, NULL, endingCount, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(countsDatasetHandle, H5T_NATIVE_UINT32, memDataspaceHandle, countsDataspaceHandle, H5P_DEFAULT, newEndingCounts));
        HDF5_EXCEPTION_CHECK(H5Sclose(countsDataspaceHandle));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(timesDatasetHandle, countsDims));
        HDF5_EXCEPTION_CALL(timesDataspaceHandle,H5Dget_space(timesDatasetHandle));
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(timesDataspaceHandle, H5S_SELECT_SET, endingStart, NULL, endingCount, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(timesDatasetHandle, H5T_NATIVE_DOUBLE, memDataspaceHandle, timesDataspaceHandle, H5P_DEFAULT, newEndingTimes));
        HDF5_EXCEPTION_CHECK(H5Sclose(timesDataspaceHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));

        // Free the buffers.
        if (newEndingCounts != NULL) delete [] newEndingCounts;
        if (newEndingTimes != NULL) delete [] newEndingTimes;
    }

    // Otherwise, we need to insert some records before and/or after the existing block so we have to copy the old data and rebuild.
    else
    {
        // Declare the new buffer.
        uint32_t combinedMinCount = min(minCount,minNewCount);
        uint32_t combinedMaxCount = max(maxCount,maxNewCount);
        size_t newSize = combinedMaxCount-combinedMinCount+1;
        uint32_t * newCounts = new uint32_t[newSize];
        double * newTimes = new double[newSize];
        for (uint i=0; i<newSize; i++)
        {
            newCounts[i] = 0;
            newTimes[i] = -1.0;
        }

        // Copy the old data into the buffer.
        hid_t memDataspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = maxCount-minCount+1;
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_UINT32, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, &newCounts[minCount-combinedMinCount]));
        HDF5_EXCEPTION_CHECK(H5Dread(timesDatasetHandle, H5T_NATIVE_DOUBLE, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, &newTimes[minCount-combinedMinCount]));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));

        // Fill in the buffer with the new data.
        for (int i=0; i<firstPassageTimes->species_count_size(); i++)
        {
            uint32_t count = firstPassageTimes->species_count(i);
            if (count < minCount || count > maxCount)
            {
                newCounts[count-combinedMinCount] = count;
                newTimes[count-combinedMinCount] = firstPassageTimes->first_passage_time(i);
            }
            else
            {
                throw InvalidArgException("firstPassageTimes", "contained duplicates of existing counts", count);
            }
        }

        // Update the dataset with the combined data.
        memDims[0] = newSize;
        countsDims[0] = newSize;
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(countsDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(countsDatasetHandle, H5T_NATIVE_UINT32, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, newCounts));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(timesDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(timesDatasetHandle, H5T_NATIVE_DOUBLE, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, newTimes));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));

        // Free the buffers.
        if (newCounts != NULL) delete [] newCounts;
        if (newTimes != NULL) delete [] newTimes;
    }

    // Close any resources.
    HDF5_EXCEPTION_CHECK(H5Dclose(countsDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Dclose(timesDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(speciesGroupHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(fptGroupHandle));
}

/*void SimulationFile::appendSpatialModelObjects(unsigned int replicate, lm::io::SpatialModel * model)
{
    if (model->sphere_xc_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");
    if (model->sphere_yc_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");
    if (model->sphere_zc_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");
    if (model->sphere_radius_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the group.
    hid_t modelHandle, spatialHandle;
    if ((modelHandle=H5Gopen2(replicateHandles->group, "Model", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(modelHandle,H5Gcreate2(replicateHandles->group, "Model", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }
    if ((spatialHandle=H5Gopen2(modelHandle, "Spatial", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(spatialHandle,H5Gcreate2(modelHandle, "Spatial", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    // Open or create the datasets.
    hid_t spheresDatasetHandle;
    if ((spheresDatasetHandle=H5Dopen2(spatialHandle, "Spheres", H5P_DEFAULT)) < 0)
    {
        // Create the datasets.
        uint RANK=2;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        dims[1] = 5;
        maxDims[0] = H5S_UNLIMITED;
        maxDims[1] = 5;
        chunkDims[0] = 100;
        chunkDims[1] = 5;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(spheresDatasetHandle,H5Dcreate2(spatialHandle, "Spheres", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Get the current size of the dataset.
    unsigned int RANK=2;
    hsize_t dims[RANK];
    hid_t dataspace_id;
    int result;
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(spheresDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

    // Extend the dataset by the number of rows in the data set.
    dims[0] += model->sphere_type_size();
    HDF5_EXCEPTION_CHECK(H5Dset_extent(spheresDatasetHandle, dims));

    // Create the memory dataset.
    hid_t memspace_id;
    hsize_t memDims[RANK];
    memDims[0] = model->sphere_type_size();
    memDims[1] = 1;
    HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

    // Write the values.
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(spheresDatasetHandle));
    hsize_t start[RANK], count[RANK];
    start[0] = dims[0]-model->sphere_type_size();
    count[0] = model->sphere_type_size();
    start[1] = 0;
    count[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_xc().data()));
    start[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_yc().data()));
    start[1] = 2;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_zc().data()));
    start[1] = 3;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_radius().data()));
    start[1] = 4;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_type().data()));

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    HDF5_EXCEPTION_CHECK(H5Dclose(spheresDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(modelHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(spatialHandle));
}

void SimulationFile::getSpatialModelObjects(unsigned int replicate, lm::io::SpatialModel * model)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open the group.
    hid_t modelHandle, spatialHandle;
    HDF5_EXCEPTION_CALL(modelHandle,H5Gopen2(replicateHandles->group, "Model", H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(spatialHandle,H5Gopen2(modelHandle, "Spatial", H5P_DEFAULT));

    // Open the datasets.
    hid_t spheresDatasetHandle;
    HDF5_EXCEPTION_CALL(spheresDatasetHandle,H5Dopen2(spatialHandle, "Spheres", H5P_DEFAULT));

    // Get the current size of the dataset.
    const uint RANK=2;
    const uint NUM_COLS=5;
    hsize_t dims[RANK];
    hid_t dataspace_id;
    int result;
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(spheresDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));

    // Create the memory dataset.
    hid_t memspace_id;
    hsize_t memDims[RANK];
    memDims[0] = TUNE_SPATIAL_MODEL_OBJECT_READ_BUFFER_SIZE;
    memDims[1] = NUM_COLS;
    HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));
    double * memData = new double[TUNE_SPATIAL_MODEL_OBJECT_READ_BUFFER_SIZE*NUM_COLS];

    // Loop over the data and read it in one section at a time.
    hsize_t start[RANK], count[RANK];
    start[0] = 0;
    start[1] = 0;
    count[0] = TUNE_SPATIAL_MODEL_OBJECT_READ_BUFFER_SIZE;
    count[1] = NUM_COLS;
    while (start[0] < dims[0])
    {
        // Figure out how many rows to read.
        if (start[0]+count[0] > dims[0])
        {
            // This is the last read, so modify its size.
            count[0] = dims[0]-start[0];

            // Modify the memory dataspace as well.
            HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
            memDims[0] = count[0];
            HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));
        }

        // Read the rows.
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, (void *)memData));

        // Add them to the model.
        for (uint index=0; index<count[0]*NUM_COLS; )
        {
            model->add_sphere_xc(memData[index++]);
            model->add_sphere_yc(memData[index++]);
            model->add_sphere_zc(memData[index++]);
            model->add_sphere_radius(memData[index++]);
            model->add_sphere_type(memData[index++]);
        }

        // Move to the next section.
        start[0] += count[0];
    }

    // Cleanup some resources.
    if (memData != NULL) delete []  memData; memData = NULL;
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    HDF5_EXCEPTION_CHECK(H5Dclose(spheresDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(modelHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(spatialHandle));
}*/


vector<double> SimulationFile::getLatticeTimes(unsigned int replicate)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Read the initial species counts.
    int RANK=1, rank;
    hsize_t dims[RANK];
    H5T_class_t type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(replicateHandles->group, "LatticeTimes", &rank));
    if (rank != RANK) throw Exception("Invalid dataset rank", filename.c_str(), "LatticeTimes");
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(replicateHandles->group, "LatticeTimes", dims, &type, &size));
    if (size != sizeof(double)) throw Exception("Invalid dataset type", filename.c_str(), "LatticeTimes");
    double * timesBuffer = new double[dims[0]];
    HDF5_EXCEPTION_CHECK(H5LTread_dataset_double(replicateHandles->group, "LatticeTimes", timesBuffer));
    vector<double> times;
    for (uint i=0; i<dims[0]; i++)
        times.push_back(timesBuffer[i]);
    delete [] timesBuffer;
    return times;
}

vector<double> SimulationFile::getSpeciesCountTimes(unsigned int replicate)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Read the initial species counts.
    int RANK=1, rank;
    hsize_t dims[RANK];
    H5T_class_t type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(replicateHandles->group, "SpeciesCountTimes", &rank));
    if (rank != RANK) throw Exception("Invalid dataset rank", filename.c_str(), "SpeciesCountTimes");
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(replicateHandles->group, "SpeciesCountTimes", dims, &type, &size));
    if (size != sizeof(double)) throw Exception("Invalid dataset type", filename.c_str(), "SpeciesCountTimes");
    double * timesBuffer = new double[dims[0]];
    HDF5_EXCEPTION_CHECK(H5LTread_dataset_double(replicateHandles->group, "SpeciesCountTimes", timesBuffer));
    vector<double> times;
    for (uint i=0; i<dims[0]; i++)
        times.push_back(timesBuffer[i]);
    delete [] timesBuffer;
    return times;
}

std::map<double, vector<int> > SimulationFile::getSpeciesCounts(unsigned int replicate)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);
	std::map<double, vector<int> > particles;

    // Read the initial species counts.
    int RANK=2, rank;
    hsize_t dims[RANK];
    H5T_class_t type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(replicateHandles->group, "SpeciesCounts", &rank));
    if (rank != RANK) throw Exception("Invalid dataset rank", filename.c_str(), "SpeciesCounts");
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(replicateHandles->group, "SpeciesCounts", dims, &type, &size));
    if (size != sizeof(int)) throw Exception("Invalid dataset type", filename.c_str(), "SpeciesCounts");
    int *counts = new int[dims[0]*dims[1]];
    HDF5_EXCEPTION_CHECK(H5LTread_dataset_int(replicateHandles->group, "SpeciesCounts", counts));

    vector<double> times=getSpeciesCountTimes(replicate);

    for (uint i=0; i<dims[0]; i++)
	{
		//printf("counts at time %g (frame %d):\n", times[i], i);
		vector<int> p(dims[1]+1);
		for (uint j=0; j<dims[1]; j++)
		{
			//printf("%d ", counts[i*dims[1] + j]);
			// expecting 1-based vector of particles
			p[j+1]=counts[i*dims[1] + j];
		}

		particles[times[i]]=p;
	}

    delete [] counts;
    return particles;
}

void SimulationFile::getLattice(unsigned int replicate, unsigned int latticeIndex, lm::rdme::Lattice * lattice)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Read the lattice data.
    char latticeDatasetName[19];
    snprintf(latticeDatasetName, sizeof(latticeDatasetName), "Lattice/%010d", latticeIndex);
    int RANK=4, rank;
    hsize_t dims[RANK];
    H5T_class_t type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(replicateHandles->group, latticeDatasetName, &rank));
    if (rank != RANK) throw Exception("Invalid dataset rank", filename.c_str(), latticeDatasetName);
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(replicateHandles->group, latticeDatasetName, dims, &type, &size));
    if (lattice->getSize().x != dims[0] || lattice->getSize().y != dims[1] || lattice->getSize().z != dims[2] || lattice->getMaxOccupancy() != dims[3]) throw Exception("Invalid lattice dimensions", filename.c_str());
    size_t particlesBufferSize = dims[0]*dims[1]*dims[2]*dims[3];

	switch(size)
	{
		case sizeof(uint8_t):
		{
			uint8_t * particlesBuffer = new uint8_t[particlesBufferSize];
			HDF5_EXCEPTION_CHECK(H5LTread_dataset(replicateHandles->group, latticeDatasetName, H5T_STD_U8LE, particlesBuffer));
			lattice->setFromRowMajorByteData(particlesBuffer, particlesBufferSize);
			delete [] particlesBuffer;
			break;
		}
		case sizeof(int):
		{
			unsigned int *particlesBuffer = new unsigned int[particlesBufferSize];
			HDF5_EXCEPTION_CHECK(H5LTread_dataset(replicateHandles->group, latticeDatasetName, H5T_STD_U32LE, particlesBuffer));
			lm::rdme::IntLattice *intlattice = (lm::rdme::IntLattice*)lattice;
			intlattice->setFromRowMajorData(particlesBuffer, particlesBufferSize);
			delete [] particlesBuffer;
			break;
		}
		default:
			throw Exception("Unsupported dataset size",size);
	}
}

// Get particle counts without actually copying to a lattice
std::map<uint32_t, uint> SimulationFile::getParticleCounts(unsigned int replicate, unsigned int latticeIndex)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Read the lattice data.
    char latticeDatasetName[19];
    snprintf(latticeDatasetName, sizeof(latticeDatasetName), "Lattice/%010d", latticeIndex);
    int RANK=4, rank;
    hsize_t dims[RANK];
    H5T_class_t type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(replicateHandles->group, latticeDatasetName, &rank));
    if (rank != RANK) throw Exception("Invalid dataset rank", filename.c_str(), latticeDatasetName);
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(replicateHandles->group, latticeDatasetName, dims, &type, &size));
    size_t particlesBufferSize = dims[0]*dims[1]*dims[2]*dims[3];

	std::map<uint32_t, uint> particles;
	switch(size)
	{
		case sizeof(uint8_t):
		{
			uint8_t * particlesBuffer = new uint8_t[particlesBufferSize];
			HDF5_EXCEPTION_CHECK(H5LTread_dataset(replicateHandles->group, latticeDatasetName, H5T_STD_U8LE, particlesBuffer));
			for(size_t i=0; i<particlesBufferSize; i++)
			{
				const particle_t p=particlesBuffer[i];
				if(p > 0)
				{
					if(particles.count(p) == 0)
						particles[p]=1;
					else
						particles[p]++;
				}
			}
			delete [] particlesBuffer;
			break;
		}
		case sizeof(int):
		{
			unsigned int *particlesBuffer = new unsigned int[particlesBufferSize];
			HDF5_EXCEPTION_CHECK(H5LTread_dataset(replicateHandles->group, latticeDatasetName, H5T_STD_U32LE, particlesBuffer));
			for(size_t i=0; i<particlesBufferSize; i++)
			{
				const particle_t p=particlesBuffer[i];
				if(p > 0)
				{
					if(particles.count(p) == 0)
						particles[p]=1;
					else
						particles[p]++;
				}
			}
			delete [] particlesBuffer;
			break;
		}
		default:
			throw Exception("Unsupported dataset size",size);
	}

	return particles;
}

void SimulationFile::closeReplicate(unsigned int replicate)
{
    map<unsigned int,ReplicateHandles *>::iterator it = openReplicates.find(replicate);
    if (it != openReplicates.end())
    {
        ReplicateHandles * handles = it->second;
        closeReplicateHandles(handles);
        delete handles;
        openReplicates.erase(it);
    }
}

void SimulationFile::closeAllReplicates()
{
    for (map<unsigned int,ReplicateHandles *>::iterator it=openReplicates.begin(); it != openReplicates.end(); it++)
    {
        ReplicateHandles * handles = it->second;
        closeReplicateHandles(handles);
        delete handles;
    }
    openReplicates.clear();
}


SimulationFile::ReplicateHandles * SimulationFile::openReplicateHandles(unsigned int replicate)
{
    // See if the replicate is already open.
    map<unsigned int,ReplicateHandles *>::iterator it = openReplicates.find(replicate);
    if (it == openReplicates.end())
    {
        ReplicateHandles * handles;

        // Construct a string representation of the replicate name.
        std::stringstream ss;
        ss.fill('0');
        ss.width(7);
        ss << replicate;
        string replicateString = ss.str();

        // Turn off exception handling again to work around odd behavior in hdf5 1.8.4.
        HDF5_EXCEPTION_CHECK(H5Eset_auto2(H5E_DEFAULT, NULL, NULL));

        // Open the replicate group.
        hid_t groupHandle;
        if ((groupHandle=H5Gopen2(simulationsGroup, replicateString.c_str(), H5P_DEFAULT)) >= 0)
        {
            // Open the rest of the handles.
            handles = new ReplicateHandles();
            handles->group = groupHandle;
            HDF5_EXCEPTION_CALL(handles->speciesCountsDataset,H5Dopen2(handles->group, "SpeciesCounts", H5P_DEFAULT));
            HDF5_EXCEPTION_CALL(handles->speciesCountTimesDataset,H5Dopen2(handles->group, "SpeciesCountTimes", H5P_DEFAULT));
        }
        else
        {
            // If we couldn't open it, try to create a new one.
            handles = createReplicateHandles(replicateString);
        }

        // Add it to the map.
        openReplicates[replicate] = handles;
        return handles;
    }
    else
    {
        // Use the already open replicate.
        return it->second;
    }

}

SimulationFile::ReplicateHandles * SimulationFile::createReplicateHandles(string replicateString)
{
    // Make sure the model is loaded, since need the number of species.
    loadModel();

    ReplicateHandles * handles = new ReplicateHandles();

    // Create the group.
    HDF5_EXCEPTION_CALL(handles->group,H5Gcreate2(simulationsGroup, replicateString.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // Create the species counts dataset
    {
        unsigned int RANK=2;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        dims[1] = numberSpecies;
        maxDims[0] = H5S_UNLIMITED;
        maxDims[1] = numberSpecies;
        chunkDims[0] = 100;
        chunkDims[1] = numberSpecies;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(handles->speciesCountsDataset,H5Dcreate2(handles->group, "SpeciesCounts", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Create the species count times dataset
    {
        unsigned int RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(handles->speciesCountTimesDataset,H5Dcreate2(handles->group, "SpeciesCountTimes", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    return handles;
}

void SimulationFile::closeReplicateHandles(ReplicateHandles * handles)
{
    HDF5_EXCEPTION_CHECK(H5Gclose(handles->group));
    HDF5_EXCEPTION_CHECK(H5Dclose(handles->speciesCountsDataset));
    HDF5_EXCEPTION_CHECK(H5Dclose(handles->speciesCountTimesDataset));
    handles->group = H5I_INVALID_HID;
    handles->speciesCountsDataset = H5I_INVALID_HID;
    handles->speciesCountTimesDataset = H5I_INVALID_HID;
}

std::map<uint, string> SimulationFile::getSpeciesNames()
{
	std::map<uint, string> names;
	string namelist=getParameter("speciesNames", "");
	if(namelist.length()==0)
		return names;

	size_t pos=0,nextpos=0;
	uint particle=1;
	while((nextpos=namelist.find(",", pos)) != string::npos)
	{
		names[particle]=namelist.substr(pos, nextpos-pos);	
	//	printf("name %d %s\n", particle, names[particle].c_str());
		particle++;
		pos=nextpos+1;
	}

	names[particle]=namelist.substr(pos);

	return names;
}

std::map<uint, string> SimulationFile::getSiteTypeNames()
{
	std::map<uint, string> names;
	string namelist=getParameter("siteTypeNames", "");
	if(namelist.length()==0)
		return names;

	size_t pos=0,nextpos=0;
	uint siteid=0;
	while((nextpos=namelist.find(",", pos)) != string::npos)
	{
		names[siteid]=namelist.substr(pos, nextpos-pos);	
	//	printf("name %d %s\n", particle, names[particle].c_str());
		siteid++;
		pos=nextpos+1;
	}

	names[siteid]=namelist.substr(pos);

	return names;
}

void
SimulationFile::arbitraryH5(byte *ptr) {
    H5MetaData* header = reinterpret_cast<H5MetaData*>(ptr);
    byte* data = ptr+sizeof(H5MetaData);
    switch (header->mode) {
        case H5MetaData::APPEND_TO_GROUP:
            appendDsToGroup(header, data);
            break;
        case H5MetaData::APPEND_TO_DATASET:
            appendToDataset(header, data);
            break;
        case H5MetaData::NEW_DATASET:
            makeNewDataset(header, data);
            break;
        default:
            throw InvalidArgException("H5MetaData::mode", "Write mode invalid");
    }
}


void
SimulationFile::appendToDataset(H5MetaData *header, byte *data) {
    hid_t replicateHandle = openReplicateHandles(header->replicate)->group;

    hsize_t dims[H5S_MAX_RANK], chunk[H5S_MAX_RANK], maxDims[H5S_MAX_RANK];
    hsize_t memDims[H5S_MAX_RANK], start[H5S_MAX_RANK], count[H5S_MAX_RANK];
    hid_t datasetHandle;
    hid_t dataspaceHandle;

    if (H5Lexists(replicateHandle, header->name, H5P_DEFAULT) > 0) {
        HDF5_EXCEPTION_CALL(datasetHandle, H5Dopen2(replicateHandle,  header->name, H5P_DEFAULT));
    } else {
        for (int i=1; i <= header->ndim; i++) {
            dims[i] = header->shape[i-1];
            chunk[i] = min((size_t)TUNE_LATTICE_GZIP_CHUNK_SIZE,header->shape[i-1]);
            maxDims[i] = header->shape[i-1];
        }

        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunk[0] = 100;

        hid_t propsHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(header->ndim+1, dims, maxDims));
        HDF5_EXCEPTION_CALL(propsHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(propsHandle, header->ndim+1, chunk));
        HDF5_EXCEPTION_CALL(datasetHandle, H5Dcreate2(replicateHandle, header->name, header->h5type, dataspaceHandle, H5P_DEFAULT, propsHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(propsHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    }

    HDF5_EXCEPTION_CALL(dataspaceHandle, H5Dget_space(datasetHandle));
    HDF5_EXCEPTION_CHECK(H5Sget_simple_extent_dims(dataspaceHandle, dims, NULL));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));

    // Extend the dataset by the number of rows in the data set.
    dims[0] += 1;
    HDF5_EXCEPTION_CHECK(H5Dset_extent(datasetHandle, dims));

    // Create the memory dataset.
    hid_t memspaceHandle;
    memcpy(memDims, dims, sizeof(dims));
    memDims[0] = 1;
    HDF5_EXCEPTION_CALL(memspaceHandle, H5Screate_simple(header->ndim+1, memDims, NULL));

    // Write the new data.
    HDF5_EXCEPTION_CALL(dataspaceHandle, H5Dget_space(datasetHandle));

    memset(start, 0, sizeof(start));
    start[0] = dims[0]-1;
    memcpy(count, memDims, sizeof(memDims));

    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, header->h5type, memspaceHandle, dataspaceHandle, H5P_DEFAULT, data));

    HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    HDF5_EXCEPTION_CHECK(H5Sclose(memspaceHandle));
    HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
}


void
SimulationFile::makeNewDataset(H5MetaData *header, byte *data) {
    hid_t replicateHandle = openReplicateHandles(header->replicate)->group;

    if (H5Lexists(replicateHandle, header->name, H5P_DEFAULT) > 0) {
        throw Exception("Dataset already exists");
    }

    hsize_t dims[H5S_MAX_RANK], chunk[H5S_MAX_RANK];

    for (int i=0; i < header->ndim; i++) {
        dims[i] = header->shape[i];
        chunk[i] = min((size_t)TUNE_LATTICE_GZIP_CHUNK_SIZE,header->shape[i]);
    }

    hid_t dataspaceHandle, dcplHandle, datasetHandle;
    HDF5_EXCEPTION_CALL(dataspaceHandle, H5Screate_simple(header->ndim, dims, NULL));
    HDF5_EXCEPTION_CALL(dcplHandle, H5Pcreate(H5P_DATASET_CREATE));
    HDF5_EXCEPTION_CHECK(H5Pset_deflate(dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
    HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, header->ndim, chunk));
    HDF5_EXCEPTION_CALL(datasetHandle, H5Dcreate2(replicateHandle, header->name, header->h5type, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, header->h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));

    HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
    HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
}

void
SimulationFile::appendDsToGroup(H5MetaData *header, byte *data) {
    hid_t replicateHandle = openReplicateHandles(header->replicate)->group;

    hid_t groupHandle;
    if (H5Lexists(replicateHandle, header->name, H5P_DEFAULT) > 0) {
        HDF5_EXCEPTION_CALL(groupHandle, H5Gopen2(replicateHandle, header->name, H5P_DEFAULT));
    } else {
        HDF5_EXCEPTION_CALL(groupHandle,H5Gcreate2(replicateHandle, header->name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    hsize_t nds;
    HDF5_EXCEPTION_CHECK(H5Gget_num_objs(groupHandle, &nds));

    char dsName[11];
    snprintf(dsName, sizeof(dsName), "%010zd", nds);

    hsize_t dims[H5S_MAX_RANK], chunk[H5S_MAX_RANK];

    for (int i=0; i < header->ndim; i++) {
        dims[i] = header->shape[i];
        chunk[i] = min((size_t)TUNE_LATTICE_GZIP_CHUNK_SIZE,header->shape[i]);
    }

    hid_t dataspaceHandle, dcplHandle, datasetHandle;
    HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(header->ndim, dims, NULL));
    HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
    HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
    HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, header->ndim, chunk));
    HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(groupHandle, dsName, header->h5type, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, header->h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));

    HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
    HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(groupHandle));
}

void
SimulationFile::arbitraryH5Lookup(byte *req, byte *&data, size_t &sz) {
    H5Lookup* header = reinterpret_cast<H5Lookup*>(req);
    switch (header->mode) {
        case H5Lookup::DATASET:
            readDataset(header, data, sz);
            break;
        case H5Lookup::ATTR:
            readAttr(header, data, sz);
            break;
        default:
            throw InvalidArgException("H5Lookup::mode", "H5 data mode invalid");
    }
}

void 
SimulationFile::readAttr(H5Lookup *req, byte *&buf, size_t &sz) {
    H5MetaData md;
    md.replicate = req->replicate;
    hsize_t dims[H5S_MAX_RANK];
    hsize_t elemSz;

    hid_t attrHandle, typeHandle, dataspaceHandle;
    HDF5_EXCEPTION_CALL(attrHandle, H5Aopen_by_name(file, req->path, req->attr, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(dataspaceHandle, H5Aget_space(attrHandle));
    HDF5_EXCEPTION_CHECK(H5Sget_simple_extent_dims(dataspaceHandle, dims, NULL));
    HDF5_EXCEPTION_CALL(md.ndim, H5Sget_simple_extent_ndims(dataspaceHandle));
    HDF5_EXCEPTION_CALL(typeHandle, H5Aget_type(attrHandle));
    HDF5_EXCEPTION_CALL(md.h5type, H5Tget_native_type(typeHandle, H5T_DIR_ASCEND));
    HDF5_EXCEPTION_CALL(elemSz, H5Tget_size(md.h5type));

    for (int i=0; i < md.ndim; i++) {
        md.shape[i] = dims[i];
    }

    md.payloadSize = elemSz*std::accumulate(md.shape, md.shape + md.ndim, 1, std::multiplies<size_t>());

    sz = sizeof(md)+md.payloadSize;
    buf = new unsigned char[sz];
    H5MetaData *hdr = reinterpret_cast<H5MetaData*>(buf);
    unsigned char *payload = buf + sizeof(md);
    *hdr = md;

    HDF5_EXCEPTION_CHECK(H5Aread(attrHandle, md.h5type, (void *)payload));

    HDF5_EXCEPTION_CHECK(H5Tclose(typeHandle));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    HDF5_EXCEPTION_CHECK(H5Aclose(attrHandle));

    if (H5Tequal(req->h5type, md.h5type) <= 0) {
        throw Exception("Requested and HDF5 data types do not match");
    }
}


void
SimulationFile::readDataset(H5Lookup *req, byte *&buf, size_t &sz) {
    H5MetaData md;
    md.replicate = req->replicate;

    hid_t datasetHandle, dataspaceHandle, typeHandle;
    H5T_class_t typeClass;
    hsize_t dims[H5S_MAX_RANK];
    hsize_t elemSz;

    HDF5_EXCEPTION_CALL(datasetHandle, H5Dopen2(file, req->path, H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(dataspaceHandle, H5Dget_space(datasetHandle));
    HDF5_EXCEPTION_CHECK(H5Sget_simple_extent_dims(dataspaceHandle, dims, NULL));
    HDF5_EXCEPTION_CALL(md.ndim, H5Sget_simple_extent_ndims(dataspaceHandle));
    HDF5_EXCEPTION_CALL(typeHandle, H5Dget_type(datasetHandle));
    HDF5_EXCEPTION_CALL(md.h5type, H5Tget_native_type(typeHandle, H5T_DIR_ASCEND));
    HDF5_EXCEPTION_CALL(elemSz, H5Tget_size(md.h5type));

    for (int i=0; i < md.ndim; i++) {
        md.shape[i] = dims[i];
    }

    md.payloadSize = elemSz*std::accumulate(md.shape, md.shape + md.ndim, 1, std::multiplies<size_t>());

    sz = sizeof(md)+md.payloadSize;
    buf = new unsigned char[sz];
    H5MetaData *hdr = reinterpret_cast<H5MetaData*>(buf);
    unsigned char *payload = buf + sizeof(md);
    *hdr = md;

    HDF5_EXCEPTION_CHECK(H5Dread(datasetHandle, md.h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)payload));
    HDF5_EXCEPTION_CHECK(H5Tclose(typeHandle));
    HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));

    if (H5Tequal(req->h5type, md.h5type) <= 0) {
        throw Exception("Requested and HDF5 data types do not match");
    }
}


}
}
}

/* Generic code to read a dataset with a hyperslab.
 * int ndims;
hsize_t dims[4];
hsize_t memdims[1];
hsize_t start[4], count[4];
hid_t dataset, dataspace, memspace;
HDF5_EXCEPTION_CALL(dataset,H5Dopen2(file, "/Model/Diffusion/Lattice", H5P_DEFAULT));
HDF5_EXCEPTION_CALL(dataspace,H5Dget_space(dataset));
ndims=H5Sget_simple_extent_ndims(dataspace);
if (ndims != 4) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
HDF5_EXCEPTION_CALL(ndims,H5Sget_simple_extent_dims(dataspace, dims, NULL));
if (dims[0]*dims[1]*dims[2]*dims[3] != latticeSize) throw InvalidArgException("lattice", "incorrect lattice size");
start[0] = 0;
start[1] = 0;
start[2] = 0;
start[3] = 0;
count[0] = dims[0];
count[1] = dims[1];
count[2] = dims[2];
count[3] = dims[3];
HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL));
memdims[0] = dims[0]*dims[1]*dims[2]*dims[3];
HDF5_EXCEPTION_CALL(memspace,H5Screate_simple(1, memdims, NULL));
HDF5_EXCEPTION_CHECK(H5Dread(dataset, H5T_NATIVE_UINT8, memspace, dataspace, H5P_DEFAULT, (void *)lattice));
HDF5_EXCEPTION_CHECK(H5Sclose(memspace));
HDF5_EXCEPTION_CHECK(H5Sclose(dataspace));
HDF5_EXCEPTION_CHECK(H5Dclose(dataset));
*/
