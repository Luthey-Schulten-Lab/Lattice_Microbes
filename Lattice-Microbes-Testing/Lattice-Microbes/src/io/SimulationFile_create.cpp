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
#include <sys/stat.h>
#include "config.h"
#include "core/Exceptions.h"
#include "io/lm_hdf5.h"
#include "io/SimulationFile.h"

using std::string;
using lm::IOException;

namespace lm {
namespace io {
namespace hdf5 {

void SimulationFile::create(const string filename)
{
    create(filename.c_str(), false);
}

void SimulationFile::create(const char * filename)
{
    create(filename, false);
}

void SimulationFile::create(const string filename, unsigned int numberSpecies)
{
    create(filename.c_str(), true, numberSpecies);
}

void SimulationFile::create(const char * filename, unsigned int numberSpecies)
{
    create(filename, true, numberSpecies);
}

void SimulationFile::create(const char * filename, bool initializeModel, unsigned int numberSpecies)
{
    // If the file exists, open it.
    struct stat fileStats;
    if (stat(filename, &fileStats) == 0) throw lm::IOException("a file with the specified filename already exists", filename);

    // Set the size of the user block.
    hid_t creationProperties;
    HDF5_EXCEPTION_CALL(creationProperties,H5Pcreate(H5P_FILE_CREATE));
    HDF5_EXCEPTION_CHECK(H5Pset_userblock(creationProperties, 512));

    // Create the file.
    hid_t file;
    HDF5_EXCEPTION_CALL(file,H5Fcreate(filename, H5F_ACC_EXCL, creationProperties, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Pclose(creationProperties));

    // Write the version.
    unsigned int version=CURRENT_VERSION;
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/", "formatVersion", &version, 1));

    // Create the main groups.
    hid_t group;
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Simulations", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));

    // See if we should initialize the model.
    if (initializeModel)
    {
        // Write the model.
        HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Reaction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Gclose(group));
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Reaction", "numberSpecies", &numberSpecies, 1));
    }

    // Close the file.
    HDF5_EXCEPTION_CHECK(H5Fclose(file));

    // Write the file magic.
    char magic[] = {'L','M','H','5'};
    FILE * fp = NULL;
    fp=fopen(filename, "rb+");
    if (fp == NULL) throw lm::IOException("the specified file could not be opened", filename);
    rewind(fp);
    fwrite (magic, 1, sizeof(magic), fp);
    if (fclose(fp) != 0)  throw lm::IOException("the specified file could not be closed", filename);
}

}
}
}


