/*
 * University of Illinois Open Source License
 * Copyright 2010 Luthey-Schulten Group, 
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
#include <map>
#include <vector>
#include "TestingHelper.h"
#include "lm/Exceptions.h"
#include "lm/io/DiffusionModel.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/ParameterValues.pb.h"
#include "lm/io/ReactionModel.pb.h"
#include "lm/io/SpatialModel.pb.h"
#include "lm/io/hdf5/HDF5.h"
#include "lm/io/hdf5/SimulationFile.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

using std::string;
using std::map;
using std::vector;
using lm::IOException;
using lm::io::DiffusionModel;
using lm::io::FirstPassageTimes;
using lm::io::ParameterValues;
using lm::io::ReactionModel;
using lm::io::SpatialModel;
using lm::io::hdf5::HDF5Exception;
using lm::io::hdf5::SimulationFile;


class SimulationFileTester : public lm::io::hdf5::SimulationFile
{
public:
    SimulationFileTester(string s):SimulationFile(s),openReplicatesP(&openReplicates),versionP(&version),numberSpeciesP(&numberSpecies){}
    virtual ReplicateHandles * openReplicateHandles(unsigned int replicate) throw(HDF5Exception) {return SimulationFile::openReplicateHandles(replicate);}
    virtual ReplicateHandles * createReplicateHandles(string replicateString) throw(HDF5Exception) {return SimulationFile::createReplicateHandles(replicateString);}
    virtual void closeReplicateHandles(ReplicateHandles * handles) throw(HDF5Exception) {SimulationFile::closeReplicateHandles(handles);}
    map<unsigned int,ReplicateHandles *> * openReplicatesP;
    unsigned int * versionP;
    unsigned int * numberSpeciesP;
};

BOOST_AUTO_TEST_SUITE(SimulationFileTest)

BOOST_AUTO_TEST_CASE(OpenNonHDF5File)
{
    SimulationFile * f = NULL;
    BOOST_CHECK_THROW(f=new SimulationFile(DATA_DIR+"/lm/io/hdf5/SimulationFile_1.h5"), IOException);
    if (f != NULL) delete f; f = NULL;
    BOOST_CHECK_THROW(f=new SimulationFile(DATA_DIR+"/lm/io/hdf5/"), IOException);
    if (f != NULL) delete f; f = NULL;
}

BOOST_AUTO_TEST_CASE(OpenGoodFile)
{
    SimulationFileTester * f = NULL;
    f=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_2.h5");

    // Make sure the file attributes were opened correctly.
    BOOST_CHECK_EQUAL(*f->versionP, 2U);

    if (f != NULL) delete f; f = NULL;
}

BOOST_AUTO_TEST_CASE(CreateFile)
{
    string filename = TMP_DIR+"/lmtest_SimulationFile_CreateFile.h5";
    remove(filename.c_str());

    // Create a new file.
    SimulationFile::create(filename);

    // Make sure we can open the file.
    SimulationFileTester * f = NULL;
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(filename));
    BOOST_CHECK_EQUAL(*f->versionP, 4);
    delete f; f = NULL;

    remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(GetParameterV2)
{
    SimulationFileTester * f = NULL;
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_2.h5"));
    BOOST_CHECK_EQUAL(*f->versionP, 2U);
    BOOST_CHECK_CLOSE(atof(f->getParameter("unused1","11.0").c_str()), 11.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test1").c_str()), 2.3, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test2").c_str()), 5.7e-16, 0.01);
    if (f != NULL) delete f; f = NULL;
}

BOOST_AUTO_TEST_CASE(GetParameterV3)
{
    SimulationFileTester * f = NULL;
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_3.h5"));
    BOOST_CHECK_EQUAL(*f->versionP, 3U);
    BOOST_CHECK_CLOSE(atof(f->getParameter("unused1","11.0").c_str()), 11.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test1").c_str()), 2.3, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test2").c_str()), 5.7e-16, 0.01);
    BOOST_CHECK_EQUAL(f->getParameter("test3"), "testing");
    BOOST_CHECK_EQUAL(f->getParameter("test4"), "");
    if (f != NULL) delete f; f = NULL;
}

BOOST_AUTO_TEST_CASE(SetParameter)
{
    // Create a new file with some parameters.
    string filename = TMP_DIR+"/lmtest_SimulationFile_SetParameter.h5";
    remove(filename.c_str());
    SimulationFile::create(filename);
    SimulationFileTester * f = NULL;
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(filename));
    BOOST_CHECK_NO_THROW(f->setParameter("test1", "23.0"));
    BOOST_CHECK_NO_THROW(f->setParameter("test2", "24.0"));
    BOOST_CHECK_NO_THROW(f->setParameter("test3", "1.7e7"));
    BOOST_CHECK_NO_THROW(f->setParameter("test4", "newtest"));
    BOOST_CHECK_NO_THROW(f->close());
    delete f; f = NULL;

    // Make sure the parameters were saved in the file.
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(filename));
    BOOST_CHECK_CLOSE(atof(f->getParameter("test1").c_str()), 23.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test2").c_str()), 24.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test3").c_str()), 1.7e7, 0.01);
    BOOST_CHECK_EQUAL(f->getParameter("test4"), "newtest");

    // Change the parameters.
    BOOST_CHECK_NO_THROW(f->setParameter("test2", "-25.0"));
    BOOST_CHECK_NO_THROW(f->setParameter("test3", "1.8e-7"));
    BOOST_CHECK_NO_THROW(f->setParameter("test4", "newtest2"));
    BOOST_CHECK_NO_THROW(f->close());
    delete f; f = NULL;

    // Make sure the changed parameters were saved in the file.
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(filename));
    BOOST_CHECK_CLOSE(atof(f->getParameter("test1").c_str()), 23.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test2").c_str()), -25.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test3").c_str()), 1.8e-7, 0.01);
    BOOST_CHECK_EQUAL(f->getParameter("test4"), "newtest2");
    delete f; f = NULL;

    remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(Checkpoint)
{
    // Create a new file with some parameters.
    string filename = TMP_DIR+"/lmtest_SimulationFile_Checkpoint.h5";
    remove(filename.c_str());
    SimulationFile::create(filename);
    SimulationFileTester * f = NULL;
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(filename));
    BOOST_CHECK_NO_THROW(f->setParameter("test1", "23.0"));
    BOOST_CHECK_NO_THROW(f->setParameter("test2", "24.0"));
    BOOST_CHECK_NO_THROW(f->setParameter("test3", "1.7e7"));

    // Checkpoint the file.
    string checkpointFilename;
    BOOST_REQUIRE_NO_THROW(checkpointFilename=f->checkpoint());
    BOOST_CHECK_EQUAL(checkpointFilename, filename+".chk");

    // Make sure we can still change the file.
    BOOST_CHECK_NO_THROW(f->setParameter("test1", "23.0"));
    BOOST_CHECK_NO_THROW(f->setParameter("test2", "-25.0"));
    BOOST_CHECK_NO_THROW(f->setParameter("test3", "1.8e-7"));
    delete f; f = NULL;

    // Make sure the changed parameters were saved in the file.
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(filename));
    BOOST_CHECK_CLOSE(atof(f->getParameter("test1").c_str()), 23.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test2").c_str()), -25.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test3").c_str()), 1.8e-7, 0.01);
    delete f; f = NULL;

    // Make sure the checkpoint file is correct.
    BOOST_REQUIRE_NO_THROW(f=new SimulationFileTester(checkpointFilename));
    BOOST_CHECK_CLOSE(atof(f->getParameter("test1").c_str()), 23.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test2").c_str()), 24.0, 0.01);
    BOOST_CHECK_CLOSE(atof(f->getParameter("test3").c_str()), 1.7e7, 0.01);
    delete f; f = NULL;

    remove(checkpointFilename.c_str());
    remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(CreateReplicateHandles)
{
    string filename = TMP_DIR+"/lmtest_SimulationFile_CreateReplicateHandles.h5";
    remove(filename.c_str());

    SimulationFileTester * f = NULL;

    // Create a new file.
    SimulationFile::create(filename, 10);

    // Create some new replicates.
    SimulationFile::ReplicateHandles * result;
    f=new SimulationFileTester(filename);
    result = f->createReplicateHandles("00");
    BOOST_CHECK_GE(result->group, 0);
    BOOST_CHECK_GE(result->speciesCountsDataset, 0);
    BOOST_CHECK_NO_THROW(f->closeReplicateHandles(result));
    BOOST_CHECK_EQUAL(result->group, H5I_INVALID_HID);
    BOOST_CHECK_EQUAL(result->speciesCountsDataset, H5I_INVALID_HID);
    if (result != NULL) delete result; result = NULL;

    result = f->createReplicateHandles("01");
    BOOST_CHECK_GE(result->group, 0);
    BOOST_CHECK_GE(result->speciesCountsDataset, 0);
    BOOST_CHECK_NO_THROW(f->closeReplicateHandles(result));
    BOOST_CHECK_EQUAL(result->group, H5I_INVALID_HID);
    BOOST_CHECK_EQUAL(result->speciesCountsDataset, H5I_INVALID_HID);
    if (result != NULL) delete result; result = NULL;

    result = f->createReplicateHandles("02");
    BOOST_CHECK_GE(result->group, 0);
    BOOST_CHECK_GE(result->speciesCountsDataset, 0);
    BOOST_CHECK_NO_THROW(f->closeReplicateHandles(result));
    BOOST_CHECK_EQUAL(result->group, H5I_INVALID_HID);
    BOOST_CHECK_EQUAL(result->speciesCountsDataset, H5I_INVALID_HID);
    if (result != NULL) delete result; result = NULL;

    result = f->createReplicateHandles("10");
    BOOST_CHECK_GE(result->group, 0);
    BOOST_CHECK_GE(result->speciesCountsDataset, 0);
    BOOST_CHECK_NO_THROW(f->closeReplicateHandles(result));
    BOOST_CHECK_EQUAL(result->group, H5I_INVALID_HID);
    BOOST_CHECK_EQUAL(result->speciesCountsDataset, H5I_INVALID_HID);
    if (result != NULL) delete result; result = NULL;

    if (f != NULL) delete f; f = NULL;

    // Make sure the replciates appeared in the file.

    remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(OpenReplicateHandles)
{
    string filename = TMP_DIR+"/lmtest_SimulationFile_OpenReplicateHandles.h5";
    remove(filename.c_str());

    SimulationFileTester * f = NULL;

    // Create a new file.
    SimulationFile::create(filename, 10);

    // Open some new replicates, making sure they are created and added to the open list.
    SimulationFile::ReplicateHandles * result, * result1;
    f=new SimulationFileTester(filename);
    BOOST_CHECK_EQUAL(f->openReplicatesP->size(), 0U);
    result = f->openReplicateHandles(0);
    BOOST_CHECK_GE(result->group, 0);
    BOOST_CHECK_GE(result->speciesCountsDataset, 0);
    BOOST_CHECK_EQUAL(f->openReplicatesP->size(), 1U);

    result1 = f->openReplicateHandles(1);
    BOOST_CHECK_GE(result1->group, 0);
    BOOST_CHECK_GE(result1->speciesCountsDataset, 0);
    BOOST_CHECK_EQUAL(f->openReplicatesP->size(), 2U);

    result = f->openReplicateHandles(2);
    BOOST_CHECK_GE(result->group, 0);
    BOOST_CHECK_GE(result->speciesCountsDataset, 0);
    BOOST_CHECK_EQUAL(f->openReplicatesP->size(), 3U);

    result = f->openReplicateHandles(10);
    BOOST_CHECK_GE(result->group, 0);
    BOOST_CHECK_GE(result->speciesCountsDataset, 0);
    BOOST_CHECK_EQUAL(f->openReplicatesP->size(), 4U);

    // Make sure opening an already open group returns the same handle and keeps the open list the same.
    result = f->openReplicateHandles(1);
    BOOST_CHECK_EQUAL(result, result1);
    BOOST_CHECK_EQUAL(result->group, result1->group);
    BOOST_CHECK_EQUAL(result->speciesCountsDataset, result1->speciesCountsDataset);
    BOOST_CHECK_EQUAL(f->openReplicatesP->size(), 4U);

    // Make sure closing the file closes all of the groups.
    f->close();
    BOOST_CHECK_EQUAL(f->openReplicatesP->size(), 0U);

    // Destroy the file object.
    if (f != NULL) delete f; f = NULL;

    remove(filename.c_str());
}

/*BOOST_AUTO_TEST_CASE(AppendSpeciesCounts)
{
    string filename = TMP_DIR+"/lmtest_SimulationFile_AppendSpeciesCounts.h5";
    remove(filename.c_str());

    SimulationFileTester * f = NULL;

    // Create a new file.
    SimulationFile::create(filename, 10);

    // Create replicate and add some species counts.
    vector<unsigned int> counts;
    for (int i=0; i<10; i++)
        counts.push_back(i*2+5);
    BOOST_CHECK_NO_THROW(f=new SimulationFileTester(filename));
    BOOST_CHECK_NO_THROW(f->appendSpeciesCounts(0, 0.0, counts));

    if (f != NULL) delete f; f = NULL;

    // Make sure the replicates appeared in the file.

    remove(filename.c_str());
}*/

BOOST_AUTO_TEST_CASE(AppendParameterValues)
{
    string filename = TMP_DIR+"/lmtest_SimulationFile_AppendParameterValues.h5";
    remove(filename.c_str());

    SimulationFileTester * f = NULL;
    ParameterValues * v = NULL;

    // Create a new file.
    SimulationFile::create(filename, 10);

    // Create replicate and add some species counts.
    BOOST_CHECK_NO_THROW(v=new ParameterValues);
    BOOST_CHECK_NO_THROW(f=new SimulationFileTester(filename));
    v->set_parameter("test1");
    for (int i=0; i<10; i++)
    {
        v->add_time((double)i);
        v->add_value(100.0-(double)i);
    }
    BOOST_CHECK_NO_THROW(f->appendParameterValues(0,v));
    v->Clear();
    v->set_parameter("test2");
    for (int i=0; i<10; i++)
    {
        v->add_time((double)i);
        v->add_value(100.0-(double)i);
    }
    BOOST_CHECK_NO_THROW(f->appendParameterValues(0,v));
    v->Clear();
    v->set_parameter("test2");
    for (int i=10; i<50; i++)
    {
        v->add_time((double)i);
        v->add_value(100.0-(double)i);
    }
    BOOST_CHECK_NO_THROW(f->appendParameterValues(0,v));

    if (f != NULL) delete f; f = NULL;
    if (v != NULL) delete v; v = NULL;

    // Make sure the replicates appeared in the file.
    //remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(SetFirstPassageTimes)
{
    string filename = TMP_DIR+"/lmtest_SimulationFile_SetFirstPassageTimes.h5";
    remove(filename.c_str());

    SimulationFileTester * f = NULL;
    FirstPassageTimes * v = NULL;

    // Create a new file.
    SimulationFile::create(filename, 10);

    // Create replicate and add some species counts.
    BOOST_CHECK_NO_THROW(v=new FirstPassageTimes);
    BOOST_CHECK_NO_THROW(f=new SimulationFileTester(filename));

    // Test add an initial set of first passage times.
    v->Clear();
    v->set_species(0);
    v->add_species_count(100); v->add_first_passage_time(0.0);
    v->add_species_count(99); v->add_first_passage_time(1.0);
    v->add_species_count(98); v->add_first_passage_time(2.0);
    v->add_species_count(101); v->add_first_passage_time(3.0);
    v->add_species_count(102); v->add_first_passage_time(4.0);
    v->add_species_count(97); v->add_first_passage_time(5.0);
    v->add_species_count(103); v->add_first_passage_time(6.0);
    BOOST_CHECK_NO_THROW(f->setFirstPassageTimes(0,v));

    // Test appending to a set of first passage times.
    v->Clear();
    v->set_species(1);
    v->add_species_count(100); v->add_first_passage_time(0.0);
    v->add_species_count(99); v->add_first_passage_time(1.0);
    v->add_species_count(98); v->add_first_passage_time(2.0);
    v->add_species_count(101); v->add_first_passage_time(3.0);
    v->add_species_count(102); v->add_first_passage_time(4.0);
    v->add_species_count(97); v->add_first_passage_time(5.0);
    v->add_species_count(103); v->add_first_passage_time(6.0);
    BOOST_CHECK_NO_THROW(f->setFirstPassageTimes(0,v));
    v->Clear();
    v->set_species(1);
    v->add_species_count(104); v->add_first_passage_time(7.0);
    v->add_species_count(105); v->add_first_passage_time(8.0);
    v->add_species_count(106); v->add_first_passage_time(9.0);
    v->add_species_count(110); v->add_first_passage_time(13.0);
    v->add_species_count(109); v->add_first_passage_time(12.0);
    v->add_species_count(108); v->add_first_passage_time(11.0);
    v->add_species_count(107); v->add_first_passage_time(10.0);
    BOOST_CHECK_NO_THROW(f->setFirstPassageTimes(0,v));
    v->Clear();
    v->set_species(1);
    v->add_species_count(111); v->add_first_passage_time(14.0);
    v->add_species_count(112); v->add_first_passage_time(15.0);
    v->add_species_count(113); v->add_first_passage_time(16.0);
    v->add_species_count(114); v->add_first_passage_time(17.0);
    v->add_species_count(115); v->add_first_passage_time(18.0);
    v->add_species_count(117); v->add_first_passage_time(20.0);
    v->add_species_count(116); v->add_first_passage_time(19.0);
    BOOST_CHECK_NO_THROW(f->setFirstPassageTimes(0,v));

    // Test inserting before and after to a set of first passage times.
    v->Clear();
    v->set_species(2);
    v->add_species_count(100); v->add_first_passage_time(0.0);
    v->add_species_count(99); v->add_first_passage_time(1.0);
    v->add_species_count(98); v->add_first_passage_time(2.0);
    v->add_species_count(101); v->add_first_passage_time(3.0);
    v->add_species_count(102); v->add_first_passage_time(4.0);
    v->add_species_count(97); v->add_first_passage_time(5.0);
    v->add_species_count(103); v->add_first_passage_time(6.0);
    BOOST_CHECK_NO_THROW(f->setFirstPassageTimes(0,v));
    v->Clear();
    v->set_species(2);
    v->add_species_count(104); v->add_first_passage_time(7.0);
    v->add_species_count(105); v->add_first_passage_time(8.0);
    v->add_species_count(106); v->add_first_passage_time(9.0);
    v->add_species_count(96); v->add_first_passage_time(9.5);
    v->add_species_count(110); v->add_first_passage_time(13.0);
    v->add_species_count(95); v->add_first_passage_time(13.1);
    v->add_species_count(109); v->add_first_passage_time(12.0);
    v->add_species_count(108); v->add_first_passage_time(11.0);
    v->add_species_count(107); v->add_first_passage_time(10.0);
    v->add_species_count(93); v->add_first_passage_time(15.0);
    v->add_species_count(94); v->add_first_passage_time(14.0);
    BOOST_CHECK_NO_THROW(f->setFirstPassageTimes(0,v));
    v->Clear();
    v->set_species(2);
    v->add_species_count(111); v->add_first_passage_time(14.0);
    v->add_species_count(92); v->add_first_passage_time(16.0);
    v->add_species_count(112); v->add_first_passage_time(15.0);
    v->add_species_count(113); v->add_first_passage_time(16.0);
    v->add_species_count(114); v->add_first_passage_time(17.0);
    v->add_species_count(91); v->add_first_passage_time(17.0);
    v->add_species_count(115); v->add_first_passage_time(18.0);
    v->add_species_count(117); v->add_first_passage_time(20.0);
    v->add_species_count(90); v->add_first_passage_time(18.0);
    v->add_species_count(116); v->add_first_passage_time(19.0);
    BOOST_CHECK_NO_THROW(f->setFirstPassageTimes(0,v));

    if (f != NULL) delete f; f = NULL;
    if (v != NULL) delete v; v = NULL;

    // Make sure the replicates appeared in the file.
    //remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(GetReactionModel)
{
    SimulationFileTester * t = NULL;
    ReactionModel * m = NULL;

    // Check getting a generic reaction model with no reactionNumber property.
    {
    m = new ReactionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_GetReactionModel_NoReactionCount.h5"));
    BOOST_CHECK_NO_THROW(t->getReactionModel(m));
    BOOST_CHECK_EQUAL(m->number_species(), 5U);
    BOOST_CHECK_EQUAL(m->number_reactions(), 0U);
    BOOST_CHECK_EQUAL(m->initial_species_count_size(), 0);
    }

    // Check getting a generic reaction model.
    {
    m = new ReactionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_GetReactionModel.h5"));
    BOOST_CHECK_NO_THROW(t->getReactionModel(m));
    BOOST_CHECK_EQUAL(m->number_species(), 5U);
    BOOST_CHECK_EQUAL(m->number_reactions(), 7U);
    BOOST_CHECK_EQUAL(m->initial_species_count_size(), 5);
    BOOST_CHECK_EQUAL(m->initial_species_count(0), 11U);
    BOOST_CHECK_EQUAL(m->initial_species_count(1), 13U);
    BOOST_CHECK_EQUAL(m->initial_species_count(2), 17U);
    BOOST_CHECK_EQUAL(m->initial_species_count(3), 19U);
    BOOST_CHECK_EQUAL(m->initial_species_count(4), 23U);
    BOOST_CHECK_EQUAL(m->reaction_size(), 7);
    BOOST_CHECK_EQUAL(m->reaction(0).type(), 1U);
    BOOST_CHECK_EQUAL(m->reaction(1).type(), 1U);
    BOOST_CHECK_EQUAL(m->reaction(2).type(), 2U);
    BOOST_CHECK_EQUAL(m->reaction(3).type(), 2U);
    BOOST_CHECK_EQUAL(m->reaction(4).type(), 91U);
    BOOST_CHECK_EQUAL(m->reaction(5).type(), 98U);
    BOOST_CHECK_EQUAL(m->reaction(6).type(), 99U);
    BOOST_CHECK_EQUAL(m->reaction(0).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(1).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(2).rate_constant_size(), 2);
    BOOST_CHECK_EQUAL(m->reaction(3).rate_constant_size(), 2);
    BOOST_CHECK_EQUAL(m->reaction(4).rate_constant_size(), 10);
    BOOST_CHECK_EQUAL(m->reaction(5).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(6).rate_constant_size(), 1);
    BOOST_CHECK_CLOSE(m->reaction(0).rate_constant(0), 1.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(0), 7.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(0), 16.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(1), 17.3, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(0), 1678.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(1), 1.0e12, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(0), 2.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(1), 3.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(2), 4.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(3), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(4), 6.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(5), 7.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(6), 8.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(7), 9.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(8), 10.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(9), 11.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(5).rate_constant(0), 5.6, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(6).rate_constant(0), 6.7, 1e-9);
    BOOST_CHECK_EQUAL(m->stoichiometric_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    int S[] = {1, -1,  0, -1,  0,  1, -1,
            -1,  1,  0,  0,  1,  1, -1,
             0,  0,  1,  1,  1,  1, -1,
             0,  0, -1,  0,  0,  1, -1,
             0,  0,  0, -1,  0,  1, -1};
    for (int i=0; i<m->stoichiometric_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->stoichiometric_matrix(i), S[i]);
    }
    BOOST_CHECK_EQUAL(m->dependency_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    uint D[] = {1,  0,  0,  1,  0,  0,  1,
                0,  1,  0,  1,  0,  0,  1,
                0,  0,  1,  0,  1,  0,  1,
                0,  1,  0,  0,  1,  1,  1,
                1,  0,  0,  0,  0,  1,  0};
    for (int i=0; i<m->dependency_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->dependency_matrix(i), D[i]);
    }
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }

    // Check getting a reaction model with noise.
    {
    m = new ReactionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_GetReactionModel_Noise.h5"));
    BOOST_CHECK_NO_THROW(t->getReactionModel(m));
    BOOST_CHECK_EQUAL(m->number_species(), 5U);
    BOOST_CHECK_EQUAL(m->number_reactions(), 7U);
    BOOST_CHECK_EQUAL(m->initial_species_count_size(), 5);
    BOOST_CHECK_EQUAL(m->initial_species_count(0), 11U);
    BOOST_CHECK_EQUAL(m->initial_species_count(1), 13U);
    BOOST_CHECK_EQUAL(m->initial_species_count(2), 17U);
    BOOST_CHECK_EQUAL(m->initial_species_count(3), 19U);
    BOOST_CHECK_EQUAL(m->initial_species_count(4), 23U);
    BOOST_CHECK_EQUAL(m->reaction_size(), 7);
    BOOST_CHECK_EQUAL(m->reaction(0).type(), 1U);
    BOOST_CHECK_EQUAL(m->reaction(1).type(), 1U);
    BOOST_CHECK_EQUAL(m->reaction(2).type(), 2U);
    BOOST_CHECK_EQUAL(m->reaction(3).type(), 2U);
    BOOST_CHECK_EQUAL(m->reaction(4).type(), 91U);
    BOOST_CHECK_EQUAL(m->reaction(5).type(), 98U);
    BOOST_CHECK_EQUAL(m->reaction(6).type(), 99U);
    BOOST_CHECK_EQUAL(m->reaction(0).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(1).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(2).rate_constant_size(), 2);
    BOOST_CHECK_EQUAL(m->reaction(3).rate_constant_size(), 2);
    BOOST_CHECK_EQUAL(m->reaction(4).rate_constant_size(), 10);
    BOOST_CHECK_EQUAL(m->reaction(5).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(6).rate_constant_size(), 1);
    BOOST_CHECK_CLOSE(m->reaction(0).rate_constant(0), 1.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(0), 7.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(0), 16.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(1), 17.3, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(0), 1678.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(1), 1.0e12, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(0), 2.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(1), 3.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(2), 4.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(3), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(4), 6.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(5), 7.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(6), 8.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(7), 9.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(8), 10.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_constant(9), 11.0, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(5).rate_constant(0), 5.6, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(6).rate_constant(0), 6.7, 1e-9);
    BOOST_CHECK_EQUAL(m->reaction(0).rate_has_noise(), false);
    BOOST_CHECK_EQUAL(m->reaction(1).rate_has_noise(), true);
    BOOST_CHECK_EQUAL(m->reaction(2).rate_has_noise(), false);
    BOOST_CHECK_EQUAL(m->reaction(3).rate_has_noise(), false);
    BOOST_CHECK_EQUAL(m->reaction(4).rate_has_noise(), true);
    BOOST_CHECK_EQUAL(m->reaction(5).rate_has_noise(), true);
    BOOST_CHECK_EQUAL(m->reaction(6).rate_has_noise(), false);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_noise_variance(), 11.11, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_noise_variance(), 4.4, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(5).rate_noise_variance(), 5.555, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_noise_tau(), 1e5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(4).rate_noise_tau(), 1e14, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(5).rate_noise_tau(), 1e15, 1e-9);
    BOOST_CHECK_EQUAL(m->stoichiometric_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    int S[] = {1, -1,  0, -1,  0,  1, -1,
            -1,  1,  0,  0,  1,  1, -1,
             0,  0,  1,  1,  1,  1, -1,
             0,  0, -1,  0,  0,  1, -1,
             0,  0,  0, -1,  0,  1, -1};
    for (int i=0; i<m->stoichiometric_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->stoichiometric_matrix(i), S[i]);
    }
    BOOST_CHECK_EQUAL(m->dependency_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    uint D[] = {1,  0,  0,  1,  0,  0,  1,
                0,  1,  0,  1,  0,  0,  1,
                0,  0,  1,  0,  1,  0,  1,
                0,  1,  0,  0,  1,  1,  1,
                1,  0,  0,  0,  0,  1,  0};
    for (int i=0; i<m->dependency_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->dependency_matrix(i), D[i]);
    }
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }
}

BOOST_AUTO_TEST_CASE(SetReactionModel)
{
    SimulationFileTester * t = NULL;
    ReactionModel * m = NULL;

    // Write out the reaction model.
    {
    m = new ReactionModel();
    string filename = TMP_DIR+"/lmtest_SimulationFile_SetReactionModel.h5";
    remove(filename.c_str());
    SimulationFile::create(filename);
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(filename));
    m->set_number_species(11U);
    m->set_number_reactions(3U);
    for (uint i=0; i<m->number_species(); i++)
        m->add_initial_species_count(i+1);
    for (uint i=0; i<m->number_reactions(); i++)
        m->add_reaction();
    m->mutable_reaction(0)->set_type(17);
    m->mutable_reaction(1)->set_type(37);
    m->mutable_reaction(2)->set_type(47);
    m->mutable_reaction(0)->add_rate_constant(1.7);
    m->mutable_reaction(1)->add_rate_constant(11.7);
    m->mutable_reaction(1)->add_rate_constant(111.7);
    m->mutable_reaction(2)->add_rate_constant(21.7);
    m->mutable_reaction(2)->add_rate_constant(211.7);
    m->mutable_reaction(2)->add_rate_constant(2111.7);
    int S[] = {0,  1,  1,
               1,  0,  1,
               0,  1,  0,
               0,  1,  0,
              -1,  0,  1,
               0,  1,  1,
               0,  1,  0,
               0,  1,  0,
               1,  0,  1,
               1,  0,  1,
               1,  0, -1};
    for (uint i=0; i<m->number_species()*m->number_reactions(); i++) m->add_stoichiometric_matrix(S[i]);
    uint D[] = {0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1};
    for (uint i=0; i<m->number_species()*m->number_reactions(); i++) m->add_dependency_matrix(D[i]);
    BOOST_CHECK_NO_THROW(t->setReactionModel(m));
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }

    // Make sure the reaction model was written correctly.
    {
    m = new ReactionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(TMP_DIR+"/lmtest_SimulationFile_SetReactionModel.h5"));
    BOOST_CHECK_NO_THROW(t->getReactionModel(m));
    BOOST_CHECK_EQUAL(m->number_species(), 11U);
    BOOST_CHECK_EQUAL(m->number_reactions(), 3U);
    BOOST_CHECK_EQUAL(m->initial_species_count_size(), 11);
    BOOST_CHECK_EQUAL(m->initial_species_count(0), 1U);
    BOOST_CHECK_EQUAL(m->initial_species_count(1), 2U);
    BOOST_CHECK_EQUAL(m->initial_species_count(2), 3U);
    BOOST_CHECK_EQUAL(m->initial_species_count(3), 4U);
    BOOST_CHECK_EQUAL(m->initial_species_count(4), 5U);
    BOOST_CHECK_EQUAL(m->initial_species_count(5), 6U);
    BOOST_CHECK_EQUAL(m->initial_species_count(6), 7U);
    BOOST_CHECK_EQUAL(m->initial_species_count(7), 8U);
    BOOST_CHECK_EQUAL(m->initial_species_count(8), 9U);
    BOOST_CHECK_EQUAL(m->initial_species_count(9), 10U);
    BOOST_CHECK_EQUAL(m->initial_species_count(10), 11U);
    BOOST_CHECK_EQUAL(m->reaction_size(), 3);
    BOOST_CHECK_EQUAL(m->reaction(0).type(), 17U);
    BOOST_CHECK_EQUAL(m->reaction(1).type(), 37U);
    BOOST_CHECK_EQUAL(m->reaction(2).type(), 47U);
    BOOST_CHECK_EQUAL(m->reaction(0).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(1).rate_constant_size(), 2);
    BOOST_CHECK_EQUAL(m->reaction(2).rate_constant_size(), 3);
    BOOST_CHECK_CLOSE(m->reaction(0).rate_constant(0), 1.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(0), 11.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(1), 111.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(0), 21.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(1), 211.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(2), 2111.7, 1e-9);
    BOOST_CHECK_EQUAL(m->stoichiometric_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    int S[] = {0,  1,  1,
               1,  0,  1,
               0,  1,  0,
               0,  1,  0,
              -1,  0,  1,
               0,  1,  1,
               0,  1,  0,
               0,  1,  0,
               1,  0,  1,
               1,  0,  1,
               1,  0, -1};
    for (int i=0; i<m->stoichiometric_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->stoichiometric_matrix(i), S[i]);
    }
    BOOST_CHECK_EQUAL(m->dependency_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    uint D[] = {0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1};
    for (int i=0; i<m->dependency_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->dependency_matrix(i), D[i]);
    }
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }

    // Make sure we can overwrite the reaction model with a new model.
    {
    m = new ReactionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(TMP_DIR+"/lmtest_SimulationFile_SetReactionModel.h5"));
    m->set_number_species(12U);
    m->set_number_reactions(4U);
    for (uint i=0; i<m->number_species(); i++)
        m->add_initial_species_count(i+3);
    for (uint i=0; i<m->number_reactions(); i++)
        m->add_reaction();
    m->mutable_reaction(0)->set_type(15);
    m->mutable_reaction(1)->set_type(30);
    m->mutable_reaction(2)->set_type(45);
    m->mutable_reaction(3)->set_type(60);
    m->mutable_reaction(0)->add_rate_constant(1.5);
    m->mutable_reaction(1)->add_rate_constant(11.5);
    m->mutable_reaction(1)->add_rate_constant(111.5);
    m->mutable_reaction(2)->add_rate_constant(21.5);
    m->mutable_reaction(2)->add_rate_constant(211.5);
    m->mutable_reaction(2)->add_rate_constant(2111.5);
    m->mutable_reaction(3)->add_rate_constant(31.5);
    m->mutable_reaction(3)->add_rate_constant(311.5);
    m->mutable_reaction(3)->add_rate_constant(3111.5);
    m->mutable_reaction(3)->add_rate_constant(31111.5);
    int S[] = {0,  1,  1,  1,
               1,  0,  1,  1,
               0,  1,  0,  1,
               0,  1,  0,  1,
              -1,  0,  1,  1,
               0,  1,  1,  1,
               0,  1,  0,  1,
               0,  1,  0,  1,
               1,  0,  1,  1,
               1,  0,  1,  1,
               1,  0,  1,  1,
               1,  0, -1,  1};
    for (uint i=0; i<m->number_species()*m->number_reactions(); i++) m->add_stoichiometric_matrix(S[i]);
    uint D[] = {0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1};
    for (uint i=0; i<m->number_species()*m->number_reactions(); i++) m->add_dependency_matrix(D[i]);
    BOOST_CHECK_NO_THROW(t->setReactionModel(m));
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }

    // Make sure the reaction model was written correctly.
    {
    m = new ReactionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(TMP_DIR+"/lmtest_SimulationFile_SetReactionModel.h5"));
    BOOST_CHECK_NO_THROW(t->getReactionModel(m));
    BOOST_CHECK_EQUAL(m->number_species(), 12U);
    BOOST_CHECK_EQUAL(m->number_reactions(), 4U);
    BOOST_CHECK_EQUAL(m->initial_species_count_size(), 12);
    BOOST_CHECK_EQUAL(m->initial_species_count(0), 3U);
    BOOST_CHECK_EQUAL(m->initial_species_count(1), 4U);
    BOOST_CHECK_EQUAL(m->initial_species_count(2), 5U);
    BOOST_CHECK_EQUAL(m->initial_species_count(3), 6U);
    BOOST_CHECK_EQUAL(m->initial_species_count(4), 7U);
    BOOST_CHECK_EQUAL(m->initial_species_count(5), 8U);
    BOOST_CHECK_EQUAL(m->initial_species_count(6), 9U);
    BOOST_CHECK_EQUAL(m->initial_species_count(7), 10U);
    BOOST_CHECK_EQUAL(m->initial_species_count(8), 11U);
    BOOST_CHECK_EQUAL(m->initial_species_count(9), 12U);
    BOOST_CHECK_EQUAL(m->initial_species_count(10), 13U);
    BOOST_CHECK_EQUAL(m->initial_species_count(11), 14U);
    BOOST_CHECK_EQUAL(m->reaction_size(), 4);
    BOOST_CHECK_EQUAL(m->reaction(0).type(), 15U);
    BOOST_CHECK_EQUAL(m->reaction(1).type(), 30U);
    BOOST_CHECK_EQUAL(m->reaction(2).type(), 45U);
    BOOST_CHECK_EQUAL(m->reaction(3).type(), 60U);
    BOOST_CHECK_EQUAL(m->reaction(0).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(1).rate_constant_size(), 2);
    BOOST_CHECK_EQUAL(m->reaction(2).rate_constant_size(), 3);
    BOOST_CHECK_EQUAL(m->reaction(3).rate_constant_size(), 4);
    BOOST_CHECK_CLOSE(m->reaction(0).rate_constant(0), 1.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(0), 11.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(1), 111.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(0), 21.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(1), 211.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(2), 2111.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(0), 31.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(1), 311.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(2), 3111.5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(3).rate_constant(3), 31111.5, 1e-9);
    int S[] = {0,  1,  1,  1,
               1,  0,  1,  1,
               0,  1,  0,  1,
               0,  1,  0,  1,
              -1,  0,  1,  1,
               0,  1,  1,  1,
               0,  1,  0,  1,
               0,  1,  0,  1,
               1,  0,  1,  1,
               1,  0,  1,  1,
               1,  0,  1,  1,
               1,  0, -1,  1};
    for (int i=0; i<m->stoichiometric_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->stoichiometric_matrix(i), S[i]);
    }
    BOOST_CHECK_EQUAL(m->dependency_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    uint D[] = {0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1,
                1,  0,  1,  1,
                1,  0,  1,  1,
                0,  1,  1,  1};
    for (int i=0; i<m->dependency_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->dependency_matrix(i), D[i]);
    }
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }

    // Write out a reaction model with noisy rates.
    {
    m = new ReactionModel();
    string filename = TMP_DIR+"/lmtest_SimulationFile_SetReactionModel_Noise.h5";
    remove(filename.c_str());
    SimulationFile::create(filename);
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(filename));
    m->set_number_species(11U);
    m->set_number_reactions(3U);
    for (uint i=0; i<m->number_species(); i++)
        m->add_initial_species_count(i+1);
    for (uint i=0; i<m->number_reactions(); i++)
        m->add_reaction();
    m->mutable_reaction(0)->set_type(17);
    m->mutable_reaction(1)->set_type(37);
    m->mutable_reaction(2)->set_type(47);
    m->mutable_reaction(0)->add_rate_constant(1.7);
    m->mutable_reaction(1)->add_rate_constant(11.7);
    m->mutable_reaction(1)->add_rate_constant(111.7);
    m->mutable_reaction(2)->add_rate_constant(21.7);
    m->mutable_reaction(2)->add_rate_constant(211.7);
    m->mutable_reaction(2)->add_rate_constant(2111.7);
    m->mutable_reaction(0)->set_rate_has_noise(true);
    m->mutable_reaction(0)->set_rate_noise_variance(17.178);
    m->mutable_reaction(0)->set_rate_noise_tau(1e5);
    m->mutable_reaction(1)->set_rate_has_noise(true);
    m->mutable_reaction(1)->set_rate_noise_variance(123.456789);
    m->mutable_reaction(1)->set_rate_noise_tau(1e9);
    int S[] = {0,  1,  1,
               1,  0,  1,
               0,  1,  0,
               0,  1,  0,
              -1,  0,  1,
               0,  1,  1,
               0,  1,  0,
               0,  1,  0,
               1,  0,  1,
               1,  0,  1,
               1,  0, -1};
    for (uint i=0; i<m->number_species()*m->number_reactions(); i++) m->add_stoichiometric_matrix(S[i]);
    uint D[] = {0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1};
    for (uint i=0; i<m->number_species()*m->number_reactions(); i++) m->add_dependency_matrix(D[i]);
    BOOST_CHECK_NO_THROW(t->setReactionModel(m));
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }

    // Make sure the reaction model was written correctly.
    {
    m = new ReactionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(TMP_DIR+"/lmtest_SimulationFile_SetReactionModel_Noise.h5"));
    BOOST_CHECK_NO_THROW(t->getReactionModel(m));
    BOOST_CHECK_EQUAL(m->number_species(), 11U);
    BOOST_CHECK_EQUAL(m->number_reactions(), 3U);
    BOOST_CHECK_EQUAL(m->initial_species_count_size(), 11);
    BOOST_CHECK_EQUAL(m->initial_species_count(0), 1U);
    BOOST_CHECK_EQUAL(m->initial_species_count(1), 2U);
    BOOST_CHECK_EQUAL(m->initial_species_count(2), 3U);
    BOOST_CHECK_EQUAL(m->initial_species_count(3), 4U);
    BOOST_CHECK_EQUAL(m->initial_species_count(4), 5U);
    BOOST_CHECK_EQUAL(m->initial_species_count(5), 6U);
    BOOST_CHECK_EQUAL(m->initial_species_count(6), 7U);
    BOOST_CHECK_EQUAL(m->initial_species_count(7), 8U);
    BOOST_CHECK_EQUAL(m->initial_species_count(8), 9U);
    BOOST_CHECK_EQUAL(m->initial_species_count(9), 10U);
    BOOST_CHECK_EQUAL(m->initial_species_count(10), 11U);
    BOOST_CHECK_EQUAL(m->reaction_size(), 3);
    BOOST_CHECK_EQUAL(m->reaction(0).type(), 17U);
    BOOST_CHECK_EQUAL(m->reaction(1).type(), 37U);
    BOOST_CHECK_EQUAL(m->reaction(2).type(), 47U);
    BOOST_CHECK_EQUAL(m->reaction(0).rate_constant_size(), 1);
    BOOST_CHECK_EQUAL(m->reaction(1).rate_constant_size(), 2);
    BOOST_CHECK_EQUAL(m->reaction(2).rate_constant_size(), 3);
    BOOST_CHECK_CLOSE(m->reaction(0).rate_constant(0), 1.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(0), 11.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_constant(1), 111.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(0), 21.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(1), 211.7, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(2).rate_constant(2), 2111.7, 1e-9);
    BOOST_CHECK_EQUAL(m->reaction(0).rate_has_noise(), true);
    BOOST_CHECK_EQUAL(m->reaction(1).rate_has_noise(), true);
    BOOST_CHECK_EQUAL(m->reaction(2).rate_has_noise(), false);
    BOOST_CHECK_CLOSE(m->reaction(0).rate_noise_variance(), 17.178, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_noise_variance(), 123.456789, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(0).rate_noise_tau(), 1e5, 1e-9);
    BOOST_CHECK_CLOSE(m->reaction(1).rate_noise_tau(), 1e9, 1e-9);
    BOOST_CHECK_EQUAL(m->stoichiometric_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    int S[] = {0,  1,  1,
               1,  0,  1,
               0,  1,  0,
               0,  1,  0,
              -1,  0,  1,
               0,  1,  1,
               0,  1,  0,
               0,  1,  0,
               1,  0,  1,
               1,  0,  1,
               1,  0, -1};
    for (int i=0; i<m->stoichiometric_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->stoichiometric_matrix(i), S[i]);
    }
    BOOST_CHECK_EQUAL(m->dependency_matrix_size(), (int)(m->number_species()*m->number_reactions()));
    uint D[] = {0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1,
                1,  0,  1,
                0,  1,  1};
    for (int i=0; i<m->dependency_matrix_size(); i++)
    {
        BOOST_CHECK_EQUAL(m->dependency_matrix(i), D[i]);
    }
    if (m != NULL) delete m; m = NULL;
    if (t != NULL) delete t; t = NULL;
    }
}

BOOST_AUTO_TEST_CASE(GetDiffusionModel)
{
    SimulationFileTester * t = NULL;
    DiffusionModel * d = NULL;

    d = new DiffusionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_GetDiffusionModel.h5"));
    BOOST_CHECK_NO_THROW(t->getDiffusionModel(d));
    BOOST_CHECK_CLOSE(d->lattice_spacing(), 16.0, 1e-9);
    BOOST_CHECK_EQUAL(d->lattice_x_size(), 3U);
    BOOST_CHECK_EQUAL(d->lattice_y_size(), 4U);
    BOOST_CHECK_EQUAL(d->lattice_z_size(), 5U);
    BOOST_CHECK_EQUAL(d->particles_per_site(), 2U);
    if (t != NULL) delete t; t = NULL;
    if (d != NULL) delete d; d = NULL;
}

/**
 * Matlab script to set some model data.

    data=zeros(3,4,5,2);
    for p=[1:size(data,4)]
        for z=[1:size(data,3)]
            for y=[1:size(data,2)]
                for x=[1:size(data,1)]
                    data(x,y,z,p)=((p-1)*20)+((z-1)*50)+((y-1)*size(data,1))+(x-1);
                end
            end
        end
    end
    h5write('testing/data/lm/io/hdf5/SimulationFile_GetDiffusionModel.h5','/Model/Diffusion/Lattice',cast(data,'uint8'));
    data=zeros(3,4,5);
    for z=[1:size(data,3)]
        for y=[1:size(data,2)]
            for x=[1:size(data,1)]
                data(x,y,z)=((z-1)*50)+((y-1)*size(data,1))+(x-1)+1;
            end
        end
    end
    h5write('testing/data/lm/io/hdf5/SimulationFile_GetDiffusionModel.h5','/Model/Diffusion/LatticeSites',cast(data,'uint8'));
 */
BOOST_AUTO_TEST_CASE(GetDiffusionModelLattice)
{
    SimulationFileTester * t = NULL;
    DiffusionModel * d = NULL;
    byte * l1 = NULL;
    uint l1Size;
    byte * l2 = NULL;
    uint l2Size;

    d = new DiffusionModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_GetDiffusionModel.h5"));
    BOOST_CHECK_NO_THROW(t->getDiffusionModel(d));
    BOOST_CHECK_CLOSE(d->lattice_spacing(), 16.0, 1e-9);
    BOOST_CHECK_EQUAL(d->lattice_x_size(), 3U);
    BOOST_CHECK_EQUAL(d->lattice_y_size(), 4U);
    BOOST_CHECK_EQUAL(d->lattice_z_size(), 5U);
    BOOST_CHECK_EQUAL(d->particles_per_site(), 2U);
    l1Size = d->lattice_x_size()*d->lattice_y_size()*d->lattice_z_size()*d->particles_per_site();
    l1 = new byte[l1Size];
    memset(l1,0xFF,l1Size);
    l2Size = d->lattice_x_size()*d->lattice_y_size()*d->lattice_z_size();
    l2 = new byte[l2Size];
    memset(l2,0xFF,l2Size);
    t->getDiffusionModelLattice(d, l1, l1Size, l2, l2Size);
    for (uint p=0; p<d->particles_per_site(); p++)
    {
        for (uint z=0; z<d->lattice_z_size(); z++)
        {
            for (uint y=0; y<d->lattice_y_size(); y++)
            {
                for (uint x=0; x<d->lattice_x_size(); x++)
                {
                    BOOST_CHECK_EQUAL(l1[x+(y*d->lattice_x_size())+(z*d->lattice_x_size()*d->lattice_y_size())+(p*d->lattice_x_size()*d->lattice_y_size()*d->lattice_z_size())], (p*20)+(z*50)+(y*d->lattice_x_size())+x);
                }
            }
        }
    }
    for (uint z=0; z<d->lattice_z_size(); z++)
    {
        for (uint y=0; y<d->lattice_y_size(); y++)
        {
            for (uint x=0; x<d->lattice_x_size(); x++)
            {
                BOOST_CHECK_EQUAL(l2[x+(y*d->lattice_x_size())+(z*d->lattice_x_size()*d->lattice_y_size())], (z*50)+(y*d->lattice_x_size())+x+1);
            }
        }
    }
    if (t != NULL) delete t; t = NULL;
    if (d != NULL) delete d; d = NULL;
    if (l1 != NULL) delete[] l1; l1 = NULL;
    if (l2 != NULL) delete[] l2; l2 = NULL;
}

/**
 * Matlab script to set some model data.

    data=zeros(3,4,5,2);
    for p=[1:size(data,4)]
        for z=[1:size(data,3)]
            for y=[1:size(data,2)]
                for x=[1:size(data,1)]
                    data(x,y,z,p)=((p-1)*20)+((z-1)*50)+((y-1)*size(data,1))+(x-1);
                end
            end
        end
    end
    h5write('testing/data/lm/io/hdf5/SimulationFile_GetDiffusionModel.h5','/Model/Diffusion/Lattice',cast(data,'uint8'));
    data=zeros(3,4,5);
    for z=[1:size(data,3)]
        for y=[1:size(data,2)]
            for x=[1:size(data,1)]
                data(x,y,z)=((z-1)*50)+((y-1)*size(data,1))+(x-1)+1;
            end
        end
    end
    h5write('testing/data/lm/io/hdf5/SimulationFile_GetDiffusionModel.h5','/Model/Diffusion/LatticeSites',cast(data,'uint8'));
 */
BOOST_AUTO_TEST_CASE(GetSpatialModelObjects)
{
    SimulationFileTester * t = NULL;
    SpatialModel * m = NULL;

    m = new SpatialModel();
    BOOST_CHECK_NO_THROW(t=new SimulationFileTester(DATA_DIR+"/lm/io/hdf5/SimulationFile_GetSpatialModelObjects.h5"));
    BOOST_CHECK_NO_THROW(t->getSpatialModel(m));
    BOOST_CHECK_EQUAL(m->sphere_xc_size(), 382);
    BOOST_CHECK_CLOSE(m->sphere_xc(0), 239.6, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_yc(0), 155.9, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_zc(0), 143.6, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_radius(0), 10.4, 1e-9);
    BOOST_CHECK_CLOSE(m->sphere_type(0), 21.0, 1e-9);
    BOOST_CHECK_CLOSE(m->sphere_xc(100), 333.0, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_yc(100), 279.8, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_zc(100), 719.6, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_radius(100), 4.0, 1e-9);
    BOOST_CHECK_CLOSE(m->sphere_type(100), 25.0, 1e-9);
    BOOST_CHECK_CLOSE(m->sphere_xc(381), 76.7, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_yc(381), 191.7, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_zc(381), 949.0, 1e-1);
    BOOST_CHECK_CLOSE(m->sphere_radius(381), 1.7, 1e-9);
    BOOST_CHECK_CLOSE(m->sphere_type(381), 32.0, 1e-9);

    if (t != NULL) delete t; t = NULL;
    if (m != NULL) delete m; m = NULL;
}

BOOST_AUTO_TEST_SUITE_END()
