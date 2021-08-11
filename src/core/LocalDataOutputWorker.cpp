/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
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

#include <queue>
#include "config.h"

#if defined(MACOSX)
#include <mach/mach_time.h>
typedef uint64_t                    timing_time_t;
#define TIMING_GET_TIME ((timing_time_t)(((double)(mach_absolute_time()))*timing_time_mult))
#define TIMING_INIT \
    mach_timebase_info_data_t timing_info;\
    mach_timebase_info(&timing_info);\
    double timing_time_mult=((double)timing_info.numer)/((double)timing_info.denom);
#elif defined(LINUX)
#include <time.h>
typedef unsigned long long          timing_time_t;
#define TIMING_GET_TIME ((clock_gettime(CLOCK_MONOTONIC, &timing_timespec)==0)?((timing_time_t)((((timing_time_t)timing_timespec.tv_sec)*1000000000ULL)+timing_timespec.tv_nsec)):((timing_time_t)0ULL))
#define TIMING_INIT \
        struct timespec timing_timespec;
#endif


#include "core/Print.h"
#include "core/DataOutputQueue.h"
#include "core/LocalDataOutputWorker.h"
#include "core/Globals.h"
#include "io/ArbitraryH5.h"
#include "FirstPassageTimes.pb.h"
#include "Lattice.pb.h"
#include "ParameterValues.pb.h"
#include "SpeciesCounts.pb.h"
#include "io/lm_hdf5.h"
#include "io/SimulationFile.h"
#include "thread/Thread.h"
#include "thread/Worker.h"
#include "rdme/ByteLattice.h"
#include "rdme/IntLattice.h"
#include "lptf/Profile.h"

using lm::Print;
using lm::io::hdf5::SimulationFile;
using lm::thread::PthreadException;

namespace lm {
namespace main {


LocalDataOutputWorker::LocalDataOutputWorker(SimulationFile * file)
:file(file),shouldCheckpoint(false)
{
    PTHREAD_EXCEPTION_CHECK(pthread_cond_init(&dataAvailable, NULL));
}

LocalDataOutputWorker::~LocalDataOutputWorker()
{
    PTHREAD_EXCEPTION_CHECK_NOTHROW(pthread_cond_destroy(&dataAvailable));
}

void LocalDataOutputWorker::pushDataSet(DataSet * dataSet)
{
    bool success=false;

    //// BEGIN CRITICAL SECTION: dataMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&dataMutex));
    if (running)
    {
        if (dataSet != NULL)
        {
            dataQueue.push(dataSet);
            PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&dataAvailable));
        }
        success = true;
    }
    else
    {
        delete dataSet;
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&dataMutex));
    //// END CRITICAL SECTION: dataMutex

    if (!success) throw lm::Exception("LocalDataOutputWorker is not running.");
}

void LocalDataOutputWorker::wake()
{
    //// BEGIN CRITICAL SECTION: dataMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&dataMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&dataAvailable));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&dataMutex));
    //// END CRITICAL SECTION: dataMutex
}

void LocalDataOutputWorker::abort()
{
    //// BEGIN CRITICAL SECTION: controlMutex,dataMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&dataMutex));
    if (running)
    {
        aborted = true;
        PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&dataAvailable));
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&dataMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex,dataMutex
}

void LocalDataOutputWorker::checkpoint()
{
    bool success=false;

    //// BEGIN CRITICAL SECTION: controlMutex,dataMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&dataMutex));
    if (running)
    {
        shouldCheckpoint = true;
        PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&dataAvailable));
        success = true;
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&dataMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex,dataMutex

    if (!success) throw lm::Exception("LocalDataOutputWorker is not running.");
}

int LocalDataOutputWorker::run()
{
    PROF_SET_THREAD(1);
    PROF_BEGIN(PROF_DATAOUTPUT_RUN);
    try
    {
        Print::printf(Print::INFO, "Data output thread running.");

        // Declare the objects used for deserializing the data sets.
        lm::io::FirstPassageTimes firstPassageTimes;
        lm::io::SpeciesCounts speciesCounts;
        lm::io::ParameterValues parameterValues;
        lm::io::Lattice lattice;

        // Start the info update clock.
        TIMING_INIT;
        timing_time_t lastUpdateTime = TIMING_GET_TIME;
        timing_time_t writingTime = 0;
        uint bytesWritten = 0;
        uint datasetsWritten = 0;
        timing_time_t startWriting;
        uint datasetsRemaining = 0;

        // Write out data until we are stopped.
        bool doAbort = false;
        while (true)
        {
            bool finished = false;
            bool doCheckpoint = false;
            DataSet * dataSet = NULL;

            //// BEGIN CRITICAL SECTION: controlMutex,dataMutex
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&dataMutex));
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));

            Print::printf(Print::VERBOSE_DEBUG, "Data output thread looping: %d data sets to write.", dataQueue.size());

            // If we need to abort, do that with the highest priority.
            if (aborted)
            {
                doAbort = true;
            }

            // If we are not running and no data is left, we are done.
            else if (!running && dataQueue.empty())
            {
                finished=true;
            }

            // If we need to write a checkpoint, do so.
            else if (shouldCheckpoint)
            {
                shouldCheckpoint = false;
                doCheckpoint = true;
            }

            // If there is data in the queue get the next entry.
            else if (!dataQueue.empty())
            {
                dataSet = dataQueue.front();
                dataQueue.pop();
                datasetsRemaining = dataQueue.size();
            }

            // Otherwise, wait for a signal that more data is available.
            else
            {
                PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
                PTHREAD_EXCEPTION_CHECK(pthread_cond_wait(&dataAvailable, &dataMutex));
                PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex)); //TODO: This is not quite deadlock-safe, refactor.
            }

            PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&dataMutex));
            //// END CRITICAL SECTION: controlMutex,dataMutex

            // If we need to abort, close the file and then exit.
            if (doAbort)
            {
                Print::printf(Print::WARNING, "Received abort signal, attempting to save data.");
                file->close();
                Print::printf(Print::WARNING, "Data saved.");
                break;
            }

            // If we need to write out a checkpoint, do so.
            if (doCheckpoint)
            {
                time_t t1 = time(NULL);
                file->checkpoint();
                time_t t2 = time(NULL);
                Print::printf(Print::INFO, "Created checkpoint file in %d seconds.", (int)(t2-t1));
            }

            // If we got a data set to output, save it.
            else if (dataSet != NULL)
            {
                PROF_BEGIN(PROF_DATAOUTPUT_WRITE_DATASET);
                // Make sure the message is valid.
                if (dataSet->size < sizeof(uint)+sizeof(uint)+sizeof(size_t)+0+sizeof(size_t)+0) throw Exception("Invalid data set, size was too small for headers", dataSet->size);
                uint type = *(uint *)&dataSet->data[0];
                uint replicate = *(uint *)&dataSet->data[sizeof(uint)];
                size_t messageSize = *(size_t *)&dataSet->data[sizeof(uint)+sizeof(uint)];
                if (messageSize > dataSet->size-(sizeof(uint)+sizeof(uint)+sizeof(size_t)+sizeof(size_t))) throw Exception("Invalid data set, message size overflows the total size", messageSize);
                byte * messageData = &dataSet->data[sizeof(uint)+sizeof(uint)+sizeof(size_t)];
                size_t payloadSize = *(size_t *)&dataSet->data[sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize];
                if (payloadSize > dataSet->size-(sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize+sizeof(size_t))) throw Exception("Invalid data set, payload size overflows the total size", payloadSize);
                byte * payloadData = &dataSet->data[sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize+sizeof(size_t)];
                if (dataSet->size != sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize+sizeof(size_t)+payloadSize) throw Exception("Invalid data set, sum of sizes was not equal to the total size", sizeof(uint)+sizeof(uint)+sizeof(size_t)+messageSize+sizeof(size_t)+payloadSize);

                Print::printf(Print::VERBOSE_DEBUG, "Saving data set from replicate %d of type %d (total %d, message %d, payload %zu)", replicate, type, dataSet->size, messageSize, payloadSize);
                datasetsWritten++;
                bytesWritten += dataSet->size;

                switch (type)
                {
                case SPECIES_COUNTS:
                    PROF_BEGIN(PROF_DESERIALIZE_COUNTS);
                    speciesCounts.Clear();
                    speciesCounts.ParseFromArray(messageData, messageSize);
                    PROF_END(PROF_DESERIALIZE_COUNTS);
                    Print::printf(Print::VERBOSE_DEBUG, "Read species count %d entries for %d species with first time %f",speciesCounts.number_entries(), speciesCounts.number_species(), speciesCounts.time(0));
                    PROF_BEGIN(PROF_DATAOUTPUT_HDF_WRITE_COUNTS);
                    startWriting = TIMING_GET_TIME;
                    file->appendSpeciesCounts(replicate, &speciesCounts);
                    writingTime += TIMING_GET_TIME-startWriting;
                    PROF_END(PROF_DATAOUTPUT_HDF_WRITE_COUNTS);
                    break;

                case FIRST_PASSAGE_TIMES:
                    PROF_BEGIN(PROF_DESERIALIZE_FPT);
                    firstPassageTimes.Clear();
                    firstPassageTimes.ParseFromArray(messageData, messageSize);
                    PROF_END(PROF_DESERIALIZE_FPT);
                    Print::printf(Print::VERBOSE_DEBUG, "Read %d fpt(s) for species %d with first count %d and time %f",firstPassageTimes.first_passage_time_size(), firstPassageTimes.species(), firstPassageTimes.species_count(0), firstPassageTimes.first_passage_time(0));
                    PROF_BEGIN(PROF_DATAOUTPUT_HDF_WRITE_FPT);
                    startWriting = TIMING_GET_TIME;
                    file->setFirstPassageTimes(replicate, &firstPassageTimes);
                    writingTime += TIMING_GET_TIME-startWriting;
                    PROF_END(PROF_DATAOUTPUT_HDF_WRITE_FPT);
                    break;

                case PARAMETER_VALUES:
                    PROF_BEGIN(PROF_DESERIALIZE_PV);
                    parameterValues.Clear();
                    parameterValues.ParseFromArray(messageData, messageSize);
                    PROF_END(PROF_DESERIALIZE_PV);
                    PROF_BEGIN(PROF_DATAOUTPUT_HDF_WRITE_PV);
                    startWriting = TIMING_GET_TIME;
                    file->appendParameterValues(replicate, &parameterValues);
                    writingTime += TIMING_GET_TIME-startWriting;
                    PROF_END(PROF_DATAOUTPUT_HDF_WRITE_PV);
                    break;

                case BYTE_LATTICE:
				{
                    PROF_BEGIN(PROF_DESERIALIZE_LATTICE);
                    lattice.Clear();
                    lattice.ParseFromArray(messageData, messageSize);
                    uint8_t *latticeData = new uint8_t[payloadSize];
                    lm::rdme::ByteLattice::copyNativeToRowMajorByte(latticeData, payloadData, lattice.lattice_x_size(),lattice.lattice_y_size(),lattice.lattice_z_size(),lattice.particles_per_site(), payloadSize);
                    PROF_END(PROF_DESERIALIZE_LATTICE);
                    Print::printf(Print::VERBOSE_DEBUG, "Read byte lattice with time %e and size %zu (%d,%d,%d,%d)",lattice.time(), payloadSize,lattice.lattice_x_size(),lattice.lattice_y_size(),lattice.lattice_z_size(),lattice.particles_per_site());
                    PROF_BEGIN(PROF_DATAOUTPUT_HDF_WRITE_LATTICE);
                    startWriting = TIMING_GET_TIME;
                    file->appendLattice(replicate, &lattice, latticeData, payloadSize);
                    writingTime += TIMING_GET_TIME-startWriting;
                    delete [] latticeData;
                    latticeData = NULL;
                    PROF_END(PROF_DATAOUTPUT_HDF_WRITE_LATTICE);
                }   break;
                case INT_LATTICE:
				{
                    PROF_BEGIN(PROF_DESERIALIZE_LATTICE);
                    lattice.Clear();
                    lattice.ParseFromArray(messageData, messageSize);
                    uint32_t *latticeData = new uint32_t[payloadSize];
                    lm::rdme::IntLattice::copyNativeToRowMajor(latticeData, payloadData, lattice.lattice_x_size(),lattice.lattice_y_size(),lattice.lattice_z_size(),lattice.particles_per_site(), payloadSize);
                    PROF_END(PROF_DESERIALIZE_LATTICE);
                    Print::printf(Print::VERBOSE_DEBUG, "Read int lattice with time %e and size %zu (%d,%d,%d,%d)",lattice.time(), payloadSize,lattice.lattice_x_size(),lattice.lattice_y_size(),lattice.lattice_z_size(),lattice.particles_per_site());
                    PROF_BEGIN(PROF_DATAOUTPUT_HDF_WRITE_LATTICE);
                    startWriting = TIMING_GET_TIME;
                    file->appendLattice_U32LE(replicate, &lattice, latticeData, payloadSize);
                    writingTime += TIMING_GET_TIME-startWriting;
                    delete [] latticeData;
                    latticeData = NULL;
                    PROF_END(PROF_DATAOUTPUT_HDF_WRITE_LATTICE);
                }   break;
                case SITE_LATTICE:
				{
                    PROF_BEGIN(PROF_DESERIALIZE_LATTICE);
                    lattice.Clear();
                    lattice.ParseFromArray(messageData, messageSize);
                    uint8_t *siteData = new uint8_t[payloadSize];
                    lm::rdme::ByteLattice::copySitesNativeToRowMajorByte(siteData, payloadData, lattice.lattice_x_size(),lattice.lattice_y_size(),lattice.lattice_z_size(), payloadSize);
                    PROF_END(PROF_DESERIALIZE_LATTICE);
                    Print::printf(Print::VERBOSE_DEBUG, "Read site lattice with time %e and size %zu (%d,%d,%d,%d)",lattice.time(), payloadSize,lattice.lattice_x_size(),lattice.lattice_y_size(),lattice.lattice_z_size(),lattice.particles_per_site());
                    PROF_BEGIN(PROF_DATAOUTPUT_HDF_WRITE_LATTICE);
                    startWriting = TIMING_GET_TIME;
                    file->appendSites(replicate, &lattice, siteData, payloadSize);
                    writingTime += TIMING_GET_TIME-startWriting;
                    delete [] siteData;
                    siteData = NULL;
                    PROF_END(PROF_DATAOUTPUT_HDF_WRITE_LATTICE);
                }   break;

                case ARBITRARY_H5:
                    file->arbitraryH5(dataSet->data);
                    break;

                case ARBITRARY_H5_READ: {
                    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&returnMutex));
                    DataSet *reqData = new DataSet(); // XXX Need to free this
                    file->arbitraryH5Lookup(dataSet->data, reqData->data, reqData->size);
                    returnList.push_front(reqData);
                    pthread_cond_broadcast(&returnAvailable);
                    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&returnMutex));
                }   break;

                default:
                    Print::printf(Print::ERROR, "Received unknown data set type: %d",type);
                }

                // We are responsible for freeing the data set.
                delete dataSet;
                dataSet = NULL;
                PROF_END(PROF_DATAOUTPUT_WRITE_DATASET);
            }

            // Otherwise if we are finished, stop the main loop.
            else if (finished)
            {
                break;
            }

            // See if we should display some stats.
            timing_time_t currentTime = TIMING_GET_TIME;
            if (currentTime-lastUpdateTime > 60*1000000000ULL)
            {
                Print::printf(Print::INFO, "Wrote %u data sets (%u bytes) in the last %0.2f seconds (%0.2f seconds writing). %u datasets queued. Flushing.",datasetsWritten,bytesWritten,((double)(currentTime-lastUpdateTime))/1000000000.0, ((double)writingTime)/1000000000.0, datasetsRemaining);
                file->flush();
                lastUpdateTime = currentTime;
                writingTime = 0;
                datasetsWritten = 0;
                bytesWritten = 0;
            }
        }
        Print::printf(Print::INFO, "Data output thread finished.");
        PROF_END(PROF_DATAOUTPUT_RUN);
        return 0;
    }
    catch (lm::io::hdf5::HDF5Exception e)
    {
        Print::printf(Print::FATAL, "HDF5 exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (lm::thread::PthreadException e)
    {
        Print::printf(Print::FATAL, "Pthread exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (lm::Exception e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (std::exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown Exception during execution (%s:%d)", __FILE__, __LINE__);
    }

    PROF_END(PROF_DATAOUTPUT_RUN);
    return -1;
}

}
}
