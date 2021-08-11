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

#include <iostream>
#include <string>
#include <map>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <ctime>
#include "config.h"
#if defined(MACOSX)
#include <sys/time.h>
#endif
#include <csignal>
#include <cerrno>
#include <unistd.h>
#include <sys/wait.h>
#include <pthread.h>
#include <google/protobuf/stubs/common.h>
#include "core/Print.h"
#include "core/Exceptions.h"
#include "core/Types.h"
#include "core/Math.h"
#include "mpi/lm_mpi.h"
#ifdef OPT_CUDA
#include "cuda/lm_cuda.h"
#endif
#include "io/lm_hdf5.h"
#include "io/SimulationFile.h"
#include "DiffusionModel.pb.h"
#include "ReactionModel.pb.h"
#include "io/SimulationParameters.h"
#include "core/CheckpointSignaler.h"
#include "core/DataOutputQueue.h"
#include "core/LocalDataOutputWorker.h"
#include "core/MPIRemoteDataOutputQueue.h"
#include "cmd/common.h"
#include "core/Globals.h"
#include "core/SignalHandler.h"
#include "core/ReplicateRunner.h"
#include "core/ResourceAllocator.h"
#include "SimulationParameters.pb.h"
#include "thread/Thread.h"
#include "thread/WorkerManager.h"
#include "lptf/Profile.h"

using std::map;
using std::list;
using lm::Print;
using lm::Exception;
using lm::main::ReplicateRunner;
using lm::main::ResourceAllocator;
using lm::me::MESolverFactory;
using lm::thread::PthreadException;

void listDevicesMPI();
void executeSimulationMPI();
void executeSimulationMPISingleMaster();
void executeSimulationMPISingleSlave();
void broadcastSimulationParameters(void * staticDataBuffer, map<string,string> & simulationParameters);
void broadcastReactionModel(void * staticDataBuffer, lm::io::ReactionModel * reactionModel);
void broadcastDiffusionModel(void * staticDataBuffer, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize);
map<string,string> receiveSimulationParameters(void * staticDataBuffer);
void receiveReactionModel(void * staticDataBuffer, lm::io::ReactionModel * reactionModel);
void receiveDiffusionModel(void * staticDataBuffer, lm::io::DiffusionModel * diffusionModel, uint8_t ** lattice, size_t * latticeSize, uint8_t ** latticeSites, size_t * latticeSitesSize);
ReplicateRunner * startReplicate(int replicate, MESolverFactory solverFactory, std::map<std::string,string> & simulationParameters, lm::io::ReactionModel * reactionModel, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize, ResourceAllocator & resourceAllocator);
ReplicateRunner * popNextFinishedReplicate(list<ReplicateRunner *> & runningReplicates, ResourceAllocator & resourceAllocator);

// Allocate the profile space.
PROF_ALLOC;

int main(int argc, char** argv)
{
    // Make sure we are using the correct protocol buffers library.
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    PROF_INIT;
    try
    {
        // Initialize the MPI library.
        lm::MPI::init(argc, argv);

        // If this is the master process, catch any problems with the command line arguments.
        int validCommandLine=0;
        if (lm::MPI::worldRank == lm::MPI::MASTER)
        {
            //Print the startup messages.
            printCopyright(argc, argv);
            lm::MPI::printCapabilities();

            // Parse the command line arguments.
            try
            {
                parseArguments(argc, argv, "lm::rdme::MpdRdmeSolver");

                if (functionOption == "help")
                {
                    // Handle help on the master process.
                    printUsage(argc, argv);
                }
                else if (functionOption == "version")
                {
                    printBuildConfig();
                }
                else if (functionOption == "devices")
                {
                    validCommandLine = 1;
                }
                else if (functionOption == "simulation")
                {
                    validCommandLine = 1;
                }
                else
                {
                    throw lm::CommandLineArgumentException("unknown function.");
                }
            }
            catch (lm::CommandLineArgumentException e)
            {
                std::cerr << "Invalid command line argument: " << e.what() << std::endl << std::endl;
                printUsage(argc, argv);
            }
        }

        // Broadcast whether the command line was valid.
        MPI_EXCEPTION_CHECK(MPI_Bcast(&validCommandLine,1,MPI_INT,lm::MPI::MASTER,MPI_COMM_WORLD));

        // If we are sure the arguments are good.
        if (validCommandLine)
        {
            // Parse them again on all processes.
            parseArguments(argc, argv, "lm::rdme::MpdRdmeSolver");

            // Perform the requested function.
            if (functionOption == "devices")
            {
                listDevicesMPI();
            }
            else if (functionOption == "simulation")
            {
                executeSimulationMPI();
            }
        }

        // Close the MPI library.
        Print::printf(Print::INFO, "Closing MPI library.");
        lm::MPI::finalize();

        Print::printf(Print::INFO, "Program execution finished.");
        PROF_WRITE;
        google::protobuf::ShutdownProtobufLibrary();
        return 0;
    }
    catch (lm::io::hdf5::HDF5Exception & e)
    {
        std::cerr << "HDF5 exception during execution: " << e.what() << std::endl;
    }
    catch (lm::MPIException & e)
    {
        std::cerr << "MPI exception during execution: " << e.what() << std::endl;
    }
    catch (lm::thread::PthreadException & e)
    {
        std::cerr << "PthreadException exception during execution: " << e.what() << std::endl;
    }
    catch (lm::Exception & e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (std::exception & e)
    {
        std::cerr << "Std exception during execution: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Unknown Exception during execution." << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD,-1);
    lm::MPI::finalize();
    PROF_WRITE;
    google::protobuf::ShutdownProtobufLibrary();
    return -1;
}

void listDevicesMPI()
{
    // Get the MPI hostname.
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int hostnameLength;
    MPI_EXCEPTION_CHECK(MPI_Get_processor_name(hostname, &hostnameLength));

    // Print the capabilities message.
    printf("Process %d running on host %s with %d/%d processor(s)", lm::MPI::worldRank, hostname, numberCpuCores, getPhysicalCpuCores());

    #ifdef OPT_CUDA
    printf(" and %d/%d CUDA device(s)", (int)cudaDevices.size(), lm::CUDA::getNumberDevices());
    #endif
    printf(".\n");

    #ifdef OPT_CUDA
    if (shouldPrintCudaCapabilities)
    {
        for (int i=0; i<(int)cudaDevices.size(); i++)
        {
            printf("  %d-%s\n", lm::MPI::worldRank, lm::CUDA::getCapabilitiesString(cudaDevices[i]).c_str());
        }
    }
    #endif
}

void executeSimulationMPI()
{
    PROF_SET_THREAD(0);
    PROF_BEGIN(PROF_SIM_RUN);

    if (lm::MPI::worldRank == lm::MPI::MASTER)
        executeSimulationMPISingleMaster();
    else
        executeSimulationMPISingleSlave();

    PROF_END(PROF_SIM_RUN);
}

void executeSimulationMPISingleMaster()
{
    Print::printf(Print::DEBUG, "MPI master process %d started.", lm::MPI::worldRank);

    // Create the resource allocator, subtract one core for the data output thread on the master.
    #ifdef OPT_CUDA
    Print::printf(Print::INFO, "Using %d processor(s) and %d CUDA device(s) per process.", numberCpuCores, (int)cudaDevices.size());
    Print::printf(Print::INFO, "Assigning %0.2f processor(s) and %0.2f CUDA device(s) per replicate.", cpuCoresPerReplicate, cudaDevicesPerReplicate);
    ResourceAllocator resourceAllocator(lm::MPI::worldRank, numberCpuCores, cpuCoresPerReplicate, cudaDevices, cudaDevicesPerReplicate);
    #else
    Print::printf(Print::INFO, "Using %d processor(s) per process.", numberCpuCores);
    Print::printf(Print::INFO, "Assigning %0.2f processor(s) per replicate.", cpuCoresPerReplicate);
    ResourceAllocator resourceAllocator(lm::MPI::worldRank, numberCpuCores, cpuCoresPerReplicate);
    #endif

    // Reserve a core for the data output thread, unless we have a flag telling us not to.
    int reservedCpuCore = 0;
    if (shouldReserveOutputCore)
    {
    	reservedCpuCore=resourceAllocator.reserveCpuCore();
    	Print::printf(Print::INFO, "Reserved CPU core %d on process %d for data output.", reservedCpuCore, lm::MPI::worldRank);
    }

    // Create a worker to handle any signals.
    lm::main::SignalHandler * signalHandler = new lm::main::SignalHandler();
    signalHandler->setAffinity(reservedCpuCore);
    signalHandler->start();

    // Create the checkpoint signaler.
    lm::main::CheckpointSignaler * checkpointSignaler = new lm::main::CheckpointSignaler();
    checkpointSignaler->setAffinity(reservedCpuCore);
    checkpointSignaler->start();
    checkpointSignaler->startCheckpointing(checkpointInterval);

    // Open the file.
    lm::io::hdf5::SimulationFile * file = new lm::io::hdf5::SimulationFile(simulationFilename);

    // Start the data output thread.
    lm::main::LocalDataOutputWorker * dataOutputWorker = new lm::main::LocalDataOutputWorker(file);
    dataOutputWorker->setAffinity(reservedCpuCore);
    dataOutputWorker->start();

    // Set the data output handler to be the worker.
    lm::main::DataOutputQueue::setInstance(dataOutputWorker);

    // Create a table for tracking the replicate assignments.
    int* assignedSimulationsTable = new int[lm::MPI::worldSize];
    for (int i=0; i<lm::MPI::worldSize; i++)
    {
        assignedSimulationsTable[i] = 0;
    }

    // Get the maximum number of simulations that can be started on each process.
    int* maxSimulationsTable = new int[lm::MPI::worldSize];
    int maxSimulations = resourceAllocator.getMaxSimultaneousReplicates();
    MPI_EXCEPTION_CHECK(MPI_Gather(&maxSimulations, 1, MPI_INT, maxSimulationsTable, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    int simultaneousReplicates=0;
    for (int i=0; i<lm::MPI::worldSize; i++) simultaneousReplicates += maxSimulationsTable[i];
    Print::printf(Print::INFO, "Number of simultaneous replicates is %d", simultaneousReplicates);
    if (simultaneousReplicates == 0) throw Exception("Invalid configuration, no replicates can be processed.");

    // Create a table for the simulation status.
    map<int,int> simulationStatusTable;
    map<int,struct timespec> simulationStartTimeTable;
    for (vector<int>::iterator it=replicates.begin(); it<replicates.end(); it++)
    {
        simulationStatusTable[*it] = 0;
        simulationStartTimeTable[*it].tv_sec = 0;
        simulationStartTimeTable[*it].tv_nsec = 0;
    }

    // MPI message variables.
    int messageWaiting;
    MPI_Status messageStatus;
    void * staticDataBuffer = NULL;
    MPI_EXCEPTION_CHECK(MPI_Alloc_mem(lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE, MPI_INFO_NULL, &staticDataBuffer));
    int finishedMessage[2];

    // Get the simulation parameters and distribute them to the slaves.
    std::map<std::string,string> simulationParameters = file->getParameters();
    broadcastSimulationParameters(staticDataBuffer, simulationParameters);

    // Get the reaction model and distribute it to the slaves.
    lm::io::ReactionModel reactionModel;
    if (solverFactory.needsReactionModel())
    {
        file->getReactionModel(&reactionModel);
        broadcastReactionModel(staticDataBuffer, &reactionModel);
    }

    // Get the diffusion model and distribute it to the slaves.
    lm::io::DiffusionModel diffusionModel;
    uint8_t * lattice=NULL, * latticeSites=NULL;
    size_t latticeSize=0, latticeSitesSize=0;
    if (solverFactory.needsDiffusionModel())
    {
        file->getDiffusionModel(&diffusionModel);
        latticeSize = diffusionModel.lattice_x_size()*diffusionModel.lattice_y_size()*diffusionModel.lattice_z_size()*diffusionModel.particles_per_site();
        lattice = new uint8_t[latticeSize];
        latticeSitesSize = diffusionModel.lattice_x_size()*diffusionModel.lattice_y_size()*diffusionModel.lattice_z_size();
        latticeSites = new uint8_t[latticeSitesSize];
        file->getDiffusionModelLattice(&diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize);
        broadcastDiffusionModel(staticDataBuffer, &diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize);
    }

    // Distribute the simulations to the processes.
    Print::printf(Print::INFO, "Starting %d replicates from file %s.", replicates.size(), simulationFilename.c_str());
    int nextProcess=-1;
    list<ReplicateRunner *> runningReplicates;
    unsigned long long noopLoopCycles=0;
    while (!globalAbort)
    {
        // Increment the loop counter.
        noopLoopCycles++;

        // Check for any messages from the slaves.
        MPI_EXCEPTION_CHECK(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &messageWaiting, &messageStatus));

        // We have an output data message waiting.
        if (messageWaiting && messageStatus.MPI_TAG == lm::MPI::MSG_OUTPUT_DATA_STATIC)
        {
            PROF_BEGIN(PROF_MASTER_READ_STATIC_MSG);
            // Read the message into the buffer.
            MPI_EXCEPTION_CHECK(MPI_Recv(staticDataBuffer, lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE, MPI_BYTE, MPI_ANY_SOURCE, lm::MPI::MSG_OUTPUT_DATA_STATIC, MPI_COMM_WORLD, &messageStatus));

            // Get the size of the message.
            int messageSize;
            MPI_EXCEPTION_CHECK(MPI_Get_count(&messageStatus, MPI_BYTE, &messageSize));
            Print::printf(Print::VERBOSE_DEBUG, "Received output data set of size %d from process %d.", messageSize, messageStatus.MPI_SOURCE);

            // Add the message to the data set queue.
            ((lm::main::DataOutputQueue *)dataOutputWorker)->pushDataSet(staticDataBuffer, (size_t)messageSize);
            noopLoopCycles = 0;
            PROF_END(PROF_MASTER_READ_STATIC_MSG);
        }

        // We have a finished simulations message waiting.
        else if (messageWaiting && messageStatus.MPI_TAG == lm::MPI::MSG_SIMULATION_FINISHED)
        {
            PROF_BEGIN(PROF_MASTER_READ_FINISHED_MSG);
            MPI_EXCEPTION_CHECK(MPI_Recv(&finishedMessage, 2, MPI_INT, MPI_ANY_SOURCE, lm::MPI::MSG_SIMULATION_FINISHED, MPI_COMM_WORLD, &messageStatus));

			struct timespec now;
			#if defined(LINUX)
			clock_gettime(CLOCK_REALTIME, &now);
			#elif defined(MACOSX)
			struct timeval now2;
			gettimeofday(&now2, NULL);
			now.tv_sec = now2.tv_sec;
			now.tv_nsec = now2.tv_usec*1000;
			#endif

            Print::printf(Print::INFO, "Replicate %d completed by process %d with exit code %d in %0.2f seconds.", finishedMessage[0], messageStatus.MPI_SOURCE, finishedMessage[1], ((double)(now.tv_sec-simulationStartTimeTable[finishedMessage[0]].tv_sec))+1e-9*((double)now.tv_nsec-simulationStartTimeTable[finishedMessage[0]].tv_nsec));
            assignedSimulationsTable[messageStatus.MPI_SOURCE]--;
            simulationStatusTable[finishedMessage[0]] = 2;
            noopLoopCycles = 0;
            PROF_END(PROF_MASTER_READ_FINISHED_MSG);
        }

        // Check for finished simulations in our process.
        ReplicateRunner * finishedReplicate;
        while ((finishedReplicate=popNextFinishedReplicate(runningReplicates, resourceAllocator)) != NULL)
        {
            PROF_BEGIN(PROF_MASTER_FINISHED_THREAD);

			struct timespec now;
			#if defined(LINUX)
			clock_gettime(CLOCK_REALTIME, &now);
			#elif defined(MACOSX)
			struct timeval now2;
			gettimeofday(&now2, NULL);
			now.tv_sec = now2.tv_sec;
			now.tv_nsec = now2.tv_usec*1000;
			#endif

            Print::printf(Print::INFO, "Replicate %d completed by process %d with exit code %d in %0.2f seconds.", finishedReplicate->getReplicate(), lm::MPI::worldRank, finishedReplicate->getReplicateExitCode(), ((double)(now.tv_sec-simulationStartTimeTable[finishedReplicate->getReplicate()].tv_sec))+1e-9*((double)now.tv_nsec-simulationStartTimeTable[finishedReplicate->getReplicate()].tv_nsec));
            assignedSimulationsTable[lm::MPI::worldRank]--;
            simulationStatusTable[finishedReplicate->getReplicate()] = 2;
            noopLoopCycles = 0;
            finishedReplicate->stop();
            delete finishedReplicate; // We are responsible for deleting the replicate runner.
            PROF_END(PROF_MASTER_FINISHED_THREAD);
        }

        // If there were no messages, see if we need to start any new simulations and then wait a while.
        if (noopLoopCycles > 1000)
        {
            // Find a simulation to perform.
            int replicate = -1;
            bool allFinished=true;
            for (vector<int>::iterator it=replicates.begin(); it<replicates.end(); it++)
            {
                if (simulationStatusTable[*it] == 0)
                {
                    replicate = *it;
                    allFinished = false;
                    break;
                }
                if (simulationStatusTable[*it] != 2) allFinished = false;
            }

            // If all of the simulations are finished, we are done.
            if (allFinished) break;

            // Find a process to perform the simulation.
            if (replicate >= 0)
            {
                int destProc = -1;
                for (int i=0; i<lm::MPI::worldSize; i++)
                {
                    if (++nextProcess >= lm::MPI::worldSize) nextProcess = 0;
                    if (assignedSimulationsTable[nextProcess] < maxSimulationsTable[nextProcess])
                    {
                        destProc = nextProcess;
                        break;
                    }
                }

                // If we found a process to perform the simulation, send it a message.
                if (destProc != -1)
                {
                    // If the process is a slave, send it a message.
                    if (destProc != lm::MPI::MASTER)
                    {
                        MPI_EXCEPTION_CHECK(MPI_Send(&replicate, 1, MPI_INT, destProc, lm::MPI::MSG_RUN_SIMULATION, MPI_COMM_WORLD));
                    }
                    else
                    {
                        // Otherwise it must be us, so start the replicate.
                        runningReplicates.push_back(startReplicate(replicate, solverFactory, simulationParameters, &reactionModel, &diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resourceAllocator));
                    }
                    assignedSimulationsTable[destProc]++;
                    simulationStatusTable[replicate] = 1;
					struct timespec now;
					#if defined(LINUX)
					clock_gettime(CLOCK_REALTIME, &now);
					#elif defined(MACOSX)
					struct timeval now2;
					gettimeofday(&now2, NULL);
					now.tv_sec = now2.tv_sec;
					now.tv_nsec = now2.tv_usec*1000;
					#endif
					simulationStartTimeTable[replicate].tv_sec = now.tv_sec;
					simulationStartTimeTable[replicate].tv_nsec = now.tv_nsec;

                    continue;
                }
            }

            PROF_BEGIN(PROF_MASTER_SLEEP);
            unsigned int sleepTime = 1000000;
            if (noopLoopCycles > 2000) sleepTime = 10000000;
            if (noopLoopCycles > 2100) sleepTime = 100000000;
            if (noopLoopCycles >= 3000 && noopLoopCycles%1000 == 0)
            {
                int replicatesRunning=0, replicatesRemaining=0;
                for (vector<int>::iterator it=replicates.begin(); it<replicates.end(); it++)
                {
                    if (simulationStatusTable[*it] == 0)
                        replicatesRemaining++;
                    else if (simulationStatusTable[*it] == 1)
                        replicatesRunning++;
                }
                Print::printf(Print::INFO, "Master sleeping, waiting for %d replicates to finish, %d left to start.",replicatesRunning,replicatesRemaining);
            }
            struct timespec requested, remainder;
            requested.tv_sec  = 0;
            requested.tv_nsec = sleepTime;
            if (nanosleep(&requested, &remainder) != 0 && errno != EINTR) throw lm::Exception("Sleep failed.");
            PROF_END(PROF_MASTER_SLEEP);
        }
    }

    Print::printf(Print::INFO, "Master shutting down.");

    // Stop checkpointing.
    checkpointSignaler->stopCheckpointing();

    // Tell all of the slave processes to stop.
    for (int destProc=0; destProc<lm::MPI::worldSize; destProc++)
    {
        int exitCode=0;
        if (destProc != lm::MPI::MASTER)
            MPI_EXCEPTION_CHECK(MPI_Send(&exitCode, 1, MPI_INT, destProc, lm::MPI::MSG_EXIT, MPI_COMM_WORLD));
    }

    // Wait for all of the processes to exit.
    MPI_EXCEPTION_CHECK(MPI_Barrier(MPI_COMM_WORLD));

    // If this was a global abort, stop the workers quickly.
    if (globalAbort)
    {
        Print::printf(Print::WARNING, "Aborting worker threads.");
        lm::thread::WorkerManager::getInstance()->abortWorkers();
    }

    // Otherwise, let them finish at their own pace.
    else
    {
        Print::printf(Print::DEBUG, "Stopping worker threads.");
        lm::thread::WorkerManager::getInstance()->stopWorkers();
    }

    // Close the simulation file.
    delete file;
    Print::printf(Print::INFO, "Simulation file closed.");

    // Cleanup any resources.
    MPI_EXCEPTION_CHECK(MPI_Free_mem(staticDataBuffer));
    delete[] assignedSimulationsTable;
    delete[] maxSimulationsTable;
    delete checkpointSignaler;
    delete signalHandler;
    delete dataOutputWorker;
    if (lattice != NULL) {delete [] lattice; lattice = NULL;}
    if (latticeSites != NULL) {delete [] latticeSites; latticeSites = NULL;}

    Print::printf(Print::DEBUG, "MPI master process %d finished.", lm::MPI::worldRank);
}

void executeSimulationMPISingleSlave()
{
    Print::printf(Print::DEBUG, "MPI slave process %d started.", lm::MPI::worldRank);

    // Create the queue to handle data output.
    lm::main::MPIRemoteDataOutputQueue * dataOutputQueue = new lm::main::MPIRemoteDataOutputQueue();
    lm::main::DataOutputQueue::setInstance(dataOutputQueue);

    // Create the resource allocator.
    #ifdef OPT_CUDA
    ResourceAllocator resourceAllocator(lm::MPI::worldRank, numberCpuCores, cpuCoresPerReplicate, cudaDevices, cudaDevicesPerReplicate);
    #else
    ResourceAllocator resourceAllocator(lm::MPI::worldRank, numberCpuCores, cpuCoresPerReplicate);
    #endif

    // Report the max simultaneous simulations to the master.
    int maxSimulations = resourceAllocator.getMaxSimultaneousReplicates();
    MPI_EXCEPTION_CHECK(MPI_Gather(&maxSimulations, 1, MPI_INT, NULL, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));

    // MPI message variables.
    int messageWaiting;
    MPI_Status messageStatus;
    void * staticDataBuffer = NULL;
    MPI_EXCEPTION_CHECK(MPI_Alloc_mem(lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE, MPI_INFO_NULL, &staticDataBuffer));
    int finishedMessage[2];

    // Read the simulation parameters.
    std::map<std::string,string> simulationParameters = receiveSimulationParameters(staticDataBuffer);

    // Read the reaction model.
    lm::io::ReactionModel reactionModel;
    if (solverFactory.needsReactionModel())
    {
        receiveReactionModel(staticDataBuffer, &reactionModel);
    }

    // Read the diffusion model.
    lm::io::DiffusionModel diffusionModel;
    uint8_t * lattice=NULL, * latticeSites=NULL;
    size_t latticeSize=0, latticeSitesSize=0;
    if (solverFactory.needsDiffusionModel())
    {
        receiveDiffusionModel(staticDataBuffer, &diffusionModel, &lattice, &latticeSize, &latticeSites, &latticeSitesSize);
    }

    // Execute simulations.
    list<ReplicateRunner *> runningReplicates;
    unsigned long long noopLoopCycles=0;
    while (true)
    {
        // Increment the loop counter.
        noopLoopCycles++;

        // See if there are any data sets to write.
        if (!dataOutputQueue->isEmpty())
        {
            Print::printf(Print::VERBOSE_DEBUG, "Sending output data set from process %d.", lm::MPI::worldRank);

            // Put the next data set into the send buffer.
            size_t messageSize=dataOutputQueue->popDataSetIntoBuffer(staticDataBuffer, lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE);

            // Send the message.
            MPI_EXCEPTION_CHECK(MPI_Send(staticDataBuffer, messageSize, MPI_BYTE, lm::MPI::MASTER, lm::MPI::MSG_OUTPUT_DATA_STATIC, MPI_COMM_WORLD));
            noopLoopCycles = 0;
            continue;
        }

        // See if there are any finished simulations.
        ReplicateRunner * finishedReplicate;
        while ((finishedReplicate=popNextFinishedReplicate(runningReplicates, resourceAllocator)) != NULL)
        {
            finishedMessage[0] = finishedReplicate->getReplicate();
            finishedMessage[1] = finishedReplicate->getReplicateExitCode();
            finishedReplicate->stop();
            delete finishedReplicate; // We are responsible for deleting the replicate runner.
            MPI_EXCEPTION_CHECK(MPI_Send(&finishedMessage, 2, MPI_INT, lm::MPI::MASTER, lm::MPI::MSG_SIMULATION_FINISHED, MPI_COMM_WORLD));
            noopLoopCycles = 0;
        }

        // See if there is a simulation to run
        MPI_EXCEPTION_CHECK(MPI_Iprobe(lm::MPI::MASTER, lm::MPI::MSG_RUN_SIMULATION, MPI_COMM_WORLD, &messageWaiting, &messageStatus));
        if (messageWaiting)
        {
            // Get the replicate number and start the simulation.
            int replicate;
            MPI_EXCEPTION_CHECK(MPI_Recv(&replicate, 1, MPI_INT, lm::MPI::MASTER, lm::MPI::MSG_RUN_SIMULATION, MPI_COMM_WORLD, &messageStatus));
            runningReplicates.push_back(startReplicate(replicate, solverFactory, simulationParameters, &reactionModel, &diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resourceAllocator));
            noopLoopCycles = 0;
            continue;
        }

        // See if there is an exit message.
        MPI_EXCEPTION_CHECK(MPI_Iprobe(lm::MPI::MASTER, lm::MPI::MSG_EXIT, MPI_COMM_WORLD, &messageWaiting, &messageStatus));
        if (messageWaiting)
        {
            int exitCode;
            MPI_EXCEPTION_CHECK(MPI_Recv(&exitCode, 1, MPI_INT, lm::MPI::MASTER, lm::MPI::MSG_EXIT, MPI_COMM_WORLD, &messageStatus));
            if (!dataOutputQueue->isEmpty()) throw lm::Exception("MPI slave process finished with queued data sets", lm::MPI::worldRank);
            noopLoopCycles = 0;
            break;
        }

        // If we haven't done anything for a while, sleep a little.
        if (noopLoopCycles > 1000)
        {
            PROF_BEGIN(PROF_SLAVE_SLEEP);
            struct timespec requested, remainder;
            requested.tv_sec  = 0;
            requested.tv_nsec = 1000000;
            if (noopLoopCycles > 2000) requested.tv_nsec = 10000000;
            if (noopLoopCycles > 2100) requested.tv_nsec = 100000000;
            if (nanosleep(&requested, &remainder) != 0 && errno != EINTR) throw lm::Exception("Sleep failed.");
            PROF_END(PROF_SLAVE_SLEEP);
        }
    }

    // Wait for all of the processes to exit.
    MPI_EXCEPTION_CHECK(MPI_Barrier(MPI_COMM_WORLD));

    // Destroy the data output queue.
    delete dataOutputQueue;
    dataOutputQueue = NULL;

    // Cleanup any resources.
    MPI_EXCEPTION_CHECK(MPI_Free_mem(staticDataBuffer));
    if (lattice != NULL) {delete [] lattice; lattice = NULL;}
    if (latticeSites != NULL) {delete [] latticeSites; latticeSites = NULL;}

    Print::printf(Print::DEBUG, "MPI slave process %d finished.", lm::MPI::worldRank);
}

void broadcastSimulationParameters(void * staticDataBuffer, map<string,string> & simulationParameters)
{
    // Send the simulation parameters to the nodes.
    lm::message::SimulationParameters msg;
    lm::io::SimulationParameters::intoMessage(msg, simulationParameters);
    int msgSize = msg.ByteSize();
    if (msgSize > lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE) throw Exception("Simulation parameter message exceeded buffer size");
    msg.SerializeToArray(staticDataBuffer, msgSize);
    MPI_EXCEPTION_CHECK(MPI_Bcast(&msgSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(staticDataBuffer, msgSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
}

void broadcastReactionModel(void * staticDataBuffer, lm::io::ReactionModel * reactionModel)
{
    // Send the reaction model to the nodes.
    int msgSize = reactionModel->ByteSize();
    if (msgSize > lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE) throw Exception("Reaction model message exceeded buffer size");
    reactionModel->SerializeToArray(staticDataBuffer, msgSize);
    MPI_EXCEPTION_CHECK(MPI_Bcast(&msgSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(staticDataBuffer, msgSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
}

void broadcastDiffusionModel(void * staticDataBuffer, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize)
{
    // Send the diffusion model to the nodes.
    int msgSize = diffusionModel->ByteSize();
    if (msgSize > lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE) throw Exception("Diffusion model message exceeded buffer size");
    diffusionModel->SerializeToArray(staticDataBuffer, msgSize);
    MPI_EXCEPTION_CHECK(MPI_Bcast(&msgSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(staticDataBuffer, msgSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(&latticeSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(lattice, latticeSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(&latticeSitesSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(latticeSites, latticeSitesSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
}

map<string,string> receiveSimulationParameters(void * staticDataBuffer)
{
    int msgSize;
    MPI_EXCEPTION_CHECK(MPI_Bcast(&msgSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(staticDataBuffer, msgSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
    lm::message::SimulationParameters msg;
    msg.ParseFromArray(staticDataBuffer, msgSize);
    std::map<std::string,string> simulationParameters = lm::io::SimulationParameters::fromMessage(msg);
    Print::printf(Print::DEBUG, "Process %d received simulation parameters: %d parameters", lm::MPI::worldRank, simulationParameters.size());
    return simulationParameters;
}

void receiveReactionModel(void * staticDataBuffer, lm::io::ReactionModel * reactionModel)
{
    int msgSize;
    MPI_EXCEPTION_CHECK(MPI_Bcast(&msgSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(staticDataBuffer, msgSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
    reactionModel->ParseFromArray(staticDataBuffer, msgSize);
    Print::printf(Print::DEBUG, "Process %d received reaction model: %d bytes", lm::MPI::worldRank, msgSize);
}

void receiveDiffusionModel(void * staticDataBuffer, lm::io::DiffusionModel * diffusionModel, uint8_t ** lattice, size_t * latticeSize, uint8_t ** latticeSites, size_t * latticeSitesSize)
{
    int msgSize;
    MPI_EXCEPTION_CHECK(MPI_Bcast(&msgSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(staticDataBuffer, msgSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
    diffusionModel->ParseFromArray(staticDataBuffer, msgSize);
    Print::printf(Print::DEBUG, "Process %d received diffusion model: %d bytes", lm::MPI::worldRank, msgSize);

    MPI_EXCEPTION_CHECK(MPI_Bcast(latticeSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    *lattice = new uint8_t[*latticeSize];
    MPI_EXCEPTION_CHECK(MPI_Bcast(*lattice, *latticeSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
    MPI_EXCEPTION_CHECK(MPI_Bcast(latticeSitesSize, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    *latticeSites = new uint8_t[*latticeSitesSize];
    MPI_EXCEPTION_CHECK(MPI_Bcast(*latticeSites, *latticeSitesSize, MPI_BYTE, lm::MPI::MASTER, MPI_COMM_WORLD));
    Print::printf(Print::DEBUG, "Process %d received diffusion model lattice: %d and %d bytes", lm::MPI::worldRank, *latticeSize, *latticeSitesSize);
}

ReplicateRunner * startReplicate(int replicate, MESolverFactory solverFactory, std::map<std::string,string> & simulationParameters, lm::io::ReactionModel * reactionModel, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize, ResourceAllocator & resourceAllocator)
{
    // Allocate resources for the replicate.
    ResourceAllocator::ComputeResources resources = resourceAllocator.assignReplicate(replicate);

    // Start a new thread for the replicate.
    Print::printf(Print::DEBUG, "Starting replicate %d in process %d (%s).", replicate, lm::MPI::worldRank, resources.toString().c_str());
    ReplicateRunner * runner = new ReplicateRunner(replicate, solverFactory, &simulationParameters, reactionModel, diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resources);
    runner->start();
    return runner;
}

ReplicateRunner * popNextFinishedReplicate(list<ReplicateRunner *> & runningReplicates, ResourceAllocator & resourceAllocator)
{
    for (list<ReplicateRunner *>::iterator it=runningReplicates.begin(); it != runningReplicates.end(); it++)
    {
        ReplicateRunner * runner = *it;
        if (runner->hasReplicateFinished())
        {
            runningReplicates.erase(it);
            resourceAllocator.removeReplicate(runner->getReplicate());
            return runner;
        }
    }
    return NULL;
}
