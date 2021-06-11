/*
 * University of Illinois Open Source License
 * Copyright 2014-2018 Luthey-Schulten Group,
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
 * Author(s): Mike Hallock
 */

#include <iostream>
#include <string>
#include <map>
#include <cstdio>
#include <cstring>
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
#include "SimulationParameters.h"
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

#include "me/MESolverFactory.h"
#include "rdme/RDMESolver.h"

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
int runReplicateNOW(int replicate, MESolverFactory solverFactory, std::map<std::string,string> & simulationParameters, lm::io::ReactionModel * reactionModel, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize, ResourceAllocator & resourceAllocator);

//ReplicateRunner * startReplicate(int replicate, MESolverFactory solverFactory, std::map<std::string,string> & simulationParameters, lm::io::ReactionModel * reactionModel, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize, ResourceAllocator & resourceAllocator);
// ReplicateRunner * popNextFinishedReplicate(list<ReplicateRunner *> & runningReplicates, ResourceAllocator & resourceAllocator);

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
                parseArguments(argc, argv, "lm::rdme::MPIMpdRdmeSolver");

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
            parseArguments(argc, argv, "lm::rdme::MPIMpdRdmeSolver");

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
        Print::printf(Print::DEBUG, "Closing MPI library.");
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
	// Print the hostname for each rank
	char host[MPI_MAX_PROCESSOR_NAME];
	int hostlen;
	MPI_Get_processor_name(host, &hostlen);
    Print::printf(Print::INFO, "MPI master process %d started on %s.", lm::MPI::worldRank, host);

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

/* TODO WTF checkpoint? 
    // Create the checkpoint signaler.
    lm::main::CheckpointSignaler * checkpointSignaler = new lm::main::CheckpointSignaler();
    checkpointSignaler->setAffinity(reservedCpuCore);
    checkpointSignaler->start();
    checkpointSignaler->startCheckpointing(checkpointInterval);
*/

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
    //int maxSimulations = resourceAllocator.getMaxSimultaneousReplicates();
    //MPI_EXCEPTION_CHECK(MPI_Gather(&maxSimulations, 1, MPI_INT, maxSimulationsTable, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
    int simultaneousReplicates=1;
    //for (int i=0; i<lm::MPI::worldSize; i++) simultaneousReplicates += maxSimulationsTable[i];
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

    void * staticDataBuffer = NULL;
    MPI_EXCEPTION_CHECK(MPI_Alloc_mem(lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE, MPI_INFO_NULL, &staticDataBuffer));

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
    list<ReplicateRunner *> runningReplicates;

	// Find a simulation to perform.
	int replicate = -1;
	for (vector<int>::iterator it=replicates.begin(); it<replicates.end(); it++)
	{
		if (simulationStatusTable[*it] == 0)
		{
			replicate = *it;
		}
		else
		{
			// Don't think this can happen.
			printf("Skip replicate %d with status %d\n", *it, simulationStatusTable[*it]);
			continue;
		}

		// bookeeping. useful?
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

		// Broadcast replicate number to start slaves
		MPI_EXCEPTION_CHECK(MPI_Bcast(&replicate, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));

		// run replicate locally
		int status = runReplicateNOW(replicate, solverFactory, simulationParameters, &reactionModel, &diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resourceAllocator);

		#if defined(LINUX)
		  clock_gettime(CLOCK_REALTIME, &now);
		#elif defined(MACOSX)
		  gettimeofday(&now2, NULL);
		  now.tv_sec = now2.tv_sec;
		  now.tv_nsec = now2.tv_usec*1000;
		#endif

		Print::printf(Print::INFO, "Replicate %d completed by process %d with exit code %d in %0.2f seconds.", replicate, lm::MPI::worldRank, status, ((double)(now.tv_sec-simulationStartTimeTable[replicate].tv_sec))+1e-9*((double)now.tv_nsec-simulationStartTimeTable[replicate].tv_nsec));
		assignedSimulationsTable[lm::MPI::worldRank]--;
		simulationStatusTable[replicate] = 2;

		// Wait for all of the processes to complete replicate.
		MPI_EXCEPTION_CHECK(MPI_Barrier(MPI_COMM_WORLD));

    }

    Print::printf(Print::INFO, "Master shutting down.");
/*
    // Stop checkpointing.
    checkpointSignaler->stopCheckpointing();
*/

	replicate=-1;
    // Tell all of the slave processes to stop.
		MPI_EXCEPTION_CHECK(MPI_Bcast(&replicate, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));
/*
    for (int destProc=0; destProc<lm::MPI::worldSize; destProc++)
    {
        int exitCode=0;
        if (destProc != lm::MPI::MASTER)
            MPI_EXCEPTION_CHECK(MPI_Send(&exitCode, 1, MPI_INT, destProc, lm::MPI::MSG_EXIT, MPI_COMM_WORLD));
    }

    // Wait for all of the processes to exit.
    MPI_EXCEPTION_CHECK(MPI_Barrier(MPI_COMM_WORLD));
*/

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
    //delete checkpointSignaler;
    delete signalHandler;
    delete dataOutputWorker;
    if (lattice != NULL) {delete [] lattice; lattice = NULL;}
    if (latticeSites != NULL) {delete [] latticeSites; latticeSites = NULL;}

    Print::printf(Print::DEBUG, "MPI master process %d finished.", lm::MPI::worldRank);
}

void executeSimulationMPISingleSlave()
{
	// Print the hostname for each rank
	char host[MPI_MAX_PROCESSOR_NAME];
	int hostlen;
	MPI_Get_processor_name(host, &hostlen);
    Print::printf(Print::INFO, "MPI slave process %d started on %s.", lm::MPI::worldRank, host);

    // Create the queue to handle data output.
    lm::main::MPIRemoteDataOutputQueue * dataOutputQueue = new lm::main::MPIRemoteDataOutputQueue();
    lm::main::DataOutputQueue::setInstance(dataOutputQueue);

    // Create the resource allocator.
    #ifdef OPT_CUDA
    ResourceAllocator resourceAllocator(lm::MPI::worldRank, numberCpuCores, cpuCoresPerReplicate, cudaDevices, cudaDevicesPerReplicate);
    #else
    ResourceAllocator resourceAllocator(lm::MPI::worldRank, numberCpuCores, cpuCoresPerReplicate);
    #endif

    void * staticDataBuffer = NULL;
    MPI_EXCEPTION_CHECK(MPI_Alloc_mem(lm::MPI::OUTPUT_DATA_STATIC_MAX_SIZE, MPI_INFO_NULL, &staticDataBuffer));

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
    while (true)
    {
		// Get the replicate number and start the simulation.
		int replicate;

		// Get Broadcast replicate number
		MPI_EXCEPTION_CHECK(MPI_Bcast(&replicate, 1, MPI_INT, lm::MPI::MASTER, MPI_COMM_WORLD));

		if(replicate == -1)
			break;

		// run replicate locally
		int status = runReplicateNOW(replicate, solverFactory, simulationParameters, &reactionModel, &diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize, resourceAllocator);

		if(!status)
			Print::printf(Print::DEBUG, "Replicate failed with status: %d.", status);

		// Wait for all of the processes to exit.
		MPI_EXCEPTION_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    }


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

int runReplicateNOW(int replicate, MESolverFactory solverFactory, std::map<std::string,string> & simulationParameters, lm::io::ReactionModel * reactionModel, lm::io::DiffusionModel * diffusionModel, uint8_t * lattice, size_t latticeSize, uint8_t * latticeSites, size_t latticeSitesSize, ResourceAllocator & resourceAllocator)
{
    PROF_BEGIN(PROF_REPLICATE_EXECUTE);

	ResourceAllocator::ComputeResources resources = resourceAllocator.assignReplicate(replicate);

    MESolver * solver = NULL;
    int status = -1;
    try
    {
        // Run the simulation using the specified solver.
        solver = solverFactory.instantiate();
        solver->initialize(replicate, &simulationParameters, &resources);
        if (solver->needsReactionModel())
        {
            ((lm::cme::CMESolver *)solver)->setReactionModel(reactionModel);
        }
        if (solver->needsDiffusionModel())
        {
            ((lm::rdme::RDMESolver *)solver)->setDiffusionModel(diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize);
        }
        solver->generateTrajectory();
        status = 0;
    }
    catch (lm::Exception e)
    {
        Print::printf(Print::ERROR, "Exception during execution of replicate: %s.", e.what());
    }
    catch (std::exception& e)
    {
        Print::printf(Print::ERROR, "Std exception during execution of replicate: %s.", e.what());
    }
    catch (...)
    {
        Print::printf(Print::ERROR, "Unknown exception during execution of replicate.");
    }

    // Free any resources.
    if (solver != NULL) {delete solver; solver = NULL;}
	
    PROF_END(PROF_REPLICATE_EXECUTE);

	return status;
}

