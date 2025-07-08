/*
 * University of Illinois Open Source License Copyright 2008-2013
 * Luthey-Schulten Group, Copyright 2012 Roberts Group, All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group University of Illinois at
 * Urbana-Champaign http://www.scs.uiuc.edu/~schulten
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the Software), to deal
 * with the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * - Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimers.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * WITH THE SOFTWARE.
 *
 * Author(s): Mike Hallock, Tyler Earnest
 *
 ***************************************************************************
 *
 * runSimulation and runSolver provide a stripped-down environment for
 * executing a single replicate.  The motivation is to have a non entry point
 * method to call for setting up a Solver and realizing the trajectory.
 * Primarily, this is to provide the python interface an easy way to run a
 * Solver without requiring setting up the output queue or dealing with
 * resource allocation.  
 */
#include "config.h"
#if defined(MACOSX)
#include <sys/sysctl.h>
#elif defined(LINUX)
#include <sys/sysinfo.h>
#endif


#include <string>
#include <map>
#include <cstdio>
#include <cstring>
#include <ctime>
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
#ifdef OPT_CUDA
#include "cuda/lm_cuda.h"
#include "rdme/MGPUMpdRdmeSolver.h"
#include "rdme/MGPUIntMpdRdmeSolver.h"
#include <cuda.h>
#endif
#include "io/lm_hdf5.h"
#include "io/SimulationFile.h"
#include "DiffusionModel.pb.h"
#include "ReactionModel.pb.h"
#include "io/SimulationParameters.h"
#include "core/CheckpointSignaler.h"
#include "core/DataOutputQueue.h"
#include "core/LocalDataOutputWorker.h"
#include "core/SignalHandler.h"
#include "core/ReplicateRunner.h"
#include "core/ResourceAllocator.h"
#include "SimulationParameters.pb.h"
#include "thread/Thread.h"
#include "thread/WorkerManager.h"
#include "lptf/Profile.h"
#include "me/MESolverFactory.h"
#include "me/MESolver.h"
#include "cme/CMESolver.h"
#include "rdme/RDMESolver.h"
#include "core/runSimulation.h"

using std::map;
using std::list;
using lm::Print;
using lm::Exception;
using lm::main::ReplicateRunner;
using lm::main::ResourceAllocator;
using lm::me::MESolverFactory;
using lm::thread::PthreadException;

extern bool globalAbort;

void abortHandler(int x)
{
	globalAbort=1;
}

void runSimulation(char *simulationFilename, int replicate, char *solverClass,
               vector<int> cudaDevices, time_t checkpointInterval)
{
    lm::me::MESolverFactory solverFactory;
    solverFactory.setSolver(solverClass);
    lm::me::MESolver * solver = NULL;
    solver = solverFactory.instantiate();

    runSolver(simulationFilename, replicate, solver, cudaDevices,checkpointInterval);

    // Free the solver we created
    if (solver != NULL) {delete solver; solver = NULL;}

}

void runSolver(char *simulationFilename, int replicate, lm::me::MESolver *solver,
               vector<int> cudaDevices, time_t checkpointInterval)
{
    PROF_SET_THREAD(0);
    PROF_BEGIN(PROF_SIM_RUN);
    double cpuCoresPerReplicate  = 1.0*cudaDevices.size();
    double cudaDevicesPerReplicate = 1.0*cudaDevices.size();
    Print::printf(Print::DEBUG, "interactive process started.");

#ifdef OPT_CUDA
    if (dynamic_cast<lm::rdme::MGPUMpdRdmeSolver*> (solver) == NULL && dynamic_cast<lm::rdme::MGPUIntMpdRdmeSolver*> (solver) == NULL && cudaDevices.size() > 1) {
        throw lm::Exception("Multiple CUDA devices specified for single GPU solver");
    }
#endif

    signal(SIGINT, abortHandler);

    // Get the number of processors. ( copied from Main.cpp )
    int numberCpuCores = 0;
#if defined(MACOSX)
    uint physicalCpuCores;
    size_t  physicalCpuCoresSize=sizeof(physicalCpuCores);
    sysctlbyname("hw.activecpu",&physicalCpuCores,&physicalCpuCoresSize,NULL,0);
    numberCpuCores = physicalCpuCores;
#elif defined(LINUX)
    #ifdef ARM
        numberCpuCores = get_nprocs_conf();
    #else
        numberCpuCores = get_nprocs();
    #endif
#else
    #error "Unsupported architecture."
#endif

    // Create the resource allocator, subtract one core for the data output thread.
#ifdef OPT_CUDA
    Print::printf(Print::INFO, "Using %d processor(s) and %zd CUDA device(s) per process.", numberCpuCores, cudaDevices.size());
    Print::printf(Print::INFO, "Assigning %0.2f processor(s) and %0.2f CUDA device(s) per replicate.", cpuCoresPerReplicate, cudaDevicesPerReplicate);
    ResourceAllocator resourceAllocator(0, numberCpuCores, cpuCoresPerReplicate, cudaDevices, cudaDevicesPerReplicate);
#else
    Print::printf(Print::INFO, "Using %d processor(s) per process.", numberCpuCores);
    Print::printf(Print::INFO, "Assigning %0.2f processor(s) per replicate.", cpuCoresPerReplicate);
    ResourceAllocator resourceAllocator(0, numberCpuCores, cpuCoresPerReplicate);
#endif

    ResourceAllocator::ComputeResources resources = resourceAllocator.assignReplicate(replicate);

    lm::main::CheckpointSignaler * checkpointSignaler = NULL;
    if (checkpointInterval>0) {
        checkpointSignaler = new lm::main::CheckpointSignaler();
        checkpointSignaler->setAffinity(0);
        checkpointSignaler->startCheckpointing(checkpointInterval);
        checkpointSignaler->start();
    }


    // Open the file.
    lm::io::hdf5::SimulationFile * file = new lm::io::hdf5::SimulationFile(simulationFilename);

    // Start the data output thread.
    lm::main::LocalDataOutputWorker * dataOutputWorker = new lm::main::LocalDataOutputWorker(file);
    //dataOutputWorker->setAffinity(reservedCpuCore);
    dataOutputWorker->start();

    // Set the data output handler to be the worker.
    lm::main::DataOutputQueue::setInstance(dataOutputWorker);

    // Get the simulation parameters.
    std::map<std::string,string> simulationParameters = file->getParameters();

    // Get the reaction model.
    lm::io::ReactionModel reactionModel;
    if (solver->needsReactionModel())
    {
        file->getReactionModel(&reactionModel);
    }

    // Get the diffusion model.
    lm::io::DiffusionModel diffusionModel;
    uint8_t * lattice=NULL, * latticeSites=NULL;
    size_t latticeSize=0, latticeSitesSize=0;
    if (solver->needsDiffusionModel())
    {
        file->getDiffusionModel(&diffusionModel);
        latticeSize = diffusionModel.lattice_x_size()*diffusionModel.lattice_y_size()*diffusionModel.lattice_z_size()*diffusionModel.particles_per_site()*diffusionModel.bytes_per_particle();
        lattice = new uint8_t[latticeSize];
        latticeSitesSize = diffusionModel.lattice_x_size()*diffusionModel.lattice_y_size()*diffusionModel.lattice_z_size();
        latticeSites = new uint8_t[latticeSitesSize];
        file->getDiffusionModelLattice(&diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize);
    }

    // Start a new thread for the replicate.
    Print::printf(Print::DEBUG, "Starting replicate %d (%s).", replicate, resources.toString().c_str());

    #if defined(OPT_CUDA)
    // Set the GPU affinity.
    if (resources.cudaDevices.size() > 0)
            lm::CUDA::setCurrentDevice(resources.cudaDevices[0]);
    #endif

    try
    {
        // Run the simulation using the specified solver.
        solver->initialize(replicate, &simulationParameters, &resources);
        if (solver->needsReactionModel())
        {
            ((lm::cme::CMESolver *)solver)->setReactionModel(&reactionModel);
        }
        if (solver->needsDiffusionModel())
        {
            ((lm::rdme::RDMESolver *)solver)->setDiffusionModel(&diffusionModel, lattice, latticeSize, latticeSites, latticeSitesSize);
        }
        solver->generateTrajectory();
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

    if (checkpointInterval>0)
        checkpointSignaler->stopCheckpointing();


    Print::printf(Print::DEBUG, "Stopping worker threads.");
    lm::thread::WorkerManager::getInstance()->stopWorkers();

    // Close the simulation file.
    delete file;
    Print::printf(Print::INFO, "Simulation file closed.");


    // Cleanup any resources.
    if (checkpointInterval>0)
        delete checkpointSignaler;
    delete dataOutputWorker;
    if (lattice != NULL) {delete [] lattice; lattice = NULL;}
    if (latticeSites != NULL) {delete [] latticeSites; latticeSites = NULL;}

    Print::printf(Print::DEBUG, "Master process finished.");


    PROF_END(PROF_SIM_RUN);
}

