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

#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <sys/stat.h>
#include "config.h"
#if defined(MACOSX)
#include <sys/sysctl.h>
#elif defined(LINUX)
#include <sys/sysinfo.h>
#endif
#ifdef OPT_CUDA
#include "cuda/lm_cuda.h"
#endif
#include "core/Types.h"
#include "core/Exceptions.h"
#include "core/Print.h"
#include "cmd/common.h"
#include "core/Globals.h"
#include "me/MESolverFactory.h"

using std::string;
using std::vector;

/**
 * Prints the copyright notice.
 */

/**
 * Gets the number of physical cpu cores on the system.
 */
int getPhysicalCpuCores()
{
    // Get the number of processors.
    #if defined(MACOSX)
    uint physicalCpuCores;
    size_t  physicalCpuCoresSize=sizeof(physicalCpuCores);
    sysctlbyname("hw.activecpu",&physicalCpuCores,&physicalCpuCoresSize,NULL,0);
    return physicalCpuCores;
    #elif defined(LINUX)
	#ifdef ARM
	return get_nprocs_conf();
	#else
	return get_nprocs();
	#endif
    #else
    #error "Unsupported architecture."
    #endif
}

/**
 * Parses the command line arguments.
 */
void parseArguments(int argc, char** argv, const char* defaultSolver)
{
    // Set any default options.
    replicates.clear();
    replicates.push_back(1);

    numberCpuCores = getPhysicalCpuCores();
    cpuCoresPerReplicate = 1.0;

    solverFactory.setSolver(defaultSolver);

    #ifdef OPT_CUDA
    cudaDevices.clear();
    for (int i=0; i<lm::CUDA::getNumberDevices(); i++)
        cudaDevices.push_back(i);
    cudaDevicesPerReplicate = 1.0;
    shouldPrintCudaCapabilities = true;

    cudaDevicesPerNode=1;
    #endif
    shouldReserveOutputCore = true;


    // Parse any arguments.
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;
        
        //See if the user is trying to get help.
        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
        	functionOption = "help";
        	break;
        }

        //See if the user is trying to get the version info.
        else if (strcmp(option, "-v") == 0 || strcmp(option, "--version") == 0) {
        	functionOption = "version";
            break;
        }



        //See if the user is trying to set the verbosity level
        else if ((strcmp(option, "-V") == 0 || strcmp(option, "--verbose") == 0) && i < (argc-1))
        {
            lm::Print::verbosityLevel(atoi(argv[++i]));
        }
        else if (strncmp(option, "--verbose=", strlen("--verbose=")) == 0)
        {
            lm::Print::verbosityLevel(atoi(option+strlen("--verbose=")));
        }


        //See if the user is trying to get the device info.
        else if (strcmp(option, "-l") == 0 || strcmp(option, "--list-devices") == 0) {
            functionOption = "devices";
        }

        //See if the user is trying to execute a simulation.
        else if (strcmp(option, "-f") == 0 || strcmp(option, "--file") == 0)
        {
            functionOption = "simulation";

            // Get the filename.
            if (i < argc-1)
                simulationFilename = argv[++i];
            else
                throw lm::CommandLineArgumentException("missing simulation filename.");
        }

        //See if the user is trying to set the replicates.
        else if ((strcmp(option, "-r") == 0 || strcmp(option, "--replicates") == 0) && i < (argc-1))
        {
            parseIntListArg(replicates, argv[++i]);
        }
        else if (strncmp(option, "--replicates=", strlen("--replicates=")) == 0)
        {
            parseIntListArg(replicates, option+strlen("--replicates="));
        }

        //See if the user is trying to set the checkpoint interval.
        else if ((strcmp(option, "-ck") == 0 || strcmp(option, "--checkpoint") == 0) && i < (argc-1))
        {
            checkpointInterval=parseTimeArg(argv[++i]);
        }
        else if (strncmp(option, "--checkpoint=", strlen("--checkpoint=")) == 0)
        {
            checkpointInterval=parseTimeArg(option+strlen("--checkpoint="));
        }

        //See if the user is trying to set the solver.
		#ifdef OPT_CUDA
        else if ((strcmp(option, "-sp") == 0 || strcmp(option, "--spatially-resolved") == 0))
        {
            solverFactory.setSolver("lm::rdme::MpdRdmeSolver");
        }
		#else
        else if ((strcmp(option, "-sp") == 0 || strcmp(option, "--spatially-resolved") == 0))
        {
            solverFactory.setSolver("lm::rdme::NextSubvolumeSolver");
        }
		#endif
        else if ((strcmp(option, "-ws") == 0 || strcmp(option, "--well-stirred") == 0))
        {
            solverFactory.setSolver("lm::cme::GillespieDSolver");
        }
        else if ((strcmp(option, "-m") == 0 || strcmp(option, "--model") == 0) && i < (argc-1))
        {
            solverFactory.setSolver(argv[++i]);
        }
        else if (strncmp(option, "--model=", strlen("--model=")) == 0)
        {
            solverFactory.setSolver(option+strlen("--model="));
        }
        else if ((strcmp(option, "-sl") == 0 || strcmp(option, "--solver") == 0) && i < (argc-1))
        {
            solverFactory.setSolver(argv[++i]);
        }
        else if (strncmp(option, "--solver=", strlen("--solver=")) == 0)
        {
            solverFactory.setSolver(option+strlen("--solver="));
        }

        //See if the user is trying to set the number of cpus.
        else if ((strcmp(option, "-c") == 0 || strcmp(option, "--cpu") == 0) && i < (argc-1))
        {
            numberCpuCores=atoi(argv[++i]);
        }
        else if (strncmp(option, "--cpu=", strlen("--cpu=")) == 0)
        {
            numberCpuCores=atoi(option+strlen("--cpu="));
        }

        //See if the user is trying to set the number of cuda devices per replicate.
         else if ((strcmp(option, "-cr") == 0 || strcmp(option, "--cpus-per-replicate") == 0) && i < (argc-1))
         {
             cpuCoresPerReplicate=parseIntReciprocalArg(argv[++i]);
         }
         else if (strncmp(option, "--cpus-per-replicate=", strlen("--cpus-per-replicate=")) == 0)
         {
             cpuCoresPerReplicate=parseIntReciprocalArg(option+strlen("--cpus-per-replicate="));
         }

#ifdef OPT_CUDA
        //See if the user is trying to set the cuda devices.
         else if ((strcmp(option, "-g") == 0 || strcmp(option, "--gpu") == 0) && i < (argc-1))
         {
             parseIntListArg(cudaDevices, argv[++i]);
         }
         else if (strncmp(option, "--gpu=", strlen("--gpu=")) == 0)
         {
             parseIntListArg(cudaDevices, option+strlen("--gpu="));
         }

        //See if the user is trying to set the number of cuda devices per replicate.
         else if ((strcmp(option, "-gr") == 0 || strcmp(option, "--gpus-per-replicate") == 0) && i < (argc-1))
         {
             cudaDevicesPerReplicate=parseIntReciprocalArg(argv[++i]);
         }
         else if (strncmp(option, "--gpus-per-replicate=", strlen("--gpus-per-replicate=")) == 0)
         {
             cudaDevicesPerReplicate=parseIntReciprocalArg(option+strlen("--gpus-per-replicate="));
         }
             
        //See if the user is trying to turn off cuda capability printing.
         else if ((strcmp(option, "-nc") == 0 || strcmp(option, "--no-capabilities") == 0))
         {
             shouldPrintCudaCapabilities = false;
         }
		// Prevent GPU Peering
		else if (strncmp(option, "--nopeer", strlen("--nopeer")) == 0)
		{
			mgpu_disablePeering=true;
		}

		// Set the number of GPUs per node for MPI RDME
         else if ((strcmp(option, "-gpn") == 0 || strcmp(option, "--gpus-per-node") == 0) && i < (argc-1))
         {
             cudaDevicesPerNode=parseIntReciprocalArg(argv[++i]);
         }

#endif
        //See if the user is trying to turn off cuda capability printing.
         else if ((strcmp(option, "-nr") == 0 || strcmp(option, "--no-reserve-core") == 0))
         {
        	 shouldReserveOutputCore = false;
         }

         else if ((strcmp(option, "-p") == 0 || strcmp(option, "--parameter") == 0))
         {
			cmdline_parameters.push_back(argv[++i]);
         }

        //This must be an invalid option.
        else {
            throw lm::CommandLineArgumentException(option);
        }
    }
}

void parseIntListArg(vector<int> & list, char * arg)
{
    list.clear();
    char * argbuf = new char[strlen(arg)+1];
    strcpy(argbuf,arg);
    char * pch = strtok(argbuf," ,;:\"");
    while (pch != NULL)
    {
        char * rangeDelimiter;
        if ((rangeDelimiter=strstr(pch,"-")) != NULL)
        {
            *rangeDelimiter='\0';
            int begin=atoi(pch);
            int end=atoi(rangeDelimiter+1);
            for (int i=begin; i<=end; i++)
                list.push_back(i);
        }
        else
        {
            if (strlen(pch) > 0) list.push_back(atoi(pch));
        }
        pch = strtok(NULL," ,;:");
    }
    delete[] argbuf;
}

time_t parseTimeArg(char * arg)
{
    char * argbuf = new char[strlen(arg)+1];
    strcpy(argbuf,arg);
    char * pch = strtok(argbuf,":");

    // Parse the arguments into tokens.
    int tokenNumber=0;
    int tokens[3];
    while (pch != NULL)
    {
        if (tokenNumber < 3)
        {
            tokens[tokenNumber++] = atoi(pch);
        }
        else
        {
            delete[] argbuf;
            throw lm::CommandLineArgumentException(arg);
        }
        pch = strtok(NULL,":");
    }
    delete[] argbuf;

    // Calculate the time from the tokens.
    time_t time=0;
    if (tokenNumber == 1)
        time = tokens[0];
    else if (tokenNumber == 2)
        time = (60*tokens[0])+tokens[1];
    else if (tokenNumber == 3)
        time = (3600*tokens[0])+(60*tokens[1])+tokens[2];

    return time;
}

float parseIntReciprocalArg(char * arg)
{
    if (strlen(arg) >= 3 && arg[0] == '1' && arg[1] == '/')
    {
        return 1.0f/(float)atoi(arg+2);
    }
    else
    {
        return (float)atoi(arg);
    }
}

/**
 * Prints the usage for the program.
 */
void printUsage(int argc, char** argv)
{
#ifndef OPT_MPI
	std::cout << "Usage: lm (-h|--help)" << std::endl;
	std::cout << "Usage: lm (-v|--version)" << std::endl;
	std::cout << "Usage: lm [OPTIONS] (-l|--list-devices)" << std::endl;
	std::cout << "Usage: lm [OPTIONS]" << std::endl;
	std::cout << "Usage: lm [OPTIONS] (-s|--script) script_filename [(-sa|--script-args) script_arguments+]" << std::endl;
    std::cout << "Usage: lm [OPTIONS] [SIM_OPTIONS] (-f|--file) simulation_filename " << std::endl;
#else
    std::cout << "Usage: mpirun lm (-h|--help)" << std::endl;
    std::cout << "Usage: mpirun lm (-v|--version)" << std::endl;
    std::cout << "Usage: mpirun lm (-l|--list-devices)" << std::endl;
    std::cout << "Usage: mpirun lm [OPTIONS] [SIM_OPTIONS] (-f|--file) simulation_filename" << std::endl;
#endif
    std::cout << std::endl;
    std::cout << "OPTIONS" << std::endl;
    std::cout << "  -c num_cpus       --cpu=num_cpus               The number of CPUs on which to execute (default all)." << std::endl;
    std::cout << "  -V                --verbose=level              Specified verbosity level. High numbers generate more noise. (default 6)" << std::endl;
    std::cout << "  -cr num           --cpus-per-replicate=num     The number of CPUs (possibly fractional) to assign per replicate, e.g. \"2\", \"1/4\" (default 1)." << std::endl;
#ifdef OPT_CUDA
    std::cout << "  -g cuda_devices   --gpu=cuda_devices           A list of cuda devices on which to execute, e.g. \"0-3\", \"0,2\" (default 0)." << std::endl;
    std::cout << "  -gr num           --gpus-per-replicate=num     The number of cuda devices (possibly fractional) to assign per replicate, e.g. \"2\", \"1/4\" (default 1)." << std::endl;
    std::cout << "  -nc               --no-capabilities            Don't print the capabilities of the CUDA devices." << std::endl;
#endif
    std::cout << "  -nr               --no-reserve-core            Don't reserve a CPU core for the output thread." << std::endl;
    std::cout << std::endl;
    std::cout << "SIM_OPTIONS" << std::endl;
    std::cout << "  -r replicates     --replicates=replicates      A list of replicates to run, e.g. \"0-9\", \"0,11,21\" (default 0)." << std::endl;
    std::cout << "  -sp               --spatially-resolved         The simulations should use the spatially resolved reaction model (default)." << std::endl;
    std::cout << "  -ws               --well-stirred               The simulations should use the well-stirred reaction model." << std::endl;
    std::cout << "  -sl solver        --solver=solver              The specific solver class to use for the simulations." << std::endl;
    std::cout << "  -ck               --checkpoint=interval        Enable checkpointing with the given interval as hh:mm:ss (default 00:00:00 -- disabled)." << std::endl;
    std::cout << "  -p name=value     --parameter name=value       Set solver parameter 'name' to 'value'" << std::endl;
}

