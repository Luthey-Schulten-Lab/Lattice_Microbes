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

#include <iostream>
#include <string>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include "core/Exceptions.h"
#include "io/SimulationFile.h"
#include "lptf/Profile.h"
#include "core/util.h"

void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);

using std::string;
using std::map;
using lm::io::hdf5::SimulationFile;

/**
 * The function being performed.
 */
string function = "";

/**
 * The file to modify.
 */
string filename = "";

/**
 * The number of species.
 */
uint numberSpecies=0;

/**
 * The parameters to set in the file.
 */
map<string,string> parameters;

// Allocate the profile space.
PROF_ALLOC;

int main(int argc, char** argv)
{	
    PROF_INIT;
	try
	{
		printCopyright(argc, argv);
		parseArguments(argc, argv);
		
		if (function == "help")
		{
			printUsage(argc, argv);
		}
		else if (function == "version")
		{
            printBuildConfig();
		}
		else if (function == "set")
		{
		    // If the file doesn't exist, create it.
		    struct stat fileStats;
		    if (stat(filename.c_str(), &fileStats) != 0)
		    {
		        if (numberSpecies > 0)
		            SimulationFile::create(filename, numberSpecies);
		        else
		            SimulationFile::create(filename);
		    }

			// Open the file.
		    SimulationFile file(filename);

		    // Set the parameters.
            printf("Setting parameters in simulation file %s:\n", filename.c_str());
		    for (std::map<std::string,string>::iterator it=parameters.begin(); it != parameters.end(); it++)
		    {
		        printf("%s=%s\n", it->first.c_str(), it->second.c_str());
		        file.setParameter(it->first, it->second);
		    }

		    // Close the file.
		    file.close();
		    printf("Done.\n");
		}
		else
		{
			throw lm::CommandLineArgumentException("unknown function.");
		}
		return 0;
	}
    catch (lm::CommandLineArgumentException e)
    {
    	std::cerr << "Invalid command line argument: " << e.what() << std::endl << std::endl;
        printUsage(argc, argv);
    }
    catch (lm::Exception e)
    {
    	std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (std::exception e)
    {
    	std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (...)
    {
    	std::cerr << "Unknown Exception during execution." << std::endl;
    }
    return -1;
}


/**
 * Parses the command line arguments.
 */
void parseArguments(int argc, char** argv)
{
    //Parse any arguments.
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;
        
        //See if the user is trying to get help.
        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
        	function = "help";
        	break;
        }
        
        //See if the user is trying to get the version info.
        else if (strcmp(option, "-v") == 0 || strcmp(option, "--version") == 0) {
        	function = "version";
        	break;
        }
            
        //See if the user is trying to specify a filename.
        else if (i == 1)
        {
        	function = "set";
    		filename = option;
        }
        
        //See if the user is trying to specify a species count.
        else if (i == 2 && strstr(option, "=") == NULL)
        {
            numberSpecies = atoi(option);
        }

        else if (strstr(option, "=") != NULL)
        {
            char * separator=strstr(option, "=");
            if (separator != NULL && separator > option && strlen(separator) > 1)
            {
                string key(option, separator-option);
                string value(separator+1);
                parameters[key] = value;
            }
        }
             
        //This must be an invalid option.
        else {
            throw lm::CommandLineArgumentException(option);
        }
    }
}

/**
 * Prints the usage for the program.
 */
void printUsage(int argc, char** argv)
{
	std::cout << "Usage: " << argv[0] << " (-h|--help)" << std::endl;
	std::cout << "Usage: " << argv[0] << " (-v|--version)" << std::endl;
	std::cout << "Usage: " << argv[0] << " filename [number_species] (key=value)+" << std::endl;
	std::cout << std::endl;
}
