/*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
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
#include <list>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include "core/Exceptions.h"
#include "ReactionModel.pb.h"
#include "io/SimulationFile.h"
#include "lptf/Profile.h"
#include "core/util.h"

void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);

using std::list;
using std::pair;
using std::string;
using std::vector;
using lm::io::ReactionModel;
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
 * The parameters to set in the file.
 */
list<pair<string,string> > parameters;

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
            bool newFile=false;
		    struct stat fileStats;
		    if (stat(filename.c_str(), &fileStats) != 0)
		    {
                SimulationFile::create(filename);
                newFile = true;
		    }

			// Open the file.
		    SimulationFile file(filename);

		    // Read the reaction model.
            ReactionModel model;
		    if (!newFile) file.getReactionModel(&model);

		    // Set the parameters.
            printf("Setting reaction model in simulation file %s:\n", filename.c_str());
		    for (list<pair<string,string> >::iterator it=parameters.begin(); it != parameters.end(); it++)
		    {
		        string key=it->first;
                string value=it->second;

                // See if the user specified a range.
                vector<int> matrixDims;
		        string range("");
		        size_t rangeStart=key.find_first_of("(");
		        if (rangeStart != string::npos)
		        {
		            range = key.substr(rangeStart);
	                key = key.substr(0, rangeStart);
		        }

                // Set the parameter.
                if (key == "numberSpecies")
                {
                	vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for the number of species.");
                    model.set_number_species((uint)lround(values[0]));
                    printf("%s=%d\n", key.c_str(), model.number_species());
                }
                else if (key == "numberReactions")
                {
                	vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for the number of reactions.");
                    model.set_number_reactions((uint)lround(values[0]));
                    for (uint i=model.reaction_size(); i<model.number_reactions(); i++)
                        model.add_reaction();
                    printf("%s=%d\n", key.c_str(), model.number_reactions());
                }
                else if (key == "InitialSpeciesCounts")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.initial_species_count_size(); i<model.number_species(); i++)
                        model.add_initial_species_count(0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_species());

    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (indices.size() >= 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_initial_species_count(indices[i],(uint)lround(values[0]));
                        }
                        printf("%s(%d...%d)=%d\n", key.c_str(), indices[0], indices[indices.size()-1], model.initial_species_count(indices[0]));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_initial_species_count(indices[i],(uint)lround(values[i]));
                            printf("%s(%d)=%d\n", key.c_str(), indices[i], model.initial_species_count(indices[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "ReactionTypes")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.reaction_size(); i<model.number_reactions(); i++)
                    	model.add_reaction();

                    // Parse the indices.
                	matrixDims.push_back(model.number_reactions());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (indices.size() >= 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.mutable_reaction(indices[i])->set_type((uint)lround(values[0]));
                        }
                        printf("%s(%d...%d)=%d\n", key.c_str(), indices[0], indices[indices.size()-1], model.reaction(indices[0]).type());
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.mutable_reaction(indices[i])->set_type((uint)lround(values[i]));
                            printf("%s(%d)=%d\n", key.c_str(), indices[i], model.reaction(indices[i]).type());
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "ReactionRateConstants")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.reaction_size(); i<model.number_reactions(); i++)
                    	model.add_reaction();

                    // Parse the indices.
                	matrixDims.push_back(model.number_reactions());
                	matrixDims.push_back(10);
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (indices.size() >= 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                        	uint rindex=indices[i]/10;
                        	uint kindex=indices[i]%10;
                        	for (uint j=model.reaction(rindex).rate_constant_size(); j<=kindex; j++)
                        		model.mutable_reaction(rindex)->add_rate_constant(0.0);
                            model.mutable_reaction(rindex)->set_rate_constant(kindex, values[0]);
                        }
                        printf("%s(%d...%d)=%e\n", key.c_str(), indices[0], indices[indices.size()-1], model.reaction(indices[0]/10).rate_constant(indices[0]%10));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                        	uint rindex=indices[i]/10;
                        	uint kindex=indices[i]%10;
                        	for (uint j=model.reaction(rindex).rate_constant_size(); j<=kindex; j++)
                        		model.mutable_reaction(rindex)->add_rate_constant(0.0);
                        	model.mutable_reaction(rindex)->set_rate_constant(kindex, values[i]);
                            printf("%s(%d)=%e\n", key.c_str(), indices[i], model.reaction(indices[i]/10).rate_constant(indices[i]%10));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "StoichiometricMatrix")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.stoichiometric_matrix_size(); i<model.number_species()*model.number_reactions(); i++)
                        model.add_stoichiometric_matrix(0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_species());
                	matrixDims.push_back(model.number_reactions());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_stoichiometric_matrix(indices[i], (int)lround(values[0]));
                        }
                        printf("%s(%d...%d)=%d\n", key.c_str(), indices[0], indices[indices.size()-1], model.stoichiometric_matrix(indices[0]));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_stoichiometric_matrix(indices[i], (int)lround(values[i]));
                            printf("%s(%d)=%d\n", key.c_str(), indices[i], model.stoichiometric_matrix(indices[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "DependencyMatrix")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.dependency_matrix_size(); i<model.number_species()*model.number_reactions(); i++)
                        model.add_dependency_matrix(0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_species());
                	matrixDims.push_back(model.number_reactions());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_dependency_matrix(indices[i], (uint)lround(values[0]));
                        }
                        printf("%s(%d...%d)=%d\n", key.c_str(), indices[0], indices[indices.size()-1], model.dependency_matrix(indices[0]));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_dependency_matrix(indices[i], (uint)lround(values[i]));
                            printf("%s(%d)=%d\n", key.c_str(), indices[i], model.dependency_matrix(indices[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "ReactionRateNoise")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.reaction_size(); i<model.number_reactions(); i++)
                    	model.add_reaction();

                    // Parse the indices.
                	matrixDims.push_back(model.number_reactions());
                	matrixDims.push_back(2);
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (indices.size() >= 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                        	uint rindex=indices[i]/2;
                        	uint kindex=indices[i]%2;
                            model.mutable_reaction(rindex)->set_rate_has_noise(true);
                            if (kindex == 0)
                            	model.mutable_reaction(rindex)->set_rate_noise_variance(values[0]);
                            else
                            	model.mutable_reaction(rindex)->set_rate_noise_tau(values[0]);
                        }
                        if (indices[0]%2 == 0)
                        	printf("%s(%d...%d)=%e\n", key.c_str(), indices[0], indices[indices.size()-1], model.reaction(indices[0]/2).rate_noise_variance());
                        else
                        	printf("%s(%d...%d)=%e\n", key.c_str(), indices[0], indices[indices.size()-1], model.reaction(indices[0]/2).rate_noise_tau());
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                        	uint rindex=indices[i]/2;
                        	uint kindex=indices[i]%2;
                            model.mutable_reaction(rindex)->set_rate_has_noise(true);
                            if (kindex == 0)
                            	model.mutable_reaction(rindex)->set_rate_noise_variance(values[i]);
                            else
                            	model.mutable_reaction(rindex)->set_rate_noise_tau(values[i]);
                            if (indices[i]%2 == 0)
                            	printf("%s(%d)=%e\n", key.c_str(), indices[i], model.reaction(indices[i]/2).rate_noise_variance());
                            else
                            	printf("%s(%d)=%e\n", key.c_str(), indices[i], model.reaction(indices[i]/2).rate_noise_tau());
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else
                {
                    printf("%s: unknown key\n",key.c_str());
                }
		    }

		    // Set the model.
		    file.setReactionModel(&model);

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
    catch (lm::CommandLineArgumentException & e)
    {
    	std::cerr << "Invalid command line argument: " << e.what() << std::endl << std::endl;
        printUsage(argc, argv);
    }
    catch (lm::Exception & e)
    {
    	std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (std::exception & e)
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
        
        //See if the user is trying to specify a key value pair.
        else if (strstr(option, "=") != NULL)
        {
            char * separator=strstr(option, "=");
            if (separator != NULL && separator > option && strlen(separator) > 1)
            {
                string key(option, separator-option);
                string value(separator+1);
                parameters.push_back(pair<string,string>(key,value));
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
	std::cout << "Usage: " << argv[0] << " filename (key=value)+" << std::endl;
	std::cout << std::endl;
}
