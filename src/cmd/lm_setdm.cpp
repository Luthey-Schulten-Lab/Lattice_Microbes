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
#include "DiffusionModel.pb.h"
#include "ReactionModel.pb.h"
#include "io/SimulationFile.h"
#include "rng/XORShift.h"
#include "lptf/Profile.h"
#include "core/util.h"

void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);

using std::list;
using std::pair;
using std::string;
using std::vector;
using lm::io::DiffusionModel;
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

/**
 * Whether we should place particles automatically.
 */
bool autoplace = true;

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
            DiffusionModel model;
		    if (!newFile) file.getDiffusionModel(&model);

		    // Set the parameters.
            printf("Setting diffusion model in simulation file %s:\n", filename.c_str());
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
                    printf("%s=%d\n", key.c_str(), model.number_reactions());
                }
                else if (key == "numberSiteTypes")
                {
                	vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for the number of site types.");
                    model.set_number_site_types((uint)lround(values[0]));
                    printf("%s=%d\n", key.c_str(), model.number_site_types());
                }
                else if (key == "latticeSize")
                {
                    vector<double> values = parseValues(value);
                    if (values.size() != 3) throw lm::Exception("Three values must be specified for the lattice size.");
                    model.set_lattice_x_size((uint)lround(values[0]));
                    model.set_lattice_y_size((uint)lround(values[1]));
                    model.set_lattice_z_size((uint)lround(values[2]));
                    printf("%s=[%d,%d,%d]\n", key.c_str(), model.lattice_x_size(), model.lattice_y_size(), model.lattice_z_size());
                }
                else if (key == "latticeSpacing")
                {
                    vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for lattice spacing.");
                    model.set_lattice_spacing(values[0]);
                    printf("%s=%e\n", key.c_str(), model.lattice_spacing());
                }
                else if (key == "particlesPerSite")
                {
                    vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for number of particles per site.");
                    model.set_particles_per_site((uint)lround(values[0]));
                    printf("%s=%d\n", key.c_str(), model.particles_per_site());
                }
                else if (key == "DiffusionMatrix")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.diffusion_matrix_size(); i<model.number_species()*model.number_site_types()*model.number_site_types(); i++)
                        model.add_diffusion_matrix(0.0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_site_types());
                	matrixDims.push_back(model.number_site_types());
                	matrixDims.push_back(model.number_species());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_diffusion_matrix(indices[i],values[0]);
                        }
                        printf("%s(%d...%d)=%e\n", key.c_str(), indices[0], indices[indices.size()-1], model.diffusion_matrix(indices[0]));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_diffusion_matrix(indices[i],values[i]);
                            printf("%s(%d)=%e\n", key.c_str(), indices[i], model.diffusion_matrix(indices[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values.", indices.size(), values.size());
                }
                else if (key == "ReactionLocationMatrix")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.reaction_location_matrix_size(); i<model.number_reactions()*model.number_site_types(); i++)
                        model.add_reaction_location_matrix(0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_reactions());
                	matrixDims.push_back(model.number_site_types());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_reaction_location_matrix(indices[i],lroundf(values[0]));
                        }
                        printf("%s(%d...%d)=%d\n", key.c_str(), indices[0], indices[indices.size()-1], model.reaction_location_matrix(indices[0]));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_reaction_location_matrix(indices[i],lroundf(values[i]));
                            printf("%s(%d)=%d\n", key.c_str(), indices[i], model.reaction_location_matrix(indices[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values.", indices.size(), values.size());
                }
                else if (key == "PlaceParticles")
                {
                	if (value == "0" || value == "false" || value == "False" || value == "FALSE")
                		autoplace = false;
                }
                else
                {
                    printf("%s: unknown key\n",key.c_str());
                }
		    }

		    // Set the model.
		    file.setDiffusionModel(&model);

		    // Create the lattice and lattice sites matrices.
		    uint8_t * data =  new uint8_t[model.lattice_x_size()*model.lattice_y_size()*model.lattice_z_size()*model.particles_per_site()];
		    for (uint i=0; i<model.lattice_x_size()*model.lattice_y_size()*model.lattice_z_size()*model.particles_per_site(); i++)
		        data[i] = 0;
		    uint8_t * sitesData =  new uint8_t[model.lattice_x_size()*model.lattice_y_size()*model.lattice_z_size()];
		    for (uint i=0; i<model.lattice_x_size()*model.lattice_y_size()*model.lattice_z_size(); i++)
		    	sitesData[i] = 0;

		    // If we can read the reaction model, fill in the initial species counts.
		    if (autoplace)
		    {
				ReactionModel reactionModel;
				file.getReactionModel(&reactionModel);
				if (reactionModel.number_species() == model.number_species() && (uint)reactionModel.initial_species_count_size() == model.number_species())
				{
					lm::rng::XORShift rng(0, 0);
					for (uint i=0; i<reactionModel.number_species(); i++)
					{
						uint numberAttempts=0;
						uint placed=0;
						while (placed < reactionModel.initial_species_count(i))
						{
							numberAttempts++;
							double randomValues[3];
							rng.getRandomDoubles(randomValues, 3);
							uint x = (uint)floor(randomValues[0]*(double)model.lattice_x_size());
							if (x == model.lattice_x_size()) x--;
							uint y = (uint)floor(randomValues[1]*(double)model.lattice_y_size());
							if (y == model.lattice_y_size()) y--;
							uint z = (uint)floor(randomValues[2]*(double)model.lattice_z_size());
							if (z == model.lattice_z_size()) z--;

							for (uint j=0; j<model.particles_per_site(); j++)
							{
								uint index=x*model.lattice_y_size()*model.lattice_z_size()*model.particles_per_site() + \
										   y*model.lattice_z_size()*model.particles_per_site() + \
										   z*model.particles_per_site() + j;
								if (data[index] == 0)
								{
									data[index] = i+1;
									placed++;
									break;
								}
							}
						}
						printf("Placed %d particles of type %d in %d attempts.\n",placed,i,numberAttempts);
					}
				}
		    }

		    // Write the lattice to the file.
		    file.setDiffusionModelLattice(&model, data, sitesData);
		    delete[] data;
		    delete[] sitesData;

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
