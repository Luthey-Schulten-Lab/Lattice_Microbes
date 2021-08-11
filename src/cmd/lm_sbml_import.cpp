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
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sbml/SBMLTypes.h>
#include "core/Exceptions.h"
#include "ReactionModel.pb.h"
#include "io/SimulationFile.h"
#include "core/Print.h"
#include "lptf/Profile.h"
#include "core/util.h"

void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);

using std::map;
using std::vector;
using std::string;
using lm::Exception;
using lm::io::ReactionModel;
using lm::io::hdf5::SimulationFile;
using lm::Print;

/**
 * The function being performed.
 */
string function = "";

/**
 * The lm file to output the model into.
 */
string outputFilename = "";

/**
 * The sbml file to import the model from.
 */
string inputFilename = "";

/**
 * Parameters specified by the user.
 */
map<string, double> userParameterValues;

void importSBMLModel(SimulationFile * lmFile, string sbmlFilename);
void importSBMLModelL3V1(ReactionModel * lmModel, Model * sbmlModel);
void importSBMLModelL3V1Kinetics(uint reactionIndex, KineticLaw * kinetics, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices, uint numberReactions, vector<string> & globalParameters, map<string,double> & globalParameterValues);
bool isFirstOrderReaction(const ASTNode * root, vector<string> & parameters);
void importFirstOrderReaction(const ASTNode * root, vector<string> & parameters, map<string,double> & parameterValues, uint reactionIndex, uint numberReactions, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices);
bool isSecondOrderReaction(const ASTNode * root, vector<string> & parameters);
void importSecondOrderReaction(const ASTNode * root, vector<string> & parameters, map<string,double> & parameterValues, uint reactionIndex, uint numberReactions, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices);
bool isSecondOrderSelfReaction(const ASTNode * root, vector<string> & parameters);
void importSecondOrderSelfReaction(const ASTNode * root, vector<string> & parameters, map<string,double> & parameterValues, uint reactionIndex, uint numberReactions, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices);

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
		else if (function == "import")
		{
            // Make sure the input file exists.
            struct stat fileStats;
            if (stat(inputFilename.c_str(), &fileStats) != 0)
            {
                throw Exception("The specified SBML file does not exist", inputFilename.c_str());
            }

		    // If the output file doesn't exist, create it.
		    if (stat(outputFilename.c_str(), &fileStats) != 0)
		    {
		        SimulationFile::create(outputFilename);
		    }

			// Open the file.
		    SimulationFile outputFile(outputFilename);

		    // Open the sbml file.
		    importSBMLModel(&outputFile, inputFilename);

		    // Close the file.
		    outputFile.close();
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
    catch (Exception e)
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

void importSBMLModel(SimulationFile * lmFile, string sbmlFilename)
{
    // Read in the SBML document.
    std::unique_ptr<SBMLDocument> sbmlDocument(readSBML(sbmlFilename.c_str()));
    if (sbmlDocument->getNumErrors() > 0)
    {
        sbmlDocument->printErrors();
        throw Exception("Error reading SBML file");
    }

    // Make sure we know how to process the document.
    if (sbmlDocument->getLevel() == 3 && sbmlDocument->getVersion() == 1)
    {
        // Build the reaction model from the SBML model.
        ReactionModel lmModel;
        Model * sbmlModel = sbmlDocument->getModel();
        importSBMLModelL3V1(&lmModel, sbmlModel);
        lmFile->setReactionModel(&lmModel);
    }
    else
        throw Exception("Unsupported SBML format", sbmlDocument->getLevel(), sbmlDocument->getVersion());
}

void importSBMLModelL3V1(ReactionModel * lmModel, Model * sbmlModel)
{
    vector<string> globalParameters;
    map<string,double> globalParameterValues;

    // Get the units for the model.
    string modelSubstanceUnits = sbmlModel->getSubstanceUnits();

    // Process any global parameters.
    for (uint i=0; i<sbmlModel->getNumParameters(); i++)
    {
        globalParameters.push_back(sbmlModel->getParameter(i)->getId());
        globalParameterValues[sbmlModel->getParameter(i)->getId()] = sbmlModel->getParameter(i)->getValue();
    }

    // Process the compartments.
    for (uint i=0; i<sbmlModel->getNumCompartments(); i++)
    {
    	globalParameters.push_back(sbmlModel->getCompartment(i)->getId());
    	globalParameterValues[sbmlModel->getCompartment(i)->getId()] = sbmlModel->getCompartment(i)->getSize();
    }

    // Process the species.
    map<string,uint> speciesIndices;
    uint numberSpecies = sbmlModel->getNumSpecies();
    lmModel->set_number_species(numberSpecies);
    for (uint i=0; i<numberSpecies; i++)
    {
        Species * species = sbmlModel->getSpecies(i);
        speciesIndices[species->getId()] = i;

        // Get the units for the amount.
        string substanceUnits = modelSubstanceUnits;
        if (species->isSetSubstanceUnits()) substanceUnits = species->getSubstanceUnits();
        if (substanceUnits != "item") throw Exception("Unsupported species substance units", substanceUnits.c_str());

        // Make sure we can process the species.
        if (!species->getHasOnlySubstanceUnits()) throw Exception("Unsupported species property", "hasOnlySubstanceUnits must be true");
        if (species->getBoundaryCondition()) throw Exception("Unsupported species property", "boundaryCondition must be false");
        if (species->getConstant()) throw Exception("Unsupported species property", "constant must be false");
        if (species->isSetConversionFactor()) throw Exception("Unsupported species property", "conversionFactor must not be set");

        // Get the initial count for the species.
        if (species->isSetInitialAmount())
        {
            lmModel->add_initial_species_count((uint)lround(species->getInitialAmount()));
        }
        else if (species->isSetInitialConcentration())
        {
            throw Exception("Unsupported species property", "initialConcentration must not be set");
        }
        else
            throw Exception("Unsupported species property", "initialAmount or initialConcentration must be set");
    }

    // Process the reactions.
    uint numberReactions = sbmlModel->getNumReactions();
    lmModel->set_number_reactions(numberReactions);
    int * S = new int[numberSpecies*numberReactions];
    uint * D = new uint[numberSpecies*numberReactions];
    for (uint i=0; i<numberSpecies*numberReactions; i++)
    {
        S[i] = 0;
        D[i] = 0;
    }
    for (uint i=0; i<numberReactions; i++)
    {
        Reaction * reaction = sbmlModel->getReaction(i);

        // Make sure we can process the reaction.
        if (reaction->getReversible()) throw Exception("Unsupported reaction property", "reversible must be false");
        if (!reaction->isSetKineticLaw()) throw Exception("Unsupported reaction property", "must have a kinetic law");

        // Go through the list of reactants.
        for (uint j=0; j<reaction->getNumReactants(); j++)
        {
            SpeciesReference * reactant = reaction->getReactant(j);
            uint speciesIndex = speciesIndices[reactant->getSpecies()];

            // Make sure we can process the reactant.
            if (!reactant->isSetStoichiometry()) throw Exception("Unsupported reaction property", "stoichiometry for reactants must be set");
            if (!reactant->getConstant()) throw Exception("Unsupported reaction property", "stoichiometry for reactants must be constant");

            // Make the proper entry in the S matrix.
            S[speciesIndex*numberReactions+i] -= reactant->getStoichiometry();
        }

        // Go through the list of products.
        for (uint j=0; j<reaction->getNumProducts(); j++)
        {
            SpeciesReference * product = reaction->getProduct(j);
            uint speciesIndex = speciesIndices[product->getSpecies()];

            // Make sure we can process the product.
            if (!product->isSetStoichiometry()) throw Exception("Unsupported reaction property", "stoichiometry for products must be set");
            if (!product->getConstant()) throw Exception("Unsupported reaction property", "stoichiometry for products must be constant");

            // Make the proper entry in the S matrix.
            S[speciesIndex*numberReactions+i] += product->getStoichiometry();
        }

        // Process the kinetic law.
        lmModel->add_reaction();
        KineticLaw * kinetics = reaction->getKineticLaw();
        importSBMLModelL3V1Kinetics(i, kinetics, lmModel, D, speciesIndices, numberReactions, globalParameters, globalParameterValues);
    }

    // Fill in the S and D matrices.
    for (uint i=0; i<numberSpecies*numberReactions; i++)
    {
        lmModel->add_stoichiometric_matrix(S[i]);
        lmModel->add_dependency_matrix(D[i]);
    }
}

void importSBMLModelL3V1Kinetics(uint reactionIndex, KineticLaw * kinetics, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices, uint numberReactions, vector<string> & globalParameters, map<string,double> & globalParameterValues)
{
    vector<string> localParameters;
    map<string,double> localParameterValues;

	// Bring the global parameters into the local scope.
    for (vector<string>::iterator it = globalParameters.begin(); it != globalParameters.end(); it++)
    {
        localParameters.push_back(*it);
        localParameterValues[*it] = globalParameterValues[*it];
    }

    // Get a list of the local parameters.
    for (uint i=0; i<kinetics->getNumLocalParameters(); i++)
    {
        LocalParameter * localParameter = kinetics->getLocalParameter(i);

        // See if the user specified this parameter.
        /*if ()
        {

        }

        // Otherwise, use the value from the file.
        else
        {*/
            if (!localParameter->isSetValue()) throw Exception("Unsupported reaction property", "value for local parameters must be set");
            localParameters.push_back(localParameter->getId());
            localParameterValues[localParameter->getId()] = localParameter->getValue();
        //}
    }

    // Figure out the reaction type.
    const ASTNode * math = kinetics->getMath();
    if (isFirstOrderReaction(math, localParameters))
        importFirstOrderReaction(math, localParameters, localParameterValues, reactionIndex, numberReactions, lmModel, D, speciesIndices);
    else if (isSecondOrderReaction(math, localParameters))
        importSecondOrderReaction(math, localParameters, localParameterValues, reactionIndex, numberReactions, lmModel, D, speciesIndices);
    else if (isSecondOrderSelfReaction(math, localParameters))
        importSecondOrderSelfReaction(math, localParameters, localParameterValues, reactionIndex, numberReactions, lmModel, D, speciesIndices);
    else
        throw Exception("Unsupported kinetic law", kinetics->getFormula().c_str());
}

void getSpeciesUsedInExpression(vector<string> & speciesUsed, const ASTNode * node, vector<string> & parameters)
{
    if (node->getType() == AST_NAME)
    {
        string name(node->getName());
        bool isParameter = false;
        for (vector<string>::iterator it = parameters.begin(); it != parameters.end(); it++)
        {
            if (*it == name)
            {
                isParameter = true;
                break;
            }
        }
        if (!isParameter)
        {
            speciesUsed.push_back(name);
        }
    }
    for (uint i=0; i<node->getNumChildren(); i++)
    {
        getSpeciesUsedInExpression(speciesUsed, node->getChild(i), parameters);
    }
}

void getOperatorsUsedInExpression(vector<string> & speciesUsed, const ASTNode * node)
{
    if (node->isOperator())
    {
        speciesUsed.push_back(string(1,node->getCharacter()));
    }
    else if (node->isName())
    {
    }
    else if (node->isNumber())
    {
    }
    else
    {
        speciesUsed.push_back(string("?"));
    }
    for (uint i=0; i<node->getNumChildren(); i++)
    {
        getOperatorsUsedInExpression(speciesUsed, node->getChild(i));
    }
}

const ASTNode * getFirstExpressionOfType(const ASTNode * node, ASTNodeType_t type)
{
    if (node->getType() == type) return node;
    for (uint i=0; i<node->getNumChildren(); i++)
    {
        const ASTNode *tmp=getFirstExpressionOfType(node->getChild(i), type);
        if (tmp != NULL) return tmp;
    }
    return NULL;
}


double calculateMultiplierInExpression(const ASTNode * node, map<string,double> & parameterValues, bool ignoreSpeciesMinusOne=false)
{
    if (node->getType() == AST_TIMES)
    {
        double value=1.0;
        for (uint i=0; i<node->getNumChildren(); i++)
        {
            value *= calculateMultiplierInExpression(node->getChild(i), parameterValues, ignoreSpeciesMinusOne);
        }
        return value;
    }
    else if (node->getType() == AST_DIVIDE && node->getNumChildren() == 2)
    {
        return calculateMultiplierInExpression(node->getChild(0), parameterValues, ignoreSpeciesMinusOne)/calculateMultiplierInExpression(node->getChild(1), parameterValues, ignoreSpeciesMinusOne);
    }
    else if (ignoreSpeciesMinusOne && node->getType() == AST_MINUS)
    {

        if (node->getNumChildren() != 2) throw Exception("Unsupported expression 1");
        if (!node->getChild(0)->isName()) throw Exception("Unsupported expression 2");
        if (parameterValues.count(node->getChild(0)->getName()) == 1) throw Exception("Unsupported expression 3");
        if (!node->getChild(1)->isInteger()) throw Exception("Unsupported expression 4");
        if (node->getChild(1)->getInteger() != 1) throw Exception("Unsupported expression 5");
        return 1.0;
    }
    else if (node->getType() == AST_INTEGER)
    {
        return (double)node->getInteger();
    }
    else if (node->getType() == AST_REAL)
    {
        return node->getReal();
    }
    else if (node->getType() == AST_NAME)
    {
        if (parameterValues.count(node->getName()) == 1)
            return parameterValues[node->getName()];
        else
            return 1.0;
    }
    else if (node->getType() == AST_NAME_AVOGADRO)
    {
    	return 6.02214179e23;
    }
    else
        throw Exception("Unsupported ast type", node->getType());
}

bool isFirstOrderReaction(const ASTNode * root, vector<string> & parameters)
{
    // Make sure the rate only depends on only one species.
    vector<string> speciesUsed;
    getSpeciesUsedInExpression(speciesUsed, root, parameters);
    if (speciesUsed.size() != 1) return false;

    // Make sure the expression only involves multiplication and division.
    vector<string> operatorsUsed;
    getOperatorsUsedInExpression(operatorsUsed, root);
    for (vector<string>::iterator it = operatorsUsed.begin(); it != operatorsUsed.end(); it++)
    {
        if (*it != "*" && *it != "/") return false;
    }

    return true;
}

void importFirstOrderReaction(const ASTNode * root, vector<string> & parameters, map<string,double> & parameterValues, uint reactionIndex, uint numberReactions, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices)
{
    // Get the index of the species.
    vector<string> speciesUsed;
    getSpeciesUsedInExpression(speciesUsed, root, parameters);
    uint speciesIndex = speciesIndices[speciesUsed[0]];

    // Set the reaction type.
    lmModel->mutable_reaction(reactionIndex)->set_type(1);

    // Set the reaction dependency.
    D[speciesIndex*numberReactions+reactionIndex] = 1;

    // Get the rate constant.
    double k=calculateMultiplierInExpression(root, parameterValues);
    lmModel->mutable_reaction(reactionIndex)->add_rate_constant(k);
}

bool isSecondOrderReaction(const ASTNode * root, vector<string> & parameters)
{
    // Make sure the rate only depends on only two species.
    vector<string> speciesUsed;
    getSpeciesUsedInExpression(speciesUsed, root, parameters);

    if (speciesUsed.size() != 2) return false;
    if (speciesUsed[0] == speciesUsed[1]) return false;

    // Make sure the expression only involves multiplication and division.
    vector<string> operatorsUsed;
    getOperatorsUsedInExpression(operatorsUsed, root);
    for (vector<string>::iterator it = operatorsUsed.begin(); it != operatorsUsed.end(); it++)
    {
        if (*it != "*" && *it != "/") return false;
    }

    return true;
}

void importSecondOrderReaction(const ASTNode * root, vector<string> & parameters, map<string,double> & parameterValues, uint reactionIndex, uint numberReactions, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices)
{
    // Get the species used.
    vector<string> speciesUsed;
    getSpeciesUsedInExpression(speciesUsed, root, parameters);

    // Set the reaction type.
    lmModel->mutable_reaction(reactionIndex)->set_type(2);

    // Set the reaction dependency.
    D[speciesIndices[speciesUsed[0]]*numberReactions+reactionIndex] = 1;
    D[speciesIndices[speciesUsed[1]]*numberReactions+reactionIndex] = 1;

    // Get the rate constant.
    double k=calculateMultiplierInExpression(root, parameterValues);
    lmModel->mutable_reaction(reactionIndex)->add_rate_constant(k);
}

bool isSecondOrderSelfReaction(const ASTNode * root, vector<string> & parameters)
{
    // Make sure the rate only depends on only two species.
    vector<string> speciesUsed;
    getSpeciesUsedInExpression(speciesUsed, root, parameters);
    if (speciesUsed.size() != 2) return false;
    if (speciesUsed[0] != speciesUsed[1]) return false;

    // Make sure the expression only involves multiplication and division.
    vector<string> operatorsUsed;
    getOperatorsUsedInExpression(operatorsUsed, root);
    uint numMinuses=0;
    for (vector<string>::iterator it = operatorsUsed.begin(); it != operatorsUsed.end(); it++)
    {
        if (*it != "*" && *it != "/" && *it != "-")
        {
        	return false;
        	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to presence of a %s operator.", *it->c_str());
        }
        if (*it == "-") numMinuses++;
    }
    if (numMinuses != 1)
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to presence of %d minus operators.", numMinuses);
    	return false;
    }

    // Make sure one of the species entries is species-1.
    const ASTNode * minusNode=getFirstExpressionOfType(root, AST_MINUS);
    if (minusNode == NULL)
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to lack of minus operator.");
    	return false;
    }
    if (minusNode->getNumChildren() != 2)
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to %d children for the minus node.", minusNode->getNumChildren());
    	return false;
    }
    if (!minusNode->getChild(0)->isName())
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to the first minus child not being a symbol.");
    	return false;
    }
    if (string(minusNode->getChild(0)->getName()) != speciesUsed[0])
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to the first child being an invalid symbol: %s.", minusNode->getChild(0)->getName());
    	return false;
    }
    if (!minusNode->getChild(1)->isInteger())
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to the second child not being an integer.");
    	return false;
    }
    if (minusNode->getChild(1)->getInteger() != 1)
    {
    	Print::printf(Print::VERBOSE_DEBUG, "Reaction was not second order self due to the second child being an invalid integer: %d.", minusNode->getChild(1)->getInteger());
    	return false;
    }

    return true;
}

void importSecondOrderSelfReaction(const ASTNode * root, vector<string> & parameters, map<string,double> & parameterValues, uint reactionIndex, uint numberReactions, ReactionModel * lmModel, uint * D, map<string,uint> & speciesIndices)
{
    // Get the species used.
    vector<string> speciesUsed;
    getSpeciesUsedInExpression(speciesUsed, root, parameters);

    // Set the reaction type.
    lmModel->mutable_reaction(reactionIndex)->set_type(3);

    // Set the reaction dependency.
    D[speciesIndices[speciesUsed[0]]*numberReactions+reactionIndex] = 1;

    // Get the rate constant.
    double k=calculateMultiplierInExpression(root, parameterValues, true);
    lmModel->mutable_reaction(reactionIndex)->add_rate_constant(k);
}

/**
 * Parses the command line arguments.
 */
void parseArguments(int argc, char** argv)
{
    // Parse any arguments.
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;
        
        // See if the user is trying to get help.
        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
        	function = "help";
        	break;
        }
        
        // See if the user is trying to get the version info.
        else if (strcmp(option, "-v") == 0 || strcmp(option, "--version") == 0) {
        	function = "version";
        	break;
        }
            
        // See if the user is trying to specify an output filename.
        else if (i == 1)
        {
        	function = "import";
    		outputFilename = option;
        }
        
        // See if the user is trying to specify an input filename
        else if (i == 2)
        {
            inputFilename = option;
        }

        //See if the user is trying to specify a key value pair.
        else if (strstr(option, "=") != NULL)
        {
            char * separator=strstr(option, "=");
            if (separator != NULL && separator > option && strlen(separator) > 1)
            {
                string key(option, separator-option);
                double value = atof(separator+1);
                userParameterValues[key] = value;
            }
        }

        // This must be an invalid option.
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
	std::cout << "Usage: " << argv[0] << " lm_filename sbml_filename (key=value)+" << std::endl;
	std::cout << std::endl;
}
