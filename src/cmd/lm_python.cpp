/*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.uiuc.edu/~schulten
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
 * Author(s): Elijah Roberts, Tyler Earnest
 */

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <google/protobuf/stubs/common.h>
#include <Python.h>
#include "core/Exceptions.h"
#include "core/Print.h"
#include "lptf/Profile.h"
#include "core/util.h"
#include "lm_module_pack.h"

using std::string;
using std::wstring;
using std::list;
using lm::Print;

wstring mkwStr(const char *str) {
    size_t n = std::strlen(str);
    wchar_t* buf = new wchar_t[n+1];
    std::mbstowcs(buf,str,n);
    buf[n] = 0;
    wstring wstr(buf);
    delete [] buf;
    return wstr;
}

#define STRING std::wstring
#define CHAR wchar_t
#define TO_PY_STR(x) mkwStr(x)
#define STRLITERAL(x) L##x
#define STRCPY(x,y) wcscpy(x,y)

#if PY_MAJOR_VERSION <3
#error "Python 2 not implemented"
#endif

/**
 * The function being performed.
 */
string function = "";

/**
 * The script filename being executed, if applicable.
 */
string scriptFilename = "";

/**
 * The arguments for the script, if applicable.
 */
list<string> scriptArguments;

/**
 * The directory containing the supporting files.
 */
STRING libDir_forPy;
string libDir;

/**
 * The directory containing the supporting files.
 */
STRING userLibDir;

void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);
void discoverEnvironment();
void initPython();
void finalizePython();
void startInterpreter();
void executeScript(string filename, list<string> arguments);
extern "C" PyObject *PyInit__lm(void);

// Allocate the profile space.
PROF_ALLOC;

int main(int argc, char** argv)
{
  setlocale(LC_ALL, "en_US.UTF-8");
    // Make sure we are using the correct protocol buffers library.
    GOOGLE_PROTOBUF_VERIFY_VERSION;

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
        else if (function == "script")
        {
            //discoverEnvironment();
            executeScript(scriptFilename, scriptArguments);
        }
        else if (function == "interpreter")
        {
            //discoverEnvironment();
            startInterpreter();
        }
        else
        {
            throw lm::CommandLineArgumentException("unknown function.");
        }
        return 0;
    }
    catch (lm::CommandLineArgumentException e)
    {
        Print::printf(Print::FATAL, "Invalid command line argument: %s\n", e.what());
        printUsage(argc, argv);
    }
    catch (lm::Exception e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s", e.what());
    }
    catch (std::exception e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s", e.what());
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown exception during execution.");
    }
    return -1;
}

/**
 * Parses the command line arguments.
 */
void parseArguments(int argc, char** argv)
{
    // Set any default options.
    function = "interpreter";
    scriptFilename = "";
    bool parsingScriptArgs = false;
    scriptArguments.clear();

    //Parse any arguments.
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;

        //See if the user is trying to get help.
        if ((!parsingScriptArgs && strcmp(option, "-h") == 0) || (strcmp(option, "--help") == 0)) {
            function = "help";
            break;
        }

        //See if the user is trying to get the version info.
        else if ((!parsingScriptArgs && strcmp(option, "-v") == 0) || (strcmp(option, "--version") == 0)) {
            function = "version";
            break;
        }

        //See if the user is trying to execute a script.
        else if ((!parsingScriptArgs && strcmp(option, "-s") == 0) || (strcmp(option, "--script") == 0))
        {
            function = "script";

            // Get the filename.
            if (i < argc-1)
                scriptFilename = argv[++i];
            else
                throw lm::CommandLineArgumentException("missing script filename.");
        }
        else if ((((!parsingScriptArgs && function == "script") && (strcmp(option, "-sa") == 0)) || (strcmp(option, "--script-args") == 0)))
        {
            parsingScriptArgs = true;
        }
        else if (parsingScriptArgs)
        {
            scriptArguments.push_back(option);
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
    std::cout << "Usage: " << argv[0] << " [OPTIONS]" << std::endl;
    std::cout << "Usage: " << argv[0] << " [OPTIONS] (-s|--script) script_filename [(-sa|--script-args) script_arguments+]" << std::endl;
    std::cout << std::endl;
}


/**
 * Figures out the environment for the program (directories, files, etc).
 */
void discoverEnvironment()
{
    char* env;
    struct stat fileStats;

    // See if we have a lib directory as an environment variable.
    if ((env=getenv("LMLIBDIR")) != NULL && stat((env+string("/lm.py")).c_str(), &fileStats) == 0 && S_ISREG(fileStats.st_mode))
    {
        libDir_forPy = TO_PY_STR(env);
        libDir = env;
    }

    // Otherwise, see if we can find the lib directory based a a few guesses.
    else if (stat("/usr/local/lib/lm/lm.py", &fileStats) == 0 && S_ISREG(fileStats.st_mode))
    {
        libDir_forPy = STRLITERAL("/usr/local/lib/lm/");
        libDir = "/usr/local/lib/lm/";
    }
    else if (stat("/usr/lib/lm/lm.py", &fileStats) == 0 && S_ISREG(fileStats.st_mode))
    {
        libDir_forPy = STRLITERAL("/usr/lib/lm/");
        libDir = "/usr/lib/lm/";
    }

    // Without a lib directory we can't continue.
    else
    {
        throw lm::Exception("Could not find installation lib directory, please set the LMLIBDIR environment variable appropriately.");
    }

    // See if we have any user paths as an environment variable.
    if ((env=getenv("LMPATH")) != NULL)
    {
        userLibDir = TO_PY_STR(env);
    }
}

/**
 * Starts the python interpreter.
 */

void initPython()
{
	PyImport_AppendInittab("_lm", PyInit__lm);
	Py_Initialize();

	PyObject *builtins = PyEval_GetBuiltins();
	PyObject *compile = PyDict_GetItemString(builtins, "compile");
	PyObject *code = PyObject_CallFunction(compile, "sss", lm_swig_wrapper, "lm.py", "exec");
	PyObject *module = PyImport_ExecCodeModule("lm", code);
	PyRun_SimpleString("import lm");
	PyRun_SimpleString("import readline");
}

void finalizePython()
{
	Py_Finalize();	
}

void startInterpreter()
{
	initPython();
	
	// Start the interpreter.
	if (PyRun_InteractiveLoop(stdout, "Lattice Microbe"))
		throw lm::Exception("Python failed");
	
	finalizePython();
}

void executeScript(string filename, std::list<string> arguments)
{
	initPython();
	Print::printf(Print::DEBUG, "Init python.");
	
	// Allocate a buffer for the arguments.
	int argc = arguments.size()+1;
	CHAR** argv = new CHAR*[argc];
	if (argv == NULL)
		throw lm::Exception("Failed to allocate memory for the script arguments.");
	
	// Copy the script name into the first argument.
	int argIndex=0;
	argv[argIndex] = new CHAR[filename.size()+1];
	if (argv[argIndex] == NULL) throw lm::Exception("Failed to allocate memory for the script arguments.");
	STRCPY(argv[argIndex], TO_PY_STR(filename.c_str()).c_str());
	argIndex++;
		
	// Copy the rest of the arguments into the buffer.
	for (std::list<string>::iterator it=arguments.begin(); it != arguments.end(); it++, argIndex++)
	{
		string arg = *it;
		argv[argIndex] = new CHAR[arg.size()+1];
		if (argv[argIndex] == NULL) throw lm::Exception("Failed to allocate memory for the script arguments.");
		STRCPY(argv[argIndex], TO_PY_STR(arg.c_str()).c_str());
	}
    Print::printf(Print::DEBUG, "Copied args.");
		
	// Set the arguments.
	PySys_SetArgv(argc, argv);
	Print::printf(Print::DEBUG, "Set args.");
	
	// Open the script.
	std::ifstream scriptFile(filename.c_str());
	scriptFile.exceptions(std::ifstream::eofbit|std::ifstream::failbit|std::ifstream::badbit);
	Print::printf(Print::DEBUG, "Opened script.");
	
	// Figure out how long the script is.
	scriptFile.seekg(0, std::ios::end);
	std::streampos length = scriptFile.tellg();
	scriptFile.seekg(0, std::ios::beg);
	Print::printf(Print::DEBUG, "Got script length %d.", (int)length);
	
	// Allocate a buffer for the script data.
	char * scriptData = new char[(int)length+1];
	memset(scriptData, 0, (int)length+1);
	if (scriptData == NULL)
	{
		scriptFile.close();
		throw lm::Exception("Could not allocate script buffer.");
	}
	Print::printf(Print::DEBUG, "Alloctaed buffer.");
	
	// Read the script data.
	scriptFile.read(scriptData, length);
	Print::printf(Print::DEBUG, "Read script.");
	
	// Close the script.
	scriptFile.close();
	Print::printf(Print::DEBUG, "Closed script.");
	
	/*// Set the replicate as a variable.
	PyObject* module = PyImport_AddModule("__main__");
	PyObject* dictionary = PyModule_GetDict(module);
	PyObject* pyReplicate = PyInt_FromLong(replicate);
	PyDict_SetItemString(dictionary, "replicate", pyReplicate);
	Py_DECREF(pyReplicate);
	*/

	// Run the script.
	if (PyRun_SimpleString(scriptData))
		throw lm::Exception("Failed to run script.");
	Print::printf(Print::DEBUG, "Finished running script.");
	
	// Free the script data buffer.
	delete[] scriptData;
	scriptData = NULL;
	
	// Stop the python environment.
	finalizePython();
	
	// Free the command line argument storage.
	if (argv != NULL)
	{
		for (int i=0; i<argc; i++)
		{
			if (argv[i] != NULL)
			{
				delete[] argv[i];
				argv[i] = NULL;
			}
		}
		delete[] argv;
		argv = NULL;
		argc = 0;
	}
}
