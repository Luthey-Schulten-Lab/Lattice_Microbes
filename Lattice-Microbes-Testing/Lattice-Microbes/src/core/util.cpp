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

#include <list>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include "config.h"
#include "core/Exceptions.h"
#include "core/util.h"

using std::list;
using std::pair;
using std::string;
using std::vector;

vector<int> parseIndices(string s, vector<int> matrixDims)
{
	vector<int> indices;

	if (s.size() == 0)
	{
		if (matrixDims.size() == 1)
		{
			for (int i=0; i<matrixDims[0]; i++)
				indices.push_back(i);
		}
		else if (matrixDims.size() == 2)
		{
			for (int i=0; i<matrixDims[0]; i++)
				for (int j=0; j<matrixDims[1]; j++)
					indices.push_back(i*matrixDims[1]+j);
		}
		else if (matrixDims.size() == 3)
		{
			for (int i=0; i<matrixDims[0]; i++)
				for (int j=0; j<matrixDims[1]; j++)
					for (int k=0; k<matrixDims[2]; k++)
						indices.push_back(i*matrixDims[1]*matrixDims[2]+j*matrixDims[2]+k);
		}
		else
			throw lm::Exception("Unsupported matrix dimension", matrixDims.size());
	}
	else
	{
		char * sbuffer = new char[s.size()+1];
		strcpy(sbuffer, s.c_str());
		int sbufferStart=0;

		// Remove any leading/trailing parentheses.
		if (sbuffer[0] == '(') sbufferStart++;
		if (sbuffer[s.size()-1] == ')') sbuffer[s.size()-1] = '\0';

		if (matrixDims.size() == 1)
		{
			pair<int,int> range = parseRange(&(sbuffer[sbufferStart]), matrixDims[0]);
			for (int i=range.first; i<=range.second; i++)
				indices.push_back(i);
		}
		else if (matrixDims.size() == 2)
		{
			char * istr = &(sbuffer[sbufferStart]);
			char * jstr = (char *)"0";
			char * split = strstr(istr, ",");
			if (split != NULL)
			{
				*split = '\0';
				jstr = split+1;
			}
			pair<int,int> irange = parseRange(istr, matrixDims[0]);
			pair<int,int> jrange = parseRange(jstr, matrixDims[1]);
			for (int i=irange.first; i<=irange.second; i++)
				for (int j=jrange.first; j<=jrange.second; j++)
					indices.push_back(i*matrixDims[1]+j);
		}
		else if (matrixDims.size() == 3)
		{
			char * istr = &(sbuffer[sbufferStart]);
			char * jstr = (char *)"0";
			char * kstr = (char *)"0";
			char * split = strstr(istr, ",");
			if (split != NULL)
			{
				*split = '\0';
				jstr = split+1;
				split = strstr(jstr, ",");
				if (split != NULL)
				{
					*split = '\0';
					kstr = split+1;
				}
			}
			pair<int,int> irange = parseRange(istr, matrixDims[0]);
			pair<int,int> jrange = parseRange(jstr, matrixDims[1]);
			pair<int,int> krange = parseRange(kstr, matrixDims[2]);
			for (int i=irange.first; i<=irange.second; i++)
				for (int j=jrange.first; j<=jrange.second; j++)
					for (int k=krange.first; k<=krange.second; k++)
						indices.push_back(i*matrixDims[1]*matrixDims[2]+j*matrixDims[2]+k);
		}
		else
			throw lm::Exception("Unsupported matrix dimension", matrixDims.size());

		delete [] sbuffer;
	}

	return indices;
}

pair<int,int> parseRange(char *rangeStr, int maxDim)
{
	if (strlen(rangeStr) == 1 && rangeStr[0] == ':')
		return pair<int,int>(0,maxDim-1);

	char * split = strstr(rangeStr, ":");
	if (split > rangeStr && split < rangeStr+strlen(rangeStr)-1)
	{
		char * start = rangeStr;
		*split = '\0';
		char * end = split+1;
		return pair<int,int>(atoi(start),atoi(end));
	}
	else
	{
		return pair<int,int>(atoi(rangeStr),atoi(rangeStr));
	}
}

vector<double> parseValues(string s)
{
    char * valueToks = new char[s.size()+1];
    strcpy(valueToks, s.c_str());
    int valueToksStart=0;

    // Remove any leading/trailing brackets.
    if (valueToks[0] == '[') valueToksStart++;
    if (valueToks[s.size()-1] == ']') valueToks[s.size()-1] = '\0';

    // Parse the values.
    vector<double> values;
    char * value = strtok(valueToks+valueToksStart, ";,");
    while (value != NULL)
    {
        double d = atof(value);
        values.push_back(d);
        value = strtok(NULL, ";,");
    }
    delete [] valueToks;
    return values;
}

/**
 * This function prints the copyright notice.
 */
void printCopyright(int argc, char** argv) {

	std::cout << argv[0] << " v" << VERSION_NUM_MAJOR << "." << VERSION_NUM_MINOR << std::endl;
    std::cout << "Copyright (C) " << COPYRIGHT_DATE << " Luthey-Schulten Group," << std::endl;
	std::cout << "University of Illinois at Urbana-Champaign." << std::endl;
	std::cout << std::endl;
}

void printBuildConfig() {
    std::cout << BUILD_CONFIGURATION << std::endl;
}
