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

#ifndef LM_EXCEPTIONS_H_
#define LM_EXCEPTIONS_H_

#include <cstdio>
#include <exception>

namespace lm
{

/// @class Exception
/// @brief Base class for exceptions.
///
/// A class for writing exceptions to a buffer that can then be written
/// to either a console or a stream.
class Exception : public std::exception
{
protected:
    static const int MAX_MESSAGE_SIZE = 1025;
    char messageBuffer[MAX_MESSAGE_SIZE];
    
public:
    /// @brief Create an Exception
	Exception(const char * message="")                                                      {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s", message);}
    /// @brief Create and Exception with one integer error code
	Exception(const char * message, const int arg)                                          {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %d", message, arg);}
    /// @brief Create and Exception with two integer error codes
    Exception(const char * message, const int arg1,    const int arg2)                      {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %d, %d", message, arg1, arg2);}
    /// @brief Create and Exception with three integer error codes
    Exception(const char * message, const int arg1,    const int arg2,    const int arg3)   {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %d, %d, %d", message, arg1, arg2, arg3);}
    /// @brief Create and Exception with one error string
	Exception(const char * message, const char * arg)                                       {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %s", message, arg);}
    /// @brief Create and Exception with two error strings
	Exception(const char * message, const char * arg1, const char* arg2)                    {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %s, %s", message, arg1, arg2);}
    /// @brief Create and Exception with three error strings
    Exception(const char * message, const char * arg1, const char* arg2,  const char* arg3) {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %s, %s, %s", message, arg1, arg2, arg3);}
    /// @brief Create and Exception with one integer error code and one error string
    Exception(const char * message, const char * arg1, const int arg2)                      {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %s, %d", message, arg1, arg2);}
    /// @brief Create and Exception with two integer error codes and one error string
    Exception(const char * message, const char * arg1, const int arg2,    const int arg3)   {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %s, %d, %d", message, arg1, arg2, arg3);}
    /// @brief Create and Exception with one integer error code, a file and a line number 
    Exception(const char * message, const int arg,     const char * file, const int line)   {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %d (%s:%d)", message, arg, file, line);}
    /// @brief Create and Exception with one error string, a file and a line number 
    Exception(const char * message, const char * arg,  const char * file, const int line)   {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %s (%s:%d)", message, arg, file, line);}
    /// @brief Destroy the Exception
	virtual ~Exception() throw() {}
    
    /// @brief Get the error string
    /// @return messageBuffer A pointer to the error message
	virtual const char * what() const throw() {return messageBuffer;}
};

/// @class CommandLineArgumentException
/// @brief Exception from a command line argument.
class CommandLineArgumentException : public Exception
{
public:
	CommandLineArgumentException(const char* message) : Exception(message) {}
    CommandLineArgumentException(const char* message, const char* arg1) : Exception(message, arg1) {}
};

/// @class InvalidArgException
/// @brief Exception for an argument to a function.
class InvalidArgException : public Exception
{
public:
	InvalidArgException(const char* argMessage) : Exception("Invalid argument", argMessage) {}
	InvalidArgException(const char* arg, const char* argMessage) : Exception("Invalid argument", arg, argMessage) {}
    InvalidArgException(const char* arg, const char* argMessage, const char * argMessageParameter) : Exception("Invalid argument", arg, argMessage, argMessageParameter) {}
    InvalidArgException(const char* arg, const char* argMessage, const int argMessageParameter) : Exception("Invalid argument", arg, argMessage, argMessageParameter) {}
    InvalidArgException(const char* arg, const char* argMessage, const int argMessageParameter1, const int argMessageParameter2) : Exception() {snprintf(messageBuffer,MAX_MESSAGE_SIZE,"%s: %s, %s (%d,%d)", "Invalid argument", arg, argMessage, argMessageParameter1, argMessageParameter2);}
};
    
/// @class IOException
/// @brief Exception when input or output fails.
class IOException : public Exception
{
public:
    IOException(const char* message, const char* arg) : Exception("IO exception", message, arg) {}
    IOException(const char* message, const int arg) : Exception("IO exception", message, arg) {}
};

}

#endif
