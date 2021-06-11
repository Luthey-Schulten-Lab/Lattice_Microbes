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
 * Author(s): Elijah Roberts, Tyler M Earnest
 */

#ifndef LM_DEBUG_H_
#define LM_DEBUG_H_

#include <string>

namespace lm {

/// @class Print
/// @brief Print messages to the console at varying levels of verbosity.
class Print
{
    static int _verbosityLevel;
public:
    // "Severity"
    static const int VERBOSE_DEBUG              = 10;
    static const int DEBUG                      =  9;
    static const int INFO                       =  4;
    static const int WARNING                    =  3;
    static const int ERROR                      =  2;
    static const int FATAL                      =  1;

    /// @brief Prints to the console at varying levels of verbosity or "Severity"
    /// @param level The level of "Severity" as indicated above
    /// @param fmt A printf style format string
    /// @param Arguments to the printf format string
    static void printf(int level, const char * fmt, ...);
    
    /// @brief Print the date and time to the console
    static void printDateTimeString();
    
    /// @brief Gets the date and time as a string
    /// @return Date and time as a string
    static std::string getDateTimeString();

    /// @brief Gets the verbosity level. Lower numbers are higher priority.
    /// @return Verbosity level
    static int verbosityLevel() {return _verbosityLevel; }

    /// @brief Sets the verbosity level. Lower numbers are higher priority.
    /// @param x Verbosity level
    static void verbosityLevel(int x) {_verbosityLevel=x; }

    /// @brief True if a call to `Print::printf(x, ...)` would print to stdout
    /// @param x Verbosity level
    /// @return Predicate
    static bool ifPrinted(int x) {return  (x <= _verbosityLevel); }
};

}

#endif

