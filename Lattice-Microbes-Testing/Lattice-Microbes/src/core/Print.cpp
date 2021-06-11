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
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdarg>
#include <unistd.h>
#include "config.h"
#ifdef OPT_PYTHON
#include "Python.h"
#endif

#include "Print.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define NO_ANSI            ""

#define TAG_FATAL    "FATAL"
#define TAG_ERROR    "ERROR"
#define TAG_WARNING  "Warning"
#define TAG_DEBUG    "Debug"
#define TAG_INFO     "Info"

int lm::Print::_verbosityLevel =  6;

void lm::Print::printf(int verbosity, const char * fmt, ...)
{
    if (verbosity <= _verbosityLevel) {

        bool color = isatty(fileno(stdout));

        char buf[2048];
        const char *ansi0, *ansi1, *tag;

        time_t now;
        time(&now);
        struct tm nowParts;
        localtime_r(&now, &nowParts);

        va_list args;

        switch (verbosity) {
            case FATAL:
                ansi0 = ANSI_COLOR_RED;
                tag = TAG_FATAL;
                break;
            case ERROR:
                ansi0 = ANSI_COLOR_RED;
                tag = TAG_ERROR;
                break;
            case WARNING:
                ansi0 = ANSI_COLOR_YELLOW;
                tag = TAG_WARNING;
                break;
            case DEBUG:
            case VERBOSE_DEBUG:
                ansi0 = ANSI_COLOR_CYAN;
                tag = TAG_DEBUG;
                break;
            default:
                ansi0 = NO_ANSI;
                tag = TAG_INFO;
                break;
        }

        ansi1 = ANSI_COLOR_RESET;

        if (!color)
            ansi0 = ansi1 = NO_ANSI;


        va_start(args,fmt);
        vsnprintf(buf, 2048,fmt, args);
        va_end(args);

#ifdef OPT_PYTHON
        if (Py_IsInitialized()) {
            PyGILState_STATE gilstate = PyGILState_Ensure(); 
            PySys_WriteStdout("%04d-%02d-%02d %02d:%02d:%02d) %s%s%s: %s\n", 
                    nowParts.tm_year+1900, nowParts.tm_mon+1, nowParts.tm_mday, 
                    nowParts.tm_hour, nowParts.tm_min, nowParts.tm_sec,
                    ansi0, tag, ansi1, buf);
            PyObject *pyStdout = PySys_GetObject("stdout");
            if (pyStdout) {
                PyObject *result = PyObject_CallMethod(pyStdout,"flush",NULL);
                Py_XDECREF(result);
            }
            PyGILState_Release(gilstate); 
        } else  {
            std::printf("%04d-%02d-%02d %02d:%02d:%02d) %s%s%s: %s\n", 
                    nowParts.tm_year+1900, nowParts.tm_mon+1, nowParts.tm_mday, 
                    nowParts.tm_hour, nowParts.tm_min, nowParts.tm_sec,
                    ansi0, tag, ansi1, buf);
            std::fflush(stdout);
        }
#else
		std::printf("%04d-%02d-%02d %02d:%02d:%02d) %s%s%s: %s\n", 
				nowParts.tm_year+1900, nowParts.tm_mon+1, nowParts.tm_mday, 
				nowParts.tm_hour, nowParts.tm_min, nowParts.tm_sec,
				ansi0, tag, ansi1, buf);
                std::fflush(stdout);
#endif
    }
}

void lm::Print::printDateTimeString()
{
    time_t now;
    time(&now);
    struct tm nowParts;
    localtime_r(&now, &nowParts);
    // Form: "year-month-day hr:min:sec"
    ::printf("%04d-%02d-%02d %02d:%02d:%02d) ", nowParts.tm_year+1900, nowParts.tm_mon+1, nowParts.tm_mday, nowParts.tm_hour, nowParts.tm_min, nowParts.tm_sec);
}

std::string lm::Print::getDateTimeString()
{
    time_t now;
    time(&now);
    struct tm nowParts;
    localtime_r(&now, &nowParts);
    std::ostringstream s;
    // Form: "year-month-day hr:min:sec"
    s << nowParts.tm_year+1900 << "-" << nowParts.tm_mon+1 << "-" << nowParts.tm_mday << " " << nowParts.tm_hour << ":" << nowParts.tm_min << ":"  << nowParts.tm_sec;
    return s.str();
}


