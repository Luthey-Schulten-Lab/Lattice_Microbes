# 
# University of Illinois Open Source License
# Copyright 2008-2018 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Michael J. Hallock and Joseph R. Peterson
# 
#

import logging

# Create the default LM Logger
LMLogger = logging.getLogger('LMLogger')

# Gracefully catch logging.NullHandler() not existing in python 2.6
try:
	nullHandlerLM = logging.NullHandler()
except AttributeError:
	class NullHandler(logging.Handler):
		def emit(self, record):
			pass
	nullHandlerLM = NullHandler()

LMLogger.addHandler(nullHandlerLM) # Set the library to initially print to nothing

# Set up the formatter
LMformatter = logging.Formatter('%(asctime)s: %(levelname)s: %(message)s')
nullHandlerLM.setFormatter(LMformatter)



def setLMLoggerLevel(level):
	"""Set the level of the logger for the application

    Args:
        level:
            The level the logger should report (e.g. logger.WARNING, logger.INFO, etc.)
    """
	LMLogger.setLevel(level)

def setLMLogFile(filename, level=logging.DEBUG):
	"""Set up file handler to print log to file

    Args:
        filename:
            The name of the file to log information
        level:
            The level of information to log
    """
	fileH = logging.FileHandler(filename, mode='w')
	fileH.setLevel(level)
	fileH.setFormatter(LMformatter)
	LMLogger.removeHandler(nullHandlerLM)
	LMLogger.addHandler(fileH)

def setLMLogConsole(level=logging.DEBUG):
	"""Set the logger to write to the console as the code is working

    Args:
        level:
            The level of information to log
    """
	consoleH = logging.StreamHandler() # Defaults to sys.stderr
	consoleH.setLevel(level)
	consoleH.setFormatter(LMformatter)
	LMLogger.removeHandler(nullHandlerLM)
	LMLogger.addHandler(consoleH)

