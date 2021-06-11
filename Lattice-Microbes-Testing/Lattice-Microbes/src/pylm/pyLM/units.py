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

import lm

def _map(f,xs):
        return [f(x) for x in xs]

############################
# Length Wrapper Functions #
############################
def angstrom(*qty):
        """Returns a representation of a number in angstroms

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(angstrom, qty)
        else:
                return float(qty[0])*1e-10

def nm(*qty):
        """Returns a representation of a number in nanometers

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(nm, qty)
        else:
                return float(qty[0])*1e-9

def micron(*qty):
        """Returns a representation of a number in micrometers

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(micron, qty)
        else:
                return float(qty[0])*1e-6

def mm(*qty):
        """Returns a representation of a number in millimeters

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(mm, qty)
        else:
                return float(qty[0])*1e-3

def cm(*qty):
        """Returns a representation of a number in centimeters

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(cm, qty)
        else:
                return float(qty[0])*1e-2

##########################
# Time Wrapper Functions #
##########################
def ns(*qty):
        """Returns a representation of a number in nanoseconds

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(ns, qty)
        else:
                return float(qty[0])*1e-9

def microsecond(*qty):
        """Returns a representation of a number in microseconds

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(microsecond, qty)
        else:
                return float(qty[0])*1e-6

def ms(*qty):
        """Returns a representation of a number in milliseconds

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(ms, qty)
        else:
                return float(qty[0])*1e-3

def second(*qty):
        """Returns a representation of a number in seconds

        Seems silly, but for completeness and ability to annotate the unit in code.

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(second, qty)
        else:
                return float(qty[0])

def minute(*qty):
        """Returns a representation of a number in minutes

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(minute, qty)
        else:
                return float(qty[0])*60.0

def hr(*qty):
        """Returns a representation of a number in hours

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(hr, qty)
        else:
                return float(qty[0])*3600.0

def day(*qty):
        """Returns a representation of a number in days

        Args:
            qty:
                A list or singleton of a number
        """
        if(len(qty) > 1):
                return _map(day, qty)
        else:
                return float(qty[0])*3600.0*24.0


