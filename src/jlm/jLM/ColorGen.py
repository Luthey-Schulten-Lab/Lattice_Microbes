# 
# University of Illinois Open Source License
# Copyright 2016-2018 Luthey-Schulten Group,
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
# Author(s): Tyler M. Earnest
# 

"""Defines ColorSeq, the color scheme generator"""
import colorsys
import itertools
import fractions

# adapted from http://stackoverflow.com/a/13781114

def _coverRange():
    yield 0
    for n in itertools.count():
        k = fractions.Fraction(1, 2**n)
        i = k.denominator # [1,2,4,8,16,...]
        for j in range(1,i,2):
            yield fractions.Fraction(j,i)

def _hsvGen(h):
    for s in [fractions.Fraction(6,10)]: # optionally use range
        for v in [fractions.Fraction(8,10),fractions.Fraction(5,10)]: # could use range too
            yield (h, s, v) # use bias for v here if you use range


class ColorGen:
    """Distinct color generator"""
    def __init__(self,hue=0.5):
        """Distinct color generator

        Args:
           hue (float):  Initial hue
        """
        self._hue = hue
        self._genHSV = itertools.chain.from_iterable(map(_hsvGen,_coverRange()))
        self._fgHex = ["#000000"]
        self._bgHex = ["#ffffff"]
        self._fgFloat = [(0.0,0.0,0.0)]
        self._bgFloat = [(1.0,1.0,1.0)]

    def _get(self,idx):
        while len(self._bgHex) <= idx:
            h,s,v = map(float, next(self._genHSV))
            r,g,b = colorsys.hsv_to_rgb((h+self._hue)%1, s, v)

            self._bgFloat.append((r,g,b))
            self._bgHex.append("#{:02x}{:02x}{:02x}".format(*map(lambda x: int(255*x), (r,g,b))))
            val =  1 - ( 0.299 *r + 0.587 * g + 0.114 * b)

            if val > 0.44:
                self._fgHex.append("#ffffff")
                self._fgFloat.append((1.0,1.0,1.0))
            else:
                self._fgHex.append("#000000")
                self._fgFloat.append((0.0,0.0,0.0))


    def floatFg(self,idx):
        """Get an appropriate foreground color for a color index as a RGB triple

        Args:
            idx (int): 
                Color index

        Returns:
            (float,float,float):
                Rgb triple
        """
        self._get(idx)
        return self._fgFloat[idx]

    def strFg(self,idx):
        """Get an appropriate foreground color for a color index as a hex string

        Args:
            idx (int):
                Color index

        Returns:
            str: 
                Hex string
        """
        self._get(idx)
        return self._fgHex[idx]

    def float(self,idx):
        """Get a RGB triple by index

        Args:
            idx (int): 
                Color index

        Returns:
            (float,float,float):
                Rgb triple
        """
        self._get(idx)
        return self._bgFloat[idx]

    def str(self,idx):
        """Get a hex string by index

        Args:
            idx (int): 
                Color index

        Returns:
            str: 
                Hex string
        """
        self._get(idx)
        return self._bgHex[idx]

ColorSeq = ColorGen()
