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

"""Auxiliary helper functions for presentation"""
import numpy as np
import re

def _htmlExponent(x):
    return "<sup>" + str(int(x)) + "</sup>"

def texUnits(s):
    """Reformat ascii number as TeX

    Args:
        s (str):
            Formatted number. See :py:func:`~jLM.DisplayUtils.numToStr`

    Returns:
        str:
            TeX math mode representation.
    """
    s = s.replace("*",r"\cdot ")
    split = re.split("\^\{([-0-9]+)\}",s)
    split[0::2] = [r"\mathrm{"+u+"}" if len(u) > 0 else '' for u in split[0::2]]
    split[1::2] = ["^{"+e+"}" if len(e) > 0 else '' for e in split[1::2]]
    return ''.join(split)

def unicodeUnits(s):
    """Reformat ascii number as unicode

    Args:
        s (str):
            Formatted number. See :py:func:`~jLM.DisplayUtils.numToStr`

    Returns:
        str:
            Unicode formatted number.
    """
    s = s.replace("*","⋅")
    split = re.split("\^\{([-0-9]+)\}",s)
    split[1::2] = [_htmlExponent(e) for e in split[1::2]]
    return ''.join(split)

def numToStr(val, unit=None, scale=1, n=3, smallest_exp=0):
    """Pretty print a number

    Units can be specified in the format, e.g. "m^{2}*s^{-1}".
    In this format, it can be converted between TeX and unicode 
    representations.

    Args:
        val (float):
            Value
        unit (str):
            Unit
        scale (float):
            Factor to scale value
        n (int):
            Precision
        smallest_exp (int):
            Smallest absolute value of exponent to write in scientific notation

    Returns:
        str: 
            Formatted number
    """
    if float("-inf") < val < float("inf"):
        val *= scale

        if abs(val) == 0:
            if n > 0:
                s = "0." + "0"*n
            else:
                s = "0"
        else:
            sgn = val/abs(val)
            x = val*sgn

            exponent = round(np.log10(x))
            mantessa = round(x*10**-exponent,n)

            if mantessa < 1:
              exponent -= 1
              mantessa = sgn*round(x*10**-exponent,n)

            if abs(exponent) < abs(smallest_exp):
              s = "{:.{n}f}".format(val, n=n)
            else:
              s = "{v:.{n:}f} × 10{exp}".format(v=mantessa,n=n, exp=_htmlExponent(exponent))

        if unit is not None:
            s = s + ' ' + unicodeUnits(unit)
    else:
        s = str(val)

    return s

def toHex(c):
    """RGB triple to hex string

    Args.
        c (float,float,float): 
            RGB tuple of values 0 <= c <= 1

    Returns:
        str:
            HTML RGB hex code.
    """
    return "#{:02x}{:02x}{:02x}".format(*map(lambda x: int(255*x), c))



def colorWheel(h):
    """RGB triple at a particular hue

    Args:
        h (float):
            Hue in range (0,1)
    
    Returns:
        (int,int,int):
            RGB triple
    """
    s=0.9
    v=0.5
    # s=0.5
    # v=0.94
    h_i = int(h*6)
    f = h*6 - h_i
    p = v * (1 - s)
    q = v * (1 - f*s)
    t = v * (1 - (1 - f) * s)
    if h_i==0:
        c =(v, t, p)
    elif h_i==1:
        c = (q, v, p)
    elif h_i==2:
        c = (p, v, t)
    elif h_i==3:
        c = (p, q, v)
    elif h_i==4:
        c = (t, p, v)
    else:
        c = (v, p, q)
    return c

def texfmt(xx,n=3,tex_escape=False,smallest_exp=3):
  """
  Format a float in LaTeXified scientific notation.

  Args:
      xx (float):
          Number

  Keyword Args:
      n (int):
          Precision
      tex_escape (bool):  
          If true, wrap result in inline math.
      smallest_exp (int): 
          Smallest absolute value of exponent to print in 
          scientific notation.

  Returns:
      str: 
          Formatted number
  """
  if abs(xx) == 0:
    s = "0.0"
  else:
    sgn = xx/abs(xx)
    x = xx*sgn

    exponent = round(np.log10(x))
    mantessa = round(x*10**-exponent,n)
    if mantessa < 1:
      exponent -= 1
      mantessa = sgn*round(x*10**-exponent,n)

    if abs(exponent) < abs(smallest_exp):
      s = "{:.{n}f}".format(xx, n=n)
    else:
      s = str(mantessa) + r'\times 10^{' + str(int(exponent)) + '}' 


  if tex_escape:
    return '$'+s+'$'
  else:
    return s


def sanitizeTex(st):
    """ Replace control characters with escaped versions

    Args:
        st (str): 
            Input string

    Returns:
        str: 
            Escaped version of `st`
    """
    s = st.replace("\\", r'\textbackslash{}')
    s = s.replace('&', r'\&')
    s = s.replace('%', r'\%')
    s = s.replace('$', r'\$')
    s = s.replace('_', r'\_')
    s = s.replace('{', r'\{')
    s = s.replace('}', r'\}')
    s = s.replace('~', r'\textasciitilde{}')
    s = s.replace('^', r'\textasciicircum{}')
    return s

