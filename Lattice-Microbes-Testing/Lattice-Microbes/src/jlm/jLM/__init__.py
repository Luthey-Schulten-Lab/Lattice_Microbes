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

"""
jLM top level.

.. py:data:: buildData

   Dictionary containing  data about build.

   buildTime
      Time/date package was built

   name
      Name of package. Should be ``jLM``

   version
      Version number of ``jLM``

   gitHash
      Git hash of the last commit

   user
      Username that built the package

   host
      Hostname of machine where the package was built
"""

import json, os

try:
    buildData = json.load(open(os.path.join(os.path.dirname(__file__),"build.json")))
except FileNotFoundError:
    import time, getpass, socket
    buildData = dict(buildTime=time.strftime("%Y-%m-%d %H:%M:%S"),
                     name='jLM',
                     version="DEVEL",
                     gitHash="<<<Development run>>>",
                     user=getpass.getuser(),
                     host=socket.gethostname())

try:
    import IPython.core.getipython as gipy
    import IPython.display as ipd
    from . import Template
    if gipy.get_ipython() is not None:
        ipd.display(Template.displayj2html("jupyterInit.html", {}))
except ImportError:
    pass

def printCompiledOptions():
    import lm
    print(lm.BUILD_CONFIGURATION)
