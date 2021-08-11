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

"""Templating engine tools"""

import os.path
try:
    import jinja2
    _j2 = jinja2.Environment(trim_blocks=True, 
                             lstrip_blocks=True,
                             loader=jinja2.FileSystemLoader(os.path.join(os.path.dirname(__file__),"templates")))
except ImportError:
    _j2 = None

def j2render(fname, ctx):
    """Renders a template given a context dictionary

    Args:
        fname (str):
            template file in the "templates" directory
        ctx (dict):
            context dictionary

    Returns:
        str: Rendered template
    """
    if _j2:
        return _j2.get_template(fname).render(ctx)
    else:
        raise ImportError("Deferred jinja2 import exception")

def displayj2html(fname, ctx):
    """Renders an HTML template and displays in Jupyter notebook

    Args:
        fname (str):
            template file in the "templates" directory
        ctx (dict):
            context dictionary
    """
    import IPython.display as ipd
    if _j2:
        return ipd.HTML(j2render(fname, ctx))
    else:
        raise ImportError("Deferred jinja2 import exception")
