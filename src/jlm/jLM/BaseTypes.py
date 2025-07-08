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

"""Simulation object base class and collections"""

import re
from . import Template

from abc import ABCMeta, abstractmethod

def _goodName(s):
    return (len(s) == len(s.encode()) and       # ensure no weird unicode stuff
            not any(x.isspace() for x in s) and # no whitespace
            ',' not in s)                       # comma used as a delimiter in HDF5

class Namespace:
    """Convienence class for TAB completing simulation objects"""
    def __init__(self, sos):
        """Convienence class for TAB completing simulation objects

        Args:
            sos (:py:class:`~jLM.BaseTypes.SimObjs`):
                Object collection
        """
        self._simobj = sos

    def __getattr__(self, name):
        obj = self._simobj._obj.get(name)

        if obj is None:
            raise AttributeError("Unknown object: {}".format(name))
        else:
            return obj

    def __dir__(self):
        return [o.name for o in self._simobj]

    def _repr_html_(self):
        return Template.j2render("objContainer.html",
                                 dict(title=self._simobj.__repr__(),
                                      annotationCol=any(o.annotation is not None for o in self._simobj),
                                      texStr=any(o._texstr is not None for o in self._simobj),
                                      objs=[dict(id=o.idx, html=o._repr_html_(), annotation=o.annotation, tex=o._texstr) 
                                               for o in self._simobj]))

class SimObjs:
    """Container to manage simulation data types"""
    def getAutoNamespace(self):
        """Get :py:class:`~jLM.BaseTypes.Namespace` object."""
        return Namespace(self)

    def matchRegex(self, regex):
        """Generator returning matches to regular expression

        Args:
            regex (str): 
                Regular expression matching the name attribute of 
                    contained objects

        Returns:
           generator:
                Generator of :py:class:`~jLM.Types.SimObj` 
        """
        for v in sorted(self._obj.values(), key=lambda x: -x.idx):
            if re.match(regex, v.name):
                yield v

    def __getitem__(self, what):
        """Get object by name or idx integer"""
        try:
            return self._obj[what]
        except KeyError:
            return self._id2obj[what]

    def __len__(self):
        return len(self._obj)

    def __iter__(self):
        return (self[i+self._idbase] for i in range(len(self._obj)))

    def __init__(self, sim, cls, idbase=0):
        self._cls = cls
        self._sim = sim
        self._idbase = idbase
        self._obj = dict()
        self._id2name = dict()
        self._id2obj = dict()
        self._frozen = False

    def freeze(self):
        """Do not create new objects from failed lookups"""
        self._frozen = True

    def unfreeze(self):
        """Create new objects from failed lookups"""
        self._frozen = False

    def is_defined(self, *args, **kwargs):
        """Check if `get` would create a new object

        Args:
            *args:
                Arguments to constructor

        Keyword Args:
            **kwargs:
                Keyword arguments to constructor

        Returns:
            bool: True if would create new object.
        """
        name = self._cls._unique_id(*args, **kwargs)
        return name in self._obj

    def get(self, *args, **kwargs):
        """Return a simulation object, create it if it does not exist

        Args:
            *args:
                Arguments to constructor

        Keyword Args:
            **kwargs:
                Keyword arguments to constructor

        Returns:
            :py:class:`~jLM.BaseTypes.SimObj`: 
                Object
        """
        name = self._cls._unique_id(*args, **kwargs)

        try:
            return self._obj[name]
        except KeyError:
            if not self._frozen:
                idx = len(self._obj)+self._idbase
                obj = self._cls(self._sim, idx, *args, **kwargs)
                self._obj[name] = obj
                self._id2name[obj.idx] = obj.name
                self._id2obj[obj.idx] = obj
                return obj
            else:
                raise



class SimObj(metaclass=ABCMeta):
    """Base class to simulation objects"""
    @classmethod
    @abstractmethod
    def _unique_id(cls, *args, **kwargs):
        pass

    def _TeXMath(self):
        if hasattr(self, '_texstr') and self._texstr is not None:
            return self._texstr 
        else:
            return self._TeX()

    def __init__(self, sim, idx, *args, **kwargs):
        if 'annotation' in kwargs:
            self.annotation = kwargs['annotation']
            del kwargs['annotation']
        else:
            self.annotation = None

        name = self._unique_id(*args, **kwargs)
        if 'texRepr' in kwargs:
            self._texstr = kwargs['texRepr']
            del kwargs['texRepr']
        else:
            self._texstr = None
        if not _goodName(name):
            raise ValueError("Invalid name: "+name)
        self.name = name
        self._sim = sim
        self.idx = idx
        self._setup(*args, **kwargs)
        self._staticAttrs = set(dir(self))

    def _dynamicAttrs(self):
        return [ x for x in set(dir(self)) - self._staticAttrs if x[0] != "_" ]

    def _setup(self, *args, **kwargs):
        pass

    def __lt__(lhs,rhs):
        return lhs.idx < rhs.idx

    def _repr_html_(self):
        try:
            return self._html()
        except AttributeError:
            return "$$"+self._TeXMath()+"$$"

