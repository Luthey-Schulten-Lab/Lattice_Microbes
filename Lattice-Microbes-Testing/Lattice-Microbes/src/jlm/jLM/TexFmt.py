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

"""Dump simulation parameters as LaTeX"""
import sys
import numpy as np
import datetime

from . import DisplayUtils as Dsp
from . import Template as Tmp

def _texHeader(net, f):
    print("% {0.name}: ({0.filename}), {1:%Y}-{1:%m}-{1:%d} {1:%H}:{1:%M}:{1:%S}".format(net, datetime.datetime.now()), file=f)
    
def regionTable(net, regions=None, file=None):
    """Write a TeX formated table of regions

    Args:
        net (:py:class:`~jLM.RDME.SpatialModel`):
            Model
        regions ([:py:class:`~jLM.Types.Region`]):
            Regions to report
        file (file-like):
            File object of output
    """
    file = file or sys.stdout
    if regions is None:
        regions = list(net.regionList) 

    ctx = dict(regions=[])
    cs = net.particleStatistics()
    annotate = False

    for r in regions:
        annotate |= r.annotation is not None
        ctx['regions'].append(
                dict( label=r._TeX(),
                      counts=cs["countByRegion"][r],
                      conc=Dsp.texfmt(cs["concByRegion"][r]*1e6),
                      volume=Dsp.texfmt(cs["regionVol"][r]*1e15),
                      siteCount=cs["regionCounts"][r],
                      annotation=r.annotation,
                      occ=Dsp.texfmt(cs["countByRegion"][r]/cs["regionCounts"][r]/net.pps)))
    ctx['annotationCol'] = annotate

    _texHeader(net,file)
    file.write(Tmp.j2render("regionTable.tex", ctx))


def parameterTable(net, params=None, file=None):
    """Write a TeX formated table of parameters

    Args:
        net (:py:class:`~jLM.RDME.SpatialModel`):
            Model
        parameters ([:py:class:`~jLM.Types.RateConst`]):
            Parameters to report
        file (file-like):
            File object of output
    """
    file = file or sys.stdout
    if params is None:
        params = list(net.rxnRateList) + list(net.diffRateList)
    ctx = dict(count=len(params), params=[])
    counts = {k:0 for k in params}

    ds = {d.idx: d for d in params if d in net.diffRateList}

    for rx in net.reactionList:
        try:
            counts[rx.rate] += 1
        except KeyError:
            pass

    for _, _, _, r in net._transitionRates:
        try:
            counts[ds[r]] += 1
        except KeyError:
            pass

    annotate = False
    for p in params:
        annotate |= p.annotation is not None
        ctx['params'].append(dict(id=p.idx, 
                                 symbol=p._TeXMath(), 
                                 annotation=p.annotation,
                                 value=Dsp.texfmt(p.value),
                                 unit=Dsp.texUnits(p._unit),
                                 count=counts[p])) 

    ctx['annotationCol'] = annotate

    _texHeader(net,file)
    file.write(Tmp.j2render("parameterTable.tex", ctx))


def reactionTable(net, rxns=None, file=None):
    """Write a TeX formated table of reactions

    Args:
        net (:py:class:`~jLM.RDME.SpatialModel`):
            Model
        rxns ([:py:class:`~jLM.Types.Reaction`]):
            Reactions to report
        file (file-like):
            File object of output
    """
    file = file or sys.stdout
    if rxns is None:
        rxns = list(net.reactionList)

    ctx = dict(rxns=[])
    annotate = False
    for r in rxns:
        annotate |= r.annotation is not None or r.rate.annotation is not None
        ctx['rxns'].append(dict(annotation=r.annotation if r.annotation is not None else r.rate.annotation,
                                idx=r.idx,
                                rxn= r._TeXMath(),
                                reg=r', '.join(net.regionList[rgidx]._TeX() for rxidx, rgidx in net._reactionLocations if rxidx==r.idx),
                                rate=Dsp.texfmt(r.rate.value),
                                rateUnit=Dsp.texUnits(r.rate._unit)))

    ctx['annotationCol'] = annotate

    _texHeader(net,file)
    file.write(Tmp.j2render("reactionTable.tex", ctx))


def modelTables(net, file=None):
    """Write a TeX formated description of the model

    Args:
        net (:py:class:`~jLM.RDME.SpatialModel`):
            Model
        file (file-like):
            File object of output
    """
    file = file or sys.stdout
    ctx = dict(name                 = net.name,
               filename             = net.filename,
               nsps                 = len(net.speciesList),
               nrxns                = len(net.reactionList),
               nsts                 = len(net.regionList),
               nrconst              = len(net.rxnRateList),
               ndconst              = len(net.diffRateList),
               siteVol              = Dsp.texfmt(1e15*net.siteV),
               volumeUnit           = r"\mathrm{fl}",
               latticeSpacing       = Dsp.texfmt(1e9*net.latticeSpacing),
               lengthUnit           = r"\mathrm{nm}",
               nPlacedParticles     = len(net._particlePlacement),
               pps                  = net.pps,
               bytesPerParticle     = net.bytesPerParticle,
               regions              = [],
               timeStepUnit         = r'\upmu\mathrm{s}',
               timeUnit             = r'\mathrm{s}',
               timeStep             = Dsp.texfmt(1e6*net.timestep),
               simTime              = Dsp.texfmt(net.simulationTime),
               hookInterval         = Dsp.texfmt(net.hookInterval),
               latticeWriteInterval = Dsp.texfmt(net.latticeWriteInterval),
               speciesWriteInterval = Dsp.texfmt(net.speciesWriteInterval),
               concUnit             = r'\upmu\mathrm{M}')

    try:
        ctx['nReplicates']= len(net.h5['Simulations'])
        ctx['replicates'] = []

        for label in net.h5['Simulations']:
            tsp = net.h5['Simulations'][label]['SpeciesCountTimes'][-1]
            tlt = net.h5['Simulations'][label]['LatticeTimes'][-1]
            nsp = net.h5['Simulations'][label]['SpeciesCountTimes'].shape[0]
            nlt = net.h5['Simulations'][label]['LatticeTimes'].shape[0]
            ct0 = np.sum(net.h5['Simulations'][label]['SpeciesCounts'][0,:])
            ct1 = np.sum(net.h5['Simulations'][label]['SpeciesCounts'][-1,:])

            ctx['replicates'].append(dict(label=label, time=Dsp.texfmt(max(tsp,tlt)),
                                          latticeEvals=nlt, speciesEvals=nsp, startCt=ct0, endCt=ct1))
    except AttributeError:
        pass


    ctx['lz'], ctx['ly'], ctx['lx'] = net.shape

    _texHeader(net,file)
    file.write(Tmp.j2render("siminfo.tex", ctx))
