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

"""SpatialModel mixins for rich display in Jupyter"""
import itertools, io, base64, random, string

import numpy as np
import warnings
from . import Template as Tmp
from . import Lattice
from . import DisplayUtils as Dsp

def _maybeJupyter(fn):
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except ImportError:
            warnings.warn("Call to {} failed due to missing module dependencies".format(fn.__name__), stacklevel=2)
            return None
    return wrapped

def _num(x, scale=1, unit=None):
    return  Dsp.numToStr(x, scale=scale, unit=unit, n=3, smallest_exp=3)

def _numFmt(fmt, *nums):
    ss = []
    for v in nums:
        try:
            x, scale, unit, n = v
        except (ValueError, TypeError):
            x, scale, unit, n = v, 1.0, None, 3
        ss.append(Dsp.numToStr(x, scale=scale, unit=unit, n=n, smallest_exp=3))
    return fmt.format(*ss)


def _maxEntropySlice(siteLattice, axis):
    maxEntropy, planeIndex = 0,0

    for i in range(siteLattice.shape[axis]):
        s = [slice(None), slice(None), slice(None)]
        s[axis] = i
        im = siteLattice[tuple(s)]
        nbins = siteLattice.max()+1
        e = sum(-x*np.log(x) for x in np.histogram(im,bins=nbins)[0]/im.size if x>0)

        if e > maxEntropy:
            maxEntropy, planeIndex = e, i

    return planeIndex

def _siteTypeImg(siteLattice, siteColors, axis0, axis1, axis2,  planeIndex):
    latticeShape = siteLattice.shape

    stRGB = np.zeros((latticeShape[axis0],latticeShape[axis1],3))
    if planeIndex is None:
        planeIndex = latticeShape[axis2]//2

    def rgb8(i):
        return np.array([int(255*c) for c in siteColors[i]], dtype=np.uint8)

    stImg = np.zeros((latticeShape[axis0],latticeShape[axis1]),dtype=np.uint8)

    for i,j in itertools.product(range(latticeShape[axis0]),range(latticeShape[axis1])):
        idx=[0,0,0]
        idx[axis0] = latticeShape[axis0]-i - 1
        idx[axis1] = j
        idx[axis2] = planeIndex
        stImg[i,j] = siteLattice[idx[0],idx[1],idx[2]]

    stRGB = np.zeros((stImg.shape[0],stImg.shape[1],3),dtype='uint8')
    for i,j in itertools.product(*map(range,stImg.shape)):
        stRGB[i,j,:] = rgb8(stImg[i,j])

    imgFile = io.BytesIO()
    import PIL.Image
    img = PIL.Image.fromarray(stRGB)
    img.save(imgFile,'PNG')
    return base64.b64encode(imgFile.getvalue()).decode()

class FileJupyterMixin:
    def _speciesJ2context(self, sps, cs): 
        ctx = super()._speciesJ2context(sps,cs)

        ns = self.speciesCounts[sps]
        ts = self.speciesCountTimes
        minIdx = np.argmin(ns)
        maxIdx = np.argmax(ns)
        trajData = dict(replicate=self.replicate,
                        duration=ts[-1],
                        startCt=ns[0],
                        endCt=ns[-1],
                        meanCt=_numFmt("{} ± {}", np.mean(ns), np.std(ns)),
                        minCt=_numFmt("{} @ t={}", (ns[minIdx], 1, None, 0), (ts[minIdx], 1.0, 's', 3)),
                        maxCt=_numFmt("{} @ t={}", (ns[maxIdx], 1, None, 0), (ts[maxIdx], 1.0, 's', 3)))

        ts = self.speciesCountTimes
        if ns.shape[0] > 256:
            stride = int(ns.shape[0]/256)
            ns = ns[::stride]
            ts = ts[::stride]

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6,4))
        ax.plot(ts,ns)

        if ts[0] == ts[-1]:
            xlim=(0,1)
        else:
            xlim=(ts[0], ts[-1])

        ns0, ns1 = int(ns.min()/1.1), int(ns.max()*1.1)
        ylim=(ns0-1,ns1+1)

        ax.set(xlim=xlim,
               ylim=ylim,
               xlabel=r"$t/\mathrm{s}$", 
               ylabel=r"$n_{\mathrm{"+sps.name+"}}$")
        fig.tight_layout()
        imgFile = io.BytesIO()
        fig.savefig(imgFile, format="SVG", transparent=True)
        plt.close(fig)
        lines = imgFile.getvalue().decode().splitlines()

        while '<svg' not in lines[0]:
            lines.pop(0)

        trajData['svg'] = '\n'.join(lines)

        ctx['trajData'] = trajData

        return ctx

    def _modelJ2context(self): 
        ctx = super()._modelJ2context()

        ctx['nReplicates']= len(self.h5['Simulations'])
        ctx['replicates'] = []

        for label in self.h5['Simulations']:
            tsp = self.h5['Simulations'][label]['SpeciesCountTimes'][-1]
            tlt = self.h5['Simulations'][label]['LatticeTimes'][-1]
            nsp = self.h5['Simulations'][label]['SpeciesCountTimes'].shape[0]
            nlt = self.h5['Simulations'][label]['LatticeTimes'].shape[0]
            ct0 = np.sum(self.h5['Simulations'][label]['SpeciesCounts'][0,:])
            ct1 = np.sum(self.h5['Simulations'][label]['SpeciesCounts'][-1,:])

            ctx['replicates'].append(dict(label=label, time=_num(max(tsp,tlt),unit='s'), 
                                          latticeEvals=nlt, speciesEvals=nsp, startCt=ct0, endCt=ct1))
        
        return ctx




class JupyterDisplayMixin:
    def _repr_html_(self):
        return Tmp.j2render("siminfo.html", self._modelJ2context())

    # Jinja2 context builders. If overidden, call super()._...J2context()

    def _speciesJ2context(self, sps, cs):
        rxns = self._reactionJ2context([r for r in self.reactionList 
                                          if sps in r.products or sps in r.reactants])

        transitionTable = np.zeros((len(self.regionList),len(self.regionList), len(self.speciesList)+1), dtype=np.int)
        transitionTable[...] = -1

        for sidx, r0idx, r1idx, rt in self._transitionRates:
            if sidx is None:
                sidx = slice(None)
            if r0idx is None:
                r0idx = slice(None)
            if r1idx is None:
                r1idx = slice(None)

            transitionTable[r0idx, r1idx, sidx] = rt

        txb = transitionTable[:,:,sps.idx]

        def style(x,y=None):
            d = dict(v=x)
            if y is not None:
                d['style'] = y
            return d

        def getDr(r0,r1):
            didx = txb[r0.idx,r1.idx]
            try:
                v = self.diffRateList[didx].value
                s = _num(v, 1e12)
            except KeyError:
                didx = -1

            if didx < 0:
                return style("undef", "undef")
            elif v < 1e-20:
                return style(s, "zero")
            elif v >= self.maxDiffusionRate():
                return style(s, "inf")
            else:
                return style(s)

        dfTable = ( [ [style('<span class="jLMnum">× µm<sup>2</sup>⋅s<sup>-1</sup></span>')] + [style(r._html()) for r in self.regionList] ] +
                    [ [style(r0._html())] + [ getDr(r0, r1) for r1 in self.regionList] for r0 in self.regionList] )

        placementList = [(x,y,z) for x,y,z,s in self._particlePlacement if s==sps.idx]

        attrs = dict()
        for attr in sps._dynamicAttrs():
            v = getattr(sps, attr)
            if isinstance(v, str):
                attrs[attr] = v
            elif np.isscalar(v):
                attrs[attr] = _num(v)
            else:
                attrs[attr] = repr(v.__class__)
                    
        ctx = dict(label=sps.name,
                   totalCount=cs['countBySpecies'][sps],
                   otherAttrs=attrs,
                   idx=sps.idx,
                   name=sps.name,
                   tex=sps._TeXMath(),
                   annotation=sps.annotation,
                   totalReactions=len(rxns['rxns']),
                   distribution=[ dict(region=st._html(), 
                                       count=cs['countBySpeciesRegion'][sps][st], 
                                       conc=_num(cs['concBySpeciesRegion'][sps][st], 1e6, "µM"))
                                  for st in self.regionList  if cs['countBySpeciesRegion'][sps][st] > 0],
                   rxns=rxns['rxns'],
                   placement=[],
                   nPlaced=len(placementList),
                   diffTbl=dfTable)

        def lookupRegion(v):
            return self.regionList[self.siteLattice[v[2],v[1],v[0]]]._html() 

        if 0 < len(placementList) <= 200:
            placement= {v:dict(count=0,x=v[0],y=v[1],z=v[2],region=lookupRegion(v)) for v in placementList}
            for v in placementList:
                placement[v]['count'] += 1

            ctx['placement'] = list(placement.values())
        else:
            pd = dict()
            for v in placementList:
                try:
                    pd[lookupRegion(v)] += 1
                except KeyError:
                    pd[lookupRegion(v)] = 1

            ctx['placement'] = [ dict(count=v, region=k) for k,v in pd.items() ]


        return ctx

    def _reactionJ2context(self, rxns):
        ctx = dict(rxns=[])
        annotate = False
        for r in rxns:
            annotate |= r.annotation is not None or r.rate.annotation is not None
            ctx['rxns'].append(dict(
                 idx=r.idx,
                 annotation=r.annotation or r.rate.annotation,
                 rxn="$" + r._TeXMath() + "$",
                 reg=r', '.join(self.regionList[rgidx]._html() for rxidx, rgidx in self._reactionLocations if rxidx==r.idx),
                 rate=_num(r.rate.value,  unit=Dsp.unicodeUnits(r.rate._unit))))
        ctx['annotationCol'] = annotate
        return ctx

    def _regionJ2context(self, plane, planeIndex, maxWidth=600, maxHeight=600):
        axis0 = 'xyz'.index(plane[0])
        axis1 = 'xyz'.index(plane[1])
        axis2 = (set((0,1,2)) - set((axis0,axis1))).pop()

        latticeShape = self.siteLattice.shape

        if planeIndex is None:
            planeIndex = latticeShape[axis2]//2

        xscl = maxWidth/latticeShape[axis0]
        yscl = maxHeight/latticeShape[axis1]

        def proc_scl(s):
            if s >= 1:
                return yscl
            else:
                inv = 1/s
                inv2 = 2**int(np.ceil(np.log2(inv)))
                return 1/inv2

        scl = proc_scl(min(xscl, yscl))

        ctx = dict(w=int(scl*latticeShape[axis0]),
                   h=int(scl*latticeShape[axis1]),
                   plane=plane,
                   reg=[r._html() for r in self.regionList],
                   b64=_siteTypeImg(self.siteLattice, 
                                   [x._floatColor() for x in self.regionList], 
                                   axis0, axis1, axis2, 
                                   planeIndex),
                   dir="xyz"[axis2],
                   planeIndex=planeIndex, 
                   renderId=_renderId())

        return ctx

    def _modelJ2context(self):
        ctx = dict()

        cs = self.particleStatistics()

        dims = np.array(self.shape)
        axis0 = np.argmax(dims)
        dims[axis0] = -1
        axis1 = np.argmax(dims)
        dims = np.array(self.shape)
        if dims[axis0] == dims[axis1]:
            axis0,axis1 = min(axis0,axis1), max(axis0,axis1)
        else:
            axis0 = axis0 if dims[axis0]>dims[axis1] else axis1
            axis1 = axis1 if dims[axis0]>dims[axis1] else axis0
        axis2 = (set((0,1,2)) - set((axis0,axis1))).pop()

        # find most interesting slice from image entropy

        planeIndex = _maxEntropySlice(self.siteLattice, axis2)

        ctx = self._regionJ2context("xyz"[axis0]+"xyz"[axis1], planeIndex, maxWidth=700, maxHeight=500)
        ctx['name'] = self.name
        ctx['filename'] = self.filename
        ctx['nsps'] = len(self.speciesList)
        ctx['nrxns'] = len(self.reactionList)
        ctx['nsts'] = len(self.regionList)
        ctx['nrconst'] = len(self.rxnRateList)
        ctx['ndconst'] = len(self.diffRateList)
        ctx['dims'] = dict(x=dims[2],y=dims[1],z=dims[0])
        ctx['dims'] = _numFmt("{} × {} × {}", (dims[2], 1, None, 0), (dims[1], 1, None, 0), (dims[0], 1, None, 0))
        ctx['siteVol'] = _num(self.siteV, 1e15, "fl")
        ctx['latticeSpacing'] = _num(self.latticeSpacing, 1e9, "nm")
        ctx['nPlacedParticles'] = len(self._particlePlacement)
        ctx['pps'] = self.pps
        ctx['bytesPerParticle'] = self.bytesPerParticle
        ctx['regions'] = []

        def maybeAttr(attr, cattr, scl, unit):
            if hasattr(self, attr) and getattr(self,attr):
                ctx[cattr] = _num(getattr(self, attr), scl, unit)
            else:
                ctx[cattr] = "undefined"

        maybeAttr("timestep", "timeStep", 1e6, "µs")
        maybeAttr("simulationTime", "simTime", 1e0, "s")
        maybeAttr("latticeWriteInterval", "latticeWriteInterval", 1e0, "s")
        maybeAttr("speciesWriteInterval", "speciesWriteInterval", 1e0, "s")
        maybeAttr("hookInterval", "hookInterval", 1e0, "s")

        for r in self.regionList:
            c = dict()
            c['label'] = r._html()
            c['counts'] = cs["countByRegion"][r]
            c['conc'] = _num(cs["concByRegion"][r], 1e6, r"µM")
            c['volume'] = _num(cs["regionVol"][r], 1e15, r"fl")
            c['siteCount'] = cs["regionCounts"][r]
            occScl = 100.0/(max(1.0, cs["regionCounts"][r])*self.pps)
            c['occ'] = r"${:.2f}\,\%$".format(occScl*cs["countByRegion"][r])
            ctx['regions'].append(c)

        return ctx

    def _x3dJ2context(self, filterFunctions=None):
        siteLattice = self.siteLattice
        if filterFunctions is None:
            filterFunctions = dict()

        latticeDims = siteLattice.shape
        scl = 7/max(latticeDims)
        centroid = 0.5*np.array(latticeDims)

        nx,ny,nz = latticeDims

        jc = dict(centroid=centroid, scl=scl, renderId=_renderId(), species=[], sites=[], modelName=self.name)

        for i,reg in enumerate(self.regionList):
            if reg in filterFunctions:
                siteMatch = (siteLattice == reg.idx)
                posMatch = filterFunctions[reg](*np.mgrid[0:nx, 0:ny, 0:nz])
                binaryLattice = np.array(siteMatch&posMatch, dtype=np.uint8)
            else:
                binaryLattice = np.array(siteLattice == reg.idx, dtype=np.uint8)
            if binaryLattice.any():
                verts, faces = Lattice.greedyMesh(binaryLattice)
                r,g,b = reg._floatColor()
                jc['sites'].append(dict(name=reg.name, 
                                        idx=i,
                                        hexColor=reg._hexColor(), 
                                        label=reg._html(),
                                        checked="" if reg.idx==0 else " checked", 
                                        faces=' '.join(str(x) for face in faces for x in face),
                                        verts=' '.join(str(x) for vert in verts for x in vert),
                                        r=r, g=g, b=b,
                                        choice="-1" if reg.idx == 0 else "0"))
                   
        particleTypes = sorted(set(s for _,_,_,s in self._particlePlacement))
        ncolors = len(particleTypes)
        particleTypes = [(s, i+len(jc['sites']), Dsp.colorWheel(i/ncolors)) for i,s in enumerate(particleTypes)]

        for s, objIdx, (r,g,b) in particleTypes:
            jc['species'].append(dict(hexColor=Dsp.toHex((r,g,b)), 
                                      r=r, g=g, b=b,
                                      label=self.speciesList[s]._html(), 
                                      idx=objIdx, 
                                      radius=0.5, 
                                      locs = [dict(x=x+0.5,y=y+0.5,z=z+0.5) 
                                               for (z,y,x,_) in 
                                                  filter(lambda x:x[3]==s, self._particlePlacement)]))

        jc['downloadX3D'] = False
        return jc

    #### Model builder context manager

    def construct(sim):
        """Track newly created model objects and display in Notebook.

        Context manager which tracks new species, reactions, etc., and displays
        a HTML summary when used in Jupyter"""
        class InteractiveConstruct:
            def _getKeys(self):
                return (set(x.name for x in sim.speciesList),
                        set(x.name for x in sim.regionList),
                        set(x.name for x in sim.rxnRateList),
                        set(x.name for x in sim.reactionList),
                        set(x.name for x in sim.diffRateList))
            def __init__(self):
                self._keys0 = self._getKeys()
            def __enter__(self):
                return self
            def __exit__(self, type, value, traceback):
                if type is None:
                    try:
                        import IPython.display as ipd
                    except ImportError:
                        warnings.warn("Construct summary not displayed", stacklevel=2)
                        return None

                    sps0, reg0, k0, rxn0, d0 = self._keys0
                    sps1, reg1, k1, rxn1, d1 = self._getKeys()
                    newSps = sorted((sim.speciesList[x] for x in sps1 - sps0), key=lambda x: x.idx)
                    newRxns = sorted((sim.reactionList[x] for x in rxn1 - rxn0), key=lambda x: x.idx)
                    newRegs = sorted((sim.regionList[x] for x in reg1 - reg0), key=lambda x: x.idx)
                    newKs = sorted((sim.rxnRateList[x] for x in k1 - k0), key=lambda x: x.idx)
                    newDs = sorted((sim.diffRateList[x] for x in d1 - d0), key=lambda x: x.idx)

                    spsCtx = dict(annotationCol=any(o.annotation is not None for o in newSps),
                                  texCol=any(o._texstr is not None for o in newSps),
                                  objs=[dict(id=o.idx, html=o._repr_html_(), annotation=o.annotation, tex=o._texstr)
                                           for o in newSps])

                    rxnsCtx = dict(annotationCol= (any(o.annotation is not None for o in newRxns)
                                                    or any(o.rate.annotation is not None for o in newRxns)),
                                   objs=[dict(id=o.idx, 
                                               html=o._repr_html_(), 
                                               annotation=(o.annotation or o.rate.annotation), 
                                               reg=r', '.join(sim.regionList[rgidx]._html() 
                                                                 for rxidx, rgidx in sim._reactionLocations if rxidx==o.idx),
                                               rate=_num(o.rate.value,  unit=Dsp.unicodeUnits(o.rate._unit)))
                                           for o in newRxns])

                    regsCtx = dict(annotationCol=any(o.annotation is not None for o in newRegs),
                                   objs=[dict(id=o.idx, html=o._repr_html_(), annotation=o.annotation)
                                            for o in newRegs])

                    ksCtx = dict(annotationCol=any(o.annotation is not None for o in newKs),
                                 texCol=any(o._texstr is not None for o in newKs),
                                 objs=[dict(id=o.idx, html=o._repr_html_(), annotation=o.annotation, tex=o._texstr,
                                             rate=_num(o.value,  unit=Dsp.unicodeUnits(o._unit)))
                                           for o in newKs])

                    dsCtx = dict(annotationCol=any(o.annotation is not None for o in newDs),
                                 texCol=any(o._texstr is not None for o in newDs),
                                 objs=[dict(id=o.idx, html=o._repr_html_(), annotation=o.annotation, tex=o._texstr,
                                             rate=_num(o.value,  unit=Dsp.unicodeUnits(o._unit)))
                                           for o in newDs])
                    ctx = dict(sps=spsCtx, rxns=rxnsCtx, regs=regsCtx, ks=ksCtx, ds=dsCtx)

                    rows = sum(map(lambda x: len(x['objs']), ctx.values()), 0)

                    if rows > 0:
                        ipd.display(ipd.HTML(Tmp.j2render("newObject.html", ctx)))

        return InteractiveConstruct()


    #### Site lattice viewer

    @_maybeJupyter
    def displayGeometry(self, filterFunctions=None, mode="widget"):
        """3-D site type lattice viewer

        The display mode can be "widget", which displays in the 
        notebook, "download_x3d", which opens a download link in the 
        notebook to the X3D scene, or "download_html", which opens a 
        download link in the notebook to a standalone HTML file. 

        To hide parts of the lattice, `filterFunctions` can be 
        specified. This option takes a list of functions which map 
        from a (x,y,z) mesh grid to a [nx,ny,nz] boolean mask 
        where only subvolumes marked True are shown. To only show
        volumes whose z coordinate are less than 32, the function

        >>> def zfilter(x,y,z):
        >>>     return z<32

        is used. Here the arguments x,y,z are 
        of type :py:class:`numpy.ndarray` and a boolean lattice is
        returned.

        Args:
            filterFunctions (dict):
                Dict of functions which take a (nx,ny,nz) mesh to a 
                bool [nx,ny,nz] mask
            mode (str):
                View mode
        """
        ctx = self._x3dJ2context(filterFunctions)

        if mode == 'widget':
            return Tmp.displayj2html("x3d.html", ctx)
        elif mode == 'download_x3d':
            ctx['downloadX3D'] = True
            xml = Tmp.j2render("structure.x3d", ctx)
            return _downloadFile(xml.encode("ascii"), self.name + ".x3d")
        elif mode == 'download_html':
            data = Tmp.j2render("standaloneX3d.html", ctx)
            return _downloadFile(data.encode("ascii"), self.name + "-3dView.html")
        else:
            raise ValueError("mode should be {widget, download_x3d, download_html}")

    @_maybeJupyter
    def showRegion(self,  plane="xz", planeIndex=None):
        """Display a slice of the site lattice

        Args:
            plane (str):
                Viewing plane, e.g. "xy"
            planeIndex (int):
                Index along the normal direction to the plane
        """
        ctx = self._regionJ2context(plane, planeIndex)
        return Tmp.displayj2html("region.html", ctx)

    @_maybeJupyter
    def showRegionStack(self, plane='xz', scl=None, maxWidth=600, maxHeight=600):
        """Display all slices of the site lattice interactively

        Args:
            plane (str):
                Viewing plane, e.g. "xy"
            scl (int):
                Scale pixels by this amount
            maxWidth (int):
                Maximum width of image
            maxHeight (int):
                Maximum height of image
        """

        siteColors = [r._floatColor() for r in self.regionList]
        htmlNames=[r._html() for r in self.regionList]
        return _showRegionStack(self.siteLattice, htmlNames, siteColors, plane=plane, scl=scl, maxWidth=maxWidth, maxHeight=maxHeight)

    # Model introspection

    @_maybeJupyter
    def showAllParameters(self):
        """Display a table of reaction rates and diffusion constants"""
        return self._paramTable(self, "All defined parameters", list(self.rxnRateList)+list(self.diffRateList))

    @_maybeJupyter
    def showRateConstants(self):
        """Display a table of reaction rates"""
        return self._paramTable(self, "Reaction rate constants", self.rxnRateList)

    @_maybeJupyter
    def showDiffusionConstants(self):
        """Display a table of diffusion constants"""
        return self._paramTable(self, "Diffusion constants", self.diffRateList)

    def _paramTable(self, title, params):
        ctx = dict(title=title, count=len(params), params=[])
        counts = {k:0 for k in params}

        ds = {d.idx: d for d in params if d in self.diffRateList}

        for rx in self.reactionList:
            try:
                counts[rx.rate] += 1
            except KeyError:
                pass

        for _, _, _, r in self._transitionRates:
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
                                     value=_num(p.value),
                                     unit=Dsp.unicodeUnits(p._unit),
                                     count=counts[p])) 
        ctx['annotationCol'] = annotate

        return Tmp.displayj2html("params.html", ctx)


    @_maybeJupyter
    def showSpecies(self,sps):
        """Show details on a species type"""
        if isinstance(sps, str):
            if self.speciesList.is_defined(sps):
                sps = self.species(sps)
            else:
                raise ValueError("unknown species: {}".format(sps))
        cs = self.particleStatistics()
        return Tmp.displayj2html("species.html", self._speciesJ2context(sps,cs))


    @_maybeJupyter
    def showAllSpecies(self):
        """Inspect all species interactively"""
        cs = self.particleStatistics()

        topctx = dict(renderId=_renderId(),
                      spData=[self._speciesJ2context(sps, cs) for sps in self.speciesList])

        return Tmp.displayj2html("allSpecies.html", topctx)

    @_maybeJupyter
    def showReactions(self, rxnList=None):
        """Display a table of all reactions"""
        ctx = self._reactionJ2context(rxnList or self.reactionList)
        return Tmp.displayj2html("reactionTable.html", ctx)



def _renderId():
    return 'jLM_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))


def _downloadFile(data, filename):
    import IPython.display as ipd
    return ipd.Javascript(Tmp.j2render("download.js", dict(filename=filename, blob=base64.b64encode(data).decode())))

def _showRegionStack(lattice, htmlNames, siteColors, plane='xz', scl=None, maxWidth=600, maxHeight=600):
    axis0 = 'xyz'.index(plane[0])
    axis1 = 'xyz'.index(plane[1])
    axis2 = (set((0,1,2)) - set((axis0,axis1))).pop()
    latticeShape = lattice.shape
    siteTypes  = np.unique(lattice)
    ntypes = len(siteTypes)

    if scl is None:
        scl = int(min(maxWidth/latticeShape[axis0], maxHeight/latticeShape[axis1]))
    else:
        scl = int(scl)

    ctx = dict(w=scl*latticeShape[axis0],
               h=scl*latticeShape[axis1],
               plane=plane,
               reg=htmlNames,
               direction="xyz"[axis2],
               dim=latticeShape[axis2],
               start=_maxEntropySlice(lattice, axis2),
               renderId=_renderId())

    pngData=[]
    for planeIndex in range(latticeShape[axis2]):
        stRGB = np.zeros((latticeShape[axis0],latticeShape[axis1],3))

        def rgb8(i):
            return np.array([int(255*c) for c in siteColors[i]], dtype=np.uint8)

        stImg = np.zeros((latticeShape[axis0],latticeShape[axis1]),dtype=np.uint8)

        for i,j in itertools.product(range(latticeShape[axis0]),range(latticeShape[axis1])):
            idx=[0,0,0]
            idx[axis0] = latticeShape[axis0]-i - 1
            idx[axis1] = j
            idx[axis2] = planeIndex
            stImg[i,j] = lattice[idx[0],idx[1],idx[2]]

        scl = int(min(maxWidth/stImg.shape[0], maxHeight/stImg.shape[1]))
        stRGB = np.zeros((stImg.shape[0],stImg.shape[1],3),dtype='uint8')
        for i,j in itertools.product(*map(range,stImg.shape)):
            stRGB[i,j,:] = rgb8(stImg[i,j])

        imgFile = io.BytesIO()
        import PIL.Image
        img = PIL.Image.fromarray(stRGB)
        img.save(imgFile,'PNG')
        pngData.append(base64.b64encode(imgFile.getvalue()).decode())

    ctx['pngData'] = '[' + ','.join('"'+x+'"' for x in pngData) + ']'
    return Tmp.displayj2html("sitestack.html", ctx)


def _showBinaryLattices(binLattices,manualColor=None, filterFunctions=None, mode="widget"):
    assert mode in ['widget', 'download_x3d', 'download_html']

    if isinstance(binLattices, np.ndarray):
        lattices = [("lattice", binLattices)]
    elif isinstance(binLattices, list):
        lattices = [ ("lattice{:02d}".format(d), l) for d,l in enumerate(binLattices) ]
    elif isinstance(binLattices, dict):
        lattices = list(sorted(binLattices.items()))
    else:
        raise TypeError

    latticeDims = lattices[0][1].shape
    scl = 7/max(latticeDims)
    centroid = 0.5*np.array(latticeDims)

    nx,ny,nz = latticeDims

    ctx = dict(centroid=centroid, scl=scl, renderId=_renderId(), sites=[], species=[])

    filterFunctions = filterFunctions or dict()

    for i,(name, binLattice) in enumerate(lattices):
        if manualColor is not None:
            r,g,b = manualColor[name]
            
        else:
            r,g,b = Dsp.colorWheel(i/len(binLattices))
        print(name, r,g,b)
        if binLattice.any():
            if name in filterFunctions:
                posMatch = filterFunctions[name](*np.mgrid[0:nx, 0:ny, 0:nz])
                binLattice = np.array(binLattice&posMatch, dtype=np.uint8)
            else:
                binLattice = np.array(binLattice, dtype=np.uint8)
            verts, faces = Lattice.greedyMesh(binLattice)
            c=Dsp.toHex((r,g,b))
            ctx['sites'].append(dict(name=name,
                                     idx=i,
                                     hexColor=c,
                                     label=r'<span class="jLMregion" style="color:white;background:{};">{}</span>'.format(c, name),
                                     checked=" checked", 
                                     faces=' '.join(str(x) for face in faces for x in face),
                                     verts=' '.join(str(x) for vert in verts for x in vert),
                                     r=r, g=g, b=b,
                                     choice="0"))
    if mode == 'widget':
        return Tmp.displayj2html("x3d.html", ctx)
    elif mode == 'download_x3d':
        ctx['downloadX3D'] = True
        xml = Tmp.j2render("structure.x3d", ctx)
        return _downloadFile(xml.encode("ascii"), "lattice.x3d")
    elif mode == 'download_html':
        data = Tmp.j2render("standaloneX3d.html", ctx)
        return _downloadFile(data.encode("ascii"), "lattice-3dView.html")


def showVolumeStack(vol, plane='xz', cmap='inferno', scl=None, maxWidth=600, maxHeight=600):
    """Display slices volumetric data interactively

    Args:
        vol (:py:class:`numpy.ndarray`):
            3-D data

    Keyword Args:
        cmap (str):
            Name of matplotlib colormap
        plane (str):
            Viewing plane, e.g. "xy"
        scl (int):
            Scale pixels by this amount
        maxWidth (int):
            Maximum width of image
        maxHeight (int):
            Maximum height of image
    """
    import matplotlib.pyplot as plt
    import PIL.Image
    axis0 = 'xyz'.index(plane[0])
    axis1 = 'xyz'.index(plane[1])
    axis2 = (set((0,1,2)) - set((axis0,axis1))).pop()
    vol = vol.transpose([axis0,axis1,axis2])
    shape = vol.shape

    if scl is None:
        scl = int(min(maxWidth/shape[0], maxHeight/shape[1]))
    else:
        scl = int(scl)

    cm = getattr(plt.cm, cmap)
    vmin = vol.min()
    vmax = vol.max()
    cbw = 16
    cbh = int(scl*shape[0]*0.6666)


    cbimg = (255*cm(np.array(cbw*[np.linspace(1,0,cbh)]).T)).astype(np.uint8)
    imgFile = io.BytesIO()
    img = PIL.Image.fromarray(cbimg)
    img.save(imgFile,'PNG')
    cbimg = base64.b64encode(imgFile.getvalue()).decode()

    ctx = dict(h=scl*shape[0],
               w=scl*shape[1],
               cbw=cbw,
               cbh=cbh,
               plane=plane,
               vmax=_num(vmax),
               vmin=_num(vmin),
               direction="xyz"[axis2],
               dim=shape[2],
               cbimg=cbimg,
               start=shape[2]//2,
               renderId=_renderId())

    sl = [slice(None),slice(None),slice(None)]

    pngData=[]
    for planeIndex in range(shape[2]):
        sl[2] = planeIndex
        im = (vol[sl[0],sl[1],sl[2]]-vmin)/(vmax-vmin)
        rgb = (255*cm(im)).astype(np.uint8)

        imgFile = io.BytesIO()
        img = PIL.Image.fromarray(rgb)
        img.save(imgFile,'PNG')
        pngData.append(base64.b64encode(imgFile.getvalue()).decode())

    ctx['pngData'] = '[' + ','.join('"'+x+'"' for x in pngData) + ']'
    return Tmp.displayj2html("datastack.html", ctx)


class _Report:
    def __init__(self):
        self._inctx = False

    def __call__(self, l, v, fmt=None):
        assert self._inctx
        if fmt:
            v = fmt.format(v)
        else:
            v = str(v)
        self._rows.append((l,v))

    def __enter__(self):
        self._inctx = True
        self._rows = []

    def __exit__(self, type, value, traceback):
        self._inctx = False
        html = ("<table class='jLMtbl'>" +
                "".join("<tr><th>{}</th><td>{}</td></tr>".format(l,v) for l,v in self._rows) +
                "</table>"
                )
        import IPython.display as ipd
        ipd.display(ipd.HTML(html))

report = _Report()