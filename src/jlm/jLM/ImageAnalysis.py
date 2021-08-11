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

import numpy as np
import cairocffi as cairo
import io
import IPython.display as ipd
import itertools
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.pyplot as plt
import PIL.Image
import webcolors


def _makeMaskPath(binarySlice):
    # collect edges of each lattice site
    # only keep boundary edges, if more than one edge traverses two points, it is internal

    if not np.any(binarySlice):
        return (), ()

    faces = dict()
    def ap(*v):
        k = tuple(int(x) for x in v)
        if k in faces:
            faces[k] = None
        else:
            faces[k] = 1

    for x,y in itertools.product(*map(range,binarySlice.shape)):
        if binarySlice[x,y]:
            ap( x,y, x+1,y )
            ap( x,y, x,y+1 )
            ap( x,y+1, x+1,y+1 )
            ap( x+1,y, x+1,y+1 )


    faces2 = set((f[:2],f[2:]) for f,v in faces.items() if v is not None )

    # build map for looking up an edge based on its verticies
    vertexTable = { x:[] for x in  set(x for x,_ in faces2) | set( x for _,x in faces2) }

    for l,r in faces2:
        vertexTable[l].append((l,r))
        vertexTable[r].append((l,r))

    # merge edges into single closed path
    face = faces2.pop()
    paths = [ [ face[0], face[1] ] ]

    def faceEq(f0,f1):
        return (f0[0] == f1[0] and f0[1] == f1[1]) or (f0[1] == f1[0] and f0[0] == f1[1]) 

    while len(faces2) > 0:
        firstFace = paths[0][-2:]
        lastFace = paths[-1][-2:]
        pl, pr = lastFace
        for face in vertexTable[pr]:
            # append next vertex, making sure it's not already visited
            if not faceEq(face, lastFace) and not faceEq(face, firstFace):
                try:
                    vertexTable[pr].remove(face)
                    faces2.remove(face)
                except KeyError:
                    break

                if pr == face[0]:
                    paths[-1].append(face[1])
                else:
                    paths[-1].append(face[0])

                break
        else: 
            # start a new path if we ran out of edges to connect
            if len(faces2) > 0:
                face = faces2.pop()
                paths.append([ face[0], face[1] ])


    newPaths = []
    for path in paths:
        # set consistent winding order
        p0,p1 = path[0:2]
        if p0[0] == p1[0]:
            if p0[1] > p1[1]:
                path = path[::-1]
        else:
            if p0[0] > p1[0]:
                path = path[::-1]

        newPaths.append(path)

    # save path in format mpl can deal with

    verts, codes = [], []

    # reverse sort by arc length
    def pathKey(x):
        return -sum((x0-x1)**2+(y0-y1)**2 for (x0,x1),(y0,y1) in zip(x,[x[-1]] + x[:-1]) )

    for pidx,path in enumerate(sorted(newPaths,key=pathKey)):
        vs, cs = [], []
        for i,vertex in enumerate(path):
            vs.append( (vertex[1], vertex[0]) )
            if i == 0:
                cs.append( "MOVETO" )
            else:
                cs.append( "LINETO" )
        
        # assume longest path is external edge and reverse so that following paths subtract fill correctly
        # if there are "islands" inside cavities, they will have the wrong winding order.
        if pidx == 0:
            vs = vs[::-1]

        vs.append( (0,0) )
        cs.append( "CLOSEPOLY" )
        verts += vs
        codes += cs
    return verts, codes


def mplBoundaryPath(binarySlice, faceColor='none', edgeColor='r', extent=None):
    """Generate a Matplotlib patch defining the boundaries of a binary image.

    Arguments:
        binarySlice:
            Binary 2D image
    Keyword Arguments:
        faceColor:
            Matplotlib color of interior (True)
        edgeColor:
            Matplotlib color of boundary
    Returns:
        matplotlib.patches.PathPatch: patches
    """
    verts, codes = _makeMaskPath(binarySlice)
    if extent:
        x0,x1,y0,y1 = extent
        dx = (x1-x0)/binarySlice.shape[0]
        dy = (y1-y0)/binarySlice.shape[1]
        v2 = []
        for x,y in verts:
            v2.append((dx*x+x0, dy*y+y0))
        verts = v2

    codes = [getattr(Path, x) for x in codes]
    return patches.PathPatch(Path(verts,codes), fc=faceColor, ec=edgeColor)



def _processImgToDisplay(img, cmap, scl, vmin, vmax, mode='rgba'):
    img = img.astype(float)
    img = (img-vmin)/(vmax-vmin)
    cm = getattr(plt.cm, cmap)
    if img.ndim == 3:
        if img.shape[2] == 4:
            rgb = img
        elif img.shape[2] == 3:
            rgb = np.zeros((img.shape[0], img.shape[1], 4))
            rgb[:,:,:3] = img
            rgb[:,:,3] = 1.0
        else:
            raise ValueError("Wrong number of channels")
    elif img.ndim == 2:
        rgb = cm(img)
    else:
        raise ValueError("Wrong dimensions for image")
        
    rgb = (255*rgb).astype(np.uint8)
    if scl is not None:
        img = np.zeros( (rgb.shape[0]*scl, rgb.shape[1]*scl, rgb.shape[2]), dtype=rgb.dtype)
        for i,j in itertools.product(range(scl), repeat=2):
            img[i::scl,j::scl, :] = rgb[...]
        rgb = img
    if mode=="bgra":
        return np.ascontiguousarray(rgb[:,:,(2,1,0,3)])
    else:
        return rgb


def imshow(img, scl=None, cmap="inferno", vmin=None, vmax=None, xlines=[], ylines=[]):
    """Quick display of array in notebook

    Arguments:
        img (array):
            A grayscale image of shape (M,N) or a rgb image of shape (M,N,3)

    Keyword Arguments:
        scl (int):
            Factor to upscale image

        cmap (str):
            Matplotlib colormap

        vmin (float):
            Minimum value

        vmax (float):
            Maximum value

        xlines (list):
            List of lines of constant x

        ylines (list):
            List of lines of constant y

    Returns:
       PIL.Image: PIL image which is displayed in Jupyter
    """
    scl = scl or 1
    vmin = vmin or img.min()
    vmax = vmax or img.max()
    rgb = _processImgToDisplay(img, cmap, scl, vmin, vmax)
    o = scl//2
    for x in xlines:
        rgb[o+x*scl,:,:] = 255
    for y in ylines:
        rgb[:,o+y*scl,:] = 255
    if rgb.shape[-1] == 4:
        return PIL.Image.fromarray(rgb,mode="RGBA")
    else:
        return PIL.Image.fromarray(rgb)


class SvgContext:
    def __init__(self, nx,ny, filename=None):
        """Context manager to draw to an SVG surface.

        Arguments:
            nx(int):
              Width in pixels
            ny(int):
              Height in pixels
        Keyword Arguments:
            filename (str):
                Write SVG to file as well as displaying in notebook
        """
        self.w = nx
        self.h = ny
        self.fileObj = io.BytesIO()
        self.surface = cairo.SVGSurface(self.fileObj,nx,ny)
        self.context = cairo.Context(self.surface)
        self.filename = filename

    def __enter__(self):
        return self.context

    def __exit__(self, type, value, traceback):
        if type is None:
            self.context.show_page()
            self.surface.finish()
            svg = self.fileObj.getvalue()
            ipd.display(ipd.SVG(svg))
            if self.filename is not None:
                open(self.filename, "wb").write(svg)



def alphaBlend(rgbaBottom, rgbaTop):
    """Alpha blend two rgba images

    Arguments:
        rgbaBottom:
            Image of dimensions (N,M,4) to use as base
        rgbaTop:
            Image of dimensions (N,M,4) to blend atop rgbaBottom
    Returns:
        Blended rgba image
    """
    rB,gB,bB,aB = np.transpose(rgbaBottom,axes=(2,0,1))
    rT,gT,bT,aT = np.transpose(rgbaTop,axes=(2,0,1))
    a1 = aT+aB*(1-aT)
    def blend(x,ax,y,ay,a):
        return np.where(a>0,(x*ax+y*ay*(1-ax))/a, 0)
    r1 = blend(rT,aT,rB,aB,a1)
    g1 = blend(gT,aT,gB,aB,a1)
    b1 = blend(bT,aT,bB,aB,a1)
    return np.transpose(np.array( [r1,g1,b1,a1] ),axes=[1,2,0])

def _interpretColor(c):
    if isinstance(c, str):
        if c[0] == '#':
            return webcolors.hex_to_rgb(c)
        else:
            return webcolors.name_to_rgb(c)
    elif hasattr(c, "__len__") and len(c) == 3:
        return c
    else:
        raise ValueError("Invalid color data")


def markup(imgIn, scl=4, cmap="inferno", 
           vmin=None,
           vmax=None,
           markerColor="blue", 
           markerWeight=1, 
           markerSize=2, 
           markers=None,
           mask=None, 
           maskColor="green", 
           maskWeight=1,
           svgFile=None):
    """Markup an image
    
    Markup generates nearest-neighbor interpolated images with binary mask outlines,
    and labeled points. Colors are specified as an RGB triple, a hex code, or a named
    color.
    
    Arguments:
        imgIn:
            Image data, either mono (N,M), rgb (N,M,3), or rgba (N,M,4)
    Keyword arguments:
        scl (int):
            Upscale image by factor
        cmap (str):
            Matplotlib colormap name:
        vmin (float):
            Minimum value
        vmax (float):
            Maximum value
        markers:
            Markers to place. Either a single x,y pair, a list of x,y pairs, or a dictionary
            where keys are the marker color and the values are the x,y pair
        markerColor:
            Color of + markers.
        markerWeight:
            Linewidth of + markers in pixels
        markerSize:
            Diameter of + markers
        mask:
            Binary mask(s) to trace. Either a single numpy array, or a dictionary with color 
            keys mapping to binary arrays. The outline is placed around the True elements
        maskColor:
            Color of mask outline
        maskWeight:
            Linewidth of mask outline in pixels
        svgFile:
            Also write the SVG output to this file
            
    Returns:
        IPython.core.display.SVG
    """
    vmin = vmin or imgIn.min()
    vmax = vmax or imgIn.max()
    rgb = _processImgToDisplay(imgIn, cmap, scl, vmin, vmax, mode="bgra")
    
    nx = imgIn.shape[1]
    ny = imgIn.shape[0]
    
    canvasX = nx*scl
    canvasY = ny*scl

    fileObj = io.BytesIO()
    surf = cairo.SVGSurface(fileObj,canvasX,canvasY)
    ctx = cairo.Context(surf)
    img = cairo.ImageSurface.create_for_data(memoryview(rgb.ravel()),cairo.FORMAT_ARGB32, canvasX, canvasY)
    ctx.set_source_surface(img)
    ctx.paint()
    
    if mask is not None:
        ctx.set_line_width(maskWeight)
        
        def mark(mask):
            for (x,y), code in zip(*_makeMaskPath(mask)):
                if code == "MOVETO":
                    ctx.move_to(scl*x,scl*y)
                elif code == "LINETO":
                    ctx.line_to(scl*x,scl*y)
                elif code == "CLOSEPOLY":
                    ctx.close_path()
                else:
                    raise RuntimeError("Bad path code")
            ctx.stroke()
            
        if isinstance(mask, dict):
            for color, mk in mask.items():
                ctx.set_source_rgb(*_interpretColor(color))
                mark(mk)
        else:
            ctx.set_source_rgb(*_interpretColor(maskColor))
            mark(mask)
            
        
    if markers is not None:
        ctx.set_line_width(markerWeight)
        def mark(x,y):
            x, y = (x+0.5)*scl, (y+0.5)*scl
            x0, y0 = (x-0.5*markerSize*scl), (y-0.5*markerSize*scl)
            x1, y1 = (x+0.5*markerSize*scl), (y+0.5*markerSize*scl)
            ctx.move_to(x, y0)
            ctx.line_to(x, y1)
            ctx.move_to(x0, y)
            ctx.line_to(x1, y)
            ctx.stroke()
            
        if isinstance(markers, dict):
            for color, pos in markers.items():
                ctx.set_source_rgb(*_interpretColor(color))
                for x,y in pos:
                    mark(x,y)
        elif hasattr(markers, "__len__"):
            ctx.set_source_rgb(*_interpretColor(markerColor))
            if hasattr(markers[0], "__len__"):
                for x,y in markers:
                    mark(x,y)
            else:
                mark(*markers)
        else:
            raise ValueError("Invalid marker data")
            
    ctx.show_page()
    surf.finish()
    v = fileObj.getvalue()
    
    if svgFile:
        with open(svgFile, "wb") as f:
            f.write(v)
    
    return ipd.SVG(v)
