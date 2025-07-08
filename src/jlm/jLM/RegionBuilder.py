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

"""Design lattice geometry"""

import numpy as np
import scipy.ndimage as ndi
import scipy.spatial as spspat
from . ColorGen import ColorSeq

from . import RDME
from . import JupyterDisplay as JD


class RegionBuilder:
    """Helper object to design site lattice geometry"""
    def __init__(self, net=None, dims=None):
        """Helper object to design site lattice geometry

        Args:
            net (:py:class:`~jLM.RDME.SpatialModel`): 
                If present, take dimensions from simulation 
            dims ([int]): 
                Dimensions of lattice
        """
        self._net = net

        if not isinstance(net, RDME.SpatialModel) and net is not None:
            raise TypeError("net must be a SpatialModel")

        if dims is not None and not hasattr(dims, "__getitem__") and not hasattr(dims, "__iter__"):
            raise TypeError("dims must be list-like")

        if (dims is None) == (net is None):
            raise ValueError("must specify dims xor net")

        if dims is None:
            self.nx, self.ny, self.nz = net.siteLattice.shape
            self.shape = net.siteLattice.shape
        else:
            self.nx, self.ny, self.nz = dims
            self.shape = dims

        self.xs = np.mgrid[0:self.nx, 0:self.ny, 0:self.nz] #: Index grid of dimensions [3, nx, ny, nz]
        self.center =np.array( [self.nx//2, self.ny//2, self.nz//2])
        self.origin =np.zeros( 3 )

    def emptyLatticeMask(self):
        """Return an empty lattice mask"""
        return np.zeros(self.shape, dtype=bool)

    @staticmethod
    @JD._maybeJupyter
    def showBinaryLattices(binLattices,manualColor=None, filterFunctions=None, mode="widget"):
        """3-D lattice mask viewer

        Lattices can be given as a single binary mask, a list of masks, or a 
        dict of masks.  The display mode can be "widget", which displays in the 
        notebook, "download_x3d", which opens a download link in the notebook to
        the X3D scene, or "download_html", which opens a download link in the 
        notebook to a standalone HTML file. 

        Partial lattices can be shown by applying a filter function.  To hide 
        parts of the lattice, `filterFunctions` can be specified.  This option 
        takes a list of functions which map from a (x,y,z) mesh grid to 
        a [nx,ny,nz] boolean mask where only subvolumes marked True are shown. 
        To only show volumes whose z coordinate are less than 32, the function
        
        >>> def zfilter(x,y,z):
        >>>     return z<32

        is used. Here the arguments x,y,z are of type :py:class:`numpy.ndarray` 
        and a boolean lattice is returned.  The filter functions are specified 
        using a dictionary where the keys correspond to the names of the 
        lattices given. These names come from the dictionary if a dict of 
        lattices is given, `latticeXY` if a list of lattices is given (XY
        being a zero-padded index into the list), or `lattice` if a single
        mask is given.

        Args:
            binLattices:
                Binary lattice or collection of binary lattices
            filterFunctions (dict(str=func)):
                Dict of functions indexed by name of lattice
            mode (str):
                View mode
        """
        return JD._showBinaryLattices(binLattices, manualColor, filterFunctions, mode)

    @staticmethod
    @JD._maybeJupyter
    def showStack(binLattices, plane='xz', scl=None, maxWidth=600, maxHeight=600):
        """Display all slices of the site lattice interactively

        Args:
            binLattices:
                Binary lattice or collection of binary lattices
                Viewing plane, e.g. "xy"
            plane (str):
                Plane to display: (xy, yz, xz)
            scl (int):
                Scale pixels by this amount
            maxWidth (int):
                Maximum width of image
            maxHeight (int):
                Maximum height of image
        """

        if isinstance(binLattices, np.ndarray):
            lattices = [("lattice", binLattices)]
        elif isinstance(binLattices, list):
            lattices = [ ("lattice{:02d}".format(d), l) for d,l in enumerate(binLattices) ]
        elif isinstance(binLattices, dict):
            lattices = list(sorted(binLattices.items()))
        else:
            raise TypeError

        lattice = np.zeros_like(lattices[0][1])
        htmlNames = [r'<span class="jLMregion" style="color:black;background:white;">BG</span>']
        siteColors = [(1,1,1)]
        for i,(n,l) in enumerate(lattices):
            lattice[l] = i+1
            args = dict(fg=ColorSeq.strFg(i+1),
                        bg=ColorSeq.str(i+1),
                        n=n)
            htmlNames.append(r'<span class="jLMregion" style="color:{fg};background:{bg};">{n}</span>'.format(**args))
            siteColors.append(ColorSeq.float(i+1))

        return JD._showRegionStack(lattice, htmlNames, siteColors, plane=plane, scl=scl, maxWidth=maxWidth, maxHeight=maxHeight)

    @staticmethod
    def transformGrid(xs, x0, alpha, beta, gamma):
        """Compute the translation/rotation of an index grid

        Args:
            xs (:py:class:`~numpy.ndarray(shape=(3,nx,ny,nz), dtype=int)`): 
                Index grid of dimensions [3, nx, ny, nz]
            x0 (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                Center of rotation
            alpha (float):
                Euler Z-rotation (radians)
            beta (float):
                Euler X-rotation (radians)
            gamma (float):
                Euler Z'-rotation (radians)

        Returns:
            :py:class:`~numpy.ndarray(shape=(3,nx,ny,nz), dtype=int)`: 
                Transformed grid
        """
        x0 = np.array(x0)
        ca, sa = np.cos(alpha), np.sin(alpha)
        cb, sb = np.cos(beta), np.sin(beta)
        cg, sg = np.cos(gamma), np.sin(gamma)
        rotZ0 = np.array([ [ca,-sa,0],[sa,ca,0], [0,0,1] ])
        rotX= np.array([[1,0,0,], [0,cb,-sb],[0,sb,cb] ])
        rotZ1 =np.array( [ [cg,-sg,0],[sg,cg,0], [0,0,1] ])
        rot = rotZ1 @ rotX @ rotZ0
        return  np.einsum("ijkl,mi->mjkl",xs  - x0[:,np.newaxis,np.newaxis,np.newaxis], rot)

    def _parseArgs(self, angles, center, xs):
        if xs is None:
            xs = self.xs
        if angles is None:
            angles = [0,0,0]
        elif np.isscalar(angles):
            raise TypeError("Need angle=(alpha, beta, gamma). ZXZ euler angles.")
        if center is None:
            center = self.center
        elif np.isscalar(center):
            raise TypeError("Need center=(x,y,z).")
        return angles, center, xs

    def ellipsoid(self, xs=None, radius=None, angles=None, center=None):
        """Construct ellipsoid mask

        Args:
            radius (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`): 
                Semiaxes of ellipse
            center (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                Centroid
            angles (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):   
                [alpha, beta, gamma] ZXZ' Euler angles
            xs (:py:class:`~numpy.ndarray(shape=(3,nx,ny,nz), dtype=int)`): 
                Index grid of dimensions [3, nx, ny, nz]

        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Lattice mask
        """
        angles, center, xs = self._parseArgs(angles, center, xs)
        if np.isscalar(radius):
            radius = [radius, radius, radius]
        elif radius is None:
            raise RuntimeError("must define radius")

        invRs =  np.array([1/x for x in radius])
        xs1 = self.transformGrid(xs,center,*angles)*invRs[:,np.newaxis,np.newaxis,np.newaxis]
        return np.sum(xs1*xs1,axis=0) <= 1

    def cylinder(self, radius, length, xs=None, angles=None, center=None):
        """Construct cylinder mask

        Args:
            radius (float):
                Radius
            length (float):
                Length
            center (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                Centroid
            angles (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                [alpha, beta, gamma] ZXZ' Euler angles
            xs (:py:class:`~numpy.ndarray(shape=(3,nx,ny,nz), dtype=int)`): 
                Index grid of dimensions [3, nx, ny, nz]

        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Lattice mask
        """
        angles, center, xs = self._parseArgs(angles, center, xs)
        xs1 = self.transformGrid(xs,center,*angles)
        m1 = np.sum(xs1[:2,...]*xs1[:2,...],axis=0) <= radius**2
        m2 = xs1[2,...] <= 0.5*length
        m3 = xs1[2,...] >= -0.5*length
        return m1&m2&m3

    def spoke(self, x0, length, spoke_radius, r, phi, theta):
        r"""Construct spoke

        For a sphere of radius, :math:`r`, centered on :math:`x_0`, a spoke is
        a cylinder centered on the surface of the sphere at the polar
        coordinates :math:`(r,\phi,\theta)`, rotated to be normal to the sphere
        surface.

        Args:
            x0 (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                Center of sphere
            length (float):
                Length of spoke. Spoke will protrude 0.5*length from the inside and outside of the sphere
            spoke_radius (float):
                Radius (thickness) of the spoke
            r (float):
                Radius of sphere
            phi (float):
                Azimuthal position of spoke with respect to the sphere center
            theta (float):
                Polar position of spoke with respect to the sphere center

        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Lattice mask
        """

        dv = np.array([
            r*np.cos(phi)*np.sin(theta),
            r*np.sin(phi)*np.sin(theta),
            r*np.cos(theta)])

        normal = dv/np.sqrt(dv@dv)
        alpha = np.arctan2(normal[0], normal[1])
        beta = np.arccos(normal[2])
        gamma = 0
        return self.cylinder(spoke_radius, length, 
                             angles=[alpha,beta,gamma],
                             center=r*normal + x0)

    def capsule(self, length, width, xs=None, angles=None, center=None):
        """Construct spherocylinder mask

        Args:
            length (float):
                Length including endcaps
            width (float):
                Diameter of cylindrical region
            center (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                Centroid
            angles (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                [alpha, beta, gamma] ZXZ' Euler angles
            xs (:py:class:`~numpy.ndarray(shape=(3,nx,ny,nz), dtype=int)`): 
                Index grid of dimensions [3, nx, ny, nz]

        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Lattice mask
        """
        angles, center, xs = self._parseArgs(angles, center, xs)
        r = width/2
        h = (length-width)/2
        hv = np.array([0,0,h])
        xs1 = self.transformGrid(xs,center,*angles)
        return (self.cylinder(r, 2*h, xs=xs1, center=self.origin) 
                | self.ellipsoid(radius=r, center=self.origin+hv, xs=xs1) 
                | self.ellipsoid(radius=r, center=self.origin-hv, xs=xs1))

    def box(self, lx, ly, lz, xs=None, angles=None, center=None):
        """Construct a rectangular cuboid mask

        Args:
            lx (float):
                Length in x-direction
            ly (float):
                Length in y-direction
            lz (float):
                Length in z-direction
            center (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                Centroid
            angles (:py:class:`~numpy.ndarray(shape=(3,), dtype=float)`):
                [alpha, beta, gamma] ZXZ Euler angles
            xs (:py:class:`~numpy.ndarray(shape=(3,nx,ny,nz), dtype=int)`):
                Index grid of dimensions [3, nx, ny, nz]

        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`:
                Lattice mask
        """
        angles, center, xs = self._parseArgs(angles, center, xs)
        xs1 = self.transformGrid(xs,center,*angles)
        return (xs1[0] < lx) & (xs1[0] >= 0) & (xs1[1] < ly) & (xs1[1] >= 0) & (xs1[2] < lz) & (xs1[2] >= 0)
    
    def compose(self, *siteSpec, net=None):
        """Compose a series of binary masks into a site lattice

        This function takes an indefinite number of (region, mask) tuples. The 
        lattice is created by setting the lattice to the index of region over 
        all masked subvolumes in the order given.  

        Args:
            *siteSpec ((:py:class:`~jLM.Types.Region`, :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`)): 
                Region, mask tuple
            net (:py:class:`~jLM.RDME.SpatialModel`): 
                If given, modify this model's site lattice
        """
        if net is not None:
            l = net.siteLattice
        elif self._net is not None:
            l = self._net.siteLattice
        else:
            raise RuntimeError("must specify SpatialModel")

        for (reg, binaryMask) in siteSpec:
            l[binaryMask] = reg.idx

    @classmethod
    def _morphApply(cls, bi, radius, se, fn):
        if se is None and radius is not None:
            se = cls.octoStructElem(radius)
        elif se is None and radius is None:
            raise ValueError("specify radius or structuring element")

        r = np.max(se.shape)//2 + 1
        pbi = np.pad(bi, r, 'constant', constant_values=0)
        mbi = fn(pbi, se)
        obi = mbi[r:-r, r:-r, r:-r]
        return obi

    @classmethod
    def _morphApplyTwice(cls, bi, radius1, se1, fn1, radius2, se2, fn2):
        if se1 is None and radius1 is not None:
            se1 = cls.octoStructElem(radius1)
        elif se1 is None and radius1 is None:
            raise ValueError("specify radius1 or structuring element1")
        if se2 is None and radius2 is not None:
            se2 = cls.octoStructElem(radius2)
        elif se2 is None and radius2 is None:
            raise ValueError("specify radius2 or structuring element2")

        r = max(np.max(se2.shape),np.max(se1.shape))//2 + 1
        pbi = np.pad(bi, r, 'constant', constant_values=0)
        mbi = fn2(fn1(pbi, se1), se2)
        obi = mbi[r:-r, r:-r, r:-r]
        return obi


    @staticmethod
    def octoStructElem(r):
        """Iterated 6-connected structuring element

        Args:
            r (int):
                Number of iterations

        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Structuring element
        """
        return ndi.iterate_structure(ndi.generate_binary_structure(3,1),int(r))

    @property
    def se6(self):
        """Structuring element connecting all 6 nearest neighbors"""
        se = np.zeros((3,3,3), dtype=bool)
        se[:,1,1] = True
        se[1,:,1] = True
        se[1,1,:] = True
        return se

    @property
    def se26(self):
        """Structuring element connecting all 26 neighbors"""
        se = np.zeros((3,3,3), dtype=bool)
        se[...] = True
        return se

    @staticmethod
    def sphereStructElem(r):
        """Spherical structuring element

        Args:
            r (int):
                radius

        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Structuring element
        """
        n = 2*r - 1
        x,y,z = np.mgrid[:n,:n,:n]
        return (x-r+1)**2 + (y-r+1)**2+(z-r+1)**2 < r**2

    @classmethod
    def dilate(cls, binaryMask, radius=None, se=None):
        """Morphological dialation of a binary mask

        Args:
            binaryMask (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Binary mask
            radius (int):   
                If provided use a 6-connected structuring element iterated `radius` times.
            se (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Structuring element
        
        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Dilated lattice
        """
        return cls._morphApply(binaryMask, radius, se, ndi.binary_dilation)

    @classmethod
    def erode(cls, binaryMask, radius=None, se=None):
        """Morphological erosion of a binary mask

        Args:
            binaryMask (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Binary mask
            radius (int):
                If provided use a 6-connected structuring element iterated 
                `radius` times.
            se (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Structuring element
        
        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Eroded lattice
        """
        return cls._morphApply(binaryMask, radius, se, ndi.binary_erosion)

    @classmethod
    def closing(cls, binaryMask, radius=None, se=None, radius1=None, se1=None):
        """Morphological closing of a binary mask

        Args:
            binaryMask (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Binary mask
            radius (int):
                If provided use a 6-connected structuring element iterated 
                `radius` times.
            radius1 (int):
                If provided use the iterated SE for the erosion (optional)
            se (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Structuring element
            se1 (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Structuring element for erosion (optional)
        
        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Closed lattice
        """
        if radius1 is not None or se1 is not None:
            return cls._morphApplyTwice(binaryMask, radius, se, ndi.binary_dilation, radius1, se1, ndi.binary_erosion)
        else:
            return cls._morphApply(binaryMask, radius, se, ndi.binary_closing)

    @classmethod
    def opening(cls, binaryMask, radius=None, se=None, radius1=None, se1=None):
        """Morphological opening of a binary mask

        Args:
            binaryMask (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Binary mask
            radius (int):   
                If provided use a 6-connected structuring element iterated 
                `radius` times.
            radius1 (int):  
                If provided use the iterated SE for the dilation (optional)
            se (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Structuring element
            se1 (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Structuring element for dilation (optional)
        
        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`:
                Opened lattice
        """
        if radius1 is not None or se1 is not None:
            return cls._morphApplyTwice(binaryMask, radius, se, ndi.binary_erosion, radius1, se1, ndi.binary_dilation)
        else:
            return cls._morphApply(binaryMask, radius, se, ndi.binary_opening)

    @classmethod
    def convexHull(cls, binaryMask):
        """Convex hull of lattice

        Args:
            binaryMask (:py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`):
                Binary mask
        
        Returns:
            :py:class:`~numpy.ndarray(shape=(nx,ny,nz), dtype=bool)`: 
                Convex hull of lattice
        """
        coords = np.argwhere(binaryMask)
        hull = spspat.ConvexHull(coords)
        triangulation = spspat.Delaunay(coords[hull.vertices])
        newMask = np.zeros_like(binaryMask)
        newMaskRavel = newMask.ravel()
        evalCoords = np.array(np.unravel_index(np.arange(binaryMask.size), binaryMask.shape)).T
        newMaskRavel[...] = triangulation.find_simplex(evalCoords) >= 0
        return newMask