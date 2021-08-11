/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 * of the Software, and to permit persons to whom the Software is furnished to 
 * do so, subject to the following conditions:
 * 
 * - Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimers.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimers in the documentation 
 * and/or other materials provided with the distribution.
 * 
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Tyler M. Earnest, Mike Hallock, Elijah Roberts, 
 *            Joseph R. Peterson
 */


%module(directors="1",threads="1", docstring="Lattice Microbes Python Bindings") lm
%include "documentation.i"
%include "typemaps.i"
%include "std_string.i"
%include "std_vector.i"
%feature("director:except") {
    if( $error != NULL ) {
        PyObject *ptype, *pvalue, *ptraceback;
        PyErr_Fetch( &ptype, &pvalue, &ptraceback );
        PyErr_Restore( ptype, pvalue, ptraceback );
        PyErr_Print();
        Py_Exit(1);
    }
}

%include "config.h" 

%{
#define SWIG_FILE_WITH_INIT
#include <vector>
#include "config.h" 
#include "core/Print.h" 
#include "core/Types.h" 
#include "core/Exceptions.h"
#include "builder/Capsule.h"
#include "builder/CapsuleShell.h"
#include "builder/Cuboid.h"
#include "builder/Hemisphere.h"
#include "builder/Shape.h"
#include "builder/Sphere.h" 
#include "builder/LatticeBuilder.h"
#include "builder/Shape.h"
#include "builder/Sphere.h"
#include "builder/Union.h"
#include "builder/UnionSet.h"
#include "builder/Difference.h"
#include "builder/Torus.h"
#include "builder/Ellipse.h"
#include "builder/Cylinder.h"
#include "builder/Cone.h"
#include "builder/Intersection.h"
#include "DiffusionModel.pb.h"
#include "ReactionModel.pb.h"
#include "SpatialModel.pb.h"
#include "io/SimulationFile.h"
#include "rdme/ByteLattice.h"
#include "rdme/IntLattice.h"
#ifdef OPT_CUDA
#include "rdme/CudaByteLattice.h"
#include "rdme/CudaIntLattice.h"
#endif
#include "rdme/Lattice.h"
#include "core/runSimulation.h"

#include "me/MESolver.h"
#include "cme/CMESolver.h"
#include "cme/GillespieDSolver.h"
#include "cme/NextReactionSolver.h"
#include "rdme/RDMESolver.h"
#ifdef OPT_CUDA
#include "rdme/MGPUMpdRdmeSolver.h"
#include "rdme/MpdRdmeSolver.h"
#include "rdme/IntMpdRdmeSolver.h"
#endif

using lm::Exception;
using lm::IOException;
using lm::InvalidArgException;
using lm::io::hdf5::HDF5Exception;
using lm::io::SpatialModel;
using lm::rdme::InvalidSiteException;
using lm::rdme::InvalidParticleException;
%}

%inline %{
typedef long time_t;
typedef unsigned char       uchar;
typedef unsigned int        uint;
typedef unsigned long       ulong;
typedef long                intv_t;
typedef unsigned long       uintv_t;
typedef unsigned char       uint8_t;
typedef unsigned int        uint32_t;
typedef uint8_t             byte;
typedef double              si_dist_t;
typedef double              si_time_t;
typedef uint32_t            lattice_size_t;
typedef uint32_t            site_size_t;
typedef uint32_t            site_t;
typedef uint32_t            particle_t;
void setVerbosityLevel(int x) { lm::Print::verbosityLevel(x); }
%}



/*
 * Typemap for returning a list of floats.
 */


%typemap(out) std::vector<double> {

    $result = PyList_New($1.size());
    int index=0;
    for (std::vector<double>::const_iterator it=$1.begin(); it != $1.end(); it++, index++)
    {
        PyList_SetItem($result, index, PyFloat_FromDouble(*it));
    }
}

/*
 * Typemap for returning a list of particle location structures.
 */
%typemap(out) std::vector<particle_loc_t> {

    $result = PyList_New($1.size());
    int index=0;
    for (std::vector<particle_loc_t>::const_iterator it=$1.begin(); it != $1.end(); it++, index++)
    {
        PyObject* resultTuple = PyTuple_New(5);
        PyTuple_SetItem(resultTuple, 0, PyInt_FromLong((*it).p));
        PyTuple_SetItem(resultTuple, 1, PyInt_FromLong((*it).x));
        PyTuple_SetItem(resultTuple, 2, PyInt_FromLong((*it).y));
        PyTuple_SetItem(resultTuple, 3, PyInt_FromLong((*it).z));
        PyTuple_SetItem(resultTuple, 4, PyInt_FromLong((*it).index));
        PyList_SetItem($result, index, resultTuple);
    }
}

/* typecheck to see if object is a tuple, suitable to be lattice_coord_t */
%typemap(typecheck) lattice_coord_t {
    $1 = PyTuple_Check($input) ?
        ((PyTuple_Size($input) == 3) ? 1 : 0) :
        0;
}

/* accept a tuple as a lattice_coord_t */
%typemap(in) lattice_coord_t {
    $1 = lattice_coord_t(
        PyInt_AsLong(PyTuple_GetItem($input, 0)),
        PyInt_AsLong(PyTuple_GetItem($input, 1)),
        PyInt_AsLong(PyTuple_GetItem($input, 2)));
}

/*
 * Type for returning a lattice coordinate.
 */
%typemap(out) lattice_coord_t {
	$result = PyTuple_New(3);
	PyTuple_SetItem($result, 0, PyInt_FromLong($1.x));
	PyTuple_SetItem($result, 1, PyInt_FromLong($1.y));
	PyTuple_SetItem($result, 2, PyInt_FromLong($1.z));
}

/*
 * Type for returning particles from CUDALattice::findAllParticles.
 */
%typemap(out) std::map<particle_t,std::vector<lattice_coord_t> > {

	$result = PyDict_New();
	for (std::map<particle_t,std::vector<lattice_coord_t> >::iterator it=$1.begin(); it != $1.end(); it++)
	{
		particle_t particle = (*it).first;
		std::vector<lattice_coord_t> coordList = (*it).second;
	
		PyObject* resultList = PyList_New(coordList.size());
		int index=0;
		for (std::vector<lattice_coord_t>::const_iterator it2=coordList.begin(); it2 != coordList.end(); it2++, index++)
		{
			lattice_coord_t value = *it2;	
			PyObject* resultTuple = PyTuple_New(3);
			PyTuple_SetItem(resultTuple, 0, PyInt_FromLong(value.x));
			PyTuple_SetItem(resultTuple, 1, PyInt_FromLong(value.y));
			PyTuple_SetItem(resultTuple, 2, PyInt_FromLong(value.z));
			PyList_SetItem(resultList, index, resultTuple);
		}
		
		PyDict_SetItem($result, PyInt_FromLong(particle), resultList);
	}
}

/*
 * Type for passing in a list of particle types.
 */
%typemap(in) std::list<particle_t> {
	if (PyList_Check($input))
	{
		for (int i=0; i<PyList_Size($input); i++)
		{
			PyObject *o = PyList_GetItem($input, i);
			if (PyInt_Check(o))
			{
				$1.push_back((particle_t)PyInt_AsLong(o));
			}
			else
			{
				PyErr_SetString(PyExc_TypeError, "must be a list of particle types");
				return NULL;
			}
		}
	}
	else
	{
		PyErr_SetString(PyExc_TypeError, "must be a list");
		return NULL;
	}	
}

%typemap(typecheck) std::list<particle_t> {
   $1 = PyList_Check($input) ? 1 : 0;
}

/*
 * Type for returning a list of particle counts.
 */
%typemap(out) std::map<particle_t,uintv_t> {
	$result = PyDict_New();
	int index=0;
	for (std::map<particle_t,uintv_t>::const_iterator it=$1.begin(); it != $1.end(); it++, index++)
	{
		PyDict_SetItem($result, PyInt_FromLong((*it).first), PyInt_FromLong((*it).second));
	}
}

%typemap(out) std::map<particle_t,uint> {
	$result = PyDict_New();
	int index=0;
	for (std::map<particle_t,uint>::const_iterator it=$1.begin(); it != $1.end(); it++, index++)
	{
		PyDict_SetItem($result, PyInt_FromLong((*it).first), PyInt_FromLong((*it).second));
	}
}

/*
 * Type for returning a list of particle names.
 */
%typemap(out) std::map<uint,std::string> {
	$result = PyDict_New();
	int index=0;
	for (std::map<uint,std::string>::const_iterator it=$1.begin(); it != $1.end(); it++, index++)
	{
        const char *n=(it->second).c_str();
		PyDict_SetItem($result, PyInt_FromLong((*it).first), PyString_FromString(n));
	}
}

%exception {
    try {
        $action
    }
    catch (HDF5Exception& e) {
        PyErr_SetString(PyExc_OSError, e.what());
        SWIG_fail;
    }
    catch (InvalidSiteException& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        SWIG_fail;
    }
    catch (InvalidParticleException& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        SWIG_fail;
    }
    catch (InvalidArgException& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        SWIG_fail;
    }
    catch (IOException& e) {
        PyErr_SetString(PyExc_OSError, e.what());
        SWIG_fail;
    }
    catch (Exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        SWIG_fail;
    }
}


%include "numpy.i"

%init %{
    import_array();
%}


%template(stdVectorInt) std::vector<int>;

%typemap(in,numinputs=0)
    (uint8_t **siteLattice, int *Nz, int *Ny, int *Nx)
    (uint8_t *data_temp, int dim1_temp, int dim2_temp,  int dim3_temp)
{
    $1 = &data_temp;
    $2 = &dim1_temp;
    $3 = &dim2_temp;
    $4 = &dim3_temp;
}

%typemap(argout)
  (uint8_t **siteLattice, int *Nz, int *Ny, int *Nx)
{
  npy_intp dims[3] = { *$2, *$3, *$4 };
  PyObject* obj = PyArray_SimpleNewFromData(3, dims, NPY_UBYTE, (void*)(*$1));
  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result,obj);
}


%typemap(in,numinputs=0)
    (uint8_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np)
    (uint8_t *data_temp, int dim1_temp, int dim2_temp,  int dim3_temp, int dim4_temp, int dim5_temp )
{
    $1 = &data_temp;
    $2 = &dim1_temp;
    $3 = &dim2_temp;
    $4 = &dim3_temp;
    $5 = &dim4_temp;
    $6 = &dim5_temp;
}

%typemap(in,numinputs=0)
    (uint32_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np)
    (uint32_t *data_temp, int dim1_temp, int dim2_temp,  int dim3_temp, int dim4_temp, int dim5_temp )
{
    $1 = &data_temp;
    $2 = &dim1_temp;
    $3 = &dim2_temp;
    $4 = &dim3_temp;
    $5 = &dim4_temp;
    $6 = &dim5_temp;
}

%typemap(argout)

    (uint8_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np)
{
    npy_intp dims[5] = { *$2, *$3, *$4, *$5, *$6 };
    PyObject* obj = PyArray_SimpleNewFromData(5, dims, NPY_UBYTE, (void*)(*$1));
  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result,obj);
}

%typemap(argout)

    (uint32_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np)
{
    npy_intp dims[5] = { *$2, *$3, *$4, *$5, *$6 };
    PyObject* obj = PyArray_SimpleNewFromData(5, dims, NPY_UINT32, (void*)(*$1));
  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result,obj);
}



namespace lm
{

namespace me
{
    class MESolver
    {
        MESolver();
        virtual ~MESolver();

        virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources)=0;
        virtual bool needsReactionModel()=0;
        virtual bool needsDiffusionModel()=0;
        //virtual void generateTrajectory()=0;
    };
}

namespace builder
{

struct point {
    point(si_dist_t x=0.0, si_dist_t y=0.0, si_dist_t z=0.0):x(x),y(y),z(z){}
    si_dist_t x;
    si_dist_t y;
    si_dist_t z;

    si_dist_t distanceSquared(const point & p2)
    {
        si_dist_t dx = p2.x - x;
        si_dist_t dy = p2.y - y;
        si_dist_t dz = p2.z - z;
        return (dx*dx + dy*dy + dz*dz);
    }

    si_dist_t distance(const point & p2) {return sqrt(distanceSquared(p2));}
};

struct bounding_box {
    bounding_box(si_dist_t x1=0.0, si_dist_t y1=0.0, si_dist_t z1=0.0, si_dist_t x2=0.0, si_dist_t y2=0.0, si_dist_t z2=0.0):min(x1,y1,z1),max(x2,y2,z2){}
    bounding_box(point min, point max):min(min),max(max){}
    point min, max;

    bounding_box joinWith(bounding_box j)
    {
        return bounding_box(::min(j.min.x,min.x),::min(j.min.y,min.y),::min(j.min.z,min.z),::max(j.max.x,max.x),::max(j.max.y,max.y),::max(j.max.z,max.z));
    }
};

struct vector {
    vector(si_dist_t x=0.0, si_dist_t y=0.0, si_dist_t z=0.0):x(x),y(y),z(z){}
    si_dist_t x;
    si_dist_t y;
    si_dist_t z;
};

class Shape
{
public:
    virtual bool boundingBoxesIntersect(Shape * query);
    virtual bool intersects(Shape * query) = 0;
    virtual bool contains(Shape * query) = 0;
    virtual bounding_box getBoundingBox();
    virtual site_t getType();
    virtual double getVolume();
};

class Sphere : public Shape
{
public:
    Sphere(point center, si_dist_t radius, site_t type);
    virtual ~Sphere();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual void setCenter(point center);
    virtual point getCenter();
    virtual si_dist_t getRadius();
};

class Hemisphere : public Shape
{
public:
    Hemisphere(point center, si_dist_t radius, vector orientation, site_t type);
    virtual ~Hemisphere();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual point getCenter();
    virtual si_dist_t getRadius();
    virtual vector getOrientation();
};

class Capsule : public Shape
{
public:
    Capsule(point p1, point p2, si_dist_t radius, site_t type);
    virtual ~Capsule();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual point getP1();
    virtual point getP2();
    virtual si_dist_t getRadius();
};
    
class Cylinder : public Shape
{
public:
    Cylinder(point p1, point p2, si_dist_t radius, site_t type);
    virtual ~Cylinder();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual point getP1();
    virtual point getP2();
    virtual si_dist_t getRadius();
};

class Cone : public Shape
{
public:
    Cone(point center, si_dist_t radius, si_dist_t height, site_t type, vector normal = vector(0.0, 0.0, 1.0));
    virtual ~Cone();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual point getCenter();
    virtual si_dist_t getRadius();
    virtual si_dist_t getHeight();
};


class CapsuleShell : public Shape
{
public:
    CapsuleShell(point p1, point p2, si_dist_t innerRadius, si_dist_t outerRadius, site_t type);
    virtual ~CapsuleShell();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual point getP1();
    virtual point getP2();
    virtual si_dist_t getInnerRadius();
    virtual si_dist_t getOuterRadius();
};

class Cuboid : public Shape
{
public:
    Cuboid(point p1, point p2, site_t type);
    Cuboid(point center, si_dist_t width, si_dist_t height, si_dist_t depth, site_t type, vector at = vector(1.0,0.0,0.0), vector up = vector(0.0,1.0,0.0));
    virtual ~Cuboid();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual point getP1();
    virtual point getP2();
};

class Union : public Shape
{
public:
    Union(Shape *s1, Shape *s2, site_t type);
    virtual ~Union();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
};

class UnionSet : public Shape
{
public:
    UnionSet(site_t type);
    virtual ~UnionSet();
    void addShape(Shape *s);
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
};

class Difference : public Shape
{
public:
    Difference(Shape *s1, Shape *s2, site_t type, bool symmetric = false);
    virtual ~Difference();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
};
    
class Intersection : public Shape
{
public:
    Intersection(Shape *s1, Shape *s2, site_t type);
    virtual ~Intersection();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
};
    
class Torus : public Shape
{
public:
    Torus(point center, si_dist_t radius1, si_dist_t radius2, site_t type, vector orientation = vector(0.0,0.0,1.0), vector up = vector(0.0,1.0,0.0));
    virtual ~Torus();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual void setCenter(point center);
    virtual point getCenter();
    virtual si_dist_t getRadius1();
    virtual si_dist_t getRadius2();
};
    
class Ellipse : public Shape
{
public:
    Ellipse(point center, si_dist_t radius1, si_dist_t radius2, si_dist_t radius3, site_t type, vector orientation1  = vector(0.0,0.0,1.0), vector orientation2 = vector(1.0, 0.0, 0.0));
    virtual ~Ellipse();
    virtual bool intersects(Shape * query);
    virtual bool contains(Shape * query);
    virtual void setCenter(point center);
    virtual point getCenter();
    virtual si_dist_t getRadius1();
    virtual si_dist_t getRadius2();
    virtual si_dist_t getRadius3();
};
    
class LatticeBuilder
{

public:
    LatticeBuilder(si_dist_t xLen, si_dist_t yLen, si_dist_t zLen, si_dist_t collisionGridSpacing, uint32_t seedTop, uint32_t seedBottom);
    virtual ~LatticeBuilder();
    virtual void addRegion(Shape * shape);
    //virtual bool placeObject(Shape * shape);
    virtual bool placeSphere(point position, si_dist_t radius, site_t type);
    virtual void removeSphere(point position, si_dist_t radius, site_t type);
    virtual uint placeRandomSphere(si_dist_t radius, site_t type, site_t region);
    virtual void placeRandomSpheres(uint count, si_dist_t radius, site_t type, site_t region);
    virtual void fillWithRandomSpheres(double volumeFraction, si_dist_t radius, site_t type, site_t region);
    virtual void getSpatialModel(lm::io::SpatialModel * model);
    virtual void addParticles(particle_t particleType, site_t siteType, uint count);
    virtual void discretizeTo(lm::rdme::Lattice * lattice, site_t obstacleSiteType, double fractionObstacleSitesOccupied);
};


}

namespace io {

class DiffusionModel
{
public:
    DiffusionModel();
    uint32_t number_species();
    uint32_t number_site_types();
    double lattice_spacing();
    uint32_t lattice_x_size();
    uint32_t lattice_y_size();
    uint32_t lattice_z_size();
    uint32_t particles_per_site();
    uint32_t bytes_per_particle();
    
    void set_lattice_spacing(double);
    void set_lattice_x_size(uint32_t);
    void set_lattice_y_size(uint32_t);
    void set_lattice_z_size(uint32_t);
    void set_particles_per_site(uint32_t);
    void set_number_reactions(uint32_t);
    void set_number_species(uint32_t);
    void set_number_site_types(uint32_t);
    void set_bytes_per_particle(uint32_t);

    double diffusion_matrix(uint32_t);
    void add_diffusion_matrix(double);
    void set_diffusion_matrix(uint32_t, double);

    uint32_t reaction_location_matrix(uint32_t);
    void add_reaction_location_matrix(uint32_t);
    void set_reaction_location_matrix(uint32_t, uint32_t);
};

class ReactionModel_Reaction
{
public:
    ReactionModel_Reaction();

    uint32_t type();
    void set_type(uint32_t);

    uint32_t rate_constant_size();
    void add_rate_constant(double);
    void set_rate_constant(uint32_t, double);
};
    
class ReactionModel
{
public:
    ReactionModel();
    uint32_t number_species();
    void set_number_species(uint32_t);

    uint32_t number_reactions();
    void set_number_reactions(uint32_t);

    uint32_t initial_species_count_size();
    void add_initial_species_count(uint32_t);
    void set_initial_species_count(uint32_t, uint32_t);

    uint32_t reaction_size();
    void add_reaction();
    ReactionModel_Reaction reaction(uint32_t);
    ReactionModel_Reaction* mutable_reaction(uint32_t);

    void add_stoichiometric_matrix(int);
    void add_dependency_matrix(uint32_t);
};
    
    

class SpatialModel
{
public:
    SpatialModel();
};

namespace hdf5 {

class SimulationFile
{
public:
    static bool isValidFile(const char * filename) throw(IOException,HDF5Exception);
    static void create(const char * filename) throw(IOException,HDF5Exception);
    
public:
    SimulationFile(const char* filename);
    virtual ~SimulationFile();
    virtual void close();
    virtual void getDiffusionModel(DiffusionModel * model);
    virtual void setDiffusionModel(lm::io::DiffusionModel * diffusionModel) throw(InvalidArgException,HDF5Exception,Exception);
    virtual void getReactionModel(ReactionModel * model);
    virtual void setReactionModel(lm::io::ReactionModel * reactionModel) throw(InvalidArgException,HDF5Exception,Exception);
    virtual void getDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::Lattice * lattice) throw(InvalidArgException,HDF5Exception,Exception);
    virtual void setDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::ByteLattice * lattice) throw(InvalidArgException,HDF5Exception,Exception);
    virtual void setDiffusionModelLattice(lm::io::DiffusionModel * m, lm::rdme::IntLattice * lattice) throw(InvalidArgException,HDF5Exception,Exception);
    virtual void getSpatialModel(SpatialModel * model);
    virtual void setSpatialModel(SpatialModel * model);
    virtual std::vector<double> getLatticeTimes(unsigned int replicate) throw(HDF5Exception,InvalidArgException);
    virtual void getLattice(unsigned int replicate, unsigned int latticeIndex, lm::rdme::Lattice * lattice) throw(HDF5Exception,InvalidArgException);
    virtual std::string getParameter(std::string key);
    virtual void setParameter(std::string k, std::string v);
    
    virtual std::vector<double> getSpeciesCountTimes(unsigned int replicate) throw(HDF5Exception,InvalidArgException);
    virtual std::map<double, std::vector<int> > getSpeciesCounts(unsigned int replicate) throw(HDF5Exception,InvalidArgException);
    std::map<uint, std::string> SimulationFile::getSpeciesNames() throw(HDF5Exception);
};

}
}

namespace cme
{

%include "numpy.i"

%init %{
    import_array();
%}

// Type maps for creating a numpy array memory view of species counts
%typemap(in,numinputs=0)
  (uint **counts, int *number)
  (uint *data_temp, int dim_temp)
{
  $1 = &data_temp;
  $2 = &dim_temp;
}
%typemap(argout)
  (uint **counts, int *number)
{
  npy_intp dims[1] = { *$2 };
  PyObject* obj = PyArray_SimpleNewFromData(1, dims, NPY_UINT, (void*)( *$1 ));
  PyArrayObject* array = (PyArrayObject*) obj;

  if(!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result, obj);
}

// Type maps for creating a numpy array view of the rates for a specific reaction
%typemap(in,numinputs=0)
  (double **rates, int *rateConstantCount)
  (double *data_temp, int dim_temp)
{
  $1 = &data_temp;
  $2 = &dim_temp;
}
%typemap(argout)
  (double **rates, int *rateConstantCount)
{
  npy_intp dims[1] = { *$2 };
  PyObject* obj = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void*)( *$1 ));
  PyArrayObject* array = (PyArrayObject*)obj;
  
  if(!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result, obj);
}



class CMESolver : public lm::me::MESolver
{
public:
    CMESolver(RandomGenerator::Distributions neededDists);
    virtual ~CMESolver();
    virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
    virtual void getSpeciesCountView(uint **counts, int *number);
    virtual void getReactionRateConstantsView(int reactionNumber, double **rates, int *rateConstantCount);
protected:
    virtual int hookSimulation(double time);
    virtual int onBeginTrajectory();
    virtual int onEndTrajectory();
};

%feature("director") GillespieDSolver;
class GillespieDSolver : public CMESolver
{
public:
    GillespieDSolver();
    virtual ~GillespieDSolver();
    virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
    virtual bool needsReactionModel();
    virtual bool needsDiffusionModel();
    //virtual void generateTrajectory();
    virtual void getSpeciesCountView(uint **counts, int *number);
    virtual void getReactionRateConstantsView(int reactionNumber, double **rates, int *rateConstantCount);
protected:
    virtual int hookSimulation(double time);
    virtual int onBeginTrajectory();
    virtual int onEndTrajectory();
};

%feature("director") NextReactionSolver;
class NextReactionSolver : public CMESolver
{
public:
    NextReactionSolver();
    NextReactionSolver(RandomGenerator::Distributions neededDists);
    virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
    virtual ~NextReactionSolver();
    virtual bool needsReactionModel();
    virtual bool needsDiffusionModel();
    //virtual void generateTrajectory();
    virtual void getSpeciesCountView(uint **counts, int *number);
    virtual void getReactionRateConstantsView(int reactionNumber, double **rates, int *rateConstantCount);

protected:
    virtual int hookSimulation(double time);
};

}

namespace rdme
{

unsigned int getCompiledLatticeMaxOccupancy();

class Lattice
{
public:
    virtual site_size_t getMaxOccupancy() const;
    virtual lattice_coord_t getSize() const;
    virtual lattice_size_t getXSize() const;
    virtual lattice_size_t getYSize() const;
    virtual lattice_size_t getZSize() const;
    virtual lattice_size_t getNumberSites() const;
    virtual si_dist_t getSpacing() const;
    
    virtual site_t getSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z)=0;
    virtual void setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t siteType)=0;
    virtual site_size_t getOccupancy(lattice_size_t x, lattice_size_t y, lattice_size_t z) const throw(InvalidSiteException)=0;
    virtual particle_t getParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_size_t particleIndex) const throw(InvalidSiteException,InvalidParticleException)=0;
    virtual std::vector<particle_loc_t> findParticles(particle_t minParticleType, particle_t maxParticleType)=0;
    
    virtual void addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle) throw(InvalidSiteException,InvalidParticleException)=0;
    virtual void removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z) throw(InvalidSiteException)=0;
    virtual std::map<particle_t,uint> getParticleCounts()=0;

};

class ByteLattice : public Lattice
{
public:
    ByteLattice(lattice_coord_t size, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    ByteLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    virtual ~ByteLattice() throw(std::bad_alloc);
    virtual site_t getSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z);
    virtual void setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t siteType);
    virtual site_size_t getOccupancy(lattice_size_t x, lattice_size_t y, lattice_size_t z) const throw(InvalidSiteException);
    virtual particle_t getParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_size_t particleIndex) const throw(InvalidSiteException,InvalidParticleException);
    virtual std::vector<particle_loc_t> findParticles(particle_t minParticleType, particle_t maxParticleType);
    
    virtual void addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle) throw(InvalidSiteException,InvalidParticleException);
	virtual void removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z) throw(InvalidSiteException);
    virtual std::map<particle_t,uint> getParticleCounts();
    virtual void getSiteLatticeView(uint8_t **siteLattice, int *Nz, int *Ny, int *Nx);
    virtual void getParticleLatticeView(uint8_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np);
};

class IntLattice : public Lattice
{
public:
    IntLattice(lattice_coord_t size, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    IntLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    virtual ~IntLattice() throw(std::bad_alloc);
    virtual site_t getSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z);
    virtual void setSiteType(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_t siteType);
    virtual site_size_t getOccupancy(lattice_size_t x, lattice_size_t y, lattice_size_t z) const throw(InvalidSiteException);
    virtual particle_t getParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, site_size_t particleIndex) const throw(InvalidSiteException,InvalidParticleException);
    virtual std::vector<particle_loc_t> findParticles(particle_t minParticleType, particle_t maxParticleType);
    
    virtual void addParticle(lattice_size_t x, lattice_size_t y, lattice_size_t z, particle_t particle) throw(InvalidSiteException,InvalidParticleException);
	virtual void removeParticles(lattice_size_t x,lattice_size_t y,lattice_size_t z) throw(InvalidSiteException);
    virtual std::map<particle_t,uint> getParticleCounts();
    virtual void getSiteLatticeView(uint8_t **siteLattice, int *Nz, int *Ny, int *Nx);
    virtual void getParticleLatticeView(uint32_t **particleLattice, int *Nw, int *Nz, int *Ny, int *Nx, int *Np);
};

#ifdef OPT_CUDA 
class CudaByteLattice: public ByteLattice
{
public:
    CudaByteLattice(lattice_coord_t size, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    CudaByteLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    virtual ~CudaByteLattice() throw(std::bad_alloc);
    virtual std::map<particle_t,uint> getParticleCounts();

    virtual void * getGPUMemorySrc();
    virtual void * getGPUMemoryDest();
    virtual void swapSrcDest();
    virtual void * getGPUMemorySiteTypes();


};
class CudaIntLattice: public IntLattice
{
public:
    CudaIntLattice(lattice_coord_t size, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    CudaIntLattice(lattice_size_t xSize, lattice_size_t ySize, lattice_size_t zSize, si_dist_t spacing, uint particlesPerSite) throw(std::bad_alloc,InvalidArgException,Exception);
    virtual ~CudaIntLattice() throw(std::bad_alloc);

    virtual void * getGPUMemorySrc();
    virtual void * getGPUMemoryDest();
    virtual void swapSrcDest();
    virtual void * getGPUMemorySiteTypes();


};


%feature("notabstract")  MpdRdmeSolver;
%feature("director") MpdRdmeSolver;
class MpdRdmeSolver : public lm::me::MESolver
{
public:
    MpdRdmeSolver();
    ~MpdRdmeSolver();
        virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
        virtual bool needsReactionModel();
        virtual bool needsDiffusionModel();
        //virtual void generateTrajectory();
protected:
    virtual int hookSimulation(double time, CudaByteLattice *lattice);
    virtual int onWriteLattice(double time, CudaByteLattice *lattice);
    virtual int onBeginTrajectory(CudaByteLattice *lattice);
    virtual int onEndTrajectory(CudaByteLattice *lattice);
    virtual void writeLatticeSites(double time, CudaByteLattice *lattice);
};


%feature("director") MGPUMpdRdmeSolver;
%feature("notabstract")  MGPUMpdRdmeSolver;
class MGPUMpdRdmeSolver : public lm::me::MESolver
{
public:
    MGPUMpdRdmeSolver();
    ~MGPUMpdRdmeSolver();
        virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
        virtual bool needsReactionModel();
        virtual bool needsDiffusionModel();
        virtual void setReactionRate(unsigned int rxid, float rate);
protected:
    virtual int hookSimulation(double time, ByteLattice *lattice);
};

%feature("director") IntMpdRdmeSolver;
%feature("notabstract")  IntMpdRdmeSolver;
class IntMpdRdmeSolver : public lm::me::MESolver
{
public:
     IntMpdRdmeSolver();
    ~IntMpdRdmeSolver();
        virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
        virtual bool needsReactionModel();
        virtual bool needsDiffusionModel();
protected:
    virtual int hookSimulation(double time, CudaIntLattice *lattice);
};

/*
%feature("director") MpdHybridSolver;
class IntMpdHybridSolver : public IntMpdRdmeSolver
{
public:
    MpdHybridSolver();
    ~MpdHybridSolver();
        virtual void initialize(unsigned int replicate, map<string,string> * parameters, ResourceAllocator::ComputeResources * resources);
        virtual bool needsReactionModel();
        virtual bool needsDiffusionModel();
        virtual void generateTrajectory();
protected:
    virtual int hookSimulation(double time, CudaIntLattice *lattice);
};
*/
#endif

}


}


void runSimulation(char *simulationFilename, int replicate, char *solverClass,
               std::vector<int> cudaDevices, time_t checkpointInterval);

void runSolver(char *simulationFilename, int replicate, lm::me::MESolver *solver,
               std::vector<int> cudaDevices, time_t checkpointInterval);


