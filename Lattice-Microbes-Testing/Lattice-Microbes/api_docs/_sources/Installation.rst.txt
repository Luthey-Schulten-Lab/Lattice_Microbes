Installation
============

Dependencies
------------

 * `CMake <https://cmake.org>`_ 3.10 or newer.

 * `Google Protocol Buffers <https://developers.google.com/protocol-buffers/>`_

 * `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_

 * `Python <https://www.python.org>`_ 3.5 or newer, and the standard scientific
   Python stack. (*Optional Python support*)

 * `Swig <http://www.swig.org>`_  (*Optional Python support*)

 * `CUDA <https://developer.nvidia.com/cuda-zone>`_, tested on 9.2, please
   report if other versions do not work (*Optional MPD-RDME support*)

 * `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_, current release tested,
   please report trouble with other versions (*Optional VMD plugin*)


Preparation
-----------

Although not required, a Python virtual environment is recommended for
interacting with Lattice Microbes using the Python bindings. We recommend using
the Anaconda python distribution since all dependencies (except CUDA) can be
installed directly using the package manager ``conda``.  Standard Python
virtual environments (i.e. ``python -m venv``) are an option as well, however
it will require that a modern c++ compiler.  Finally, a system-wide install is
possible, however doing so is discouraged.

Conda (recommended)
~~~~~~~~~~~~~~~~~~~

Install CUDA following the `documentation
<https://docs.nvidia.com/cuda/cuda-quick-start-guide/index.html>`_.

Install `Miniconda <https://conda.io/miniconda.html>`_ (or alternatively
`Anaconda <https://www.anaconda.com/distribution/>`_) and dependencies

.. code-block:: bash

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p $HOME/miniconda
    rm /tmp/miniconda.sh
    $HOME/miniconda/bin/conda create -n lattice_microbes
    source $HOME/miniconda/bin/activate lattice_microbes
    conda install -y protobuf hdf5 numpy cython cmake swig jupyter matplotlib \
        ipywidgets h5py tqdm pillow jinja2 scipy
    conda install -y -c anaconda git gcc_linux-64 gxx_linux-64 binutils_linux-64

Activation of the virtual environment will change your environment variables
(e.g. PATH, CC, CXX) so that programs will be searched for in your virtual
environment first. To activate the virtual environment created above,
``lattice_microbes``,  run

.. code-block:: bash

    source $HOME/miniconda/bin/activate lattice_microbes

When your are done using the virtual environment either run

.. code-block:: bash

    source deactivate

or simply close the terminal window.


Python Virtualenv or Global Install (advanced)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install CUDA, CMake, Google Protocol buffers, HDF5, and Swig. Your Linux
distribution package manager (or MacOS homebrew) should be sufficient, however
you will need a modern C++ compiler which supports C++11.

You will need the following Python modules:

 * h5py

 * NumPy

 * Cython

 * Jupyter

 * Matplotlib

 * IPywidgets

 * tqdm

 * Pillow

 * Jinja2

 * SciPy

They can be installed using your Linux distribution's package manager, or
through ``pip``.


Build
-----

Create a working directory for the build. This is a temporary directory can be
placed anywhere in the filesystem. Then generate the build files

.. code-block:: bash

    mkdir build
    cd build
    cmake ../Lattice-Microbes/src # or wherever the source was extracted to

Then,

.. code-block:: bash

    make

Install
-------

Simply run

.. code-block:: bash

    make install

This will install the python modules and executables to your virtual
environment prefix.  The `make install` should have also picked up the VMD
executable from your path and installed the molfile plugin properly.  Finally,
start Jupyter and get started with LM.

Advanced Build Configuration
============================

The CMake build system provides a GUI for specifying build options.  To
configure the build, use ``ccmake`` for a curses-based interface, or
``cmake-gui`` for a GUI. If your prefer the command line, run 
``cmake -D<OPTION_NAME>=<VALUE> .`` in a configured build directory to set
build options. The standard flags are:

    **OPT_CUDA**
        Enable CUDA support

    **OPT_MPI**
        Enable MPI support. *Not currently working*

    **OPT_PYTHON**
        Enable Python bindings.

    **OPT_VMD**
        Build VMD plugin

There are also advanced options whose visibility can be toggled in the GUI. The
most common are

    **MPD_GLOBAL_S_MATRIX**
        Declare the stoichiometric matrix as ``__global__`` instead of
        ``__constant__``.

    **MPD_GLOBAL_T_MATRIX**
        Declare the transition matrix as ``__global__`` instead of
        ``__constant__``.

    **MPD_LATTICE_MAX_OCCUPANCY**
        Specify the maximum number of particles per subvolume. Currently only 8
        or 16 particles are valid choices.

    **MPD_MAX_REACTION_TABLE_ENTRIES**
        Specify the size of the reaction table.

    **MPD_MAX_RL_MATRIX_ENTRIES**
        Specify the size of the reaction location matrix.

    **MPD_MAX_S_MATRIX_ENTRIES**
        Specify the size of the stoichiometric matrix if it has been declared in
        ``__constant__`` memory.

    **MPD_MAX_TRANSITION_TABLE_ENTRIES**
        Specify the size of the transition table if it has been declared in
        ``__constant__`` memory.

