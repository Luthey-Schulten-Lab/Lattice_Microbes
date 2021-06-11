pyLM
----

pyLM is a Problem Solving Environment (PSE) for biological simulations
[PETE2013]_.  Written in Python, it wraps and extends the highly optimized
multi-GPU `Lattice Microbes`_ stochastic simulation software
[ROBE2009]_, [ROBE2013]_.  The PSE is comprised of a base set of functionality
to set up, monitor and modify simulations, as well as a set of standard
post-processing routines that interface to other Python packages, including
NumPy, SciPy, H5py, iGraph to name a few. See [PETE2013]_ for additional
information as well as the user guide on the main `Lattice Microbes`_ website.
If you use pyLM in your simulations, please cite references [PETE2013]_ and
[ROBE2013]_.

.. _`Lattice Microbes`: http://www.scs.illinois.edu/schulten/lm

.. autosummary::
    :toctree: _autosummary

    pyLM
    pyLM.CME
    pyLM.ipyInterface
    pyLM.LMLogger
    pyLM.RDME
    pyLM.units

