Introduction
============

Studying cellular processes, which show inherently noisy non-deterministic
behavior, with single molecule resolution on timescales of biological relevance
such as the lifetime of a cell, requires considerable computational effort
[EARN2018]_.  Lattice Microbes [ROBE2009]_ [ROBE2013]_ is software developed to
sample realizations of the spatially homogeneous and heterogeneous stochastic
Master equations, with thousands of reactions among hundreds of molecular
species.  The software uses graphics processing units (GPUs) to exploit the
natural parallelism afforded by the Master equations to access timescales
orders-of-magnitude larger than other particle- and grid-based software for
sampling stochastic cellular processes.  While Lattice Microbes was originally
designed for simulating single *E. coli* cells on one GPU, the desire to
simulate cellular consortia and larger species like yeast, drove the
development of a new version of Lattice Microbes that utilizes multiple GPUs to
share the work [HALL2014]_. In addition to larger simulations, multiple GPUs
allow small simulations to be completed more quickly.

Capabilities
------------

Simulation modes
~~~~~~~~~~~~~~~~

Lattice Microbes can be used to simulate trajectories following either a
chemical master equation (CME),

.. math::
    \frac{\mathrm{d}P(\mathbf{x},t)}{\mathrm{d}t} 
    = \sum_{r}^{R} [-a_r({{\mathbf{x}}}) P({{\mathbf{x}}},t) 
    + a_r({{\mathbf{x}}}_\nu-\mathbf{S_r}) P({{\mathbf{x}}}-\mathbf{S_r},t)],

or a reaction-diffusion master equation (RDME),

.. math::
    \frac{\mathrm{d} P(\mathbf{x},t)}{\mathrm{d}t} 
    =& \sum_{\nu}^{V}\sum_{r}^{R} [-a_r({{\mathbf{x}}}_\nu) P({{\mathbf{x}}}_\nu,t) 
    + a_r({{\mathbf{x}}}_\nu-\mathbf{S_r}) P({{\mathbf{x}}}_\nu-\mathbf{S_r},t)]\\
    &+ \sum_{\nu}^{V}\sum_{\xi}^{\pm\hat{i},\hat{j},\hat{k}}\sum_{\alpha}^{N}
    [-d^{\alpha} x_{\nu}^{\alpha} P({{\mathbf{x}}},t) + 
    d^{\alpha} (x_{\nu+\xi}^{\alpha}+1) 
    P({{\mathbf{x}}}+1_{\nu+\xi}^{\alpha}-1_{\nu}^{\alpha},t)]

using a variety of methods.  Both CME and RDME simulations support a number of
different reaction types including zeroth, first, second, and second order self
reactions. Additionally, the Python interfaces to LM allow for interactive
construction of the simulation geometry in RDME simulations.

Through the use of custom solvers, different simulation methods can be combined
to simulate cellular processes spanning different length and concentration scales.

Execution modes
~~~~~~~~~~~~~~~

TODO

Command line

Python: pyLM

Python: jLM


