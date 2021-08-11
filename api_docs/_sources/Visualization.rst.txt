Visualization using VMD
~~~~~~~~~~~~~~~~~~~~~~~

`Visual Molecular Dynamics <http://www.ks.uiuc.edu/Research/vmd/>`_ can be used
to visualize RDME simulation output. If LM was built and installed with VMD
support, the required molfile plugin to load LM trajectories will be installed
automatically. VMD has excellent `documentation
<http://www.ks.uiuc.edu/Research/vmd/current/docs.html>`_,  so only LM specific
features of the plugin will be discussed here.

After loading a trajectory, both the site type lattice and the RDME particles
will be displayed as points. Since this will be cluttered and hard to
interpret, the ``atomselect`` syntax used in the Representations dialog can be
used to select which components of the simulation to display. The lattice
geometry can be visualized using the selector ``segname SITE``, where as the
chemical particles can be visualized using ``segname RDME``. For example, for a
cell :math:`1\;\mu\mathrm{m}` thick with its  longest dimension aligned with
the :math:`z` axis, the ``atomselect`` syntax to show only the bottom part of
the cell is

.. code-block:: none

    segname SITE and name membrane and y < 5000

where the name of the membrane site type is ``membrane``, and we display only
the bottom :math:`500\;\mathrm{nm}` of the  cell. Note that units in VMD are in
Ångstroms, not nanometers. 

Consider the system

.. math::
    \mathrm{A} + \mathrm{B} &\to \mathrm{AB}\\
    \mathrm{A} + \mathrm{C} &\to \mathrm{AC}

To view only the dimers in the simulation volume, one could write

.. code-block:: none

    segname RDME and (name AB or name AC)

or instead use a regular expression

.. code-block:: none

    segname RDME and name "A[BC]"

By combining multiple representations with different ``atomselect`` statements
and drawing methods, it is possible to quickly create interactively
visualizations in VMD which allow the dynamics of interest to be investigated.

