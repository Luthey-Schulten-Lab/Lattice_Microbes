# Lattice Microbes

**Lattice Microbes** is a program used to simulate stochastic reaction-diffusions dynamics.  It consists of implementations of both the well-stirred chemical master equation (CME) and spatailly-resolved reaction-diffusion master equation (RDME).

## Modules

**jLM** is a module used for constructing cell geomtries and defining reaction and diffusion rules for a RDME simualtion with a python/jupyter notebook interface.  Experimental data, such as cell organization from cryo-electron tomograms, can be incorporated into simulation cell geometry also using region builder functions.

**pyLM** can define well-stirred reaction systems without diffusion to be simulated using the CME with a python interface.

## Solvers

*GillespieDSolver*

*MpdRdmeSolver*

*IntMpdRdmeSolver*

*MGPUMpdRdmeSolver*

## Resources

- [Installation guide](http://faculty.scs.illinois.edu/schulten/software_manuals/LM_2.3_INSTALL_annotations.pdf): Instructions for installing **Lattice Microbes**.
- [User Guide](http://faculty.scs.illinois.edu/schulten/software_manuals/UsersGuide.pdf): Tutorial, full Python API description, and usage information.
- [Citing Lattice Microbes](http://faculty.scs.illinois.edu/schulten/Software2.0.html): How to cite the code.
