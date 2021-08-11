# Lattice Microbes

**Lattice Microbes** is a program used to simulate stochastic reaction-diffusions dynamics.  It consists of implementations of both the well-stirred chemical master equation (CME) and spatailly-resolved reaction-diffusion master equation (RDME).

## Modules

**jLM** is a python/jupyter notebook interface used for constructing cell geomtries and defining reaction and diffusion rules for a RDME simualtion.  Experimental data, such as cell organization from cryo-electron tomograms, can be incorporated into simulation cell geometry also using region builder functions.

**pyLM** is a python interface used to define well-stirred reaction systems without diffusion to be simulated using the CME.

## Solvers

**GillespieDSolver** is the Gillespie direct chemical master equation solver for well-stirred simulations

**MpdRdmeSolver** is the base version of the multi-particle diffusion reaction-diffusion master equation solver for spatial simulations

**IntMpdRdmeSolver** is the MpdRdmeSolver, but particle indexes are stored as 32-bit numbers rather than 8-bit.  This is used for systems with many types of molecules and intermediates.

**MGPUMpdRdmeSolver** is the multi-GPU MpdRdmeSolver which is used to divide the simulation space for one spatial simulation across multiple GPUs.

## Resources

- [Installation guide](http://faculty.scs.illinois.edu/schulten/software_manuals/LM_2.3_INSTALL_annotations.pdf): Instructions for installing **Lattice Microbes**.
- [User Guide](http://faculty.scs.illinois.edu/schulten/software_manuals/UsersGuide.pdf): Tutorial, full Python API description, and usage information.
- [Citing Lattice Microbes](http://faculty.scs.illinois.edu/schulten/Software2.0.html): How to cite the code.
