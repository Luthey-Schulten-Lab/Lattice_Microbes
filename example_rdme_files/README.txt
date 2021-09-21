Example executable rdme simulation files.

The tutorial contains a simple E. coli model with just a few reactions.

The MinCell tar file contains a simulation of just the localized reactions in the minimal cell.  Extract the simulation file by running:

tar -xzvf MinCell.tar.gz 

To run the example rdme simulation files after extracting from the tar, run the following line:

lm -g 0 -r 1 -sl lm::rdme::IntMpdRdmeSolver -f MinCell/MinCell_5s_testing.lm
or
lm -g 0 -r 1 -sl lm::rdme::MpdRdmeSolver -f LM_tutorial.lm

-g the index of the GPU to use
-r the index of the replicate to run
-sl the Reaction-Diffusion Master equation Solver to use
-f the simulation filename

LM_tutorial.lm can be simulated using only MpdRdmeSolver
MinCell_5s_testing.lm will only run correctly using IntMpdRdmeSolver
