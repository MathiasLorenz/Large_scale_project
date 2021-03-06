# Large scale project
This large scale project is an implementation of the Jacobi method for solving 
steady-state partial differential equations. The implementation is based on 
OpenMP, openMPI and CUDA to test the performance gains of using parallel
technology for finding the steady-state solution to a partial differential
equation.
In our implementation we have tested the program with a simple heat equation
using a sinusodial function as the source term. 
The repository contain all source code for the implementation along with shell
and matlab files which can be used to test it.

## Folder structure
The repository contain two folders, 
- Poisson 
	- 	Contain all files needed for compiling the program including 
		a makefile and a readme describing the functionality of the 
		program in depth.
- Testing
	- 	Containg shell scripts and matlab code used for testing and visualizing
		the results of the program.
- Technical report describing the process of developing the code and the theory 
	the behind problem and method.
- Additional files used by us for streamlining the development process of the 
  project.

### Prerequisites
In order to compile and run the code the following libraries must be available: 
(the versions specified in paranthesis is the tested versions).
- GCC	( 6.4.0 built with OpenMP suport )
- MPI   ( openMPI 3.1.1 )
- CUDA  ( 9.2 )

### Installing
The poisson folder should be readily usable by the make utility. Please examine
the makefile before compiling to check for compliation flags and compilers. 

## Running the tests
For direct usage of the jacobiSolver in the commandline please see the 
readme located in the Poisson folder. This describes the functionality in 
depth.
The test scripts located in Testing/shell is designed to run on the LSF based
cluster at the Technical University of Denmark. Please modify these in order to
use them elsewhere. 
No additional packages are needed in order to run the matlab code used in the 
visualization. Matlab version R2017a was used for testing.

## Authors
* **Mathias Lorenz** -  [MathiasLorenz](https://github.com/MathiasLorenz)
* **Tim Felle Olsen** -  [TimFelle](https://github.com/TimFelle)

## Project supervisors
* **Bernd Dammann**
* **Hans Henrik Brandenborg Sørensen**
