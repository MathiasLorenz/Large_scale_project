README for jacobiSolver
-------------------------------------------------------------------------------
The Jacobi solver have been implemented into the jacobiSolver executable.

Default call to the function will be as bellow, please be adviced mpiexec or 
mpirun must be used for all mpi versions.
	> ./jacobiSolver METHOD NX NY NZ

where:
	METHOD :
		omp2d	: Is the Jacobi method using OpenMP for a 2D problem.
		omp3d	: Is the Jacobi method using OpenMP for a 3D problem.
		mpi3d_1	: Is the Jacobi method using MPI with a single split along Z.
		mpi3d_2	: Is the Jacobi method using MPI with multiple split along Z.
		cuda_1	: First implementation using cuda for a single GPU.
		mixed_1	: First implementation using both MPI and CUDA. Very naive.
		mixed_2	: Second implementation using both MPI and CUDA. Copies only 
				  needed data.

	NX NY NZ : Integer numbers.
		The number of points in the problem, fx. 7 by 7. Nz is the size in the 
		3rd dimension. NY and NZ can be omitted for cubic domains.
			B B B B B B B
			B X X X X X B
			B X X X X X B
			B X X X X X B
			B X X X X X B
			B X X X X X B
			B B B B B B B
		Where the boundary is set during the function and the elemnets in X is
		computed.

-------------------------------------------------------------------------------
The solver will read for several environmental variables and change behavior
based on these. The environmental variables supported are listed bellow.

OpenMP: All default OpenMP controls such as setting number of threads
		are present and jomp have been implemented such that the schedule
		can be set during Runtime using OMP_SCHEDULE.

Extra environmental variables: ( > ENV_NAME=value ./jacobiSolver ...)
	PROBLEM_NAME: [default: sin]
		Defines the problem that should to be solved.

		sin : 	The well defined example where the Boundary condition and
				source is defined as:
				(2D version):
					f(x,y)   = 2*pi^2*sin(pi*x)*sin(pi*y)
					u(x,y)   = 0, for (x,y) on boundary.
				(3D version):
					f(x,y,z) = 3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)
					u(x,y,z) = 0, for (x,y,z) on boundary.

	OUTPUT_INFO: [default: timing]
		Defines the output type.

		matrix_slize: 
			Defines that the result matrix should be printed as the
			output. This will print the slize at the center of z for,
			3 dimensional problems.
		matrix_full:
			Defines that the result matrix should be printed as the
			output for the full 3 dimensional problem. The format of 
			the output will be:
				Nx Ny Nz
				U(0,0,0) U(0,0,1) ....
		error:
			Prints an output of the dimension sizes in terms of number of
			gridpoints in each and the maximal absolute difference between 
			the analytical solution and the solution calculated in the program.	
		timing:
			Defines that the output should be the timing info. Giving
			the memory footprint (kB), Mflops and Walltime for the main loop.

	TOLERANCE:	[default: 1e-6]
		Sets the tolerance for the program to terminate when the maximal 
		absolute error falls under this number.

	USE_TOLERANCE: [default: on] 
		Defines if the tolerance should be used or not.

		on :	Use the tolerance as stop criterion.
		off:	Do not use the tolerance. Force to do MAXITER iterations.

	MAX_ITER: [default: 1e5]
		The maximal number of iterations done by the solver.


% = = EOF = = %
