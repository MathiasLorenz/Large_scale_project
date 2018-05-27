// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "tests.h"
#include "jacobi_util.h"
#include "tests_cuda.h"

double MFLOP=0.0;
double MEMORY=0.0;
double TIME_SPENT=0.0;

// ============================================================================
// MAIN FUNCTION
int main(int argc, char *argv[])
{
	// ------------------------------------------------------------------------
	// Handle input and help.
	if (argc == 1){
		printf(
		"\nThis function accepts the following arguments,\n"
		"	testName,	name of the test that should be performed\n"
		"	Nx			size of the first axis in the grid.\n"
		"	Ny			size of the second axis in the grid.\n"
		"	Nz			size of the third axis in the grid.\n\n"
		"See README.txt for additional information.\n\n");
		return EXIT_SUCCESS;
	} else if (argc < 3){
		fprintf(stderr,
			"\nNot enough input arguments\n"
			"Solution:\n"
			"  Run the function again with no inputs to see help.\n"
			"  See README.txt for full documentation.\n\n");
		return EXIT_FAILURE;
	}

	// ------------------------------------------------------------------------
	// Initalize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int size, rank;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	// ------------------------------------------------------------------------
	// Reading the inputs

	char const *T = argv[1];
	int Nx, Ny, Nz;

	if (argc < 4)
		Nx = Ny = Nz = atoi(argv[2]);
	else if (argc < 5){
		Nx = atoi(argv[2]);
		Ny = atoi(argv[3]);
		Nz = 3;
	} else {
		Nx = atoi(argv[2]);
		Ny = atoi(argv[3]);
		Nz = atoi(argv[4]);
	}
	if (Nx <= 2 || Ny <= 2 || Nz <= 2) {
		fprintf(stderr,
		"\nInvalid dimension inputs. Minimal value is 3.\n"
		"	Nx: %i, Ny: %i, Nz: %i.\n\n",Nx,Ny,Nz);
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	// ------------------------------------------------------------------------
	// Setup information structure
	
	Information information;
	write_information(&information,Nx,Ny,Nz,rank,size);
	
	// ------------------------------------------------------------------------
	// Handle Enviromental values

	char *problem_name, *output_info, *use_tol, *tol_env, *maxiter;

	if ( (problem_name = getenv("PROBLEM_NAME")) == NULL )
		putenv("PROBLEM_NAME=sin");
	if ( (output_info = getenv("OUTPUT_INFO")) == NULL )
		putenv("OUTPUT_INFO=timing");
	if ( (use_tol = getenv("USE_TOLERANCE")) == NULL )
		putenv("USE_TOLERANCE=on");
	if ( (maxiter = getenv("MAX_ITER")) == NULL)
		putenv( "MAX_ITER=10000" );
	if ( (tol_env = getenv("TOLERANCE")) == NULL)
		putenv( "TOLERANCE=1e-6" );
	
	// ------------------------------------------------------------------------
	// Make the call for the desired test

	if (strcmp(T, "omp2d") == 0)
		test_jacobi_2D(Nx, Ny);

	else if (strcmp(T, "omp3d") == 0)
		test_jacobi_3D(Nx, Ny, Nz);

	else if (strcmp(T, "mpi3d_1") == 0)
		test_jacobi_mpi3D_1(&information);

	else if (strcmp(T, "mpi3d_2") == 0)
		test_jacobi_mpi3D_2(&information);

	else if (strcmp(T, "cuda_1") == 0)
		test_cuda_1(&information);

	else if (strcmp(T, "mixed_1") == 0)
		test_jacobi_mixed_1(&information);
	else if (strcmp(T, "mixed_2") == 0)
		test_jacobi_mixed_2(&information);

	else {
		fprintf(stderr, 
			"\nInvalid test name: %s\n"
			"   Accepts: omp2d, omp3d, mpi3d_1, mpi3d_2, cuda_1, mixed_1\n\n",T);
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	// ------------------------------------------------------------------------
	// Handling the printing of statistics and data.
	
	if (strcmp("timing",getenv("OUTPUT_INFO")) == 0 && rank == 0){
		printf("Memory: %10.4f ", MEMORY);
		printf("Mflops: %10.4f ", MFLOP/TIME_SPENT );
		printf("Walltime: %10.4f\n", TIME_SPENT);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	free_information_arrays(&information);
    MPI_Finalize();
	return EXIT_SUCCESS;
}
