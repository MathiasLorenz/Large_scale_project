// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#include "tests.h"
#include "tests_cuda.h"

double FLOP=0.0;

// ============================================================================
// MAIN FUNCTION
int main(int argc, char const *argv[])
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
	if (Nx <= 2 || Ny <= 2){
		fprintf(stderr,
		"\nInvalid dimension inputs. Minimal value is 3.\n"
		"	Nx: %i, Ny: %i, Nz: %i.\n\n",Nx,Ny,Nz);
		return EXIT_FAILURE;
	}
	
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
	double t = omp_get_wtime();

	if (strcmp(T,"omp2d") == 0)
		test_jacobi_2D(Nx, Ny);

	else if (strcmp(T,"omp3d") == 0)
		test_jacobi_3D(Nx, Ny, Nz);

	else if (strcmp(T,"cuda") == 0)
		test_cuda(Nx, Ny, Nz);

	else {
		fprintf(stderr, 
			"\nInvalid test name\n"
			"   Accepts: omp2d, omp3d, cuda\n\n");
		return EXIT_FAILURE;
	}

	// ------------------------------------------------------------------------
	// Handling the printing of statistics and data.
	
	double timespent = omp_get_wtime() - t;
	if (strcmp("timing",getenv("OUTPUT_INFO")) == 0){
		printf("Mflops: %10.4f ", FLOP/timespent*1e-6 );
		printf("Walltime: %10.4f\n", timespent);
	}

	return EXIT_SUCCESS;
}