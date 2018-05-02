// ============================================================================
// INCLUDES
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <omp.h>

#include "init_data.h"
#include "jacobi_util.h"
#include "tests_cuda.h"
#include "poisson_cuda.h"
#include "matrix_routines.h"
#include "cuda_routines.h"

extern double MEMORY;
extern double TIME_SPENT;

// ============================================================================
// CUDA TEST
void test_cuda_1(Information *information)
{
	// Read the information structure
	int rank = information->rank;
	int Nx = information->global_Nx;
	int Ny = information->global_Ny;
	int Nz = information->global_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];

	// Allocation
	double *U, *F, *Unew, *A = NULL;
	U = dmalloc_3d_l(loc_Nx, loc_Ny, loc_Nz);
	F = dmalloc_3d_l(loc_Nx, loc_Ny, loc_Nz);
	Unew = dmalloc_3d_l(loc_Nx, loc_Ny, loc_Nz);
	if (!U || !F || !Unew) {
		// Consider fail cases if one thread ies and the others do not
		fprintf(stderr, "Error in malloc, pointer is NULL.\n");
		return;
	}

	// Array for true solution if requested
	if (strcmp("error", getenv("OUTPUT_INFO")) == 0)
	{
		A = dmalloc_3d_l(loc_Nx, loc_Ny, loc_Nz);
		if (!A) { fprintf(stderr, "Error in malloc, pointer is NULL.\n"); return; }
		generate_true_solution(A, information);
	}

	// Initialise the boundary values
	if (strcmp("sin", getenv("PROBLEM_NAME")) == 0)
	init_sin_3D(U, F, Unew, Nx, Ny, Nz);
	// else if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
	//    init_rad_2D(U, F, Unew, Nx, Ny);
	else {
		fprintf(stderr, "Problem type is not supported.\n");
		return;
	}

	// Handle the environmental variables
	int maxiter = atoi(getenv("MAX_ITER"));
	double tol = atof(getenv("TOLERANCE"));

	// Main computation and time
	double t = omp_get_wtime();
	jacobi_cuda_1(Nx,Ny,Nz,maxiter,tol,U,F,Unew);

	// Save global variables
	TIME_SPENT = omp_get_wtime() - t;
	MEMORY = 3.0 * Nx * Ny * Nz * 8.0 / 1024.0;

	// Print the needed information
	if (strcmp("matrix", getenv("OUTPUT_INFO")) == 0)
		array_print_3d_slice(U, Nx, Ny, Nz, Nz / 2, "%10g ");
	else if (strcmp("full_matrix", getenv("OUTPUT_INFO")) == 0)
		array_print_3d(U, Nx, Ny, Nz, "%10g ");
	else if (strcmp("full_matrix_mpi_z_slice", getenv("OUTPUT_INFO")) == 0)
		print_jacobi3d_z_sliced(U, information, "%10g ");
	else if (strcmp("error", getenv("OUTPUT_INFO")) == 0)
	{
		// Compute absolute error
		double abs_err = 0.0;
		check_true_solution(A, U, &abs_err, information);
		//printf("I'm rank %d, local error: %f\n", rank, loc_abs_err);
		if (rank == 0)
			printf("Grid: %d %d %d, error: %.10f\n", Nx, Ny, Nz, abs_err);
	}

	// Free the arrays created for the computation
	free(U); free(F); free(A); free(Unew);	
}

