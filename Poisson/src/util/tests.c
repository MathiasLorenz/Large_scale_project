// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#include "tests.h"
#include "poisson.h"
#include "init_data.h"
#include "jacobi_util.h"
#include "matrix_routines.h"

extern double MEMORY;
extern double TIME_SPENT;


// ============================================================================
// 2D TEST FUNCTION

void test_jacobi_2D(int Nx, int Ny)
{
	// Allocation
	double *U = dmalloc_2d_l(Nx, Ny);
	double *F = dmalloc_2d_l(Nx, Ny);
	double *Unew = dmalloc_2d_l(Nx, Ny);
	if (!U || !F || !Unew)
	{
		fprintf(stderr, "Error in malloc, Pointer is NULL.\n");
		return;
	}

	// Initialise the boundary values
	if (strcmp("sin", getenv("PROBLEM_NAME")) == 0)
		init_sin_2D(U, F, Unew, Nx, Ny);
	else if (strcmp("rad", getenv("PROBLEM_NAME")) == 0)
		init_rad_2D(U, F, Unew, Nx, Ny);
	else
	{
		fprintf(stderr, "Problem type is not supported.\n");
		return;
	}

	// Handle the environmental variables
	int maxiter = atoi(getenv("MAX_ITER"));
	double tol = atof(getenv("TOLERANCE"));

	// Main computation and time
	double t = omp_get_wtime();

	jacobi_openmp_2D(Nx, Ny, maxiter, tol, U, F, Unew);

	// Save global variables
	TIME_SPENT = omp_get_wtime() - t;
	MEMORY = 3.0 * Nx * Ny * 8.0 / 1024.0;

	// Print the needed information
	if (strcmp("matrix_full", getenv("OUTPUT_INFO")) == 0)
		array_print_2d(U, Nx, Ny, "%10g ");

	// Free the arrays created for the computation
	free(U);
	free(F);
	free(Unew);
}

// ============================================================================
// 3D TEST FUNCTION

void test_jacobi	(Information *information, char const *solver)
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

	if (strcmp(solver, "omp3d") == 0)
		jacobi_openmp_3D(Nx, Ny, Nz, maxiter, tol, U, F, Unew);
	else if (strcmp(solver, "cuda_1") == 0)
		jacobi_cuda_1(information,maxiter,tol,U,F,Unew);

	// Save global variables
	
	TIME_SPENT = omp_get_wtime() - t;
	MEMORY = 3.0 * Nx * Ny * Nz * 8.0 / 1024.0;

	// Print the needed information
	if (strcmp("matrix_slice", getenv("OUTPUT_INFO")) == 0)
		array_print_3d_slice(U, Nx, Ny, Nz, Nz / 2, "%10g ");
	else if (strcmp("matrix_full", getenv("OUTPUT_INFO")) == 0)
		array_print_3d(U, Nx, Ny, Nz, "%10g ");
	else if (strcmp("error", getenv("OUTPUT_INFO")) == 0)
		compute_global_error(information, A, U);

	// Free the arrays created for the computation
	free(U); free(F); free(A); free(Unew);		
}

// ============================================================================
// MPI TEST FUNCTION

void test_jacobi_mpi(Information *information, char const *solver)
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

	// Initialise the problem
	if (strcmp("sin", getenv("PROBLEM_NAME")) == 0)
		init_sin_mpi3D_2(U, F, Unew, information);
	else {
		fprintf(stderr, "Problem type is not supported.\n");
		return;
	}

	// Handle the environmental variables
	int maxiter = atoi(getenv("MAX_ITER"));
	double tol = atof(getenv("TOLERANCE"));
	
	// Main computation and time
	double t = omp_get_wtime();

	if (     strcmp(solver, "mpi3d_1") == 0)
		jacobi_mpi3D_1(information, maxiter, tol, U, F, Unew);
	else if (strcmp(solver, "mpi3d_2") == 0)
		jacobi_mpi3D_2(information, maxiter, tol, U, F, Unew);

	else if (strcmp(solver, "mixed_1") == 0)
		jacobi_mixed_1(information, maxiter, tol, U, F, Unew);
	else if (strcmp(solver, "mixed_2") == 0)
		jacobi_mixed_2(information, maxiter, tol, U, F, Unew);

	MPI_Barrier(MPI_COMM_WORLD);
	// Save global variables
	TIME_SPENT = omp_get_wtime() - t;
	MEMORY = 3.0 * Nx * Ny * Nz * 8.0 / 1024.0;

	// Print the needed information
	if (     strcmp("matrix_slice", getenv("OUTPUT_INFO")) == 0)
		array_print_3d_slice(U, Nx, Ny, Nz, Nz / 2, "%10g ");
	else if (strcmp("matrix_full", getenv("OUTPUT_INFO")) == 0)
		print_jacobi3d_z_sliced(U, information, "%10g ");
	else if (strcmp("error", getenv("OUTPUT_INFO")) == 0)
		compute_global_error(information, A, U);
	
	MPI_Barrier(MPI_COMM_WORLD);
	// Free the arrays created for the computation
	free(U); free(F); free(A); free(Unew);		
}
