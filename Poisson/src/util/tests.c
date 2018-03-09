// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "tests.h"
#include "matrix_routines.h"
#include "poisson.h"
#include "init_data.h"

extern double MEMORY;
extern double TIME_SPENT;

// ============================================================================
// JACOBI 2D TEST

void test_jacobi_2D(int Nx, int Ny)
{
	// Allocation
    double *U = dmalloc_2d_l(Nx, Ny);
    double *F = dmalloc_2d_l(Nx, Ny);
    double *Unew = dmalloc_2d_l(Nx, Ny);
    if(!U || !F || !Unew)
        {fprintf(stderr,"Error in malloc, Pointer is NULL.\n"); return;}

	// Initialise the boundary values
    if (strcmp("sin",getenv("PROBLEM_NAME")) == 0)
        init_sin_2D(U, F, Unew, Nx, Ny);
    else if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
        init_rad_2D(U, F, Unew, Nx, Ny);
    else
        {fprintf(stderr,"Error in problem specification.\n"); return;}
    
	// Handle the environmental variables
    int maxiter = atoi(getenv("MAX_ITER"));
    double tol  = atof(getenv("TOLERANCE"));

	// Main computation and time
	double t = omp_get_wtime();

    jacobi_openmp_2D(Nx, Ny, maxiter, tol, U, F, Unew);

	// Save global variables
	TIME_SPENT 	= omp_get_wtime() - t;
	MEMORY		= 3.0*Nx*Ny*8.0/1024.0;

	// Print the needed information
	if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0)
        array_print_2d(U, Nx, Ny, "%10g ");

	// Free the arrays created for the computation
    free(U);
    free(F);
    free(Unew);
}

// ============================================================================
// JACOBI 3D TEST

void test_jacobi_3D(int Nx, int Ny, int Nz)
{
	// Allocation
    double *U = dmalloc_3d_l(Nx, Ny, Nz);
    double *F = dmalloc_3d_l(Nx, Ny, Nz);
    double *Unew = dmalloc_3d_l(Nx, Ny, Nz);
    if(!U || !F || !Unew)
        {fprintf(stderr,"Error in malloc, Pointer is NULL.\n"); return;}

	// Initialise the boundary values
    if (strcmp("sin",getenv("PROBLEM_NAME")) == 0)
        init_sin_3D(U, F, Unew, Nx, Ny, Nz);
    // else if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
    //    init_rad_2D(U, F, Unew, Nx, Ny);
    else
        {fprintf(stderr,"Error in problem specification.\n"); return;}
    
	// Handle the environmental variables
    int maxiter = atoi(getenv("MAX_ITER"));
    double tol  = atof(getenv("TOLERANCE"));

	// Main computation and time
	double t = omp_get_wtime();
    jacobi_openmp_3D(Nx, Ny, Nz, maxiter, tol, U, F, Unew);

	// Save global variables
	TIME_SPENT 	= omp_get_wtime() - t;
	MEMORY		= 3.0*Nx*Ny*Nz*8.0/1024.0;

	// Print the needed information
    if (strcmp("matrix", getenv("OUTPUT_INFO")) == 0)
        array_print_3d_slice(U, Nx, Ny, Nz, Nz/2, "%10g ");
    else if (strcmp("full_matrix", getenv("OUTPUT_INFO")) == 0)
        array_print_3d(U, Nx, Ny, Nz, "%10g ");

	// Free the arrays created for the computation
    free(U);
    free(F);
    free(Unew);
}

// ============================================================================
// 3D MPI TEST 1

void test_jacobi_mpi3D_1(int Nx, int Ny, int Nz)
{
	MPI_Comm comm = MPI_COMM_WORLD;
	int size, rank;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	// Main computation and time
	double t = omp_get_wtime();
    printf("Hi, I'm %d out of %d\n", rank, size);

	// Save global variables
	TIME_SPENT 	= omp_get_wtime() - t;
	MEMORY		= 3.0*Nx*Ny*Nz*8.0/1024.0;
}