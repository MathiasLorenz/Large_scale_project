// ============================================================================
// INCLUDES
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "tests_cuda.h"
#include "cuda_routines.h"

// ============================================================================
// CUDA TEST
void test_cuda(int Nx, int Ny, int Nz)
{
	// Allocation
	/*
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
    */

	// Handle the environmental variables
    int maxiter = atoi(getenv("MAX_ITER"));
    double tol  = atof(getenv("TOLERANCE"));

	// Main computation and time
	double t = omp_get_wtime();
    printf("We are now testing the cuda function\n");
	cu_test<<<1,1>>>();
	checkCudaErrors(cudaDeviceSynchronize());

	// Save global variables
	TIME_SPENT 	= omp_get_wtime() - t;
	MEMORY		= 3.0*Nx*Ny*Nz*8.0/1024.0;

	/*
	// Print the needed information
    if (strcmp("matrix", getenv("OUTPUT_INFO")) == 0)
        array_print_3d_slice(U, Nx, Ny, Nz, Nz/2, "%10g ");
    else if (strcmp("full_matrix", getenv("OUTPUT_INFO")) == 0)
        array_print_3d(U, Nx, Ny, Nz, "%10g ");

	// Free the arrays created for the computation
    free(U);
    free(F);
    free(Unew);
	*/
}

