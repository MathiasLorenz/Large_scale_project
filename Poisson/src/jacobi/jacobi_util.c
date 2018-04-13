// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

#include "matrix_routines.h"
#include "jacobi_util.h"

// ============================================================================
// FUNCTIONS TO HANDLE THE INFORMATION STRUCTURE

void write_information(Information *information, int Nx, int Ny, int Nz,
	int rank, int size)
{
	information->size = size;
	information->rank = rank;
	information->glo_Nx = Nx;
	information->glo_Ny = Ny;
	information->glo_Nz = Nz;
	information->loc_Nx = malloc(size*sizeof(int));
	information->loc_Ny = malloc(size*sizeof(int));
	information->loc_Nz = malloc(size*sizeof(int));
	int loc_Nz;
	for (int r = 0; r < size; r++)
	{
		if 		(r == 0) 				{loc_Nz = floor(Nz/size) + 1; }
		else if (0 < r && r < (size-1)) {loc_Nz = floor(Nz/size) + 2; }
		else 							{loc_Nz = Nz - (size - 1)*floor(Nz/size) + 1;}

		information->loc_Nx[r] = Nx;
		information->loc_Ny[r] = Ny;
		information->loc_Nz[r] = loc_Nz;
	}
}

void free_information_arrays(Information *information)
{
	free(information->loc_Nx);
	free(information->loc_Ny);
	free(information->loc_Nz);
}

// ============================================================================
// FUNCTION FOR GENERATING AND CHECKING SOLUTION AGAINST TRUE SOLUTION
void generate_true_solution(double *A, Information *information)
{
	if (!A || !information) { fprintf(stderr, "Pointer is NULL.\n"); return; }

	// Read information
	int size = information->size;
	int rank = information->rank;
	int Nx = information->glo_Nx;
	int Ny = information->glo_Ny;
	int Nz = information->glo_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];

	// Rewritting to C style coordinates
	int K = loc_Nx, J = loc_Ny, I = loc_Nz;

	// Setting up steps and variables 
	double hi = 2.0 / (Nz - 1.0); // Global Nz is intended
	double hj = 2.0 / (J - 1.0);
	double hk = 2.0 / (K - 1.0);

	// This is based on an offset where z is split once.
	double z = -1.0;
	for (int r = 1; r <= rank; r++)
		z += hi*(information->loc_Nz[r-1]-2.0);

	for (int i = 0; i < I; i++)
	{
		double y = -1.0;
		for (int j = 0; j < J; j++)
		{
			double x = -1.0;
			for (int k = 0; k < K; k++)
			{
				A[IND_3D(i, j, k, I, J, K)] =
					sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
				x += hk;
			}
			y += hj;
		}
		z += hi;
	}
}

// Compute absolute error
void check_true_solution(double *A, double *U, double *abs_err,
	Information *information)
{
	if (!A || !U || !abs_err || !information)
	{ fprintf(stderr, "Pointer is NULL.\n"); return; }
	*abs_err = 0.0;

	// Extract problem dimensions
	int rank = information->rank;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];
	// Rewritting to C style coordinates
	int K = loc_Nx, J = loc_Ny, I = loc_Nz;

	// Consider not including boundary?
	for (int i = 0; i < I; ++i)
		for (int j = 0; j < J; ++j)
			for (int k = 0; k < K; ++k) {
				int ijk = IND_3D(i, j, k, I, J, K);
				double u = U[ijk];
				double a = A[ijk];
				double err = fabs(u - a);
				if (err > *abs_err) *abs_err = err;
			}
}



// ============================================================================
// JACOBI ITERATION

void jacobi_iteration(int I, int J, int K, int rank, int global_Nz,
					double *U, double *F, double *Unew)
{
	// Setting up steps
    double hi = 2.0/(global_Nz-1.0);
	double hj = 2.0/(J-1.0);
	double hk = 2.0/(K-1.0);
    double stepi = hi*hi;
	double stepj = hj*hj;
	double stepk = hk*hk;
	double f3 = 1.0/3.0;
	double f6 = 1.0/6.0;

	// Loop over all interior points
	for (int i = 1; i < I - 1; i++) {
		for (int j = 1; j < J - 1; j++) {
			for (int k = 1; k < K - 1; k++) {
				// Save i, j, k index once
				int ijk = IND_3D(i, j, k, I, J, K);

				// Linear indexing with macro
				double ui = U[IND_3D(i - 1, j, k, I, J, K)] 
					+ U[IND_3D(i + 1, j, k, I, J, K)] 
					+ f3 * stepi * F[ijk];
				double uj = U[IND_3D(i, j - 1, k, I, J, K)] 
					+ U[IND_3D(i, j + 1, k, I, J, K)] 
					+ f3 * stepj * F[ijk];
				double uk = U[IND_3D(i, j, k - 1, I, J, K)] 
					+ U[IND_3D(i, j, k + 1, I, J, K)] 
					+ f3 * stepk * F[ijk];

				// Collect terms
				Unew[ijk] = f6 * (ui + uj + uk);

				// Compute the stop criterion
				// Remember to implement tolerance
				/*
				double u = U[IND_3D(i, j, k, I, J, K)];
				double unew = Unew[IND_3D(i, j, k, I, J, K)];
				norm_diff  += (u - unew)*(u - unew);
				*/
			}
		}
	}
}
// END OF FILE
// ============================================================================