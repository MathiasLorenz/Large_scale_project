// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#include "matrix_routines.h"
#include "jacobi_util.h"

// ============================================================================
// FUNCTIONS TO HANDLE THE INFORMATION STRUCTURE

void write_information(Information *information, int Nx, int Ny, int Nz,
	int rank, int size)
{
	information->size = size;
	information->rank = rank;
	information->global_Nx = Nx;
	information->global_Ny = Ny;
	information->global_Nz = Nz;
	information->loc_Nx = malloc(size*sizeof(int));
	information->loc_Ny = malloc(size*sizeof(int));
	information->loc_Nz = malloc(size*sizeof(int));
	if (!information->loc_Nx || !information->loc_Ny || !information->loc_Nz)
	{ fprintf(stderr, "Error in malloc, pointer is NULL.\n"); return; }

	if (size == 1)
	{
		information->loc_Nx[0] = Nx;
		information->loc_Ny[0] = Ny;
		information->loc_Nz[0] = Nz;
	} else {
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

	// Handle the environmental variables
	information->maxit	= atoi(getenv("MAX_ITER"));
	information->iter	= -1;
	information->tol	= atof(getenv("TOLERANCE"));

	// Variables for tolerance calculations
	if (strcmp("on",getenv("USE_TOLERANCE")) == 0) 
		information->use_tol = true;
	else
		information->use_tol = false;

	information->local_frobenius = 0.0;
	information->frobenius_error = 10.0;
}

void free_information_arrays(Information *information)
{
	free(information->loc_Nx);
	free(information->loc_Ny);
	free(information->loc_Nz);
}

// ============================================================================
// Function to compute neighbors for mpi threads.
void compute_neighbors(Information *information,
	int *neighbour_1, int *neighbour_2)
{
	int rank   = information->rank;
	int size   = information->size;
	if (rank == 0) {
		*neighbour_1 = 1;
	} else if (rank == size - 1) {
		*neighbour_1 = size - 2;
	} else {
		*neighbour_1 = rank - 1; 
		*neighbour_2 = rank + 1;
	}
}

// ============================================================================
// JACOBI ITERATION

void jacobi_iteration(Information *information,
					double *U, double *F, double *Unew)
{
	int rank = information->rank;
	int Nx = information->global_Nx;
	int Ny = information->global_Ny;
	int Nz = information->global_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];

    int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	// Setting up steps
    double hi = 2.0/(Nz-1.0);
	double hj = 2.0/(Ny-1.0);
	double hk = 2.0/(Nx-1.0);
    double stepi = hi*hi;
	double stepj = hj*hj;
	double stepk = hk*hk;
	double f3 = 1.0/3.0;
	double f6 = 1.0/6.0;

	// For relative error stopping
	information->local_frobenius = 0.0;

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

				// Tolerance criterion
				// For small problems it is more efficient to put outside the loop,
				// however for large problems (as we wish to solve) looping only once
				// is more efficient.				
				if (information->use_tol)
				{
					double uij    = U[ijk];
					double unewij = Unew[ijk];
					information->local_frobenius += (uij - unewij)*(uij - unewij);
				}
			}
		}
	}
}

// ============================================================================
// JACOBI ITERATION with interior and boundary split

void jacobi_iteration_separate(Information *information,
					double *U, double *F, double *Unew, const char *ver)
{
	int rank = information->rank;
	int Nx = information->global_Nx;
	int Ny = information->global_Ny;
	int Nz = information->global_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];

    int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	// Setting up steps
    double hi = 2.0/(Nz-1.0);
	double hj = 2.0/(Ny-1.0);
	double hk = 2.0/(Nx-1.0);
    double stepi = hi*hi;
	double stepj = hj*hj;
	double stepk = hk*hk;
	double f3 = 1.0/3.0;
	double f6 = 1.0/6.0;

	// Loop over points. Either interior or boundary
	if (strcmp(ver, "i") == 0) // interior
	{
		// Loop over all interior points
		for (int i = 2; i < I - 2; i++) {
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

					// Tolerance criterion
					// For small problems it is more efficient to put outside the loop,
					// however for large problems (as we wish to solve) looping only once
					// is more efficient.				
					if (information->use_tol)
					{
						double uij    = U[ijk];
						double unewij = Unew[ijk];
						information->local_frobenius += (uij - unewij)*(uij - unewij);
					}
				}
			}
		}
	}
	if (strcmp(ver, "b") == 0) // boundary
	{
		for (int b = 0; b < 2; b++)
		{
			// Determine which boundary we are on
			int i = ( b == 0) ? 1 : I - 2;

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

					// Tolerance criterion
					// For small problems it is more efficient to put outside the loop,
					// however for large problems (as we wish to solve) looping only once
					// is more efficient.				
					if (information->use_tol)
					{
						double uij    = U[ijk];
						double unewij = Unew[ijk];
						information->local_frobenius += (uij - unewij)*(uij - unewij);
					}
				}
			}
		}
	}
}

// ============================================================================
// COMPUTE MAXIMAL ABSOLUTE DIFFERENCE ERROR

// Function to print error for OUTPUT_INFO=error
void print_error(Information *information, double *U)
{
	int rank = information->rank;
	int Nx = information->global_Nx;
	int Ny = information->global_Nx;
	int Nz = information->global_Nx;

	int iter = information->iter;
	double frobenius_error = information->frobenius_error;

	double global_error = 0.0;
	compute_global_max_error(information, U, &global_error);

	if (rank == 0)
	{
		printf("Grid: %d %d %d, ",Nx, Ny, Nz);
		if (information->use_tol)
			printf("iter: %d, frobenius error: %.7e ", iter, frobenius_error);
		printf("Max absolute error: %.7e\n",global_error);
		
	}
}

// Compute max absolute error
void compute_local_max_error(Information *information, double *U, double *local_error)
{
	if (!U || !local_error || !information)
	{ fprintf(stderr, "Pointer is NULL.\n"); return; }
	*local_error = 0.0;

	// Extract problem dimensions
	int rank = information->rank;
	int Nx = information->global_Nx;
	int Ny = information->global_Ny;
	int Nz = information->global_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];
	// Rewritting to C style coordinates
	int K = loc_Nx, J = loc_Ny, I = loc_Nz;

	// Setting up steps and variables 
	double hi = 2.0 / (Nz - 1.0);
	double hj = 2.0 / (Ny - 1.0);
	double hk = 2.0 / (Nx - 1.0);

	// Determine how far we are in the z-direction
	double z = -1.0;
	for (int r = 0; r < rank; r++)
		z += hi*(information->loc_Nz[r]-2.0);

	// Loop over all interior points
	double x, y, true_sol, err;
	for (int i = 0; i < I; i++)
	{
		y = -1.0;
		for (int j = 0; j < J; j++)
		{
			x = -1.0;
			for (int k = 0; k < K; k++)
			{
				int ijk = IND_3D(i, j, k, I, J, K);
				true_sol = sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
				err = fabs(U[ijk] - true_sol);
				if (err > *local_error) *local_error = err; // Save largest element
				x += hk;
			}
			y += hj;
		}
		z += hi;
	}
}
// Function to compute global error on solution.
void compute_global_max_error(Information *information, double *U,
	double *global_error)
{
	*global_error = 0.0;
	double local_error = 0.0;

	compute_local_max_error(information, U, &local_error);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&local_error, global_error, 1, MPI_DOUBLE, MPI_MAX,
		0, MPI_COMM_WORLD);
}


// ============================================================================
// COMPUTE STOP CRITERION BASED ON FROBENIUS

// Calculate if norm criterion is reached
bool norm_early_stop(Information *information)
{
	int size = information->size;
	
	// Reduce and send global norm diff to all threads
	if (size > 1)
		MPI_Allreduce(
			&information->local_frobenius, &information->frobenius_error,
			1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	else
		information->frobenius_error = information->local_frobenius;

	information->local_frobenius = 0.0;
	information->frobenius_error = sqrt(information->frobenius_error);
	
	// Returning boolean for stopping
	return (information->frobenius_error < information->tol);
}


// END OF FILE
// ============================================================================