// ============================================================================
// INCLUDES

#include <stdio.h>
#include <math.h>

#include "matrix_routines.h"
#include "init_data.h"

// ============================================================================
// SINUSOIDAL PROBLEM IN 2D

void init_sin_2D(double *U, double *F, double *Unew, int Nx, int Ny)
{
	if (!U || !F || !Unew)
	{
		fprintf(stderr, "Pointer is NULL.\n");
		return;
	}

	// Rewritting to C style coordinates
	int I = Ny, J = Nx;

	// Setting up steps and variables
	double hi = 2.0 / (I - 1.0);
	double hj = 2.0 / (J - 1.0);
	double M_PI_sq = M_PI*M_PI;

	double y = -1.0;

	// Initialising arrays
	for (int i = 0; i < I; i++)
	{
		double x = -1.0;
		for (int j = 0; j < J; j++)
		{
			F[IND_2D(i, j, I, J)] = 2.0 * M_PI_sq * sin(M_PI * x) * sin(M_PI * y);
			x += hj;
		}
		y += hi;
	}

	// Set U and Unew to 0
	for (int i = 0; i < I*J; i++)
	{
		U[i] = 0.0;
		Unew[i] = 0.0;
	}
}

// ============================================================================
// RADIATOR PROBLEM

void init_rad_2D(double *U, double *F, double *Unew, int Nx, int Ny)
{
	if (!U || !F || !Unew) { fprintf(stderr, "Pointer is NULL.\n"); return; }

	int I = Ny, J = Nx;

	// Setting up steps and variables
	double hi = 2.0 / (I - 1.0);
	double hj = 2.0 / (J - 1.0);

	// Set up F to be the radiator problem.
	for (int i = 0; i < I*J; i++) 
		F[i] = 0.0;
	for (int i = lround(1.0 / (hi * 3.0)); i <= lround(2.0 / (hi * 3.0)); i++)
		for (int j = lround(1.0 / hj); j <= lround(4.0 / (hj * 3.0)); j++)
		{
			F[IND_2D(i, j, I, J)] = 200;
		}

	// Set U to zero
	for (int i = 0; i < I*J; i++) 
		{ U[i] = 0.0; Unew[i] = 0.0; }

	// Set BC of U
	for (int i = 0; i < I; i++)
	{
		U[IND_2D(i,0,I,J)] = Unew[IND_2D(i,0,I,J)] = 20;
		U[IND_2D(i,J-1,I,J)] = Unew[IND_2D(i,J-1,I,J)] = 20;
	}
	for (int j = 0; j < J; j++)
		U[IND_2D(I-1,j,I,J)] = Unew[IND_2D(I-1,j,I,J)] = 20;

}

// ============================================================================
// SINUSOIDAL PROBLEM IN 3D

void init_sin_3D(double *U, double *F, double *Unew, int Nx, int Ny, int Nz)
{
	if (!U || !F || !Unew) { fprintf(stderr, "Pointer is NULL.\n"); return; }

	// Rewritting to C style coordinates
	int K = Nx, J = Ny, I = Nz;

	// Setting up steps and variables
	double hi = 2.0 / (I - 1.0);
	double hj = 2.0 / (J - 1.0);
	double hk = 2.0 / (K - 1.0);
	double M_PI_sq = M_PI*M_PI;
	
	// Initialising arrays
	// #pragma omp parallel for private(i, j, x, y) shared(M_PI_sq, F, U, Unew)
	double z = -1.0;
	for (int i = 0; i < I; i++)
	{
		double y = -1.0;
		for (int j = 0; j < J; j++)
		{
			double x = -1.0;
			for (int k = 0; k < K; k++)
			{
				F[IND_3D(i, j, k, I, J, K)] = 3.0 * M_PI_sq * 
					sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
				x += hk;
			}
			y += hj;
		}
		z += hi;
	}

	// Set U and Unew to 0
	for (int i = 0; i < I*J*K; i++)
	{
		U[i] = 0.0;
		Unew[i] = 0.0;
	}
}

// END OF FILE
// ============================================================================
