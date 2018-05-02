#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "matrix_routines.h"
#include "cuda_routines.h"

__global__ void jacobi_cuda_iteration(int I, int J, int K,
	double *U, double *F, double *Unew)
{
	// Setting up steps
	double hi = 2.0/(I-1.0);
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