// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

#include "matrix_routines.h"
#include "poisson.h"

extern double FLOP;

// ============================================================================
// JACOBI 3D SOLVER

void jacobi_openmp_3D(int Nx, int Ny, int Nz, int maxit, double threshold,
    double *U, double *F, double *Unew)
{

	// ------------------------------------------------------------------------
	// Preparation

    if(!U || !F || !Unew) {fprintf(stderr,"Pointer is NULL.\n"); return;}

	// Rewritting to C style coordinates
    int i, j, k, I, J, K;
	I = Nz; J = Ny; K = Nx;

	// Setting up steps
    double hi = 2.0/(I-1.0);
	double hj = 2.0/(J-1.0);
	double hk = 2.0/(K-1.0);
    double stepi = hi*hi;
	double stepj = hj*hj;
	double stepk = hk*hk;
	double f3 = 1.0/3.0;
	double f6 = 1.0/6.0;

    // Prepare stop criterion
    bool use_tol = false; int iter = 0; 
	double norm_diff = 10.0;

    if (strcmp("on",getenv("USE_TOLERANCE")) == 0) { use_tol = true; }

	// ------------------------------------------------------------------------
	// Run the iterative method
    for(iter = 0; iter < maxit ; iter++){
        norm_diff = 0.0;

		// Run a single iteration in OpenMP parallel.
        #pragma omp parallel for\
            private(i,j,k) reduction(+: norm_diff) schedule(runtime)
        for(i=1; i < I-1; i++){
            for(j=1; j < J-1; j++){
				for(k=1; k < K-1; k++){
					// Linear indexing with macro
					double ui =
							U[IND_3D(i-1, j, k, I, J, K)] 
						+ U[IND_3D(i+1, j, k, I, J, K)]
						+ f3*stepi*F[IND_3D(i, j, k, I, J, K)];
					double uj =
							U[IND_3D(i, j-1, k, I, J, K)] 
						+ U[IND_3D(i, j+1, k, I, J, K)]
						+ f3*stepj*F[IND_3D(i, j, k, I, J, K)];
					double uk =
							U[IND_3D(i, j, k-1, I, J, K)] 
						+ U[IND_3D(i, j, k+1, I, J, K)]
						+ f3*stepk*F[IND_3D(i, j, k, I, J, K)];
					
					// Collect the terms
					Unew[IND_3D(i, j, k, I, J, K)] = f6*( ui + uj + uk );

					// Compute the stop criterion
					double u    = U[IND_3D(i, j, k, I, J, K)];
					double unew = Unew[IND_3D(i, j, k, I, J, K)];
					norm_diff  += (u - unew)*(u - unew);
				}
            }
        }

        // Swap the arrays and check for convergence
        swap_array( &U, &Unew );
        norm_diff = sqrt(norm_diff);

        if (use_tol && (norm_diff < threshold))
            break;
    }

	// ------------------------------------------------------------------------
	// Finalise

    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0){
        // Print exit cause
        if(norm_diff < threshold && use_tol)
            fprintf(stdout, "Exited because norm < threshold\n");
        else
            fprintf(stdout, "Exited because iter = maxit\n");
    }

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0){
        FLOP = (19.0*I*J*K + 4.0)*iter;
    }
}
// END OF FILE
// ============================================================================