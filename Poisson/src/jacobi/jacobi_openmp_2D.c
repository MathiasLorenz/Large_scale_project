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
// JACOBI 2D SOLVER

void jacobi_openmp_2D(int Nx, int Ny, int maxit, double threshold,
    double *U, double *F, double *Unew)
{

	// ------------------------------------------------------------------------
	// Preparation

    if(!U || !F || !Unew) {fprintf(stderr,"Pointer is NULL.\n"); return;}

	// Rewritting to C style coordinates
    int i, j;
	int I,J;
	I = Ny;
	J = Nx;

	// Setting up steps
    double hi = 2.0/(I-1.0);
	double hj = 2.0/(J-1.0);
    double stepi = hi*hi;
	double stepj = hj*hj;

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
            private(i,j) reduction(+: norm_diff) schedule(runtime)
        for(i=1; i < I-1; i++){
            for(j=1; j < J-1; j++){

                // Linear indexing with macro
				double ui =
					  U[IND_2D(i-1, j, I, J)] 
					+ U[IND_2D(i+1, j, I, J)]
					+ 0.5*stepi*F[IND_2D(i, j, I, J)];
				double uj =
					  U[IND_2D(i, j-1, I, J)] 
					+ U[IND_2D(i, j+1, I, J)]
					+ 0.5*stepj*F[IND_2D(i, j, I, J)];
				// Collect the terms
                Unew[IND_2D(i, j, I, J)] = 0.25*( ui + uj );

				// Compute the stop criterion
                double uij   = U[IND_2D(i, j, I, J)];
                double unewij = Unew[IND_2D(i, j, I, J)];
                norm_diff += (uij - unewij)*(uij - unewij);
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
        FLOP = (14.0*I*J + 4.0)*iter;
    }
}
// END OF FILE
// ============================================================================