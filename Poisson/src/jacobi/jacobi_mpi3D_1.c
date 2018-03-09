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
#include "poisson.h"

extern double MFLOP;

// ============================================================================
// JACOBI 3D SOLVER

void jacobi_mpi3D_1(int loc_Nx, int loc_Ny, int loc_Nz, int maxit, double threshold, int rank,
    double *U, double *F, double *Unew)
{

	// ------------------------------------------------------------------------
	// Preparation

    if(!U || !F || !Unew) {fprintf(stderr,"Pointer is NULL.\n"); return;}

	// Rewritting to C style coordinates
    int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	// Remember to implement tolerance
	/*
    // Prepare stop criterion
    bool use_tol = false;
	double norm_diff = 10.0;
	
	if (strcmp("on",getenv("USE_TOLERANCE")) == 0) { use_tol = true; }
	*/
	
	// Prepare stop criterion
	int iter = 0; 

    

	// ------------------------------------------------------------------------
	// Run the iterative method
    for(iter = 0; iter < maxit ; iter++){
		// Remember to implement tolerance
		/*
        norm_diff = 0.0;
		*/

		// Compute the iteration of the jacobi method
        jacobi_iteration(I, J, K, rank, U, F, Unew);

        // Swap the arrays and check for convergence
        swap_array( &U, &Unew );

		// Remember to implement tolerance
		/*
        norm_diff = sqrt(norm_diff);

        if (use_tol && (norm_diff < threshold))
            break;
		*/
    }

	// ------------------------------------------------------------------------
	// Finalise
	
	MFLOP = 1e-6*(19.0*I*J*K + 4.0)*iter;

	// Print the information requested
    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0){
		fprintf(stdout, "Exited because iter = maxit\n");

		/*
		// Remember to implement tolerance
        if(norm_diff < threshold && use_tol)
            fprintf(stdout, "Exited because norm < threshold\n");
        else
            fprintf(stdout, "Exited because iter = maxit\n");
		*/
    }
}
// END OF FILE
// ============================================================================