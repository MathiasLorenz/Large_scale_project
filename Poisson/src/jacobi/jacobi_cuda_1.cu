#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "poisson_cuda.h"
#include "matrix_routines.h"
#include "cuda_routines.h"
#include "jacobi_util_cuda.h"

extern double MFLOP;

void jacobi_cuda_1(int Nx, int Ny, int Nz, int maxit, double threshold, double *U, double *F, double *Unew)
{
    if(!U || !F || !Unew) {
		fprintf(stderr,"Pointer is NULL.\n");
		return;
	}

	// Rewritting to C style coordinates
    int I, J, K;
	I = Nz; J = Ny; K = Nx;

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
        
		jacobi_cuda_iteration<<<1,1>>>(I, J, K, U, F, Unew);
		checkCudaErrors(cudaDeviceSynchronize());

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