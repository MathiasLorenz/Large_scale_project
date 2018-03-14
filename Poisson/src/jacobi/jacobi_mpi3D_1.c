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
#include "poisson.h"

extern double MFLOP;

// ============================================================================
// JACOBI 3D SOLVER

void jacobi_mpi3D_1(int loc_Nx, int loc_Ny, int loc_Nz, int maxit, double threshold, int rank,
    double *U, double *F, double *Unew)
{
	// ------------------------------------------------------------------------
	// Preparation
	// Send and receive buffers for MPI
	double *s_buf = dmalloc_2d_l(loc_Nx, loc_Ny); // Buffer to send with MPI
	double *r_buf = dmalloc_2d_l(loc_Nx, loc_Ny); // Buffer to send with MPI
    if(!U || !F || !Unew || !s_buf || !r_buf) {
		fprintf(stderr,"Pointer is NULL.\n");
		return;
	}

	// Set buffer size
	int N_buffer = loc_Nx*loc_Ny;

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

		// Extract boundaries
		double *U_ptr;
		if (rank == 0) {
			U_ptr = &U[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];
			memcpy(s_buf, U_ptr, N_buffer);
		} else { // rank == 1 now
			U_ptr = &U[IND_3D(1, 0, 0, I, J, K)];
			memcpy(s_buf, U_ptr, N_buffer);
		}

		// Determine source and destination
		int src, dest;
		if (rank == 0) {src = 0; dest = 1;}
		else {src = 1; dest = 0;}

		// Send boundaries
		printf("I'm rank %d before send.\n"
				"N_buffer = %d, s_buf size = %zu\n"
				"src = %d, dest = %d\n",
				rank, N_buffer, sizeof(s_buf), src, dest);
		MPI_Sendrecv(s_buf, N_buffer, MPI_DOUBLE, dest, 0, r_buf, N_buffer,
			MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("I'm rank %d after send.\n", rank);
		MPI_Barrier(MPI_COMM_WORLD);

		// Insert received boundaries
		memcpy(U_ptr, r_buf, N_buffer);

		// Remember to implement tolerance
		/*
        norm_diff = sqrt(norm_diff);

        if (use_tol && (norm_diff < threshold))
            break;
		*/
    }

	// ------------------------------------------------------------------------
	// Finalise

	free(s_buf);
	free(r_buf);
	
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