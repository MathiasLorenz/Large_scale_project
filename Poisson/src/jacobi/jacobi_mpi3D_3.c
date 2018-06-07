// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "matrix_routines.h"
#include "jacobi_util.h"
#include "poisson.h"

extern double MFLOP;

// ============================================================================
// JACOBI 3D SOLVER
// This version (intents) to hide latency by computing interior points while
// the boundary is send between processes.

void jacobi_mpi3D_3(Information *information, double *U, double *F, double *Unew)
{
	// Read the information structure
	int rank   = information->rank;
	int size   = information->size;
	int Nx	   = information->global_Nx;
	int Ny	   = information->global_Ny;
	int Nz	   = information->global_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];
	int maxit  = information->maxit;
	
	// Must run on more than one core
	if (size < 2)
	{
		fprintf(stderr, "Version mpi3d_3 must run on multiple cores. Exiting.\n");
		return;
	}

	// Number of requests for MPI send/recv. 
	int num_req = (rank > 1 && rank < (size - 1)) ? 4 : 2;
	MPI_Request req[num_req];

	// ------------------------------------------------------------------------
	// Preparation
	// Send and receive buffers for MPI
	double *s_buf1 = dmalloc_2d_l(loc_Nx, loc_Ny); // Buffer 1 to send with MPI
	double *r_buf1 = dmalloc_2d_l(loc_Nx, loc_Ny); // Buffer 1 to send with MPI
	double *s_buf2 = dmalloc_2d_l(loc_Nx, loc_Ny); // Buffer 2 to send with MPI
	double *r_buf2 = dmalloc_2d_l(loc_Nx, loc_Ny); // Buffer 2 to send with MPI
    if(!U || !F || !Unew || !s_buf1 || !r_buf1 || !s_buf2 || !r_buf2) {
		fprintf(stderr,"Pointer is NULL.\n");
		return;
	}

	// Set buffer size
	int N_buffer = loc_Nx*loc_Ny;

	// Rewritting to C style coordinates
    int J, K;
	J = loc_Ny; K = loc_Nx;
	
	// Prepare stop criterion
	int iter = 0;

	// ------------------------------------------------------------------------
	// Run the iterative method
    for(iter = 0; iter < maxit ; iter++){

		// Compute the boundary first so we can send it while computing int
		jacobi_iteration_separate(information, U, F, Unew, "b");

		// Extract boundaries
		double *U_ptr_s1, *U_ptr_s2;
		double *U_ptr_r1, *U_ptr_r2;
		if (rank == 0) {
			// First rank needs the last updated index
			U_ptr_s1 = &Unew[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];
			U_ptr_r1 = &Unew[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			memcpy(s_buf1, U_ptr_s1, N_buffer*sizeof(double));
		} else if (rank == (size - 1)){
			// Last rank needs the first updated index
			U_ptr_s1 = &Unew[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_r1 = &Unew[IND_3D(0, 0, 0, I, J, K)];

			memcpy(s_buf1, U_ptr_s1, N_buffer*sizeof(double));
		} else {
			// All other ranks needs the first and last updated index
			U_ptr_s1 = &Unew[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_s2 = &Unew[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];

			U_ptr_r1 = &Unew[IND_3D(0, 0, 0, I, J, K)];
			U_ptr_r2 = &Unew[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			memcpy(s_buf1, U_ptr_s1, N_buffer*sizeof(double));
			memcpy(s_buf2, U_ptr_s2, N_buffer*sizeof(double));
		}
		
		// Determine source and destination
		int neighbour_1, neighbour_2;
		compute_neighbors(information, &neighbour_1, &neighbour_2);

		// Send boundaries and receive boundaries
		MPI_Isend(s_buf1, N_buffer, MPI_DOUBLE, neighbour_1, 0, MPI_COMM_WORLD,
				&req[0]);
		if ( rank != 0 && rank != (size - 1) )
			MPI_Isend(s_buf2, N_buffer, MPI_DOUBLE, neighbour_2, 0, MPI_COMM_WORLD,
					&req[2]);

		MPI_Irecv(r_buf1, N_buffer, MPI_DOUBLE, neighbour_1, 0,
					MPI_COMM_WORLD, &req[1]);
		if ( rank != 0 && rank != (size - 1) )
			MPI_Irecv(r_buf2, N_buffer, MPI_DOUBLE, neighbour_2, 0,
					MPI_COMM_WORLD, &req[3]);

		// Compute interior while boundary is being sent
		jacobi_iteration_separate(information, U, F, Unew, "i");

		MPI_Waitall(num_req, req, MPI_STATUS_IGNORE);

		// Synchronize and copy buffers
		memcpy(U_ptr_r1, r_buf1, N_buffer*sizeof(double));
		if (rank > 0 && rank < (size - 1) )
			memcpy(U_ptr_r2, r_buf2, N_buffer*sizeof(double));

		// Swap the arrays
        swap_array( &U, &Unew );

		// Stop early if relative error is used
		MPI_Allreduce(&information->norm_diff, &information->global_norm_diff,
			1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		if ( (information->use_tol) && (information->global_norm_diff < information->tol) )
			break;
    }

	information->iter = iter;

	// ------------------------------------------------------------------------
	// Finalise
	MPI_Barrier(MPI_COMM_WORLD); // needed?
	free(s_buf1); free(s_buf2);
	free(r_buf1); free(r_buf2);
	
	// Flop Counts:
	// jacobi_iteration: (iter)
	//		Constants: 
	//			Simple:		6
	//			Divisions:	5
	// 		Update:
	//			Simple:		15
	MFLOP += 1e-6*( (6.0 + 5.0*4.0 ) + 15.0*Nx*Ny*Nz)*iter;
}
// END OF FILE
// ============================================================================