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

void jacobi_mpi3D_1(Information *information, double *U, double *F, double *Unew)
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

	// Can only run on two cores.
	if (size != 2)
	{
		fprintf(stderr, "Version mpi3d_1 can only run on 2 cores. Exiting.\n");
		return;
	}

	MPI_Request req;

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
    int J, K;
	J = loc_Ny; K = loc_Nx;
	
	// Prepare stop criterion
	int iter = 0; 
    
	// ------------------------------------------------------------------------
	// Run the iterative method
    for(iter = 0; iter < maxit ; iter++){

		// Compute the iteration of the jacobi method
        jacobi_iteration(information, U, F, Unew);

        // Swap the arrays and check for convergence
        swap_array( &U, &Unew );

		// Extract boundaries
		double *U_ptr_s, *U_ptr_r;
		if (rank == 0) {
			U_ptr_s = &U[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];
			U_ptr_r = &U[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];
			memcpy(s_buf, U_ptr_s, N_buffer*sizeof(double));
		} else { // rank == 1 now
			U_ptr_s = &U[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_r = &U[IND_3D(0, 0, 0, I, J, K)];
			memcpy(s_buf, U_ptr_s, N_buffer*sizeof(double));
		}

		// Determine destination
		int dst;
		if (rank == 0) 	{ dst = 1;}
		else 			{ dst = 0;}

		// Send and recieve boundaries
		MPI_Isend(s_buf, N_buffer, MPI_DOUBLE, dst, 0, MPI_COMM_WORLD, &req);
		MPI_Recv(r_buf, N_buffer, MPI_DOUBLE, dst, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// Synchronize and swap
		MPI_Barrier(MPI_COMM_WORLD);
		memcpy(U_ptr_r, r_buf, N_buffer*sizeof(double));

		// Stop early if relative error is used
		if (information->use_tol)
		{
			MPI_Allreduce(&information->norm_diff, &information->global_norm_diff,
				1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			if (information->global_norm_diff < information->tol)
				break;
		}
    }

	// ------------------------------------------------------------------------
	// Finalise
	MPI_Barrier(MPI_COMM_WORLD);
	free(s_buf);
	free(r_buf);

	information->iter = iter;
	
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