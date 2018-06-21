// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "matrix_routines.h"
#include "jacobi_util.h"
#include "poisson.h"

#include "cuda_routines.h"
#include "jacobi_util_cuda.h"

extern double MFLOP;

// ============================================================================
// JACOBI 3D SOLVER

void jacobi_mixed_1(Information *information, double *U, double *F, double *Unew)
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

	// Check if program can run
	int Ndevices;
	cuda_get_device_count(&Ndevices);
	if ( Ndevices < 1)
	{
		fprintf(stderr, "Version mixed_1 must run on a CUDA device.\n");
		return;
	}
	if ( size < 2)
	{
		fprintf(stderr, "Version mixed_1 must run on multiple mpi ranks.\n");
		return;
	}

	// Set the CUDA device
	cuda_set_device(rank);

	Information *information_cuda;
	cuda_malloc((void**) &information_cuda,sizeof(Information));
	copy_information_cuda(information_cuda,information);
	
	MPI_Request req;

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
    int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	// Allocate device arrays
	double *U_cuda, *F_cuda, *Unew_cuda;
	int arraySizes = I*J*K*sizeof(double);
	cuda_malloc((void**)&U_cuda,    arraySizes);
	cuda_malloc((void**)&F_cuda,    arraySizes);
	cuda_malloc((void**)&Unew_cuda, arraySizes);
	
	// Copy over F and Unew to the device
	copy_to_device(F,   arraySizes,F_cuda   );
	copy_to_device(Unew,arraySizes,Unew_cuda);
	
	// Prepare stop criterion
	int iter = 0;

	// ------------------------------------------------------------------------
	// Run the iterative method
    for(iter = 0; iter < maxit ; iter++){

		// Copy over the result of the last itteration
		copy_to_device(U,   arraySizes,U_cuda   );

		// Compute the iteration of the jacobi method on the GPU
        jacobi_iteration_cuda(
			information, information_cuda, U_cuda, F_cuda, Unew_cuda
		);

		// Copy back the result to the CPU
		copy_from_device(Unew,arraySizes,Unew_cuda);
		
        // Swap the arrays and check for convergence
        swap_array( &U, &Unew );

		// Extract boundaries
		double *U_ptr_s1, *U_ptr_s2;
		double *U_ptr_r1, *U_ptr_r2;
		if (rank == 0) {
			// First rank needs the last updated index
			U_ptr_s1 = &U[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];
			U_ptr_r1 = &U[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			memcpy(s_buf1, U_ptr_s1, N_buffer*sizeof(double));
		} else if (rank == (size - 1)){
			// Last rank needs the first updated index
			U_ptr_s1 = &U[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_r1 = &U[IND_3D(0, 0, 0, I, J, K)];

			memcpy(s_buf1, U_ptr_s1, N_buffer*sizeof(double));
		} else {
			// All other ranks needs the first and last updated index
			U_ptr_s1 = &U[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_s2 = &U[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];

			U_ptr_r1 = &U[IND_3D(0, 0, 0, I, J, K)];
			U_ptr_r2 = &U[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			memcpy(s_buf1, U_ptr_s1, N_buffer*sizeof(double));
			memcpy(s_buf2, U_ptr_s2, N_buffer*sizeof(double));
		}
		
		// Determine source and destination
		int neighbour_1, neighbour_2;
		compute_neighbors(information, &neighbour_1, &neighbour_2);

		MPI_Barrier(MPI_COMM_WORLD); // Maybe remove?

		// Send boundaries and receive boundaries
		MPI_Isend(s_buf1, N_buffer, MPI_DOUBLE, neighbour_1, iter, MPI_COMM_WORLD, &req);
		if ( rank != 0 && rank != (size - 1) )
			MPI_Isend(s_buf2, N_buffer, MPI_DOUBLE, neighbour_2, iter, MPI_COMM_WORLD, &req);

		MPI_Recv(r_buf1, N_buffer, MPI_DOUBLE, neighbour_1, iter, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		if ( rank != 0 && rank != (size - 1) )
			MPI_Recv(r_buf2, N_buffer, MPI_DOUBLE, neighbour_2, iter, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		// Synchronize and swap
		MPI_Barrier(MPI_COMM_WORLD);
		memcpy(U_ptr_r1, r_buf1, N_buffer*sizeof(double));
		if (rank > 0 && rank < (size - 1) )
			memcpy(U_ptr_r2, r_buf2, N_buffer*sizeof(double));

		// Stop early if relative error is used
		/*
		if (information->use_tol)
		{

			// Copy norm_diff^2 from device to host and sqrt it
			copy_from_device(&information->norm_diff, 1*sizeof(double),
				&information_cuda->norm_diff);
			//information->norm_diff = sqrt(information->norm_diff);
			//Information tmp; // ONLY to get tmp.norm_diff, pointers are invalid
			//copy_from_device_void(&tmp, sizeof(Information), information_cuda);
			//information->norm_diff = sqrt(tmp.norm_diff);

			MPI_Allreduce(&information->norm_diff, &information->global_norm_diff,
				1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			if ( (information->global_norm_diff < information->tol) )
				break;
		}
		*/
    }

	// Save number of iterations
	information->iter = iter;

	// ------------------------------------------------------------------------
	// Finalise
	MPI_Barrier(MPI_COMM_WORLD);
	free(s_buf1); free(s_buf2);
	free(r_buf1); free(r_buf2);
	cuda_free(U_cuda   );
	cuda_free(F_cuda   );
	cuda_free(Unew_cuda);
	free_information_cuda(information_cuda);
	
	// Flop Counts:
	// jacobi_iteration_cuda: (iter)
	// 		Simple: 	21
	//		Divisions: 	5
	MFLOP += 1e-6*(21.0 + 5.0*4.0 )*iter*Nx*Ny*Nz;
}
// END OF FILE
// ============================================================================