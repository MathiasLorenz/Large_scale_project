// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

#include "matrix_routines.h"
#include "jacobi_util.h"
#include "poisson.h"

#include "cuda_routines.h"
#include "jacobi_util_cuda.h"

extern double MFLOP;

// ============================================================================
// JACOBI 3D SOLVER

void jacobi_mixed_2(Information *information, double *U, double *F, double *Unew)
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

	// Set the CUDA device
	setCudaDevice(rank);
	
	Information *information_cuda;
	cuda_malloc((void**)&information_cuda,sizeof(Information));
	copy_information_cuda(information_cuda, information);
	
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

	// Allocate and initialize device arrays
	double *U_cuda, *F_cuda, *Unew_cuda;
	int arraySizes = I*J*K*sizeof(double);
	cuda_malloc((void**)&U_cuda,    arraySizes);
	cuda_malloc((void**)&F_cuda,    arraySizes);
	cuda_malloc((void**)&Unew_cuda, arraySizes);

	copy_to_device_async(U,   arraySizes,U_cuda   );
	copy_to_device_async(F,   arraySizes,F_cuda   );
	copy_to_device_async(Unew,arraySizes,Unew_cuda);
	cuda_synchronize();

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
        jacobi_iteration_cuda(
			information, information_cuda, U_cuda, F_cuda, Unew_cuda
		);

        // Swap the arrays
		cuda_synchronize();
        swap_array( &U_cuda, &Unew_cuda );

		// Extract boundaries
		double *U_ptr_s1, *U_ptr_s2;
		double *U_ptr_r1, *U_ptr_r2;
		if (rank == 0) {
			// First rank needs the last updated index
			U_ptr_s1 = &U_cuda[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];
			U_ptr_r1 = &U_cuda[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			copy_from_device_async(s_buf1, N_buffer*sizeof(double), U_ptr_s1);
		} else if (rank == (size - 1)){
			// Last rank needs the first updated index
			U_ptr_s1 = &U_cuda[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_r1 = &U_cuda[IND_3D(0, 0, 0, I, J, K)];

			copy_from_device_async(s_buf1, N_buffer*sizeof(double), U_ptr_s1);
		} else {
			// All other ranks needs the first and last updated index
			U_ptr_s1 = &U_cuda[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_s2 = &U_cuda[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];

			U_ptr_r1 = &U_cuda[IND_3D(0, 0, 0, I, J, K)];
			U_ptr_r2 = &U_cuda[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			copy_from_device_async(s_buf1, N_buffer*sizeof(double), U_ptr_s1);
			copy_from_device_async(s_buf2, N_buffer*sizeof(double), U_ptr_s2);
		}
		
		// Determine source and destination
		int neighbour_1, neighbour_2;
		compute_neighbors(information, &neighbour_1, &neighbour_2);

		cuda_synchronize(); // Maybe remove?

		// Send boundaries and receive boundaries
		MPI_Isend(s_buf1, N_buffer, MPI_DOUBLE, neighbour_1, 0, MPI_COMM_WORLD, &req);
		if ( rank != 0 && rank != (size - 1) )
			MPI_Isend(s_buf2, N_buffer, MPI_DOUBLE, neighbour_2, 0, MPI_COMM_WORLD, &req);

		MPI_Recv(r_buf1, N_buffer, MPI_DOUBLE, neighbour_1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		if ( rank != 0 && rank != (size - 1) )
			MPI_Recv(r_buf2, N_buffer, MPI_DOUBLE, neighbour_2, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		// Synchronize and copy back
		MPI_Barrier(MPI_COMM_WORLD);
		copy_to_device_async(r_buf1, N_buffer*sizeof(double), U_ptr_r1);
		if (rank > 0 && rank < (size - 1) )
			copy_to_device_async(r_buf2, N_buffer*sizeof(double), U_ptr_r2);
		
		cuda_synchronize();
		// Remember to implement tolerance
		/*
        norm_diff = sqrt(norm_diff);

        if (use_tol && (norm_diff < threshold))
            break;
		*/
    }

	// Save number of iterations
	information->iter = iter;

	// ------------------------------------------------------------------------
	// Finalise
	MPI_Barrier(MPI_COMM_WORLD);

	// Copy back the result
	copy_from_device_async(U,   arraySizes,U_cuda   );
	cuda_synchronize();

	// Free the arrays
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