// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "matrix_routines.h"
#include "jacobi_util.h"
#include "poisson.h"

#include "cuda_routines.h"
#include "jacobi_util_cuda.h"

extern double MFLOP;

// ============================================================================
// JACOBI 3D SOLVER

void jacobi_mixed_3(Information *information, double *U, double *F, double *Unew)
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
	
	// Number of requests for MPI send/recv. 
	int num_req = (rank > 1 && rank < (size - 1)) ? 4 : 2;
	MPI_Request req[num_req];

	// ------------------------------------------------------------------------
	// Preparation
	// Send and receive buffers for MPI
	int loc_size = loc_Nx*loc_Ny*sizeof(double);
	double *s_buf1, *r_buf1, *s_buf2, *r_buf2;
	cuda_malloc_host((void**)&s_buf1, loc_size); // Buffer 1 to send with MPI
	cuda_malloc_host((void**)&r_buf1, loc_size); // Buffer 1 to receive with MPI
	cuda_malloc_host((void**)&s_buf2, loc_size); // Buffer 2 to send with MPI
	cuda_malloc_host((void**)&r_buf2, loc_size); // Buffer 2 to receive with MPI
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

	copy_to_device(U,   arraySizes,U_cuda   );
	copy_to_device(F,   arraySizes,F_cuda   );
	copy_to_device(Unew,arraySizes,Unew_cuda);
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
		
		// Compute the iteration of the jacobi method. First boundary
		jacobi_iteration_cuda_separate(
			information, information_cuda, U_cuda, F_cuda, Unew_cuda, "b");
		cuda_synchronize();

		// Extract boundaries
		double *U_ptr_s1, *U_ptr_s2;
		double *U_ptr_r1, *U_ptr_r2;
		if (rank == 0) {
			// First rank needs the last updated index
			U_ptr_s1 = &Unew_cuda[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];
			U_ptr_r1 = &Unew_cuda[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			copy_from_device(s_buf1, N_buffer*sizeof(double), U_ptr_s1);
		} else if (rank == (size - 1)){
			// Last rank needs the first updated index
			U_ptr_s1 = &Unew_cuda[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_r1 = &Unew_cuda[IND_3D(0, 0, 0, I, J, K)];

			copy_from_device(s_buf1, N_buffer*sizeof(double), U_ptr_s1);
		} else {
			// All other ranks needs the first and last updated index
			U_ptr_s1 = &Unew_cuda[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_s2 = &Unew_cuda[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];

			U_ptr_r1 = &Unew_cuda[IND_3D(0, 0, 0, I, J, K)];
			U_ptr_r2 = &Unew_cuda[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];

			copy_from_device(s_buf1, N_buffer*sizeof(double), U_ptr_s1);
			copy_from_device(s_buf2, N_buffer*sizeof(double), U_ptr_s2);
		}
		
		// Determine source and destination
		int neighbour_1, neighbour_2;
		compute_neighbors(information, &neighbour_1, &neighbour_2);

		cuda_synchronize();

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
		jacobi_iteration_cuda_separate(
			information, information_cuda, U_cuda, F_cuda, Unew_cuda, "i");

		// Wait for boundaries to be completed
		MPI_Waitall(num_req, req, MPI_STATUS_IGNORE);

		// Synchronize and copy back
		copy_to_device(r_buf1, N_buffer*sizeof(double), U_ptr_r1);
		if (rank > 0 && rank < (size - 1) )
			copy_to_device(r_buf2, N_buffer*sizeof(double), U_ptr_r2);
		
		cuda_synchronize();

		// Swap the arrays
        swap_array( &U_cuda, &Unew_cuda );
		
    }

	// Save number of iterations
	information->iter = iter;

	// ------------------------------------------------------------------------
	// Finalise
	MPI_Barrier(MPI_COMM_WORLD);

	// Copy back the result
	copy_from_device(U,   arraySizes,U_cuda   );
	cuda_synchronize();

	// Free the arrays
	cuda_host_free(s_buf1); cuda_host_free(s_buf2);
	cuda_host_free(r_buf1); cuda_host_free(r_buf2);
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