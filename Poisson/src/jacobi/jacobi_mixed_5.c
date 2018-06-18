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

void jacobi_mixed_5(Information *information, double *U, double *F, double *Unew)
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
		fprintf(stderr, "Version mixed_5 must run on a CUDA device.\n");
		return;
	}
	if ( size < 2)
	{
		fprintf(stderr, "Version mixed_5 must run on multiple mpi ranks.\n");
		return;
	}

	// Set the CUDA device
	cuda_set_device(rank);
	
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

	copy_to_device_async(U,   arraySizes,U_cuda   );
	copy_to_device_async(F,   arraySizes,F_cuda   );
	copy_to_device_async(Unew,arraySizes,Unew_cuda);
	cuda_synchronize();
	
	// Prepare stop criterion
	int iter = 0;

	// ------------------------------------------------------------------------
	// Run the iterative method
    for(iter = 0; iter < maxit ; iter++){
		
		// Compute the iteration of the jacobi method. First boundary
		jacobi_iteration_cuda_separate(
			information, information_cuda, U_cuda, F_cuda, Unew_cuda, "b");
		
		// Compute interior while boundary is being sent
		jacobi_iteration_cuda_separate(
			information, information_cuda, U_cuda, F_cuda, Unew_cuda, "i");

		cuda_wait_boundary();
		// Extract boundaries
		double *U_ptr_s1, *U_ptr_s2;
		double *U_ptr_r1, *U_ptr_r2;
		if (rank == 0) {
			// First rank needs the last updated index
			U_ptr_s1 = &Unew_cuda[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];
			U_ptr_r1 = &Unew_cuda[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];
		} else if (rank == (size - 1)){
			// Last rank needs the first updated index
			U_ptr_s1 = &Unew_cuda[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_r1 = &Unew_cuda[IND_3D(0, 0, 0, I, J, K)];
		} else {
			// All other ranks needs the first and last updated index
			U_ptr_s1 = &Unew_cuda[IND_3D(1, 0, 0, I, J, K)];
			U_ptr_s2 = &Unew_cuda[IND_3D(loc_Nz - 2, 0, 0, I, J, K)];

			U_ptr_r1 = &Unew_cuda[IND_3D(0, 0, 0, I, J, K)];
			U_ptr_r2 = &Unew_cuda[IND_3D(loc_Nz - 1, 0, 0, I, J, K)];
		}
		
		// Determine source and destination
		int neighbour_1, neighbour_2;
		compute_neighbors(information, &neighbour_1, &neighbour_2);

		printf("Okay, Ready to send\n");
		printf("%f\n",U_ptr_s1[0]);
		printf("ITS A POINTER\n");
		// Send boundaries and receive boundaries
		MPI_Isend(U_ptr_s1, N_buffer, MPI_DOUBLE, neighbour_1, 0, MPI_COMM_WORLD,
				&req[0]);
		printf("Okay, Ready to send time 2\n");
		if ( rank != 0 && rank != (size - 1) )
			MPI_Isend(U_ptr_s2, N_buffer, MPI_DOUBLE, neighbour_2, 0, MPI_COMM_WORLD,
					&req[2]);

		printf("Okay, Ready to receive\n");
		MPI_Irecv(U_ptr_r1, N_buffer, MPI_DOUBLE, neighbour_1, 0,
					MPI_COMM_WORLD, &req[1]);
		if ( rank != 0 && rank != (size - 1) )
			MPI_Irecv(U_ptr_r2, N_buffer, MPI_DOUBLE, neighbour_2, 0,
					MPI_COMM_WORLD, &req[3]);
		printf("Okay, Receive completed\n");
		// Wait for boundaries to be completed
		MPI_Waitall(num_req, req, MPI_STATUS_IGNORE);

		// Synchronize and copy back
		cuda_synchronize();

		// Swap the arrays
        swap_array( &U_cuda, &Unew_cuda );
		
    }

	// ------------------------------------------------------------------------
	// Finalise
	MPI_Barrier(MPI_COMM_WORLD);

	// Copy back the result
	copy_from_device_async(U,   arraySizes,U_cuda   );
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
}
// END OF FILE
// ============================================================================