// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "poisson.h"
#include "matrix_routines.h"
#include "cuda_routines.h"
#include "jacobi_util_cuda.h"

extern double MFLOP;

// ============================================================================
// JACOBI 3D SOLVER USING CUDA WITH PEER ACCESS BETWEEN 2 GPUS

void jacobi_cuda_2(Information *information, double *U, double *F, double *Unew)
{
    // Check for errors in the input
    if(!information || !U || !F || !Unew)
    { fprintf(stderr,"Pointer is NULL.\n"); return; }

    //printf("I'm %s on line %d.\n", __FILE__, __LINE__);
    
    // Read the information structure
    int rank   = information->rank;
    int size   = information->size;
	int Nx 	   = information->global_Nx;
	int Ny 	   = information->global_Ny;
	int Nz 	   = information->global_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];
    int maxit  = information->maxit;

    // Check if we run with the correct settings
	int Ndevices;
	cudaGetDeviceCount(&Ndevices);
	if ( (Ndevices != 2) && (size != 1) )
	{
		fprintf(stderr, "Version cuda_2 must run on exactly 1 core and 2 GPUs.\n");
		return;
    }

    //printf("I'm %s on line %d.\n", __FILE__, __LINE__);

    // Rewritting to C style coordinates and determine split in z-dim
    int J, K;
    J = loc_Ny; K = loc_Nx;
    int loc_Nz_split[2], I[2];
    loc_Nz_split[0] = I[0] = floor(Nz/2) + 1;
    loc_Nz_split[1] = I[1] = floor(Nz/2) + 2;
    
    // Make copies of the information structure for the GPUs
    Information info_to_gpu[2];
    info_to_gpu[0] = info_to_gpu[1] = *information;

    // Insert the GPU splits into information
    info_to_gpu[0].loc_Nz[rank] = loc_Nz_split[0];
    info_to_gpu[1].loc_Nz[rank] = loc_Nz_split[1];

    // Size of array that has to be copied
    int buffer_size = loc_Nz*loc_Ny*sizeof(double);

    // Determine array sizes
    int array_size[2];
    array_size[0] = loc_Nz_split[0]*K*J*sizeof(double);
    array_size[1] = loc_Nz_split[1]*K*J*sizeof(double);

    // Pointer offset for U on the second gpu
    int ptr_offset_gpu1 = loc_Nz_split[0]*K*J;

    cudaSetDevice(0); // Just because

    // Enable peer access between GPUs
    cuda_enable_peer_access();

    //printf("I'm %s on line %d.\n", __FILE__, __LINE__);

    // Define the information structure on the device
    Information *information_cuda[2];

    // Arrays, one for each gpu
    double *U_cuda[2], *F_cuda[2], *Unew_cuda[2];

	// ========================================================================
    // Preparation. Do for both GPUs. Need pointers to arrays on both 
    for (int gpu = 0; gpu < 2; ++gpu)
    {
        cudaSetDevice(gpu);

        // Copy information structure to the current device
        cuda_malloc( (void**)&information_cuda[gpu], sizeof(Information) );
        copy_information_cuda(information_cuda[gpu], &info_to_gpu[gpu]);

        // Allocat on the current gpu
        cuda_malloc((void**)&U_cuda[gpu],    array_size[gpu]);
        cuda_malloc((void**)&F_cuda[gpu],    array_size[gpu]);
        cuda_malloc((void**)&Unew_cuda[gpu], array_size[gpu]);

        // Copy data to the GPU
        // rewrite with copy_to_device_async
        checkCudaErrors(cudaMemcpyAsync(U_cuda[gpu]   , U+gpu*ptr_offset_gpu1,
            array_size[gpu], cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpyAsync(F_cuda[gpu]   , F+gpu*ptr_offset_gpu1,
            array_size[gpu], cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpyAsync(Unew_cuda[gpu], Unew+gpu*ptr_offset_gpu1,
            array_size[gpu], cudaMemcpyHostToDevice));
    }

    //printf("I'm %s on line %d.\n", __FILE__, __LINE__);



    // Pointers to send between gpus
    double *U_ptr_send[2], *U_ptr_recv[2];
    // First rank needs the last updated index
    U_ptr_send[0] = &Unew_cuda[0][IND_3D(I[0] - 2, 0, 0, I[0], J, K)];
    U_ptr_recv[0] = &Unew_cuda[0][IND_3D(I[0] - 1, 0, 0, I[0], J, K)];

    // Last rank needs the first updated index
    U_ptr_send[1] = &Unew_cuda[1][IND_3D(1, 0, 0, I[1], J, K)];
    U_ptr_recv[1] = &Unew_cuda[1][IND_3D(0, 0, 0, I[1], J, K)];

    // ------------------------------------------------------------------------
    // Setup blocks for the GPU
    dim3 BlockSize = dim3(32,32,32);
    dim3 BlockAmount = dim3( K/BlockSize.x + 1, J/BlockSize.y + 1, I[1]/BlockSize.z + 1 );
    
	// ------------------------------------------------------------------------
	// Prepare stop criterion
	int iter = 0;

	// ========================================================================
	// Run the iterative method
	cuda_peer_device_sync();
    for(iter = 0; iter < maxit ; iter++){
        // Start computations
        for (int gpu = 0; gpu < 2; ++gpu)
        {
            cudaSetDevice(gpu);
            jacobi_iteration_kernel<<<BlockSize,BlockAmount>>>
                (information_cuda[gpu], U_cuda[gpu], F_cuda[gpu], Unew_cuda[gpu]);
        }

        // Wait for both kernels
        cuda_peer_device_sync();
        cudaSetDevice(0);

        // Send and recv boundary arrays
        // From 0 to 1
        checkCudaErrors(cudaMemcpyPeerAsync(U_ptr_recv[1], 1,
            U_ptr_send[0], 0, buffer_size ));

        // From 1 to 0
        checkCudaErrors(cudaMemcpyPeerAsync(U_ptr_recv[0], 0,
            U_ptr_send[1], 1, buffer_size ));

        // Wait for both kernels
        cuda_peer_device_sync();

        // Swap arrays
        swap_array(&U_cuda[0], &Unew_cuda[0]);
        swap_array(&U_cuda[1], &Unew_cuda[1]);
	}
	
	// Save number of iterations
    information->iter = iter;
    
    //printf("I'm %s on line %d.\n", __FILE__, __LINE__);

	// ========================================================================
	// Finalise
    // Copy data from the GPU
    cudaSetDevice(0);
    checkCudaErrors(cudaMemcpyAsync(U, U_cuda[0] , array_size[0], cudaMemcpyDeviceToHost));
    cudaSetDevice(1);
    checkCudaErrors(cudaMemcpyAsync(U+ptr_offset_gpu1, U_cuda[1] , array_size[1], cudaMemcpyDeviceToHost));
    cuda_peer_device_sync();
    cudaSetDevice(0);

    //printf("I'm %s on line %d.\n", __FILE__, __LINE__);
	
	// Flop Counts:
	// jacobi_iteration_cuda: (iter)
	// 		Simple: 	21
	//		Divisions: 	5
	MFLOP += 1e-6*(21.0 + 5.0*4.0 )*iter*Nx*Ny*Nz;

	// ------------------------------------------------------------------------
	// Free the arrays
	cudaFree(U_cuda[0]);
	cudaFree(F_cuda[0]);
	cudaFree(Unew_cuda[0]);
    free_information_cuda(information_cuda[0]);
    cudaFree(U_cuda[1]);
	cudaFree(F_cuda[1]);
	cudaFree(Unew_cuda[1]);
    free_information_cuda(information_cuda[1]);
    
    //printf("I'm %s on line %d.\n", __FILE__, __LINE__);
}
// END OF FILE
// ============================================================================