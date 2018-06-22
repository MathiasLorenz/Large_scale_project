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
// JACOBI 3D SOLVER USING CUDA

void jacobi_cuda_2(Information *information, double *U, double *F, double *Unew)
{
    // Check for errors in the input
    if(!information || !U || !F || !Unew)
    { fprintf(stderr,"Pointer is NULL.\n"); return; }
    
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

    // Rewritting to C style coordinates
    int I, J, K;
    I = loc_Nz; J = loc_Ny; K = loc_Nx;

    // ------------------------------------------------------------------------
    // Now define the arrays which hold the problem on the device
    int arraySizes = I*J*K*sizeof(double);
    
    // Check if a device is available
	int Ndevices;
	cudaGetDeviceCount(&Ndevices);
	if ( (Ndevices != 2) && (size != 1) )
	{
		fprintf(stderr, "Version cuda_2 must run on exactly 1 core and 2 GPUs.\n");
		return;
    }

    cudaSetDevice(0); // Just because
    
    // Enable peer access between GPUs
    cuda_enable_peer_access();

    // Define the information structure on the device and allocate
    Information *information_cuda;
    checkCudaErrors(cudaMalloc( (void**)&information_cuda, sizeof(Information) ));

    // Arrays, one for each gpu
    double *U_cuda[2], *F_cuda[2], *Unew_cuda[2];

	// ========================================================================
    // Preparation. Do for both GPUs. Need pointers to arrays on both 
    for (int gpu = 0; gpu < 2; gpu++)
    {
        cudaSetDevice(gpu);
        
        // Copy information structure to the current device
        copy_information_cuda(information_cuda,information);

        // Allocat on the current gpu
        checkCudaErrors(cudaMalloc((void**)&U_cuda[gpu],    arraySizes));
        checkCudaErrors(cudaMalloc((void**)&F_cuda[gpu],    arraySizes));
        checkCudaErrors(cudaMalloc((void**)&Unew_cuda[gpu], arraySizes));

        // Copy data to the GPU
        checkCudaErrors(cudaMemcpyAsync(U_cuda[gpu]   , U   , arraySizes, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpyAsync(F_cuda[gpu]   , F   , arraySizes, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpyAsync(Unew_cuda[gpu], Unew, arraySizes, cudaMemcpyHostToDevice));
    }


    // ------------------------------------------------------------------------
    // Setup blocks for the GPU
    dim3 BlockSize = dim3(32,32,32);
    dim3 BlockAmount = dim3( K/BlockSize.x + 1, J/BlockSize.y + 1, I/BlockSize.z + 1 );
    

	// ------------------------------------------------------------------------
	// Prepare stop criterion
	int iter = 0;

	// ========================================================================
	// Run the iterative method
	checkCudaErrors(cudaDeviceSynchronize());
    for(iter = 0; iter < maxit ; iter++){
		// Compute the iteration of the jacobi method
		jacobi_iteration_kernel<<<BlockSize,BlockAmount>>>(information_cuda, U_cuda, F_cuda, Unew_cuda);
		
		// Wait for kernel to complete
		checkCudaErrors(cudaDeviceSynchronize());

        // Swap the arrays
        swap_array( &U_cuda, &Unew_cuda );
	}
	
	// Save number of iterations
	information->iter = iter;

	// ========================================================================
	// Finalise
	// Copy data from the GPU
	checkCudaErrors(cudaMemcpy(U   , U_cuda   , arraySizes, cudaMemcpyDeviceToHost));
	
	// Flop Counts:
	// jacobi_iteration_cuda: (iter)
	// 		Simple: 	21
	//		Divisions: 	5
	MFLOP += 1e-6*(21.0 + 5.0*4.0 )*iter*Nx*Ny*Nz;

	// ------------------------------------------------------------------------
	// Free the arrays
	cudaFree(U_cuda);
	cudaFree(F_cuda);
	cudaFree(Unew_cuda);
	free_information_cuda(information_cuda);
}
// END OF FILE
// ============================================================================