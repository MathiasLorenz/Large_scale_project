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

void jacobi_cuda_1(Information *information, double *U, double *F, double *Unew)
{
	// Check if a device is available
	int Ndevices;
	cudaGetDeviceCount(&Ndevices);
	if ( Ndevices < 1)
	{
		fprintf(stderr, "Version cuda_1 must run on a CUDA device.\n");
		return;
	}
	cudaSetDevice(0);
	// Check for errors in the input
	if(!U || !F || !Unew) { fprintf(stderr,"Pointer is NULL.\n"); return; }

	// ========================================================================
	// Preparation

	// Firstly lets define the information structure on the device
	Information *information_cuda;
	checkCudaErrors(cudaMalloc( (void**)&information_cuda, sizeof(Information) ));
	copy_information_cuda(information_cuda,information);

	// Read the information structure
	int rank   = information->rank;
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
	double *U_cuda, *F_cuda, *Unew_cuda;

	checkCudaErrors(cudaMalloc((void**)&U_cuda,    arraySizes));
	checkCudaErrors(cudaMalloc((void**)&F_cuda,    arraySizes));
	checkCudaErrors(cudaMalloc((void**)&Unew_cuda, arraySizes));

	// Copy data to the GPU
	checkCudaErrors(cudaMemcpyAsync(U_cuda   , U   , arraySizes, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpyAsync(F_cuda   , F   , arraySizes, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpyAsync(Unew_cuda, Unew, arraySizes, cudaMemcpyHostToDevice));
	
	// ------------------------------------------------------------------------
	// Setup blocks for the GPU
	dim3 BlockSize = dim3(16,16,4);
	dim3 BlockAmount = dim3( K/BlockSize.x + 1, J/BlockSize.y + 1, I/BlockSize.z + 1 );

	// ------------------------------------------------------------------------
	// Prepare stop criterion
	int iter = 0;

	// ========================================================================
	// Run the iterative method
	checkCudaErrors(cudaDeviceSynchronize());
    for(iter = 0; iter < maxit ; iter++){
		// Compute the iteration of the jacobi method
		jacobi_iteration_kernel<<<BlockAmount,BlockSize>>>(information_cuda, U_cuda, F_cuda, Unew_cuda);
		
		// Wait for kernel to complete
		checkCudaErrors(cudaDeviceSynchronize());

		// Swap the arrays
		swap_array( &U_cuda, &Unew_cuda );
		if (information->use_tol)
			compute_relative_norm_cuda(information,information_cuda,U_cuda,Unew_cuda);
		// Stop early if relative error is used.
		// Second operand is only evaluated if the first is true
		if (information->use_tol && norm_early_stop(information))
			{iter++; break;}
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