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

void jacobi_cuda_1(Information *information, int maxit,
	double threshold, double *U, double *F, double *Unew)
{
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
	int rank = information->rank;
	int Nx 	 = information->global_Nx;
	int Ny 	 = information->global_Ny;
	int Nz 	 = information->global_Nz;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];

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
	dim3 BlockSize = dim3(32,32,32);
	dim3 BlockAmount = dim3( I/BlockSize.x + 1, J/BlockSize.y + 1, K/BlockSize.z + 1 );

	//printf("BS.x: %3d, BS.y: %3d, BS.z: %3d\n",BlockSize.x,BlockSize.y,BlockSize.z);
	//printf("BA.x: %3d, BA.y: %3d, BA.z: %3d\n",BlockAmount.x,BlockAmount.y,BlockAmount.z);

	// ------------------------------------------------------------------------
	// Prepare stop criterion
	int iter = 0;

	// Remember to implement tolerance
	/*
    // Prepare stop criterion
    bool use_tol = false;
	double norm_diff = 10.0;
	
	if (strcmp("on",getenv("USE_TOLERANCE")) == 0) { use_tol = true; }
	*/

	// ========================================================================
	// Run the iterative method
	checkCudaErrors(cudaDeviceSynchronize());
    for(iter = 0; iter < maxit ; iter++){
		// Remember to implement tolerance
		/*
        norm_diff = 0.0;
		*/

		// Compute the iteration of the jacobi method
		jacobi_cuda_iteration<<<BlockSize,BlockAmount>>>(information_cuda, U_cuda, F_cuda, Unew_cuda);
		
		checkCudaErrors(cudaDeviceSynchronize());

        // Swap the arrays and check for convergence
        swap_array( &U_cuda, &Unew_cuda );


		// Compute the stop criterion
		// Remember to implement tolerance
		/*
		double u = U_cuda[IND_3D(i, j, k, I, J, K)];
		double unew = Unew_cuda[IND_3D(i, j, k, I, J, K)];
		norm_diff  += (u - unew)*(u - unew);
		*/
		// Remember to implement tolerance
		/*
        norm_diff = sqrt(norm_diff);

        if (use_tol && (norm_diff < threshold))
            break;
		*/
    }

	// ========================================================================
	// Finalise
	// Copy data from the GPU
	checkCudaErrors(cudaMemcpy(U   , U_cuda   , arraySizes, cudaMemcpyDeviceToHost));
	
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
	


	// ------------------------------------------------------------------------
	// Free the arrays
	cudaFree(U_cuda);
	cudaFree(F_cuda);
	cudaFree(Unew_cuda);
	free_information_cuda(information_cuda);
}
// END OF FILE
// ============================================================================