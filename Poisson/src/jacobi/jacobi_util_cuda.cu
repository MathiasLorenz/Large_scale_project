#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "matrix_routines.h"
#include "cuda_routines.h"

// ============================================================================
// FUNCTIONS TO HANDLE THE INFORMATION STRUCTURE ON THE DEVICE

void copy_information_cuda(Information *information_cuda, Information *information)
{
	// Consider to get all the memory pinned before copying

	int size = information->size;

	cudaMemcpyAsync( &information_cuda->size, &information->size,
		sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpyAsync( &information_cuda->rank, &information->rank,
		sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpyAsync( &information_cuda->global_Nx, &information->global_Nx,
		sizeof(int),cudaMemcpyHostToDevice );
	cudaMemcpyAsync( &information_cuda->global_Ny, &information->global_Ny,
		sizeof(int),cudaMemcpyHostToDevice );
	cudaMemcpyAsync( &information_cuda->global_Nz, &information->global_Nz,
		sizeof(int),cudaMemcpyHostToDevice );

	cudaDeviceSynchronize();

	printf("Starting the Malloc procedures\n");

	// So the current errors are placed in these three lines. We need to find out how to 
	// allocate these arrays correctly.
	checkCudaErrors(cudaMalloc( (void**) &information_cuda->loc_Nx, size*sizeof(int)));
	checkCudaErrors(cudaMalloc( (void**) &information_cuda->loc_Ny, size*sizeof(int)));
	checkCudaErrors(cudaMalloc( (void**) &information_cuda->loc_Nz, size*sizeof(int)));

	printf("Completed malloc, moving on to memcpy\n");

	checkCudaErrors(cudaMemcpyAsync(&information_cuda->loc_Nx, &information->loc_Nx,
		information->size*sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpyAsync(&information_cuda->loc_Ny, &information->loc_Ny,
		information->size*sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpyAsync(&information_cuda->loc_Nz, &information->loc_Nz,
		information->size*sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaDeviceSynchronize());
}

void free_information_arrays_cuda(Information *information_cuda)
{
	cudaFree(information_cuda->loc_Nx);
	cudaFree(information_cuda->loc_Ny);
	cudaFree(information_cuda->loc_Nz);
}

// ============================================================================
// CUDA VERSION OF THE ITTERATIVE CORE

__global__ void jacobi_cuda_iteration(Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda)
{

	// Determine where the thread is located
	int j = threadIdx.x + blockDim.x*blockIdx.x;
	int i = threadIdx.y + blockDim.y*blockIdx.y;
	int k = threadIdx.z + blockDim.z*blockIdx.z;

	// Read the needed data from the information structure
	int rank = information_cuda->rank;
	int Nx = information_cuda->global_Nx;
	int Ny = information_cuda->global_Ny;
	int Nz = information_cuda->global_Nz;
	int loc_Nx = information_cuda->loc_Nx[rank];
	int loc_Ny = information_cuda->loc_Ny[rank];
	int loc_Nz = information_cuda->loc_Nz[rank];

	printf("The rank is: %d, Hello i am, i: %d, j: %d, k: %d\n",rank,i,j,k);
    int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	// Setting up steps
	double hi = 2.0/(Nz-1.0);
	double hj = 2.0/(Ny-1.0);
	double hk = 2.0/(Nx-1.0);
	double stepi = hi*hi;
	double stepj = hj*hj;
	double stepk = hk*hk;
	double f3 = 1.0/3.0;
	double f6 = 1.0/6.0;

	// Compute new value
	if ( (i > 0 && j > 0 && k > 0) && (i < (I-1) && j < (J-1) && k < (K-1) ) )
	{
		// Save i, j, k index once
		int ijk = IND_3D(i, j, k, I, J, K);

		// Linear indexing with macro
		double ui = U_cuda[IND_3D(i - 1, j, k, I, J, K)] 
			+ U_cuda[IND_3D(i + 1, j, k, I, J, K)] 
			+ f3 * stepi * F_cuda[ijk];
		double uj = U_cuda[IND_3D(i, j - 1, k, I, J, K)] 
			+ U_cuda[IND_3D(i, j + 1, k, I, J, K)] 
			+ f3 * stepj * F_cuda[ijk];
		double uk = U_cuda[IND_3D(i, j, k - 1, I, J, K)] 
			+ U_cuda[IND_3D(i, j, k + 1, I, J, K)] 
			+ f3 * stepk * F_cuda[ijk];

		// Collect terms
		Unew_cuda[ijk] = f6 * (ui + uj + uk) + 10;
		printf("%f\n",f6 * (ui + uj + uk) + 10);
	}
}

