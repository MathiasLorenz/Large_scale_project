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
	// Solution to use temporary structures was found at:
	// https://stackoverflow.com/questions/31133522/simple-operation-on-structure-in-cuda-segmentation-fault
	// Consider to get all the memory pinned before copying
	cudaSetDevice(0);
	int size = information->size;

	// Allocate the temporary information structure
	Information Temp;

	// Simple Structure Elements
	Temp.size = information->size;
	Temp.rank = information->rank;
	Temp.global_Nx = information->global_Nx;
	Temp.global_Ny = information->global_Ny;
	Temp.global_Nz = information->global_Nz;

	// Allocate and copy the Arrays
	checkCudaErrors(cudaMalloc( (void**) &Temp.loc_Nx, size*sizeof(int)));
	checkCudaErrors(cudaMalloc( (void**) &Temp.loc_Ny, size*sizeof(int)));
	checkCudaErrors(cudaMalloc( (void**) &Temp.loc_Nz, size*sizeof(int)));

	checkCudaErrors(cudaMemcpyAsync(
		Temp.loc_Nx, 
		information->loc_Nx,
		information->size*sizeof(int), 
		cudaMemcpyHostToDevice
	));
	checkCudaErrors(cudaMemcpyAsync(
		Temp.loc_Ny, 
		information->loc_Ny,
		information->size*sizeof(int), 
		cudaMemcpyHostToDevice
	));
	checkCudaErrors(cudaMemcpyAsync(
		Temp.loc_Nz, 
		information->loc_Nz,
		information->size*sizeof(int), 
		cudaMemcpyHostToDevice
	));

	checkCudaErrors(cudaDeviceSynchronize());

	// Copy over the information structure
	checkCudaErrors(cudaMemcpy(
		information_cuda, 
		&Temp,
		sizeof(Information), 
		cudaMemcpyHostToDevice
	));
	
	checkCudaErrors(cudaDeviceSynchronize());
}

__global__ void free_information_arrays_cuda(Information *information_cuda)
{
	free(information_cuda->loc_Nx);
	free(information_cuda->loc_Ny);
	free(information_cuda->loc_Nz);
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

    int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	if ( 
		( (i > 0) && (j > 0) && (k > 0)) 
		&& ((i < (I-1)) && (j < (J-1)) && (k < (K-1) )) ) 
	{
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
		Unew_cuda[ijk] = f6 * (ui + uj + uk);
	}
}

