#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "jacobi_util.h"
#include "jacobi_util_cuda.h"
#include "matrix_routines.h"
#include "cuda_routines.h"

// ============================================================================
// FUNCTIONS TO HANDLE THE INFORMATION STRUCTURE ON THE DEVICE

void copy_information_cuda(Information *information_cuda, Information *information)
{
	// Solution to use temporary structures was found at:
	// https://stackoverflow.com/questions/31133522/simple-operation-on-structure-in-cuda-segmentation-fault
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
	
	Temp.maxit				= information->maxit;
	Temp.iter				= information->iter;
	Temp.tol				= information->tol;
	Temp.use_tol			= information->use_tol;
	Temp.local_frobenius	= information->local_frobenius;
	Temp.frobenius_error	= information->frobenius_error;

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

void free_information_cuda(Information *information_cuda){
	free_information_arrays_cuda<<<1,1>>>(information_cuda);
	cudaFree(information_cuda);
}

// As we these free's are launched in a kernel we de a normal free and
// NOT a cudaFree.
__global__ void free_information_arrays_cuda(Information *information_cuda)
{
	free(information_cuda->loc_Nx);
	free(information_cuda->loc_Ny);
	free(information_cuda->loc_Nz);
}

// ============================================================================
// CUDA VERSION OF THE ITERATIVE CORE

void jacobi_iteration_cuda(Information *information, Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda)
{
	int rank = information->rank;
	int K = information->loc_Nx[rank];
	int J = information->loc_Ny[rank];
	int I = information->loc_Nz[rank];

	// Blocks
	dim3 BlockSize = dim3(16,16,4);
	dim3 BlockAmount = dim3( K/BlockSize.x + 3, J/BlockSize.y + 3, I/BlockSize.z + 3 );

	jacobi_iteration_kernel<<<BlockAmount,BlockSize>>>
		(information_cuda, U_cuda, F_cuda, Unew_cuda);
}

__global__ void jacobi_iteration_kernel(Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda)
{
	// Determine where the thread is located
	int k = threadIdx.x + blockDim.x*blockIdx.x;
	int j = threadIdx.y + blockDim.y*blockIdx.y;
	int i = threadIdx.z + blockDim.z*blockIdx.z;

	// Read the needed data from the information structure
	int rank = information_cuda->rank;
	int Nx = information_cuda->global_Nx;
	int Ny = information_cuda->global_Ny;
	int Nz = information_cuda->global_Nz;
	int loc_Nx = information_cuda->loc_Nx[rank];
	int loc_Ny = information_cuda->loc_Ny[rank];
	int loc_Nz = information_cuda->loc_Nz[rank];

	// For relative error stopping
	information_cuda->local_frobenius = 0.0;

    int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	// Compute grid points
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

// ============================================================================
// CUDA VERSION OF THE ITTERATIVE CORE
// This version handles interior and boundary separately

// Consider changing the BlockSize and BlockAmount for boundary call
void jacobi_iteration_cuda_separate(Information *information, Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda, const char *ver)
{
	int rank = information->rank;
	int K = information->loc_Nx[rank];
	int J = information->loc_Ny[rank];
	int I = information->loc_Nz[rank];

	// interior or boundary
	if (strcmp(ver, "i") == 0)
	{
		dim3 BlockSize = dim3(16, 16, 4);
		dim3 BlockAmount = dim3( K/BlockSize.x + 3, J/BlockSize.y + 3, I/BlockSize.z + 3 );
		jacobi_iteration_kernel_interior<<<BlockAmount,BlockSize>>>
				(information_cuda, U_cuda, F_cuda, Unew_cuda);
	}
	if (strcmp(ver, "b") == 0)   // boundary
	{
		dim3 BlockSize = dim3(32, 16, 2);
		dim3 BlockAmount = dim3( K/BlockSize.x + 3, J/BlockSize.y + 3, 1 );
		jacobi_iteration_kernel_boundary<<<BlockAmount,BlockSize>>>
				(information_cuda, U_cuda, F_cuda, Unew_cuda);
	}
}
void jacobi_iteration_cuda_separate_stream(Information *information, Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda, const char *ver,void *stream)
{
	cudaStream_t *s;
	s = (cudaStream_t*)stream;
	int rank = information->rank;
	int K = information->loc_Nx[rank];
	int J = information->loc_Ny[rank];
	int I = information->loc_Nz[rank];

	// interior or boundary
	if (strcmp(ver, "i") == 0)
	{
		dim3 BlockSize = dim3(16, 16, 4);
		dim3 BlockAmount = dim3( K/BlockSize.x + 3, J/BlockSize.y + 3, I/BlockSize.z + 3 );
		jacobi_iteration_kernel_interior<<<BlockAmount,BlockSize,0,*s>>>
				(information_cuda, U_cuda, F_cuda, Unew_cuda);
	}
	if (strcmp(ver, "b") == 0)   // boundary
	{
		dim3 BlockSize = dim3(32, 16, 2);
		dim3 BlockAmount = dim3( K/BlockSize.x + 3, J/BlockSize.y + 3, 1 );
		jacobi_iteration_kernel_boundary<<<BlockAmount,BlockSize,0,*s>>>
				(information_cuda, U_cuda, F_cuda, Unew_cuda);
	}
}


// Kernel for interior points. Starts being used in mixed_3
__global__ void jacobi_iteration_kernel_interior(Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda)
{
	// Determine where the thread is located
	int k = threadIdx.x + blockDim.x*blockIdx.x;
	int j = threadIdx.y + blockDim.y*blockIdx.y;
	int i = threadIdx.z + blockDim.z*blockIdx.z;

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

	// Setting up steps
	double hi = 2.0/(Nz-1.0);
	double hj = 2.0/(Ny-1.0);
	double hk = 2.0/(Nx-1.0);
	double stepi = hi*hi;
	double stepj = hj*hj;
	double stepk = hk*hk;
	double f3 = 1.0/3.0;
	double f6 = 1.0/6.0;

	// Interior case
	// Consider optimizing this as kernels for i \in {0, 1, K-2, K-1} are not used
	if ( 
		( (i > 1) && (j > 0) && (k > 0)) 
		&& ((i < (I-2)) && (j < (J-1)) && (k < (K-1) )) ) 
	{
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

// Kernel for boundary points. Starts being used in mixed_3
__global__ void jacobi_iteration_kernel_boundary(Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda)
{
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

	// Setting up steps
	double hi = 2.0/(Nz-1.0);
	double hj = 2.0/(Ny-1.0);
	double hk = 2.0/(Nx-1.0);
	double stepi = hi*hi;
	double stepj = hj*hj;
	double stepk = hk*hk;
	double f3 = 1.0/3.0;
	double f6 = 1.0/6.0;

	
	// Determine where the thread is located
	int k = threadIdx.x + blockDim.x*blockIdx.x;
	int j = threadIdx.y + blockDim.y*blockIdx.y;
	int i = ( (threadIdx.z + blockDim.z*blockIdx.z) == 0) ? 1 : I - 2;

	// Boundary case
	if ( 
		( (j > 0) && (k > 0)) 
		&& ( (j < (J-1)) && (k < (K-1) )) ) 
	{
		// Compute new value
		// Do for first boundary (i = 1)
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

// ============================================================================
// Function to compute the frobenious norm 

void compute_relative_norm_cuda(
	Information *information, Information *information_cuda,
	double *U_cuda, double *Unew_cuda)
{
	int N = 1;
	frobenious_kernel<<<1,N>>>(information_cuda,U_cuda,Unew_cuda);

	copy_from_device(
		&information->local_frobenius, sizeof(double),
		&information_cuda->local_frobenius
	);
}
__global__ void frobenious_kernel(
	Information *information_cuda, double *U_cuda, double *Unew_cuda)
{
	double local_frobenius=0.0;
	int rank = information_cuda->rank;
	int loc_Nx = information_cuda->loc_Nx[rank];
	int loc_Ny = information_cuda->loc_Ny[rank];
	int loc_Nz = information_cuda->loc_Nz[rank];
	
	int I, J, K;
	I = loc_Nz; J = loc_Ny; K = loc_Nx;

	//local_frobenius[threadIdx.x] = 0.0;
	// Loop over all interior points
	for (int ijk = threadIdx.x; ijk < I*J*K; ijk += blockDim.x) {
		//int ijk = IND_3D(i, j, k, I, J, K);
		double uij    = U_cuda[ijk];
		double unewij = Unew_cuda[ijk];
		local_frobenius += (uij - unewij)*(uij - unewij);
	}
	information_cuda->local_frobenius += local_frobenius;
}