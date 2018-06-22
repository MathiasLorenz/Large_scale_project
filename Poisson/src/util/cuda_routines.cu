#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "matrix_routines.h"
#include "cuda_routines.h"
#include "jacobi_util.h"

// ============================================================================
// ALLOCATE AND DEALLOCATE DATA ON THE DEVICE

void cuda_malloc(void **device_array, int N_bytes)
{
	checkCudaErrors(cudaMalloc(device_array , N_bytes));
}
void cuda_malloc_host(void **host_array, int N_bytes)
{
	checkCudaErrors(cudaMallocHost(host_array, N_bytes));
}

void cuda_free(double *device_array)
{
	checkCudaErrors(cudaFree(device_array));
}
void cuda_host_free(double *host_array)
{
	checkCudaErrors(cudaFreeHost(host_array));
}

// ============================================================================
// COPY DATA TO AND FROM DEVICE

void copy_to_device_async(double *host, int N_bytes, double *device)
{
	checkCudaErrors(cudaMemcpyAsync(device, host, N_bytes, cudaMemcpyHostToDevice));
}
void copy_from_device_async(double *host, int N_bytes, double *device)
{
	checkCudaErrors(cudaMemcpyAsync(host, device, N_bytes, cudaMemcpyDeviceToHost));
}
void copy_to_device(double *host, int N_bytes, double *device)
{
	checkCudaErrors(cudaMemcpy(device, host, N_bytes, cudaMemcpyHostToDevice));
}
void copy_from_device(double *host, int N_bytes, double *device)
{
	checkCudaErrors(cudaMemcpy(host, device, N_bytes, cudaMemcpyDeviceToHost));
}
/*
void copy_from_device_void(void *host, int N_bytes, void *device)
{
	checkCudaErrors(cudaMemcpy(host, device, N_bytes, cudaMemcpyDeviceToHost));
}
*/

// ============================================================================
// UTILITY FUNCTIONS

void cuda_synchronize(){
	checkCudaErrors(cudaDeviceSynchronize());
}
void cuda_get_device_count(int *count)
{
	cudaGetDeviceCount( count );
}
void cuda_set_device(int rank)
{
	int Ndevices, device;
	cudaGetDeviceCount( &Ndevices );
	if (Ndevices > 1)
		device = (rank % 2 == 0) ? 0 : 1;
	else 
		device = 0;
	cudaSetDevice(device);
}


// Enable peer access
void cuda_enable_peer_access(const int access_from, const int access_to)
{
	cudaSetDevice(access_from);
	checkCudaErrors(cudaDeviceEnablePeerAccess(access_to, 0));
}