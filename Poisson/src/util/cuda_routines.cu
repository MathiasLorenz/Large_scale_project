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

void cuda_malloc(void** device_array, int N){
	checkCudaErrors(cudaMalloc( device_array ,N));
}
void cuda_free(double *device_array){
	checkCudaErrors(cudaFree(device_array));
}

// ============================================================================
// COPY DATA TO AND FROM DEVICE

void copy_to_device(double *host, int N, double *device){
	checkCudaErrors(cudaMemcpyAsync(device, host, N, cudaMemcpyHostToDevice));
}
void copy_from_device(double *host, int N, double *device){
	checkCudaErrors(cudaMemcpyAsync(host, device, N, cudaMemcpyDeviceToHost));
}

// ============================================================================
// UTILITY FUNCTIONS

void cuda_synchronize(){
	checkCudaErrors(cudaDeviceSynchronize());
}