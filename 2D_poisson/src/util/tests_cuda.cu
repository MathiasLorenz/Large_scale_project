// ============================================================================
// INCLUDES
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "tests_cuda.h"
#include "cuda_routines.h"

// ============================================================================
// CUDA TEST
void test_cuda(int Nx, int Ny, int Nz)
{
	printf("We are now testing the cuda function\n");
	cu_test<<<1,1>>>();
	checkCudaErrors(cudaDeviceSynchronize());
}

