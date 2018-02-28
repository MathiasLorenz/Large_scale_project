#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "cuda_routines.h"

__global__ void cu_test()
{
	int localThread = threadIdx.x;
	int block 	= blockIdx.x;
	int threadNoBlk = blockDim.x;
	int globalThread= localThread+block*threadNoBlk;
	int threadNoGlb = gridDim.x*threadNoBlk;
	if (globalThread == 100)
		{ int *a = (int*) 0x10000; *a = 0; }
	printf("This is a test\n I am thread %2d out of %2d in block %d.\n "
		"My global thread id is %3d out of %3d\n",
		localThread,threadNoBlk,block,globalThread,threadNoGlb);
	return;
}
