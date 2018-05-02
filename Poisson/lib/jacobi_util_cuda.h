#ifndef __JACOBI_UTIL_CUDA__
#define __JACOBI_UTIL_CUDA__

__global__ void jacobi_cuda_iteration(int I, int J, int K,
	double *U, double *F, double *Unew);


#endif // __JACOBI_UTIL_CUDA__