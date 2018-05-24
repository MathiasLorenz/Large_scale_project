#ifndef __JACOBI_UTIL_CUDA__
#define __JACOBI_UTIL_CUDA__

#include "jacobi_util.h"
void copy_information_cuda(Information *information_cuda, Information *information);
__global__ void free_information_arrays_cuda(Information *information_cuda);
__global__ void 
	jacobi_cuda_iteration(Information *information,double *U_cuda, double *F_cuda, double *Unew_cuda);


#endif // __JACOBI_UTIL_CUDA__