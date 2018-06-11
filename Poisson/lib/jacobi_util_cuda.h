#ifndef __JACOBI_UTIL_CUDA__
#define __JACOBI_UTIL_CUDA__

#include "jacobi_util.h"
#ifdef __cplusplus
extern "C"{
#endif //__cplusplus

void copy_information_cuda(Information *information_cuda, Information *information);
void free_information_cuda(Information *information_cuda);
void jacobi_iteration_cuda(
	Information *information, 
	Information *information_cuda, 
	double *U_cuda, double *F_cuda, double *Unew_cuda
);
void jacobi_iteration_cuda_separate(
	Information *information,
	Information *information_cuda,
	double *U_cuda, double *F_cuda, double *Unew_cuda,
	const char *ver
);
void cuda_wait_boundary();

#ifdef __cplusplus
__global__ void free_information_arrays_cuda(Information *information_cuda);
__global__ void jacobi_iteration_kernel(
	Information *information_cuda, 
	double *U_cuda, double *F_cuda, double *Unew_cuda
);
__global__ void jacobi_iteration_kernel_tol(
	Information *information_cuda, 
	double *U_cuda, double *F_cuda, double *Unew_cuda
);
__global__ void jacobi_iteration_kernel_interior(
	Information *information_cuda, 
	double *U_cuda, double *F_cuda, double *Unew_cuda
);
__global__ void jacobi_iteration_kernel_boundary(
	Information *information_cuda, 
	double *U_cuda, double *F_cuda, double *Unew_cuda
);
#endif //__cplusplus


#ifdef __cplusplus
}
#endif //__cplusplus

#endif // __JACOBI_UTIL_CUDA__