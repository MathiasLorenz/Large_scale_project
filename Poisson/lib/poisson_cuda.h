#ifndef __POISSON_CUDA_H
#define __POISSON_CUDA_H

#include "jacobi_util.h"

// Define the cuda methods that work on the GPU
void jacobi_cuda_1(Information *information, int maxit, double threshold, 
	double *U, double *F, double *Unew);

#endif
