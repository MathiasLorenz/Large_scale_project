#ifndef __POISSON_CUDA_H
#define __POISSON_CUDA_H

// Define the cuda methods that work on the GPU

void jacobi_cuda_1(int Nx, int Ny, int Nz, int maxit, double threshold, double *U, double *F, double *Unew);

#endif
