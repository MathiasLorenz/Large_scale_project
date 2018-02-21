#ifndef __TESTS_H
#define __TESTS_H

// This function contains the prototypes for all tests that does not include
// any CUDA code.

void test_jacobi_2D(int Nx, int Ny, double tol, int maxiter);
void test_jacobi_3D(int Nx, int Ny, int Nz, double tol, int maxiter);

#endif // __TESTS_H
