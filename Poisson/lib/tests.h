#ifndef __TESTS_H
#define __TESTS_H

// This function contains the prototypes for all tests that does not include
// any CUDA code.

void test_jacobi_2D(int Nx, int Ny);
void test_jacobi_3D(int Nx, int Ny, int Nz);

#endif // __TESTS_H