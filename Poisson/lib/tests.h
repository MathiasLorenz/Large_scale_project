#ifndef __TESTS_H
#define __TESTS_H
// This function contains the prototypes for all tests that does not include
// any CUDA code.

#include "jacobi_util.h"

void test_jacobi_2D(int Nx, int Ny);
void test_jacobi_3D(int Nx, int Ny, int Nz);

void test_jacobi_mpi3D_1(int Nx, int Ny, int Nz);
void test_jacobi_mpi3D_2(Information *information);

#endif // __TESTS_H
