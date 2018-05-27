#ifndef __TESTS_H
#define __TESTS_H
// This function contains the prototypes for all tests that does not include
// any CUDA code.

#include "jacobi_util.h"

// Makes sure that linker can understand the code
#ifdef __cplusplus
extern "C"{
#endif //__cplusplus

void test_jacobi_2D(int Nx, int Ny);
void test_jacobi_3D(int Nx, int Ny, int Nz);

void test_jacobi_mpi3D_1(Information *information);
void test_jacobi_mpi3D_2(Information *information);

void test_jacobi_mixed_1(Information *information);
void test_jacobi_mixed_2(Information *information);


#ifdef __cplusplus
} // Extern "C"
#endif //__cplusplus

#endif // __TESTS_H
