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
void test_jacobi_3D(Information *information, char const *solver);

#ifdef __cplusplus
} // Extern "C"
#endif //__cplusplus

#endif // __TESTS_H
