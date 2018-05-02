#ifndef __TESTS_CUDA__
#define __TESTS_CUDA__

#include "jacobi_util.h"

// Makes sure that linker can understand the code
#ifdef __cplusplus
extern "C"{
#endif //__cplusplus

// Define all test functions which include CUDA code.
void test_cuda_1(Information *information);


#ifdef __cplusplus
} // Extern "C"
#endif //__cplusplus

#endif //__TESTS_CUDA__
