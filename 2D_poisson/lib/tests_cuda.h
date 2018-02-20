#ifndef __TESTS_CUDA__
#define __TESTS_CUDA__

// Makes sure that linker can understand the code
#ifdef __cplusplus
extern "C"{
#endif

// Define all test functions which include CUDA code.
void test_cuda(int N, double tol, int maxiter);

#ifdef __cplusplus
} // Extern "C"
#endif

#endif
