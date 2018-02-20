#ifndef __CUDA_TESTS__
#define __CUDA_TESTS__

// Define all the functions which must run on the GPU
#ifdef __cplusplus
extern "C"
{	
	void cuda_test();
	__global__ void cu_test();
}
#else
	void cuda_test();
	void cu_test();
#endif


#endif
