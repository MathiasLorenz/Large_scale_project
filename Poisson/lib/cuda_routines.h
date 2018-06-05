#ifndef __CUDA_ROUTINES__
#define __CUDA_ROUTINES__

#include "jacobi_util.h"
#ifdef __cplusplus
extern "C"{
#endif //__cplusplus

// ============================================================================
// ALLOCATE AND DEALLOCATE DATA ON THE DEVICE

void cuda_malloc(void **device_array, int N);
void cuda_malloc_host(void **host_array, int N);
void cuda_free(double *device_array);
void cuda_host_free(double *host_array);

// ============================================================================
// COPY DATA TO AND FROM DEVICE

void copy_to_device(double *host, int N, double *device);
void copy_from_device(double *host, int N, double *device);

// ============================================================================
// UTILITY FUNCTIONS

void cuda_synchronize();
void setCudaDevice(int rank);

#ifdef __cplusplus
}
#endif //__cplusplus
#endif