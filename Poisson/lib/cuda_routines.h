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

void copy_to_device_async(double *host, int N, double *device, void *stream);
void copy_from_device_async(double *host, int N, double *device, void *stream);

//void copy_from_device_void(void *host, int N_bytes, void *device);

// ============================================================================
// UTILITY FUNCTIONS

void cuda_synchronize();
void cuda_stream_synchronize(void *stream);

void cuda_create_stream(void **stream);
void cuda_destroy_stream(void  *stream);

void cuda_set_device(int rank);
void cuda_get_device_count(int *count);

// ============================================================================
// PEER-TO-PEER FUNCTIONS

void cuda_enable_peer_access();
void cuda_peer_device_sync();

#ifdef __cplusplus
}
#endif //__cplusplus
#endif