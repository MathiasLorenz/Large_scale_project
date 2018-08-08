#ifndef __JACOBI_UTIL__
#define __JACOBI_UTIL__

#include <stdbool.h> // For bool

// Makes sure that linker can understand the code
#ifdef __cplusplus
extern "C"{
#endif //__cplusplus
// ============================================================================
// Handle the information structure
typedef struct Information {
	// MPI information of size and local rank.
	int rank;
	int size;

	// Global dimensions of the problem
	int global_Nx;
	int global_Ny;
	int global_Nz;

	// Arrays holding the local number of gridpoints for each rank.
	int *loc_Nx;	
	int *loc_Ny;
	int *loc_Nz;

	// Variables for looping
	int 	maxit;
	int		iter;

	// For relative stop criterion
	bool 	use_tol;
	double  tol;
	double  local_frobenius;
	double  frobenius_error;

} Information;

void write_information(Information *information, int Nx, int Ny, int Nz, int rank, int size);
void free_information_arrays(Information *information);

// Function to determine neighbour ranks
void compute_neighbors(
	Information *information, int *neighbour_1, int *neighbour_2
);
// ============================================================================
// Functions for computing jacobi iterations
void jacobi_iteration(
	Information *information,double *U,double *F,double *Unew
);
void jacobi_iteration_separate(
	Information *information,double *U,double *F,double *Unew,const char *ver
);

// ============================================================================
// Functions to handle maximal absolute error
void print_error(Information *information, double *U);

void compute_local_max_error(
	Information *information, double *U, double *local_error);
void compute_global_max_error(
	Information *information, double *U, double *global_error);

// ============================================================================
// Function to handle frobenius stop criterion
bool norm_early_stop(Information *information);

#ifdef __cplusplus
} // Extern "C"
#endif //__cplusplus

#endif // __JACOBI_UTIL__
