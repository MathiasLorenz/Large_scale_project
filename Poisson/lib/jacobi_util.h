#ifndef __JACOBI_UTIL__
#define __JACOBI_UTIL__

// Makes sure that linker can understand the code
#ifdef __cplusplus
extern "C"{
#endif //__cplusplus

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
} Information;

void write_information(Information *information, int Nx, int Ny, int Nz, int rank, int size);
void free_information_arrays(Information *information);
void jacobi_iteration(int I, int J, int K, int rank, int global_Nz,
	double *U, double *F, double *Unew);
void generate_true_solution(double *A, Information *information);
<<<<<<< HEAD
void check_true_solution(double *A, double *U, double *norm_diff,
	Information *information);

#ifdef __cplusplus
} // Extern "C"
#endif //__cplusplus

=======
void compute_max_error(Information *information, double *A, double *U, double *local_error);
void compute_global_error(Information *information, double *A, double *U);
>>>>>>> 3ffed517998227385e269cf68a693722327445e8
#endif // __JACOBI_UTIL__
