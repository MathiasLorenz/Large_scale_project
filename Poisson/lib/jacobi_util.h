#ifndef __JACOBI_UTIL__
#define __JACOBI_UTIL__


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
void check_true_solution(double *A, double *U, double *norm_diff,
	Information *information);

#endif // __JACOBI_UTIL__
