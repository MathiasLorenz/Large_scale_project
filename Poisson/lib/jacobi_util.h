#ifndef __JACOBI_UTIL__
#define __JACOBI_UTIL__


typedef struct Information {
	int rank;
	int size;

	int global_Nx;
	int global_Ny;
	int global_Nz;

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
