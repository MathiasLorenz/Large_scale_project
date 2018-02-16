#ifndef __MY_MATRIX_ROUTINES__
#define __MY_MATRIX_ROUTINES__

// Macros:
// Min
#define MIN(a, b) ((a) < (b) ? (a) : (b))


double ** dmalloc_2d(int m, int n);
void dfree_2d(double **A);
double * dmalloc_1d(int m);
void dvector_print(double *v, const int m);
void dmatrix_print_2d(double **A, const int m, const int n, const char* set);
void matmult_nat(int m, int n, int k, double **A, double **B, double **C);
void swap_matrix(double ***A, double ***B);
void swap_double(double *A, double *B);


double frobenius_difference(double **A, double **B, int N_total);
double frobenius_norm(double **A, int N_total);

#endif // __MY_MATRIX_ROUTINES__
