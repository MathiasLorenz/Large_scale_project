#ifndef __MY_MATRIX_ROUTINES__
#define __MY_MATRIX_ROUTINES__

// Macros:
// Min
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Macro for accessing 2D array A with linear indexing
// This is called with first the Matrix then the array and finally the element.
#define IND_2D(i, j, I, J) ((i)*(J) + (j))
#define IND_3D(i, j, k, I, J, K) ((i)*(J)*(K) + (j)*(K) + (k))


// Allocation and Deallocation of matrices and vectors
double *  dmalloc_1d	(int m);
double *  dmalloc_2d_l	(int m, int n);
double ** dmalloc_2d	(int m, int n);
void dfree_2d	(double **A);

// Printing of matrices and vectors
void dvector_print(double *v, const int m);
void dmatrix_print_2d(double **A, const int m, const int n, const char* set);
void array_print_2d(double *A, const int m, const int n, const char* set);

// Matrix-Matrix opperations
void matmult_nat(int m, int n, int k, double **A, double **B, double **C);

// Swapping of arrays
void swap_double(double *A, double *B);
void swap_array (double **A, double **B);
void swap_matrix(double ***A, double ***B);

// Norms and other measures
double frobenius_difference(double **A, double **B, int N_total);
double frobenius_norm(double **A, int N_total);

#endif // __MY_MATRIX_ROUTINES__
