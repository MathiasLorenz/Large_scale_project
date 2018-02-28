// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_routines.h"

/*
#if defined(__MACH__) && defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif
*/

// ============================================================================
// ALLOCATION AND DEALLOCATION

// Allocate vector.
double * dmalloc_1d(int m)
{
    if (m <= 0) return NULL;
    double * arr = malloc(m*sizeof(*arr));
    if (arr)
        return arr;
    else
        return NULL;
}

// Allocate a linear indexed m x n matrix
double * dmalloc_2d_l(int m, int n)
{
    if (m <= 0 || n <= 0) return NULL;
    double *A = malloc(m * n * sizeof(double *));
    if (A == NULL) return NULL;
	A[0] = 0.0;
    for (int i = 1; i < m; i++)
        A[i] = A[0] + i * n;

    return A;
}

// Allocate a double-prec m x n matrix
double ** dmalloc_2d(int m, int n)
{
    if (m <= 0 || n <= 0) return NULL;
    double **A = malloc(m * sizeof(double *));
    if (A == NULL) return NULL;
    A[0] = malloc(m*n*sizeof(double));
    if (A[0] == NULL) {
        free(A);
        return NULL;
    }
    for (int i = 1; i < m; i++)
        A[i] = A[0] + i * n;

    return A;
}


void dfree_2d(double **A)
{
    free(A[0]);
    free(A);
}

// ============================================================================
// PRINTING

// Print 2d array
void array_print_2d(double *A, const int m, const int n,const char* fmt)
{
    if (!A) {
        fprintf(stderr, "Pointer is NULL\n");
        return;
    }
	fprintf(stdout, "\n  ");
    for(size_t i = 0; i < m; i++) {
        for(size_t j = 0; j < n; j++) {
            fprintf(stdout, fmt, A[IND_2D(i,j,m,n)]);
        }
        fprintf(stdout, "\n  ");
    }
    printf("\n");
}

// Print 2d matrix
void dmatrix_print_2d(double **A, const int m, const int n,const char* fmt)
{
    if (!A) {
        fprintf(stderr, "Pointer is NULL\n");
        return;
    }
	fprintf(stdout, "\n  ");
    for(size_t i = 0; i < m; i++) {
        for(size_t j = 0; j < n; j++) {
            fprintf(stdout, fmt, A[i][j]);
        }
        fprintf(stdout, "\n  ");
    }
    printf("\n");
}

// Print vector
void dvector_print(double *v, const int m)
{
	if (!v) {
        fprintf(stderr, "Pointer is NULL\n");
        return;
    }
	fprintf(stdout, "\n  ");
    for(size_t i = 0; i < m; i++)
        fprintf(stdout, "  %5.3f\n", v[i]);
    fprintf(stdout, "\n");
}

// ============================================================================
// MATRIX-MATRIX OPERATIONS

// Multiply matrices A and B into C.
void matmult_nat(int m, int n, int k, 
	double** A, double** B, double** C)
{

    if (!A || !B || !C) {
        fprintf(stderr, "Pointer is NULL\n");
        return;
    }

    // Initialize C matrix
    for (size_t i = 0; i < m*n; i++)
        C[0][i] *= 0;

    // Calculate matrix-matrix product
    for(size_t i = 0; i < m; i++) {
        for(size_t j = 0; j < n; j++) {
            // Every element, do vector-vector multiplication
            for(size_t kk = 0; kk < k; kk++) {
                C[i][j] += A[i][kk]*B[kk][j];
            }
        }
    }
}

// ============================================================================
// SWAPPING OPERATIONS
void swap_array(double ** A, double ** B)
{
    if(!A || !B) {fprintf(stderr,"Pointer is NULL.\n"); return;}
    double * temp;
    temp = *A;
    *A = *B;
    *B = temp;
}

void swap_matrix(double ***A, double ***B)
{
    if(!A || !B) {fprintf(stderr,"Pointer is NULL.\n"); return;}
    double **temp;
    temp = *A;
    *A = *B;
    *B = temp;
}

void swap_double(double *A, double *B)
{
    if(!A || !B) {fprintf(stderr,"Pointer is NULL.\n"); return;}
    double temp;
    temp = *A;
    *A = *B;
    *B = temp;
}

// ============================================================================
// NORMS AND MEASUREMENTS

double frobenius_difference(double **A, double **B, int N_total)
{
    if(!A || !B) {fprintf(stderr,"Pointer is NULL.\n"); return -1.0;}
    int i, j;
    double sum = 0.0;
    for(i = 0; i < N_total; i++)
        for(j = 0; j < N_total; j++)
            sum += (A[i][j]-B[i][j])*(A[i][j]-B[i][j]);

    sum = sqrt(sum);
    return sum;
}

double frobenius_norm(double **A, int N_total)
{
    if(!A) {fprintf(stderr,"Pointer is NULL.\n"); return -1.0;}
    int i, j;
    double sum = 0.0;
    for(i = 0; i < N_total; i++)
        for(j = 0; j < N_total; j++)
            sum += A[i][j]*A[i][j];

    sum = sqrt(sum);
    return sum;
}
