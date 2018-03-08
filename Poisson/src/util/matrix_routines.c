// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_routines.h"

// ============================================================================
// ALLOCATION AND DEALLOCATION

// Allocate vector.
double * dmalloc_1d(int Nx)
{
    if (Nx <= 0) return NULL;
    double * arr = malloc(Nx*sizeof(*arr));
    if (arr)
        return arr;
    else
        return NULL;
}

// Allocate a double-prec Nx x n matrix
double ** dmalloc_2d(int Nx, int Ny)
{
    if (Nx <= 0 || Ny <= 0) return NULL;
    double **A = malloc(Nx * sizeof(double *));
    if (A == NULL) return NULL;
    A[0] = malloc(Nx*Ny*sizeof(double));
    if (A[0] == NULL) {
        free(A);
        return NULL;
    }
    for (int i = 1; i < Nx; i++)
        A[i] = A[0] + i * Ny;

    return A;
}

// Allocate a linear indexed Nx x n matrix
double * dmalloc_2d_l(int Nx, int Ny)
{
	int S = Nx*Ny;
    if (S <= 0) return NULL;
    double * arr = malloc(S*sizeof(*arr));
    if (arr)
        return arr;
    else
        return NULL;
}

// Allocate a linear indexed Nx x n x Nz matrix
double * dmalloc_3d_l(int Nx, int Ny, int Nz)
{
	int S = Nx*Ny*Nz;
    if (S <= 0) return NULL;
    double * arr = malloc(S*sizeof(*arr));
    if (arr)
        return arr;
    else
        return NULL;
}

// Free 2 dimentional array
void free_2d(double **A)
{
    free(A[0]);
    free(A);
}

// ============================================================================
// PRINTING

// Print vector
void dvector_print(double *v, const int Nx)
{
	if (!v) {
        fprintf(stderr, "Pointer is NULL\n");
        return;
    }
	fprintf(stdout, "\n  ");
    for(size_t i = 0; i < Nx; i++)
        fprintf(stdout, "  %5.3f\n", v[i]);
    fprintf(stdout, "\n");
}

// Print 2d array
void array_print_2d(double *A, const int Nx, const int Ny, const char* fmt)
{
    int I = Ny;
    int J = Nx;
    if (!A) {
        fprintf(stderr, "Pointer is NULL\n");
        return;
    }
	fprintf(stdout, "\n  ");
    for(size_t i = 0; i < I; i++) {
        for(size_t j = 0; j < J; j++) {
            fprintf(stdout, fmt, A[IND_2D(i,j,I,J)]);
        }
        fprintf(stdout, "\n  ");
    }
    printf("\n");
}

// Print 2d matrix
void dmatrix_print_2d(double **A, const int Nx, const int Ny,const char* fmt)
{
	int I = Ny;
    int J = Nx;

    if (!A) {
        fprintf(stderr, "Pointer is NULL\n");
        return;
    }
	fprintf(stdout, "\n  ");
    for(size_t i = 0; i < I; i++) {
        for(size_t j = 0; j < J; j++) {
            fprintf(stdout, fmt, A[i][j]);
        }
        fprintf(stdout, "\n  ");
    }
    printf("\n");
}

// Print 3d array slice
void array_print_3d_slice(double *A, const int Nx, const int Ny, const int Nz, 
	const int slice, const char* fmt)
{
    if (!A) { fprintf(stderr, "Pointer is NULL\n"); return; }
	if (slice > Nz) 
	{ 
		fprintf(stderr, "Slice is not contained in the data\n"); 
		return; 
	}
    int I = Nz;
    int J = Ny;
	int K = Nx;

	int i = slice;
	fprintf(stdout, "\n  ");
    for(size_t j = 0; j < J; j++) {
        for(size_t k = 0; k < K; k++) {
            fprintf(stdout, fmt, A[IND_3D(i,j,k,I,J,K)]);
        }
        fprintf(stdout, "\n  ");
    }
    printf("\n");
}

// Print 3d array 
void array_print_3d(double *A, const int Nx, const int Ny, const int Nz,
	const char* fmt)
{
    if (!A) { fprintf(stderr, "Pointer is NULL\n"); return; }
    int I = Nz;
    int J = Ny;
	int K = Nx;

	fprintf(stdout, "%d %d %d\n",Nx,Ny,Nz);
    for(size_t i = 0; i < I; i++)
        for(size_t j = 0; j < J; j++)
			for(size_t k = 0; k < K; k++)
            	fprintf(stdout, fmt, A[IND_3D(i,j,k,I,J,K)]);
	
    printf("\n");
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
