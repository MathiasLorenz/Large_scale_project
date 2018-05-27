// ============================================================================
// INCLUDES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "jacobi_util.h"
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
    double *arr = malloc(S*sizeof(double));
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

// Print 3d array sliced in the z-axis with MPI
// All ranks send their data to rank 0, which prints to stdout
void print_jacobi3d_z_sliced(const double *U, Information *information, const char* fmt)
{
	int rank = information->rank;
	int size = information->size;
	int loc_Nx = information->loc_Nx[rank];
	int loc_Ny = information->loc_Ny[rank];
	int loc_Nz = information->loc_Nz[rank];
	int global_Nz = information->global_Nz;

    // Initialise
    if (!U) { fprintf(stderr, "Pointer is NULL\n"); return; }
    int I = loc_Nz;
    int J = loc_Ny;
	int K = loc_Nx;

    // Get MPI communicator size
	MPI_Request req;

	// Make ready for recieving the data
	int N_buffer_r = loc_Nx*loc_Ny*max_array_int(information->loc_Nz, size);
	double *r_buf = malloc( N_buffer_r*sizeof(double) );
	if (!r_buf) { fprintf(stderr, "Pointer is NULL\n"); return; }

	// Handle transfer of data from processors with rank different from 0
    if ( rank > 0 )
    {
		// Number of elements to send
		int N_buffer = loc_Nx*loc_Ny*loc_Nz;

        // Send to rank 0
        MPI_Isend(U, N_buffer, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&req);

    } else // If we are rank 0, collect all data and print it
	{ 
		// Print the global dimensions.
		fprintf(stdout, "%d %d %d\n",loc_Nx,loc_Ny,global_Nz);

        // Print rank 0 without the last slice.
        for(size_t i = 0; i < (I-1); i++) {
            for(size_t j = 0; j < J; j++)
                for(size_t k = 0; k < K; k++)
                    fprintf(stdout, fmt, U[IND_3D(i,j,k,I,J,K)]); }


        // For rank 1..end-1 receive the data and print it
        for (int r = 1; r < (size - 1); r++)
        {
            // Read the local dimensions for rank r.
			int loc_Nz_r = information->loc_Nz[r];

            // Get data with blocking receive
            MPI_Recv(r_buf, N_buffer_r, MPI_DOUBLE, r, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Print from slice 1 to the second last
            for(size_t i = 1; i < (loc_Nz_r-1); i++) {
                for(size_t j = 0; j < J; j++)
                    for(size_t k = 0; k < K; k++)
                        fprintf(stdout, fmt, r_buf[IND_3D(i,j,k,loc_Nz_r,J,K)]);}

        }

		// Get the local size of last rank.
		int loc_Nz_r = information->loc_Nz[size-1];

		// Recieve the data.
        MPI_Recv(r_buf, N_buffer_r, MPI_DOUBLE, size-1, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Print from slice 1 to the last.
        for(size_t i = 1; i < loc_Nz_r; i++)
		{
            for(size_t j = 0; j < J; j++)
                for(size_t k = 0; k < K; k++)
                    fprintf(stdout, fmt, r_buf[IND_3D(i,j,k,loc_Nz_r+1,J,K)]);
		}
		printf("\n");
    }

	// Make sure all processors are synchronized before exiting.
    MPI_Barrier(MPI_COMM_WORLD);

	// Free the buffer.
	free(r_buf);
    return;
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
// MAX AND MIN FUNCTIONS

int max_array_int(int *array, int size)
{
	int old = array[0];
	for (int i = 1; i < size; i++)
		if (old < array[i]) old = array[i];
	return old;
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
