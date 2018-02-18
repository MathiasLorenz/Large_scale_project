#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "tests.h"
#include "matrix_routines.h"
#include "poisson.h"
#include "init_data.h"

/*
void test_jacobi_sequential(int N,double tol,int maxiter){

    int N_total = N+2;              // Total grid points
    double h    = 2.0/(N + 1.0);    // stepsize

    double **U = dmalloc_2d(N_total,N_total);
    double **f = dmalloc_2d(N_total,N_total);
    double **R = dmalloc_2d(N_total,N_total);
    if(!U || !f || !R) {fprintf(stderr,"Error in malloc, Pointer is NULL.\n");
        return;}

    if (strcmp("sin",getenv("PROBLEM_NAME")) == 0)
        init_sin(U, f, R, N_total, h);
    if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
        init_rad(U, f, R, N_total, h);

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0)
        printf("Memory: %10.4f ", 3.0*N_total*N_total*8/1024);


    jacobi_sequential(N_total, maxiter, tol, U, f, R);

    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0)
        dmatrix_print_2d(U, N_total, N_total, "%10g ");

    dfree_2d(U);
    dfree_2d(f);
    dfree_2d(R);

}
*/
void test_jacobi(int N,double tol,int maxiter){

    int N_total = N+2;              // Total grid points
    double h    = 2.0/(N + 1.0);    // stepsize

    double **U = dmalloc_2d(N_total,N_total);
    double **f = dmalloc_2d(N_total,N_total);
    double **R = dmalloc_2d(N_total,N_total);
    if(!U || !f || !R) {fprintf(stderr,"Error in malloc, Pointer is NULL.\n");
        return;}

    if (strcmp("sin",getenv("PROBLEM_NAME")) == 0)
        init_sin(U, f, R, N_total, h);
    if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
        init_rad(U, f, R, N_total, h);

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0)
        printf("Memory: %10.4f ", 3.0*N_total*N_total*8/1024);

    jacobi_openmp(N_total, maxiter, tol, *U, *f, *R);

    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0)
        dmatrix_print_2d(U, N_total, N_total, "%10g ");

    dfree_2d(U);
    dfree_2d(f);
    dfree_2d(R);
}
/*
void test_gs_sequential(int N,double tol,int maxiter){
    int N_total = N+2;              // Total grid points
    double h    = 2.0/(N + 1.0);    // stepsize

    double **U = dmalloc_2d(N_total,N_total);
    double **f = dmalloc_2d(N_total,N_total);
    if(!U || !f) {fprintf(stderr,"Error in malloc, Pointer is NULL.\n");
        return;}

    if (strcmp("sin",getenv("PROBLEM_NAME")) == 0)
        init_sin_gs(U, f, N_total, h);
    if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
        init_rad_gs(U, f, N_total, h);

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0)
        printf("Memory: %10.4f ", 2.0*N_total*N_total*8/1024);

    gauss_seidel(N_total,maxiter,tol,U,f);

    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0)
        dmatrix_print_2d(U, N_total, N_total, "%10g ");

    dfree_2d(U);
    dfree_2d(f);
}
void test_jacobi_openmpv2(int N,double tol,int maxiter){

    int N_total = N+2;              // Total grid points
    double h    = 2.0/(N + 1.0);    // stepsize

    double **U = dmalloc_2d(N_total,N_total);
    double **f = dmalloc_2d(N_total,N_total);
    double **R = dmalloc_2d(N_total,N_total);
    if(!U || !f || !R) {fprintf(stderr,"Error in malloc, Pointer is NULL.\n");
        return;}

    if (strcmp("sin",getenv("PROBLEM_NAME")) == 0)
        init_sin(U, f, R, N_total, h);
    if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
        init_rad(U, f, R, N_total, h);

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0)
        printf("Memory: %10.4f ", 3.0*N_total*N_total*8/1024);

    jacobi_openmpv2(N_total, maxiter, tol, U, f, R);



    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0)
        dmatrix_print_2d(U, N_total, N_total, "%10g ");

    dfree_2d(U);
    dfree_2d(f);
    dfree_2d(R);
}
void test_jacobi_openmpft(int N,double tol,int maxiter){

    int N_total = N+2;              // Total grid points
    double h    = 2.0/(N + 1.0);    // stepsize

    double **U = dmalloc_2d(N_total,N_total);
    double **f = dmalloc_2d(N_total,N_total);
    double **R = dmalloc_2d(N_total,N_total);
    if(!U || !f || !R) {fprintf(stderr,"Error in malloc, Pointer is NULL.\n");
        return;}

    if (strcmp("sin",getenv("PROBLEM_NAME")) == 0)
        init_sin_omp(U, f, R, N_total, h);
    if (strcmp("rad",getenv("PROBLEM_NAME")) == 0)
        init_rad_omp(U, f, R, N_total, h);

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0)
        printf("Memory: %10.4f ", 3.0*N_total*N_total*8/1024);

    jacobi_openmp(N_total, maxiter, tol, U, f, R);



    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0)
        dmatrix_print_2d(U, N_total, N_total, "%10g ");

    dfree_2d(U);
    dfree_2d(f);
    dfree_2d(R);
}
*/
