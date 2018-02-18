#ifndef __POISSON_H
#define __POISSON_H

//void jacobi_sequential(int N_total, int maxit, double threshold,
//        double **u, double **f, double **tmp);
void jacobi_openmp(int N_total, int maxit, double threshold,
        double *u, double *f, double *tmp);

/*
void gauss_seidel(int N_total, int maxit, double threshold,
        double **u, double **f);
void jacobi_openmpv2(int N_total, int maxit, double threshold,
        double **u, double **f, double **tmp);
void jacobi_vectorize(int N_total, int maxit, double threshold,
        double **u, double **f, double **tmp);
*/

#endif // __POISSON_H
