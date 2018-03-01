#ifndef __POISSON_H
#define __POISSON_H

void jacobi_openmp_2D(int Nx, int Ny, int maxit, double threshold,
        double *u, double *f, double *tmp);
void jacobi_openmp_3D(int Nx, int Ny, int Nz, int maxit, double threshold,
        double *u, double *f, double *tmp);

#endif // __POISSON_H
