#ifndef __POISSON_H
#define __POISSON_H

// Makes sure that linker can understand the code
#ifdef __cplusplus
extern "C"{
#endif //__cplusplus

void jacobi_openmp_2D(int Nx, int Ny, int maxit, double threshold,
        double *u, double *f, double *tmp);
void jacobi_openmp_3D(int Nx, int Ny, int Nz, int maxit, double threshold,
        double *u, double *f, double *tmp);
void jacobi_mpi3D_1(int loc_Nx, int loc_Ny, int loc_Nz, int maxit, double threshold, int rank, int global_Nz,
    double *U, double *F, double *Unew);
void jacobi_mpi3D_2(Information *information, int maxit, double threshold,
    double *U, double *F, double *Unew);


#ifdef __cplusplus
} // Extern "C"
#endif //__cplusplus

#endif // __POISSON_H
