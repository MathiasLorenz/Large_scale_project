#ifndef __INIT_DATA_H
#define __INIT_DATA_H

#include "jacobi_util.h"

void init_sin_2D(double *U, double *F, double *Unew, int Nx, int Ny);
void init_rad_2D(double *U, double *F, double *Unew, int Nx, int Ny);
void init_sin_3D(double *U, double *F, double *Unew, int Nx, int Ny, int Nz);
void init_sin_mpi3D_1(double *U, double *F, double *Unew, 
			int loc_Nx, int loc_Ny, int loc_Nz, int rank, int global_Nz);
void init_sin_mpi3D_2(double *U, double *F, double *Unew, struct Information *information);

#endif // __INIT_DATA_H
