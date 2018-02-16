#ifndef __INIT_DATA_H
#define __INIT_DATA_H

void init_sin(double **u, double **f, double **R, int N_total, double h);
void init_rad(double **u, double **f, double **R, int N_total, double h);
void init_sin_gs(double **u, double **f, int N_total, double h);
void init_rad_gs(double **u, double **f, int N_total, double h);
void init_sin_omp(double **u, double **f, double **R, int N_total, double h);
void init_rad_omp(double **u, double **f, double **R, int N_total, double h);
#endif // __INIT_DATA_H
