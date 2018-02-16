#include <stdio.h>
#include <math.h>
#include "init_data.h"

void init_sin(double **u, double **f, double **R, int N_total, double h){
    if(!f || !u || !R) {fprintf(stderr,"Pointer is NULL.\n"); return;}

    double M_PI_sq = M_PI*M_PI;
    double x = -1.0, y = -1.0, res = 0.0;
    for(int i = 0; i < N_total; i++) {
        x = -1.0;
        for(int j = 0; j < N_total; j++) {
            res = 2.0*M_PI_sq*sin(M_PI*x)*sin(M_PI*y);
            f[i][j] = res;
            x += h;
        }
        y += h;
    }

    for(int i = 0; i < N_total*N_total; i++) {
        u[0][i] = 0.0;
        R[0][i] = 0.0;
    }
}
void init_rad(double **u, double **f, double **R, int N_total, double h){

    if(!f || !u || !R) {fprintf(stderr,"Pointer is NULL.\n"); return;}


    for(int i = 0; i < N_total*N_total; i++)
        f[0][i] = 0.0;

    for(int j = lround(1.0/h); j <= lround(4.0/(h*3.0)); j++) {
        for(int i = lround(1.0/(h*3.0)); i <= lround(2.0/(h*3.0)); i++) {
            f[i][j] = 200;
        }
    }

    // Set u to zeros
    for(int i = 0; i < N_total*N_total; i++) {
        u[0][i] = 0.0;
        R[0][i] = 0.0;
    }

    // Set BC of u
    for(int i = 0; i < N_total; i++) {
        u[i][N_total-1] = R[i][N_total-1]   = 20;
        u[N_total-1][i] = R[N_total-1][i]   = 20;
        u[i][0]         = R[i][0]           = 20;
    }

}


void init_sin_gs(double **u, double **f, int N_total, double h){
    if(!f || !u) {fprintf(stderr,"Pointer is NULL.\n"); return;}

    double M_PI_sq = M_PI*M_PI;
    double x = -1.0, y = -1.0, res = 0.0;
    for(int i = 0; i < N_total; i++) {
        x = -1.0;
        for(int j = 0; j < N_total; j++) {
            res = 2.0*M_PI_sq*sin(M_PI*x)*sin(M_PI*y);
            f[i][j] = res;
            x += h;
        }
        y += h;
    }

    for(int i = 0; i < N_total*N_total; i++) {
        u[0][i] = 0.0;
    }
}
void init_rad_gs(double **u, double **f, int N_total, double h){

    if(!f || !u) {fprintf(stderr,"Pointer is NULL.\n"); return;}


    for(int i = 0; i < N_total*N_total; i++)
        f[0][i] = 0.0;

    for(int i = lround(1.0/h+0.5); i <= lround(4.0/(h*3.0)-0.5); i++) {
        for(int j = lround(1.0/(h*3.0)+0.5); j <= lround(2.0/(h*3.0)-0.5); j++) {
            f[j][i] = 200;
        }
    }

    // Set u to zeros
    for(int i = 0; i < N_total*N_total; i++)
        u[0][i] = 0.0;

    // Set BC of u
    for(int i = 0; i < N_total; i++) {
        u[i][N_total-1] = 20;
        u[N_total-1][i] = 20;
        u[i][0]         = 20;
    }

}

void init_sin_omp(double **u, double **f, double **R, int N_total, double h){
    if(!f || !u || !R) {fprintf(stderr,"Pointer is NULL.\n"); return;}


    double M_PI_sq = M_PI*M_PI;
    double x = 0.0, y = 0.0;

    int i, j;
    #pragma omp parallel for private(i,j,x,y) shared(M_PI_sq, f, u, R)
    for(i = 1; i < N_total-1; i++) {
        x = -1.0;
        y = -1.0 + i*h;
        for(j = 1; j < N_total-1; j++) {
            x += h;
            f[i][j] = 2.0*M_PI_sq*sin(M_PI*x)*sin(M_PI*y);
            u[i][j] = 0.0;
            R[i][j] = 0.0;
        }
    } // End of parallel

    // Set BC
    for(int i = 0; i < N_total; i++) {
        //u[0][i] = 0.0;
        //R[0][i] = 0.0;
        u[0][i] = R[0][i] = f[0][i] = 0;
        u[i][0] = R[i][0] = f[i][0] = 0;
        u[N_total-1][i] = R[N_total-1][i] = f[N_total-1][i] = 0;
        u[i][N_total-1] = R[i][N_total-1] = f[i][N_total-1] = 0;
    }
}
void init_rad_omp(double **u, double **f, double **R, int N_total, double h){

    if(!f || !u || !R) {fprintf(stderr,"Pointer is NULL.\n"); return;}

    // Initialize interior points with first touch
    int i, j;
    #pragma omp parallel for private(i,j) shared(f, u, R)
    for(i = 1; i < N_total-1; i++) {
        for(j = 1; j < N_total-1; j++) {
            f[i][j] = 0.0;
            u[i][j] = 0.0;
            R[i][j] = 0.0;
        }
    } // End parallel


    // Update interior points for source
    for(int j = lround(1.0/h); j <= lround(4.0/(h*3.0)); j++) {
        for(int i = lround(1.0/(h*3.0)); i <= lround(2.0/(h*3.0)); i++) {
            f[i][j] = 200;
        }
    }


    // Set BC of all
    for(int i = 0; i < N_total; i++) {
        u[i][N_total-1] = R[i][N_total-1]   = 20;
        u[N_total-1][i] = R[N_total-1][i]   = 20;
        u[i][0]         = R[i][0]           = 20;
        u[0][i]         = R[0][i]           = 0.0;
        f[i][0] = f[0][i] = f[i][N_total-1] = f[N_total-1][i] = 0.0;
    }

}
