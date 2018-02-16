#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "matrix_routines.h"
#include "poisson.h"

#define ACCESS_2D(A, r, c, rC) (A[(r)*(rC) + c])

extern double MFLOP;

void jacobi_openmp(int N_total, int maxit, double threshold,
    double **u, double **f, double **tmp)
{
    if(!u || !f || !tmp) {fprintf(stderr,"Pointer is NULL.\n"); return;}
    int N = N_total-2;
    int i, j, k = 0;
    double frac = 1.0/4.0;
    double h = 2.0/(N+1.0);
    double step = h*h;
    double norm_diff = 10.0;



    // Jacobi iteration
    bool use_tol = false;
    if (strcmp("on",getenv("USE_TOLERANCE")) == 0)
        use_tol = true;

    for(k = 0; k < maxit ; k++){
        norm_diff = 0.0;
        // Linear indexing
        double * u_lin = *u;

        #pragma omp parallel for\
            private(i,j) reduction(+: norm_diff) schedule(runtime)
        for(i=1; i<N+1; i++){
            for(j=1; j<N+1; j++){
                //tmp[i][j] = frac*(u[i-1][j]+u[i+1][j]+
                //    u[i][j-1]+u[i][j+1]+step*f[i][j]);

                // With linear indexing
                //tmp[i][j] -> tmp[i*N + j]
                //tmp[i][j] = frac*(u_lin[(i-1)*N + j] + u_lin[(i+1)*N + j]
                //                + u_lin[i*N + j-1] + u_lin[i*N + j]
                //                + step*f[i][j]);

                // Linear indexing and macro
                tmp[i][j] = frac*( ACCESS_2D(u_lin, i-1, j, N)
                                 + ACCESS_2D(u_lin, i+1, j, N)
                                 + ACCESS_2D(u_lin, i, j-1, N)
                                 + ACCESS_2D(u_lin, i, j+1, N)
                                 + step*f[i][j]);

                norm_diff += (u[i][j]-tmp[i][j])*(u[i][j]-tmp[i][j]);
            }
        }
        swap_matrix(&tmp,&u);
        norm_diff = sqrt(norm_diff);

        if (use_tol && (norm_diff < threshold))
            break;
    }

    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0){
        // Print exit cause
        if(norm_diff < threshold && use_tol)
            fprintf(stdout, "Exited because norm < threshold\n");
        else
            fprintf(stdout, "Exited because iter = maxit\n");
    }

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0){
        MFLOP = (10.0*N*N + 4.0)*k;
    }
}
