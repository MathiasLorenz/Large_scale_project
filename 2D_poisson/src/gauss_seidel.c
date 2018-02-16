#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "matrix_routines.h"
#include "poisson.h"

extern double MFLOP;

void gauss_seidel(int N_total, int maxit, double threshold,
        double **u, double **f)
{
    if(!u || !f) {fprintf(stderr,"Pointer is NULL.\n"); return;}
    int N = N_total-2;
    int i, j, k = 0;
    double frac = 1.0/4.0;
    double h = 2.0/(N+1.0);
    double step = h*h;
    double norm0 = 0.0, norm1 = 0.0;
    double norm_diff = 10.0;

    // Jacobi iteration

    bool use_tol = false;
    if (strcmp("on",getenv("USE_TOLERANCE")) == 0)
        use_tol = true;

    for(k = 0; k < maxit ; k++){
        for(i=1; i<N+1; i++){
            for(j=1; j<N+1; j++){
                u[i][j] = frac*(u[i-1][j]+u[i+1][j]+
                    u[i][j-1]+u[i][j+1]+step*f[i][j]);
                    norm1 += u[i][j]*u[i][j];
                }
            }
            // Calclulate frobenius norm
            norm1 = sqrt(norm1);

            norm_diff = fabs(norm0 - norm1);
            swap_double(&norm0, &norm1);

        if (use_tol && (norm_diff < threshold))
            break;
    }

    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0){
        // Print exit cause
        if(norm_diff < threshold && use_tol)
            fprintf(stdout, "Exited because norm > threshold\n");
        else
            fprintf(stdout, "Exited because iter > maxit\n");
    }

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0){
        MFLOP = (8.0*N*N + 4.0)*k;
    }
}
