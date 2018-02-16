#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "matrix_routines.h"
#include "poisson.h"

extern double MFLOP;

void jacobi_openmpv2(int N_total, int maxit, double threshold,
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

        k = 0;
  #pragma omp parallel if (N > 200)
  {
    for( ; k < maxit ;){
        #pragma omp single
        {
          norm_diff = 0.0;
        }
        #pragma omp for\
            private(i,j) reduction(+: norm_diff) schedule(runtime)
        for(i=1; i<N+1; i++){
            for(j=1; j<N+1; j++){
                tmp[i][j] = frac*(u[i-1][j]+u[i+1][j]+
                    u[i][j-1]+u[i][j+1]+step*f[i][j]);
                norm_diff += (u[i][j]-tmp[i][j])*(u[i][j]-tmp[i][j]);
            }
        }

        #pragma omp single
        {
          swap_matrix(&tmp,&u);
          norm_diff = sqrt(norm_diff);
          k++;
        }

        //if (use_tol && (norm_diff < threshold))
        //    break;
    }
  }

    if (strcmp("matrix",getenv("OUTPUT_INFO")) == 0){
        if(norm_diff < threshold && use_tol)
        // Print exit cause
            fprintf(stdout, "Exited because norm < threshold\n");
        else
            fprintf(stdout, "Exited because iter = maxit\n");
            //fprintf(stdout, "norm: %lf\n",norm_diff);
    }
    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0){
        MFLOP = (10.0*N*N + 4.0)*k;
    }
}
