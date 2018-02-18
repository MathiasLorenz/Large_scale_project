#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "matrix_routines.h"
#include "poisson.h"

extern double MFLOP;

void jacobi_openmp(int N_total, int maxit, double threshold,
    double *u, double *f, double *tmp)
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


        #pragma omp parallel for\
            private(i,j) reduction(+: norm_diff) schedule(runtime)
        for(i=1; i<N+1; i++){
            for(j=1; j<N+1; j++){
                // Linear indexing with macro
                ACCESS_2D(tmp, i, j, N_total) =
                        frac*( ACCESS_2D(u, i-1, j, N_total)
                             + ACCESS_2D(u, i+1, j, N_total)
                             + ACCESS_2D(u, i, j-1, N_total)
                             + ACCESS_2D(u, i, j+1, N_total)
                             + step*ACCESS_2D(f, i, j, N_total));

                double uij   = ACCESS_2D(u, i, j, N_total);
                double tmpij = ACCESS_2D(tmp, i, j, N_total);
                norm_diff += (uij-tmpij)*(uij-tmpij);
            }
        }
        //swap_matrix(&tmp,&u);
        swap_array(&tmp,&u);
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
