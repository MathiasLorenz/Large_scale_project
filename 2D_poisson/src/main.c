#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include "tests.h"

double MFLOP=0.0;

int main(int argc, char const *argv[])
{
    // Handle inputs and default values
    if (argc < 3){
        fprintf(stderr, "Not enough inputs\n");
        return EXIT_FAILURE;
    }
    char const * T = argv[1];
    int N = atoi(argv[2]);

    double tol = 1e-6;
    int maxiter= 1e5;

    if (argc >= 4)  maxiter = atoi(argv[3]);
    if (argc >= 5)  tol = atof(argv[4]);

    // Handle Enviromental values
    char *problem_name, *output_info, *use_tol;

    if ( (problem_name = getenv("PROBLEM_NAME")) == NULL )
        putenv("PROBLEM_NAME=rad");
    if ( (output_info = getenv("OUTPUT_INFO")) == NULL )
        putenv("OUTPUT_INFO=timing");
    if ( (use_tol = getenv("USE_TOLERANCE")) == NULL )
        putenv("USE_TOLERANCE=on");

    // Make the call for the desired function
    double t = omp_get_wtime();
    if      (strcmp(T,"jseq") == 0) {
        fprintf(stderr, "Not available in this build. Only jomp is available.\n");
        return EXIT_FAILURE;
    }
    else if (strcmp(T,"jomp") == 0)
        test_jacobi(N, tol, maxiter);
    else if (strcmp(T,"gsseq") == 0) {
        fprintf(stderr, "Not available in this build. Only jomp is available.\n");
        return EXIT_FAILURE;
    }
    else if (strcmp(T,"jomp2") == 0) {
        fprintf(stderr, "Not available in this build. Only jomp is available.\n");
        return EXIT_FAILURE;
    }
    else if (strcmp(T,"jompft") == 0) {
        fprintf(stderr, "Not available in this build. Only jomp is available.\n");
        return EXIT_FAILURE;
    }
    else  {
        fprintf(stderr, "Invalid input\n   Accepts: jomp");
        return EXIT_FAILURE;
    }
    double timespent = omp_get_wtime() - t;

    if (strcmp("timing",getenv("OUTPUT_INFO")) == 0){
        printf("Mflops: %10.4f ", MFLOP/timespent*1e-6 );
        printf("Walltime: %10.4f\n", timespent);
    }
  return EXIT_SUCCESS;
}
