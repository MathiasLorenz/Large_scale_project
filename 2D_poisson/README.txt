README for jacobiSolver
-----------------------
The Jacobi solver have been implemented into the jacobiSolver executable.

Default call to the function will be
    > ./jacobiSolver METHOD GRIDSIZE [MAXITER [TOLERANCE]]
where:
    METHOD  : jseq, jomp, gsseq.
        jseq  : Is the sequential Jacobi method.
        jomp  : Is the parallelized Jacobi method using OpenMP.
        gsseq : Is the sequential Gauss-Seidel method.
        jomp2 : Attempt to improve the parallel version.
        jompft: Parallelized version using First Touch.

    GRIDSIZE : Integer number.
        The number of inner points in the problem, fx. 5
            B B B B B B B B
            B X X X X X X B
            B X X X X X X B
            B X X X X X X B
            B X X X X X X B
            B X X X X X X B
            B B B B B B B B
        Where the boundary is set during the function and the elemnets in X is
        computed.

    MAXITER[optional] : Integer number.
        The maximal number of iterations done by the solver.
        Default is: 10^(5)

    TOLERANCE[optional] : Floating point number
        Sets the tolerance for the methods. The methods will terminate depending
        on the method:
            jseq  : When the change norm_F(U_new - U_old) falls below the
                tolerance. The norm used is Frobenius norm.
            jomp  : When the change norm_F(U_new - U_old) falls below the
                tolerance. The norm used is Frobenius norm.
            gsseq : When the change ( norm_F(U_new) - norm_F(U_old) ) falls below
                the tolerance. The norm used is Frobenius norm.

The solver will read for several environmental variables and change behavior
based on these. The environmental variables supported is:

    OpenMP: All default OpenMP controls such as setting number of threads
            are present and jomp have been implemented such that the schedule
            can be set during Runtime using OMP_SCHEDULE.

    Extra environmental variables: ( > ./jacobiSolver ... ENV_NAME=value )
        PROBLEM_NAME: [default: rad]
            Defines the problem that should to be solved.
            sin : The well defined example where the Boundary condition and
                  source is defined as:
                        f(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)
                        u(x,y) = 0, for (x,y) on boundary.
            rad : The radiator problem with Boundary and source is defined as:
                        f(x,y) = 200 for x in [ 0, 1/3 ] and y in [ -2/3, -1/3 ]
                        f(x,y) = 0 otherwise.
                        u(x,1) = 20, u(x ,-1) =  0,
                        u(1,y) = 20, u(-1, y) = 20.

        OUTPUT_INFO: [default: timing]
            Defines the output type.
            matrix: Defines that the result matrix should be printed as the
                    output.
            timing: Defines that the output should be the timing info. Giving
                    the memory footprint, Mflops and Walltime for the main loop.

        USE_TOLERANCE : Defines if the tolerance should be used or not.
            on :    Use the tolerance as stop criterion.
            off:    Do not use the tolerance. Force to do MAXITER iterations.

% = = EOF = = %
