# put all things related to the compiler here
CC       = suncc
OPTIM    = -g -fast -xrestrict
WARN     =
LIBS     =
CHIP     =
PARA_A	 = -xautopar -xreduction -xloopinfo -xvpara -xipo
PARA_M	 = -xopenmp -xloopinfo -xvpara
LFLAGS	 = -xopenmp
CFLAGS   = $(OPTIM) $(WARN) $(LIBS) $(CHIP)


# -- File definitions and macros
EXEC     = jacobiSolver
MOBJECTS = main.o matrix_routines.o tests.o init_data.o
ALGOS    = jacobi_sequential.o jacobi_openmp.o gauss_seidel.o \
		   jacobi_vectorize.o jacobi_openmpv2.o

OBJECTS  = $(MOBJECTS) $(ALGOS)
# -- Compile main program
$(EXEC) : $(OBJECTS)
	$(CC) $(CFLAGS) $(LFLAGS) -o $(EXEC) $(OBJECTS)


jacobi_sequential.o :
init_data.o : CFLAGS += $(PARA_M)
jacobi_openmp.o : CFLAGS += $(PARA_M)
jacobi_openmpv2.o : CFLAGS += $(PARA_M)

# -- Utility commands
clean :
	@rm -f *.o core
realclean : clean
	@rm -f $(EXEC)

# -- Compile dependecies
# DO NOT DELETE

gauss_seidel.o: matrix_routines.h poisson.h
init_data.o: init_data.h
jacobi_openmp.o: matrix_routines.h poisson.h
jacobi_openmpv2.o: matrix_routines.h poisson.h
jacobi_sequential.o: matrix_routines.h poisson.h
jacobi_vectorize.o: matrix_routines.h poisson.h
main.o: tests.h
matrix_routines.o: matrix_routines.h
tests.o: tests.h matrix_routines.h poisson.h init_data.h
