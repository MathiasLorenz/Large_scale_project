#==============================================================================
# -- Definitions

# - Compilers
CC	= mpicc
CU	= nvcc

# - C Flags
OPTIM   =-O3
WARN    =-Wall
LIBS    =-std=gnu99
CHIP    =
PARA    =-fopenmp
LFLAGS  =-fopenmp -lm

CFLAGS  = $(OPTIM) $(WARN) $(LIBS) $(CHIP) $(PARA)

# - CU Flags
CUOPT 	=-O3
CUOPTS	=
CUWARN	=-Wall
CULIBS	=-lcudart -L/appl/cuda/9.1/lib64
CUARCH	=-arch=sm_30
CUPARA	=

# - CU Library setup
CUDA_PATH ?=/appl/cuda/9.1
CUINCLUDES =-I$(CUDA_PATH)/include -I$(CUDA_PATH)/samples/common/inc

CUFLAGS=-Xcompiler "$(CUOPT) $(CUPIC) $(CUOMP)" $(CUARCH) $(CUOPTS)
LUFLAGS=$(CUINCLUDES) $(CULIBS)


# -- End of Definitions
#==============================================================================
# -- File definitions and macros

# - Excecutable
EXEC    = jacobiSolver.bin

# - Source files
UOBJS   = main tests init_data matrix_routines
JOBJS   = jacobi_openmp_2D jacobi_openmp_3D jacobi_mpi3D_1 jacobi_mpi3D_2 jacobi_util
CUUOBJS = tests_cuda cuda_routines
CUJOBJS = jacobi_cuda_1 jacobi_util_cuda

# - Directories
SRCDIR	= src
LIBDIR	= lib
OBJDIR	= obj
JDIR 	= $(SRCDIR)/jacobi
UDIR 	= $(SRCDIR)/util

# - Combining files to their correct directory
SOURCES	 = $(addprefix $(JDIR)/, $(addsuffix .c, $(JOBJS)))
SOURCES	+= $(addprefix $(JDIR)/, $(addsuffix .cu, $(CUJOBJS)))
SOURCES	+= $(addprefix $(UDIR)/, $(addsuffix .c, $(UOBJS)))
SOURCES	+= $(addprefix $(UDIR)/, $(addsuffix .cu, $(CUUOBJS)))

OBJECTS	 = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(UOBJS) $(JOBJS)))
OBJECTS	+= $(addprefix $(OBJDIR)/, $(addsuffix .o, $(CUUOBJS) $(CUJOBJS)))

# -- End of File definitions
#==============================================================================
# -- Compilations                             (Avoid changes beyond this point)

# - Main
$(EXEC) : $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) -I$(LIBDIR) $(LFLAGS) $(LUFLAGS)
	@echo -e "\nmake: '$(EXEC)' was build sucessfully."

# - C Files
.SUFFIXES: .c .cu
$(OBJDIR)/%.o: $(JDIR)/%.c | $(OBJDIR)
	$(CC) -o $@ -c $< -I$(LIBDIR) $(CFLAGS)
$(OBJDIR)/%.o: $(UDIR)/%.c | $(OBJDIR)
	$(CC) -o $@ -c $< -I$(LIBDIR) $(CFLAGS)
$(OBJDIR)/%.o: $(JDIR)/%.cu | $(OBJDIR)
	$(CU) -o $@ -c $< -I$(LIBDIR) $(CUINCLUDES) $(CUFLAGS)
$(OBJDIR)/%.o: $(UDIR)/%.cu | $(OBJDIR)
	$(CU) -o $@ -c $< -I$(LIBDIR) $(CUINCLUDES) $(CUFLAGS)

# -- End of compilations
#==============================================================================
# -- Utility commands

$(OBJDIR):
	@mkdir -p $(OBJDIR)
clean :
	@rm -fr $(OBJDIR) core
realclean : clean
	@rm -f $(EXEC)
depend :
	makedepend -Y$(LIBDIR) $(SOURCES)
	@rm -f Makefile.bak MDP.err

# -- End of utility commands
#==============================================================================
# -- Compile dependecies
# DO NOT DELETE

src/jacobi/jacobi_openmp_2D.o: lib/matrix_routines.h lib/jacobi_util.h
src/jacobi/jacobi_openmp_2D.o: lib/poisson.h
src/jacobi/jacobi_openmp_3D.o: lib/matrix_routines.h lib/jacobi_util.h
src/jacobi/jacobi_openmp_3D.o: lib/poisson.h
src/jacobi/jacobi_mpi3D_1.o: lib/matrix_routines.h lib/jacobi_util.h
src/jacobi/jacobi_mpi3D_1.o: lib/poisson.h
src/jacobi/jacobi_mpi3D_2.o: lib/matrix_routines.h lib/jacobi_util.h
src/jacobi/jacobi_mpi3D_2.o: lib/poisson.h
src/jacobi/jacobi_util.o: lib/matrix_routines.h lib/jacobi_util.h
src/jacobi/jacobi_cuda_1.o: lib/poisson_cuda.h lib/matrix_routines.h
src/jacobi/jacobi_cuda_1.o: lib/jacobi_util.h lib/cuda_routines.h
src/jacobi/jacobi_cuda_1.o: lib/jacobi_util_cuda.h
src/jacobi/jacobi_util_cuda.o: lib/matrix_routines.h lib/jacobi_util.h
src/jacobi/jacobi_util_cuda.o: lib/cuda_routines.h
src/util/main.o: lib/tests.h lib/jacobi_util.h lib/tests_cuda.h
src/util/tests.o: lib/tests.h lib/jacobi_util.h lib/poisson.h lib/init_data.h
src/util/tests.o: lib/matrix_routines.h
src/util/matrix_routines.o: lib/jacobi_util.h lib/matrix_routines.h
src/util/tests_cuda.o: lib/init_data.h lib/jacobi_util.h lib/tests_cuda.h
src/util/tests_cuda.o: lib/poisson_cuda.h lib/matrix_routines.h
src/util/tests_cuda.o: lib/cuda_routines.h
