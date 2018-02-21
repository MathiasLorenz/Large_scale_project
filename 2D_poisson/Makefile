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
EXEC    = jacobiSolver

# - Source files
COBJS   = main tests init_data matrix_routines
COBJS  += jacobi_openmp_2D
CUOBJS	= tests_cuda cuda_routines

# - Directories
SRCDIR	= src
LIBDIR  = lib
OBJDIR  = obj

# - Combining files to their correct directory
SOURCES  = $(addprefix $(SRCDIR)/, $(addsuffix .c, $(COBJS)))
SOURCES	+= $(addprefix $(SRCDIR)/, $(addsuffix .cu,$(CUOBJS)))
OBJECTS  = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(COBJS)))
OBJECTS += $(addprefix $(OBJDIR)/, $(addsuffix .o, $(CUOBJS)))

# -- End of File definitions
#==============================================================================
# -- Compilations                             (Avoid changes beyond this point)

# - Main
$(EXEC) : $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) -I$(LIBDIR) $(LFLAGS) $(LUFLAGS)
	@echo -e "\nmake: '$(EXEC)' was build sucessfully."

# - C Files
.SUFFIXES: .c .cu
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) -o $@ -c $< -I$(LIBDIR) $(CFLAGS)
$(OBJDIR)/%.o: $(SRCDIR)/%.cu | $(OBJDIR)
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

src/main.o: lib/tests.h lib/cuda_tests.h
src/tests.o: lib/tests.h lib/matrix_routines.h lib/poisson.h lib/init_data.h
src/init_data.o: lib/init_data.h
src/matrix_routines.o: lib/matrix_routines.h
src/jacobi_openmp.o: lib/matrix_routines.h lib/poisson.h
src/cuda_tests.o: lib/cuda_tests.h lib/cuda_routines.h
src/cuda_routines.o: lib/cuda_routines.h
