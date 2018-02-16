#==============================================================================
# -- Definitions

# - Compilers
CC	=gcc

# - C flags
OPTIM    =-O3
WARN     =-Wall
LIBS     =-std=gnu99 -lm
CHIP     =
PARA     =-fopenmp
LFLAGS   =-fopenmp

CFLAGS   = $(OPTIM) $(WARN) $(LIBS) $(CHIP) $(PARA)

# -- End of Definitions
#==============================================================================
# -- File definitions and macros

EXEC     = jacobiSolver
COBJS    = main tests init_data matrix_routines gauss_seidel
COBJS   += jacobi_sequential jacobi_openmp jacobi_openmpv2 jacobi_vectorize

SRCDIR	 = src
LIBDIR   = lib
OBJDIR   = obj

SOURCES  = $(addprefix $(SRCDIR)/, $(addsuffix .c, $(COBJS)))
OBJECTS  = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(COBJS)))

# -- End of File definitions
#==============================================================================
# -- Compilations                             (Avoid changes beyond this point)

# - Main
$(EXEC) : $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) -I$(LIBDIR) $(CFLAGS) $(LFLAGS)
	@echo -e "\nmake: '$(EXEC)' was build sucessfully."

# - C Files
.SUFFIXES: .c
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) -o $@ -c $< -I$(LIBDIR) $(CFLAGS)

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

src/main.o: lib/tests.h
src/tests.o: lib/tests.h lib/matrix_routines.h lib/poisson.h lib/init_data.h
src/init_data.o: lib/init_data.h
src/matrix_routines.o: lib/matrix_routines.h
src/gauss_seidel.o: lib/matrix_routines.h lib/poisson.h
src/jacobi_sequential.o: lib/matrix_routines.h lib/poisson.h
src/jacobi_openmp.o: lib/matrix_routines.h lib/poisson.h
src/jacobi_openmpv2.o: lib/matrix_routines.h lib/poisson.h
src/jacobi_vectorize.o: lib/matrix_routines.h lib/poisson.h
