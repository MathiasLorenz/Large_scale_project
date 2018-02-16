#==============================================================================
# -- Definitions

# - Compilers
CC	= gcc

# - C flags
OPTIM    = -O3
WARN     = -Wall
LIBS     = -std=gnu99
CHIP     =

CFLAGS   = $(OPTIM) $(WARN) $(LIBS) $(CHIP)

# -- End of Definitions
#==============================================================================
# -- File definitions and macros

EXEC     = test
COBJS    = main cfile

SRCDIR	 = src
LIBDIR   = lib
OBJDIR   = objects

SOURCES  = $(addprefix $(SRCDIR)/, $(addsuffix .c, $(COBJS)))
OBJECTS  = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(COBJS)))

# -- End of File definitions
#==============================================================================
# -- Compilations                             (Avoid changes beyond this point)

# - Main
$(EXEC) : $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) -I$(LIBDIR) $(CFLAGS)
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
	@MakeDepend -Y$(LIBDIR) $(SOURCES)
	@rm -f makefile.bak MDP.err

# -- End of utility commands
#==============================================================================
# -- Compile dependecies
# DO NOT DELETE

src/main.o: lib/hfile.h
src/cfile.o: lib/hfile.h
