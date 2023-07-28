# Source files
PROG = prmtop.f90
MODS = $(wildcard mod_*.f90)
OBJS = $(patsubst %.f90,%.o,$(MODS))

# Compiler
FC = gfortran

# Fortran prod and debug flags:
FCFLAGS = -O2 -fbacktrace -Wall -Wextra -Wno-tabs -fimplicit-none
DBFLAGS = -O0 -fbacktrace -fcheck-bounds

# NetCDF compilation flags
NC_INCL = $(shell nc-config --includedir)
NC_LIBS = $(shell nc-config --libs)
NC_ADDT = -lnetcdff
NCFLAG1 = -I$(NC_INCL)
NCFLAG2 =   $(NC_LIBS) $(NC_ADDT)

# Target
PROGRAM = prmtop.x

default: $(PROGRAM)
#	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(NCFLAG1) -o $@ $(PROG) $^ $(FCFLAGS) $(NCFLAG2)

$(OBJS): %.o : %.f90
	$(FC) $(NCFLAG1) -c $< $(FCFLAGS) $(NCFLAG2)

mod_interfaces.o mod_readprmtop.o mod_readnetcdf.o : mod_procedures.o
mod_readprmtop.o mod_readnetcdf.o : mod_interfaces.o

clean:
	rm $(PROGRAM) $(OBJS) $(wildcard *.mod)

debug:
	@echo "Modules sources: $(MODS)"
	@echo "Modules objects: $(OBJS)"
	@echo "NetCDF flags: $(NCFLAG1)"
	@echo "NetCDF flags: $(NCFLAG2)"
	@echo "Clean command: rm $(PROGRAM) $(OBJS) $(wildcard *.mod)"

.PHONY: default clean debug
