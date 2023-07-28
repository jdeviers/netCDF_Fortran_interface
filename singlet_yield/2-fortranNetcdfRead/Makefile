# Source files
PROG = $(wildcard prog_*.f90)

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
PROGRAM = openNetCDF.x

default: $(PROGRAM)

$(PROGRAM):
	$(FC) $(NCFLAG1) -o $@ $(PROG) $(FCFLAGS) $(NCFLAG2)


clean:
	rm $(PROGRAM) $(OBJS) $(wildcard *.mod)

debug:
	@echo "NetCDF flags: $(NCFLAG1)"
	@echo "NetCDF flags: $(NCFLAG2)"
	@echo "Clean command: rm $(PROGRAM)"

.PHONY: default clean debug
