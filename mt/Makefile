# Makefile
default: all
include ../Makefile.config

#
srcs=mt19937_mod.f90

# suffix
.SUFFIXES: .f90 .f .o

%.o: %.f90 Makefile
	$(F90) $(F90FLAGS) $(OPTF90) -module $(MODDIR) -c $<

%.o: %.f Makefile
	$(F77) $(F77FLAGS) $(OPTF77) -c $<

all: $(srcs:.f90=.o)

clean:
	 $(RM) *.o *.mod
