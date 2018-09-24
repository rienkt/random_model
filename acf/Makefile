# Makefile for random model generation
default: all
include ../Makefile.config

# Categories of source files
srcs=acf.f90
objs=$(srcs:.f90=.o)

F90INCLUDE=-I$(NRDir)

# suffix
.SUFFIXES: .f90 .f .o

%.o: %.f90 Makefile
	$(F90) $(F90FLAGS) $(OPTF90) $(F90INCLUDE) -module $(MODDIR) -c $<

%.o: %.f Makefile
	$(F77) $(F77FLAGS) $(OPTF77) -c $<

# target
all: $(objs)

clean:
	$(RM) *.o *.mod
