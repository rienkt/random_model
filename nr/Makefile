# Makefile for random model generation
default: all
include ../Makefile.config

# Categories of source files
NRType=nrtype.f90
NRCommon=nr.f90 nrutil.f90
NRSubs=sort.f90 static.f90 sp_func.f90 fnd_root.f90

# suffix
.SUFFIXES: .f90 .f .o .a
#%.o: %.f90 Makefile

# target
all: $(NRType:.f90=.o) $(NRCommon:.f90=.o) $(NRSubs:.f90=.o) ../Makefile.config

.f90.o:
	$(F90) $(F90FLAGS) $(OPTF90) -module $(MODDIR) -c $<

.f77.o:
	$(F77) $(F77FLAGS) $(OPTF77) -c $<

.o.a:
	ar -r $@ $<


clean:
	$(RM) *.o *.mod
