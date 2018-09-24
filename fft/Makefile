# Makefile for random model generation
default: all
include ../Makefile.config

# Categories of source files
fft1dsrcs=fft1d/
fft2dsrcs=fft2d/fftsg2d.f fft2d/fftsg.f

# suffix
#.SUFFIXES: .f90 .f .o

%.o: %.f90 Makefile
	$(F90) $(F90FLAGS) $(OPTF90) -module $(MODDIR) -c $<

%.o: %.f Makefile
	$(F77) $(F77FLAGS) $(OPTF77) -c $<

# target
all: 

fft1d: $(fft1dsrcs:.f=.o)

fft2d: $(fft2dsrcs:.f=.o)

clean:
	$(RM) *.o *.mod
