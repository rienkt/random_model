# Makefile for random model generation
default: all

# Configuration (better use configure... but I'm not sure how to use :-p)
include ../Makefile.config

FINCLUDE=-I$(NRDir) -I$(FFTDir) -I$(MTDir) -I$(ACFDir) -I$(Rnd2dDir)
#FLIBS=-L$(NRDir) -L$(FFTDir) -L$(MTDir) -L$(ACFDir)

convPS=convPS.f90
extractModel=extractModel.f90
nr     =nrtype.f90 nr.f90 nrutil.f90 sort.f90 static.f90 sp_func.f90 fnd_root.f90
hdrs  =header2d.f90

# main
convPSsrc=$(addprefix $(MiscDir), $(convPS))
extractModelsrc==$(addprefix $(MiscDir), $(extractModel))

# common subroutines
nrsrc=$(addprefix $(NRDir), $(nr))
#rndsrc=$(addprefix $(MTDir), $(rnd))
hdrssrc=$(addprefix $(Main2dDir), $(hdrs))

# object files
convPSobj=$(convPS:.f90=.o)
extractModelobj=$(extractModel:.f90=.o)

hdrsobj=$(addprefix $(Rnd2dDir)/, $(hdrs:.f90=.o))
nrobj=$(addprefix $(NRDir)/, $(nr:.f90=.o))
#libsobj=$(addprefix $(MTDir)/,$(rnd:.f90=.o)) $(addprefix $(FFTDir)/,$(fft:.f=.o)) $(addprefix $(ACFDir)/, $(acf:.f90=.o))
objects=$(hdrssobj) $(nrobj) $(libobj)


.SUFFIXES: .f90 .f .o
%.o: %.f90 Makefile
	$(F90) $(F90FLAGS) $(OPTF90) $(FINCLUDE) -module $(MODDIR) -c $<

%.o: %.f Makefile
	$(F77) $(F77FLAGS) $(OPTF77) -module $(MODDir) -c $<

# Dependencies

all: convPS extractModel

convPS: $(convPSobj) $(hdrsobj)
	@$(F90) -o $(BINDIR)/convPS \
	$(LDFLAGS) $(FINCLUDE) -module $(MODDIR)\
	$(FLIB) $(convPSobj) $(nrobj) $(hdrsobj) 


extractModel: $(extractModelobj) $(hdrsobj)
	@$(F90) -o $(BINDIR)/extractModel \
	$(LDFLAGS) $(FINCLUDE) -module $(MODDIR)\
	$(FLIB) $(extractModelobj) $(nrobj) $(hdrsobj) 

clean:
	$(RM) *.o *.mod
