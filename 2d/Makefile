# Makefile for random model generation
default: all

# Configuration (better use configure... but I'm not sure how to use :-p)
include ../Makefile.config

FINCLUDE=-I$(NRDir) -I$(FFTDir) -I$(MTDir) -I$(ACFDir)
#FLIBS=-L$(NRDir) -L$(FFTDir) -L$(MTDir) -L$(ACFDir)

est    =estpsd2d.f90
gen    =rnd2d.f90
genagu =rnd2d_agu.f90

hdrs   =header2d.f90
subs   =fft2dmod.f90 subrnd2d.f90
nr     =nrtype.f90 nr.f90 nrutil.f90 sort.f90 static.f90 sp_func.f90 fnd_root.f90
rnd    =mt19937_mod.f90
fft    =fftsg.f fftsg2d.f
acf    =acf.f90

# main
gensrc=$(addprefix $(Rnd2dDir), $(gen))
estsrc=$(addprefix $(Rnd2dDir), $(est))
hdrssrc=$(addprefix $(Rnd2dDir), $(hdrs))
subssrc=$(addprefix $(Rnd2dDir), $(subs))

# common subroutines
nrsrc=$(addprefix $(NRDir), $(nr))
rndsrc=$(addprefix $(MTDir), $(rnd))

genaguobj=$(genagu:.f90=.o)
genobj=$(gen:.f90=.o)
estobj=$(est:.f90=.o)
hdrsobj=$(hdrs:.f90=.o)
subsobj=$(subs:.f90=.o)
nrobj=$(addprefix $(NRDir)/, $(nr:.f90=.o))
libsobj=$(addprefix $(MTDir)/,$(rnd:.f90=.o)) $(addprefix $(FFTDir)/,$(fft:.f=.o)) $(addprefix $(ACFDir)/, $(acf:.f90=.o))
objects=$(subsobj) $(nrobj) $(libobj)


.SUFFIXES: .f90 .f .o
%.o: %.f90 Makefile
	$(F90) $(F90FLAGS) $(OPTF90) $(FINCLUDE) -module $(MODDIR) -c $<

%.o: %.f Makefile
	$(F77) $(F77FLAGS) $(OPTF77) -module $(MODDir) -c $<


# Dependencies

all: gen est

gen: $(hdrsobj) $(subsobj) $(genobj)
	@$(F90) -o $(BINDIR)/rndgen2d \
	$(LDFLAGS) $(FINCLUDE) -module $(MODDIR)\
	$(FLIB) $(genobj) $(nrobj) $(libsobj) $(hdrsobj) $(subsobj) 

est: $(hdrsobj) $(subsobj)  $(estobj)
	@$(F90) -o $(BINDIR)/rndest2d \
	$(LDFLAGS) $(FINCLUDE) -module $(MODDIR)\
	$(FLIB) $(estobj) $(nrobj) $(libsobj) $(hdrsobj) $(subsobj) 

genagu: $(genaguobj) $(hdrsobj) $(subsobj) 
	@$(F90) -o $(BINDIR)/rndgenagu2d \
	$(LDFLAGS) $(FINCLUDE) -module $(MODDIR)\
	$(FLIB) $(genaguobj) $(nrobj) $(libsobj) $(hdrsobj) $(subsobj) 

clean:
	$(RM) *.o *.mod
