# Makefile for random model generation
default: all

# Configuration (better use configure... but I'm not sure how to use :-p)
include Makefile.config

# Categories of source files
EstRnd2d=estpsd2d.f90 
GenRnd2d=rnd2d.f90
Rnd2dSubs=header2d.f90 fft2dmod.f90 
ConvPS=convPS.f90



# A function to make object file names out of source file names
# Use it like "$(call Objects, CommonSource GenexralPurpose)" 
# - Note missing comma!
Objects=$(sort $(foreach SourceFileList,$(1),\
        $(foreach SourceFile,$($(SourceFileList)),$(patsubst %.f,%.o,\
		$(patsubst %.f90,%.o,$(SourceFile))))))

all: rndest2d rndgen2d convPS


rndest2d: 
	cd $(NRDir); $(MAKE)
	cd $(FFTDir); $(MAKE) fft2d
	cd $(MTDir); $(MAKE) 
	cd $(ACFDir); $(MAKE)
	cd $(Rnd2dDir); $(MAKE) est

rndgen2d: 
	cd $(NRDir); $(MAKE) 
	cd $(FFTDir); $(MAKE) fft2d
	cd $(MTDir); $(MAKE) 
	cd $(ACFDir); $(MAKE)
	cd $(Rnd2dDir); $(MAKE) gen

convPS:
	cd $(NRDir); $(MAKE)
	cd $(Rnd2dDir); $(MAKE) gen
	cd $(MISCDir); $(MAKE)

clean:
	cd $(NRDir); $(MAKE) clean
	cd $(FFTDir); $(MAKE) clean
	cd $(MTDir); $(MAKE) clean
	cd $(ACFDir); $(MAKE) clean
	cd $(Rnd2dDir); $(MAKE) clean
	cd $(MISCDir); $(MAKE) clean
