.PHONY: clean link rm-link

FC=ifort
#FFLAGS=-fpp -O3 -xHost -ipo
#FFLAGS=-C -O0 -g -traceback -fpe0 -ftrapuv
LIBS=-lfftw3 -lm
LDFLAGS=

# Switches for debug mode
ifneq ($(DEBUG),)
	CPPFLAGS:=-D_VERBOSE $(CPPFLAGS)
	FFLAGS=-O1 -g -traceback -fpe0 -ftrapuv -fpp
else
	CPPFLAGS=
	FFLAGS=-fpp -O3 -xHost -ipo
endif

# Debug just FFT
ifeq (fft,$(findstring $(DEBUG),fft))
	CPPFLAGS:=-D_VERBOSE_FFT $(CPPFLAGS)
	FFLAGS=-O1 -g -traceback -fpe0 -ftrapuv -fpp
endif

# Create a dependency list using makedepf90.  All files  
# that needs to be compiled to build the program,  
# i.e all source files except include files, should  
# be given on the command line to makedepf90.   
# 
# The argument to the '-o' option will be the name of the 
# resulting program when running 'make', in this case  
# 'foobar' 
.depend: 
	makedepf90 -W -o Surface.x *.f90 > .depend

# Include the dependency-list created by makedepf90 below 
-include .depend

%.o: %.f90 
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $<


# target 'clean' for deleting object- *.mod- and other  
# unwanted files 
clean: 
	rm -f *.o *.mod core .depend

link:
	ln -fs ${PWD}/Surface.x ${HOME}/bin/

rm-link:
	rm -f ${HOME}/bin/Surface.x


