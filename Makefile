.PHONY: clean link rm-link

FC=ifort
FFLAGS=-fast -ftrapuv
#FFLAGS=-O0 -g -traceback -fpe0 -ftrapuv
LIBS=-lfftw3 -lm
LDFLAGS=

# Include the dependency-list created by makedepf90 below 
include .depend


%.o: %.f90 
	$(FC) -c $(FFLAGS) -o $@ $<


# target 'clean' for deleting object- *.mod- and other  
# unwanted files 
clean: 
	rm -f *.o *.mod core .depend

link:
	ln -fs ${PWD}/Surface.x ${HOME}/bin/

rm-link:
	rm -f ${HOME}/bin/Surface.x


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
