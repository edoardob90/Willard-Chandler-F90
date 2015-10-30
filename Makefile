.PHONY: clean

FC=ifort
FFLAGS=-C -g -O0
LIBS=

# Suffix-rules:  Begin by throwing away all old suffix- 
# rules, and then create new ones for compiling  
# *.f90-files. 
.SUFFIXES: 
.SUFFIXES: .f90 .o


.f90.o: 
	$(FC) -c $(FFLAGS) $<

# Include the dependency-list created by makedepf90 below 
include .depend


# target 'clean' for deleting object- *.mod- and other  
# unwanted files 
clean: 
	rm -f *.o *.mod core


# Create a dependency list using makedepf90.  All files  
# that needs to be compiled to build the program,  
# i.e all source files except include files, should  
# be given on the command line to makedepf90.   
# 
# The argument to the '-o' option will be the name of the 
# resulting program when running 'make', in this case  
# 'foobar' 
depend .depend: 
	makedepf90 -W -o Surface.x *.f90 > .depend
