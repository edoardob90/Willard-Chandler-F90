.PHONY: all clean

F90=gfortran
FOPTS=-g -O2 
OBJECTS=mod_surf.o mod_interp.o inout.o calculations.o surface.o

all: surface.x

surface.x: $(OBJECTS)
	$(F90) $(FOPTS) -o surface.x $(OBJECTS)

%.o: %.f90 
	$(F90) $(FOPTS) -c $< 

clean: 
	rm -f $(OBJECTS)
