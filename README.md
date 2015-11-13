# Willard-Chandler-F90
This code has two main purposes:
  
  1- Given a MD trajectory with an interface between two phases (according to some kind of local order parameter), it computes and draws the location of the Willard-Chandler instataneous dividing surface.
  
  2- (STILL IN DEV) It calculates the "interface profile function" (i.e. the height of the surface in each point of the grid) and Fourier-trasforms it.
  
  From the average of the Fourier coefficients, within the framework of the Capillary Fluctuation Method, one can then get the stiffness of the interface and the interfacial free energy for that particular surface. This is accomplished by means of a fitting procedure of 1/A(k)^2 vs k^2, where k are the wave vector of the Fourier decomposition.


There are a couple of LIMITATIONS:
  
  i) one should know "a priori" the reference value of the order parameter, which recognizes the system in phase A and phase B.
  
  ii) the codes DOES NOT compute the fit of the Fourier coefficients (at least, up to now).



VERY BRIEF USAGE INSTRUCTIONS:
  
  The code can be compiled directly with a simple "make" command in the main directory. You must have installed "makedepf90", which builds for you all the dependecies of the code. BE AWARE that the Makefile is configured to use the Intel(R) Fortran Compiler, which gives the best performances on Intel CPU.
  If one wants to use the GNU Fortran compiler:
    - modify the compiler variable: FC=gfortran
    - modify the flags: FFLAGS=-03
  
  The only external tool needed is the FFTW library, freely available and installable/compilable. The Makefile will link the executable by looking up FFTW libraries in the standard location. If they are in a particular path, the user should add the option "-L/path/to/fftw3" to the LDFLAGS variable.
