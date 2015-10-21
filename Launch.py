#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#	Launch.py
#  
#	Launcher for the Willard-Chandler code.
#	
#	Functions: read input arguments, check, prepare input file for surface.F90 and launch it.

import sys,math,argparse


def main():
	
	# Very first thing: read in command line arguments and options
	parser = argparse.ArgumentParser(description="""Python code to analyze an interface trajectory and compute the Willard-Chandler dividing surface (wrapper for surface.F90)""")
	parser.add_argument('input_traj', action="store", help="Trajectory input file. MUST have the OP as fourth column.")
	parser.add_argument('out_files', nargs=2, action="store", help="Output files: wrapped traj, surface plot")
	parser.add_argument('-s','--stride',dest="stride", action="store", help="Stride for reading trajectory")
	parser.add_argument('-b','--box', nargs=3, type=float, action="store", dest="box_lengths", help="Simulation box size", required=True)
	parser.add_argument('--opref', nargs=2, type=float, action="store", dest="op_ref", help="The two reference values of the OP for the different phases. By default, the constant is set to the arithmetic average of these two values", required=True)
	parser.add_argument('--xi', type=float, action="store", dest="xi", help="Sigma value for Gaussian functions")

	args = parser.parse_args() # Run the parser and store arguments
	
	# Second: check the consistency of arguments
	#	'argparse' automatically deals with compulsory and optional arguments. We only need to check:
	#		1) OP reference values are non-zero
	#		2) the variance parameter for the Gaussian function is given, otherwise define automatically one
	#		3) stride? Default = 1 (analyze each frame)
	#
	# 1 check
	if ( any(args.op_ref[i] == 0.0e0 for i in range(len(args.op_ref))) ):
		print "OP reference values must be non-zero! Abort."
		sys.exit(1)
	# 2 check
	if ( args.xi is None ):
		xi = 0.5e0 # MUST BE CHECKED !!!
	else:
		xi = args.xi
	# 3 check
	if (args.stride is None):
		stride = 1
	else:
		stride = args.stride
	
	
	return 0

# This is what is run if the file is run as a script. Simple :-)
if __name__ == '__main__':
	main()

