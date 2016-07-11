#!/usr/bin/env python

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import SU2
import util

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
	
	# Command Line Options
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
	                  help="read config from FILE", metavar="FILE")
	parser.add_option("-n", "--partitions", dest="partitions", default=0,
	                  help="number of PARTITIONS", metavar="PARTITIONS")
	
	parser.add_option("-l", "--flevel", dest="flevel", default=-1,
	                  help="fidelity level to run", metavar="FLEVEL")	
	
	(options, args)=parser.parse_args()
	
	
	options.partitions = int( options.partitions )
	options.flevel     = int( options.flevel )
	
	
	if options.flevel < 0 :
		print "  ## ERROR : Please choose a fidelity level to run (option -l or --flevel)"
		sys.exit(1);
	
	# Check and print options
	#CheckOptimizationOptions ( options.filename, options.flevel );
	
	util.options.CheckOptimizationOptions( options.filename, options.flevel );
	
	x = [1,2,3,4]
	y = [2,2,2,2]
	util.CFD.NozzleGeoFile(x, y);
	
	
# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
