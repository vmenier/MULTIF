#!/usr/bin/env python

import os, sys

from optparse import OptionParser

import multif
from multif import SU2

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
	
	output = 'verbose';
	
	if output == 'verbose':         
	    sys.stdout.write('-' * 90);
	    sys.stdout.write('\n');
	    sys.stdout.write('\t __  __ _   _ _  _____ ___         ___ \t\n') ;
	    sys.stdout.write('\t|  \/  | | | | ||_   _|_ _|  ___  | __|\t\n');
	    sys.stdout.write('\t| |\/| | |_| | |__| |  | |  |___| | _| \t\t Dev. : R. Fenrich & V. Menier\n');
	    sys.stdout.write('\t|_|  |_|\___/|____|_| |___|       |_|  \t\t        Stanford University\n\n');
	    sys.stdout.write('-' * 90);
	    sys.stdout.write('\n\n');
	
	# Command Line Options
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
	                  help="read config from FILE", metavar="FILE")
	parser.add_option("-n", "--partitions", dest="partitions", default=1,
	                  help="number of PARTITIONS", metavar="PARTITIONS")
	
	parser.add_option("-l", "--flevel", dest="flevel", default=0,
	                  help="fidelity level to run", metavar="FLEVEL")   
	
	parser.add_option("-v", "--verif", dest="verification", default=False,
	                  help="Run verification test case?")
	                  
	parser.add_option("-d", "--deform",
	                  dest="deform", default=False,
	                  help="Use mesh deformation?")
	
	parser.add_option("-s", "--sample",
	                  dest="samples_filename", default=False,
	                  help="Run samples from file?")
	
	parser.add_option("", "--from",
	                  dest="samples_from", default=False,
	                  help="ID of first sample to run?")
	
	parser.add_option("", "--to",
	                  dest="samples_to", default=False,
	                  help="ID of last sample to run?")
	
	parser.add_option("", "--onebyone",
	                  dest="onebyone", default=False,
	                  help="Run samples one by one?")
	
	(options, args)=parser.parse_args()
	
	options.partitions = int( options.partitions )
	options.flevel     = int( options.flevel )
	
	options.samples_from = int(options.samples_from)
	options.samples_to   = int(options.samples_to)
	
	if options.flevel < 0:
	    sys.stderr.write("  ## ERROR : Please choose a fidelity level to run (option -l or --flevel)");
	    sys.exit(0);
	
	if int(options.verification):
		config = multif.nozzle.SetupVerificationConfig();
		options.flevel = 0;		
	else:		
		if not os.path.isfile(options.filename) :
		    sys.stderr.write("  ## ERROR : could not find configuration file %s\n\ns" % options.filename);
		    sys.exit(0);		
		config = SU2.io.Config(options.filename)
		
	nozzle = multif.nozzle.NozzleSetup(config, options.flevel, output);
	nozzle.verification = int(options.verification);
	nozzle.meshDeformationFlag = int(options.deform);
	nozzle.partitions = int(options.partitions);
	
	nozzle.onebyone = options.onebyone;
	
	if not options.samples_filename:
		sys.stderr.write("  ## ERROR : No sample file name was provided (option --sample).\n");
		sys.exit(0);
	
	
	if nozzle.NbrDVTot > 0 :
		nozzle.UpdateDV(output);
		nozzle.SetupWall(output);
	
	#multif.samples.RunOneSample(777, nozzle);
	multif.samples.RunSamples(nozzle, options.samples_filename, options.samples_from, options.samples_to);	
	

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------


if __name__ == '__main__':
    main()