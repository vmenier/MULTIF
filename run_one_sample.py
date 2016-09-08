#!/usr/bin/env python

import os, time, sys, shutil, copy



from optparse import OptionParser
import textwrap

import multif
from multif import SU2



# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
	 
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
	parser.add_option("-n", "--partitions", dest="partitions", default=0,
	                  help="number of PARTITIONS", metavar="PARTITIONS")
	
	parser.add_option("-l", "--flevel", dest="flevel", default=-1,
	                  help="fidelity level to run", metavar="FLEVEL")	
	
	parser.add_option("-i", "--id", dest="id_sample", default=-1,
	                  help="ID of the sample to be run", metavar="IDSAMPLE")
	
	parser.add_option("-s", "--sample", dest="sample_filename", default="TPMC_10k.txt",
	                  help="name of the sample file", metavar="samplefilename")
	
	
	
	(options, args)=parser.parse_args()
	
	options.partitions = int( options.partitions )
	options.flevel     = int( options.flevel )
	
	options.id_sample = int ( options.id_sample )
	
	if options.flevel < 0 :
		sys.stderr.write("  ## ERROR : Please choose a fidelity level to run (option -l or --flevel)\n\n");
		sys.exit(0);
		
	if options.id_sample <= 0 :
		sys.stderr.write("  ## ERROR : Please choose a sample to run\n\n");
		sys.exit(0);
	
		
	nozzle = multif.nozzle.NozzleSetup( options.filename, options.flevel, 'quiet');
	
	### --- Read sample file and get DV
	
	#options.sample_filename = 'TPMC_10k.txt';
	
	try:
		fil = open(options.sample_filename, 'r');
	except:
		sys.stderr.write("  ## ERROR : Could not open %s\n" % options.sample_filename);
		sys.exit(0);
	
	for ilin, line in enumerate(fil):
		if ilin == options.id_sample:
			break;
	
	#print line.split(",");
	dv_list = line.split(",") ;
	
	#for i in range(len(dv_list)):
	#	print "%d : %lf" % (i, float(dv_list[i]));
	
	NbrDV = len(dv_list);
	
	if NbrDV != len(nozzle.DV_List):
		sys.stderr.write("  ## ERROR : Inconsistent numbers of DV\n");
		sys.exit(0);
		
	for i in range(NbrDV):
		#print "%d : %lf -> %lf" % (i, nozzle.DV_List[i], float(dv_list[i]));
		nozzle.DV_List[i] = float(dv_list[i]);

	
	config = SU2.io.Config(options.filename);
	
	print nozzle.coefs;
	
	### --- Update DV + compute inner wall
	
	nozzle.UpdateDV(config,'verbose');
	
	nozzle.SetupInnerWall(config);

	nozzle.runDir = "run_%d" % options.id_sample;
	
	newpath = nozzle.runDir; 
	print newpath
	if not os.path.exists(newpath):
	    os.makedirs(newpath)
	multif.MEDIUMF.Run(nozzle);
		
	#nozzle = multif.nozzle.NozzleSetup( options.filename, options.flevel );
	#
	#### HACK
	##multif.MEDIUMF.AEROSPostProcessing(nozzle);
	##sys.exit(1);
	#
	#if nozzle.method == 'NONIDEALNOZZLE' :
	#	multif.LOWF.Run(nozzle);
	#elif nozzle.method == 'EULER' or nozzle.method == 'RANS':
	#	multif.MEDIUMF.Run(nozzle);
	#	
	## --- Output functions 
	#
	##nozzle.WriteOutputFunctions_Plain ();
	#nozzle.WriteOutputFunctions_Dakota ();
	#
	#sys.stdout.write('\n');
	#
	## --- Print warning in case the wrong SU2 version was run
	#if nozzle.method != 'NONIDEALNOZZLE' and nozzle.SU2Version != 'OK':
	#	sys.stdout.write('\n');
	#	sys.stdout.write('#' * 90);
	#	sys.stdout.write('\n  ## WARNING : You are not using the right version of SU2. This may have caused robustness issues.\n');
	#	sys.stdout.write('#' * 90);
	#	sys.stdout.write('\n\n');
	#
# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
