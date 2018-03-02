#!/usr/bin/env python

import os, sys

from optparse import OptionParser

import multif

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
    parser.add_option("-n", "--ntasks", dest="nTasks", default=1,
                      help="number of tasks", metavar="NTASKS")    
    parser.add_option("-c", "--cpus-per-task", dest="cpusPerTask",
                      default=1, help="cpus requested per task",
                      metavar="CPUS_PER_TASK")
    parser.add_option("-l", "--flevel", dest="flevel", default=0,
                      help="fidelity level to run", metavar="FLEVEL")                  
    parser.add_option("-d", "--deform",
                      dest="deform", default=False,
                      help="Use mesh deformation?")
                      
    parser.add_option("-g", "--postpro",dest="postpro", default=False, action="store_true",help="Run post-processing functions only?")
    
    parser.add_option("-s", "--skipaero",dest="skipaero", default=False, action="store_true",help="Skip aero analysis")
    
    (options, args)=parser.parse_args()
    
    options.nTasks = int( options.nTasks )
    options.cpusPerTask = int( options.cpusPerTask )
    options.flevel     = int( options.flevel )
    
    if options.flevel < 0:
        sys.stderr.write("  ## ERROR : Please choose a fidelity level to run (option -l or --flevel)");
        sys.exit(0);
    
    if not os.path.isfile(options.filename) :
        sys.stderr.write("  ## ERROR : could not find configuration file %s\n\ns" % options.filename);
        sys.exit(0);        
    
    config = multif.SU2.io.Config(options.filename)        
    nozzle = multif.nozzle.NozzleSetup(config, options.flevel, output);
    nozzle.nTasks = int(options.nTasks);
    nozzle.cpusPerTask = int(options.cpusPerTask);

    # --- Check parallel computing specifications
    ncores = nozzle.nTasks*nozzle.cpusPerTask # requested by user
    if output == 'verbose':
        sys.stdout.write("User has requested %i total cores for calculation.\n" % ncores);
            
    postpro = 0;
    if options.postpro :
        postpro = 1;
        
    skipaero = 0;
    if options.skipaero :
        skipaero = 1;
        
    if nozzle.method == 'NONIDEALNOZZLE' :
        multif.LOWF.Run(nozzle, output=output);
    elif nozzle.dim == '2D':
        multif.MEDIUMF.Run(nozzle, output=output, postpro=postpro);
    elif nozzle.dim == '3D':
        multif.HIGHF.Run(nozzle,output=output, postpro=postpro, skipAero=skipaero);
    
    # --- Print warning in case the wrong SU2 version was run
    if nozzle.method != 'NONIDEALNOZZLE' and nozzle.cfd.su2_version != 'OK':
        sys.stdout.write('\n');
        sys.stdout.write('#' * 90);
        sys.stdout.write('\n  ## WARNING : You are not using the right version of SU2. This may have caused robustness issues.\n');
        sys.stdout.write('#' * 90);
        sys.stdout.write('\n\n');

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
