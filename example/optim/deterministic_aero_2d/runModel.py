#!/usr/bin/env python

import os, sys
import subprocess 

from optparse import OptionParser

import multif

if __name__ == '__main__':

    output = 'verbose'
    
    if output == 'verbose':         
        print('-' * 90)
        print('')
        print('\t __  __ _   _ _  _____ ___         ___ \t') 
        print('\t|  \/  | | | | ||_   _|_ _|  ___  | __|\t')
        print('\t| |\/| | |_| | |__| |  | |  |___| | _| \t\t Dev. : R. Fenrich & V. Menier')
        print('\t|_|  |_|\___/|____|_| |___|       |_|  \t\t        Stanford University\n')
        print('-' * 90)
        print('\n')
    
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
    parser.add_option("-g", "--postpro", dest="postpro", default=False, 
                      action="store_true", 
                      help="Run post-processing functions only?")
    parser.add_option("-s", "--skipaero", dest="skipaero", default=False, 
                      action="store_true", help="Skip aero analysis")
    
    (options, args)=parser.parse_args()
    
    options.nTasks = int( options.nTasks )
    options.cpusPerTask = int( options.cpusPerTask )
    options.flevel     = int( options.flevel )
    
    if options.flevel < 0:
        raise RuntimeError("Please choose a fidelity level to run (option -l or --flevel)")
    
    if not os.path.isfile(options.filename):
        raise RuntimeError("Could not find configuration file %s" % options.filename)       
    
    config = multif.models.aero.SU2.io.Config(options.filename)        
    nozzle = multif.nozzle.nozzle.NozzleSetup(config, options.flevel, output)
    nozzle.nTasks = int(options.nTasks)
    nozzle.cpusPerTask = int(options.cpusPerTask)

    # --- Check parallel computing specifications
    ncores = nozzle.nTasks*nozzle.cpusPerTask # requested by user
    if output == 'verbose':
        print("User has requested %i total cores for calculation." % ncores)
            
    postpro = 0
    if options.postpro:
        postpro = 1
        
    skipaero = 0
    if options.skipaero:
        skipaero = 1
    
    # --- Run analysis
    multif.analysis.run(nozzle, output=output, post_process_only=postpro,
        skip_aero=skipaero)
    
    # --- Print warning in case the wrong SU2 version was run
    # if nozzle.method != 'NONIDEALNOZZLE' and nozzle.cfd.su2_version != 'OK':
    #     print('')
    #     print('#' * 90)
    #     print('WARNING: You are not using the right version of SU2. This may have caused robustness issues.')
    #     print('#' * 90)
    #     print('')
        
    
