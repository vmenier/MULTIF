#!/usr/bin/env python

import os, sys

from optparse import OptionParser
import multiprocessing

import multif


def runSample_wrap(sample):
    sucess, val_out = sample.RunSample();
    return [sample.run_id, sucess, val_out];
    
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
    parser.add_option("-f", "--file", dest="filename", default=False,
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS used to run each individual ", metavar="PARTITIONS")  
    parser.add_option("-p", "--poolpartitions", dest="poolpartitions", default=1,
                      help="number of partitions used for the whole set of runs", metavar="PARTITIONS")                   
    parser.add_option("-l", "--flevel", dest="flevel", default=-1,
                      help="fidelity level to run", metavar="FLEVEL")                  
    parser.add_option("-d", "--deform",
                      dest="deform", default=False,
                      help="Use mesh deformation?")
                      
    parser.add_option("-s", "--filesamples", dest="samples_filename", default=False,
                      help="Samples file", metavar="FILE")                      
                      
    parser.add_option("-b", "--beg", dest="run_beg", default=-1,
                      help="Id of first sample (from 0 to N-1, bounds are included)", metavar="FLEVEL")   
                      
    parser.add_option("-e", "--end", dest="run_end", default=-1,
                      help="Id of last sample (from 0 to N-1, bounds are included)", metavar="FLEVEL")   
                          
    (options, args)=parser.parse_args()
    
    options.partitions     = int( options.partitions )
    options.poolpartitions = int( options.poolpartitions )
    options.flevel         = int( options.flevel )
    
    if options.flevel < 0:
        sys.stderr.write("  ## ERROR : Please choose a fidelity level to run (option -l or --flevel)\n\n");
        sys.exit(0);
    
    if options.filename != False:
        if not os.path.isfile(options.filename) :
            sys.stderr.write("  ## ERROR : could not find configuration file %s\n\ns" % options.filename);
            sys.exit(0);        
    else :         
        sys.stderr.write("  ## ERROR : Please provide a configuration file. (option -f)\n");
        sys.exit(0);    
            
    if options.samples_filename != False:
        if not os.path.isfile(options.samples_filename) :
            sys.stderr.write("  ## ERROR : could not find configuration file %s\n\ns" % options.samples_filename);
            sys.exit(0);
    else :         
        sys.stderr.write("  ## ERROR : Please provide a file containing samples. (option --filesamples)\n");
        sys.exit(0);    

    #--- Check number of samples in file
    
    Nbs = sum(1 for line in open(options.samples_filename))
    
    sys.stdout.write("-- %s opened.\n   Number of samples: %d\n" % (options.samples_filename,Nbs));
    
    if Nbs <= 0 :
        sys.stderr.write("  ## ERROR : Sample file is empty.\n\n");
        sys.exit(0);           
        
    #--- Check parallel options
    
    sys.stdout.write("\nChecking parallel options:\n");
    
    NbrCpu = multiprocessing.cpu_count();
    sys.stdout.write("\t%d cpus available total.\n" % NbrCpu);
    sys.stdout.write("\t%d parallel pool instances max requested.\n" % options.poolpartitions);
    sys.stdout.write("\t%d cpus used for each MULTI-F run.\n" % options.partitions);
    
    MaxCpuUse = options.poolpartitions*options.partitions;
    if MaxCpuUse > NbrCpu :
        sys.stderr.write("\n  ## ERROR parallel options: %d processes will be used simutaneously, but only %d cpus are available on this machine.\n\n" \
        % (MaxCpuUse, NbrCpu));
        sys.exit(0);
    
    #--- Check run id bounds
    
    run_beg = int( options.run_beg );
    run_end = int( options.run_end ); 
       
    if run_beg >= 0 and run_end >= 0:
        if run_beg > run_end or run_beg > Nbs-1 or run_end > Nbs-1:
            sys.stderr.write("  ## ERROR : Invalid sample bounds: [%d, %d] (%d samples in file)\n\ns" % (run_beg,run_end, Nbs));
            sys.exit(0);                   
    elif run_beg < 0 and run_end < 0:
        #  No bounds provided, we run the whole sample file
        sys.stdout.write("   ## WARNING : No sample bounds were provided. All the samples will be run.\n");
        run_beg = 0;
        run_end = Nbs-1;
    else :
        # Only one valid bound provided? This is a bug
        sys.stderr.write("  ## ERROR : Invalid sample bounds: [%d, %d]\n\ns" % (run_beg,run_end));
        sys.exit(0);
    
    NbrRun = run_end-run_beg+1;
    
    sys.stdout.write("-- Info : Running samples %d to %d (%d run(s) total).\n\n" % (run_beg, run_end, NbrRun));
    
    #--- Fill samples data structure    
    
    samples_tab = [];
    
    for iRun in range(run_beg, run_end+1):
        
        samples_tab.append(multif.samples.Sample(iRun));
                
        samples_tab[-1].multif_dir      = os.path.dirname(os.path.abspath(__file__));
        samples_tab[-1].working_rootdir = os.getcwd();
        samples_tab[-1].samples_file    = options.samples_filename;
        samples_tab[-1].stdout          = 'stdout.job'; 
        samples_tab[-1].stderr          = 'stderr.job';      
        samples_tab[-1].cfg_file        = options.filename;
        samples_tab[-1].input_file      = "inputDV.in";
        samples_tab[-1].fidelity        = options.flevel;
        samples_tab[-1].partitions      = options.partitions;
    
    #--- Create ./runs folder
    
    runs_dirNam = "runs"; # wrapping folder containing all local run dirs
    
    if not os.path.isdir(runs_dirNam):
        os.mkdir(runs_dirNam);
    
    #--- Start python's multiprocessing pool
    
    if options.poolpartitions > 1:
        sys.stdout.write("Initializing multiprocessing pool using %d processors.\n" % options.poolpartitions);
        try:
    	    pool = multiprocessing.Pool(processes=options.poolpartitions, maxtasksperchild=1);
            #pool = multiprocessing.Pool(processes=options.poolpartitions);
        except:
            sys.stderr.write("  ## ERROR : Parallel initialization failed.\n\n");
            sys.exit(0);
    
    if options.poolpartitions == 1:
        for i in range(NbrRun):
            samples_tab[i].RunSample();
    else:
        mEval = [];
        rEval = [];
    	for i in range(NbrRun):
            mEval.append(-1);
            rEval.append(-1);
            
        for i in range(NbrRun):
            iRun = samples_tab[i].run_id;    
            mEval[i] = pool.apply_async(runSample_wrap,(samples_tab[i],));
    	
        pool.close();
        pool.join();
    
        for i in range(len(mEval)):
            rEval[i] = mEval[i].get();
        
        for i in range(len(mEval)):
            print rEval[i];
        
        
        
# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
