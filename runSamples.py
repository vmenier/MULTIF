 #!/usr/bin/env python

import os, sys

from optparse import OptionParser
import multiprocessing

import multif
import time

def runSample_wrap(sample):
    sucess, val_out = sample.RunSample();
    return [sample.run_id, sucess, val_out];
    
def runSamplePostpro_wrap(sample):
    sucess, val_out = sample.RunSamplePostpro();
    return [sample.run_id, sucess, val_out];
    
def runSampleSkipAero_wrap(sample):
    sucess, val_out = sample.RunSampleSkipAero();
    return [sample.run_id, sucess, val_out];
    
def runSampleVisu_wrap(sample):
    return sample.RunSampleVisu();




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
    # parser.add_option("-n", "--ntasks", dest="nTasks", default=1,
    #                   help="number of tasks", metavar="NTASKS")    
    # parser.add_option("-c", "--cpus-per-task", dest="cpusPerTask",
    #                   default=1, help="cpus requested per task",
    #                   metavar="CPUS_PER_TASK")    
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
    parser.add_option("-o", "--output", dest="outfile", default="data.out",
                      help="Output file", metavar="OUTPUT")
                      
    parser.add_option("-g", "--postpro",
                      dest="postpro", default=False, action="store_true",
                      help="Run post-processing functions only?")
                     
    parser.add_option("-k", "--skipaero",
                      dest="skipaero", default=False, action="store_true",
                      help="Skip aero analysis?")
                           
    parser.add_option("-v", "--visu",
                      dest="visu", default=False, action="store_true",
                      help="Run visualization functions only?")
                          
    (options, args)=parser.parse_args()
    
    options.partitions     = int( options.partitions )
    # options.nTasks = int( options.nTasks )
    # options.cpusPerTask = int( options.cpusPerTask )    
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

    if options.postpro :
        sys.stdout.write("  -- Info : Running postprocessing functions only.\n\n");
        
    if options.skipaero :
        sys.stdout.write("  -- Info : Skipping aero analysis.\n\n");
    
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
        #samples_tab[-1].partitions      = options.partitions;
        samples_tab[-1].cpusPerTask     = options.partitions;
        samples_tab[-1].postpro         = options.postpro;
        samples_tab[-1].skipaero        = options.skipaero;
    
    #--- Create ./runs folder
    
    runs_dirNam = "runs_%s_%d" % (options.samples_filename, options.flevel); # wrapping folder containing all local run dirs
    
    ## Save it if it already exists
    #if os.path.isdir(runs_dirNam):
    #    try : 
    #        runs_dirNam_save = "%s_%s" % (runs_dirNam, time.strftime("%Y%m%d_%I:%M:%S"));
    #        os.rename(runs_dirNam, runs_dirNam_save);
    #        sys.stdout.write("Renamed previous %s folder into %s.\n" % (runs_dirNam,runs_dirNam_save));
    #    except:
    #        sys.stderr.write("  ## ERROR : Unable to save previous folder %s. Rename/delete it and try again.\n\n" % (runs_dirNam));
    #        sys.exit(0);
    #
    #os.mkdir(runs_dirNam);
    
    
    if not os.path.isdir(runs_dirNam):
        os.mkdir(runs_dirNam);
        
    #--- Create visu folder if needed
    visu_dirNam = "visu"
    if options.visu:
        if not os.path.isdir(visu_dirNam):
            os.mkdir(visu_dirNam);
    
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
        
        rEval = [];
    	for i in range(NbrRun):
            rEval.append(-1);
        
        for i in range(NbrRun):
            
            
            rEval[i] = runSample_wrap(samples_tab[i]);
            #if options.visu :
            #    rEval[i] = runSampleVisu_wrap(samples_tab[i]); 
            #elif options.postpro :
            #    rEval[i] = runSamplePostpro_wrap(samples_tab[i]);             
            #elif options.skipaero :
            #    rEval[i] = runSampleSkipAero_wrap(samples_tab[i]); 
            #else:
            #    rEval[i] = runSample_wrap(samples_tab[i]);
    else:
        mEval = [];
        rEval = [];
    	for i in range(NbrRun):
            mEval.append(-1);
            rEval.append(-1);
            
        for i in range(NbrRun):
            iRun = samples_tab[i].run_id;   
            
            if options.visu :
                mEval[i] = pool.apply_async(runSampleVisu_wrap,(samples_tab[i],));
            elif  options.postpro :
                mEval[i] = pool.apply_async(runSamplePostpro_wrap,(samples_tab[i],));
            elif  options.skipaero :
                mEval[i] = pool.apply_async(runSampleSkipAero_wrap,(samples_tab[i],));
            else :
                mEval[i] = pool.apply_async(runSample_wrap,(samples_tab[i],));
    	
        pool.close();
        pool.join();
        
        for i in range(len(mEval)):
            rEval[i] = mEval[i].get();
    
    for i in range(len(rEval)):
        print rEval[i];
    
    if not options.visu:
        
        # Write results to file
        f = open(options.outfile,'w');
        for i in range(len(rEval)):
            for j in range(len(rEval[i][2])-1):
                f.write('%0.16f, ' % rEval[i][2][j]);
            f.write('%0.16f\n' % rEval[i][2][-1]);
        f.close();
        sys.stdout.write("-- %s written with sample data.\n" % options.outfile);
    
        
        
# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
