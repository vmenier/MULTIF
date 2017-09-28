import time, os, shutil, subprocess, datetime, sys
import numpy as np

import multif

class Sample:
    
    def __init__(self,run_id):
        
        self.run_id = run_id;        
        self.multif_dir = os.path.dirname(os.path.abspath(__file__));
        self.working_rootdir =  os.getcwd();
        
        self.samples_file = "";
        
        self.stdout = "stdout.job";
        self.stderr = "stderr.job";
                
        self.cfg_file   = "general.cfg";
        self.input_file = "inputDV.in";
        
        self.fidelity = 0;
        
        self.partitions = 1;
        
    def RunSample(self):
        
        sys.stdout.write('-- Running sample %d \n' % self.run_id);
        
        run_id = self.run_id;
        
        #--- Go to root working dir
        
        os.chdir(self.working_rootdir);
        
        #--- Create run working dir
        
        runs_dirNam = "runs"; # wrapping folder containing all local run dirs
        
        if not os.path.isdir(runs_dirNam):
            os.mkdir(runs_dirNam);
            os.chdir(runs_dirNam);
        else:
            os.chdir(runs_dirNam);
        
        dirNam = "run_%d" % run_id; # local run dir
        
        if os.path.isdir(dirNam):
            shutil.rmtree(dirNam);
        os.mkdir(dirNam);
        os.chdir(dirNam);
        
        #--- Open log files
        
        stdout_hdl = open(self.stdout,'w'); # new targets
        stderr_hdl = open(self.stderr,'w');
        
        success = False;
        val_out = [False];
        
        try: # run with redirected outputs
            
            sav_stdout, sys.stdout = sys.stdout, stdout_hdl; 
            sav_stderr, sys.stderr = sys.stderr, stderr_hdl; 
                        
            #--- Copy cfg file
            shutil.copyfile(os.path.join(self.working_rootdir,self.cfg_file),self.cfg_file);
            
            #--- Create DV file
            try:
                self.FormatDVFile();
            except:
                sys.exit(0);
                    
            #--- Setup nozzle data structure
    	    config = multif.SU2.io.Config(self.cfg_file);
            config.INPUT_DV_NAME = self.input_file;
    	    nozzle = multif.nozzle.NozzleSetup(config, self.fidelity);
    	    nozzle.partitions = int(self.partitions);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            # Hack to debug error handling: make it diverge
            #if run_id == 2:
            #    nozzle.mission.mach = 1e6;
            
            #--- Run analysis
            
            output = 'verbose';
            
            nozzle.cfd.su2_max_iterations = 300;
            
            if nozzle.method == 'NONIDEALNOZZLE' :
                multif.LOWF.Run(nozzle, output);
            elif nozzle.dim == '2D':
                multif.MEDIUMF.Run(nozzle, output);
            elif nozzle.dim == '3D':
                multif.HIGHF.Run(nozzle, output);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            success = True;
            
        except:
            sys.stdout = sav_stdout;
            sys.stderr = sav_stderr;
            sys.stderr.write("## Error : Run %d failed.\n" % run_id);
            return success, val_out;
        
        sys.stdout = sav_stdout;
        sys.stderr = sav_stderr;
        
        return success, val_out;
    
    def FormatDVFile(self):
        
        run_id = self.run_id;
        
        samples_filename = os.path.join(self.working_rootdir,self.samples_file);
        
        try:
            hdl = np.loadtxt(samples_filename);
        except:
            sys.stderr.write("  ## ERROR : Unable to open samples file %s. It might be invalid.\n" % (samples_filename));
            sys.exit(0);
        
        if run_id > len(hdl) or run_id < 0:
            sys.stderr.write("  ## ERROR : Invalid run_id=%d (%d samples in %s)\n" % (run_id, len(hdl), samples_filename));
            sys.exit(0);
        
        # --- Write 
        
        try:  
            fil = open(self.input_file, "w");
            for i in range(len(hdl[run_id])):
                fil.write("%.16le\n" % hdl[run_id][i]);
            fil.close(); 
        except:
            sys.stderr.write("  ## ERROR : run_id=%d : Unable to write DV file. \n" % (run_id));
            sys.exit(0);
        
        return;
        
            
        
        
        
        
    