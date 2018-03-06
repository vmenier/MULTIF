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
        
        #self.partitions = 1;
        self.nTasks = 1;
        self.cpusPerTask = 1;
        
    def RunSample(self):
        
        sys.stdout.write('-- Running sample %d \n' % self.run_id);
        
        redirect = True;
        
        run_id = self.run_id;
        
        #--- Go to root working dir
        
        os.chdir(self.working_rootdir);
        
        #--- Create run working dir
        
        runs_dirNam = "runs_%s_%d" % (self.samples_file, self.fidelity); # wrapping folder containing all local run dirs
        
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
            
            
            if redirect:
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
            config.OUTPUT_GRADIENTS= 'NO'
    	    nozzle = multif.nozzle.NozzleSetup(config, self.fidelity);
    	    #nozzle.partitions = int(self.partitions);
            nozzle.nTasks = int(self.nTasks);
            nozzle.cpusPerTask = int(self.cpusPerTask);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            # Hack to debug error handling: make it diverge
            #if run_id == 2:
            #    nozzle.mission.mach = 1e6;
            
            #--- Run analysis
            
            output = 'verbose';
            
            #nozzle.cfd.su2_max_iterations = 2; # hack
            
            if nozzle.method == 'NONIDEALNOZZLE' :
                multif.LOWF.Run(nozzle, output=output);
            elif nozzle.dim == '2D':
                multif.MEDIUMF.Run(nozzle, output=output);
            elif nozzle.dim == '3D':
                multif.HIGHF.Run(nozzle, output=output);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            success = True;
            
        except:
            if redirect:
                sys.stdout = sav_stdout;
                sys.stderr = sav_stderr;
            sys.stderr.write("## Error : Run %d failed.\n" % run_id);
            #raise;
            return success, val_out;
        
        if redirect:
            sys.stdout = sav_stdout;
            sys.stderr = sav_stderr;
        
        return success, val_out;


    def RunSamplePostpro(self):
        
        sys.stdout.write('-- Running sample %d postpro\n' % self.run_id);
        
        
        success = False;
        val_out = [False];
        
        run_id = self.run_id;
                        
        #--- Go to working dir
        
        runs_dirNam = "runs_%s_%d" % (self.samples_file, self.fidelity); # wrapping folder containing all local run dirs
        
        dirNam = "run_%d" % run_id; # local run dir
        locDir = os.path.join(self.working_rootdir,runs_dirNam,dirNam);
        
        if not os.path.isdir(locDir):
            sys.stderr.write ("## ERROR : %s does not exit. Skip.\n" % locDir);
            sys.exit(0);
        else:
            os.chdir(locDir);
            
            
        #--- Copy cfg file
        
        shutil.copyfile(os.path.join(self.working_rootdir,self.cfg_file),self.cfg_file);
        
        if not self.cfg_file:
            sys.stderr.write ("## ERROR run %d: configuration file %s does not exit. Skip.\n" % (run_id, self.cfg_file));
            return success, val_out;
        
        if not self.input_file:
            sys.stderr.write ("## ERROR run %d: input file %s does not exit. Skip.\n" % (run_id, self.input_file));
            return success, val_out;
            
        #--- Open log files
        
        stdout_hdl = open("postpro_%s" % self.stdout,'w'); # new targets
        stderr_hdl = open("postpro_%s" % self.stderr,'w');
        
        try: # run with redirected outputs
            
            sav_stdout, sys.stdout = sys.stdout, stdout_hdl; 
            sav_stderr, sys.stderr = sys.stderr, stderr_hdl; 
            
            #--- Setup nozzle data structure
            
    	    config = multif.SU2.io.Config(self.cfg_file);
            config.INPUT_DV_NAME = self.input_file;
    	    nozzle = multif.nozzle.NozzleSetup(config, self.fidelity);
    	    #nozzle.partitions = int(self.partitions);
            nozzle.nTasks = int(self.nTasks);
            nozzle.cpusPerTask = int(self.cpusPerTask);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            output = 'verbose';
            
            if nozzle.method == 'NONIDEALNOZZLE' :
                sys.stderr.write("## ERROR : Postprocessing only not available for 1D.\n");
                sys.exit(0);
            elif nozzle.dim == '2D':
                multif.MEDIUMF.Run(nozzle, output=output, postpro=1);
            elif nozzle.dim == '3D':
                multif.HIGHF.Run(nozzle, output=output, postpro=1);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            success = True;
            
        except:
            sys.stdout = sav_stdout;
            sys.stderr = sav_stderr;
            sys.stderr.write("## Error : Run %d failed.\n" % run_id);
            #raise;
            return success, val_out;
        
        sys.stdout = sav_stdout;
        sys.stderr = sav_stderr;
        
        return success, val_out;
        
    def RunSampleSkipAero(self):
        
        sys.stdout.write('-- Running sample %d skipaero\n' % self.run_id);
        
        run_id = self.run_id;
                        
        #--- Go to working dir
        
        runs_dirNam = "runs_%s_%d" % (self.samples_file, self.fidelity); # wrapping folder containing all local run dirs
        
        dirNam = "run_%d" % run_id; # local run dir
        locDir = os.path.join(self.working_rootdir,runs_dirNam,dirNam);
        
        
        if not os.path.isdir(locDir):
            sys.stderr.write ("## ERROR : %s does not exit. Skip.\n" % locDir);
            sys.exit(0);
        else:
            os.chdir(locDir);
            
            
        #--- Copy cfg file
        shutil.copyfile(os.path.join(self.working_rootdir,self.cfg_file),self.cfg_file);
        
        if not self.cfg_file:
            sys.stderr.write ("## ERROR run %d: configuration file %s does not exit. Skip.\n" % (run_id, self.cfg_file));
            sys.exit(0);
        
        if not self.input_file:
            sys.stderr.write ("## ERROR run %d: input file %s does not exit. Skip.\n" % (run_id, self.input_file));
            sys.exit(0);
            
        #--- Open log files
        
        stdout_hdl = open("skipaero_%s" % self.stdout,'w'); # new targets
        stderr_hdl = open("skipaero_%s" % self.stderr,'w');
        
        success = False;
        val_out = [False];
        
        try: # run with redirected outputs
            
            sav_stdout, sys.stdout = sys.stdout, stdout_hdl; 
            sav_stderr, sys.stderr = sys.stderr, stderr_hdl; 
            
            #--- Setup nozzle data structure
            
    	    config = multif.SU2.io.Config(self.cfg_file);
            config.INPUT_DV_NAME = self.input_file;
    	    nozzle = multif.nozzle.NozzleSetup(config, self.fidelity);
    	    #nozzle.partitions = int(self.partitions);
            nozzle.nTasks = int(self.nTasks);
            nozzle.cpusPerTask = int(self.cpusPerTask);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            output = 'verbose';
            
            if nozzle.method == 'NONIDEALNOZZLE' :
                sys.stderr.write("## ERROR : Skip aero only not available for 1D.\n");
                sys.exit(0);
            elif nozzle.dim == '2D':
                multif.MEDIUMF.Run(nozzle, output=output, skipAero=1);
            elif nozzle.dim == '3D':
                multif.HIGHF.Run(nozzle, output=output, skipAero=1);
            
            tag_out, val_out, gra_out, gratag_out = nozzle.GetOutputFunctions();
            
            success = True;
            
        except:
            sys.stdout = sav_stdout;
            sys.stderr = sav_stderr;
            sys.stderr.write("## Error : Run %d failed.\n" % run_id);
            #raise;
            return success, val_out;
        
        sys.stdout = sav_stdout;
        sys.stderr = sav_stderr;
        
        return success, val_out;
 
    
    def RunSampleVisu(self):
        sys.stdout.write('-- Running sample %d visualization\n' % self.run_id);
        run_id = self.run_id;
                        
        #--- Go to working dir
        
        runs_dirNam = "runs_%s_%d" % (self.samples_file, self.fidelity); # wrapping folder containing all local run dirs
        
        dirNam = "run_%d" % run_id; # local run dir
        locDir = os.path.join(self.working_rootdir,runs_dirNam,dirNam);
        
        if not os.path.isdir(locDir):
            sys.stderr.write ("## ERROR : %s does not exit. Skip.\n" % locDir);
            sys.exit(0);
        else:
            os.chdir(locDir);
        
        if not self.cfg_file:
            sys.stderr.write ("## ERROR run %d: configuration file %s does not exit. Skip.\n" % (run_id, self.cfg_file));
            sys.exit(0);
        
        if not self.input_file:
            sys.stderr.write ("## ERROR run %d: input file %s does not exit. Skip.\n" % (run_id, self.input_file));
            sys.exit(0);
            
        #--- Open log files
        
        stdout_hdl = open("visu_%s" % self.stdout,'w'); # new targets
        stderr_hdl = open("visu_%s" % self.stderr,'w');
        
        success = False;
        val_out = [False];
        
        #--- Setup nozzle data structure
        
        config = multif.SU2.io.Config(self.cfg_file);
        config.INPUT_DV_NAME = self.input_file;
        nozzle = multif.nozzle.NozzleSetup(config, self.fidelity);
        #nozzle.partitions = int(self.partitions);
        nozzle.nTasks = int(self.nTasks);
        nozzle.cpusPerTask = int(self.cpusPerTask);
        
        visu_dirNam = os.path.join(self.working_rootdir,"visu");
        visu_prefix = "run%d_" % run_id; 
        
        multif.visu.Viz_NozzleVisu(nozzle, visu_path=visu_dirNam, visu_prefix=visu_prefix);
        
        #try: # run with redirected outputs
        #    
        #    sav_stdout, sys.stdout = sys.stdout, stdout_hdl; 
        #    sav_stderr, sys.stderr = sys.stderr, stderr_hdl; 
        #    
        #    #--- Setup nozzle data structure
        #    
    	#    config = multif.SU2.io.Config(self.cfg_file);
        #    config.INPUT_DV_NAME = self.input_file;
    	#    nozzle = multif.nozzle.NozzleSetup(config, self.fidelity);
    	#    nozzle.partitions = int(self.partitions);
        #    
        #    visu_dirNam = os.path.join(self.working_rootdir,"visu");
        #    visu_prefix = "run%d_" % run_id; 
        #    
        #    multif.visu.Viz_NozzleVisu(nozzle, visu_path=visu_dirNam, visu_prefix=visu_prefix);
        #    
        #    success = True;
        #    
        #except:
        #    sys.stdout = sav_stdout;
        #    sys.stderr = sav_stderr;
        #    sys.stderr.write("## Error : Run %d failed.\n" % run_id);
        #    return success, val_out;
        #
        #sys.stdout = sav_stdout;
        #sys.stderr = sav_stderr;
        #
        return 1;
    
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
        
            
        
        
        
        
    