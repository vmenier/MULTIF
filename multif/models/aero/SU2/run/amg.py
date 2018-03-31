import os, sys, shutil, copy, time

from .. import io   as su2io
from .. import amginria as su2amg
from interface import CFD as SU2_CFD
from ..amginria import _amgio as amgio

def amg ( config , kind='' ):
    
    sys.stdout.write("Run Anisotropic Mesh Adaptation\n");
    
    #--- Check options
    
    adap_options = ['ADAP_SIZES', 'ADAP_SUBITE', 'ADAP_SENSOR', \
    'ADAP_BACK', 'ADAP_HMAX', 'ADAP_HMIN', 'ADAP_HGRAD', 'ADAP_RESIDUAL_REDUCTION', 'ADAP_EXT_ITER']
    required_options = ['ADAP_SIZES', 'ADAP_SUBITE', \
    'ADAP_SENSOR', 'MESH_FILENAME', 'RESTART_SOL', 'MESH_OUT_FILENAME']
    
    if not all (opt in config for opt in required_options):
        err = '\n\n## ERROR : Missing options: \n'
        for opt in required_options:
            if not opt in config:
                err += opt + '\n'
        raise RuntimeError , err
    
    # Print adap options
    sys.stdout.write(su2amg.print_adap_options(config, adap_options));
    
    #--- How many iterative loops? Using what prescribed mesh sizes? 
        
    mesh_sizes   = su2amg.get_mesh_sizes(config);
    sub_iter     = su2amg.get_sub_iterations(config);
    
    # solver iterations/ residual reduction param for each size level
    adap_ext_iter = su2amg.get_ext_iter(config);
    adap_res = su2amg.get_residual_reduction(config);

    adap_sensor = config.ADAP_SENSOR;
    sensor_avail = ['MACH', 'PRES', 'MACH_PRES'];
    
    if adap_sensor not in sensor_avail:
        raise RuntimeError , 'Unknown adaptation sensor (ADAP_SENSOR option)\n'
        
    if len(mesh_sizes) != len(sub_iter):
        raise RuntimeError , 'Inconsistent number of mesh sizes and sub-iterations'
    
    #--- Change current directory
    
    warn = False;
    adap_dir = './ADAP';
    cwd = os.getcwd();
        
    if os.path.exists(adap_dir):
        sys.stdout.write('./ADAP exists. Removing old mesh adaptation in 10s.\n')
        sys.stdout.flush();
        if warn : time.sleep(10);
        shutil.rmtree(adap_dir);
    
    os.makedirs(adap_dir);
    os.chdir(adap_dir);
    sys.stdout.write('The %s folder was deleted\n' % adap_dir);
    
    os.symlink(os.path.join(cwd, config.MESH_FILENAME), config.MESH_FILENAME);
    os.symlink(os.path.join(cwd, config.SOLUTION_FLOW_FILENAME), config.SOLUTION_FLOW_FILENAME);
    
    
    #--- Compute initial solution if needed, else link current files
    
    config_cfd = copy.deepcopy(config);
    for opt in adap_options:
        config_cfd.pop(opt, None)
    config_cfd.LOW_MEMORY_OUTPUT = "YES";
    
    # HACK
    #config_cfd.EXT_ITER = 2;
    
    current_mesh     = "Initial_mesh";
    current_solution = "Initial_solution";
        
    if config['RESTART_SOL'] == 'NO':
        
        stdout_hdl = open('ini.stdout','w'); # new targets
        stderr_hdl = open('ini.stderr','w');
        
        success = False;
        val_out = [False];
        
        sys.stdout.write('Running initial solution. Log file : %s\n' % 'ini.stdout')
        
        try: # run with redirected outputs
            
            sav_stdout, sys.stdout = sys.stdout, stdout_hdl; 
            sav_stderr, sys.stderr = sys.stderr, stderr_hdl;
        
            current_mesh     = config['MESH_FILENAME'];
            current_solution = "ini_restart_flow.dat"
            
            config_cfd.CONV_FILENAME         = "ini_history"
            config_cfd.RESTART_FLOW_FILENAME = current_solution;
            
            SU2_CFD(config_cfd);
            
        except:
            sys.stdout = sav_stdout;
            sys.stderr = sav_stderr;
            raise;
        
        sys.stdout = sav_stdout;
        sys.stderr = sav_stderr;
        
    else:
        required_options=['SOLUTION_FLOW_FILENAME']
        if not all (opt in config for opt in required_options):
            err = '\n\n## ERROR : RESTART_SOL is set to YES, but the solution is missing:\n'
            for opt in required_options:
                if not opt in config:
                    err += opt + '\n'
            raise RuntimeError , err
        
        current_mesh     = config['MESH_FILENAME'];
        current_solution = config['SOLUTION_FLOW_FILENAME'];
        
    #--- Check existence of initial mesh, solution
    
    required_files = [current_mesh,current_solution];
    
    if not all (os.path.exists(fil) for fil in required_files):
        err = '\n\n## ERROR : Can\'t find:\n'
        for fil in required_files:
            if not os.path.exists(fil):
                err += fil + '\n'
        raise RuntimeError , err
    
    #--- Start looping
    
    # Get mesh dimension
    dim = su2amg.get_su2_dim(current_mesh);
    if ( dim != 2 and dim != 3 ):
        raise RuntimeError , "Wrong dimension number\n"
    
    #--- AMG parameters
    
    config_amg = dict();
    
    config_amg['hgrad']       = float(config['ADAP_HGRAD'])
    config_amg['hmax']        = float(config['ADAP_HMAX'])
    config_amg['hmin']        = float(config['ADAP_HMIN'])
    config_amg['mesh_in']     = 'current.meshb'
    config_amg['mesh_out']    = 'current.new.meshb';
    config_amg['metric_in']   = ''
    config_amg['sol_in']      = 'current_sensor.solb'
    config_amg['itp_sol_in']  = 'current.solb'
    config_amg['adap_source'] = ''
    
    if 'ADAP_BACK' in config:
        config_amg['adap_back'] = config['ADAP_BACK'];
    else:
        config_amg['adap_back'] = config['MESH_FILENAME'];
        
    back_name, back_extension = os.path.splitext(config_amg['adap_back']);
    
    if not os.path.exists(config_amg['adap_back']):
        raise RuntimeError , "\n\n##ERROR : Can't find back mesh: %s.\n\n" % config_amg['adap_back']
    
    if back_extension == ".su2":
        amgio.py_ConvertSU2toInria(config_amg['adap_back'], "", "amg_back")
        config_amg['adap_back'] = "amg_back.meshb"
        
    
    if 'ADAP_SOURCE' in config:
        config_amg['adap_source'] = config['ADAP_SOURCE'];
    
    global_iter = 0;
    
    for iSiz in range(len(mesh_sizes)):
        
        mesh_size   = int(mesh_sizes[iSiz]);
        nSub        = int(sub_iter[iSiz]);
        
        for iSub in range(nSub):
            
            print "Global iter %d : Size %d, sub_ite %d" % (global_iter, mesh_size, iSub)
            
            #--- Convert current mesh/solution to inria format
            
            amgio.py_ConvertSU2toInria(current_mesh, current_solution, "current")
            
            if not os.path.exists("current.solb"):
                raise RuntimeError , "\n##ERROR : Can't find solution.\n"
            if not os.path.exists("current.meshb"):
                raise RuntimeError , "\n##ERROR : Can't find mesh.\n"
            
            #--- Get sensor
            
            amgio.py_SplitSolution(current_solution, dim, "current", adap_sensor);
            
            if not os.path.exists("current_sensor.solb"):
                raise RuntimeError , "\n##ERROR : Can't find adap sensor.\n"
                        
            #--- Run amg
            
            config_amg['size']        = mesh_size;
            config_amg['amg_log']     = 'ite%d.amg.stdout' % (global_iter);
            
            sys.stdout.write("Running amg. Log : %s\n" % config_amg['amg_log']);
            
            if os.path.exists("current.itp.solb"):
                os.remove("current.itp.solb");
                        
            try :
                su2amg.amg_call(config_amg);
            except:
                raise RuntimeError , "\n##ERROR : Call to AMG failed.\n"
            
            if not os.path.exists(config_amg['mesh_out']):
                raise RuntimeError , "\n##ERROR : Mesh adaptation failed.\n"
            
            if not os.path.exists("current.itp.solb"):
                raise RuntimeError , "\n##ERROR AMG: Solution interpolation failed.\n"            
            
            #--- Convert output

            current_mesh = "ite%d.su2" % global_iter
            current_solution = "ite%d.dat" % global_iter            
            amgio.py_ConvertInriatoSU2(config_amg['mesh_out'], "current.itp.solb", "ite%d" % global_iter)
                        
            if not os.path.exists(current_mesh) or not os.path.exists(current_solution) :
                raise RuntimeError , "\n##ERROR : Conversion to SU2 failed.\n"
                                    
            #--- Run su2
            
            log = 'ite%d.SU2'%global_iter;
            stdout_hdl = open('%sstdout'%log,'w'); # new targets
            stderr_hdl = open('%sstderr'%log,'w');
            
            success = False;
            val_out = [False];
            
            sys.stdout.write('Running SU2_CFD. Log file : %s.std[out/err]\n' % log)
        
            try: # run with redirected outputs
            
                sav_stdout, sys.stdout = sys.stdout, stdout_hdl; 
                sav_stderr, sys.stderr = sys.stderr, stderr_hdl;
                
                current_solution_ini = "ite%d_ini.dat" % global_iter
                os.rename(current_solution, current_solution_ini);
                
                config_cfd.MESH_FILENAME          = current_mesh;
                config_cfd.CONV_FILENAME          = "ite%d_history" % global_iter
                config_cfd.SOLUTION_FLOW_FILENAME = current_solution_ini;
                config_cfd.RESTART_FLOW_FILENAME  = current_solution;
                
                config_cfd.RESIDUAL_REDUCTION = float(adap_res[iSiz]);
                config_cfd.EXT_ITER = int(adap_ext_iter[iSiz]);
                
                SU2_CFD(config_cfd);
                
                if not os.path.exists(current_solution) :
                    raise RuntimeError , "\n##ERROR : SU2_CFD Failed.\n"
            
            except:
                sys.stdout = sav_stdout;
                sys.stderr = sav_stderr;
                raise;
            
            sys.stdout = sav_stdout;
            sys.stderr = sav_stderr;
            
            to_remove = ["current.itp.solb", config_amg['mesh_in'], config_amg['mesh_out'], config_amg['sol_in'],config_amg['itp_sol_in']];
            for fil in to_remove:
                if os.path.exists(fil) : os.remove(fil);
            
            global_iter += 1;
    
    os.rename(current_solution,os.path.join(cwd,config.RESTART_FLOW_FILENAME));
    os.rename(current_mesh,os.path.join(cwd,config.MESH_OUT_FILENAME));
    
    sys.stdout.write("\nMesh adaptation successfully ended. Results files:\n");
    sys.stdout.write("%s\n%s\n\n" % (config.MESH_OUT_FILENAME,config.RESTART_FLOW_FILENAME));
    