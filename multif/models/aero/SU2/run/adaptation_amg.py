import os, sys, shutil, copy, time



from .. import io   as su2io
from .. import mesh as su2mesh

import interface as su2run




def amg_Inisol (config_adap, config):
    
    # TO DO HERE : DEAL WITH MESH FILE FORMAT
    #   If mesh format = inria and a solution is provided, continue
    #   If mesh format = SU2 and a solution is provided, call amg_Inisol to convert to the Inria format
    #   If a solution is not provided, call amg_Inisol which converts if needed
    
    rootDir    = config_adap.rootDir;
    adapDir    = config_adap.adapDir;
    partitions = config_adap.partitions;
    
    # --- Outputs: current mesh, restart, and sensor
    
    restart_name_dest  = "%s/current_restart.solb" % adapDir
    sensor_name_dest   = "%s/current_sensor.solb" % adapDir;  
    mesh_name_dest     = "%s/current.meshb" % adapDir;
    
    if config['ADAP_RESTART'] == 'NO':
        config_ini = copy.deepcopy(config);
                
        #state = su2io.State()
        #state.find_files(config_ini)

        # Set the number of partitions for parallel computations
        config_ini.NUMBER_PART = partitions
        config_ini.CONSOLE = 'QUIET'

        config_name = config_ini._filename;
        config_name = "adap.cfg"
        print "CONFIG_NAME %s" % config_name
        #sys.exit(1)
        
        folder_inisol='initial_solution';
        #pull="../%s" % config_name;
        pull=""
        link="../%s" % config_ini['MESH_FILENAME'];

        config_ini.RESTART_FLOW_FILENAME = 'initial_sol.solb';
        config_ini.MESH_FORMAT = 'INRIA';
        config_ini.RESTART_SOL = 'NO';
        
        #config_ini.EXT_ITER = 1;
        
        jobNam = "SU2_ini.job"
        config_ini.OUTPUT_LOG = jobNam;
        print " Running SU2 on initial mesh \n Log: %s\n" % (jobNam) ;
                
        with su2io.redirect_folder(folder_inisol,pull,link):
            info = su2run.CFD(config_ini)
            shutil.copyfile (config_ini.RESTART_FLOW_FILENAME, restart_name_dest);  # restart solution
            shutil.copyfile ('mach.solb', sensor_name_dest );  # the solution used for error estimation (migh be different than the restart, ex: Mach)
            shutil.copyfile (config_ini.MESH_FILENAME        , mesh_name_dest   );
        print "...Done.\n\n";
    
    else :
        
        #--- An initial restart solution and an initial sensor are provided.
        #    -> Check their existence and copy them
        
        restart_name = "%s/%s" % (rootDir,config.ADAP_INI_RESTART_FILE);
        sensor_name  = "%s/%s" % (rootDir,config.ADAP_INI_SENSOR_FILE );
        mesh_name    = "%s/%s" % (rootDir,config.ADAP_INI_MESH_FILE   );
        
        error=0;
        
        if not os.path.isfile(restart_name):
            sys.stdout.write(' ## ERROR : No restart solution file was found.\n\n')
            error=1;
        
        if not os.path.isfile(sensor_name):
            sys.stdout.write(' ## ERROR : No restart sensor file was found.\n\n')
            error=1;
        
        if not os.path.isfile(mesh_name):
            sys.stdout.write(' ## ERROR : No restart mesh file was found.\n\n')
            error=1;
            
        if error == 1:
            sys.exit(1);
        
        shutil.copyfile (restart_name, restart_name_dest);   # restart solution
        shutil.copyfile (sensor_name , sensor_name_dest );   # the solution used for error estimation (migh be different than the restart, ex: Mach)
        shutil.copyfile (mesh_name   , mesh_name_dest   );
        
    
#: def amg_Inisol

def Parse_Adap_Options (config, config_adap):
    
    # --- Get the number of desired mesh adaptation loops
    #     and the mesh complexity for each of them
    
    config_adap.adap_complex = config['ADAP_COMPLEXITIES'].strip('()')
    config_adap.adap_complex = config_adap.adap_complex.split(",")
    
    config_adap.NbrIte = len(config_adap.adap_complex);
    
    if config_adap.NbrIte < 1 :
        print "  ## ERROR : Invalid number of iterations."
        sys.exit(1)
    
    # --- Get the number of sub-iterations for each mesh complexity
    
    config_adap.adap_subite = config['ADAP_SUBITE'].strip('()')
    config_adap.adap_subite = config_adap.adap_subite.split(",")
    
    if config_adap.NbrIte != len(config_adap.adap_subite) :
        print "  ## ERROR mesh adaptation: the number of mesh complexities (%d) and the number of sub-iterations (%d) don't match.\n" % (config_adap.NbrIte, config_adap.adap_subite)
        sys.exit(1)
    
    config_adap.NbrIteGlo = 0;
    for i in range(config_adap.NbrIte) :
        config_adap.adap_complex[i] = int(config_adap.adap_complex[i])
        config_adap.adap_subite[i]  = int(config_adap.adap_subite[i])
        config_adap.NbrIteGlo += config_adap.adap_subite[i];
    
    if not ('SU2_RUN' in os.environ ):
        os.environ['SU2_RUN'] = "";
    
    if 'ADAP_PATH' in config:
        path = "%s/" % config['ADAP_PATH'];
        os.environ['SU2_RUN'] = "%s:%s" % (os.environ['SU2_RUN'], path);
        os.environ['PATH'] = "%s:%s" % (os.environ['PATH'], path);
        config_adap.ADAP_PATH = path;

    else:
        config_adap.ADAP_PATH = "";
    
    if 'ADAP_BACK' in config:
        print " config['ADAP_BACK'] = %s\n" % config['ADAP_BACK'];
        if config['ADAP_BACK']=="YES":
            config_adap.ADAP_BACK = 1;
            
            if not 'ADAP_BACK_NAME' in config :
                print " ## ERROR : a back mesh name must be provided.\n ADAP_BACK option is ignored.\n";
                config_adap.ADAP_BACK = 0;
            else :
                config_adap.ADAP_BACK_NAME = config['ADAP_BACK_NAME'];
                
                #if not os.path.isfile(config_adap.ADAP_BACK_NAME):
                #    print " ## ERROR : Back mesh %s NOT FOUND!\n ADAP_BACK option is ignored.\n" %    config_adap.ADAP_BACK_NAME;
                #    config_adap.ADAP_BACK = 0;
        else:
            config_adap.ADAP_BACK = 0;
    else:
        config_adap.ADAP_BACK = 0;
            
    if (config_adap.ADAP_BACK):
        print "  -- Info : Use %s as a back mesh." % config_adap.ADAP_BACK_NAME;
    else:
        print "DONT GET BACK MESH";
    
    if 'ADAP_HIN' in config:
        config_adap.HMIN = float(config['ADAP_HMIN']);
    else:
        config_adap.HMIN = 0;
    
    if 'ADAP_HMAX' in config :
        config_adap.HMAX = float(config['ADAP_HMAX']);
    else:
        config_adap.HMAX = 1e300;
        
    if 'ADAP_HGRAD' in config:
        config_adap.HGRAD = float(config['ADAP_HGRAD']);
    else :
        config_adap.HGRAD = 3.0;
        
    config_adap.adap_ref = '';
    config_adap.RANS = 'NO';
    if config.PHYSICAL_PROBLEM == "NAVIER_STOKES":
        config_adap.RANS = 'YES';
        
        if not 'MARKER_HEATFLUX' in config:
            print "  ## ERROR RANS mesh adaptation: MARKER_HEATFLUX must be provided\n";
            sys.exit(1)
        
        adap_ref = [];
        
        if 'MARKER_FAR' in config:
            
            for i in range(len(config['MARKER_FAR'])):
                adap_ref.append(config['MARKER_FAR'][i]);
                
        if 'MARKER_SYM' in config:
        
            for i in range(len(config['MARKER_SYM'])):
                adap_ref.append(config['MARKER_SYM'][i]);
        
        
        config_adap.adap_ref = adap_ref;
    
    
    # --- Print summary
    
    #print "Mesh adaptation summary: \n"
    #
    #for i in range(config_adap.NbrIte) :
    
# def: Parse_Adap_Options

def Call_AMG (config_adap, config):
    
    ite_glo = config_adap.ite_glo;
    ite_cpx = config_adap.ite_cpx;
    ite_sub = config_adap.ite_sub;
    
    outNam = "current.new.meshb";
    itpNam = "current_restart.itp.solb";
    
    rootDir    = config_adap.rootDir;
    
    cpx = config_adap.adap_complex[ite_cpx];
    
    hmin  = config_adap.HMIN;
    hmax  = config_adap.HMAX;
    hgrad = config_adap.HGRAD;
    
    jobNam = "AMG.%d.%.0f.job" % (ite_glo,cpx);
    
    if config_adap.ADAP_PATH != "":
        path = config_adap.ADAP_PATH;
    else:
        path = "";
    
    if os.path.isfile(outNam):
        os.remove(outNam);
        
    if os.path.isfile(itpNam):
        os.remove(itpNam);
        
    if (config_adap.ADAP_BACK):
        back = " -back %s/%s " % (rootDir, config_adap.ADAP_BACK_NAME);
    else:
        back = "";
    
    adapsrc = "%s/adap.source" % rootDir;
    if not os.path.exists(adapsrc):
        adapsrc = "./adap.source";
        open(adapsrc, 'a').close();
    
    amg_cmd = "%samg -in current.meshb -sol current_sensor.solb -source %s -p 2 \
     -c %f -hgrad %.2f -hmin %le -hmax %le -out current.new.meshb \
    -itp current_restart.solb %s  -nordg " \
    % (path, adapsrc, cpx, hgrad, hmin, hmax, back)
    
    surf_ids = ""
    
    if ( config_adap.RANS == 'YES' ):
        # Freeze the boundary layer and adap rest of the domain    
        adap_ref = config_adap.adap_ref;
        
        if len(adap_ref) < 1:
            sys.stdout.write("  ## ERROR mesh adaptation: no surface IDs were provided.\n");
            sys.exit(1);
        
        surf_ids = "%d" % int(adap_ref[0]);
        for i in range(1,len(adap_ref)):
            surf_ids =  "%s,%d" % (surf_ids,int(adap_ref[i]))
        
        surf_ids = "-adap-domn-ids  0  -adap-surf-ids %s" % (surf_ids)
        
        source = ""
    else :
        source = "-source %s" % adapsrc;
    
    amg_cmd = "%s %s %s > %s" % (amg_cmd, surf_ids, source, jobNam)
        
    print " Running AMG \n Log: %s\n" % (jobNam);
    
    os.system(amg_cmd);
    
    print "...Done.\n\n";
    
    if not os.path.isfile(outNam):
        msg = "  ## ERROR : AMG failed at global iteration %d (see %s)\n" % (ite_glo, jobNam);
        msg = "     Command : %s\n\n" % amg_cmd;
        sys.stdout.write(msg);
        sys.exit(1);
        
    if not os.path.isfile(itpNam):
        msg = "  ## ERROR AMG : solution interpolation failed  at global iteration %d (see %s)\n" % (ite_glo, jobNam);
        msg = "     Command : %s\n\n" % amg_cmd;
        sys.stdout.write(msg);
        sys.exit(1);
    
    os.rename(itpNam, "current.new_ini.solb");
    
#: def amg_RunIte

def Call_SU2 (config_adap, config):
    
    outNam = 'current.new_restart.solb';
    
    ite_glo = config_adap.ite_glo;
    ite_cpx = config_adap.ite_cpx;
    
    rootDir    = config_adap.rootDir;
    adapDir    = config_adap.adapDir;
    partitions = config_adap.partitions;
        
    cpx = config_adap.adap_complex[ite_cpx];
    
    config_cur = copy.deepcopy(config);

    state = su2io.State()
    state.find_files(config_cur)
    
    jobNam = "SU2.%d.%.0f.job" % (ite_glo,cpx);
    
    config_cur.NUMBER_PART = partitions
    config_cur.CONSOLE = 'QUIET'
    config_cur.OUTPUT_LOG = jobNam;
    
    config_cur.RESTART_SOL            = 'YES'
    config_cur.MESH_FILENAME          = 'current.new.meshb'
    config_cur.SOLUTION_FLOW_FILENAME = 'current.new_ini.solb';
    config_cur.RESTART_FLOW_FILENAME  = outNam;
    config_cur.MESH_FORMAT            = 'INRIA';

    config_cur.VOLUME_FLOW_FILENAME = "flow.%d.%.0f" % (ite_glo,cpx);
    config_cur.SURFACE_FLOW_FILENAME = "surface_flow.%d.%.0f" % (ite_glo,cpx);

    #config_cur.EXT_ITER = 1;
    
    print " Running SU2_CFD \n Log: %s\n" % (jobNam) ;
    
    if os.path.isfile(outNam):
        os.remove(outNam);
    
    info = su2run.CFD(config_cur)
    
    print "...Done.\n\n";
    
    if not os.path.isfile(outNam):
        msg = "  ## ERROR : SU2 failed at global iteration %d\n" % (ite_glo);
        sys.stdout.write(msg);
        sys.exit(1);
        
    # Solution merging
    
    jobNam = "SU2_merge.%d.%.0f.job" % (ite_glo,cpx);
    config_cur.SOLUTION_FLOW_FILENAME = 'current.new_restart.dat';
    config_cur.OUTPUT_LOG = jobNam;
    
    print " Merging solutions \n Log: %s\n" % (jobNam) ;
    
    info = su2run.merge(config_cur)
    
    state.update(info)
    
    shutil.copyfile ("mach.solb", "current.new_sensor.solb");
    
#: Call_SU2

def adaptation_amg ( config, partitions=1,  kind='' ):
    
    # --- Warning before deleting old files?    
    warn = 1;    
    
    config_adap = su2io.Config();
    Parse_Adap_Options(config, config_adap);
    
    config_adap.partitions = partitions;
    
    NbrIte       = config_adap.NbrIte;
    adap_complex = config_adap.adap_complex;
    adap_subite  = config_adap.adap_subite;
    
    # ---------------------------------------------
    # --- Create adaptation folder
    # ---------------------------------------------
    
    #for i in range(NbrIte) :
    #    print "Iteration %d : Complexity = %d, %d sub-iterations.\n" % (i,adap_complex[i], adap_subite[i])
        
    rootDir = os.getcwd();
    adapDir = "%s/ADAP" % rootDir;
         
    config_adap.rootDir = rootDir;
    config_adap.adapDir = adapDir;
    
    if os.path.exists(adapDir):
        sys.stdout.write('./ADAP exists. Removing old mesh adaptation in 10s.')
        sys.stdout.flush();
        if warn : time.sleep(10);
        shutil.rmtree(adapDir);
        sys.stdout.write(' Done!\n\n')
  
    os.makedirs(adapDir);
    os.chdir(adapDir);
    
    # ---------------------------------------------
    # --- Initial solution
    # ---------------------------------------------
    
    amg_Inisol (config_adap, config);
    # Out : current.meshb
    #       current_restart.solb
    #       current_sensor.solb
    
    # ---------------------------------------------
    # --- Run mesh adaptation loop
    # ---------------------------------------------
    
    ite_glo = 1;
    for ite_cpx in range(NbrIte):
                
        NbrSub = adap_subite [ite_cpx];
        Cpx    = adap_complex[ite_cpx];
        
        config_adap.ite_cpx = ite_cpx;
        
        plu=""
        if NbrSub>1 : plu="s"
        
        print " -- Mesh complexity %.0f: %d sub-iteration%s.\n" % (adap_complex[ite_cpx], NbrSub, plu)
        
        for ite_sub in range(NbrSub) :
            
            config_adap.ite_sub = ite_sub; 
            
            print "    -- Sub-iteration %d\n" % ite_sub
            
            config_adap.ite_glo = ite_glo;
            
            #--- Mesh adaptation and solution interpolation
            
            # In : current.meshb
            #      current_restart.solb
            #      current_sensor.solb
            Call_AMG(config_adap, config);
            # Out : current.new.meshb
            #       current.new_ini.solb
            
            #--- CFD computation
            
            # In : current.new.meshb
            #      current.new_ini.solb
            Call_SU2(config_adap, config);
            # Out : current.new_restart.solb
            #       current.new_sensor.solb
            
            
            #--- Save files
            shutil.copyfile ("current.new.meshb", "ite.%d.%.0f.meshb"%(ite_glo, Cpx));
            shutil.copyfile ("current.new_sensor.solb", "ite.%d.%.0f.solb"%(ite_glo, Cpx));
            shutil.copyfile ("current.new_restart.solb", "ite.%d.%.0f_restart.solb"%(ite_glo, Cpx));
            
            historyNam = "history.dat"
            if os.path.exists(historyNam):
                shutil.copyfile (historyNam, "ite.%d.%.0f_history.dat"%(ite_glo, Cpx));
                
            #--- Rename files for the next iteration
            shutil.move("current.new.meshb"       , "current.meshb");
            shutil.move("current.new_restart.solb", "current_restart.solb");
            shutil.move("current.new_sensor.solb" , "current_sensor.solb");        
            
            ite_glo = ite_glo+1;
                        
        #: End sub-iteration
    
    #: ite_cpx
        
    sys.exit(1)
    
    
    