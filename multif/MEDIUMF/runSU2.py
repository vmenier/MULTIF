# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import multif

from meshgeneration import *
from .. import SU2

class Solver_Options:
    def __init__(self):
        pass


def CheckSU2Version(nozzle):
    import subprocess;
    su2_exe = '%s/SU2_CFD' % nozzle.SU2_RUN;
        
    #sys.path.append("/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/");
    #sys.pythonpath.append("/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/");
    #os.environ['PATH'] = ':'.join('/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/')
    
    nozzle.SU2Version = '';
    
    print "EXE = %s" % su2_exe;
    
    
    try :
        cmd = [su2_exe];
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None);
    except subprocess.CalledProcessError as err: 
        if ( 'DARPA' in err.output ):
            sys.stdout.write('Check SU2 version : OK\n');
            nozzle.LocalRelax = 'YES';
            nozzle.SU2Version = 'OK';
        else:
            sys.stdout.write('\n');
            sys.stdout.write('#' * 90);
            sys.stdout.write('\n  ## WARNING : You are not using the right version of SU2. This may cause robustness issues.\n');
            sys.stdout.write('#' * 90);
            sys.stdout.write('\n\n');
            nozzle.LocalRelax = 'NO';
            nozzle.SU2Version = 'NOT_OK';
            


def SetupConfig_old (solver_options):
    
    config = SU2.io.Config();
    
    # --- Options
    
    Mach = solver_options.Mach;
    Pres = solver_options.Pres;
    Temp = solver_options.Temp;
    
    InletPstag = solver_options.InletPstag;
    InletTstag = solver_options.InletTstag;
    
    LocalRelax = solver_options.LocalRelax;
    
    NbrIte = solver_options.NbrIte;
    
    mesh_name = solver_options.mesh_name;
    restart_name = solver_options.restart_name;
    
    convergence_order = solver_options.convergence_order;
    
    partitions = solver_options.partitions;

    MaxCFL = solver_options.MaxCFL;
    
    # --- SU2_RUN
    
    config.SU2_RUN = solver_options.SU2_RUN;
    
    # --- Governing
    config.PHYSICAL_PROBLEM= 'EULER';
    config.MATH_PROBLEM= 'DIRECT';
    config.RESTART_SOL= 'NO';
    config.SYSTEM_MEASUREMENTS= 'SI';
    
    # --- Free stream
    
    config.MACH_NUMBER='%lf' % Mach;
    config.FREESTREAM_PRESSURE='%lf' % Pres;
    config.FREESTREAM_TEMPERATURE='%lf' % Temp;
    config.REF_DIMENSIONALIZATION= 'DIMENSIONAL';
    
    # --- Boundary conditions
    
    config.MARKER_EULER= '( 9, 10, 11, 12 )';
    config.MARKER_INLET= '( 13, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag);
    config.MARKER_FAR= '( 6, 7, 8 )';
    config.MARKER_SYM= '( 1, 2 )';
    config.MARKER_OUTLET= '( 3, %lf,  4, %lf,  5, %lf)' % (Pres, Pres, Pres);
    
    # --- Numerical method
    
    config.NUM_METHOD_GRAD= 'WEIGHTED_LEAST_SQUARES';
    config.CFL_NUMBER= '%lf' % MaxCFL;
    config.CFL_ADAPT= 'NO';
    config.MAX_DELTA_TIME= '1E6';
    config.EXT_ITER= NbrIte;
    config.LINEAR_SOLVER= 'FGMRES';
    config.LINEAR_SOLVER_ERROR= '1E-6';
    config.LINEAR_SOLVER_ITER= '3';
    
    # --- Slope limiter
    
    config.REF_ELEM_LENGTH= '0.005 ';
    config.LIMITER_COEFF= '0.3';
    config.SHARP_EDGES_COEFF= '3.0';
    config.LIMITER_ITER= '150';
    config.REF_SHARP_EDGES= '3.0';
    config.SENS_REMOVE_SHARP= 'YES';
    
    # --- Multigrid
    
    config.MGLEVEL= '3';
    config.MGCYCLE= 'V_CYCLE';
    config.MG_PRE_SMOOTH= '( 1, 2, 3, 3 )';
    config.MG_POST_SMOOTH= '( 0, 0, 0, 0 )';
    config.MG_CORRECTION_SMOOTH= '( 0, 0, 0, 0 )';
    config.MG_DAMP_RESTRICTION= '0.75';
    config.MG_DAMP_PROLONGATION= '0.75';
    
    # --- Flow numerical method
    
    config.CONV_NUM_METHOD_FLOW= 'ROE';
    config.SPATIAL_ORDER_FLOW= '2ND_ORDER_LIMITER';
    config.SLOPE_LIMITER_FLOW= 'VENKATAKRISHNAN';
    config.AD_COEFF_FLOW= '( 0.15, 0.5, 0.05 )';
    config.TIME_DISCRE_FLOW= 'EULER_IMPLICIT';
    
    # --- Convergence parameters
    
    config.CONV_CRITERIA= 'RESIDUAL';
    config.RESIDUAL_REDUCTION= convergence_order;
    config.RESIDUAL_MINVAL= '-12';
    config.STARTCONV_ITER= '25';
    
    # --- Input/Output
    
    config.MESH_FILENAME= mesh_name;
    config.OUTPUT_FORMAT= 'TECPLOT';
    config.CONV_FILENAME= 'history';
    config.RESTART_FLOW_FILENAME= restart_name;
    config.WRT_SOL_FREQ= '1000';
    config.WRT_CON_FREQ= '1';
    
    # --- Misc
    config.NUMBER_PART= partitions;
    
    # --- Local relaxation / CFL
    #     Note: these options are only available in a custom version of su2:
    #                     https://github.com/vmenier/SU2/tree/darpa
    if (LocalRelax == 'YES') :
        config.RELAXATION_LOCAL= 'YES';
        config.CFL_ADAPT_LOCAL= 'YES';
        config.HARD_LIMITING_PARAM= '(0.15, 1e-5)';
        config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, %lf )' % MaxCFL;
        config.RESIDUAL_MAXVAL= 2;
                
    return config;
    
def runSU2_old ( nozzle ):
    
    solver_options = Solver_Options();
    
    solver_options.Mach = nozzle.mission.mach;
    solver_options.Pres = nozzle.environment.P;
    solver_options.Temp = nozzle.environment.T;
    
    solver_options.InletPstag = nozzle.inlet.Pstag;
    solver_options.InletTstag = nozzle.inlet.Tstag;
    
    solver_options.LocalRelax = nozzle.LocalRelax;
    
    solver_options.NbrIte = int(nozzle.su2_max_iterations);
    
    solver_options.SU2_RUN = nozzle.SU2_RUN;
    
    solver_options.mesh_name    = nozzle.mesh_name;
    solver_options.restart_name = nozzle.restart_name;
    
    solver_options.convergence_order = nozzle.su2_convergence_order;
    
    solver_options.MaxCFL = nozzle.MaxCFL;
    
    solver_options.Dimension = 2;
    
    
        
    GenerateNozzleMesh(nozzle);
    
    config = SetupConfig(solver_options);
    
    nozzle.OUTPUT_FORMAT = config['OUTPUT_FORMAT'];
    nozzle.CONV_FILENAME = config['CONV_FILENAME'];
    
    info = SU2.run.CFD(config);
    
    #return info;



def SetupConfig (solver_options):

    config = SU2.io.Config();

    # --- Options

    Mach = solver_options.Mach;
    Pres = solver_options.Pres;
    Temp = solver_options.Temp;

    InletPstag = solver_options.InletPstag;
    InletTstag = solver_options.InletTstag;

    LocalRelax = solver_options.LocalRelax;

    NbrIte = solver_options.NbrIte;

    mesh_name = solver_options.mesh_name;
    restart_name = solver_options.restart_name;

    convergence_order = solver_options.convergence_order;

    Reynolds = solver_options.Reynolds;
    Reynolds_length = solver_options.Reynolds_length;

    method = solver_options.Method;

    Dim = solver_options.Dimension;

    Pt = solver_options.Pt;
    Tt = solver_options.Tt;

	
    wall_temp = solver_options.wall_temp;
    wall_temp_values = solver_options.wall_temp_values;
	
    # --- SU2_RUN

    config.SU2_RUN = solver_options.SU2_RUN;

    config.NUMBER_PART =  solver_options.nproc;

    if Dim == '2D' :
        config.AXISYMMETRIC= 'YES';

        
    # --- Governing

    if method == 'EULER':
        config.PHYSICAL_PROBLEM= 'EULER';

        # --- Numerical method

        config.NUM_METHOD_GRAD= 'WEIGHTED_LEAST_SQUARES';
        config.CFL_NUMBER= '5';
        config.CFL_ADAPT= 'NO';
        config.MAX_DELTA_TIME= '1E6';
        config.LINEAR_SOLVER= 'FGMRES';
        config.LINEAR_SOLVER_ERROR= '1E-6';
        config.LINEAR_SOLVER_ITER= '3';

        config.LIMITER_ITER= '150';

    elif method == 'RANS':
        config.PHYSICAL_PROBLEM= 'NAVIER_STOKES';
        config.KIND_TURB_MODEL= 'SST'
        config.REYNOLDS_NUMBER= '%lf' % Reynolds;
        config.REYNOLDS_LENGTH= '%lf' % Reynolds_length;
        config.VISCOSITY_MODEL= 'SUTHERLAND';
        config.MU_CONSTANT= 1.716E-5;
        config.MU_REF= 1.716E-5;
        config.MU_T_REF= 273.15;

        config.NUM_METHOD_GRAD= 'GREEN_GAUSS';

        config.CFL_NUMBER= '5';
        config.CFL_ADAPT= 'NO';

        config.LINEAR_SOLVER= 'FGMRES';
        config.LINEAR_SOLVER_PREC= 'LU_SGS';
        config.LINEAR_SOLVER_ERROR= '1E-4';
        config.LINEAR_SOLVER_ITER= '3';


    config.MATH_PROBLEM= 'DIRECT';
    config.RESTART_SOL= 'NO';
    config.SYSTEM_MEASUREMENTS= 'SI';
    config.REGIME_TYPE= 'COMPRESSIBLE';

    config.EXT_ITER= NbrIte;

    config.RK_ALPHA_COEFF= "( 0.66667, 0.66667, 1.000000 )";

    # --- Free stream

    config.MACH_NUMBER='%lf' % Mach;


    config.FREESTREAM_PRESSURE='%lf' % Pres;
    config.FREESTREAM_TEMPERATURE='%lf' % Temp;
    config.REF_DIMENSIONALIZATION= 'DIMENSIONAL';

    # --- Boundary conditions

    if Dim == '2D':
        if method == 'EULER':
            config.MARKER_EULER= '( PhysicalLine1, PhysicalLine2, PhysicalLine3 )';
        elif method == 'RANS':
            config.MARKER_HEATFLUX= '( PhysicalLine1, 0.0, PhysicalLine2, 0.0, PhysicalLine3, 0.0 )';
        config.MARKER_INLET= '( PhysicalLine8, %lf, %lf, 1.0, 0.0, 0.0, PhysicalLine4,  %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag,Tt, Pt);
        config.MARKER_FAR= '( PhysicalLine5 )';
        config.MARKER_SYM= '( PhysicalLine7 )';
        config.MARKER_OUTLET= '( PhysicalLine6, %lf)' % (Pres);
    else:
        config.MARKER_EULER= '( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, \
        PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, \
        PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )';
        config.MARKER_INLET= '( PhysicalSurface15, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag);
        config.MARKER_FAR= '( PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )';
        config.MARKER_SYM= '( PhysicalSurface19, PhysicalSurface20 )';
        config.MARKER_OUTLET= '( PhysicalSurface22, %lf)' % (Pres);

        #MARKER_EULER= ( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )
        #MARKER_INLET= ( PhysicalSurface16,225.423, 8624.06, 1.0, 0.0, 0.0, PhysicalSurface15, 601, 275000, 1.0, 0.0, 0.0 )
        #MARKER_FAR= ( PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )
        #MARKER_SYM= ( PhysicalSurface19, PhysicalSurface20 )
        #MARKER_OUTLET= ( PhysicalSurface22, 7505.2400000)



    # --- Slope limiter

    config.REF_ELEM_LENGTH= '0.01 ';
    config.LIMITER_COEFF= '0.3';
    config.SHARP_EDGES_COEFF= '3.0';
    config.REF_SHARP_EDGES= '3.0';
    config.SENS_REMOVE_SHARP= 'NO';

    # --- Multigrid

    config.MGLEVEL= '3';
    config.MGCYCLE= 'V_CYCLE';
    config.MG_PRE_SMOOTH= '( 1, 2, 3, 3 )';
    config.MG_POST_SMOOTH= '( 0, 0, 0, 0 )';
    config.MG_CORRECTION_SMOOTH= '( 0, 0, 0, 0 )';
    config.MG_DAMP_RESTRICTION= '0.75';
    config.MG_DAMP_PROLONGATION= '0.75';

    # --- Flow numerical method
    if method == 'EULER':
        config.CONV_NUM_METHOD_FLOW= 'ROE';
        config.SPATIAL_ORDER_FLOW= '2ND_ORDER_LIMITER';
        config.SLOPE_LIMITER_FLOW= 'VENKATAKRISHNAN';
        config.AD_COEFF_FLOW= '( 0.15, 0.5, 0.05 )';
        config.TIME_DISCRE_FLOW= 'EULER_IMPLICIT';
    else :
        config.CONV_NUM_METHOD_FLOW= 'JST';
        config.SPATIAL_ORDER_FLOW= '2ND_ORDER_LIMITER';
        config.SLOPE_LIMITER_FLOW= 'VENKATAKRISHNAN';
        config.AD_COEFF_FLOW= '( 0.15, 0.5, 0.05 )';
        config.TIME_DISCRE_FLOW= 'EULER_IMPLICIT';
        config.ENTROPY_FIX_COEFF= 0.0;
        config.AD_COEFF_FLOW= "( 0.15, 0.5, 0.02 )";

        config.CONV_NUM_METHOD_TURB= 'SCALAR_UPWIND'
        config.SPATIAL_ORDER_TURB= '2ND_ORDER_LIMITER'
        config.SLOPE_LIMITER_TURB= 'VENKATAKRISHNAN'
        config.VISCOUS_LIMITER_TURB= 'NO'
        config.TIME_DISCRE_TURB= 'EULER_IMPLICIT'
        config.CFL_REDUCTION_TURB= '0.5'
        config.RELAXATION_FACTOR_TURB= '1.0'


    # --- Convergence parameters

    config.CONV_CRITERIA= 'RESIDUAL';
    config.RESIDUAL_REDUCTION= convergence_order;
    config.RESIDUAL_MINVAL= '-12';
    config.STARTCONV_ITER= '25';

    # --- Input/Output

    config.MESH_FILENAME= mesh_name;
    config.OUTPUT_FORMAT= solver_options.output_format;
    config.CONV_FILENAME= 'history';
    config.RESTART_FLOW_FILENAME= restart_name;
    config.WRT_SOL_FREQ= '1000';
    config.WRT_CON_FREQ= '1';
	
    # --- Local relaxation / CFL
    #     Note: these options are only available in a custom version of su2:
    #                     https://github.com/vmenier/SU2/tree/darpa
    if (LocalRelax == 'YES' and method == 'EULER') :
        config.RELAXATION_LOCAL= 'YES';
        config.CFL_ADAPT_LOCAL= 'YES';
        config.HARD_LIMITING_PARAM= '(0.15, 1e-5)';
        config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 30.0 )';
        config.RESIDUAL_MAXVAL= 2;
		
    # --- Setup wall temp distribution
    
    if ( wall_temp == 1 ) :
    	
    	nbv = len(wall_temp_values);
    	print "NBV %d" % nbv
    	print wall_temp_values
    	
    	temp_kwd = "%lf, %lf" % (wall_temp_values[0][0], wall_temp_values[0][1]);
    	
    	for i in range(1,nbv):
    		temp_kwd = "%s, %lf, %lf" % (temp_kwd, wall_temp_values[i][0], wall_temp_values[i][1]);
    
    	temp_kwd = "(%s)" % temp_kwd;
    	
    	config.MARKER_WALL_TEMP= "( PhysicalLine1 )";
    	config.WALL_TEMP_DEFINITION = temp_kwd;
    	
    return config;



def runSU2 ( nozzle ):

    solver_options = Solver_Options();

    solver_options.Method = nozzle.method;

    solver_options.Mach = nozzle.mission.mach;
    solver_options.Pres = nozzle.environment.P;
    solver_options.Temp = nozzle.environment.T;

    solver_options.InletPstag = nozzle.inlet.Pstag;
    solver_options.InletTstag = nozzle.inlet.Tstag;

    solver_options.LocalRelax = nozzle.LocalRelax;

    solver_options.NbrIte = int(nozzle.su2_max_iterations);
 
    solver_options.output_format = nozzle.OUTPUT_FORMAT;

    solver_options.SU2_RUN = nozzle.SU2_RUN;

    solver_options.mesh_name    = nozzle.mesh_name;
    solver_options.restart_name = nozzle.restart_name;

    solver_options.convergence_order = nozzle.su2_convergence_order;

    gam   = 1.4;
    R     = 287.06;
    Cv    = 717.645;
    Su    = 110.4;

    M      = nozzle.mission.mach;
    Ps     = nozzle.environment.P;
    Ts     = nozzle.environment.T;
    D      = nozzle.wall.geometry.radius(nozzle.wall.geometry.length);

    mu     = 1.716e-5*((Ts/273.15)**1.5)*(273.15 + Su)/(Ts + Su);      # Sutherland law 
    rho    = Ps / ( (gam-1.) * Cv * Ts )                               # density
    c      = np.sqrt( gam * (Ps/rho));                                 # speed of sound
    U      = M*c                                                       # velocity
    Rey    = rho*U*D/mu;                                               # Reynolds number

    solver_options.Reynolds_length = D;
    solver_options.Reynolds        = Rey;

    
    solver_options.nproc = nozzle.partitions;
    
    solver_options.Pt = Ps + 0.5*rho*U*U;
    solver_options.Tt = Ts*(1.+0.5*(gam-1.)*M*M);

    # --- Setup wall temperature distribution
    
    solver_options.wall_temp = 0;
    solver_options.wall_temp_values = [];

    if ( nozzle.wall_temp == 1 ) :
    	if ( nozzle.method != 'RANS' ):
    		sys.stderr.write('  ## ERROR : Wall temperature distribution only available for RANS.\n');
    		sys.exit(1);
    	solver_options.wall_temp = nozzle.wall_temp;
    	solver_options.wall_temp_values = nozzle.wall.temperature.thicknessNodes;

    #print "Rey %lf mu %lf rho %lf  T %lf  P %lf  D %lf" % (Rey, mu, rho, Ts, Ps,  D)
    #sys.exit(1)
    solver_options.Dimension = '2D';

    GenerateNozzleMesh(nozzle);

    config = SetupConfig(solver_options);

    nozzle.OUTPUT_FORMAT = config['OUTPUT_FORMAT'];
    nozzle.CONV_FILENAME = config['CONV_FILENAME'];

    info = SU2.run.CFD(config);

    #return info;
