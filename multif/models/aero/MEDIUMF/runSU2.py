"""
Run SU2 for 2D axisymmetric nozzle using Euler or RANS.
"""

import os, sys, shutil
import numpy as np

import multif
from meshgeneration import GenerateNozzleMesh, GenerateNozzleMesh_Deform
from meshgeneration import GenerateNozzleExitMesh
from .. import SU2

class Solver_Options:
    def __init__(self):
        pass


def CheckSU2Version(nozzle):
    import subprocess

    su2_exe = '%s/SU2_CFD' % nozzle.cfd.su2_run

    #sys.path.append("/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/")
    #sys.pythonpath.append("/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/")
    #os.environ['PATH'] = ':'.join('/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/')
    
    nozzle.cfd.su2_version = ''
    
    try :
        cmd = [su2_exe]
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
    except subprocess.CalledProcessError as err: 
        if ( 'DARPA' in err.output ):
            sys.stdout.write('Check SU2 version : OK\n')
            # Local relaxation has a few bugs and is found to not be necessary
            # for good convergence as of 10/3/17. Turn off for now.
            #nozzle.cfd.local_relax = 'YES'
            nozzle.cfd.local_relax = 'NO'
            nozzle.cfd.su2_version = 'OK'
        else:
            sys.stdout.write('\n')
            sys.stdout.write('#' * 90)
            sys.stdout.write('\n  ## WARNING : You are not using the right version of SU2. This may cause robustness issues.\n')
            sys.stdout.write('#' * 90)
            sys.stdout.write('\n\n')
            nozzle.cfd.local_relax = 'NO'
            nozzle.cfd.su2_version = 'NOT_OK'


def CheckSU2Convergence ( history_filename, field_name ) :

    #plot_format      = con.OUTPUT_FORMAT
    #plot_extension   = SU2.io.get_extension(plot_format)
    #history_filename = nozzle.cfd.conv_filename + plot_extension
    #special_cases    = SU2.io.get_specialCases(config)

    history      = SU2.io.read_history( history_filename )
    
    plot = SU2.io.read_plot(history_filename)

    Res = history[field_name]    
    
    NbrIte = len(Res)
    
    if ( NbrIte < 1 ):
        IniRes = 0
        FinRes = 0
    else :
        IniRes = Res[0]
        FinRes = Res[NbrIte-1]

    #print "Initial res = %le, Final res = %lf, DIFF = %lf\n" % (IniRes, FinRes, ResDif)
    return IniRes, FinRes
            
def SetupConfig (solver_options):
    
    config = SU2.io.Config()
    
    # --- Options
    
    Mach = solver_options.Mach
    Pres = solver_options.Pres
    Temp = solver_options.Temp
    
    mu_ref = solver_options.mu
    
    InletPstag = solver_options.InletPstag
    InletTstag = solver_options.InletTstag
    
    LocalRelax = solver_options.LocalRelax
    
    NbrIte = solver_options.NbrIte
    
    mesh_name = solver_options.mesh_name
    restart_name = solver_options.restart_name
    
    convergence_order = solver_options.convergence_order
    
    Reynolds = solver_options.Reynolds
    Reynolds_length = solver_options.Reynolds_length
    
    method = solver_options.Method
    
    Dim = solver_options.Dimension
    
    Pt = solver_options.Pt
    Tt = solver_options.Tt
    
    if hasattr(solver_options,'wall_temp'):
        wall_temp = solver_options.wall_temp
        wall_temp_values = solver_options.wall_temp_values
    else:
        wall_temp = 0
    
    # --- SU2_RUN
    
    config.SU2_RUN = solver_options.SU2_RUN
    
    config.NUMBER_PART =  solver_options.nproc
    
    if Dim == '2D' :
        config.AXISYMMETRIC= 'YES'
    
    # --- Governing
    
    if method == 'EULER':
        
        config.PHYSICAL_PROBLEM= 'EULER'
    
        # --- Numerical method
    
        config.NUM_METHOD_GRAD= 'WEIGHTED_LEAST_SQUARES'
        config.CFL_NUMBER= '5'
        config.CFL_ADAPT= 'NO'
        config.CFL_ADAPT_PARAM= '( 1.5, 0.5, 1.25, 50.0 )'
        config.MAX_DELTA_TIME= '1E6'
        config.LINEAR_SOLVER= 'FGMRES'
        config.LINEAR_SOLVER_ERROR= '1E-6'
        config.LINEAR_SOLVER_ITER= '3'
        
        config.LIMITER_ITER= '1000'
    
    elif method == 'RANS':
        config.PHYSICAL_PROBLEM= 'NAVIER_STOKES'
        config.KIND_TURB_MODEL= 'SST'
        config.REYNOLDS_NUMBER= '%lf' % Reynolds
        config.REYNOLDS_LENGTH= '%lf' % Reynolds_length
        config.VISCOSITY_MODEL= 'SUTHERLAND'
        config.MU_CONSTANT= 1.716E-5
        config.MU_REF= mu_ref
        config.MU_T_REF= Temp
        
        config.NUM_METHOD_GRAD= 'GREEN_GAUSS'
        
        config.CFL_NUMBER= '2'
        config.CFL_ADAPT= 'YES'
        config.CFL_ADAPT_PARAM= '( 1.5, 0.5, 1.25, 50.0 )'
        
        config.LINEAR_SOLVER= 'FGMRES'
        config.LINEAR_SOLVER_PREC= 'LU_SGS'
        config.LINEAR_SOLVER_ERROR= '1E-4'
        config.LINEAR_SOLVER_ITER= '3'
        
        config.LIMITER_ITER= '1000'
        
    config.MATH_PROBLEM= 'DIRECT'
    config.RESTART_SOL= 'NO'
    config.SYSTEM_MEASUREMENTS= 'SI'
    config.REGIME_TYPE= 'COMPRESSIBLE'
    
    config.EXT_ITER=  NbrIte
    
    config.RK_ALPHA_COEFF= "( 0.66667, 0.66667, 1.000000 )"
    
    # --- Free stream
    
    config.MACH_NUMBER='%lf' % Mach
    
    config.FREESTREAM_PRESSURE='%lf' % Pres
    config.FREESTREAM_TEMPERATURE='%lf' % Temp
    config.REF_DIMENSIONALIZATION= 'DIMENSIONAL'
    
    # --- Boundary conditions
    
    if Dim == '2D':
        if method == 'EULER':
            config.MARKER_EULER= '( 1, 2, 3 )'
        elif method == 'RANS':
            config.MARKER_HEATFLUX= '( 1, 0.0, 2, 0.0, 3, 0.0 )'
        config.MARKER_INLET= '( 8, %lf, %lf, 1.0, 0.0, 0.0)' % (InletTstag,InletPstag)
        config.MARKER_FAR= '( 4, 5, 6 )'
        config.MARKER_SYM= '( 7 )'
        #config.MARKER_OUTLET= '( 6, %lf)' % (Pres)
        #config.MARKER_THRUST= '( 9 )'
        config.MARKER_INTERNAL= '( 9 )'
    else:
        config.MARKER_EULER= '( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, \
        PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, \
        PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )'
        config.MARKER_INLET= '( PhysicalSurface15, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag)
        config.MARKER_FAR= '( PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )'
        config.MARKER_SYM= '( PhysicalSurface19, PhysicalSurface20 )'
        config.MARKER_OUTLET= '( PhysicalSurface22, %lf)' % (Pres)
    
        #MARKER_EULER= ( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )
        #MARKER_INLET= ( PhysicalSurface16,225.423, 8624.06, 1.0, 0.0, 0.0, PhysicalSurface15, 601, 275000, 1.0, 0.0, 0.0 )
        #MARKER_FAR= ( PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )
        #MARKER_SYM= ( PhysicalSurface19, PhysicalSurface20 )
        #MARKER_OUTLET= ( PhysicalSurface22, 7505.2400000)
    
    # --- Slope limiter
    
    config.REF_SHARP_EDGES= '3.0'
    config.SENS_REMOVE_SHARP= 'NO'
    
    # --- Multigrid
    
    config.MGLEVEL= '3'
    config.MGCYCLE= 'V_CYCLE'
    config.MG_PRE_SMOOTH= '( 1, 2, 3, 3 )'
    config.MG_POST_SMOOTH= '( 0, 0, 0, 0 )'
    config.MG_CORRECTION_SMOOTH= '( 0, 0, 0, 0 )'
    config.MG_DAMP_RESTRICTION= '0.75'
    config.MG_DAMP_PROLONGATION= '0.75'
    
    config.MUSCL_FLOW= 'YES'
    config.VENKAT_LIMITER_COEFF= 0.05
    config.JST_SENSOR_COEFF= '( 0.5, 0.02 )'
    
    # --- Flow numerical method
    if method == 'EULER':
        config.CONV_NUM_METHOD_FLOW= 'JST'
        config.SLOPE_LIMITER_FLOW= 'VENKATAKRISHNAN'
        config.TIME_DISCRE_FLOW= 'EULER_IMPLICIT'
    else :
        config.CONV_NUM_METHOD_FLOW= 'JST'
        config.SLOPE_LIMITER_FLOW= 'VENKATAKRISHNAN'
        config.TIME_DISCRE_FLOW= 'EULER_IMPLICIT'
        config.ENTROPY_FIX_COEFF= 0.0
    
        #config.CONV_NUM_METHOD_TURB= 'SCALAR_UPWIND'
        #config.SPATIAL_ORDER_TURB= '2ND_ORDER_LIMITER'
        #config.SLOPE_LIMITER_TURB= 'VENKATAKRISHNAN'
        #config.VISCOUS_LIMITER_TURB= 'NO'
        #config.TIME_DISCRE_TURB= 'EULER_IMPLICIT'
        #config.CFL_REDUCTION_TURB= '0.5'
        #config.RELAXATION_FACTOR_TURB= '1.0'
    
        config.CONV_NUM_METHOD_TURB= 'SCALAR_UPWIND'
        config.SLOPE_LIMITER_TURB= 'VENKATAKRISHNAN'
        config.TIME_DISCRE_TURB= 'EULER_IMPLICIT'
        config.CFL_REDUCTION_TURB= '0.5'
        config.RELAXATION_FACTOR_TURB= '1.0'
    
    # --- Convergence parameters
    
    if solver_options.gradients == 'ADJOINT':
        convergence_order = max(convergence_order, 8)
    
    config.CONV_CRITERIA= 'RESIDUAL'
    config.RESIDUAL_REDUCTION= convergence_order
    config.RESIDUAL_MINVAL= '-12'
    config.STARTCONV_ITER= '25'
    
    # --- Input/Output
        
    config.MESH_FILENAME= mesh_name
    config.OUTPUT_FORMAT= solver_options.output_format
    config.CONV_FILENAME= 'history'
    config.RESTART_FLOW_FILENAME= restart_name
    config.WRT_SOL_FREQ= '1000'
    config.WRT_CON_FREQ= '1'
    
    # --- Local relaxation / CFL
    #     Note: these options are only available in a custom version of su2:
    #                     https://github.com/vmenier/SU2/tree/darpa
    if (LocalRelax == 'YES' and method == 'EULER') :
        config.RELAXATION_LOCAL= 'YES'
        config.CFL_ADAPT_LOCAL= 'YES'
        config.HARD_LIMITING_PARAM= '(0.15, 1e-5)'
        #config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 30.0 )'
        config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 20.0 )'
        config.RESIDUAL_MAXVAL= 2
        
    # --- Setup wall temp distribution
    
    if ( wall_temp == 1 ) :
        
        nbv = len(wall_temp_values)
        
        temp_kwd = "%lf, %lf" % (wall_temp_values[0][0], wall_temp_values[0][1])
        
        for i in range(1,nbv):
            temp_kwd = "%s, %lf, %lf" % (temp_kwd, wall_temp_values[i][0], wall_temp_values[i][1])
    
        temp_kwd = "(%s)" % temp_kwd
        
        config.MARKER_WALL_TEMP= "( 1 )"
        config.WALL_TEMP_DEFINITION = temp_kwd
    
    #--- I/O
    
    config.WRT_BINARY_RESTART = 'NO'
    config.READ_BINARY_RESTART= 'NO'

    # --- Objective function
    config.OBJECTIVE_FUNCTION= 'THRUST_NOZZLE'
    
    #--- Deprecated options
    #config.REF_ELEM_LENGTH= '0.01 '
    #config.LIMITER_COEFF= '0.3'
    #config.SHARP_EDGES_COEFF= '3.0'
    #config.AD_COEFF_FLOW= "( 0.15, 0.5, 0.02 )"
    #config.THRUST_FILENAME= "thrust_nodef.dat"
    #config.SPATIAL_ORDER_FLOW= '2ND_ORDER_LIMITER'
    #config.AD_COEFF_FLOW= '( 0.15, 0.5, 0.05 )'
    
    return config
    
def checkResidual(config=[]):
    
    idRes = 13
    
    if( os.path.isfile('history.csv') ):
        history = np.loadtxt('history.csv',skiprows=1,delimiter=',')
        finalResidual = history[-1,11]
        residualReduction = max(history[:,idRes]) - history[-1,idRes]        
    elif( os.path.isfile('history.dat') ):
        history = np.loadtxt('history.dat',skiprows=3,delimiter=',')
        finalResidual = history[-1,idRes]
        residualReduction = max(history[:,idRes]) - history[-1,idRes]
    elif( os.path.isfile('history.vtk') ):
        history = np.loadtxt('history.vtk',skiprows=1,delimiter=',')
        
        filres=open('history.vtk','r')
        header=filres.readline()
        kwd=filres.readline().split(",")
        
        idRes=-1
        for i in range(len(kwd)):
            if kwd[i] == '"Res_Flow[0]"':
                idRes = i
                break
        
        filres.close()
        
        finalResidual = history[-1,idRes]
        residualReduction = max(history[:,idRes]) - history[-1,idRes]
    else: # bypass solution checking
        history = -1
        finalResidual = -1
        
        if "RESIDUAL_REDUCTION" in config:
            residualReduction = config.RESIDUAL_REDUCTION -1
        sys.stderr.write("  ## ERROR No history file found.\n")
        sys.exit()
                
    return history, finalResidual, residualReduction


def callSU2(config):
    """
    Helper function for dispatching SU2, detecting divergence failure and
    checking residual decrease.

    Inputs:
    config: dictionary of SU2 configuration keyword --> value pairs

    Returns:
    history: Numpy array containing data from SU2's history file
    finalResidual: final log of density residual
    residualReduction: orders of magnitude decrease in density residual
    """

    try:
        info = SU2.run.CFD(config)
        print(info)
    except SU2.DivergenceFailure as e:   
        su2history = open('about.txt','a')
        su2history.write('SU2 calculation with baseline params diverged.\n')
        su2history.close()
        # raise e
    # any other failure should stop program and be reported

    # --- Check SU2 residual here    
    history, finalResidual, residualReduction = checkResidual(config)       
    su2history = open('about.txt','a')
    su2history.write('Final residual: %0.16f\n' % finalResidual)
    su2history.write('Residual reduction: %0.16f\n' % residualReduction)
    su2history.close()

    return history, finalResidual, residualReduction


def runSU2(nozzle, sst_perturbation=None, output='verbose'):
    """
    Provided a nozzle instance, dispatch SU2 to determine requested aerodynamic
    quantities of interest. Medium-fidelity (2D axisymmetric Euler or RANS).

    Arguments:
    sst_perturbation: either 'SSTC1', 'SSTC2', 'SSTC3', 'SSTP1C1', 'SSTP1C2', or None
        specifying the type of SST perturbation that should be used, if any
    """

    # --- Setup solver options
    solver_options = Solver_Options()
    
    solver_options.Method = nozzle.method
    
    solver_options.Mach = nozzle.mission.mach
    solver_options.Pres = nozzle.environment.P
    solver_options.Temp = nozzle.environment.T
    
    solver_options.InletPstag = nozzle.inlet.Pstag
    solver_options.InletTstag = nozzle.inlet.Tstag
    
    solver_options.LocalRelax = nozzle.cfd.local_relax
    solver_options.NbrIte = int(nozzle.cfd.su2_max_iterations)
    solver_options.output_format = nozzle.cfd.output_format
    solver_options.SU2_RUN = nozzle.cfd.su2_run
    
    solver_options.mesh_name    = nozzle.cfd.mesh_name
    solver_options.restart_name = nozzle.cfd.restart_name
    solver_options.convergence_order = nozzle.cfd.su2_convergence_order
    
    solver_options.dv_coefs = []
    
    # --- Specify wall variable data when adjoint gradients are requested
    if( nozzle.gradientsMethod == 'ADJOINT' ):
    
        iTag = -1
        for i in range(len(nozzle.DV_Tags)):
            Tag = nozzle.DV_Tags[i]
            if (Tag == "WALL"):
                iTag = i
                break
        
        if ( iTag < 0 ):
            sys.stderr.write("  ## ERROR SU2 adjoint computation: Wall parameterization not specified.\n")
            sys.exit()

        nbr_dv = max(nozzle.wall.dv)+1

        for i in range(nbr_dv):
            id_dv = nozzle.DV_Head[iTag] + i
            solver_options.dv_coefs.append(nozzle.dvList[id_dv])
    
        #for iCoef in range(len(nozzle.wall.dv)):
        #    id_dv = nozzle.DV_Head[iTag] + nozzle.wall.dv[iCoef]  
        #    if id_dv >= nozzle.DV_Head[iTag]:
        #        print "id_dv %d val %lf" % (id_dv, nozzle.dvList[id_dv])
        #        solver_options.dv_coefs.append(nozzle.dvList[id_dv])
        
        solver_options.wall_coefs_dv = nozzle.wall.dv
        solver_options.gradients     = nozzle.gradientsMethod
        
        solver_options.param = nozzle.param
        
        if nozzle.param == '2D':
            solver_options.wall_coefs    = nozzle.wall.coefs
        else:
            solver_options.gradients        = nozzle.gradientsMethod
            solver_options.centerline_coefs = nozzle.wall.centerline.coefs
            solver_options.majoraxis_coefs  = nozzle.wall.majoraxis.coefs
            solver_options.minoraxis_coefs  = nozzle.wall.minoraxis.coefs
            
    else:
    
        solver_options.gradients     = nozzle.gradientsMethod
    
    gam   = 1.4
    R     = 287.06
    Cv    = 717.645
    Su    = 110.4
    
    M      = nozzle.mission.mach
    Ps     = nozzle.environment.P
    Ts     = nozzle.environment.T
    D      = nozzle.wall.geometry.radius(nozzle.wall.geometry.length)
    
    mu     = 1.716e-5*((Ts/273.15)**1.5)*(273.15 + Su)/(Ts + Su)      # Sutherland law 
    rho    = Ps / ( (gam-1.) * Cv * Ts )                               # density
    c      = np.sqrt( gam * (Ps/rho))                                 # speed of sound
    U      = M*c                                                       # velocity
    Rey    = rho*U*D/mu                                               # Reynolds number
    
    solver_options.Reynolds_length = D
    solver_options.Reynolds        = Rey
        
    solver_options.mu = mu

    solver_options.nproc = nozzle.cpusPerTask
    
    solver_options.Pt = Ps + 0.5*rho*U*U
    solver_options.Tt = Ts*(1.+0.5*(gam-1.)*M*M)
    
    # --- Setup wall temperature distribution
    
    solver_options.wall_temp = 0
    solver_options.wall_temp_values = []
    
    if ( nozzle.wallTempFlag == 1 ) :
        if ( nozzle.method != 'RANS' ):
            raise RuntimeError('Wall temperature distribution only available for RANS.')
        solver_options.wall_temp = nozzle.wallTempFlag
        solver_options.wall_temp_values = nozzle.wall.temperature.thicknessNodes

    # --- Setup meshes
    solver_options.Dimension = '2D'
    
    if ( nozzle.meshDeformationFlag ):
        GenerateNozzleMesh_Deform(nozzle)
    else:
        GenerateNozzleMesh(nozzle)
        
    if nozzle.method == "RANS":
        #GenerateNozzleExitMesh(nozzle)
        solver_options.NbrIte = max(solver_options.NbrIte,5000)
    
    # --- Setup config file options
    config = SetupConfig(solver_options)
    
    nozzle.cfd.output_format = config['OUTPUT_FORMAT']
    nozzle.cfd.conv_filename = config['CONV_FILENAME']

    # --- Add config file options for SST RANS perturbations if necessary

    if sst_perturbation is not None:
        # Update config file parameters common to all perturbations
        config['RESTART_SOL'] = 'NO'
        config['USING_UQ'] = 'YES'
        config['BETA_DELTA'] = 1
        config['URLX'] = 0.1

        # Specific parameters for each perturbation
        componentality = [1, 2, 3, 1, 2]
        permute = ['NO', 'NO', 'NO', 'YES', 'YES']

        # nozzle.cfd.sst_tags is ['SSTC1', 'SSTC2', 'SSTC3', 'SSTP1C1', 'SSTP1C2']
        i = nozzle.cfd.sst_tags.index(sst_perturbation)
        config['COMPONENTALITY'] = componentality[i]
        config['PERMUTE'] = permute[i]  

        # XXX Efficiency can be increasing by linking previously obtained
        # SU2 mesh instead of redeforming baseline mesh 
        # linkFiles = ['nozzle.su2']

        # XXX Does this need to be changed for the adjoint as well?
    
    # --- Remove thrust file if it exists
    
    if os.path.exists("thrust_nodef.dat"): os.remove("thrust_nodef.dat")
    
    # --- Setup file and flags to record SU2's progress
    
    su2history = open('about.txt','w')
    su2history.close()
    # savefiles = ['config_CFD.cfg','history.csv','history.dat','history.vtk',
    # 'nozzle.su2','nozzle.dat'] # list of files to save in case of SU2 failure
    # dirname = 'su2_run1' # local subdirectory to store files in case of SU2 failure

    # --- Run SU2
    callSU2(config)

    # # --- Rerun SU2 if necessary (EULER only)
    # if( nozzle.method == 'EULER' ):
    #     if( finalResidual > 0 ): # SU2 diverged

    #         # Copy all original files into a directory for troubleshooting
    #         if not os.path.isdir(dirname):
    #             os.mkdir(dirname)
    #         for f in savefiles:
    #             if( os.path.isfile(f) ):
    #                 shutil.copyfile(f,os.path.join(dirname,f))

    #         # Notify reason for 2nd run
    #         su2history = open('about.txt','a')
    #         su2history.write('Rerunning SU2 with more conservative parameters.\n')
    #         su2history.write('Failed run files are stored in %s.\n' % dirname)
    #         su2history.close() 

    #         # Change SU2 parameters
    #         config = SetupConfig(solver_options)
    #         if( solver_options.LocalRelax == 'YES' ):
    #             config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 10.0 )'
    #             config.LIMITER_ITER= int(config.LIMITER_ITER)*2. 
    #             config.EXT_ITER = int(config.EXT_ITER)*3.           
    #         else:
    #             config.CFL_NUMBER = float(config.CFL_NUMBER)/2.
    #             config.LIMITER_ITER = int(config.LIMITER_ITER)*2.
    #             config.EXT_ITER = int(config.EXT_ITER)*3.

    #         # Rerun SU2      
    #         info = SU2.run.CFD(config)

    #         # Check residual
    #         history, finalResidual, residualReduction = checkResidual(config)
    #         su2history = open('about.txt','a')
    #         su2history.write('Final residual: %0.16f\n' % finalResidual)
    #         su2history.write('Residual reduction: %0.16f\n' % residualReduction)
    #         su2history.close() 

    #     elif( residualReduction < config.RESIDUAL_REDUCTION ): # residual not reduced enough

    #         # Copy all original files into a directory for troubleshooting
    #         if not os.path.isdir(dirname):
    #             os.mkdir(dirname)
    #         for f in savefiles:
    #             if( os.path.isfile(f) ):
    #                 shutil.copyfile(f,os.path.join(dirname,f))

    #         # Notify reason for 2nd run
    #         su2history = open('about.txt','a')
    #         su2history.write('Rerunning SU2 since residual reduction is inadequate.\n')
    #         su2history.write('Initial run files are stored in %s.\n' % dirname)
    #         su2history.close()

    #         # Change SU2 parameters
    #         config = SetupConfig(solver_options)
    #         if( os.path.isfile('nozzle.dat') ): # implement restart

    #             config.RESTART_SOL= 'YES' # restart from previous solution
    #             os.rename('nozzle.dat','solution_flow.dat')
                
    #             # Gauge progress made in last N iter
    #             N = min(300,history[:,0].size)
    #             modVarLast300Iters = np.mean(np.abs(history[-N:-1,11] - np.mean(history[-N:-1,11])))
    #             if( modVarLast300Iters > 0.5 ): # sufficient progress made, make no changes
    #                 pass # originally 0.2 was used
    #             elif( solver_options.LocalRelax == 'YES' ): # no progress made, tighten CFL
    #                 su2history = open('about.txt','a')
    #                 su2history.write('More conservative parameters are being used.\n')
    #                 su2history.close() 
    #                 config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 10.0 )'
    #                 config.LIMITER_ITER= int(config.LIMITER_ITER)*2 
    #                 config.EXT_ITER = int(config.EXT_ITER)*2              
    #             else:
    #                 su2history = open('about.txt','a')
    #                 su2history.write('More conservative parameters are being used.\n')
    #                 su2history.close() 
    #                 config.CFL_NUMBER = float(config.CFL_NUMBER)/2
    #                 config.LIMITER_ITER = int(config.LIMITER_ITER)*2
    #                 config.EXT_ITER = int(config.EXT_ITER)*2

    #             config.RESIDUAL_REDUCTION = float(config.RESIDUAL_REDUCTION) - residualReduction

    #         else: # possibly a different residual diverged, etc.

    #             su2history = open('about.txt','a')
    #             su2history.write('More conservative parameters are being used.')
    #             su2history.close() 

    #             if( solver_options.LocalRelax == 'YES' ):
    #                 config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 10.0 )'
    #                 config.LIMITER_ITER= int(config.LIMITER_ITER)*2 
    #                 config.EXT_ITER = int(config.EXT_ITER)*3           
    #             else:
    #                 config.CFL_NUMBER = float(config.CFL_NUMBER)/2
    #                 config.LIMITER_ITER = int(config.LIMITER_ITER)*2
    #                 config.EXT_ITER = int(config.EXT_ITER)*3

    #         # Rerun SU2
    #         info = SU2.run.CFD(config)

    #         # Check residual
    #         history, finalResidual, residualReduction = checkResidual(config)
    #         su2history = open('about.txt','a')
    #         su2history.write('Final residual: %0.16f\n' % finalResidual)
    #         su2history.write('Residual reduction: %0.16f\n' % residualReduction)
    #         su2history.close() 
        
    # --- Adjoint computation (if required)
    
    if nozzle.gradientsFlag and nozzle.qoi.getGradient('THRUST') is not None:
        
        if nozzle.gradientsMethod == 'ADJOINT':
            
            # --- AD            
            config_AD = setupConfig_AD (solver_options)
            info = SU2.run.CFD(config_AD)
            
            # --- DOT            
            config_DOT = setupConfig_DOT (solver_options)            
            sys.stdout.write("  -- Running SU2_DOT\n")            
            info = SU2.run.DOT(config_DOT)
            
            # --- Check convergence            
            # IniRes, FinRes = CheckSU2Convergence("history_adj.dat", "Res_AdjFlow[0]")
            
            # if IniRes < FinRes:
            #     raise RuntimeWarning("Discrete adjoint solution is NOT converged. Finite differences will be run instead.")                
            # else :
            #     sys.stdout.write("  -- Info : Discrete adjoint solution converged.\n")
            #     # Assumes gradient of thrust w.r.t. all design variables not governing nozzle wall are 0
            #     gtmp = np.zeros((len(nozzle.derivativesDV),))
            #     gtmp[0:max(nozzle.wall.dv)+1] = Read_Gradients_AD(nozzle)
            #     nozzle.qoi.setGradient('THRUST', gtmp)
                        
        else:
            
            # Estimate gradients via finite difference after return call
            print('WARNING: Thrust gradients are desired, but adjoint will not be used. Are you sure?')
          
    return


def setupConfig_AD (solver_options):
    # ---
    
    config = SU2.io.Config()
    
    # --- Options
    
    Mach = solver_options.Mach
    Pres = solver_options.Pres
    Temp = solver_options.Temp
    
    InletPstag = solver_options.InletPstag
    InletTstag = solver_options.InletTstag
    
    LocalRelax = solver_options.LocalRelax
    
    NbrIte = solver_options.NbrIte
    
    mesh_name = solver_options.mesh_name
    restart_name = solver_options.restart_name
    
    convergence_order = solver_options.convergence_order
    
    partitions = solver_options.nproc
    
    Reynolds = solver_options.Reynolds
    Reynolds_length = solver_options.Reynolds_length
    
    method = solver_options.Method
    
    Dim = solver_options.Dimension
    
    Pt = solver_options.Pt
    Tt = solver_options.Tt
    
    if hasattr(solver_options,'wall_temp'):
        wall_temp = solver_options.wall_temp
        wall_temp_values = solver_options.wall_temp_values
    else:
        wall_temp = 0
    
    
    config.NUMBER_PART= partitions
    
    config.SU2_RUN = solver_options.SU2_RUN
    
    
    # --- Governing
    
    if method == 'EULER':
        config.PHYSICAL_PROBLEM= 'EULER'
    
        # --- Numerical method
    
        config.NUM_METHOD_GRAD= 'WEIGHTED_LEAST_SQUARES'
        config.CFL_NUMBER= '15'
        config.CFL_ADAPT= 'NO'
        config.MAX_DELTA_TIME= '1E6'
        config.LINEAR_SOLVER= 'FGMRES'
        config.LINEAR_SOLVER_ERROR= '1E-6'
        config.LINEAR_SOLVER_ITER= '3'
    
        config.LIMITER_ITER= '200'
    
    elif method == 'RANS':
        config.PHYSICAL_PROBLEM= 'NAVIER_STOKES'
        config.KIND_TURB_MODEL= 'SST'
        config.REYNOLDS_NUMBER= '%lf' % Reynolds
        config.REYNOLDS_LENGTH= '%lf' % Reynolds_length
        config.VISCOSITY_MODEL= 'SUTHERLAND'
        config.MU_CONSTANT= 1.716E-5
        config.MU_REF= 1.716E-5
        config.MU_T_REF= 273.15
    
        config.NUM_METHOD_GRAD= 'GREEN_GAUSS'
    
        config.CFL_NUMBER= '5'
        config.CFL_ADAPT= 'NO'
        
        config.LINEAR_SOLVER= 'FGMRES'
        config.LINEAR_SOLVER_PREC= 'LU_SGS'
        config.LINEAR_SOLVER_ERROR= '1E-4'
        config.LINEAR_SOLVER_ITER= '3'
    
    config.SYSTEM_MEASUREMENTS= 'SI'
    config.REGIME_TYPE= 'COMPRESSIBLE'
    
    config.EXT_ITER= NbrIte
    
    config.RK_ALPHA_COEFF= "( 0.66667, 0.66667, 1.000000 )"
    
    # ---
    
    config.MATH_PROBLEM= 'DISCRETE_ADJOINT'
    
    config.AXISYMMETRIC= 'YES'
    config.CFL_ADAPT= 'NO'
    
    #config.LINEAR_SOLVER= 'FGMRES'
    #config.LINEAR_SOLVER_ERROR= 1E-6
    #config.LINEAR_SOLVER_ITER= 10
    
    config.RESTART_SOL= 'NO'
    config.EXT_ITER= 1000
    config.RK_ALPHA_COEFF= '( 0.66667, 0.66667, 1.000000 )'
    config.MACH_NUMBER= 0.511000
    config.FREESTREAM_PRESSURE= 18754.000000
    config.FREESTREAM_TEMPERATURE= 216.700000
    config.REF_DIMENSIONALIZATION= 'DIMENSIONAL'
    
    #config.MARKER_EULER= '( ( PhysicalLine1, PhysicalLine2, PhysicalLine3 ) )'
    #config.MARKER_INLET= '( PhysicalLine8, 955.000000, 97585.000000, 1.0, 0.0, 0.0, PhysicalLine4,  228.016984, 22181.944264, 1.0, 0.0, 0.0 )'
    #config.MARKER_FAR= '( ( PhysicalLine5 ) )'
    #config.MARKER_SYM= '( ( PhysicalLine7 ) )'
    #config.MARKER_OUTLET= '( PhysicalLine6, 18754.000000)'
    
    # --- Boundary conditions
    
    if Dim == '2D':
        if method == 'EULER':
            config.MARKER_EULER= '( 1, 2, 3, 10 )'
        elif method == 'RANS':
            config.MARKER_HEATFLUX= '( 1, 0.0, 2, 0.0, 3, 0.0, 10, 0.0 )'
        #config.MARKER_INLET= '( 8, %lf, %lf, 1.0, 0.0, 0.0, 4,  %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag,Tt, Pt)
        config.MARKER_INLET= '( 8, %lf, %lf, 1.0, 0.0, 0.0)' % (InletTstag,InletPstag)
        config.MARKER_FAR= '( 4, 5, 6 )'
        config.MARKER_SYM= '( 7 )'
        #config.MARKER_OUTLET= '( 6, %lf)' % (Pres)
        #config.MARKER_THRUST= '( 9 )'
        config.MARKER_INTERNAL= '( 9 )'
    else:
        config.MARKER_EULER= '( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, \
        PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, \
        PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )'
        config.MARKER_INLET= '( PhysicalSurface15, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag)
        config.MARKER_FAR= '( PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )'
        config.MARKER_SYM= '( PhysicalSurface19, PhysicalSurface20 )'
        config.MARKER_OUTLET= '( PhysicalSurface22, %lf)' % (Pres)

    #config.REF_ELEM_LENGTH= 0.01  # deprecated
    #config.LIMITER_COEFF= 0.3
    #config.SHARP_EDGES_COEFF= 3.0
    config.REF_SHARP_EDGES= 3.0
    config.SENS_REMOVE_SHARP= 'YES'
    config.MGLEVEL= 3
    config.MGCYCLE= 'V_CYCLE'
    config.MG_PRE_SMOOTH= '( 1, 2, 3, 3 )'
    config.MG_POST_SMOOTH= '( 0, 0, 0, 0 )'
    config.MG_CORRECTION_SMOOTH= '( 0, 0, 0, 0 )'
    config.MG_DAMP_RESTRICTION= 0.75
    config.MG_DAMP_PROLONGATION= 0.75
    config.CONV_NUM_METHOD_FLOW= 'JST'
    #config.SPATIAL_ORDER_FLOW= '2ND_ORDER_LIMITER'
    config.SLOPE_LIMITER_FLOW= 'VENKATAKRISHNAN'
    #config.AD_COEFF_FLOW= '( 0.15, 0.5, 0.02 )'
    config.TIME_DISCRE_FLOW= 'EULER_IMPLICIT'
    config.CONV_CRITERIA= 'RESIDUAL'
    config.RESIDUAL_REDUCTION= 10
    config.RESIDUAL_MINVAL= -200
    config.STARTCONV_ITER= 25
    config.OUTPUT_FORMAT= 'TECPLOT'
    config.CONV_FILENAME= 'history_adj'
    config.RESTART_ADJ_FILENAME= 'nozzle_adj.dat'
    config.WRT_SOL_FREQ= 100
    config.WRT_CON_FREQ= 1
    
    config.MESH_FILENAME= 'nozzle.su2'
    config.SOLUTION_FLOW_FILENAME= 'nozzle.dat'
    
    config.OBJECTIVE_FUNCTION= 'THRUST_NOZZLE'
    
    
    #--- ADJ options
    
    config.OBJECTIVE_WEIGHT = 1.0
    config.REF_SHARP_EDGES= 3.0
    config.SENS_REMOVE_SHARP= 'YES'
    config.CONV_NUM_METHOD_ADJFLOW= 'JST'
    config.MUSCL_ADJFLOW= 'YES'
    config.SLOPE_LIMITER_ADJFLOW= 'VENKATAKRISHNAN'
    config.ADJ_SHARP_LIMITER_COEFF= 3.0
    config.ADJ_JST_SENSOR_COEFF= '( 0.0, 0.00 )'
    config.CFL_REDUCTION_ADJFLOW= 0.75
    config.TIME_DISCRE_ADJFLOW= 'EULER_IMPLICIT'
    config.LIMIT_ADJFLOW= '1E6'
    config.MG_ADJFLOW= 'YES'    
    
    config.MARKER_MONITORING= '( 9 )'


    #--- I/O
    
    config.WRT_BINARY_RESTART = 'NO'
    config.READ_BINARY_RESTART= 'NO'
    
    return config
    
def setupConfig_DOT (solver_options):
    
    # ---

    #config = SU2.io.Config('dot.cfg')
    #return config
    
    config = SU2.io.Config()
    
    # --- Options
    
    Mach = solver_options.Mach
    Pres = solver_options.Pres
    Temp = solver_options.Temp
    
    InletPstag = solver_options.InletPstag
    InletTstag = solver_options.InletTstag
    
    LocalRelax = solver_options.LocalRelax
    
    NbrIte = solver_options.NbrIte
    
    mesh_name = solver_options.mesh_name
    restart_name = solver_options.restart_name
    
    convergence_order = solver_options.convergence_order
    
    partitions = 1#solver_options.nproc
    
    Reynolds = solver_options.Reynolds
    Reynolds_length = solver_options.Reynolds_length
    
    method = solver_options.Method
    
    Dim = solver_options.Dimension
    
    Pt = solver_options.Pt
    Tt = solver_options.Tt
    
    if hasattr(solver_options,'wall_temp'):
        wall_temp = solver_options.wall_temp
        wall_temp_values = solver_options.wall_temp_values
    else:
        wall_temp = 0
    
    
    config.NUMBER_PART= partitions
    
    config.SU2_RUN = solver_options.SU2_RUN
    
    dv_coefs = solver_options.dv_coefs
    
    param = solver_options.param
    
    if param == '2D':
        wall_coefs    = solver_options.wall_coefs    
    else :
        centerline_coefs = solver_options.centerline_coefs
        majoraxis_coefs  = solver_options.majoraxis_coefs 
        minoraxis_coefs  = solver_options.minoraxis_coefs 
        
    wall_coefs_dv = solver_options.wall_coefs_dv 
    
    
    
    # --- Governing
    
    if method == 'EULER':
        config.PHYSICAL_PROBLEM= 'EULER'
    
        # --- Numerical method
    
        config.NUM_METHOD_GRAD= 'WEIGHTED_LEAST_SQUARES'
        config.CFL_NUMBER= '15'
        config.CFL_ADAPT= 'YES'
        config.CFL_ADAPT_PARAM= '( 1.5, 0.5, 1.25, 50.0 )'
        config.MAX_DELTA_TIME= '1E6'
        config.LINEAR_SOLVER= 'FGMRES'
        config.LINEAR_SOLVER_ERROR= '1E-6'
        config.LINEAR_SOLVER_ITER= '3'
    
        config.LIMITER_ITER= '200'
    
    elif method == 'RANS':
        config.PHYSICAL_PROBLEM= 'NAVIER_STOKES'
        config.KIND_TURB_MODEL= 'SST'
        config.REYNOLDS_NUMBER= '%lf' % Reynolds
        config.REYNOLDS_LENGTH= '%lf' % Reynolds_length
        config.VISCOSITY_MODEL= 'SUTHERLAND'
        config.MU_CONSTANT= 1.716E-5
        config.MU_REF= 1.716E-5
        config.MU_T_REF= 273.15
    
        config.NUM_METHOD_GRAD= 'GREEN_GAUSS'
    
        config.CFL_NUMBER= '5'
        
        config.CFL_ADAPT= 'YES'
        config.CFL_ADAPT_PARAM= '( 1.5, 0.5, 1.25, 50.0 )'
    
        config.LINEAR_SOLVER= 'FGMRES'
        config.LINEAR_SOLVER_PREC= 'LU_SGS'
        config.LINEAR_SOLVER_ERROR= '1E-4'
        config.LINEAR_SOLVER_ITER= '3'
    
    config.SYSTEM_MEASUREMENTS= 'SI'
    config.REGIME_TYPE= 'COMPRESSIBLE'
    
    config.EXT_ITER= NbrIte
    
    config.RK_ALPHA_COEFF= "( 0.66667, 0.66667, 1.000000 )"
    
    # -------
    
    
    config.WRT_BINARY_RESTART = 'NO'
    config.READ_BINARY_RESTART= 'NO'
    
    
    config.MATH_PROBLEM= 'DISCRETE_ADJOINT'
    
    config.AXISYMMETRIC= 'YES'
    config.CFL_ADAPT= 'NO'
    
    #config.LINEAR_SOLVER= 'FGMRES'
    #config.LINEAR_SOLVER_ERROR= 1E-6
    #config.LINEAR_SOLVER_ITER= 10
    
    config.RESTART_SOL= 'NO'
    config.EXT_ITER= 1000
    config.RK_ALPHA_COEFF= '( 0.66667, 0.66667, 1.000000 )'
    config.MACH_NUMBER= 0.511000
    config.FREESTREAM_PRESSURE= 18754.000000
    config.FREESTREAM_TEMPERATURE= 216.700000
    config.REF_DIMENSIONALIZATION= 'DIMENSIONAL'
    
    
    
    #config.MARKER_EULER= '( ( PhysicalLine1, PhysicalLine2, PhysicalLine3 ) )'
    #config.MARKER_INLET= '( PhysicalLine8, 955.000000, 97585.000000, 1.0, 0.0, 0.0, PhysicalLine4,  228.016984, 22181.944264, 1.0, 0.0, 0.0 )'
    #config.MARKER_FAR= '( ( PhysicalLine5 ) )'
    #config.MARKER_SYM= '( ( PhysicalLine7 ) )'
    #config.MARKER_OUTLET= '( PhysicalLine6, 18754.000000)'
    
    # --- Boundary conditions
    
    if Dim == '2D':
        if method == 'EULER':
            config.MARKER_EULER= '( 1, 2, 3, 10 )'
        elif method == 'RANS':
            config.MARKER_HEATFLUX= '( 1, 0.0, 2, 0.0, 3, 0.0, 10, 0.0 )'
        config.MARKER_INLET= '( 8, %lf, %lf, 1.0, 0.0, 0.0, 4,  %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag,Tt, Pt)
        config.MARKER_FAR= '( 5 )'
        config.MARKER_SYM= '( 7 )'
        config.MARKER_OUTLET= '( 6, %lf)' % (Pres)
        #config.MARKER_THRUST= '( 9 )'
        config.MARKER_INTERNAL= '( 9 )'
    else:
        config.MARKER_EULER= '( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, \
        PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, \
        PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )'
        config.MARKER_INLET= '( PhysicalSurface15, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag)
        config.MARKER_FAR= '( PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )'
        config.MARKER_SYM= '( PhysicalSurface19, PhysicalSurface20 )'
        config.MARKER_OUTLET= '( PhysicalSurface22, %lf)' % (Pres)
    
    #config.REF_ELEM_LENGTH= 0.01 
    #config.LIMITER_COEFF= 0.3
    #config.SHARP_EDGES_COEFF= 3.0
    config.REF_SHARP_EDGES= 3.0
    config.SENS_REMOVE_SHARP= 'YES'
    config.MGLEVEL= 3
    config.MGCYCLE= 'V_CYCLE'
    config.MG_PRE_SMOOTH= '( 1, 2, 3, 3 )'
    config.MG_POST_SMOOTH= '( 0, 0, 0, 0 )'
    config.MG_CORRECTION_SMOOTH= '( 0, 0, 0, 0 )'
    config.MG_DAMP_RESTRICTION= 0.75
    config.MG_DAMP_PROLONGATION= 0.75
    config.CONV_NUM_METHOD_FLOW= 'JST'
    #config.SPATIAL_ORDER_FLOW= '2ND_ORDER_LIMITER'
    config.SLOPE_LIMITER_FLOW= 'VENKATAKRISHNAN'
    #config.AD_COEFF_FLOW= '( 0.15, 0.5, 0.02 )'
    config.TIME_DISCRE_FLOW= 'EULER_IMPLICIT'
    config.CONV_CRITERIA= 'RESIDUAL'
    config.RESIDUAL_REDUCTION= 10
    config.RESIDUAL_MINVAL= -200
    config.STARTCONV_ITER= 25
    config.OUTPUT_FORMAT= 'TECPLOT'
    config.CONV_FILENAME= 'history_adj'
    config.RESTART_ADJ_FILENAME= 'nozzle_adj.dat'
    config.WRT_SOL_FREQ= 100
    config.WRT_CON_FREQ= 1
    
    config.MESH_FILENAME= 'nozzle.su2'
    config.SOLUTION_FLOW_FILENAME= 'nozzle.dat'
    
    config.OBJECTIVE_FUNCTION= 'THRUST_NOZZLE'
    
    config.SOLUTION_ADJ_FILENAME= 'nozzle_adj.dat'
    
    
    # ------------- DOT PARAMETERS -------------

    config.GEO_MODE= 'FUNCTION'
    config.GEO_MARKER= '( 1 )'
    config.DV_MARKER= '( 1 )'
    
    

    
    if param == '2D':
        
        # --- SETUP DV_PARAM
    
        dv_Parameters = []
        dv_FFDTag     = []
        dv_Size       = []
        
        dv_kind  = ""
        dv_value = []
    
        for i in range(len(dv_coefs)):
        
            dv_Parameters.append([float(dv_coefs[i])])
            dv_Size.append(1)
            dv_FFDTag.append([])
    
            dv_kind = dv_kind + "BSPLINECOEF "
            dv_value.append(0.001)
    
        config.DV_KIND = dv_kind    
        config.DV_PARAM = { 'FFDTAG' : dv_FFDTag     ,
                           'PARAM'  : dv_Parameters ,
                           'SIZE'   : dv_Size}
        config.DV_VALUE = dv_value
        
        bsplinecoefs    = "%lf" % wall_coefs[0]
        bsplinecoefs_dv = "%d" % (wall_coefs_dv[0]+1)
        
        for i in range(1, len(wall_coefs)):
            bsplinecoefs    = bsplinecoefs + ", %lf" % wall_coefs[i]
            bsplinecoefs_dv = bsplinecoefs_dv + ", %d" % (wall_coefs_dv[i]+1)
        bsplinecoefs    = "(%s)" % bsplinecoefs
        bsplinecoefs_dv = "(%s)" % bsplinecoefs_dv
        
        config.BSPLINECOEFS    = bsplinecoefs
        config.BSPLINECOEFS_DV = bsplinecoefs_dv
    
    else:
        
        print "3D PARAM"
        
        idx=0
        
        wrt_center_coefs = "%lf" % centerline_coefs[0]
        wrt_center_dv    = "%d" % (wall_coefs_dv[idx]+1)
        idx+=1
        for i in range(1, len(centerline_coefs)):
            wrt_center_coefs = wrt_center_coefs + ", %lf" % centerline_coefs[i]
            wrt_center_dv    = wrt_center_dv + ", %d" % (wall_coefs_dv[idx]+1)
            idx+=1
            
        wrt_r1_coefs = "%lf" % majoraxis_coefs[0]
        wrt_r1_dv    = "%d" % (wall_coefs_dv[idx]+1)
        idx+=1
        for i in range(1, len(majoraxis_coefs)):
            wrt_r1_coefs = wrt_r1_coefs + ", %lf" % majoraxis_coefs[i]
            wrt_r1_dv    = wrt_r1_dv + ", %d" % (wall_coefs_dv[idx]+1)
            idx+=1
            
        wrt_r2_coefs = "%lf" % minoraxis_coefs[0]
        wrt_r2_dv    = "%d" % (wall_coefs_dv[idx]+1)
        idx+=1
        for i in range(1, len(minoraxis_coefs)):
            wrt_r2_coefs = wrt_r2_coefs + ", %lf" % minoraxis_coefs[i]
            wrt_r2_dv    = wrt_r2_dv + ", %d" % (wall_coefs_dv[idx]+1)
            idx+=1
        
        #--- DV kind, param, value
            
        dv_Parameters = []
        dv_FFDTag     = []
        dv_Size       = []
        
        dv_kind  = ""
        dv_value = []
    
        for i in range(len(dv_coefs)):
        
            dv_Parameters.append([float(dv_coefs[i])])
            dv_Size.append(1)
            dv_FFDTag.append([])
    
            dv_kind = dv_kind + "NOZZLE3DCOEF "
            dv_value.append(0.001)
    
        config.DV_KIND = dv_kind    
        config.DV_PARAM = { 'FFDTAG' : dv_FFDTag     ,
                           'PARAM'  : dv_Parameters ,
                           'SIZE'   : dv_Size}
        config.DV_VALUE = dv_value
        
        
        #wrt_dv_kind  = []
        #wrt_dv_param = []
        #wrt_dv_value = []
        #for i in range(len(wall_coefs_dv)):
        #     wrt_dv_kind.append("NOZZLE3DCOEF")
        #     wrt_dv_param.append(1)
        #     wrt_dv_value.append(0.001)
        #
        #config.DV_KIND  = wrt_dv_kind  
        #config.DV_PARAM = wrt_dv_param
        #config.DV_VALUE = wrt_dv_value
        
        config.NOZZLE3DCOEFS1    = '(%s)' % wrt_center_coefs
        config.NOZZLE3DCOEFS1_DV = '(%s)' % wrt_center_dv
        config.NOZZLE3DCOEFS2    = '(%s)' % wrt_r1_coefs
        config.NOZZLE3DCOEFS2_DV = '(%s)' % wrt_r1_dv
        config.NOZZLE3DCOEFS3    = '(%s)' % wrt_r2_coefs
        config.NOZZLE3DCOEFS3_DV = '(%s)' % wrt_r2_dv
        
    
    config.DEFORM_LINEAR_ITER= 100
    config.DEFORM_NONLINEAR_ITER= 50
    config.DEFORM_CONSOLE_OUTPUT= 'YES'
    config.DEFORM_TOL_FACTOR= 0.0001
    config.DEFORM_STIFFNESS_TYPE= 'WALL_DISTANCE'
    
    #config.SAVE_DEF_FILE= "YES"
    
    return config


def SetupConfig_DEF (solver_options):
    
    config    = SU2.io.Config()

    # --- Options

    Mach = solver_options.Mach
    Pres = solver_options.Pres
    Temp = solver_options.Temp

    InletPstag = solver_options.InletPstag
    InletTstag = solver_options.InletTstag

    partitions = solver_options.nproc

    Reynolds = solver_options.Reynolds
    Reynolds_length = solver_options.Reynolds_length

    method = solver_options.Method

    Dim = solver_options.Dimension

    Pt = solver_options.Pt
    Tt = solver_options.Tt
    
    Dim = solver_options.Dimension
    
    config.SU2_RUN = solver_options.SU2_RUN

    # --- 

    config.MESH_FILENAME= 'nozzle.su2'
    config.DV_KIND= 'SURFACE_FILE'
    
    config.DV_PARAM = { 'FFDTAG' : [[]]     ,
                       'PARAM'  : [[1]] ,
                       'SIZE'   : [1]}
    config.DV_VALUE = [0.001]
    
    
    config.DV_MARKER= '( 1 )'
    config.MOTION_FILENAME= 'mesh_motion.dat'
    config.DEFORM_LINEAR_SOLVER= 'FGMRES'
    config.DEFORM_LINEAR_ITER= 200
    config.DEFORM_NONLINEAR_ITER= 5
    config.DEFORM_CONSOLE_OUTPUT= 'YES'
    config.DEFORM_TOL_FACTOR= 1e-6
    config.DEFORM_STIFFNESS_TYPE= 'WALL_DISTANCE'
    config.HOLD_GRID_FIXED= 'NO'
    config.HOLD_GRID_FIXED_COORD= '(-1e6,-1e6,-1e6,1e6,1e6,1e6)'
    config.VISUALIZE_DEFORMATION= 'YES'
    config.MARKER_MOVING= '( 1 )'
    config.NUMBER_PART= 1

    # --- Boundary conditions

    if Dim == '2D':
        if method == 'EULER':
            config.MARKER_EULER= '( 1, 2, 3, 10 )'
            config.DEFORM_LINEAR_ITER= 200
        elif method == 'RANS':
            config.MARKER_HEATFLUX= '( 1, 0.0, 2, 0.0, 3, 0.0, 10, 0.0 )'
            config.DEFORM_LINEAR_ITER= 500
        config.MARKER_INLET= '( 8, %lf, %lf, 1.0, 0.0, 0.0, 4,  %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag,Tt, Pt)
        config.MARKER_FAR= '( 5 )'
        config.MARKER_SYM= '( 7 )'
        config.MARKER_OUTLET= '( 6, %lf)' % (Pres)
        #config.MARKER_THRUST= '( 9 )'
        config.MARKER_INTERNAL= '( 9 )'
    else:
        config.MARKER_EULER= '( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, \
        PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, \
        PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )'
        config.MARKER_INLET= '( PhysicalSurface15, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag)
        config.MARKER_FAR= '( PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )'
        config.MARKER_SYM= '( PhysicalSurface19, PhysicalSurface20 )'
        config.MARKER_OUTLET= '( PhysicalSurface22, %lf)' % (Pres)
    
    
    return config

def Compute_Thrust_Gradients_FD (nozzle):
    
    nbr_dv = max(nozzle.wall.dv)+1
    
    thrust_nodef = Get_Thrust_File(nozzle)
    
    if thrust_nodef < 0 :
        sys.stderr.write("  ## ERROR Compute_Thrust_Gradients_FD : No baseline thrust value was found.\n")
        sys.exit(1)
        thrust_nodef = 0.0
    
    thrust_grad = np.zeros(nbr_dv)
    
    ## --- Load deformation files 
    ##     (They are outputs from custom SU2_CFD_AD)
    #
    #hdl_def = []
    #
    #for i in range(nbr_dv):
    #    filNam = "wall_%d.dat" % i
    #    try:
    #        hdl_def.append(np.loadtxt(filNam))
    #    except:
    #        sys.stderr.write("  ## ERROR Compute_Thrust_Gradients_FD : %s not found. Abort.\n" % filNam)
    #        return thrust_grad
    #
    #    sys.stdout.write("%s loaded.\n" % filNam)
    
    # --- Call deformation + CFD for each point
    
    
    solver_options = Solver_Options()
    
    solver_options.Method = nozzle.method
    
    solver_options.Mach = nozzle.mission.mach
    solver_options.Pres = nozzle.environment.P
    solver_options.Temp = nozzle.environment.T
    
    solver_options.InletPstag = nozzle.inlet.Pstag
    solver_options.InletTstag = nozzle.inlet.Tstag
    
    solver_options.LocalRelax = nozzle.cfd.local_relax
    
    solver_options.NbrIte = int(nozzle.cfd.su2_max_iterations)
    
    solver_options.output_format = nozzle.cfd.output_format
    
    solver_options.SU2_RUN = nozzle.cfd.su2_run
    
    solver_options.mesh_name    = nozzle.cfd.mesh_name
    solver_options.restart_name = nozzle.cfd.restart_name
    
    solver_options.convergence_order = nozzle.cfd.su2_convergence_order
    
    solver_options.dv_coefs = []
    
    iTag = -1
    for i in range(len(nozzle.DV_Tags)):
        Tag = nozzle.DV_Tags[i]
        if (Tag == "WALL"):
            iTag = i
            break
        
    if ( iTag < 0 ):
        sys.stderr.write("  ## ERROR SU2 adjoint computation: Wall parameterization not specified.\n")
        sys.exit()
    
    nbr_dv = max(nozzle.wall.dv)+1
    
    for i in range(nbr_dv):
        id_dv = nozzle.DV_Head[iTag] + i
        print "id_dv %d val %lf" % (id_dv, nozzle.dvList[id_dv])
        solver_options.dv_coefs.append(nozzle.dvList[id_dv])
    
    
    #for iCoef in range(len(nozzle.wall.dv)):
    #    id_dv = nozzle.DV_Head[iTag] + nozzle.wall.dv[iCoef]  
    #    if id_dv >= nozzle.DV_Head[iTag]:
    #        print "id_dv %d val %lf" % (id_dv, nozzle.dvList[id_dv])
    #        solver_options.dv_coefs.append(nozzle.dvList[id_dv])
    
    solver_options.gradients     = nozzle.gradientsMethod
    solver_options.wall_coefs    = nozzle.wall.coefs
    solver_options.wall_coefs_dv = nozzle.wall.dv
    
    gam   = 1.4
    R     = 287.06
    Cv    = 717.645
    Su    = 110.4
    
    M      = nozzle.mission.mach
    Ps     = nozzle.environment.P
    Ts     = nozzle.environment.T
    D      = nozzle.wall.geometry.radius(nozzle.wall.geometry.length)
    
    mu     = 1.716e-5*((Ts/273.15)**1.5)*(273.15 + Su)/(Ts + Su)      # Sutherland law 
    rho    = Ps / ( (gam-1.) * Cv * Ts )                               # density
    c      = np.sqrt( gam * (Ps/rho))                                 # speed of sound
    U      = M*c                                                       # velocity
    Rey    = rho*U*D/mu                                               # Reynolds number
    
    solver_options.Reynolds_length = D
    solver_options.Reynolds        = Rey
    
    #solver_options.nproc = nozzle.partitions
    solver_options.nproc = nozzle.cpusPerTask
    
    solver_options.Pt = Ps + 0.5*rho*U*U
    solver_options.Tt = Ts*(1.+0.5*(gam-1.)*M*M)
    
    # --- Setup wall temperature distribution
    
    solver_options.wall_temp = 0
    solver_options.wall_temp_values = []
    
    if ( nozzle.wallTempFlag == 1 ) :
        if ( nozzle.method != 'RANS' ):
            sys.stderr.write('  ## ERROR : Wall temperature distribution only available for RANS.\n')
            sys.exit(1)
        solver_options.wall_temp = nozzle.wallTempFlag
        solver_options.wall_temp_values = nozzle.wall.temperature.thicknessNodes
    
    #print "Rey %lf mu %lf rho %lf  T %lf  P %lf  D %lf" % (Rey, mu, rho, Ts, Ps,  D)
    #sys.exit(1)
    solver_options.Dimension = '2D'
    
    GenerateNozzleMesh(nozzle)
    
    config = SetupConfig(solver_options)
    
    nozzle.cfd.output_format = config['OUTPUT_FORMAT']
    nozzle.cfd.conv_filename = config['CONV_FILENAME']
    
    config_DEF = SetupConfig_DEF (solver_options)
    
    # --- Check whether mesh deformation input files (wall_*.dat) exist
    #     If not, generate them using SU2
    
    flag=0
    for idv in range(nbr_dv):
        if not os.path.exists(config_DEF.MOTION_FILENAME):
            flag=1
            break
    
    if flag == 1:
        # --- Generate wall_*.dat files (inputs to SU2_DEF)
        config_AD = setupConfig_AD (solver_options)
        config_AD.EXT_ITER=0
        info = SU2.run.CFD(config_AD)
        config_DOT = setupConfig_DOT (solver_options)
        info = SU2.run.DOT(config_DOT)
    
    for idv in range(nbr_dv):
    
        # --- Call def
        
        mesh_out_filename = "nozzle_%d.su2" % idv
        
        config_DEF.MOTION_FILENAME   = "wall_%d.dat" % idv
        
        if not os.path.exists(config_DEF.MOTION_FILENAME):
            sys.stderr.write("  ## ERROR FD gradients: %s not found.\n" % config_DEF.MOTION_FILENAME)
            return thrust_grad
        
        config_DEF.MESH_OUT_FILENAME = mesh_out_filename
        config_DEF.SU2_RUN = solver_options.SU2_RUN
        
        if os.path.exists(mesh_out_filename): os.remove(mesh_out_filename)
        
        sys.stdout.write("  -- Running SU2_DEF\n")
        
        SU2.run.DEF(config_DEF)
        
        # --- Call CFD
        
        thrust_filename = 'thrust_%d.dat' % idv # output from SU2 containing thrust
        
        if os.path.exists(thrust_filename): os.remove(thrust_filename)
        
        config_CFD = SetupConfig(solver_options)
        
        config_CFD.OBJECTIVE_FUNCTION= 'THRUST_NOZZLE'
        #config_CFD.MARKER_THRUST= '( 9 ) '
        config.MARKER_INTERNAL= '( 9 )'
        config_CFD.THRUST_FILENAME= thrust_filename
        config_CFD.MESH_FILENAME= config_DEF.MESH_OUT_FILENAME
        config_CFD.RESTART_FLOW_FILENAME= "nozzle_%d.dat" % idv
        #config_CFD.EXT_ITER= 1
        
        info = SU2.run.CFD(config_CFD)
        
        if not os.path.exists(thrust_filename): 
            sys.stderr.write("  ## ERROR Compute_Thrust_Gradients_FD : output thrust from SU2 not found.\n \
            Are you using the right SU2 version?\n")
            return thrust_grad
        
        thrust = np.loadtxt(thrust_filename)
        thrust_grad[idv] = thrust-thrust_nodef
                
    return thrust_grad


def Read_Gradients_AD (nozzle):

    nbr_dv = max(nozzle.wall.dv)+1
    
    hdl_grad = np.loadtxt('./of_grad.dat', skiprows=1)

    thrust_grad = np.zeros(len(nozzle.dvList))

    iTag = -1
    for i in range(len(nozzle.DV_Tags)):
        Tag = nozzle.DV_Tags[i]
        if (Tag == "WALL"):
            iTag = i
            break

    if ( iTag < 0 ):
        sys.stderr.write("  ## ERROR SU2 adjoint computation: Wall parameterization not specified.\n")
        sys.exit()

    nbr_dv = max(nozzle.wall.dv)+1

    for i in range(nbr_dv):
        id_dv = nozzle.DV_Head[iTag] + i
        thrust_grad[id_dv] = hdl_grad[i]
        #print "id_dv %d val %lf" % (id_dv, nozzle.dvList[id_dv])
        
    return thrust_grad
