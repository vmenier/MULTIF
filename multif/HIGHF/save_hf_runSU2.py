# -*- coding: utf-8 -*-

import sys

from hf_meshgeneration import *
from .. import SU2

class Solver_Options:
	def __init__(self):
		pass


def CheckSU2Version(nozzle):
	import subprocess;
	su2_exe = '%s/SU2_CFD' % nozzle.cfd.su2_run;
		
	#sys.path.append("/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/");
	#sys.pythonpath.append("/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/");
	#os.environ['PATH'] = ':'.join('/Users/menier/codes/SU2_DARPA/SU2_CFD/bin/')
	
	nozzle.cfd.su2_version = '';
	
	print "EXE = %s" % su2_exe;
	
	try :
		cmd = [su2_exe];
		out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None);
	except subprocess.CalledProcessError as err: 
		if ( 'DARPA' in err.output ):
			sys.stdout.write('Check SU2 version : OK\n');
			#nozzle.cfd.local_relax = 'YES';
			nozzle.cfd.local_relax = 'NO';
			nozzle.cfd.su2_version = 'OK';
		else:
			sys.stdout.write('\n');
			sys.stdout.write('#' * 90);
			sys.stdout.write('\n  ## WARNING : You are not using the right version of SU2. This may cause robustness issues.\n');
			sys.stdout.write('#' * 90);
			sys.stdout.write('\n\n');
			nozzle.cfd.local_relax = 'NO';
			nozzle.cfd.su2_version = 'NOT_OK';			

	
def SetupConfig (solver_options):
	
	config = SU2.io.Config();
	
	# --- Options
	
	Mach = solver_options.Mach;
	Pres = solver_options.Pres;
	Temp = solver_options.Temp;
	
	InletPstag = solver_options.InletPstag;
	InletTstag = solver_options.InletTstag;
	
	Pt = solver_options.Pt;
	Tt = solver_options.Tt;
	
	LocalRelax = solver_options.LocalRelax;
	
	NbrIte = solver_options.NbrIte;
	
	mesh_name = solver_options.mesh_name;
	restart_name = solver_options.restart_name;
	
	convergence_order = solver_options.convergence_order;
	
	Reynolds = solver_options.Reynolds;
	Reynolds_length = solver_options.Reynolds_length;
	
	method = solver_options.Method;
	
	Dim = solver_options.Dimension;
	
	# --- SU2_RUN
	
	config.SU2_RUN = solver_options.SU2_RUN;
	
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
			config.MARKER_EULER= '( 6 )';
		elif method == 'RANS':
			config.MARKER_HEATFLUX= '( 6, 0.0 )';
		config.MARKER_INLET= '( 1, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag);
		config.MARKER_FAR= '( 5, 4 )';
		config.MARKER_SYM= '( 2 )';
		config.MARKER_OUTLET= '( 3, %lf)' % (Pres);
	else:
		config.MARKER_EULER= '( 1, 2, 3, 4, \
		5, 6, 7, 8, 9, 10, \
		11, 12, 13, 14 )';
		config.MARKER_INLET= '( 16, %lf, %lf, 1.0, 0.0, 0.0 , 15, %lf, %lf, 1.0, 0.0, 0.0 )' % (Tt, Pt, InletTstag,InletPstag);
		config.MARKER_FAR= '(  17, 18, 21 )';
		config.MARKER_SYM= '( 19, 20 )';
		config.MARKER_OUTLET= '( 22, %lf)' % (Pres);
    
	#if Dim == '2D':
	#	if method == 'EULER':
	#		config.MARKER_EULER= '( PhysicalLine6 )';
	#	elif method == 'RANS':
	#		config.MARKER_HEATFLUX= '( PhysicalLine6, 0.0 )';
	#	config.MARKER_INLET= '( PhysicalLine1, %lf, %lf, 1.0, 0.0, 0.0 )' % (InletTstag,InletPstag);
	#	config.MARKER_FAR= '( PhysicalLine5, PhysicalLine4 )';
	#	config.MARKER_SYM= '( PhysicalLine2 )';
	#	config.MARKER_OUTLET= '( PhysicalLine3, %lf)' % (Pres);
	#else:
	#	config.MARKER_EULER= '( PhysicalSurface1, PhysicalSurface2, PhysicalSurface3, PhysicalSurface4, \
	#	PhysicalSurface5, PhysicalSurface6, PhysicalSurface7, PhysicalSurface8, PhysicalSurface9, PhysicalSurface10, \
	#	PhysicalSurface11, PhysicalSurface12, PhysicalSurface13, PhysicalSurface14 )';
	#	config.MARKER_INLET= '( PhysicalSurface16, %lf, %lf, 1.0, 0.0, 0.0 , PhysicalSurface15, %lf, %lf, 1.0, 0.0, 0.0 )' % (Tt, Pt, InletTstag,InletPstag);
	#	config.MARKER_FAR= '(  PhysicalSurface17, PhysicalSurface18, PhysicalSurface21 )';
	#	config.MARKER_SYM= '( PhysicalSurface19, PhysicalSurface20 )';
	#	config.MARKER_OUTLET= '( PhysicalSurface22, %lf)' % (Pres);
		
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
	config.OUTPUT_FORMAT= 'TECPLOT';
	config.CONV_FILENAME= 'history';
	config.RESTART_FLOW_FILENAME= restart_name;
	config.WRT_SOL_FREQ= '1000';
	config.WRT_CON_FREQ= '1';
	
	# --- Misc
	config.NUMBER_PART= 1;
	
	# --- Local relaxation / CFL
	#     Note: these options are only available in a custom version of su2:
	#     				https://github.com/vmenier/SU2/tree/darpa
	if (LocalRelax == 'YES' and method == 'EULER') :
		config.RELAXATION_LOCAL= 'YES';
		config.CFL_ADAPT_LOCAL= 'YES';
		config.HARD_LIMITING_PARAM= '(0.15, 1e-5)';
		config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 30.0 )';
		config.RESIDUAL_MAXVAL= 2;
				
	return config;


def HF_runSU2 ( nozzle ):
	
	solver_options = Solver_Options();
	
	solver_options.Method = nozzle.method;
	
	solver_options.Mach = nozzle.mission.mach;
	solver_options.Pres = nozzle.environment.P;
	solver_options.Temp = nozzle.environment.T;
	
	solver_options.InletPstag = nozzle.inlet.Pstag;
	solver_options.InletTstag = nozzle.inlet.Tstag;
	
	solver_options.LocalRelax = nozzle.cfd.local_relax;
	
	solver_options.Dimension = nozzle.dim;
	
	if nozzle.method == 'EULER':
		solver_options.NbrIte = 300;
	else:
		solver_options.NbrIte = 2000;
		
	solver_options.SU2_RUN = nozzle.cfd.su2_run;
	
	solver_options.mesh_name    = nozzle.cfd.mesh_name;
	solver_options.restart_name = nozzle.cfd.restart_name;
	
	solver_options.convergence_order = nozzle.cfd.su2_convergence_order;
	
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
	
	solver_options.Pt = Ps + 0.5*rho*U*U;
	solver_options.Tt = Ts*(1.+0.5*(gam-1.)*M*M);
	
	
	solver_options.Reynolds_length = D;
	solver_options.Reynolds        = Rey;
	
	config = SetupConfig(solver_options);
	
	nozzle.cfd.output_format = config['OUTPUT_FORMAT'];
	nozzle.cfd.conv_filename = config['CONV_FILENAME'];
	
	info = SU2.run.CFD(config);
	
	#return info;
	

	
