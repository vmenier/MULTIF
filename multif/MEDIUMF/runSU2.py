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
	config.CFL_NUMBER= '25';
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
	config.NUMBER_PART= 1;
	
	# --- Local relaxation / CFL
	#     Note: these options are only available in a custom version of su2:
	#     				https://github.com/vmenier/SU2/tree/darpa
	if (LocalRelax == 'YES') :
		config.RELAXATION_LOCAL= 'YES';
		config.CFL_ADAPT_LOCAL= 'YES';
		config.HARD_LIMITING_PARAM= '(0.15, 1e-5)';
		config.CFL_ADAPT_LOCAL_PARAM= '( 0.1, 1.5, 1e-12, 30.0 )';
		config.RESIDUAL_MAXVAL= 2;
				
	return config;
	
def runSU2 ( nozzle ):
	
	solver_options = Solver_Options();
	
	solver_options.Mach = nozzle.mission.mach;
	solver_options.Pres = nozzle.environment.P;
	solver_options.Temp = nozzle.environment.T;
	
	solver_options.InletPstag = nozzle.inlet.Pstag;
	solver_options.InletTstag = nozzle.inlet.Tstag;
	
	solver_options.LocalRelax = nozzle.LocalRelax;
	
	solver_options.NbrIte = 2; #300;
	
	solver_options.SU2_RUN = nozzle.SU2_RUN;
	
	solver_options.mesh_name    = nozzle.mesh_name;
	solver_options.restart_name = nozzle.restart_name;
	
	solver_options.convergence_order = nozzle.su2_convergence_order;
		
	GenerateNozzleMesh(nozzle);
	
	config = SetupConfig(solver_options);
	
	nozzle.OUTPUT_FORMAT = config['OUTPUT_FORMAT'];
	nozzle.CONV_FILENAME = config['CONV_FILENAME'];
	
	info = SU2.run.CFD(config);
	
	#return info;
	
	