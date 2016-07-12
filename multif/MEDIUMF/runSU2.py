# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import multif
from .. import SU2

def SetupConfig (nozzle):
	
	config = SU2.io.Config();
	
	# --- Governing
	
	config.PHYSICAL_PROBLEM= 'EULER';
	config.MATH_PROBLEM= 'DIRECT';
	config.RESTART_SOL= 'NO';
	config.SYSTEM_MEASUREMENTS= 'SI';
	
	# --- Free stream
	
	config.MACH_NUMBER='0.9';
	config.FREESTREAM_PRESSURE='23842.787765';
	config.FREESTREAM_TEMPERATURE='218.807923';
	config.REF_DIMENSIONALIZATION= 'DIMENSIONAL';
	
	# --- Boundary conditions
	
	config.MARKER_EULER= '( 9, 10, 11, 12 )';
	config.MARKER_INLET= '( 13, 1021.500000, 144925.000000, 1.0, 0.0, 0.0 )';
	config.MARKER_FAR= '( 6, 7, 8 )';
	config.MARKER_SYM= '( 1, 2 )';
	config.MARKER_OUTLET= '( 3, 23842.787765,  4, 23842.787765,  5, 23842.787765)';
	
	# --- Numerical method
	
	config.NUM_METHOD_GRAD= 'WEIGHTED_LEAST_SQUARES';
	config.CFL_NUMBER= '25';
	config.CFL_ADAPT= 'NO';
	config.MAX_DELTA_TIME= '1E6';
	config.EXT_ITER= 500;
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
	config.RESIDUAL_REDUCTION= '3';
	config.RESIDUAL_MINVAL= '-12';
	config.STARTCONV_ITER= '25';
	
	# --- Input/Output
	
	config.MESH_FILENAME= 'axinoz.su2';
	config.OUTPUT_FORMAT= 'TECPLOT';
	config.CONV_FILENAME= 'history';
	config.RESTART_FLOW_FILENAME= 'restart_flow.dat';
	config.WRT_SOL_FREQ= '1000';
	config.WRT_CON_FREQ= '1';
	
	# --- Misc
	config.NUMBER_PART= 1;
	
	# --- Local relaxation / CFL
	
	config.RELAXATION_LOCAL= 'YES';
	
	return config;
	
	