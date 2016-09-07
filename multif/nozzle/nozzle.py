# -*- coding: utf-8 -*-


"""

R. Fenrich & V. Menier, July 2016

"""

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
from .. import SU2
import multif

import material
import component
import inlet
import environment
import fluid
import mission
import tolerance
import geometry

from .. import _meshutils_module
import ctypes
import numpy as np

from parserDV import *

class Nozzle:
	def __init__(self):
		pass
		
	def SetupDV (self, config):
		nozzle = self;
		
		coefs_size = nozzle.coefs_size;

		NbrDVTot = 0;		
		
				
		if 'DV_LIST' in config:

			dv_keys = ('WALL', 'INLET_TSTAG', 'INLET_PSTAG', 'THERMAL_CONDUCTIVITY', \
			'THERMAL_DIFFUSIVITY', 'ELASTIC_MODULUS', 'ATM_PRES', 'ATM_TEMP', \
			'LOWER_WALL_THICKNESS', 'UPPER_WALL_THICKNESS');

			hdl = config['DV_LIST'].strip('()');
			hdl = hdl.split(",");

			dv_keys_size = len(hdl)/2;
			
			nozzle.DV_Tags = []; # tags: e.g. wall, tstag etc.
			nozzle.DV_Head = []; # correspondance btw DV_Tags and DV_List
			nozzle.wall_dv = []; # 
			
			for i in range(0,2*dv_keys_size,2):
				key = hdl[i].strip();
				NbrDV = int(hdl[i+1]);
								
				if key == 'WALL' :
					nozzle.DV_Tags.append(key);
					nozzle.DV_Head.append(NbrDVTot);

					if 'GEOM_WALL_COEFS_DV' in config:

						count_wall_dv = 0;

						wall_dv = config['GEOM_WALL_COEFS_DV'].strip('()');
						wall_dv = wall_dv.split(",");

						wall_dv_size = len(wall_dv);
						
						if  wall_dv_size != coefs_size  :
							sys.stderr.write('  ## ERROR : Inconsistent number of design variables for the inner wall: %d coefs provided, %d design variables.\n',coefs_size,  wall_dv_size);
							sys.stderr.exit(0);
						
						tag = np.zeros(NbrDV);
						for iw in range(0,wall_dv_size):

							val = int(wall_dv[iw]);

							if ( val > NbrDV or val < 0 ):
								sys.stderr.write('  ## ERROR : Inconsistent design variables were provided for the inner wall (idx=%d).\n' % val);
								sys.exit(0);

							nozzle.wall_dv.append(val-1);

							if val-1 >= 0 :
								tag[val-1] = 1;
							
						for iw in range(0,NbrDV):
							if tag[iw] != 1 :
								sys.stderr.write('  ## ERROR : Design variable %d of inner wall is not defined.\n' % (iw+1));
								sys.exit(0);
					
					else :
						sys.stderr.write('  ## ERROR : Expected GEOM_WALL_COEFS_DV.\n');
						sys.exit(0);
				
				elif ( key == 'LOWER_WALL_THICKNESS' ):
					
					if NbrDV != len(nozzle.lower_wall_thickness) :
						sys.stderr.write('  ## ERROR : Inconsistent number of DV for lower wall thickness definition.\n');
						sys.exit(0);
					
					nozzle.DV_Tags.append(key);
					nozzle.DV_Head.append(NbrDVTot);

				elif ( key == 'UPPER_WALL_THICKNESS' ):
					
					if NbrDV != len(nozzle.upper_wall_thickness) :
						sys.stderr.write('  ## ERROR : Inconsistent number of DV for upper wall thickness definition.\n');
						sys.exit(0);
					
					nozzle.DV_Tags.append(key);
					nozzle.DV_Head.append(NbrDVTot);
				
				elif (
				key == 'INLET_TSTAG' or key == 'INLET_PSTAG' 
				or key == 'THERMAL_CONDUCTIVITY' or key == 'THERMAL_DIFFUSIVITY' 
				or key == 'ELASTIC_MODULUS' or key == 'ATM_PRES'
				or key == 'ATM_TEMP'
				):
					
					if NbrDV != 1 :
						sys.stderr.write('  ## ERROR : Only one design variable expected for %s (%d given).\n' % (key, NbrDV));
						sys.exit(0);
					
					nozzle.DV_Tags.append(key);
					nozzle.DV_Head.append(NbrDVTot);
					
				elif key == 'INLET_PSTAG':
					
					nozzle.DV_Tags.append(key);
					nozzle.DV_Head.append(NbrDVTot);
					
				
										
				else :
					str = '';
					for k in dv_keys :
						str = "%s %s " % (str,k);
					sys.stderr.write('  ## ERROR : Unknown design variable key : %s\n' % key);
					sys.stderr.write('             Expected = %s\n\n' % str);
					sys.exit(0);
					
				NbrDVTot = NbrDVTot + NbrDV;
				
			# for i in keys
		
		else :
			sys.stdout.write('\n  -- Info : No design variable set was defined. Running the baseline parameters.\n\n');
		
		
		nozzle.NbrDVTot = NbrDVTot;
		
		if NbrDVTot > 0 :
			nozzle.DV_Head.append(NbrDVTot);

	def SetupFidelityLevels (self, config, flevel, output='verbose'):
		
		nozzle = self;
		
		fidelity_tags = config['FIDELITY_LEVELS_TAGS'].strip('()');
		fidelity_tags = fidelity_tags.split(",");

		NbrFidLev = len(fidelity_tags);

		if NbrFidLev < 1 :
			print "  ## ERROR : No fidelity level was defined.\n"
			sys.exit(0)

		if output == 'verbose':
		  sys.stdout.write('\n%d fidelity level(s) defined. Summary :\n' % (NbrFidLev));

		  sys.stdout.write('-' * 90);
		  sys.stdout.write('\n%s | %s | %s\n' % ("Level #".ljust(10), "Tag".ljust(10),"Description.".ljust(70)));
		  sys.stdout.write('-' * 90);
		  sys.stdout.write('\n'); 
		elif output == 'quiet':
		  pass
		else:
		  raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')

		for i in range(NbrFidLev) :
			tag = fidelity_tags[i];
			kwd = "DEF_%s" % tag;

			if kwd not in config :
				sys.stderr.write("\n  ## ERROR : The fidelity level tagged %s is not defined.\n" % kwd);
				sys.exit(0);

			cfgLvl = config[kwd].strip('()');
			cfgLvl = cfgLvl.split(",");

			method = cfgLvl[0];

			description = "";

			if i == flevel :
				nozzle.method = method;

			if method == 'NONIDEALNOZZLE' :

				tol = float(cfgLvl[1]);
				if tol < 1e-30 :
					sys.stderr.write("  ## ERROR : Wrong tolerance for fidelity level %d (tagged %s)\n" % (i,tag));
					sys.exit(0);

				if i == flevel :
					nozzle.tolerance = tolerance.Tolerance();
					nozzle.tolerance.setRelTol(tol);
					nozzle.tolerance.setAbsTol(tol);
					#nozzle.tolerance.exitTempPercentError = tol;
				description = "ODE solver relative and absolute tolerance set to %le." % (tol);	

			elif method == 'RANS' or method == 'EULER':
				dim = cfgLvl[1];
				if dim != '2D' and dim != '3D' :
					sys.stderr.write("  ## ERROR : Wrong dimension for fidelity level %d (tagged %s) : only 2D or 3D simulations" % (i,tag));
					sys.exit(0);
				
				nozzle.Dim = dim;
				
				meshsize = cfgLvl[2];	
				if meshsize != 'COARSE' and meshsize != 'MEDIUM' and meshsize != 'FINE' :
					sys.stderr.write("  ## ERROR : Wrong mesh level for fidelity level %d (tagged %s) : must be set to either COARSE, MEDIUM or FINE" % (i,tag));
					sys.exit(0);
				description = "%s %s CFD simulation using the %s mesh level." % (dim, method, meshsize);

				if i == flevel :
					nozzle.meshsize = meshsize;
					
					nozzle.bl_ds        = 0.000007;
					nozzle.bl_ratio     = 1.3; 
					nozzle.bl_thickness = 0.1;
					
					if meshsize == 'COARSE':
						nozzle.meshhl = [0.1, 0.07, 0.06, 0.006, 0.0108];
					elif meshsize == 'MEDIUM':
						nozzle.meshhl = [0.1, 0.07, 0.06, 0.006, 0.0108];
					elif meshsize == 'FINE':
						nozzle.meshhl = [0.1, 0.07, 0.06, 0.006, 0.0108];

			else :
				sys.stderr.write("  ## ERROR : Unknown governing method (%s) for fidelity level %s." % (method, tag));
				sys.stderr.write("  Note: it must be either NONIDEALNOZZLE, EULER, or RANS\n");
				sys.exit(0);

			if output == 'verbose':
				sys.stdout.write("   %s | %s | %s \n" % ( ("%d" % i).ljust(7), tag.ljust(10), textwrap.fill(description,70, subsequent_indent="".ljust(26))) );
			elif output == 'quiet':
				pass
			else:
				raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')

		if output == 'verbose':
			sys.stdout.write('-' * 90);
			sys.stdout.write('\n\n');
		elif output == 'quiet':
			pass
		else:
			raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')

		if flevel >= NbrFidLev :
			sys.stderr.write("\n ## ERROR : Level %d not defined !! \n\n" % flevel);
			sys.exit(0);

		if output == 'verbose':
			sys.stdout.write('  -- Info : Fidelity level to be run : %d\n\n' % flevel);
		elif output == 'quiet':
			pass
		else:
			raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')
		
	def SetupMission(self, config):
	
		nozzle = self;
	
		mission_id = int(config["MISSION"]);

		if(mission_id == 1): # static sea-level thrust case
			altitude = 0.
			mach     = 0.01
			inletTs  = 888.3658
			inletPs  = 3.0550e5
		elif(mission_id == 2): # intermediate case
			altitude = 15000.
			mach     = 0.5
			inletTs  = 942.9857
			inletPs  = 2.3227e5
		elif(mission_id == 3): # high speed, high altitude case
			altitude = 35000.
			mach     = 0.9
			inletTs  = 1021.5
			inletPs  = 1.44925e5
		elif(mission_id == 4): # case with shock in nozzle
			altitude = 0.
			mach     = 0.01
			inletTs  = 900.
			inletPs  = 1.3e5
		elif(mission_id == 5): # subsonic flow
			altitude = 0.
			mach     = 0.01
			inletTs  = 900.
			inletPs  = 1.1e5
		else : 
			sys.stderr.write("\n ## ERROR : UNKNOWN MISSION ID %d !! \n\n" % mission);
			sys.exit(0);

		nozzle.mission = mission.Mission(mission_id);
		nozzle.mission.setMach(mach);
		nozzle.inlet   = inlet.Inlet(inletPs,inletTs);

		hInf = 25; # W/m^2/K, heat transfer coeff. from external wall to env.
		nozzle.environment = environment.Environment(altitude,hInf);

		# --- Setup fluid

		heatRatio = 1.4;
		gasCst    = 287.06;
		nozzle.fluid = fluid.Fluid(heatRatio, gasCst);
		
		# --- Setup convergence parameter
		
		if 'SU2_CONVERGENCE_ORDER' in config:
			nozzle.su2_convergence_order = config['SU2_CONVERGENCE_ORDER'];
		else:
			nozzle.su2_convergence_order = 3;
		
	def SetupBSplineCoefs(self, config):
		
		nozzle = self;
		
		wall_keys = ('GEOM_WALL_PARAM','GEOM_WALL_KNOTS','GEOM_WALL_COEFS');

		if all (key in config for key in wall_keys):

			# --- Get coefs

			hdl = config['GEOM_WALL_COEFS'].strip('()');
			hdl = hdl.split(",");

			coefs_size = len(hdl);

			coefs = [];

			for i in range(0,coefs_size):
				coefs.append(float(hdl[i]));

			# --- Get knots

			hdl = config['GEOM_WALL_KNOTS'].strip('()');
			hdl = hdl.split(",");

			knots_size = len(hdl);

			knots = [];

			for i in range(0,knots_size):
				knots.append(float(hdl[i]));


			#print "%d coefs, %d knots\n" % (coefs_size, knots_size);

		else:
			str = '';
			for key in wall_keys :
				str = "%s %s " % (str,key);
			sys.stderr.write('\n  ## ERROR : NO INNER WALL DEFINITION IS PROVIDED.\n');
			sys.stderr.write('             Expected = %s\n\n' % str);
			sys.exit(0);
		
		nozzle.coefs =  coefs;
		nozzle.knots =  knots;
		nozzle.coefs_size = coefs_size;
		
		thermalConductivity   = 8.6;    # W/m*K, thermal conductivity of wall
		coeffThermalExpansion = 2.3e-6; # 1/K, coeff. of thermal expansion
		elasticModulus        = 80e9;   # Pa, elastic modulus
		poissonRatio          = 0.3;    # Poisson's ratio
		
		nozzle.wall = component.AxisymmetricWall();
		nozzle.wall.material = material.Material(thermalConductivity,coeffThermalExpansion,elasticModulus,poissonRatio);
	
	def ParseThickness(self, config, key):
		nozzle = self;
		
		loc_name = '%s_WALL_THICKNESS_LOCATIONS' % key;
		val_name = '%s_WALL_THICKNESS_VALUES' % key;
		
		wall_keys = (loc_name, val_name);
		
		if all (key in config for key in wall_keys):
    	
			hdl = config[loc_name].strip('()');
			hdl = hdl.split(",");
  
			size_loc = len(hdl);
  		
			wall_thickness = [[0 for i in range(2)] for j in range(size_loc)];
						
			for i in range(0,size_loc):
				wall_thickness[i][0] = float(hdl[i]);
				if wall_thickness[i][0] < 0.0 or wall_thickness[i][0] > 1.0:
					sys.stderr.write('\n ## ERROR : Invalid wall thickness definition (%s).\n' % key);
					sys.stderr.write('              All values must be between 0 and 1.\n\n');
					sys.exit(0);
			
			hdl = config[val_name].strip('()');
			hdl = hdl.split(",");
			size_val = len(hdl);
			
			if size_val != size_loc :
				sys.stderr.write('\n ## ERROR : Inconsistent wall thickness definition (%s).\n' % key);
				sys.stderr.write('              Same number of locations and values required.\n\n');
				sys.exit(0);
			
			for i in range(0,size_val):
				wall_thickness[i][1] = float(hdl[i]);
			
			if wall_thickness[0][0] != 0.0 or wall_thickness[size_val-1][0] != 1.0:
				sys.stderr.write('\n ## ERROR : Invalid wall thickness definition (%s).\n' % key);
				sys.stderr.write('              First and last loc values must be 0 and 1 resp.\n\n');
				sys.exit(0);
				
			return wall_thickness;
		
		else:
			raise; 
		
	
	def SetupWallThickness(self, config):
  
		nozzle = self;
  	
		nozzle.upper_wall_thickness = [[0.0,0.01], [1.0, 0.01]];
		nozzle.lower_wall_thickness = [[0.0,0.01], [1.0, 0.01]];

		try :
			nozzle.lower_wall_thickness = nozzle.ParseThickness(config, "LOWER");
			#print lower_wall_thickness
		except:
			nozzle.lower_wall_thickness = [[0.0,0.01], [1.0, 0.01]];
			
		#print nozzle.lower_wall_thickness;
		
		try:
			nozzle.upper_wall_thickness = nozzle.ParseThickness(config, "UPPER");
		except:
			nozzle.lower_wall_thickness = [[0.0,0.01], [1.0, 0.01]];
		#print nozzle.upper_wall_thickness
		

	def SetupInnerWall (self, config):
		
		nozzle = self;
			
		coefs = nozzle.coefs;
		knots = nozzle.knots;
		
		coefs_size = nozzle.coefs_size;
		
		
		nozzle.height = coefs[coefs_size-1];
		nozzle.length = coefs[coefs_size/2-1];
		

		if nozzle.method == 'RANS' or nozzle.method == 'EULER':
			x    = [];
			y    = [];

			nx = 100;
			_meshutils_module.py_BSplineGeo3 (knots, coefs, x, y, nx);

			nozzle.xwall = x;
			nozzle.ywall = y;
			
		#coefs = np.array(([0.0000, 0.0000, 0.1500, 0.1700, 
	  #    0.1900, 0.2124, 0.2269, 0.2734, 0.3218, 0.3218, 0.3230, 0.3343, 0.3474, 
	  #    0.4392, 0.4828, 0.5673, 0.6700, 0.6700],[0.3255, 0.3255, 0.3255, 
	  #    0.3255, 0.3255, 0.3238, 0.2981, 0.2817, 0.2787, 0.2787, 0.2787, 0.2797, 
	  #    0.2807, 0.2936, 0.2978, 0.3049, 0.3048, 0.3048]));
    
		coefsnp = np.empty([2, coefs_size/2]);
    
		for i in range(0, coefs_size/2):
			coefsnp[0][i] = coefs[i];
			coefsnp[1][i] = coefs[i+coefs_size/2];
    
		#thicknessNodeArray = np.array(([0., 0.33, 0.67],[0.01, 0.01, 0.01]));
		
		
		#nozzle.wall.geometry = geometry.Bspline
		nozzle.wall.geometry = geometry.Bspline(coefsnp);
		
		#print thicknessNodeArray;
		
		# --- LOWER WALL THICKNESS
		
		size_thickness = len(nozzle.lower_wall_thickness);
		thicknessNodeArray = np.zeros(shape=(2,size_thickness))
		
		for i in range(size_thickness):
			thicknessNodeArray[0][i] = nozzle.lower_wall_thickness[i][0]*nozzle.length;
			thicknessNodeArray[1][i] = nozzle.lower_wall_thickness[i][1];
		
		nozzle.wall.lower_thickness = geometry.PiecewiseLinear(thicknessNodeArray);

		# --- UPPER WALL THICKNESS
		
		size_thickness = len(nozzle.upper_wall_thickness);
		thicknessNodeArray = np.zeros(shape=(2,size_thickness))
		
		for i in range(size_thickness):
			thicknessNodeArray[0][i] = nozzle.upper_wall_thickness[i][0]*nozzle.length;
			thicknessNodeArray[1][i] = nozzle.upper_wall_thickness[i][1];
		
		nozzle.wall.upper_thickness = geometry.PiecewiseLinear(thicknessNodeArray);
		
		
		nozzle.wall.thickness = nozzle.wall.lower_thickness;
			

		#sys.exit(1);
		#	
		#if nozzle.method == 'RANS' or nozzle.method == 'EULER':
		#	
		#	coefs = [0.0000, 0.0000, 0.1500, 0.1700, 0.1900, 0.2124, 0.2269, 0.2734, 0.3218, 0.3218, 0.3230, 0.3343, 0.3474, 0.4392, 0.4828, 0.5673, 0.6700, 0.6700, 0.3255, 0.3255, 0.3255, 0.3255, 0.3255, 0.3238, 0.2981, 0.2817, 0.2787, 0.2787, 0.2787, 0.2797, 0.2807, 0.2936, 0.2978, 0.3049, 0.3048, 0.3048];
		#	knots = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 15.0, 15.0, 15.0];
	  #
		#	x    = [];
		#	y    = [];
		#	
		#	nx = 100;
		#	_meshutils_module.py_BSplineGeo3 (knots, coefs, x, y, nx);
		#	
		#	nozzle.xwall = x;
		#	nozzle.ywall = y;
		#	
		#else:
		#	coefs = np.array(([0.0000, 0.0000, 0.1500, 0.1700, 
	  #	    0.1900, 0.2124, 0.2269, 0.2734, 0.3218, 0.3218, 0.3230, 0.3343, 0.3474, 
	  #	    0.4392, 0.4828, 0.5673, 0.6700, 0.6700],[0.3255, 0.3255, 0.3255, 
	  #	    0.3255, 0.3255, 0.3238, 0.2981, 0.2817, 0.2787, 0.2787, 0.2787, 0.2797, 
	  #	    0.2807, 0.2936, 0.2978, 0.3049, 0.3048, 0.3048]));
	  #	
		#	thicknessNodeArray = np.array(([0., 0.33, 0.67],[0.01, 0.01, 0.01]));
	  #	
		#	#nozzle.wall.geometry = geometry.Bspline
		#	nozzle.wall.geometry = geometry.Bspline(coefs);
		#	nozzle.wall.thickness = geometry.PiecewiseLinear(thicknessNodeArray);
		#	
		#	thermalConductivity   = 8.6;    # W/m*K, thermal conductivity of wall
		#	coeffThermalExpansion = 2.3e-6; # 1/K, coeff. of thermal expansion
		#	elasticModulus        = 80e9;   # Pa, elastic modulus
		#	poissonRatio          = 0.3;    # Poisson's ratio
		#	
		#	nozzle.wall.material = material.Material(thermalConductivity,coeffThermalExpansion,elasticModulus,poissonRatio);

	def ParseDV (self, config):
		
		nozzle = self;
		
		if 'INPUT_DV_FORMAT' in config:
			format = config['INPUT_DV_FORMAT'];
		else :
			sys.stderr.write("  ## ERROR : Input DV file format not specified. (INPUT_DV_FORMAT expected: PLAIN or DAKOTA)\n");
			sys.exit(0);
		
				
		if 'INPUT_DV_NAME' in config:
			filename = config['INPUT_DV_NAME'];
		else :
			sys.stderr.write("  ## ERROR : Input DV file name not specified. (INPUT_DV_NAME expected)\n");
			sys.exit(0);
				
		if format == 'PLAIN':
			DV_List, OutputCode, Derivatives_DV = ParseDesignVariables_Plain(filename);	
			NbrDV = len(DV_List);		
			
		elif format == 'DAKOTA' :
			DV_List, OutputCode, Derivatives_DV = ParseDesignVariables_Dakota(filename);	
			NbrDV = len(DV_List);
		else:
			sys.stderr.write("  ## ERROR : Unknown DV input file format %s\n" % format);
			sys.exit(0);
		
		if NbrDV != nozzle.NbrDVTot : 
			sys.stderr.write("  ## Error : Inconsistent number of design variables are given in %s\n" % filename);
			sys.stderr.write("             %d given, %d expected\n" % (NbrDV, nozzle.NbrDVTot ));
			sys.exit(0);
	
		nozzle.DV_List = DV_List;
	
	def UpdateDV(self, config, output='verbose'):
		
		nozzle = self;
		
		NbrTags = len(nozzle.DV_Tags);
		
		prt_name = [];
		prt_basval = [];
		prt_newval = [];	
		
		NbrChanged = 0; # Count the total number of changed parameters
		                # Note: different from the number of DV, because one DV might correspond to more BSP coefs
		
		for iTag in range(0,NbrTags):
			Tag = nozzle.DV_Tags[iTag];
			NbrDV = nozzle.DV_Head[iTag+1] - nozzle.DV_Head[iTag]; 
			
			if Tag == 'WALL':
				
				for iCoef in range(0,len(nozzle.wall_dv)):
					id_dv = nozzle.DV_Head[iTag] + nozzle.wall_dv[iCoef];
					
					# --- Update coef iCoef if required
					if id_dv >= 0 :
						#sys.stdout.write("Update Bspline coef %d: %lf -> %lf\n" % (iCoef, nozzle.coefs[iCoef], nozzle.DV_List[id_dv]));
						prt_name.append('Bspline coef #%d' % (iCoef+1));
						prt_basval.append('%.4lf'% nozzle.coefs[iCoef]);
						prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
						nozzle.coefs[iCoef] = nozzle.DV_List[id_dv];
						NbrChanged = NbrChanged+1;
				
			elif Tag == 'INLET_TSTAG':
				
				id_dv = nozzle.DV_Head[iTag];
				
				#sys.stdout.write("Update inlet temperature: %lf -> %lf\n" % (nozzle.inlet.Tstag,nozzle.DV_List[id_dv]));
				
				prt_name.append('Inlet temperature');
				prt_basval.append('%.2lf'% nozzle.inlet.Tstag);
				prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
				
				nozzle.inlet.Tstag = nozzle.DV_List[id_dv];
				
			elif Tag == 'INLET_PSTAG':
				
				id_dv = nozzle.DV_Head[iTag];
				#sys.stdout.write("Update inlet pressure: %lf -> %lf\n" % (nozzle.inlet.Pstag, nozzle.DV_List[id_dv]));
				
				prt_name.append('Inlet pressure');
				prt_basval.append('%.2lf'% nozzle.inlet.Pstag);
				prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
				
				nozzle.inlet.Pstag = nozzle.DV_List[id_dv];
			
			elif Tag == 'THERMAL_CONDUCTIVITY':
				
				id_dv = nozzle.DV_Head[iTag];
				#sys.stdout.write("Update thermal conductivity: %lf -> %lf\n" % (nozzle.wall.material.k, nozzle.DV_List[id_dv]));
				
				prt_name.append('Thermal conductivity');
				prt_basval.append('%.2lf'% nozzle.wall.material.k);
				prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
				
				nozzle.wall.material.k = nozzle.DV_List[id_dv];
				
			elif Tag == 'THERMAL_DIFFUSIVITY':
				
				id_dv = nozzle.DV_Head[iTag];
				#sys.stdout.write("Update thermal diffusivity: %lf -> %lf\n" % (nozzle.wall.material.k, nozzle.DV_List[id_dv]));
				
				prt_name.append('Thermal diffusivity');
				prt_basval.append('%.2le'% nozzle.wall.material.alpha);
				prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
				
				nozzle.wall.material.alpha = nozzle.DV_List[id_dv];
			
			elif Tag == 'ELASTIC_MODULUS':
				
				id_dv = nozzle.DV_Head[iTag];
				#sys.stdout.write("Update elastic modulus: %lf -> %lf\n" % (nozzle.wall.material.E, nozzle.DV_List[id_dv]));
				
				prt_name.append('Elastic modulus');
				prt_basval.append('%.2le'% nozzle.wall.material.E);
				prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
				
				nozzle.wall.material.E = nozzle.DV_List[id_dv];
			
			elif Tag == 'ATM_PRES':
				
				id_dv = nozzle.DV_Head[iTag];
				#sys.stdout.write("Update atmospheric pressure: %lf -> %lf\n" % (nozzle.environment.P, nozzle.DV_List[id_dv]));
				
				prt_name.append('Atmospheric pressure');
				prt_basval.append('%.2lf'% nozzle.environment.P);
				prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
				
				nozzle.environment.P = nozzle.DV_List[id_dv];
				
			elif Tag == 'ATM_TEMP':
				
				id_dv = nozzle.DV_Head[iTag];
				#sys.stdout.write("Update atmospheric temperature: %lf -> %lf\n" % (nozzle.environment.T, nozzle.DV_List[id_dv]));
				
				prt_name.append('Atmospheric temperature');
				prt_basval.append('%.2lf'% nozzle.environment.T);
				prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
				
				nozzle.environment.T = nozzle.DV_List[id_dv];
				
			elif Tag == 'LOWER_WALL_THICKNESS':
				
				for i in range(NbrDV):
					id_dv = nozzle.DV_Head[iTag]+i;
					
					prt_name.append('Lower wall thickness t=%.3lf' % nozzle.lower_wall_thickness[i][0]);
					prt_basval.append('%.4lf'% nozzle.lower_wall_thickness[i][1]);
					prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
					
					nozzle.lower_wall_thickness[i][1] = nozzle.DV_List[id_dv];
					
			elif Tag == 'UPPER_WALL_THICKNESS':
      
				for i in range(NbrDV):
					id_dv = nozzle.DV_Head[iTag]+i;
      
					prt_name.append('Upper wall thickness t=%.3lf' % nozzle.upper_wall_thickness[i][0]);
					prt_basval.append('%.4lf'% nozzle.upper_wall_thickness[i][1]);
					prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
      
					nozzle.upper_wall_thickness[i][1] = nozzle.DV_List[id_dv];
			
			else :
				sys.stderr.write("  ## Error : Unknown tag name %s\n" % Tag);
				sys.exit(0);
				
			if Tag != 'WALL':
				NbrChanged = NbrChanged+1;
		
		# --- Print summary
		if output == 'verbose':
		  sys.stdout.write('\n%d parameter(s) updated according to %d design variable(s). Summary:\n' % (NbrChanged, nozzle.NbrDVTot));
      
		  sys.stdout.write('-' * 79);
		  sys.stdout.write('\n%s | %s | %s\n' % ("DV name".ljust(30), "Baseline value".ljust(20),"Updated value".ljust(20)));
		  sys.stdout.write('-' * 79);
		  sys.stdout.write('\n');
		  for i in range(0,len(prt_name)):
			  sys.stdout.write('%s | %s | %s\n' % (prt_name[i].ljust(30), prt_basval[i].ljust(20),prt_newval[i].ljust(20)));
		  sys.stdout.write('-' * 79);	
		  sys.stdout.write('\n\n');
		elif output == 'quiet':
		  pass
		else:
		  raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')
		
	def SetupOutputFunctions (self, config):
		
		nozzle = self;
		
		nozzle.Output_Tags = [];
		
		nozzle.Output_Volume = 0 ;
		nozzle.Output_Thrust = 0 ;
		
		nozzle.thermal_stress    = 0.0;
		nozzle.mechanical_stress = 0.0;
		
		nozzle.GetOutput = dict();
		
		nozzle.GetOutput['VOLUME'] = 0;
		nozzle.GetOutput['THRUST'] = 0;
		
		if 'OUTPUT_NAME' in config:
			nozzle.Output_Name = config['OUTPUT_NAME'];
		else :
			sys.stderr.write("  ## ERROR : Output function file name not specified in config file. (OUTPUT_NAME expected)\n");
			sys.exit(0);
		
				
		# --- Initialize outputs
		nozzle.Thrust = -1;
		nozzle.Volume = -1;
		
		if 'OUTPUT_FUNCTIONS' in config:

			dv_keys = ('VOLUME', 'THRUST');

			hdl = config['OUTPUT_FUNCTIONS'].strip('()');
			hdl = hdl.split(",");
			
			for i in range(0,len(hdl)):
				
				key = hdl[i].strip();
						
				if (
				key == 'VOLUME' or key == 'THRUST' or 'MECHANICAL_STRESS' or 'THERMAL_STRESS'
				):
					nozzle.GetOutput[key] = 1;
					nozzle.Output_Tags.append(key);
				
				else :
					str = '';
					for k in dv_keys :
						str = "%s %s " % (str,k);
					sys.stderr.write('  ## ERROR : Unknown output function name : %s\n' % key);
					sys.stderr.write('             Expected = %s\n\n' % str);
					sys.exit(0);
				
			# for i in keys
		
		if len(nozzle.Output_Tags) == 0 :
			sys.stderr.write("  ## Error : No output function was given.\n\n");
			sys.exit(0);


	def WriteOutputFunctions_Plain (self, output='verbose'):
  	
		nozzle = self;

		filename = nozzle.Output_Name;
		
		if output == 'verbose':
			sys.stdout.write('\n');
			str = " Post-processing ";
			nch = (60-len(str))/2;
			sys.stdout.write('-' * nch);
			sys.stdout.write(str);
			sys.stdout.write('-' * nch);
			sys.stdout.write('\n\n');
		elif output == 'quiet':
			pass
		else:
			raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')
		
		try:
			fil = open(filename, 'w');
		except:
			sys.stderr.write("  ## ERROR : Could not open output file %s\n" % filename);
			sys.exit(0);
  
		if output == 'verbose':
			sys.stdout.write('  -- Info : Output functions file : %s\n' % filename);
		elif output == 'quiet':
			pass
		else:
			raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')

		for i in range(0, len(nozzle.Output_Tags)):
  		
			tag = nozzle.Output_Tags[i];
  		
			if tag == 'THRUST':
				fil.write('%0.16f\n' % nozzle.Thrust);
				if output == 'verbose':
					sys.stdout.write('      Thrust = %0.16f\n' % nozzle.Thrust);
  
			if tag == 'VOLUME':
				fil.write('%0.16f\n' % nozzle.Volume);
				if output == 'verbose':
					sys.stdout.write('      Volume = %0.16f\n' % nozzle.Volume);
					
			if tag == 'MECHANICAL_STRESS':
				fil.write('%0.16f mechanical_stress\n' % nozzle.mechanical_stress);
			
			if tag == 'THERMAL_STRESS':
				fil.write('%0.16f thermal_stress\n' % nozzle.thermal_stress);					
  	
		if output == 'verbose':
			sys.stdout.write('\n');
		fil.close();
		
	def WriteOutputFunctions_Dakota (self):
	
		nozzle = self;
	
		filename = nozzle.Output_Name;
	
		try:
			fil = open(filename, 'w');
		except:
			sys.stderr.write("  ## ERROR : Could not open output file %s\n" % filename);
			sys.exit(0);
	
		sys.stdout.write('  -- Info : Output functions file : %s\n' % filename);
	
		for i in range(0, len(nozzle.Output_Tags)):
	
			tag = nozzle.Output_Tags[i];
	
			if tag == 'THRUST':
				fil.write('%0.16f thrust\n' % nozzle.Thrust);
				#sys.stdout.write('%0.16f thrust\n' % nozzle.Thrust);
	
			if tag == 'VOLUME':
				fil.write('%0.16f volume\n' % nozzle.Volume);
				#sys.stdout.write(' %0.16f \n' % nozzle.Volume);
			
			if tag == 'MECHANICAL_STRESS':
				fil.write('%0.16f mechanical_stress\n' % nozzle.mechanical_stress);
			
			if tag == 'THERMAL_STRESS':
				fil.write('%0.16f thermal_stress\n' % nozzle.thermal_stress);

		sys.stdout.write('\n');
		fil.close();
		

def NozzleSetup( config_name, flevel, output='verbose' ):
	import tempfile
		
	if not os.path.isfile(config_name) :
		sys.stderr.write("  ## ERROR : could not find configuration file %s\n\ns" % config_name);
		sys.exit(0);
	
	nozzle = multif.nozzle.nozzle.Nozzle()
	
	config = SU2.io.Config(config_name)
	
	# --- File names
	
 	#nozzle.mesh_name    =	 'blabla.su2'
 	#nozzle.restart_name =  'blabla.dat'

	#hdl, toto = tempfile.mkstemp(suffix='.su2');
 	#
	#print "TMPNAME = %s" %  (toto);

	nozzle.mesh_name    =  'nozzle.su2'; #tempfile.mkstemp(suffix='.su2');
	nozzle.restart_name =  'nozzle.dat'; #tempfile.mkstemp(suffix='.dat');
	
	if 'TEMP_RUN_DIR' in config:
		if config['TEMP_RUN_DIR'] == 'YES':
			nozzle.runDir       =  tempfile.mkdtemp();	
		else:
			nozzle.runDir = '';
	else:
		nozzle.runDir = '';
	
	# --- Path to SU2 exe
	
	if 'SU2_RUN' in config:
		nozzle.SU2_RUN = config['SU2_RUN'];
	else:
		nozzle.SU2_RUN = os.environ['SU2_RUN'];
	
	# --- Parse fidelity levels
	
	nozzle.SetupFidelityLevels(config, flevel, output);
	
	# --- Set flight regime + fluid
	
	nozzle.SetupMission(config);
	
	# --- Setup Bspline coefficients 
	
	nozzle.SetupBSplineCoefs(config);
	
	# --- Setup wall thickness
	
	nozzle.SetupWallThickness(config);
	
	# --- Setup DV definition : 
		
	nozzle.SetupDV(config);
	
	# --- If input DV are expected, parse them and update nozzle
	
	if nozzle.NbrDVTot > 0 :	
		
		# Parse DV from input DV file (plain or dakota format)
		nozzle.ParseDV(config);
		
		# Update DV using values provided in input DV file
		nozzle.UpdateDV(config,output);
	
	# --- Compute inner wall's bspline from coefs and knots
	#     It is done now because the coefs might have been updated according to the input DV file
	
	nozzle.SetupInnerWall(config);
	
	# --- Get output functions to be returned
	
	nozzle.SetupOutputFunctions(config);
		
	#sys.exit(1);
	return nozzle;	
