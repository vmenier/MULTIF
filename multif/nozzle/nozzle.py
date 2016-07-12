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
import pylab

class Nozzle:
	def __init__(self):
		pass


def NozzleSetup( config_name, flevel ):
	
	if not os.path.isfile(config_name) :
		msg = "  ## ERROR : could not find configuration file %s\n\ns" % config_name;
		sys.stdout.write(msg);
		sys.exit(1);
	
	nozzle = multif.nozzle.nozzle.Nozzle()
	
	config = SU2.io.Config(config_name)
	
	fidelity_tags = config['FIDELITY_LEVELS_TAGS'].strip('()');
	fidelity_tags = fidelity_tags.split(",");
	
	NbrFidLev = len(fidelity_tags);
	
	if NbrFidLev < 1 :
		print "  ## ERROR : No fidelity level was defined.\n"
		sys.exit(1)
	
	sys.stdout.write('\n%d fidelity level(s) defined. Summary :\n' % (NbrFidLev));
	
	sys.stdout.write('-' * 90);
	sys.stdout.write('\n%s | %s | %s\n' % ("Level #".ljust(10), "Tag".ljust(10),"Description.".ljust(70)));
	sys.stdout.write('-' * 90);
	sys.stdout.write('\n'); 
	
	for i in range(NbrFidLev) :
		tag = fidelity_tags[i];
		kwd = "DEF_%s" % tag;
	
		if kwd not in config :
			print "\n  ## ERROR : The fidelity level tagged %s is not defined.\n" % kwd;
			sys.exit(1)
	
		cfgLvl = config[kwd].strip('()');
		cfgLvl = cfgLvl.split(",");
	
		method = cfgLvl[0];
	
		description = "";
		
		if i == flevel :
			nozzle.method = method;
	
		if method == 'NONIDEALNOZZLE' :
						
			tol = float(cfgLvl[1]);
			if tol < 1e-30 :
				print "  ## ERROR : Wrong tolerance for fidelity level %d (tagged %s)\n" % (i,tag);
			description = "Quasi 1D non ideal nozzle with tolerance set to %lf." % (tol);	
			
			if i == flevel :
				nozzle.tolerance = tolerance.Tolerance();
				nozzle.tolerance.setTol(tol);
			
		elif method == 'RANS' or method == 'EULER':
			dim = cfgLvl[1];
			if dim != '2D' and dim != '3D' :
				print "  ## ERROR : Wrong dimension for fidelity level %d (tagged %s) : only 2D or 3D simulations" % (i,tag);
			meshsize = cfgLvl[2];	
			if meshsize != 'COARSE' and meshsize != 'MEDIUM' and meshsize != 'FINE' :
				print "  ## ERROR : Wrong mesh level for fidelity level %d (tagged %s) : must be set to either COARSE, MEDIUM or FINE" % (i,tag);
			description = "%s %s CFD simulation using the %s mesh level." % (dim, method, meshsize);
			
			if i == flevel :
				nozzle.meshsize = meshsize;
				if meshsize == 'COARSE':
					nozzle.meshhl = [0.1, 0.07, 0.06, 0.006, 0.0108];
				elif meshsize == 'MEDIUM':
					nozzle.meshhl = [0.1, 0.07, 0.06, 0.006, 0.0108];
				elif meshsize == 'FINE':
					nozzle.meshhl = [0.1, 0.07, 0.06, 0.006, 0.0108];
		
		else :
			print "  ## ERROR : Unknown governing method (%s) for fidelity level %s." % (method, tag);
			print "  Note: it must be either NONIDEALNOZZLE, EULER, or RANS\n";
	
			
		sys.stdout.write("   %s | %s | %s \n" % ( ("%d" % i).ljust(7), tag.ljust(10), textwrap.fill(description,70, subsequent_indent="".ljust(26))) );
	
	sys.stdout.write('-' * 90);
	sys.stdout.write('\n\n');
	
	if flevel >= NbrFidLev :
		sys.stderr.write("\n ## ERROR : Level %d not defined !! \n\n" % flevel);
		sys.exit(0);
	
	sys.stdout.write('  -- Info : Fidelity level to be run : %d\n\n' % flevel);
	
	
	# --- Set flight regime
	
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
	
	# --- Setup nozzle inner wall
	
	nozzle.wall = component.AxisymmetricWall();
	
	if nozzle.method == 'RANS' or nozzle.method == 'EULER':
		
		coefs = [0.0000, 0.0000, 0.1500, 0.1700, 0.1900, 0.2124, 0.2269, 0.2734, 0.3218, 0.3218, 0.3230, 0.3343, 0.3474, 0.4392, 0.4828, 0.5673, 0.6700, 0.6700, 0.3255, 0.3255, 0.3255, 0.3255, 0.3255, 0.3238, 0.2981, 0.2817, 0.2787, 0.2787, 0.2787, 0.2797, 0.2807, 0.2936, 0.2978, 0.3049, 0.3048, 0.3048];
		knots = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 15.0, 15.0, 15.0];

		x    = [];
		y    = [];
		
		nx = 100;
		_meshutils_module.py_BSplineGeo3 (knots, coefs, x, y, nx);
		
		nozzle.xwall = x;
		nozzle.ywall = y;
		
	else:
		coefs = np.array(([0.0000, 0.0000, 0.1500, 0.1700, 
  	    0.1900, 0.2124, 0.2269, 0.2734, 0.3218, 0.3218, 0.3230, 0.3343, 0.3474, 
  	    0.4392, 0.4828, 0.5673, 0.6700, 0.6700],[0.3255, 0.3255, 0.3255, 
  	    0.3255, 0.3255, 0.3238, 0.2981, 0.2817, 0.2787, 0.2787, 0.2787, 0.2797, 
  	    0.2807, 0.2936, 0.2978, 0.3049, 0.3048, 0.3048]));
  	
		thicknessNodeArray = np.array(([0., 0.33, 0.67],[0.01, 0.01, 0.01]));
  	
		#nozzle.wall.geometry = geometry.Bspline
		nozzle.wall.geometry = geometry.Bspline(coefs);
		nozzle.wall.thickness = geometry.PiecewiseLinear(thicknessNodeArray);
		
		thermalConductivity   = 8.6;    # W/m*K, thermal conductivity of wall
		coeffThermalExpansion = 2.3e-6; # 1/K, coeff. of thermal expansion
		elasticModulus        = 80e9;   # Pa, elastic modulus
		poissonRatio          = 0.3;    # Poisson's ratio
		
		nozzle.wall.material = material.Material(thermalConductivity,coeffThermalExpansion,elasticModulus,poissonRatio);
		
		
		
	return nozzle;
	