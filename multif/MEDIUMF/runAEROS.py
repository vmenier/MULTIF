# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import multif
import ctypes
import numpy as np
import pylab
from .. import _nozzle_module

from postprocessing import *

def runAEROS ( nozzle ):
	
	# --- Get the CFD solution at the inner wall
	# SolExtract : Solution (x, y, sol1, sol2, etc.)
	# Size : [NbrVer, SolSiz]
	# idHeader : id of each solution field, e.g. mach_ver89 = SolExtract[89][idHeader['MACH']]
	
	SolExtract, Size, idHeader  = ExtractSolutionAtWall(nozzle);
	
	# --- Material properties
	
	#self.k = k # W/m*K, thermal conductivity of wall
	#self.alpha = alpha # 1/K, coeff. of thermal expansion
	#self.E = E # Pa, elastic modulus
	#self.v = v # Poisson's ratio
	
	Ek   = nozzle.wall.material.E;     # (float > 0) Young’s modulus of k­th material.                    
	Nuk  = nozzle.wall.material.v;     # (float < 0.5) Poisson’s ratio of k­th material.                  
	Rhok = 0.0;                        # (float > 0) Density of k­th material (mass per unit volume).     -> Victorien : ?
	Tk   = 0.005;                      # (float > 0) Shell thickness of k­th material.                    -> Victorien : ?
	Ak   = nozzle.wall.material.alpha; # (float >= 0) Coefficient of thermal expansion of k­th material.  
	Rk   = 0.0;                        # (float > 0) Reference temperature of k­th material.              -> Victorien : ?
	
	Lbd  = nozzle.wall.material.k;     # W/m*K, thermal conductivity of wall                              
	
	iPres = idHeader['Pressure'];
	iTemp = idHeader['Temperature'];
	
	#print "\n\n#####################################";
	#print "###  INTERFACE WITH AERO-S GOES HERE" ;
	#print "#####################################\n\n";
	
	## --- How to get x, y, P, T :
	#for i in range(0,Size[0]):
	#	print "VER %d : (x,y) = (%lf, %lf) , Pres = %lf, Temp = %lf" % (i, SolExtract[i][0], SolExtract[i][1], SolExtract[i][iPres], SolExtract[i][iTemp]);
		
	f1 = open("NOZZLE.txt", 'w');
	print >> f1, "%d %d %d %f" % (Size[0], 2, 1, 0.02);
	for i in range(0,Size[0]):
		print >> f1, "%lf %lf" % (SolExtract[i][0], SolExtract[i][1]);
	print >> f1, "%d 0.0 -1" % (0);
	print >> f1, "%d 0.0 -1" % (Size[0]-1);
	print >> f1, "%d 0 0.0 -1" %(0);
	print >> f1, "%lf %lf %lf %lf %lf %lf" % (Ek, Nuk, Rhok, Tk, Ak, Rk);
	f1.close();

	f2 = open("BOUNDARY.txt", 'w');
	print >> f2, "%d" % (Size[0]);
	for i in range(0,Size[0]):
		print >> f2, "%lf %lf %lf" % (SolExtract[i][0], SolExtract[i][iPres], SolExtract[i][iTemp]);
	f2.close();

	_nozzle_module.generate();
	
	os.system("aeros nozzle.aeros")
	
	AEROSPostProcessing ( nozzle );
	
def AEROSPostProcessing ( nozzle ):
	
	# --- Open MECHANICAL_STRESS
	
	try:
		fil = open("MECHANICAL_STRESS", "r" );
	except IOError:
		sys.stderr.write('## ERROR : UNABLE TO OPEN MECHANICAL_STRESS FILE. RETURN 0.\n');
		nozzle.mechanical_stress = 0;
		return;
	
	lines = [line.split() for line in fil];
	
	max_mech = 0.0;
	for i in range(2,len(lines)):
		max_mech = max(float(lines[i][0]), max_mech);
		
	nozzle.mechanical_stress = max_mech;
	
	fil.close();
	
	# --- Open THERMAL_STRESS
	
	try:
		fil = open("THERMAL_STRESS", "r" );
	except IOError:
		sys.stderr.write('## ERROR : UNABLE TO OPEN THERMAL_STRESS FILE. RETURN 0.\n');
		nozzle.thermal_stress = 0;
		return;
	
	lines = [line.split() for line in fil];
	
	max_therm = 0.0;
	for i in range(2,len(lines)):
		max_therm = max(float(lines[i][0]), max_therm);
	
	nozzle.thermal_stress = max_therm;
	
	fil.close();
	
	
	
	
	
