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
	
	print "\n\n#####################################";
	print "###  INTERFACE WITH AERO-S GOES HERE" ;
	print "#####################################\n\n";
	
	# --- How to get x, y, P, T :
	for i in range(0,Size[0]):
		print "VER %d : (x,y) = (%lf, %lf) , Pres = %lf, Temp = %lf" % (i, SolExtract[i][0], SolExtract[i][1], SolExtract[i][iPres], SolExtract[i][iTemp]);
		
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
