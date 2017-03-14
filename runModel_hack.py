#!/usr/bin/env python

import os, time, sys, shutil, copy

from optparse import OptionParser
import textwrap

import multif
from multif import SU2

from scipy.interpolate import splev, splrep
from scipy import interpolate as itp

import HIGHF
import numpy as np

from multif.nozzle import geometry

def UpdateCrossSections(DVin):
	
	update_r1 = np.loadtxt("r1dv.dat");
	update_r2 = np.loadtxt("r2dv.dat");
	
	inpData = np.loadtxt('inputsection.dat');
	
	NbrSec = len(inpData);
	
	x  = inpData[:,0]
	r1 = inpData[:,1]
	r2 = inpData[:,2]
	
	#x = np.array([ 0.,   0.4,  0.8,  1.2,  1.6,  2.,   2.4,  2.8,  3.2,  3.6,  4. ]);
	#r1 = np.zeros(len(x));
	#r2 = np.zeros(len(x));
	
	# --- R1
	
	for i in range(len(update_r1)):
		id_r = update_r1[i][0];
		id_dv = update_r1[i][1];
		print "r1 %d (DV %d): %lf -> %lf" % (id_r, id_dv, r1[id_r], DVin[id_dv])
		r1[id_r] = DVin[id_dv];
		
	# --- R2
	
	for i in range(len(update_r2)):
		id_r  = update_r2[i][0];
		id_dv = update_r2[i][1];	
		print "r2 %d (DV %d): %lf -> %lf" % (id_r, id_dv, r2[id_r], DVin[id_dv])
		r2[id_r] = DVin[id_dv];
		
		
	for i in range(len(x)):
		if (x[i] <= 2.0 + 1e-12 ):
			r2[i] = r1[i];
	
	return [x, r1, r2]
	
	
	


def GenCAD(idSam, hdlSam, sizes, params):
	
	cadNam	= "cad_%d.geo"%idSam;
	cadNamAxi = "cad_axi_%d.geo"%idSam;
	
	mshNam	= "cad_%d.su2"%idSam;
	mshNamAxi = "cad_axi_%d.su2"%idSam;
	
	# --- Get DVin
	
	DVin = hdlSam[idSam];
	
	# --- Parse centerline bspline
	
	inpBsp = np.loadtxt('inputbspline.dat');
	k = 3;
	Nbc = len(inpBsp)/2;
	[t,c] = np.split(inpBsp,2)
	tck = [t,c,3];
	
	# --- Update centerline coefs

	#tckc =  UpdateCenterlineDV(tck, DVin);
	tckc = tck;
	def fz(x):
		return splev(x, tckc);
	
	# --- Update cross sections
	
	[x, r1, r2] = UpdateCrossSections(DVin)
	
	print x
	print r1 
	print r2
	
	tck1 = splrep(x, r1, xb=x[0], xe=x[-1], k=3, s=0);
	def fr1(x):
		return splev(x, tck1);
		
	tck2 = splrep(x, r2, xb=x[0], xe=x[-1], k=3, s=0);
	def fr2(x):
		return splev(x, tck2);
	
	
	#HIGHF.HF_DefineCAD(cadNam,x,fr1,fr2, fz, sizes, params);	
	#
	#HIGHF.MF_DefineAxiSymCAD ("axi.geo", x, fr1, fr2, fz, sizes, params);
	
	#gmsh_executable = 'gmsh';
	#try :
	#	cmd = [gmsh_executable, '-3', cadNam, '-o', mshNam];
	#	out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
	#except:
	#	raise;
	
	NbrLnk = 5000;
	
	xthroat = params[3];
	xexit   = params[4];
	xmax    = params[5];
	ymax    = params[6];
	zmax    = params[7];
	
	hl0 = sizes[3];  # Edge size in the 'small' plume area
	hl1 = sizes[4];  # Edge size in the vicinity of the aircraft
	
	xx = np.linspace(0, xexit, NbrLnk);
	yy = HIGHF.MF_GetRadius (xx, x, fr1, fr2, fz, params);
	
	return [xx,yy];

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

	output = 'verbose';
	
	if output == 'verbose':		 
		sys.stdout.write('-' * 90);
		sys.stdout.write('\n');
		sys.stdout.write('\t __  __ _   _ _  _____ ___		 ___ \t\n') ;
		sys.stdout.write('\t|  \/  | | | | ||_   _|_ _|  ___  | __|\t\n');
		sys.stdout.write('\t| |\/| | |_| | |__| |  | |  |___| | _| \t\t Dev. : R. Fenrich & V. Menier\n');
		sys.stdout.write('\t|_|  |_|\___/|____|_| |___|	   |_|  \t\t		Stanford University\n\n');
		sys.stdout.write('-' * 90);
		sys.stdout.write('\n\n');
	
	# Command Line Options
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="read config from FILE", metavar="FILE")
	parser.add_option("-n", "--partitions", dest="partitions", default=0,
					  help="number of PARTITIONS", metavar="PARTITIONS")
	
	parser.add_option("-l", "--flevel", dest="flevel", default=-1,
					  help="fidelity level to run", metavar="FLEVEL")	
	
	(options, args)=parser.parse_args()
	
	options.partitions = int( options.partitions )
	options.flevel	 = int( options.flevel )
	
	if options.flevel < 0 :
		sys.stderr.write("  ## ERROR : Please choose a fidelity level to run (option -l or --flevel)");
		sys.exit(0);
	
	nozzle = multif.nozzle.NozzleSetup( options.filename, options.flevel, output, options.partitions );

	
	# --- Update geometry according to sample file
	
	idSam = 2;
	
	FilNam = "samples.dat"
	hdlSam = np.loadtxt(FilNam);
	
	idSam = 6;
	
	params = np.zeros(100)
	params[0] = -0.0739843; # z crd of the cut at the throat
	params[1] = -0.458624;  # z crd of the flat exit (bottom of the shovel)
	params[2] = -0.15;	  # z coordinate of the top of the shovel
	params[3] = 2.0;		# x_throat 
	params[4] = 4.0;		# x_exit 
	#params[5] = 10;   		  # xmax
	params[5] = 7;   		  # xmax
	params[6] = 8;			  # ymax
	params[7] = 10;   		  # zmax
	params[8] = 0.04;	   # thickness of shovel + exit

	# --- Size parameters : 
	# 		- inside nozzle
	# 		- surface aircraft
	# 		- farfield
	# 		- (box) plume region
	# 		- (box) aircraft vicinity
	#sizes = [0.02, 0.1, 1,  0.03, 0.08 ];
	#sizes = [0.035, 0.1, 1,  0.06, 0.1 ];
	sizes = [0.035, 0.2, 1,  0.09, 0.15 ];
	
	
	[x,y] = GenCAD(idSam, hdlSam, sizes, params)
	
	radiusTab = np.zeros(shape=(2,len(x)));
	for j in range(len(x)):
		radiusTab[0][j] = x[j];
		radiusTab[1][j] = y[j];
	
	# --- plot geometry before/after
	
	#import matplotlib.pyplot as plt
	#
	#plt.plot(x, y, label="NEW GEOM")
	#plt.legend()
	#plt.show()
	
	

	nozzle.wall.geometry = geometry.PiecewiseLinear(radiusTab);
	
	print nozzle.wall.geometry.length
	
   # if nozzle.method == 'NONIDEALNOZZLE' :
   #	 #multif.LOWF.Run(nozzle, output);
   # elif nozzle.method == 'EULER' or nozzle.method == 'RANS':
   #	 #multif.MEDIUMF.Run(nozzle, output);
	
	
	
	

	if nozzle.method == 'NONIDEALNOZZLE' :
		multif.LOWF.Run(nozzle, output);
	elif nozzle.method == 'EULER' or nozzle.method == 'RANS':
		multif.MEDIUMF.Run(nozzle, output);
	
	# --- Print warning in case the wrong SU2 version was run
	if nozzle.method != 'NONIDEALNOZZLE' and nozzle.SU2Version != 'OK':
		sys.stdout.write('\n');
		sys.stdout.write('#' * 90);
		sys.stdout.write('\n  ## WARNING : You are not using the right version of SU2. This may have caused robustness issues.\n');
		sys.stdout.write('#' * 90);
		sys.stdout.write('\n\n');

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
	main()
