#!/usr/bin/env python

import os, time, sys, shutil, copy

from optparse import OptionParser
import textwrap

import multif
from multif import SU2

from multif import nozzle

from multif.MEDIUMF import *

from decimal import Decimal

from multif import _meshutils_module

class optionsmesh:
	def __init__(self):
		pass

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
	
	nozzle = multif.nozzle.nozzle.Nozzle()
	config = SU2.io.Config('optim.cfg')
	
	nozzle.SetupBSplineCoefs(config);
	nozzle.SetupDV(config);
	
	if nozzle.NbrDVTot > 0 :	
		
		# Parse DV from input DV file (plain or dakota format)
		nozzle.ParseDV(config);
		
		# Update DV using values provided in input DV file
		nozzle.UpdateDV(config);
	
	nozzle.method = 'RANS';
	
	nozzle.SetupInnerWall(config);
	
	nozzle.meshhl = [0.1, 0.07, 0.06, 0.006, 0.0108];
	
	#nozzle.bl_ds        = Decimal(config.BLDS);
	#nozzle.bl_ratio     = Decimal(config.BLRATIO); 
	#nozzle.bl_thickness = Decimal(config.BLTHICKNESS);
	
	nozzle.bl_ds        = Decimal(config.BLDS);
	nozzle.bl_ratio     = Decimal(config.BLRATIO); 
	nozzle.bl_thickness = Decimal(config.BLTHICKNESS);
	
	print "BOUNDARY LAYER PARAMETERS: ds=%le, ratio=%f, thickness=%f" % (nozzle.bl_ds, nozzle.bl_ratio, nozzle.bl_thickness);
	
	nozzle.mesh_name = 'nozzle.su2';
	
	GenerateNozzleMesh(nozzle);
	
	
	_meshutils_module.py_ConvertSU2toGMSH("nozzle.su2", "", "nozzle")
	
	#mesh_options = optionsmesh();
	#mesh_options.xwall  = nozzle.xwall;
	#mesh_options.ywall  = nozzle.ywall;
	#mesh_options.hl     = nozzle.meshhl; 
	#mesh_options.method = nozzle.method; # Euler or RANS
	
# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
