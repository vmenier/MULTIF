# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import multif
import ctypes
import numpy as np
import pylab

from postprocessing import *

def runAEROS ( nozzle ):
	
	# --- Get the CFD solution at the inner wall
	# SolExtract : Solution (x, y, sol1, sol2, etc.)
	# Size : [NbrVer, SolSiz]
	# idHeader : id of each solution field, e.g. mach_ver89 = SolExtract[89][idHeader['MACH']]
	
	SolExtract, Size, idHeader  = ExtractSolutionAtWall(nozzle);
	
	iPres = idHeader['Pressure'];
	iTemp = idHeader['Temperature'];
	
	print "\n\nINTERFACE WITH AEROS GOES HERE\n\n" ;
	
	# --- How to get x, y, P, T :
	#for i in range(0,10):
	#	print "VER %d : (x,y) = (%lf, %lf) , Pres = %lf, Temp = %lf" % (i, SolExtract[i][0], SolExtract[i][1], SolExtract[i][iPres], SolExtract[i][iTemp]);
		
		
	
	
	
	
	
		
	
	