# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import util

def GenerateNozzleMesh ():
	# 
	# 
	print "TOTO";
	
def NozzleGeoFile( xwall, ywall):
	
	# Domain definition
	
	nx = len(xwall);
	
	length = xwall[nx-1];
	
	CrdBox = [[0 for x in range(2)] for y in range(9)] 
	
	CrdBox[1][0] = 0;          CrdBox[1][1] = 0;
	CrdBox[2][0] = length;     CrdBox[2][1] = 0;
	CrdBox[3][0] = 1.5;        CrdBox[3][1] = 0;
	CrdBox[4][0] = 1.5;        CrdBox[4][1] = 2.5;
	CrdBox[5][0] = -0.67;      CrdBox[5][1] = 2.5;
	CrdBox[6][0] = -0.67;      CrdBox[6][1] = 0.4244;
	CrdBox[7][0] = 0.1548;     CrdBox[7][1] = 0.4244;
	CrdBox[8][0] = length;     CrdBox[8][1] = ywall[nx-1]+0.012;
		
	fil = open('nozzle.geo', 'w');
	
	fil.write('toto');
	
	fil.close();
	

