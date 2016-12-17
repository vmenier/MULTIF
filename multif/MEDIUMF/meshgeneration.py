# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy, math
from optparse import OptionParser
import textwrap
import ctypes
import numpy as np

import multif




class optionsmesh:
	def __init__(self):
		pass


def GenerateNozzleMesh (nozzle):
	import tempfile
	
	#hdl, nozzle.tmpGeoNam = tempfile.mkstemp(suffix='.geo');
	#hdl, nozzle.tmpMshNam = tempfile.mkstemp(suffix='.mesh');
	
	nozzle.tmpGeoNam = 'nozzle_tmp.geo';
	nozzle.tmpMshNam = 'nozzle_tmp.mesh';
	
	# --- Write geo file
	
	mesh_options = optionsmesh();
	mesh_options.xwall  = nozzle.xwall;
	mesh_options.ywall  = nozzle.ywall;
	mesh_options.hl     = nozzle.meshhl; 
	mesh_options.method = nozzle.method; # Euler or RANS
	
	mesh_options.ds          = nozzle.bl_ds;
	mesh_options.ratio       = nozzle.bl_ratio 
	mesh_options.thickness   = nozzle.bl_thickness; 
	
	mesh_options.x_thrust = nozzle.x_thrust;  	
	
	NozzleGeoFile(nozzle.tmpGeoNam, mesh_options);
	
	# --- Call Gmsh
	
	try :   
		CallGmsh(nozzle);
	except:
		print "\n  ## ERROR : Mesh generation failed.\n";
		sys.exit(0);
	
	# --- Mesh preprocessing
	try :
		MeshPrepro(nozzle);
	except :
		print "\n  ## ERROR : Mesh preprocessing failed.\n";
		sys.exit(0);
	

def CallGmsh (nozzle):
	import subprocess
	gmsh_executable = 'gmsh';
	try :
		cmd = [gmsh_executable, '-2', nozzle.tmpGeoNam, '-o', nozzle.tmpMshNam];
		out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
	except:
		raise;
	
def MeshPrepro (nozzle):
	from .. import _meshutils_module
	
	out = _meshutils_module.py_MeshPrepro2D (nozzle.tmpMshNam, nozzle.mesh_name);
	
	if ( out == 0 ) :
		raise;
	
def NozzleGeoFile_old(FilNam, Mesh_options):
	
	# --- Options
		
	xwall  = Mesh_options.xwall;
	ywall  = Mesh_options.ywall;
 	hl     = Mesh_options.hl;
	method = Mesh_options.method;
	
	
	ds        =  Mesh_options.ds;       
	ratio     =  Mesh_options.ratio;   
	thickness =  Mesh_options.thickness;
	
	
	# --- Domain definition
	
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
	
	try:
		fil = open(FilNam, 'w');
	except:
		sys.stderr.write("  ## ERROR : Could not open %s\n" % FilNam);
		sys.exit(0);
	
	sys.stdout.write("%s OPENED.\n" % FilNam);
	
	# --- These sizes are not in use anymore -> See background field definitions
	sizWal = 0.3;
	sizFar = 0.3;
	sizSym = 0.3;
	
	fil.write('Point(1) = {%lf, %lf, 0, %lf};\n' % (CrdBox[1][0], CrdBox[1][1], sizWal));
	fil.write('Point(2) = {%lf, %lf, 0, %lf};\n' % (CrdBox[2][0], CrdBox[2][1], sizWal));
	fil.write('Point(3) = {%lf, %lf, 0, %lf};\n' % (CrdBox[3][0], CrdBox[3][1], sizSym));
	fil.write('Point(4) = {%lf, %lf, 0, %lf};\n' % (CrdBox[4][0], CrdBox[4][1], sizFar));
	fil.write('Point(5) = {%lf, %lf, 0, %lf};\n' % (CrdBox[5][0], CrdBox[5][1], sizFar));
	fil.write('Point(6) = {%lf, %lf, 0, %lf};\n' % (CrdBox[6][0], CrdBox[6][1], sizWal));
	fil.write('Point(7) = {%lf, %lf, 0, %lf};\n' % (CrdBox[7][0], CrdBox[7][1], sizWal));
	fil.write('Point(8) = {%lf, %lf, 0, %lf};\n' % (CrdBox[8][0], CrdBox[8][1], sizWal));
	fil.write('Point(9)  = {%lf, %lf, 0, %lf};\n'% (CrdBox[3][0], CrdBox[8][1], sizWal));
	fil.write('Point(10) = {%lf, %lf, 0, %lf};\n'% (CrdBox[3][0], CrdBox[7][1]+0.25*CrdBox[8][1], sizWal));
	fil.write('Point(11) = {%lf, %lf, 0, %lf};\n'% (CrdBox[6][0], CrdBox[7][1]+0.25*CrdBox[8][1], sizWal));
	
	# --- Add B-Spline control points
	
	vid = 11;
	for i in range(0,nx) :
		vid=vid+1;
		fil.write('Point(%d)  = {%lf, %lf, 0, %lf};\n'% (vid, xwall[nx-i-1], ywall[nx-i-1], sizWal));
	
	# --- Add domain lines
	
	fil.write('Line(1)  = {1, 2};\n');
	fil.write('Line(2)  = {2, 3};\n');
	fil.write('Line(3)  = {3, 9};\n');
	fil.write('Line(4)  = {9, 10};\n');
	fil.write('Line(5)  = {10, 4};\n');
	fil.write('Line(6)  = {4, 5};\n');
	fil.write('Line(7)  = {5, 11};\n');
	fil.write('Line(8)  = {11, 6};\n');
	fil.write('Line(9)  = {6, 7};\n');
	fil.write('Line(10) = {7, 8};\n');
	fil.write('Line(11) = {8, 12};\n');
	
	# --- Define B-Spline
	
	fil.write('BSpline(12) = { 12');
	eid=13;
	
	for i in range(eid, eid+nx-1):
		fil.write(', %d' % i);
	
	fil.write('};\n');
	
	fil.write('Line(13) = {%d, 1};\n' % (i));
	
	# --- Plane surface
	
	fil.write('Line Loop(14) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};\n');
	fil.write('Plane Surface(14) = {14};\n');
		
	hl1 = hl[0];
	hl2 = hl[1];
	hl3 = hl[2];
	hl4 = hl[3];
	hl5 = hl[4];
	
	#print "HL = %lf %lf %lf %lf %lf\n" % (hl1, hl2, hl3, hl4, hl5);
	
	fields = [[-100, 100, -10, 0.8,    hl3], \
	          [-100, 100, -10, 0.5,    hl1],\
	          [-100, 100, -10, 0.5,    hl5],\
	          [-100, 0.3, -10, 0.37,   hl4],\
						[-100, 100, -10, 0.8,    hl2],\
	          [-100, length+0.04, -10, 0.5, hl4]];
	
	NbrFld = len(fields);
	
	for i in range(0,NbrFld):
		fil.write('Field[%d] = Box;\n'     % (i+1             ));
		fil.write('Field[%d].VIn = %lf;\n' % (i+1,fields[i][4]));
		fil.write('Field[%d].VOut = %lf;\n'% (i+1,sizFar      ));
		fil.write('Field[%d].XMin = %lf;\n'% (i+1,fields[i][0]));
		fil.write('Field[%d].XMax = %lf;\n'% (i+1,fields[i][1]));
		fil.write('Field[%d].YMin = %lf;\n'% (i+1,fields[i][2]));
		fil.write('Field[%d].YMax = %lf;\n'% (i+1,fields[i][3]));	
	
	fil.write('Field[%d] = Min;\n'        % (NbrFld+1));
	fil.write('Field[%d].FieldsList = { 1'% (NbrFld+1));	
	
	for i in range(2,NbrFld+1) :
		fil.write(', %d ' % i);
	
	fil.write('};\n');
	
	fil.write('Background Field = %d;\n' % (NbrFld+1));
	
	# --- Add boundary layer definition
	if method == 'RANS':
		fil.write('Field[%d] = BoundaryLayer;          \n' % (NbrFld+2));
		fil.write('Field[%d].EdgesList = {9,10,11,12}; \n' % (NbrFld+2));
		fil.write('Field[%d].NodesList = {6,111};      \n' % (NbrFld+2));
		fil.write('Field[%d].hfar = 1;                 \n' % (NbrFld+2));
		fil.write('Field[%d].hwall_n = %le;       \n' % ((NbrFld+2), ds));
		fil.write('Field[%d].hwall_t = 0.300000;       \n' % (NbrFld+2));
		fil.write('Field[%d].ratio = %le;              \n' % ((NbrFld+2),ratio));
		fil.write('Field[%d].thickness = %le;         \n' % ((NbrFld+2), thickness));
		fil.write('BoundaryLayer Field = %d;           \n' % (NbrFld+2));
	
	
	else :
		
		fil.write('Delete {                                             \n');
		fil.write('  Surface{14};                                       \n');
		fil.write('}                                                    \n');
		fil.write('Line(15) = {12, 2};                                  \n');
		fil.write('Line Loop(16) = {10, 11, 15, 2, 3, 4, 5, 6, 7, 8, 9};\n');
		fil.write('Plane Surface(17) = {16};                            \n');
		fil.write('Line Loop(18) = {12, 13, 1, -15};                    \n');
		fil.write('Plane Surface(19) = {18};                            \n');
		fil.write('Physical Surface(20) = {19, 17};                     \n');
		fil.write('Physical Line(1)  = {1};                             \n');
		fil.write('Physical Line(2)  = {2};                             \n');
		fil.write('Physical Line(3)  = {3};                             \n');
		fil.write('Physical Line(4)  = {4};                             \n');
		fil.write('Physical Line(5)  = {5};                             \n');
		fil.write('Physical Line(6)  = {6};                             \n');
		fil.write('Physical Line(7)  = {7};                             \n');
		fil.write('Physical Line(8)  = {8};                             \n');
		fil.write('Physical Line(9)  = {9};                             \n');
		fil.write('Physical Line(10) = {10};                            \n');
		fil.write('Physical Line(11) = {11};                            \n');
		fil.write('Physical Line(12) = {12};                            \n');
		fil.write('Physical Line(13) = {13};                            \n');
		fil.write('Physical Line(14) = {14};                            \n');
	
	
	fil.close();

	
def NozzleGeoFile(FilNam, Mesh_options):

	# --- Options

	xwall  = Mesh_options.xwall;
	ywall  = Mesh_options.ywall;
 	hl     = Mesh_options.hl;
	method = Mesh_options.method;


	ds        =  Mesh_options.ds;       
	ratio     =  Mesh_options.ratio;   
	thickness =  Mesh_options.thickness;

	x_thrust  = Mesh_options.x_thrust;

	# --- Domain definition

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

	try:
		fil = open(FilNam, 'w');
	except:
		sys.stderr.write("  ## ERROR : Could not open %s\n" % FilNam);
		sys.exit(0);

	sys.stdout.write("%s OPENED.\n" % FilNam);

	# --- These sizes are not in use anymore -> See background field definitions
	sizWal = 0.3;
	sizFar = 0.3;
	sizSym = 0.3;
		
	exit_vid = -1;
	
	fil.write('Point(1) = {%lf, %lf, 0, %lf};\n' % (CrdBox[1][0], CrdBox[1][1], sizWal));
	fil.write('Point(2) = {%lf, %lf, 0, %lf};\n' % (x_thrust, CrdBox[2][1], sizWal));
	fil.write('Point(3) = {%lf, %lf, 0, %lf};\n' % (CrdBox[3][0], CrdBox[3][1], sizSym));
	fil.write('Point(4) = {%lf, %lf, 0, %lf};\n' % (CrdBox[4][0], CrdBox[4][1], sizFar));
	fil.write('Point(5) = {%lf, %lf, 0, %lf};\n' % (CrdBox[5][0], CrdBox[5][1], sizFar));
	fil.write('Point(6) = {%lf, %lf, 0, %lf};\n' % (CrdBox[6][0], CrdBox[6][1], sizWal));
	fil.write('Point(7) = {%lf, %lf, 0, %lf};\n' % (CrdBox[7][0], CrdBox[7][1], sizWal));
	fil.write('Point(8) = {%lf, %lf, 0, %lf};\n' % (CrdBox[8][0], CrdBox[8][1], sizWal));
	fil.write('Point(9)  = {%lf, %lf, 0, %lf};\n'% (CrdBox[3][0], CrdBox[8][1], sizWal));
	fil.write('Point(10) = {%lf, %lf, 0, %lf};\n'% (CrdBox[3][0], CrdBox[7][1]+0.25*CrdBox[8][1], sizWal));
	fil.write('Point(11) = {%lf, %lf, 0, %lf};\n'% (CrdBox[6][0], CrdBox[7][1]+0.25*CrdBox[8][1], sizWal));

	# --- Add B-Spline control points

	vid = 11;
	for i in range(0,nx) :
		vid=vid+1;
		fil.write('Point(%d)  = {%lf, %lf, 0, %lf};\n'% (vid, xwall[nx-i-1], ywall[nx-i-1], sizWal));
		
		if ( math.fabs(x_thrust-xwall[nx-i-1]) < 1e-6 ):
			exit_vid = vid;
			print "x %lf exit_vid = %d" % (xwall[nx-i-1], exit_vid)

	if ( exit_vid == -1 ):
		print " ## ERROR : Coordinates for the thrust computation don't match.";
		sys.exit(1);
	
	# --- Add domain lines

	fil.write('Line(1)  = {1, 2};\n');
	fil.write('Line(2)  = {2, 3};\n');
	fil.write('Line(3)  = {3, 9};\n');
	fil.write('Line(4)  = {9, 10};\n');
	fil.write('Line(5)  = {10, 4};\n');
	fil.write('Line(6)  = {4, 5};\n');
	fil.write('Line(7)  = {5, 11};\n');
	fil.write('Line(8)  = {11, 6};\n');
	fil.write('Line(9)  = {6, 7};\n');
	fil.write('Line(10) = {7, 8};\n');
	fil.write('Line(11) = {8, 12};\n');

	# --- Define B-Spline

	fil.write('BSpline(12) = { 12');
	eid=13;

	for i in range(eid, exit_vid+1):
		fil.write(', %d' % i);

	fil.write('};\n');
	
	fil.write('BSpline(15) = { %d'%(exit_vid));
	
	for i in range(exit_vid+1, eid+nx-1):
		fil.write(', %d' % i);

	fil.write('};\n');

	fil.write('Line(13) = {%d, 1};\n' % (i));

	# --- Plane surface

	fil.write('Line Loop(14) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,15, 13};\n');
	fil.write('Plane Surface(14) = {14};\n');

	hl1 = hl[0];
	hl2 = hl[1];
	hl3 = hl[2];
	hl4 = hl[3];
	hl5 = hl[4];

	#print "HL = %lf %lf %lf %lf %lf\n" % (hl1, hl2, hl3, hl4, hl5);

	fields = [[-100, 100, -10, 0.8,    hl3], \
	          [-100, 100, -10, 0.5,    hl1],\
	          [-100, 100, -10, 0.5,    hl5],\
	          [-100, 0.3, -10, 0.37,   hl4],\
			  [-100, 100, -10, 0.8,    hl2],\
	          [-100, length+0.04, -10, 0.5, hl4]];

	NbrFld = len(fields);

	for i in range(0,NbrFld):
		fil.write('Field[%d] = Box;\n'     % (i+1             ));
		fil.write('Field[%d].VIn = %lf;\n' % (i+1,fields[i][4]));
		fil.write('Field[%d].VOut = %lf;\n'% (i+1,sizFar      ));
		fil.write('Field[%d].XMin = %lf;\n'% (i+1,fields[i][0]));
		fil.write('Field[%d].XMax = %lf;\n'% (i+1,fields[i][1]));
		fil.write('Field[%d].YMin = %lf;\n'% (i+1,fields[i][2]));
		fil.write('Field[%d].YMax = %lf;\n'% (i+1,fields[i][3]));	

	fil.write('Field[%d] = Min;\n'        % (NbrFld+1));
	fil.write('Field[%d].FieldsList = { 1'% (NbrFld+1));	

	for i in range(2,NbrFld+1) :
		fil.write(', %d ' % i);

	fil.write('};\n');

	fil.write('Background Field = %d;\n' % (NbrFld+1));

	# --- Add boundary layer definition
	if method == 'RANS':
		fil.write('Field[%d] = BoundaryLayer;          \n' % (NbrFld+2));
		fil.write('Field[%d].EdgesList = {9,10,11,12}; \n' % (NbrFld+2));
		fil.write('Field[%d].NodesList = {6,111};      \n' % (NbrFld+2));
		fil.write('Field[%d].hfar = 1;                 \n' % (NbrFld+2));
		fil.write('Field[%d].hwall_n = %le;       \n' % ((NbrFld+2), ds));
		fil.write('Field[%d].hwall_t = 0.300000;       \n' % (NbrFld+2));
		fil.write('Field[%d].ratio = %le;              \n' % ((NbrFld+2),ratio));
		fil.write('Field[%d].thickness = %le;         \n' % ((NbrFld+2), thickness));
		fil.write('BoundaryLayer Field = %d;           \n' % (NbrFld+2));


	else :
    
		fil.write('Delete {                                                  \n');
		fil.write('  Surface{14};                                            \n');
		fil.write('}                                                         \n');
		fil.write('Line(16) = {%d, 2};                                       \n'%exit_vid);
		fil.write('Line Loop(17) = {16, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}; \n');
		fil.write('Plane Surface(18) = {17};                                 \n');
		fil.write('Line Loop(19) = {1, -16, 15, 13};                         \n');
		fil.write('Plane Surface(20) = {19};                                 \n');
		fil.write('Physical Surface(21) = {20, 18};                          \n');
		
		
		fil.write('Physical Line(1)  = {1};                             \n');
		fil.write('Physical Line(2)  = {2};                             \n');
		fil.write('Physical Line(3)  = {3};                             \n');
		fil.write('Physical Line(4)  = {4};                             \n');
		fil.write('Physical Line(5)  = {5};                             \n');
		fil.write('Physical Line(6)  = {6};                             \n');
		fil.write('Physical Line(7)  = {7};                             \n');
		fil.write('Physical Line(8)  = {8};                             \n');
		fil.write('Physical Line(9)  = {9};                             \n');
		fil.write('Physical Line(10) = {10};                            \n');
		fil.write('Physical Line(11) = {11};                            \n');
		fil.write('Physical Line(12) = {12,15};                         \n');
		fil.write('Physical Line(13) = {13};                            \n');
		fil.write('Physical Line(14) = {14};                            \n');


	fil.close();

	fil = open("%s.opt"%FilNam,'w');
	fil.write("Mesh.SaveElementTagType = 2;\n");
	fil.close();