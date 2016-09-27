# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import multif
import ctypes
import numpy as np
from .. import _nozzle_module

from postprocessing import *

def runAEROS ( nozzle ):
    
    # --- Get the CFD solution at the inner wall
    # SolExtract : Solution (x, y, sol1, sol2, etc.)
    # Size : [NbrVer, SolSiz]
    # idHeader : id of each solution field, e.g. mach_ver89 = SolExtract[89][idHeader['MACH']]
    
    SolExtract, Size, idHeader  = ExtractSolutionAtWall(nozzle);
    
    # --- Wall thicknesses
    
    xtab = np.linspace(0, nozzle.length, num=10);
    for i in range(len(xtab)):
        x = xtab[i];
        hl = nozzle.wall.thermal_layer.thickness.radius(x);
        hu = nozzle.wall.load_layer.thickness.radius(x);
        print "x = %lf lower thickness = %lf upper thickness = %lf " % (x, hl, hu);
    
    # --- Material properties
    
    # Thermal layer
    k_axial = nozzle.wall.thermal_layer.material.getThermalConductivity('axial');
    k_radial = nozzle.wall.thermal_layer.material.getThermalConductivity('radial');
    Lbd_therm = (k_axial + k_radial)/2.;
    
    # Load layer
    
    # Young’s modulus of k­th material in axial direction. (float > 0)
    E_axial = nozzle.wall.load_layer.material.getElasticModulus('axial');
    E_radial = nozzle.wall.load_layer.material.getElasticModulus('radial');
    Ek   = (E_axial + E_radial)/2.;     
    # Poisson's ratio of kth material. (float < 0.5)
    Nuk  = nozzle.wall.load_layer.material.getPoissonRatio();          
    # Density of kth material. (float > 0)
    Rhok = nozzle.wall.load_layer.material.getDensity();
    # Shell thickness of kth material. (float > 0)             -> Victorien : ?
    Tk   = 0.005;
    # Coefficient of thermal expansion of k­th material. (float >= 0)  
    A_axial = nozzle.wall.load_layer.material.getThermalExpansionCoef('axial');
    A_radial = nozzle.wall.load_layer.material.getThermalExpansionCoef('radial');
    Ak   = (A_axial + A_radial)/2.; 
    # Reference temperature of k­th material. (float > 0) 
    Rk   = nozzle.environment.T;
    # Thermal conductivity of wall (W/m*K)
    k_axial = nozzle.wall.load_layer.material.getThermalConductivity('axial');
    k_radial = nozzle.wall.load_layer.material.getThermalConductivity('radial');
    Lbd  = (k_axial + k_radial)/2.;     
    
    # Heat transfer coefficient to environment (W/m^2-K)         
    Hk   = nozzle.environment.hInf;    
    
    iPres = idHeader['Pressure'];
    iTemp = idHeader['Temperature'];

    # --- Mesh parameters
    Nn   = 21; # Number of nodes in longitudinal direction
    Mn   = 21; # Number of nodes in circumferential direction 
    Tn   = 4;  # Number of nodes through thickness of thermal insulating layer
    
    ## --- How to get x, y, P, T :
    #for i in range(0,Size[0]):
    #    print "VER %d : (x,y) = (%lf, %lf) , Pres = %lf, Temp = %lf" % (i, SolExtract[i][0], SolExtract[i][1], SolExtract[i][iPres], SolExtract[i][iTemp]);
    
    f1 = open("NOZZLE.txt", 'w');
    print >> f1, "%d %d %d %f" % (Size[0], 2, 2, 0.02);
    for i in range(0,Size[0]):
        print >> f1, "%lf %lf" % (SolExtract[i][0], SolExtract[i][1]);
    print >> f1, "%d 0.0 -1 0 %lf %lf" % (0, nozzle.wall.load_layer.thickness.radius(0), nozzle.wall.thermal_layer.thickness.radius(0));
    print >> f1, "%d 0.0 -1 0 %lf %lf" % (Size[0]-1, nozzle.wall.load_layer.thickness.radius(nozzle.length), nozzle.wall.thermal_layer.thickness.radius(nozzle.length));
    print >> f1, "0 2 0.0 -1 %d %d 0 1 %d" % (Nn, Mn, Tn);
    print >> f1, "%lf %lf %lf %lf %lf %lf %lf" % (Ek, Nuk, Rhok, Tk, Ak, Lbd, Hk); # material properties for the upper layer
    print >> f1, "0 0 0 0 0 %lf 0" % (Lbd_therm); # material properties for the lower layer
    f1.close();

    f2 = open("BOUNDARY.txt", 'w');
    print >> f2, "%d" % (Size[0]);
    for i in range(0,Size[0]):
        print >> f2, "%lf %lf %lf %lf" % (SolExtract[i][0], SolExtract[i][iPres], SolExtract[i][iTemp], Rk); # XXX
    f2.close();
	
	## --- How to get x, y, P, T :
	#for i in range(0,Size[0]):
	#	print "VER %d : (x,y) = (%lf, %lf) , Pres = %lf, Temp = %lf" % (i, SolExtract[i][0], SolExtract[i][1], SolExtract[i][iPres], SolExtract[i][iTemp]);

    _nozzle_module.generate();       # generate the meshes for thermal and structural analyses

    os.system("aeros nozzle.aeroh"); # execute the thermal analysis
    _nozzle_module.convert();        # convert temperature output from thermal analysis to input for structural analysis
    os.system("aeros nozzle.aeros"); # execute the structural analysis
    
    AEROSPostProcessing ( nozzle );
    
def AEROSPostProcessing ( nozzle ):
    
        
    # --- Open MECHANICAL_STRESS
    
    try:
        fil = open("MECHANICAL_STRESS", "r" );
    except IOError:
        sys.stderr.write('\n ## ERROR : UNABLE TO OPEN MECHANICAL_STRESS FILE. RETURN 0.\n\n');
        nozzle.max_mechanical_stress = 0;
        return;
    
    lines = [line.split() for line in fil];
    
    max_mech = 0.0;
    for i in range(2,len(lines)):
        max_mech = max(float(lines[i][0]), max_mech);
        
    nozzle.max_mechanical_stress = max_mech;
    
    fil.close();
    
    # --- Open THERMAL_STRESS
    
    try:
        fil = open("THERMAL_STRESS", "r" );
    except IOError:
        sys.stderr.write('\n ## ERROR : UNABLE TO OPEN THERMAL_STRESS FILE. RETURN 0.\n\n');
        nozzle.max_thermal_stress = 0;
        return;
    
    lines = [line.split() for line in fil];
    
    max_therm = 0.0;
    for i in range(2,len(lines)):
        max_therm = max(float(lines[i][0]), max_therm);
    
    nozzle.max_thermal_stress = max_therm;
    
    fil.close();
    
    
    
    
    
