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
    
    # Start of example for accessing nozzle properties
    
    print '\nEntered runAEROS\n'
    
    x = np.linspace(0,nozzle.length,10)
    
    # First print all information related to geometry
    print '--- Wall:'
    print 'B-spline coefs (x-r coordinates): {}'.format(nozzle.wall.coefs)
    print 'B-spline knots: {}'.format(nozzle.wall.knots)
    print ''
    for i in range(len(nozzle.wall.layer)):
        print '--- %s:' % nozzle.wall.layer[i].name
        print 'material: %s' % nozzle.wall.layer[i].material.name
        print 'thickness node x-coordinate: {} m'.format(nozzle.wall.layer[i].thickness.nodes[0,:])
        print 'thickness node local n-coordinate: {} m'.format(nozzle.wall.layer[i].thickness.nodes[1,:])
        print '\n'
    print '--- Baffles:'
    print 'number: %i' % nozzle.baffles.n
    print 'location (x-coordinate): {} m'.format(nozzle.baffles.location)
    print 'thickness: {} m'.format(nozzle.baffles.thickness)
    print 'height: {} m'.format(nozzle.baffles.height)
    print ''
    print '--- Stringers:'
    print 'number: %i' % nozzle.stringers.n
    print 'thickness: {} m'.format(nozzle.stringers.thickness)
    print 'height: {} m'.format(nozzle.stringers.height)
    print ''
    
    # aside: to easily calculate wall shape or thickness as a function of the
    # x-coordinate in an item's local coordinates, use the radius method
    r_innerwall = nozzle.wall.geometry.radius(x) # returns shape of inner wall at x in r-coordinate
    t_thermal_layer = nozzle.wall.geometry.radius(x) # returns thickness of thermal layer at x in local n-coordinate
    
    # Next, print all information related to materials
    for k in nozzle.materials:
        print '--- %s material:' % nozzle.materials[k].name
        print 'type: %s' % nozzle.materials[k].type
        print 'density: %10.4f kg/m^3' % nozzle.materials[k].getDensity()
        if nozzle.materials[k].name == 'TI-HC' or nozzle.materials[k].name == 'GR-BMI':
            print 'elastic modulus: {} Pa'.format(nozzle.materials[k].getElasticModulus())
            print 'shear modulus: %e Pa' % nozzle.materials[k].getShearModulus()
            print 'poisson ratio: %1.4f' % nozzle.materials[k].getPoissonRatio()
        if nozzle.materials[k].name == 'GR-BMI':            
            print 'mutual influence coefs: {}'.format(nozzle.materials[k].getMutualInfluenceCoefs())
        print 'thermal conductivity: {} W/m-K'.format(nozzle.materials[k].getThermalConductivity())
        if nozzle.materials[k].name == 'TI-HC' or nozzle.materials[k].name == 'GR-BMI':
            print 'thermal expansion coef: {} 1/K'.format(nozzle.materials[k].getThermalExpansionCoef())
        print ''
        
    # So for example to obtain the Young's Moduli for the innermost Gr-BMI layer
    [E1, E2] = nozzle.wall.layer[1].material.getElasticModulus()
    # Example, obtain all 3 coefs of thermal expansion for outermost Gr-BMI layer
    [alpha1, alpha2, alpha12] = nozzle.wall.layer[3].material.getThermalExpansionCoef()
    # Example, obtain thermal conductivities for innermost Gr-BMI layer (there
    # are 3 values: (1st in-plane direction, 2nd in-plane direction, out-of-plane direction))
    [k1, k2, k3] = nozzle.wall.layer[1].material.getThermalConductivity()
    
    # End of example
    
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

    boundaryFlag = 0; # 0: inlet fixed, 1: baffles fixed, 2: both inlet and baffles fixed
    thermalFlag = 1;  # 0: structural analysis only, 1: both thermal and structural analyses
    
    ## --- How to get x, y, P, T :
    #for i in range(0,Size[0]):
    #    print "VER %d : (x,y) = (%lf, %lf) , Pres = %lf, Temp = %lf" % (i, SolExtract[i][0], SolExtract[i][1], SolExtract[i][iPres], SolExtract[i][iTemp]);
    
    f1 = open("NOZZLE.txt", 'w');
    print >> f1, "%d %d %d %f %d %d" % (Size[0], 2, 2, 0.02, boundaryFlag, thermalFlag);
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

    if thermalFlag > 0:
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
    
    
    
    
    
