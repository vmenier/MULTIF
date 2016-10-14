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
    
    nozzle.wall.load_layer    = nozzle.wall.layer[1];
    nozzle.wall.thermal_layer = nozzle.wall.layer[0];
    
    print '\nEntered runAEROS\n'
    
    # Start of example for accessing nozzle properties
    
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
    #r_innerwall = nozzle.wall.geometry.radius(x) # returns shape of inner wall at x in r-coordinate
    #t_thermal_layer = nozzle.wall.layer[0].thickness.radius(x) # returns thickness of thermal layer at x in local n-coordinate
    
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
    #[E1, E2] = nozzle.wall.layer[1].material.getElasticModulus()
    # Example, obtain all 3 coefs of thermal expansion for outermost Gr-BMI layer
    #[alpha1, alpha2, alpha12] = nozzle.wall.layer[3].material.getThermalExpansionCoef()
    # Example, obtain thermal conductivities for innermost Gr-BMI layer (there
    # are 3 values: (1st in-plane direction, 2nd in-plane direction, out-of-plane direction))
    #[k1, k2, k3] = nozzle.wall.layer[1].material.getThermalConductivity()
    # End of example
    
    SolExtract, Size, idHeader  = ExtractSolutionAtWall(nozzle);
    
    # Average thermal conductivity of load layer (W/m*K)
    [k11, k12, k13] = nozzle.wall.layer[1].material.getThermalConductivity();
    k2 = nozzle.wall.layer[2].material.getThermalConductivity();
    [k31, k32, k33] = nozzle.wall.layer[3].material.getThermalConductivity();
    Lbd = (k11+k12+k13+3*k2+k31+k32+k33)/9
    
    iPres = idHeader['Pressure'];
    iTemp = idHeader['Temperature'];
    
    # --- Mesh parameters
    Nn   = 21; # Number of nodes in longitudinal direction
    Mn   = 21; # Number of nodes in circumferential direction 
    Tn   = 4;  # Number of nodes through thickness of thermal insulating layer
    
    boundaryFlag = 0; # 0: inlet fixed, 1: baffles fixed, 2: both inlet and baffles fixed
    thermalFlag = 1;  # 0: structural analysis only, 1: both thermal and structural analyses
    
    f1 = open("NOZZLE.txt", 'w');
    print >> f1, "%d %d %d %f %d %d %d" % (Size[0], 2, 4, 0.02, boundaryFlag, thermalFlag, 3);
    for i in range(0,Size[0]):
        print >> f1, "%lf %lf" % (SolExtract[i][0], SolExtract[i][1]);
    print >> f1, "%d 0.0 -1 0 %lf %lf %lf %lf" % (0, nozzle.wall.layer[1].thickness.radius(0),
             nozzle.wall.layer[2].thickness.radius(0), nozzle.wall.layer[3].thickness.radius(0),
             nozzle.wall.layer[0].thickness.radius(0));
    print >> f1, "%d 0.0 -1 0 %lf %lf %lf %lf" % (Size[0]-1, nozzle.wall.layer[1].thickness.radius(nozzle.length),
             nozzle.wall.layer[2].thickness.radius(nozzle.length), nozzle.wall.layer[3].thickness.radius(nozzle.length),
             nozzle.wall.layer[0].thickness.radius(nozzle.length));
    print >> f1, "0 1 2 3 0.0 -1 %d %d 0 1 %d" % (Nn, Mn, Tn);
    # material properties for the load layer
    for i in range(1,4):
      if nozzle.wall.layer[i].material.type == 'ISOTROPIC':
        print >> f1, "ISOTROPIC %lf %lf %lf 0 %lf %lf %lf" % \
                (nozzle.wall.layer[i].material.getElasticModulus(), nozzle.wall.layer[i].material.getPoissonRatio(),
                 nozzle.wall.layer[i].material.getDensity(), nozzle.wall.layer[i].material.getThermalExpansionCoef(),
                 Lbd, nozzle.environment.hInf);
      else:
        [E1, E2] = nozzle.wall.layer[i].material.getElasticModulus()
        [mu1, mu2] = nozzle.wall.layer[i].material.getMutualInfluenceCoefs()
        [alpha1, alpha2, alpha12] = nozzle.wall.layer[i].material.getThermalExpansionCoef()
        print >> f1, "ANISOTROPIC %lf %lf %lf %lf %lf %lf %lf 0 %lf %lf %lf %lf %lf" % \
                (E1, E2, nozzle.wall.layer[i].material.getPoissonRatio(), nozzle.wall.layer[i].material.getShearModulus(),
                 mu1, mu2, nozzle.wall.layer[i].material.getDensity(), alpha1, alpha2, alpha12, Lbd, nozzle.environment.hInf);
    # material properties for the thermal layer
    print >> f1, "0 0 0 0 0 %lf 0" % (nozzle.wall.layer[0].material.getThermalConductivity());
    f1.close();
    
    f2 = open("BOUNDARY.txt", 'w');
    print >> f2, "%d" % (Size[0]);
    for i in range(0,Size[0]):
        print >> f2, "%lf %lf %lf %lf" % (SolExtract[i][0], SolExtract[i][iPres], SolExtract[i][iTemp], nozzle.environment.T);
    f2.close();
    
    print "ENTER GENERATE"
    _nozzle_module.generate();       # generate the meshes for thermal and structural analyses
    print "EXIT GENERATE"
    
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
    
    
    
    
    
