# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy, math
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
    print 'material: %s' % nozzle.baffles.material.name
    print 'location (x-coordinate): {} m'.format(nozzle.baffles.location)
    print 'thickness: {} m'.format(nozzle.baffles.thickness)
    print 'height: {} m'.format(nozzle.baffles.height)
    print ''
    print '--- Stringers:'
    print 'number: %i' % nozzle.stringers.n
    print 'material: %s' % nozzle.stringers.material.name
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

    # merge lists
    vertices = sorted(list(set(list(nozzle.wall.layer[0].thickness.nodes[0,:])+list(nozzle.wall.layer[1].thickness.nodes[0,:])+
                               list(nozzle.wall.layer[2].thickness.nodes[0,:])+list(nozzle.wall.layer[3].thickness.nodes[0,:])+
                               list(nozzle.baffles.location))))
    #for k in range(len(vertices)):
    #    print ' k = %d, vertices[k] = %f' % (k, vertices[k])

    points = sorted(list(set(vertices+[item[0] for item in SolExtract])))
    #for k in range(len(points)):
    #    print ' k = %d, points[k] = %f' % (k, points[k])
    
    # Average thermal conductivity of load layer (W/m*K)
    #[k11, k12, k13] = nozzle.wall.layer[1].material.getThermalConductivity();
    #k2 = nozzle.wall.layer[2].material.getThermalConductivity();
    #[k31, k32, k33] = nozzle.wall.layer[3].material.getThermalConductivity();
    #Lbd = (k11+k12+k13+3*k2+k31+k32+k33)/9
    
    iPres = idHeader['Pressure'];
    iTemp = idHeader['Temperature'];
    
    # --- Mesh parameters
    lc   = 0.02; # Characteristic length (i.e. element size)
    Ns   = max(2,nozzle.stringers.n); # number of panels (i.e. circumferential subdivisions of the mesh)
    Mn   = max(2,(2*math.pi*nozzle.wall.geometry.radius(points[0])/Ns)/lc+1); # Number of nodes in circumferential direction per panel
    Tn   = 4;  # Number of nodes through thickness of thermal insulating layer

    boundaryFlag = 1 if len(nozzle.baffles.location) > 0 else 0; # 0: inlet fixed, 1: baffles fixed, 2: both inlet and baffles fixed
    thermalFlag = 1;  # 0: structural analysis only, 1: both thermal and structural analyses

    # XXX stringer height and thickness (currently assumed to be constant)
    Ws = nozzle.stringers.height[0] if nozzle.stringers.n > 0 else 0;
    Ts = nozzle.stringers.thickness[0] if nozzle.stringers.n > 0 else 0;

    materialNames = [nozzle.materials[k].name for k in nozzle.materials]
    # material ids of the thermal and load layers
    M = [materialNames.index(nozzle.wall.layer[i].material.name) for i in range(len(nozzle.wall.layer))]
    # material id of baffles
    Mb = materialNames.index(nozzle.baffles.material.name) if len(nozzle.baffles.location) > 0 else -1
    # material id of stringers
    Ms = materialNames.index(nozzle.stringers.material.name) if nozzle.stringers.n > 0 else -1
    
    f1 = open("NOZZLE.txt", 'w');
    print >> f1, "%d %d %d %f %d %d %d" % (len(points), len(vertices), len(nozzle.materials), lc, boundaryFlag, thermalFlag, 3);
    # points
    for i in range(len(points)):
        print >> f1, "%lf %lf" % (points[i], nozzle.wall.geometry.radius(points[i]));
    # vertices
    for i in range(len(vertices)):  
        Wb = nozzle.baffles.height[nozzle.baffles.location.index(vertices[i])] if vertices[i] in nozzle.baffles.location else 0 # height of baffle
        Nb = max((Wb-Ws)/lc+1,2); # number of nodes on radial edge of baffle (not including overlap with stringer)
        Tb = nozzle.baffles.thickness[nozzle.baffles.location.index(vertices[i])] if vertices[i] in nozzle.baffles.location else 0 # thickness of baffle
        print >> f1, "%d %f %d %d %lf %lf %lf %lf %lf %lf" % (points.index(vertices[i]), Wb, Mb, Nb,
                 nozzle.wall.layer[1].thickness.radius(vertices[i]), nozzle.wall.layer[2].thickness.radius(vertices[i]),
                 nozzle.wall.layer[3].thickness.radius(vertices[i]), nozzle.wall.layer[0].thickness.radius(vertices[i]),
                 Tb, Ts);
    # panels
    for i in range(1,len(vertices)):
        Nn = max(2,(vertices[i]-vertices[i-1])/lc+1); # number of nodes on longitudial edge
        Sn = max(Ws/lc+1,2); # number of nodes on radial edge of stringers
        print >> f1, "%d %d %d %d %f %d %d %d %d %d %d" % (M[1], M[2], M[3], Ns, Ws, Ms, Nn, Mn, Sn, M[0], Tn);
    # material properties
    for k in nozzle.materials:
      if nozzle.materials[k].type == 'ISOTROPIC':
        if nozzle.materials[k].name == 'CMC':
            print >> f1, "ISOTROPIC 0 0 %lf 0 %lf 0" % (nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalConductivity())
        else:
            print >> f1, "ISOTROPIC %lf %lf %lf %lf %lf %lf" % (nozzle.materials[k].getElasticModulus(),
                     nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalExpansionCoef(), nozzle.materials[k].getThermalConductivity(),
                     nozzle.environment.hInf);
      else:
        [E1, E2] = nozzle.materials[k].getElasticModulus()
        [mu1, mu2] = nozzle.materials[k].getMutualInfluenceCoefs()
        [alpha1, alpha2, alpha12] = nozzle.materials[k].getThermalExpansionCoef()
        [k11, k12, k13] = nozzle.materials[k].getThermalConductivity();
        print >> f1, "ANISOTROPIC %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf" % \
                (E1, E2, nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getShearModulus(), mu1, mu2,
                 nozzle.materials[k].getDensity(), alpha1, alpha2, alpha12, (k11+k12+k13)/3, nozzle.environment.hInf);
    f1.close();
    
    f2 = open("BOUNDARY.txt", 'w');
    print >> f2, "%d" % (Size[0]);
    for i in range(0,Size[0]):
        print >> f2, "%lf %lf %lf %lf" % (SolExtract[i][0], SolExtract[i][iPres], SolExtract[i][iTemp], nozzle.environment.T);
    f2.close();
    
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
    
    
    
    
    
