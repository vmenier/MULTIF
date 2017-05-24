# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy, math
from optparse import OptionParser
import textwrap
import multif
import ctypes
import numpy as np
from .. import _nozzle_module

from AEROSpostprocessing import *

def runAEROS ( nozzle, output='verbose' ):

    # --- Extract flow solution at the wall
    
    # SolExtract : Solution (x, y, sol1, sol2, etc.)
    # Size : [NbrVer, SolSiz]
    # idHeader : id of each solution field, 
    #            e.g. mach_ver89 = SolExtract[89][idHeader['MACH']]
    if nozzle.method == 'NONIDEALNOZZLE':
        SolExtract = nozzle.wallResults;
        Size = [x for x in nozzle.wallResults.shape];
        idHeader = {'Temperature': 1, 'Pressure': 2};
    else:
        SolExtract, Size, idHeader  = ExtractSolutionAtWall(nozzle);
        
    ## How to get x, y, P, T :
    #for i in range(0,Size[0]):
    #    print "VER %d : (x,y) = (%.*lf, %.*lf) , Pres = %lf, Temp = %lf" % \
    #                    (i, 16, SolExtract[i][0], 16, SolExtract[i][1], \
    #                    SolExtract[i][iPres], SolExtract[i][iTemp]);

    iPres = idHeader['Pressure'];
    iTemp = idHeader['Temperature'];
    
    # --- Assemble important points to align mesh with
    
    # merge lists
    ends = [min(list(nozzle.wall.layer[0].thickness.nodes[0,:])), \
            max(list(nozzle.wall.layer[0].thickness.nodes[0,:]))]
    #for i in range(len(ends)):
    #    print ' i = %d, ends[i] = %f' % (i, ends[i])

    vertices = sorted(list(set(ends+list(nozzle.baffles.location))))
    #for j in range(len(vertices)):
    #    print ' j = %d, vertices[j] = %f' % (j, vertices[j])

    points = sorted(list(set(vertices+[item[0] for item in SolExtract]+ \
                    list(nozzle.wall.layer[0].thickness.nodes[0,:])+ \
                    list(nozzle.wall.layer[2].thickness.nodes[0,:])+ \
                    list(nozzle.wall.layer[3].thickness.nodes[0,:])+ \
                    list(nozzle.wall.layer[4].thickness.nodes[0,:])+ \
                    list(nozzle.stringers.thickness.nodes[0,:])+ \
                    list(nozzle.stringers.height.nodes[0,:]))))
    #for k in range(len(points)):
    #    print ' k = %d, points[k] = %f' % (k, points[k])
    
    # --- Set important flags
    
    # Determine how stringer height is defined:
    #     0: stringer height is defined w.r.t. nozzle wall 
    #        (specifically, the center of the load layers)
    #     1: stringer height is defined w.r.t. the axis of 
    #        rotation of the nozzle, i.e. total y coordinate
    stringerFlag = nozzle.stringers.absoluteRadialCoord;   
    
    # Determine where nozzle boundary is fixed:
    #     0: inlet fixed
    #     1: baffle edges fixed
    #     2: both inlet and baffle edges fixed
    boundaryFlag = 1 if len(nozzle.baffles.location) > 0 else 0;
    
    # Determine whether to perform structural and/or thermal analyses:
    #     0: structural analysis only
    #     1: both thermal and structural analyses
    thermalFlag = 1 if nozzle.thermalFlag == 1 else 0;

    # Determine structural analysis type:
    #     0: nonlinear structural analysis
    #     1: linear structural analysis
    linearFlag = nozzle.linearStructuralAnalysis;
    
    # --- Mesh parameters
    
    # Assign a continuous range of mesh parameters based on user-specified thermostructuralFidelityLevel float. This
    # float is between 0 and 1. Choosing 0.5 gives the baseline mesh parameters.
    
    # lc: characteristic length (i.e. element size)
    # Tn1: number of nodes through thickness of thermal insulating layer
    # Tn2: number of nodes through thickness of gap between thermal and load layers
    # Ln: number of nodes through each half of the thickness of the load layer (thermal model)
    # Mn: number of nodes in circumferential direction per baffle
    if nozzle.thermostructuralFidelityLevel <= 0.5:
        lc = np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[0.04,0.02]);
        Tn1 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[2,4]));
        Tn2 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[1,2]));
        Ln = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[1,2]));
        Mn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[25,40]));
        Sn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[5,8]));
    else: # nozzle.thermostructuralFidelityLevel betwen 0.5 and 1
        lc = np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[0.02,0.01]);
        Tn1 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[4,8]));
        Tn2 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[2,4]));
        Ln = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[2,4]));
        Mn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[40,55]));
        Sn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[8,11]));    
    # Ns: number of panels (i.e. circumferential subdivisions of the mesh)        
    Ns   = max(2,nozzle.stringers.n); 
    
    #Mn   = max(2,(2*math.pi*nozzle.wall.geometry.radius(points[0])/Ns)/lc+1); # Number of nodes in circumferential direction per panel
    #if stringerFlag == 0: # height defined w.r.t. nozzle exterior wall
    #    Sn = max(nozzle.stringers.height.radius(0)/lc+1,2) if nozzle.stringers.n > 0 else 0; # number of nodes on radial edge of stringers
    #else:
    #    Sn = max((nozzle.stringers.height.radius(0)-nozzle.wall.geometry.radius(points[0]))/lc+1,2) if nozzle.stringers.n > 0 else 0;
    
    # --- Material properties
    
    materialNames = [nozzle.materials[k].name for k in nozzle.materials]
    # material ids of the thermal and load layers
    M = list();
    for i in range(len(nozzle.wall.layer)):
        if i != 1: # skip the air gap layer
            M.append(materialNames.index(nozzle.wall.layer[i].material.name));
    M = [materialNames.index(nozzle.wall.layer[i].material.name) \
         for i in range(len(nozzle.wall.layer))]
    # material id of baffles
    Mb = materialNames.index(nozzle.baffles.material.name) if \
         len(nozzle.baffles.location) > 0 else -1
    # material id of stringers
    Ms = materialNames.index(nozzle.stringers.material.name) if \
         nozzle.stringers.n > 0 else -1

    # --- Write nozzle definition to file for Aero-S
    
    f1 = open("NOZZLE.txt", 'w');
    verboseFlag = 0.;
    if output == 'verbose':
        verboseFlag = 1.;
    print >> f1, "%d %d %d %f %d %d %d %d %d %d %d" % (len(points), \
        len(vertices), len(nozzle.materials), lc, boundaryFlag, \
        thermalFlag, 3, 2, linearFlag, verboseFlag, stringerFlag);
    
    # points
    for i in range(len(points)):
        # Tg: thickness of gap between thermal and load layers
        Tg = nozzle.wall.layer[1].thickness.radius(0.) 
        # Ts: thickness of stringers
        Ts = nozzle.stringers.thickness.radius(points[i]) if nozzle.stringers.n > 0 else 0; 
        # Ws: # height of stringer
        Ws = nozzle.stringers.height.radius(points[i]) if nozzle.stringers.n > 0 else 0; 
        print >> f1, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e" % (points[i],nozzle.wall.geometry.radius(points[i]),
             nozzle.wall.geometry.radiusGradient(points[i]),
             nozzle.wall.layer[2].thickness.radius(points[i]), 
             nozzle.wall.layer[3].thickness.radius(points[i]),
             nozzle.wall.layer[4].thickness.radius(points[i]), 
             nozzle.wall.layer[0].thickness.radius(points[i]), Tg, Ts, Ws);
        
    # vertices
    for i in range(len(vertices)): 
        # Wb: height of baffle 
        Wb = nozzle.baffles.height[nozzle.baffles.location.index(vertices[i])] \
             if vertices[i] in nozzle.baffles.location else 0
        # Nb: number of nodes on radial edge of baffle (not including overlap with stringer)
        #Nb = max(Wb/lc-Ns+1,2); 
        Nb = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,1],[2,4]));
        Tb = nozzle.baffles.thickness[nozzle.baffles.location.index(vertices[i])] \
             if vertices[i] in nozzle.baffles.location else 0 # thickness of baffle
        print >> f1, "%d %0.16e %d %d %0.16e" % (points.index(vertices[i]), Wb, Mb, Nb, Tb);
                 
    # panels
    for i in range(1,len(vertices)):
        #Nn = max(2,(vertices[i]-vertices[i-1])/lc+1); # number of nodes on longitudial edge
        Nn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,1],[7,15]));
        print >> f1, "%d %d %d %d %d %d %d %d %d %d %d %d %d" % (M[2], M[3], M[4], Ns, Ms, Nn, Mn, Sn, M[0], M[1], Tn1, Tn2, Ln);
    # material properties
    for k in nozzle.materials:
      if nozzle.materials[k].type == 'ISOTROPIC':
        if nozzle.materials[k].name == 'CMC':
            print >> f1, "ISOTROPIC %0.16e %0.16e %0.16e %0.16e %0.16e 0" % (nozzle.materials[k].getElasticModulus(),
                     nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalExpansionCoef(), nozzle.materials[k].getThermalConductivity())
        elif nozzle.materials[k].name == 'AIR':
            print >> f1, "ISOTROPIC 0 0 %0.16e 0 %0.16e 0" % (nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalConductivity())   
        else:
            print >> f1, "ISOTROPIC %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e" % (nozzle.materials[k].getElasticModulus(),
                     nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalExpansionCoef(), nozzle.materials[k].getThermalConductivity(),
                     nozzle.environment.hInf);
      else:
        [E1, E2] = nozzle.materials[k].getElasticModulus()
        [mu1, mu2] = nozzle.materials[k].getMutualInfluenceCoefs()
        [alpha1, alpha2, alpha12] = nozzle.materials[k].getThermalExpansionCoef()
        [k11, k12, k13] = nozzle.materials[k].getThermalConductivity();
        print >> f1, "ANISOTROPIC %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e" % \
                (E1, E2, nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getShearModulus(), mu1, mu2,
                 nozzle.materials[k].getDensity(), alpha1, alpha2, alpha12, (k11+k12+k13)/3, nozzle.environment.hInf);
    f1.close();
    
    f2 = open("BOUNDARY.txt", 'w');
    print >> f2, "%d" % (Size[0]);
    for i in range(0,Size[0]):
        print >> f2, "%0.16e %0.16e %0.16e %0.16e" % (SolExtract[i][0], SolExtract[i][iPres], SolExtract[i][iTemp], nozzle.environment.T);
    f2.close();
    
    _nozzle_module.generate();       # generate the meshes for thermal and structural analyses

    if thermalFlag > 0:
      os.system("aeros nozzle.aeroh"); # execute the thermal analysis
      _nozzle_module.convert();        # convert temperature output from thermal analysis to input for structural analysis
    os.system("aeros nozzle.aeros"); # execute the structural analysis
    os.system("aeros nozzle.aeros.cmc"); # execute the structural analysis of the cmc layer
    
    #nozzle.wallTemp = SolExtract[:,iTemp];
    
    return 0;
