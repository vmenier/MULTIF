# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy, math
from optparse import OptionParser
import textwrap
import multif
import ctypes
import numpy as np
from .. import _nozzle_module

from AEROSpostprocessing import *
from SU2postprocessing import ExtractSolutionAtWall

def runAEROS ( nozzle, output='verbose' ):
    
    if nozzle.dim == '3D':
            
        # Introduction to this function and available nozzle data
        print
        print '*****'
        print
                
        print 'Exterior environment information:'
        print nozzle.environment.__dict__
        print
            
        print 'NOZZLE INTERIOR WALL:'
        print 'The interior nozzle wall is parameterized with 3 B-splines.'
        print 'One B-spline parameterizes the centerline, and two more '
        print 'parameterize the major and minor axes of an ellipse centered '
        print 'on the centerline, drawn in the vertical plane. Two other '
        print 'parameters, the shovel exit height and shovel inlet start '
        print 'angle are used to define the flattening out of the geometry '
        print 'near the exit. The nozzle is axisymmetric up until the throat.'
        print 'Nozzle wall centerline geometry:'
        print nozzle.wall.centerline.geometry.__dict__
        print 'Nozzle wall major axis geometry:'
        print nozzle.wall.majoraxis.geometry.__dict__
        print 'Nozzle wall minor axis geometry:'
        print nozzle.wall.minoraxis.geometry.__dict__
        print 'Nozzle length: %f' % nozzle.length
        print 'Nozzle wall struct:'
        print nozzle.wall.__dict__
        print 'Nozzle wall shovel exit height: %f' % nozzle.wall.shovel_height
        print 'Nozzle wall shovel inlet start angle: %f' % nozzle.wall.shovel_start_angle
        print
        
        print 'WALL LAYERS:'
        print 'Each wall layer is defined using a bilinear distribution, '
        print 'where thickness is a function of global axial coordinate x '
        print 'and global angular coordinate theta. Data is arranged in a '
        print 'Numpy array in a rectilinear grid as follows:'
        print 'nozzle.wall.layer[i].thicknessNodes = '
        print '  [[x1, y1, t1],'
        print '   [x2, y1, t2],'
        print '   [x3, y1, t3],'
        print '   [x4, y1, t4],'
        print '   [x1, y2, t5],'
        print '       ... '
        print '        etc.  ]'
        for i in range(len(nozzle.wall.layer)):
            print 'Nozzle wall layer: %s' % nozzle.wall.layer[i].name
            print nozzle.wall.layer[i].__dict__
        print
        print 'Material properties are accessed the same way as before'
        print
        
        print 'STRINGERS:'
        print 'Since the nozzle is nonaxisymmetric, stringers are defined by '
        print 'angular position around the centerline. Thicknesses are a '
        print 'function of global axial X-coordinate and stringer angle. '
        print 'A stringer is assumed to travel along the same angular coord.'
        print nozzle.stringers.__dict__
        print 'Example access of stringers height definition for 2nd stringer:'
        print nozzle.stringers.height[1].__dict__
        print 'Example access of stringers thickness definition for 2nd stringer:'
        print nozzle.stringers.thickness[1].__dict__
        print
        
        print 'BAFFLES:'
        print 'Baffles are defined as before.'
        print nozzle.baffles.__dict__
        print
        
        print 'EXTERIOR:'
        'The nozzle exterior is complex and fully defined in CAD, not here.'
        print nozzle.exterior.__dict__
        print
        
        print 'Entire NOZZLE struct:'
        for k in nozzle.__dict__:
            print k
            print nozzle.__dict__[k]
            print        
        
    else: # assume 2D axisymmetric nozzle

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
        ends = [min(list(nozzle.wall.layer[0].thicknessNodes[:,0])), \
                max(list(nozzle.wall.layer[0].thicknessNodes[:,0]))]
        #for i in range(len(ends)):
        #    print ' i = %d, ends[i] = %f' % (i, ends[i])
    
        vertices = sorted(list(set(ends+list(nozzle.baffles.location))))
        #for j in range(len(vertices)):
        #    print ' j = %d, vertices[j] = %f' % (j, vertices[j])
    
        points = sorted(list(set(vertices+[item[0] for item in SolExtract]+ \
                        list(nozzle.wall.layer[0].thicknessNodes[:,0])+ \
                        list(nozzle.wall.layer[2].thicknessNodes[:,0])+ \
                        list(nozzle.wall.layer[3].thicknessNodes[:,0])+ \
                        list(nozzle.wall.layer[4].thicknessNodes[:,0])+ \
                        list(nozzle.stringers.thicknessNodes[:,0])+ \
                        list(nozzle.stringers.heightNodes[:,0]))))
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
        linearFlag = nozzle.linearStructuralAnalysisFlag;
        
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
        if nozzle.wallTempFlag == 1: # wall temperature is assigned by user
            for i in range(0,Size[0]):
                Temp = nozzle.wall.temperature.geometry.radius(SolExtract[i][0])
                print >> f2, "%0.16e %0.16e %0.16e %0.16e" % (SolExtract[i][0], SolExtract[i][iPres], Temp, nozzle.environment.T);
            f2.close();        
        else: # wall temperature is extracted from flow
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
