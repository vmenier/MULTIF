# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy, math
from optparse import OptionParser
import textwrap
import multif
import ctypes
import numpy as np
from .. import _nozzle_module

from postprocessing import *

def runAEROS ( nozzle, output='verbose' ):
    
    # --- Get the CFD solution at the inner wall
    # SolExtract : Solution (x, y, sol1, sol2, etc.)
    # Size : [NbrVer, SolSiz]
    # idHeader : id of each solution field, e.g. mach_ver89 = SolExtract[89][idHeader['MACH']]

    # Start of example for accessing nozzle properties
    
    # First print all information related to geometry
    #print '--- Wall:'
    #print 'B-spline coefs (x-r coordinates): {:{prec}}'.format(nozzle.wall.coefs,prec=16)
    #print 'B-spline knots: {:{prec}}'.format(nozzle.wall.knots,prec=16)
    #print ''
    #for i in range(len(nozzle.wall.layer)):
    #    print '--- %s:' % nozzle.wall.layer[i].name
    #    print 'material: %s' % nozzle.wall.layer[i].material.name
    #    print 'thickness node x-coordinate: {} m'.format(nozzle.wall.layer[i].thickness.nodes[0,:])
    #    print 'thickness node local n-coordinate: {} m'.format(nozzle.wall.layer[i].thickness.nodes[1,:])
    #    print '\n'
    #print '--- Baffles:'
    #print 'number: %i' % nozzle.baffles.n
    #print 'material: %s' % nozzle.baffles.material.name
    #print 'layer 1: %s (ratio: %f)' % (nozzle.baffles.material.layer[0].material.name,nozzle.baffles.material.layer[0].ratio)
    #print 'layer 3: %s (ratio: %f)' % (nozzle.baffles.material.layer[1].material.name,nozzle.baffles.material.layer[1].ratio)
    #print 'layer 4: %s (ratio: %f)' % (nozzle.baffles.material.layer[2].material.name,nozzle.baffles.material.layer[2].ratio)
    #print 'location (x-coordinate): {} m'.format(nozzle.baffles.location)
    #print 'thickness: {} m'.format(nozzle.baffles.thickness)
    #print 'height: {} m'.format(nozzle.baffles.height)
    #print '\nFor example to access material properties of layer in baffle:'
    #print 'elastic modulus layer 1: {} Pa'.format(nozzle.baffles.material.layer[0].material.getElasticModulus())
    #print ''
    #print '--- Stringers:'
    #print 'number: %i' % nozzle.stringers.n
    #print 'material: %s' % nozzle.stringers.material.name
    #print 'thickness node x-coordinate: {} m'.format(nozzle.stringers.thickness.nodes[0,:])
    #print 'thickness node local n-coordinate: {} m'.format(nozzle.stringers.thickness.nodes[1,:])
    #print 'height node x-coordinate: {} m'.format(nozzle.stringers.height.nodes[0,:])
    #print 'height node local n-coordinate: {} m'.format(nozzle.stringers.height.nodes[1,:])    
    #print '\n' 
    
    # aside: to easily calculate wall shape or thickness as a function of the
    # x-coordinate in an item's local coordinates, use the radius method
    #r_innerwall = nozzle.wall.geometry.radius(x) # returns shape of inner wall at x in r-coordinate
    #t_thermal_layer = nozzle.wall.layer[0].thickness.radius(x) # returns thickness of thermal layer at x in local n-coordinate
    
    # Next, print all information related to materials
    #for k in nozzle.materials:
    #    print '--- %s material:' % nozzle.materials[k].name
    #    print 'type: %s' % nozzle.materials[k].type
    #    print 'density: %10.4f kg/m^3' % nozzle.materials[k].getDensity()
    #    if nozzle.materials[k].name == 'TI-HC' or nozzle.materials[k].name == 'GR-BMI':
    #        print 'elastic modulus: {} Pa'.format(nozzle.materials[k].getElasticModulus())
    #        print 'shear modulus: %e Pa' % nozzle.materials[k].getShearModulus()
    #        print 'poisson ratio: %1.4f' % nozzle.materials[k].getPoissonRatio()
    #    if nozzle.materials[k].name == 'GR-BMI':            
    #        print 'mutual influence coefs: {}'.format(nozzle.materials[k].getMutualInfluenceCoefs())
    #    print 'thermal conductivity: {} W/m-K'.format(nozzle.materials[k].getThermalConductivity())
    #    if nozzle.materials[k].name == 'TI-HC' or nozzle.materials[k].name == 'GR-BMI':
    #        print 'thermal expansion coef: {} 1/K'.format(nozzle.materials[k].getThermalExpansionCoef())
    #    print ''
        
    # So for example to obtain the Young's Moduli for the innermost Gr-BMI layer
    #[E1, E2] = nozzle.wall.layer[1].material.getElasticModulus()
    # Example, obtain all 3 coefs of thermal expansion for outermost Gr-BMI layer
    #[alpha1, alpha2, alpha12] = nozzle.wall.layer[3].material.getThermalExpansionCoef()
    # Example, obtain thermal conductivities for innermost Gr-BMI layer (there
    # are 3 values: (1st in-plane direction, 2nd in-plane direction, out-of-plane direction))
    #[k1, k2, k3] = nozzle.wall.layer[1].material.getThermalConductivity()
    # End of example

    if nozzle.method == 'NONIDEALNOZZLE':
        SolExtract = nozzle.wallResults;
        Size = [x for x in nozzle.wallResults.shape];
        idHeader = {'Temperature': 1, 'Pressure': 2};
    else:
        SolExtract, Size, idHeader  = ExtractSolutionAtWall(nozzle);

    # merge lists
    vertices = sorted(list(set(list(nozzle.wall.layer[0].thickness.nodes[0,:])+list(nozzle.wall.layer[2].thickness.nodes[0,:])+
                               list(nozzle.wall.layer[3].thickness.nodes[0,:])+list(nozzle.wall.layer[4].thickness.nodes[0,:])+
                               list(nozzle.baffles.location)+list(nozzle.stringers.thickness.nodes[0,:])+
                               list(nozzle.stringers.height.nodes[0,:]))))

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
    # lc: characteristic length (i.e. element size)
    # Tn1: number of nodes through thickness of thermal insulating layer
    # Tn2: number of nodes through thickness of gap between thermal and load layers
    # Ln: number of nodes through each half of the thickness of the load layer (thermal model)
    if nozzle.thermostructuralFidelityLevel <= 0.5:
        lc = np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[0.04,0.02]);
        Tn1 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[2,4]));
        Tn2 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[1,2]));
        Ln = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[1,2]));
    else: # nozzle.thermostructuralFidelityLevel betwen 0.5 and 1
        lc = np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[0.02,0.01]);
        Tn1 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[4,8]));
        Tn2 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[2,4]));
        Ln = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[2,4]));
    Ns   = max(2,nozzle.stringers.n); # number of panels (i.e. circumferential subdivisions of the mesh)
    Mn   = max(2,(2*math.pi*nozzle.wall.geometry.radius(points[0])/Ns)/lc+1); # Number of nodes in circumferential direction per panel
    Sn   = max(nozzle.stringers.height.radius(0)/lc+1,2) if nozzle.stringers.n > 0 else 0; # number of nodes on radial edge of stringers
    
    ## --- How to get x, y, P, T :
    #for i in range(0,Size[0]):
    #    print "VER %d : (x,y) = (%.*lf, %.*lf) , Pres = %lf, Temp = %lf" % (i, 16, SolExtract[i][0], 16, SolExtract[i][1], SolExtract[i][iPres], SolExtract[i][iTemp]);

    boundaryFlag = 1 if len(nozzle.baffles.location) > 0 else 0; # 0: inlet fixed, 1: baffles fixed, 2: both inlet and baffles fixed
    if nozzle.thermalFlag == 1:
        thermalFlag = 1;  # 0: structural analysis only, 1: both thermal and structural analyses
    else: # only perform structural analysis
        thermalFlag = 0;
    linearFlag = nozzle.linearStructuralAnalysis; # 0: nonlinear structural analysis, 1: linear structural analysis

    materialNames = [nozzle.materials[k].name for k in nozzle.materials]
    # material ids of the thermal and load layers
    M = list();
    for i in range(len(nozzle.wall.layer)):
        if i != 1: # skip the air gap layer
            M.append(materialNames.index(nozzle.wall.layer[i].material.name));
    M = [materialNames.index(nozzle.wall.layer[i].material.name) for i in range(len(nozzle.wall.layer))]
    # material id of baffles
    Mb = materialNames.index(nozzle.baffles.material.name) if len(nozzle.baffles.location) > 0 else -1
    # material id of stringers
    Ms = materialNames.index(nozzle.stringers.material.name) if nozzle.stringers.n > 0 else -1
    
    f1 = open("NOZZLE.txt", 'w');
    print >> f1, "%d %d %d %f %d %d %d %d %d" % (len(points), len(vertices), len(nozzle.materials), lc, boundaryFlag, thermalFlag, 3, 2, linearFlag);
    # points
    for i in range(len(points)):
        print >> f1, "%.*lf %.*lf %.*lf" % (16, points[i], 16, nozzle.wall.geometry.radius(points[i]), 16, nozzle.wall.geometry.radiusGradient(points[i]));
        
    # vertices
    for i in range(len(vertices)):  
        Wb = nozzle.baffles.height[nozzle.baffles.location.index(vertices[i])] if vertices[i] in nozzle.baffles.location else 0 # height of baffle
        Ws = nozzle.stringers.height.radius(vertices[i]) if nozzle.stringers.n > 0 else 0; # height of stringer
        Nb = max((Wb-Ws)/lc+1,2); # number of nodes on radial edge of baffle (not including overlap with stringer)
        Tg = nozzle.wall.layer[1].thickness.radius(0.) # thickness of gap between thermal and load layers
        Tb = nozzle.baffles.thickness[nozzle.baffles.location.index(vertices[i])] if vertices[i] in nozzle.baffles.location else 0 # thickness of baffle
        Ts = nozzle.stringers.thickness.radius(vertices[i]) if nozzle.stringers.n > 0 else 0; # thickness of stringers
        print >> f1, "%d %0.16f %d %d %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f" % (points.index(vertices[i]), Wb, Mb, Nb,
                 nozzle.wall.layer[2].thickness.radius(vertices[i]), nozzle.wall.layer[3].thickness.radius(vertices[i]),
                 nozzle.wall.layer[4].thickness.radius(vertices[i]), nozzle.wall.layer[0].thickness.radius(vertices[i]),
                 Tg, Tb, Ts, Ws);
                 
    # panels
    for i in range(1,len(vertices)):
        Nn = max(2,(vertices[i]-vertices[i-1])/lc+1); # number of nodes on longitudial edge
        print >> f1, "%d %d %d %d %d %d %d %d %d %d %d %d %d" % (M[2], M[3], M[4], Ns, Ms, Nn, Mn, Sn, M[0], M[1], Tn1, Tn2, Ln);
    # material properties
    for k in nozzle.materials:
      if nozzle.materials[k].type == 'ISOTROPIC':
        if nozzle.materials[k].name == 'CMC':
            print >> f1, "ISOTROPIC %lf %lf %lf %0.12f %lf 0" % (nozzle.materials[k].getElasticModulus(),
                     nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalExpansionCoef(), nozzle.materials[k].getThermalConductivity())
        elif nozzle.materials[k].name == 'AIR':
            print >> f1, "ISOTROPIC 0 0 %lf 0 %lf 0" % (nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalConductivity())   
        else:
            print >> f1, "ISOTROPIC %lf %lf %lf %0.12f %lf %lf" % (nozzle.materials[k].getElasticModulus(),
                     nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getDensity(),
                     nozzle.materials[k].getThermalExpansionCoef(), nozzle.materials[k].getThermalConductivity(),
                     nozzle.environment.hInf);
      else:
        [E1, E2] = nozzle.materials[k].getElasticModulus()
        [mu1, mu2] = nozzle.materials[k].getMutualInfluenceCoefs()
        [alpha1, alpha2, alpha12] = nozzle.materials[k].getThermalExpansionCoef()
        [k11, k12, k13] = nozzle.materials[k].getThermalConductivity();
        print >> f1, "ANISOTROPIC %lf %lf %lf %lf %lf %lf %lf %0.12f %0.12f %0.12f %lf %lf" % \
                (E1, E2, nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getShearModulus(), mu1, mu2,
                 nozzle.materials[k].getDensity(), alpha1, alpha2, alpha12, (k11+k12+k13)/3, nozzle.environment.hInf);
    f1.close();
    
    f2 = open("BOUNDARY.txt", 'w');
    print >> f2, "%d" % (Size[0]);
    for i in range(0,Size[0]):
        print >> f2, "%0.8f %0.8f %0.8f %0.8f" % (SolExtract[i][0], SolExtract[i][iPres], SolExtract[i][iTemp], nozzle.environment.T);
    f2.close();
    
    _nozzle_module.generate();       # generate the meshes for thermal and structural analyses

    if thermalFlag > 0:
      os.system("aeros nozzle.aeroh"); # execute the thermal analysis
      _nozzle_module.convert();        # convert temperature output from thermal analysis to input for structural analysis
    os.system("aeros nozzle.aeros"); # execute the structural analysis
    os.system("aeros nozzle.aeros.cmc"); # execute the structural analysis of the cmc layer
    
    nozzle.wallTemp = SolExtract[:,iTemp];
    
    AEROSPostProcessing ( nozzle, output );

# Kreselmeier-Steinhauser function
def ksFunction(x,p):
    return (1./p)*np.log(np.sum(np.exp(p*x)));
    
# Modified P-norm function
def pnFunction(x,p):
    return ((1./len(x))*np.sum(np.power(x,p)))**(1./p)

# Assign stresses
def assignStress(nozzle, index, output='verbose'):
    
    # Mapping from MULTI-F internal ordering to AERO-S stress file suffixes
    suffix = [0, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

    # KS and PN params for stresses
    ks_param1 = 50.;
    pn_param1 = 10.;

    # Read stress data
    filename = 'STRESS.' + str(suffix[index]);
    try:
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
    except:
        if output=='verbose':
            sys.stdout.write('WARNING: could not open STRESS file for component %i\n' % index);
        return 0;
        
    stress = data[:,-1];
    
    # Assign stresses
    stemp = np.mean(stress);
    nozzle.max_total_stress[index] = np.max(stress);
    nozzle.ks_total_stress[index] = ksFunction(stress/stemp,ks_param1)*stemp;
    nozzle.pn_total_stress[index] = pnFunction(stress/stemp,pn_param1)*stemp;

    return 0;

def assignFailureCriteria(nozzle, index, material, output='verbose'):
    
    # Mapping from MULTI-F internal ordering to AERO-S stress file suffixes
    suffix = [0, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
    
    # KS and PN params for failure criteria
    ks_param2 = 50.;
    pn_param2 = 10.;
    
    # Assign failure criteria
    if not hasattr(material,'failureType'):
        if output=='verbose':
            sys.stdout.write('WARNING: no failure type for component %i\n' % index);
        return 0;
    elif material.failureType == 'VON_MISES':
        # Von Mises is read in through the STRESS files
        filename = 'STRESS.' + str(suffix[index]);
        try: 
            data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        except IOError:
            if output=='verbose':
                sys.stdout.write('WARNING: could not open STRESS file for component %i\n' % index);
            return 0;
        stress = data[:,-1];        
        # Assign failure criterion
        yieldStress = material.yieldStress;        
        nozzle.max_failure_criteria[index] = np.max(stress/yieldStress);
        nozzle.ks_failure_criteria[index] = ksFunction(stress/yieldStress,ks_param2);
        nozzle.pn_failure_criteria[index] = pnFunction(stress/yieldStress,pn_param2);
    elif material.failureType == 'PRINCIPLE_FAILURE_STRAIN':
        sys.stderr.write('\n ## ERROR: NOT IMPLEMENTED: principle failure strain failure type\n\n');
        sys.exit(0);
    elif material.failureType == 'LOCAL_FAILURE_STRAIN':
        # Local strains are read in through the STRAINXX and STRAINYY files
        filename = 'STRAINXX.' + str(suffix[index]);
        try:
            data = np.loadtxt(filename,dtype=float,skiprows=3); # strains in 4th column (0-indexed)
        except IOError:
            if output=='verbose':
                sys.stdout.write('WARNING: could not open STRAINXX file for component %i\n' % index);
            return 0;        
        strainxx = data[:,-1];
        try:
            filename = 'STRAINYY.' + str(suffix[index]);
        except IOError:
            if output=='verbose':
                sys.stdout.write('WARNING: could not open STRAINYY file for component %i\n' % index);
            return 0;            
        data = np.loadtxt(filename,dtype=float,skiprows=3); # strains in 4th column (0-indexed)
        strainyy = data[:,-1];
        # Assign failure criterion
        failureStrain = material.getFailureLimit();
        failxx = np.empty((strainxx.size));
        failyy = np.empty((strainyy.size));
        for i in range(len(failxx)):
            if strainxx[i] >= 0:
                failxx = failureStrain[0];
            else:
                failxx = failureStrain[1];
            if strainyy[i] >= 0:
                failyy = failureStrain[2];
            else:
                failyy = failureStrain[3];
        strain = np.vstack((strainxx,strainyy));
        failStrain = np.vstack((failxx,failyy));
        nozzle.max_failure_criteria[index] = np.max(strain/failStrain);
        nozzle.ks_failure_criteria[index] = ksFunction(strain/failStrain,ks_param2);
        nozzle.pn_failure_criteria[index] = pnFunction(strain/failStrain,pn_param2);        
    else:
        sys.stderr.write('\n ## ERROR: failure type not accepted.\n\n');
        sys.exit(0);    
    
    return 0;
    
def AEROSPostProcessing ( nozzle, output='verbose' ):
    
    # ---- KS and modified P-norm parameters
    ks_param = 50.;
    pn_param = 10.;
    
    # ---- Indexing arrays
    mat = [nozzle.wall.layer[q].material for q in range(5)];
    mat.append(nozzle.stringers.material);
    for i in range(nozzle.baffles.n):
        mat.append(nozzle.baffles.material);
        
    # ---- Load and assign total stress results if necessary 
    if ( sum(nozzle.GetOutput['MAX_TOTAL_STRESS']) > 0 or 
         sum(nozzle.GetOutput['KS_TOTAL_STRESS']) > 0 or
         sum(nozzle.GetOutput['PN_TOTAL_STRESS']) > 0):
             
        for i in range(len(nozzle.GetOutput['MAX_TOTAL_STRESS'])):
            
            if ( nozzle.GetOutput['MAX_TOTAL_STRESS'][i] >= 0 or
                 nozzle.GetOutput['KS_TOTAL_STRESS'][i] >= 0 or
                 nozzle.GetOutput['PN_TOTAL_STRESS'][i] >= 0):
                     
                     assignStress(nozzle, i, output);
                     
    # ---- Load and assign failure criteria results if necessary        
    if ( sum( nozzle.GetOutput['MAX_FAILURE_CRITERIA']) > 0 or
         sum( nozzle.GetOutput['KS_FAILURE_CRITERIA']) > 0 or
         sum( nozzle.GetOutput['PN_FAILURE_CRITERIA']) > 0):
             
        for i in range(len(nozzle.GetOutput['MAX_FAILURE_CRITERIA'])):
            
            if ( nozzle.GetOutput['MAX_FAILURE_CRITERIA'][i] >= 0 or
                 nozzle.GetOutput['KS_FAILURE_CRITERIA'][i] >= 0 or
                 nozzle.GetOutput['PN_FAILURE_CRITERIA'][i] >= 0):
                     
                     assignFailureCriteria(nozzle, i, mat[i], output);
    
    # ---- Load and assign mechanical stress results if necessary
    if sum(nozzle.GetOutput['MAX_MECHANICAL_STRESS']) > 0:
        
        # Thermal layer
        filename = 'MECHANICAL_STRESS.0';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_mechanical_stress[0] = np.max(data[:,-1]);
        #nozzle.ks_mechanical_stress[0] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_mechanical_stress[0] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
        
        
        # Inner load layer
        filename = 'MECHANICAL_STRESS.1';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_mechanical_stress[2] = np.max(data[:,-1]);
        #nozzle.ks_mechanical_stress[2] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_mechanical_stress[2] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
        
        # Middle load layer
        filename = 'MECHANICAL_STRESS.2';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_mechanical_stress[3] = np.max(data[:,-1]);
        #nozzle.ks_mechanical_stress[3] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_mechanical_stress[3] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;    
        
        # Upper load layer
        filename = 'MECHANICAL_STRESS.3';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_mechanical_stress[4] = np.max(data[:,-1]);
        #nozzle.ks_mechanical_stress[4] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_mechanical_stress[4] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;  
        
        # Stringers
        filename = 'MECHANICAL_STRESS.4';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_mechanical_stress[5] = np.max(data[:,-1]);
        #nozzle.ks_mechanical_stress[5] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_mechanical_stress[5] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
    
        # Each baffle
        for i in range(7,nozzle.baffles.n+7):
            filename = 'MECHANICAL_STRESS.' + str(i-2);
            data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
            stemp = np.mean(data[:,-1]);
            nozzle.max_mechanical_stress[i-1] = np.max(data[:,-1]);
            #nozzle.ks_mechanical_stress[i-1] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
            #nozzle.pn_mechanical_stress[i-1] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
    
    # ---- Load and assign thermal stress results if necessary
    if sum(nozzle.GetOutput['MAX_THERMAL_STRESS']) > 0:
        
        # Thermal layer
        filename = 'THERMAL_STRESS.0';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_thermal_stress[0] = np.max(data[:,-1]);
        #nozzle.ks_mechanical_stress[0] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_mechanical_stress[0] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;

        # Inner load layer
        filename = 'THERMAL_STRESS.1';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_thermal_stress[2] = np.max(data[:,-1]);
        #nozzle.ks_thermal_stress[2] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_thermal_stress[2] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
        
        # Middle load layer
        filename = 'THERMAL_STRESS.2';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_thermal_stress[3] = np.max(data[:,-1]);
        #nozzle.ks_thermal_stress[3] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_thermal_stress[3] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;    
        
        # Upper load layer
        filename = 'THERMAL_STRESS.3';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_thermal_stress[4] = np.max(data[:,-1]);
        #nozzle.ks_thermal_stress[4] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_thermal_stress[4] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;  
        
        # Stringers
        filename = 'THERMAL_STRESS.4';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_thermal_stress[5] = np.max(data[:,-1]);
        #nozzle.ks_thermal_stress[5] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        #nozzle.pn_thermal_stress[5] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
    
        # Each baffle
        for i in range(7,nozzle.baffles.n+7):
            filename = 'THERMAL_STRESS.' + str(i-2);
            data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
            stemp = np.mean(data[:,-1]);
            nozzle.max_thermal_stress[i-1] = np.max(data[:,-1]);
            #nozzle.ks_thermal_stress[i-1] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
            #nozzle.pn_thermal_stress[i-1] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;

    # ---- Load and assign temperature results if necessary, AND
    # load and assign temperature ratio results if necessary
    if ( sum(nozzle.GetOutput['MAX_TEMPERATURE']) > 0 or 
         sum(nozzle.GetOutput['KS_TEMPERATURE']) > 0 or
         sum(nozzle.GetOutput['PN_TEMPERATURE']) > 0 or
         sum(nozzle.GetOutput['MAX_TEMP_RATIO']) > 0 or
         sum(nozzle.GetOutput['KS_TEMP_RATIO']) > 0 or
         sum(nozzle.GetOutput['PN_TEMP_RATIO']) > 0):    

        # Thermal layer
        data = nozzle.wallTemp;
        nozzle.max_temperature[0] = np.max(data);
        nozzle.max_temp_ratio[0] = nozzle.max_temperature[0]/nozzle.wall.layer[0].material.Tmax;
        stemp = np.mean(data);
        nozzle.ks_temperature[0] = ksFunction(data/stemp,ks_param)*stemp;
        nozzle.pn_temperature[0] = pnFunction(data/stemp,pn_param)*stemp;
        nozzle.ks_temp_ratio[0] = ksFunction(data/nozzle.wall.layer[0].material.Tmax,ks_param);
        nozzle.pn_temp_ratio[0] = pnFunction(data/nozzle.wall.layer[0].material.Tmax,pn_param);
    
        # Inner load layer
        filename = 'TEMP.1';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_temperature[2] = np.max(data[:,-1]);
        nozzle.max_temp_ratio[2] = nozzle.max_temperature[2]/nozzle.wall.layer[2].material.Tmax;        
        nozzle.ks_temperature[2] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        nozzle.pn_temperature[2] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
        nozzle.ks_temp_ratio[2] = ksFunction(data[:,-1]/nozzle.wall.layer[2].material.Tmax,ks_param);
        nozzle.pn_temp_ratio[2] = pnFunction(data[:,-1]/nozzle.wall.layer[2].material.Tmax,pn_param);
        
        # Inner load layer
        filename = 'TEMP.2';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_temperature[3] = np.max(data[:,-1]);
        nozzle.max_temp_ratio[3] = nozzle.max_temperature[3]/nozzle.wall.layer[3].material.Tmax;
        nozzle.ks_temperature[3] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        nozzle.pn_temperature[3] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;
        nozzle.ks_temp_ratio[3] = ksFunction(data[:,-1]/nozzle.wall.layer[3].material.Tmax,ks_param);
        nozzle.pn_temp_ratio[3] = pnFunction(data[:,-1]/nozzle.wall.layer[3].material.Tmax,pn_param);
        
        # Inner load layer
        filename = 'TEMP.3';
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        stemp = np.mean(data[:,-1]);
        nozzle.max_temperature[4] = np.max(data[:,-1]);
        nozzle.max_temp_ratio[4] = nozzle.max_temperature[4]/nozzle.wall.layer[4].material.Tmax;
        nozzle.ks_temperature[4] = ksFunction(data[:,-1]/stemp,ks_param)*stemp;
        nozzle.pn_temperature[4] = pnFunction(data[:,-1]/stemp,pn_param)*stemp;    
        nozzle.ks_temp_ratio[4] = ksFunction(data[:,-1]/nozzle.wall.layer[4].material.Tmax,ks_param);
        nozzle.pn_temp_ratio[4] = pnFunction(data[:,-1]/nozzle.wall.layer[4].material.Tmax,pn_param);
        
