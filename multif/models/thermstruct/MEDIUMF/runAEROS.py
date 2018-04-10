import os
import numpy as np
import subprocess

import multif
from ....models import _nozzle_module

from ....models.aero.MEDIUMF.postprocessing import ExtractSolutionAtWall
from ....models.aero.HIGHF.aeros import hf_FluidStructureInterpolation


def getMass(nozzle, run_analysis=True, output='verbose'):
    """
    Build mesh if necessary, run AERO-S mass calculation, and return total mass
    of nozzle, and the wall mass of the nozzle (excluding stringers & baffles).
    """

    # Check that runAEROS has previously been called to generate the mesh and
    # AeroS input files. If not, generate the mesh.
    if not os.path.exists('nozzle.aeros.mass'):
        prepareAEROS(nozzle, run_analysis=run_analysis, output=output)

    if run_analysis:
        # Calculate mass of thermal layer
        if 'SLURM_NTASKS' in os.environ:
            subprocess.call(['srun','-n','1','aeros','nozzle.aeros.cmc.mass'])
        else:
            subprocess.call(['aeros','nozzle.aeros.cmc.mass'])

        # Calculate mass of thermal model (thermal layers + approximate load layers)
        # Inaccurate since load layers are averaged.
        #os.system("aeros nozzle.aeroh.mass") 
        #m2 = float(np.loadtxt("MASS.txt.thermal")) # mass of thermal model

        # Calculate mass of load layers and stringers and baffles
        if 'SLURM_NTASKS' in os.environ:
            subprocess.call(['srun','-n','1','aeros','nozzle.aeros.mass'])
        else:
            subprocess.call(['aeros','nozzle.aeros.mass'])

    # Post-process

    # Get mass of thermal layer
    m1 = float(np.loadtxt("MASS.txt.cmc")) # mass of CMC structural model

    # Get mass of load layers and stringers and baffles
    m3a = float(np.loadtxt("MASS.txt")) # mass of structural model
    m3b = float(np.loadtxt("MASS.txt.2")) # mass of load layer in structural model

    total_mass = m1 + m3a
    wall_mass = m1 + m3b

    return total_mass, wall_mass


def writeGeometry(nozzle, mesh_params=None, output='verbose'):
    """
    Write geometry files used for thermal/structural analyses and mass 
    calculations
    """

    # --- Set important flags

    # Determine how stringer height is defined:
    #     0: stringer height is defined w.r.t. nozzle wall 
    #        (specifically, the center of the load layers)
    #     1: stringer height is defined w.r.t. the axis of 
    #        rotation of the nozzle, i.e. total z coordinate
    stringerFlag = nozzle.stringers.absoluteRadialCoord   
    
    # Determine where nozzle boundary is fixed:
    #     0: inlet fixed
    #     1: baffle edges fixed
    #     2: both inlet and baffle edges fixed
    boundaryFlag = 1 if len(nozzle.baffles.location) > 0 else 0
    
    # Determine whether to perform structural and/or thermal analyses:
    #     0: structural analysis only
    #     1: both thermal and structural analyses
    thermalFlag = 1 if nozzle.thermalFlag == 1 else 0

    # Determine whether thermal model needs to be built so mass can be
    # calculated. A thermal analysis is not necessarily run unless needed.
    if 'MASS' in nozzle.qoi.names or 'MASS_WALL_ONLY' in nozzle.qoi.names:
        thermalFlagForMass = 1
    elif thermalFlag == 1:
        thermalFlagForMass = 1
    else:
        thermalFlagForMass = 0

    # Determine structural analysis type:
    #     0: nonlinear structural analysis
    #     1: linear structural analysis
    linearFlag = nozzle.linearStructuralAnalysisFlag

    # --- Mesh parameters
    
    # Assign a continuous range of mesh parameters based on user-specified thermostructuralFidelityLevel float. This
    # float is between 0 and 1. Choosing 0.5 gives the baseline mesh parameters.
    
    # lc: characteristic length (i.e. element size)
    # Tn1: number of nodes through thickness of thermal insulating layer
    # Tn2: number of nodes through thickness of gap between thermal and load layers
    # Ln: number of nodes through each half of the thickness of the load layer (thermal model)
    # Mn: number of nodes in circumferential direction per baffle (note: for transfinite meshing of baffles, set Mn such that (Mn-1)%3 == 0)
    # Ns: number of panels (i.e. circumferential subdivisions of mesh)
    # Nn: number of nodes on longitudinal edge
    # Nb: number of nodes on radial edge of baffle (not including overlap with
    #     stringer). If set to -1, unstructured meshing is used for baffle
    # Sn: number of nodes on radial edges of stringers. Should be the same as
    #     Nb. If set to -1, unstructured meshing is used for stringers.
    if mesh_params is not None:

        lc = mesh_params['lc']
        Tn1 = mesh_params['Tn1']
        Tn2 = mesh_params['Tn2']
        Ln = mesh_params['Ln']
        Mn = mesh_params['Mn']
        #Sn = mesh_params['Sn']
        Ns = mesh_params['Ns']
        Nn = mesh_params['Nn']
        Nb = mesh_params['Nb']

    else:
        if nozzle.thermostructuralFidelityLevel <= 0.5:
            lc = np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[0.32,0.16])
            Tn1 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[2,4]))
            Tn2 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[1,2]))
            Ln = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[1,2]))
            Mn = 3*np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[10,30]))+1
            #Sn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,0.5],[5,8]))
        else: # nozzle.thermostructuralFidelityLevel betwen 0.5 and 1
            lc = np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[0.16,0.08])
            Tn1 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[4,8]))
            Tn2 = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[2,4]))
            Ln = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[2,4]))
            Mn = 3*np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[30,50]))+1
            #Sn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0.5,1],[8,11]))    
        # Ns: number of panels (i.e. circumferential subdivisions of the mesh)        
        Ns   = max(2,nozzle.stringers.n)     
        Nn = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,1],[10,40]))
        Nb = np.round(np.interp(nozzle.thermostructuralFidelityLevel,[0,1],[20,50]))
    Sn = Nb

    # --- Assemble important points to align mesh with
    
    ends = [nozzle.wall.geometry.xstart, nozzle.wall.geometry.xend]
    vertices = sorted(list(set(ends+list(nozzle.baffles.location))))
    points = sorted(list(set(vertices+ \
                    list(nozzle.wall.layer[0].thicknessNodes[:,0])+ \
                    list(nozzle.wall.layer[2].thicknessNodes[:,0])+ \
                    list(nozzle.wall.layer[3].thicknessNodes[:,0])+ \
                    list(nozzle.wall.layer[4].thicknessNodes[:,0])+ \
                    list(nozzle.stringers.thicknessNodes[:,0])+ \
                    list(nozzle.stringers.heightNodes[:,0]))))
    
    # --- Material properties
    
    materialNames = [nozzle.materials[k].name for k in nozzle.materials]
    # material ids of the thermal and load layers
    M = list()
    for i in range(len(nozzle.wall.layer)):
        if i != 1: # skip the air gap layer
            M.append(materialNames.index(nozzle.wall.layer[i].material.name))
    M = [materialNames.index(nozzle.wall.layer[i].material.name) \
         for i in range(len(nozzle.wall.layer))]
    # material id of baffles
    Mb = materialNames.index(nozzle.baffles.material.name) if \
         len(nozzle.baffles.location) > 0 else -1
    # material id of stringers
    Ms = materialNames.index(nozzle.stringers.material.name) if \
         nozzle.stringers.n > 0 else -1

    # --- Set inner wall shape

    # Start and end angle of shovel geometry
    theta_in = nozzle.wall.shovel_start_angle
    theta_out = nozzle.wall.shovel_end_angle

    # Write nozzle inner wall shape to file using 3d dim irregardless of 
    # specified parameterization and dimension
    if nozzle.dim == '3D':

        nc  = len(nozzle.wall.centerline.geometry.coefs[0])
        nm1 = len(nozzle.wall.majoraxis.geometry.coefs[0])
        nm2 = len(nozzle.wall.minoraxis.geometry.coefs[0])
        kc  = len(nozzle.wall.centerline.geometry.knots)
        km1 = len(nozzle.wall.majoraxis.geometry.knots)
        km2 = len(nozzle.wall.minoraxis.geometry.knots)

        xcentercoefs = nozzle.wall.centerline.geometry.coefs[0][:]
        ycentercoefs = nozzle.wall.centerline.geometry.coefs[1][:]
        xmajorcoefs  = nozzle.wall.majoraxis.geometry.coefs[0][:]
        ymajorcoefs  = nozzle.wall.majoraxis.geometry.coefs[1][:]
        xminorcoefs  = nozzle.wall.minoraxis.geometry.coefs[0][:]
        yminorcoefs  = nozzle.wall.minoraxis.geometry.coefs[1][:]
        centerknots  = nozzle.wall.centerline.geometry.knots[:]
        majorknots   = nozzle.wall.majoraxis.geometry.knots[:]
        minorknots   = nozzle.wall.minoraxis.geometry.knots[:]

    else: # assume 2D dimension given, upscale to 3D

        # Assign centerline here (constant)
        xi = nozzle.wall.geometry.xstart # inlet x-coord
        xe = nozzle.wall.geometry.xend   # outlet x-coord
        nc = 4
        kc = nc + 4
        xcentercoefs = [xi, xi, xe, xe]
        ycentercoefs = [0., 0., 0., 0.]
        centerknots  = [0,0,0,0] + range(1,nc-3) + [kc,kc,kc,kc]

        # Assign major axis here
        if nozzle.param == '2D':   

            nm1 = len(nozzle.wall.geometry.coefs[0])
            km1 = len(nozzle.wall.knots)
            print nozzle.wall.geometry.coefs
            xmajorcoefs  = nozzle.wall.geometry.coefs[0][:]
            ymajorcoefs  = nozzle.wall.geometry.coefs[1][:]
            majorknots   = nozzle.wall.geometry.knots[:]

        else: # assume 3D parameterization downscaled to 2d dimension

            # In this situation, the 2d axisymmetric geometry is defined by a 
            # piecewise linear function. Let each break node in the piecewise
            # linear function be a B-spline coefficient.

            nm1 = nozzle.wall.geometry.nx
            km1 = nozzle.wall.geometry.nx + 4 # 3rd degree B-spline
            xmajorcoefs  = nozzle.wall.geometry.nodes[:,0]
            ymajorcoefs  = nozzle.wall.geometry.nodes[:,1]
            majorknots   = [0,0,0,0] + range(1,nm1-3) + [km1,km1,km1,km1]

        # Assign minor axis here (axisymmetric, so same as major axis)
        nm2 = nm1
        km2 = km1
        xminorcoefs = xmajorcoefs
        yminorcoefs = ymajorcoefs
        minorknots  = majorknots

    # Write BSPLINE.txt file
    f0 = open("BSPLINE.txt", 'w')
    if nozzle.dim == '3D':
        print >> f0, "%d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d" % (
                nc, kc, nm1, km1, nm2, km2, 
                theta_in, theta_out, nozzle.baffles.halfWidth,
                nozzle.exterior.geometry['top'].angle, nozzle.exterior.geometry['top'].offset,
                nozzle.exterior.geometry['top'].a, nozzle.exterior.geometry['top'].b,
                nozzle.exterior.geometry['bottom'].angle, nozzle.exterior.geometry['bottom'].offset,
                nozzle.exterior.geometry['bottom'].a, nozzle.exterior.geometry['bottom'].b,
                len(points),
                nozzle.wall.layer[2].nAngularBreaks, nozzle.wall.layer[3].nAngularBreaks,
                nozzle.wall.layer[4].nAngularBreaks, nozzle.wall.layer[0].nAngularBreaks)
    else:
        print >> f0, "%d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d" % (
                nc, kc, nm1, km1, nm2, km2, 
                theta_in, theta_out, nozzle.baffles.halfWidth,
                nozzle.exterior.geometry['top'].angle, nozzle.exterior.geometry['top'].offset,
                nozzle.exterior.geometry['top'].a, nozzle.exterior.geometry['top'].b,
                nozzle.exterior.geometry['bottom'].angle, nozzle.exterior.geometry['bottom'].offset,
                nozzle.exterior.geometry['bottom'].a, nozzle.exterior.geometry['bottom'].b,
                0, 0, 0, 0, 0)        

    # center line control points
    for i in range(nc):
        print >> f0, "%0.16e %0.16e" % (xcentercoefs[i], ycentercoefs[i])

    # major axis knot values
    for i in range(kc):
        print >> f0, "%0.16e" % (centerknots[i])

    # major axis control points
    for i in range(nm1):
        print >> f0, "%0.16e %0.16e" % (xmajorcoefs[i], ymajorcoefs[i])

    # major axis knot values
    for i in range(km1):
        print >> f0, "%0.16e" % (majorknots[i])

    # minor axis control points
    for i in range(nm2):
        print >> f0, "%0.16e %0.16e" % (xminorcoefs[i], yminorcoefs[i])

    # minor axis knot values
    for i in range(km2):
        print >> f0, "%0.16e" % (minorknots[i])

    if nozzle.dim == '3D':

        # thickness of inner load layer
        for i in range(len(points)):
            for j in range(nozzle.wall.layer[2].nAngularBreaks):
                angle = nozzle.wall.layer[2].thicknessNodes[j*nozzle.wall.layer[2].nAxialBreaks,1]
                print >> f0, "%0.16e %0.16e %0.16e" % \
                        (points[i], angle, nozzle.wall.layer[2].thickness.height(points[i],angle))

        # thickness of middle load layer
        for i in range(len(points)):
            for j in range(nozzle.wall.layer[3].nAngularBreaks):
                angle = nozzle.wall.layer[3].thicknessNodes[j*nozzle.wall.layer[3].nAxialBreaks,1]
                print >> f0, "%0.16e %0.16e %0.16e" % \
                        (points[i], angle, nozzle.wall.layer[3].thickness.height(points[i],angle))

        # thickness of outer load layer
        for i in range(len(points)):
            for j in range(nozzle.wall.layer[4].nAngularBreaks):
                angle = nozzle.wall.layer[4].thicknessNodes[j*nozzle.wall.layer[4].nAxialBreaks,1]
                print >> f0, "%0.16e %0.16e %0.16e" % \
                        (points[i], angle, nozzle.wall.layer[4].thickness.height(points[i],angle))

        # thickness of thermal layer
        for i in range(len(points)):
            for j in range(nozzle.wall.layer[0].nAngularBreaks):
                angle = nozzle.wall.layer[0].thicknessNodes[j*nozzle.wall.layer[0].nAxialBreaks,1]
                print >> f0, "%0.16e %0.16e %0.16e" % \
                        (points[i], angle, nozzle.wall.layer[0].thickness.height(points[i],angle))

    f0.close()
    
    # --- Write nozzle definition to file for Aero-S

    if nozzle.dim == '3D':

        f1 = open("NOZZLE.txt", 'w')
        verboseFlag = 0.
        if output == 'verbose':
            verboseFlag = 1.
        print >> f1, "%d %d %d %f %d %d %d %d %d %d %d" % (len(points), \
            len(vertices), len(nozzle.materials), lc, boundaryFlag, \
            thermalFlagForMass, 3, 2, linearFlag, verboseFlag, stringerFlag)
        
        # points
        for i in range(len(points)):
            # Tg: thickness of gap between thermal and load layers
            Tg = nozzle.wall.layer[1].thickness.radius(0.) 
            # Ts: thickness of stringers
            Ts1 = nozzle.stringers.thickness[0].radius(points[i]) if nozzle.stringers.n > 0 else 0 
            Ts2 = nozzle.stringers.thickness[1].radius(points[i]) if nozzle.stringers.n > 0 else 0
            # Ws: # height of stringer
            print >> f1, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e" % \
                (points[i],
                -1, #nozzle.wall.geometry.radius(points[i]),
                -1, #nozzle.wall.geometry.radiusGradient(points[i]),
                nozzle.wall.layer[2].thickness.height(points[i],0.), # thickness at theta = 0
                nozzle.wall.layer[3].thickness.height(points[i],0.),
                nozzle.wall.layer[4].thickness.height(points[i],0.),
                nozzle.wall.layer[0].thickness.height(points[i],0.),
                Tg, Ts1, Ts2) #Ts, Ws)

        # vertices
        for i in range(len(vertices)): 
            # Wb: height of baffle 
            Wb = 1 if vertices[i] in nozzle.baffles.location else 0
            Tb = nozzle.baffles.thickness[nozzle.baffles.location.index(vertices[i])] \
                if vertices[i] in nozzle.baffles.location else 0 # thickness of baffle
            print >> f1, "%d %0.16e %d %d %0.16e" % (points.index(vertices[i]), Wb, Mb, Nb, Tb)
        
        # panels
        for i in range(1,len(vertices)):
            #Nn = max(2,(vertices[i]-vertices[i-1])/lc+1) # number of nodes on longitudial edge
#            print >> f1, "%d %d %d %d %d %d %d %d %d %d %d %d %d" % (M[2], M[3], M[4], Ns, Ms, Nn, Mn, Sn, M[0], M[1], Tn1, Tn2, Ln)
            print >> f1, "%d %d %d %d %d %d %d %d %d %d %d %d %d" % (M[2], M[3], M[4], Ns, Ms, Nn, Mn, Sn, M[0], M[1], Tn1, Tn2, Ln)

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
                            nozzle.environment.hInf)
            else:
                [E1, E2] = nozzle.materials[k].getElasticModulus()
                [mu1, mu2] = nozzle.materials[k].getMutualInfluenceCoefs()
                [alpha1, alpha2, alpha12] = nozzle.materials[k].getThermalExpansionCoef()
                [k11, k12, k13] = nozzle.materials[k].getThermalConductivity()
                print >> f1, "ANISOTROPIC %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e" % \
                        (E1, E2, nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getShearModulus(), mu1, mu2,
                        nozzle.materials[k].getDensity(), alpha1, alpha2, alpha12, (k11+k12+k13)/3, nozzle.environment.hInf)
        
        f1.close()
    
    else: # assume 2d dimension of nozzle
    
        f1 = open("NOZZLE.txt", 'w')
        verboseFlag = 0.
        if output == 'verbose':
            verboseFlag = 1.
        print >> f1, "%d %d %d %f %d %d %d %d %d %d %d" % (len(points), \
            len(vertices), len(nozzle.materials), lc, boundaryFlag, \
            thermalFlagForMass, 3, 2, linearFlag, verboseFlag, stringerFlag)
        
        # points
        for i in range(len(points)):
            # Tg: thickness of gap between thermal and load layers
            Tg = nozzle.wall.layer[1].thickness.radius(0.) 
            if nozzle.stringers.heightDefinition == 'EXTERIOR' or \
               nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT':
                if isinstance(nozzle.stringers.thickness,list):
                    Ts = nozzle.stringers.thickness[0].radius(points[i]) if nozzle.stringers.n > 0 else 0
                    Ws = nozzle.stringers.thickness[1].radius(points[i]) if nozzle.stringers.n > 0 else 0
                else:
                    Ts = nozzle.stringers.thickness.radius(points[i]) if nozzle.stringers.n > 0 else 0 
                    Ws = Ts
            else:
                # Ts: thickness of stringers
                Ts = nozzle.stringers.thickness.radius(points[i]) if nozzle.stringers.n > 0 else 0 
                # Ws: # height of stringer
                Ws = nozzle.stringers.height.radius(points[i]) if nozzle.stringers.n > 0 else 0 
            print >> f1, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e" % (points[i],nozzle.wall.geometry.radius(points[i]),
                nozzle.wall.geometry.radiusGradient(points[i]),
                nozzle.wall.layer[2].thickness.radius(points[i]), 
                nozzle.wall.layer[3].thickness.radius(points[i]),
                nozzle.wall.layer[4].thickness.radius(points[i]), 
                nozzle.wall.layer[0].thickness.radius(points[i]),
                Tg, Ts, Ws)

        # vertices
        for i in range(len(vertices)): 
            # Wb: height of baffle 
            # Wb = nozzle.baffles.height[nozzle.baffles.location.index(vertices[i])] \
            #     if vertices[i] in nozzle.baffles.location else 0
            Wb = 1 if vertices[i] in nozzle.baffles.location else 0
            Tb = nozzle.baffles.thickness[nozzle.baffles.location.index(vertices[i])] \
                if vertices[i] in nozzle.baffles.location else 0 # thickness of baffle
            print >> f1, "%d %0.16e %d %d %0.16e" % (points.index(vertices[i]), Wb, Mb, Nb, Tb)
                    
        # panels
        for i in range(1,len(vertices)):
            #Nn = max(2,(vertices[i]-vertices[i-1])/lc+1) # number of nodes on longitudial edge
            print >> f1, "%d %d %d %d %d %d %d %d %d %d %d %d %d" % (M[2], M[3], M[4], Ns, Ms, Nn, Mn, Sn, M[0], M[1], Tn1, Tn2, Ln)

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
                            nozzle.environment.hInf)
            else:
                [E1, E2] = nozzle.materials[k].getElasticModulus()
                [mu1, mu2] = nozzle.materials[k].getMutualInfluenceCoefs()
                [alpha1, alpha2, alpha12] = nozzle.materials[k].getThermalExpansionCoef()
                [k11, k12, k13] = nozzle.materials[k].getThermalConductivity()
                print >> f1, "ANISOTROPIC %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e" % \
                        (E1, E2, nozzle.materials[k].getPoissonRatio(), nozzle.materials[k].getShearModulus(), mu1, mu2,
                        nozzle.materials[k].getDensity(), alpha1, alpha2, alpha12, (k11+k12+k13)/3, nozzle.environment.hInf)
        
        f1.close()

    # END NOZZLE.txt writing

    return


def writeBoundaryConditions2D(nozzle, run_analysis=True, output='verbose'):
    """
    Write BOUNDARY.txt boundary condition file for Aero-S input.
    """

    # Determine whether to perform structural and/or thermal analyses:
    #     0: structural analysis only
    #     1: both thermal and structural analyses
    thermalFlag = 1 if nozzle.thermalFlag == 1 else 0

    # --- Extract flow solution at the wall
    
    try:
        # SolExtract : Solution (x, y, sol1, sol2, etc.)
        # Size : [NbrVer, SolSiz]
        # idHeader : id of each solution field, 
        #            e.g. mach_ver89 = SolExtract[89][idHeader['MACH']]
        if nozzle.dim != '3D':
            if nozzle.method == 'NONIDEALNOZZLE':
                SolExtract = nozzle.wallResults
                Size = [x for x in nozzle.wallResults.shape]
                idHeader = {'Temperature': 1, 'Pressure': 2}
            else:
                SolExtract, Size, idHeader  = ExtractSolutionAtWall(nozzle)

            iPres = idHeader['Pressure']
            iTemp = idHeader['Temperature']

        ## How to get x, y, P, T :
        #for i in range(0,Size[0]):
        #    print "VER %d : (x,y) = (%.*lf, %.*lf) , Pres = %lf, Temp = %lf" % \
        #                    (i, 16, SolExtract[i][0], 16, SolExtract[i][1], \
        #                    SolExtract[i][iPres], SolExtract[i][iTemp])
        
        f2 = open("BOUNDARY.txt", 'w')
        if nozzle.dim != '3D':
            print >> f2, "%d" % (Size[0])
            if nozzle.wallTempFlag == 1: # wall temperature is assigned by user
                for i in range(0,Size[0]):
                    Temp = nozzle.wall.temperature.geometry.radius(SolExtract[i][0])
                    print >> f2, "%0.16e %0.16e %0.16e %0.16e" % (SolExtract[i][0], SolExtract[i][iPres], Temp, nozzle.environment.T)
            else: # wall temperature is extracted from flow
                for i in range(0,Size[0]):
                    print >> f2, "%0.16e %0.16e %0.16e %0.16e" % (SolExtract[i][0], SolExtract[i][iPres], SolExtract[i][iTemp], nozzle.environment.T)
        else:
            print >> f2, "2"
            print >> f2, "%0.16e 0 0 %0.16e" % (nozzle.wall.geometry.xstart, nozzle.environment.T)
            print >> f2, "%0.16e 0 0 %0.16e" % (nozzle.wall.geometry.xend, nozzle.environment.T)
        f2.close()
    except:
        f2 = open("BOUNDARY.txt", 'w')
        print >> f2, "2"
        print >> f2, "%0.16e 0 0 %0.16e" % (nozzle.wall.geometry.xstart, nozzle.environment.T)
        print >> f2, "%0.16e 0 0 %0.16e" % (nozzle.wall.geometry.xend, nozzle.environment.T)
        f2.close()
        print "ERROR: BOUNDARY.txt could not be written"
        print

    return


def writeBoundaryConditions3D(nozzle, run_analysis=True, output='verbose'):
    """
    Write TEMPERATURES.txt.thermal, TEMPERATURES.txt, and PRESSURES.txt if
    nozzle geometry is 3D nonaxisymmetric.
    """

    # Determine whether to perform structural and/or thermal analyses:
    #     0: structural analysis only
    #     1: both thermal and structural analyses
    thermalFlag = 1 if nozzle.thermalFlag == 1 else 0

    try:
        if nozzle.dim == '3D' and run_analysis:
            #--- Get solution from fluid calculation

            print "Interface AEROS"

            MshNam_str = "nozzle.mesh"
            MshNam_cfd = "nozzle.su2"
            SolNam_cfd = "nozzle.dat"

            Crd, Tri, Pres, Temp = multif.models.aero.HIGHF.aeros.hf_FluidStructureInterpolation(MshNam_str, MshNam_cfd, SolNam_cfd)

            if thermalFlag > 0:
                # temperatures for the thermal model
                f0 = open("TEMPERATURES.txt.thermal", 'r')
                f0.readline()
                f1 = open("TEMPERATURES.txt.thermal.3d", 'w')
                print >> f1, "TEMPERATURE"
                for line in f0:
                    nodeId = int(line.split()[0])
                    print >> f1, "%d %0.16e" % (nodeId, Temp[nodeId-1])
                f0.close()
                f1.close()
                os.rename("TEMPERATURES.txt.thermal.3d", "TEMPERATURES.txt.thermal")
            else:
                # temperatures for the structural model
                f0 = open("TEMPERATURES.txt", 'r')
                f0.readline()
                f1 = open("TEMPERATURES.txt.3d", 'w')
                print >> f1, "TEMPERATURE"
                for line in f0:
                    nodeId = int(line.split()[0])
                    print >> f1, "%d %0.16e" % (nodeId, Temp[nodeId-1])
                f0.close()
                f1.close()
                os.rename("TEMPERATURES.txt.3d", "TEMPERATURES.txt")
            
            # pressures for the structural model
            f0 = open("PRESSURES.txt", 'r')
            f0.readline()
            f1 = open("PRESSURES.txt.3d", 'w')
            print >> f1, "PRESSURE"
            i = 0
            for line in f0:
                elemId = int(line.split()[0])
                avgPres = (Pres[Tri[i][0]]+Pres[Tri[i][1]]+Pres[Tri[i][2]])/3
                print >> f1, "%d %0.16e" % (elemId, avgPres)
                i = i+1
            f0.close()
            f1.close()
            os.rename("PRESSURES.txt.3d", "PRESSURES.txt")

    except:
        print("WARNING: 3D AERO-S boundary condition files are filled with fluff.")
        # temperatures for the structural model
        f1 = open("TEMPERATURES.txt", 'w')
        print >> f1, "TEMPERATURE"
        print >> f1, "1 300"
        f1.close()           
        # pressures for the structural model
        f1 = open("PRESSURES.txt", 'w')
        print >> f1, "PRESSURE"
        print >> f1, "1 50000"
        f1.close()

    return


def prepareAEROS(nozzle, run_analysis=True, output='verbose'):
    """
    Generate mesh and AERO-S input files.
    """

    writeGeometry(nozzle, output=output)
    writeBoundaryConditions2D(nozzle, run_analysis=run_analysis, output=output)
    _nozzle_module.generate()
    # 3D boundary conditions depend on mesh generation
    writeBoundaryConditions3D(nozzle, run_analysis=run_analysis, output=output)

    return


def callThermalAnalysis():

    # Thermal analysis
    if 'SLURM_NTASKS' in os.environ:
        subprocess.call(['srun','-n','1','aeros','nozzle.aeroh'])
    else:
        subprocess.call(['aeros','nozzle.aeroh'])

    # Convert temp. output from thermal analysis to input for structural analysis
    _nozzle_module.convert()

    # Structural analysis of CMC layer
    if 'SLURM_NTASKS' in os.environ:
        subprocess.call(['srun','-n','1','aeros','nozzle.aeros.cmc'])
    else:
        subprocess.call(['aeros','nozzle.aeros.cmc'])

    return


def callStructuralAnalysis():

    # Structural analysis of load layers + baffles and stringers
    if 'SLURM_NTASKS' in os.environ:
        subprocess.call(['srun','-n','1','aeros','nozzle.aeros'])
    else:
        subprocess.call(['aeros','nozzle.aeros'])

    return


def runAEROS(nozzle, output='verbose', run_analysis=True, mesh_params=None):
    """
    All-in-one function that generates mesh, AERO-S input files, and runs
    AERO-S.
    """
    
    # Determine whether to perform structural and/or thermal analyses:
    #     0: structural analysis only
    #     1: both thermal and structural analyses
    thermalFlag = 1 if nozzle.thermalFlag == 1 else 0

    # Determine whether to perform structural analysis:
    #     0: no
    #     1: yes
    structuralFlag = 1 if nozzle.structuralFlag == 1 else 0

    # Generate meshes and AERO-S input files
    prepareAEROS(nozzle, run_analysis=run_analysis, output=output)

    # Execute analyses
    if run_analysis and thermalFlag:
        callThermalAnalysis()
    if run_analysis and structuralFlag:
        callStructuralAnalysis()

    return 0
