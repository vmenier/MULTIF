# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy, math
from optparse import OptionParser
import textwrap
import ctypes
import numpy as np

import multif

from scipy import interpolate as itp
from scipy.interpolate import splev, splrep

from .. import SU2
from .. import nozzle as noz

class optionsmesh:
    def __init__(self):
        pass



def ComputeDs(nozzle):
    yplus = nozzle.cfd.bl_yplus;
    
    gam   = 1.4;
    R     = 287.06;
    Cv    = 717.645;
    Su    = 110.4;
    
    M      = nozzle.mission.mach;
    Ps     = nozzle.environment.P;
    Ts     = nozzle.environment.T;
    D      = nozzle.wall.geometry.radius(nozzle.wall.geometry.length);
    
    mu     = 1.716e-5*((Ts/273.15)**1.5)*(273.15 + Su)/(Ts + Su);      # Sutherland law 
    rho    = Ps / ( (gam-1.) * Cv * Ts )                               # density
    c      = np.sqrt( gam * (Ps/rho));                                 # speed of sound
    U      = M*c                                                       # velocity
    Rey    = rho*U*D/mu;                                               # Reynolds number
    
    Re_x = rho * U * D / mu;
    Cf = 0.026/math.pow(Re_x,1./7.);
    
    tau_w = 0.5*Cf*rho*U*U;
    Ufric = np.sqrt(tau_w/rho);
    
    ds = yplus*mu/(Ufric*rho);
    
    print "yplus %lf Re_x %lf " % (yplus, Re_x)
    
    return ds;
    

def GenerateNozzleMesh (nozzle):
    import tempfile
    
    #hdl, nozzle.tmpGeoNam = tempfile.mkstemp(suffix='.geo');
    #hdl, nozzle.tmpMshNam = tempfile.mkstemp(suffix='.mesh');
    
    nozzle.tmpGeoNam = 'nozzle_tmp.geo';
    nozzle.tmpMshNam = 'nozzle.su2';
    
    # --- Write geo file
    
    mesh_options = optionsmesh();
    mesh_options.xwall  = nozzle.cfd.x_wall;
    mesh_options.ywall  = nozzle.cfd.y_wall;
    mesh_options.hl     = nozzle.cfd.meshhl;
    mesh_options.method = nozzle.method; # Euler or RANS
    
    # --- Compute ds based on y+ value
    
    mesh_options.ds = ComputeDs(nozzle);
    #mesh_options.ds          = nozzle.cfd.bl_ds;
    
    mesh_options.ratio       = nozzle.cfd.bl_ratio 
    mesh_options.thickness   = nozzle.cfd.bl_thickness; 
    
    mesh_options.x_thrust = nozzle.cfd.x_thrust;      
    
    NozzleGeoFile(nozzle.tmpGeoNam, mesh_options);
    #NozzleGeoFileRoundedEdges(nozzle.tmpGeoNam, mesh_options);
    
    # --- Call Gmsh
    
    try :   
        CallGmsh(nozzle);
    except:
        print "\n  ## ERROR : Mesh generation failed.\n";
        sys.exit(0);
    
    if nozzle.method == "RANS":
        import subprocess
        gmsh_executable = 'gmsh';
        try :
            cmd = [gmsh_executable, '-2', "exit.geo", '-o', "exit.mesh"];
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
        except:
            raise;
    
    ## --- Mesh preprocessing
    #try :
    #    MeshPrepro(nozzle);
    #except :
    #    print "\n  ## ERROR : Mesh preprocessing failed.\n";
    #    sys.exit(0);
    

def GenerateNozzleExitMesh(nozzle):
    
    mesh_options = optionsmesh();
    mesh_options.xwall  = nozzle.cfd.x_wall;
    mesh_options.ywall  = nozzle.cfd.y_wall;
    mesh_options.hl     = nozzle.cfd.meshhl;
    mesh_options.method = nozzle.method; # Euler or RANS
    
    # --- Compute ds based on y+ value
    
    mesh_options.ds = ComputeDs(nozzle);
    #mesh_options.ds          = nozzle.cfd.bl_ds;
    
    mesh_options.ratio       = nozzle.cfd.bl_ratio 
    mesh_options.thickness   = nozzle.cfd.bl_thickness; 
    
    mesh_options.x_thrust = nozzle.cfd.x_thrust;      
    
    ExitGeoFile("exit.geo", mesh_options);
    
    import subprocess
    gmsh_executable = 'gmsh';
    try :
        cmd = [gmsh_executable, '-2', "exit.geo", '-o', "exit.mesh"];
        print "%s" % cmd
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
    except:
        raise;
    
        


def Project_Wall (nozzle, motion):
    
    Nbv = len(motion);
    
    # --- Project baseline surface onto the new CAD model
    
    x=[];
    for i in range(0,Nbv):    
        x.append(float(motion[i][1]))
    
    ynew = [];
    dydx = [];
    
    # --- Get baseline length
    
    lenBas = max(x)-min(x);
    #print "LENGTH BASELINE : %lf" % lenBas
    
    lenCur = nozzle.wall.geometry.length;
    
    scax = lenCur/lenBas;
    
    xnew = [];    
    for i in range(0,Nbv):
        xnew.append(0.9999*scax*x[i])
    
    xwall  = nozzle.xwall;
    ywall  = nozzle.ywall;
    
    #print "MAX xwall %lf, max xnew %lf" % (max(xwall), max(xnew))
    
    #from scipy.interpolate import interp1d
    #f = interp1d(xwall, ywall)
    #xx = np.linspace(0, max(xwall), num=41, endpoint=True)
    #import matplotlib.pyplot as plt
    #plt.plot(xwall, ywall, '-', xx, f(xx), 'o', xnew, f(xnew), 'o')
    #plt.show()
    
    #sys.exit(1)
    
    _meshutils_module.py_BSplineGeo3LowF (nozzle.wall.knots, nozzle.wall.coefs, xnew, ynew, dydx);
    
    motion_new = [];
    for i in range(0,Nbv):
        #motion_new.append([motion[i][0],xnew[i],f(xnew[i])])
        #print "%lf %lf" % (xnew[i],f(xnew[i]))
        motion_new.append([motion[i][0],xnew[i],ynew[i]])
    
    return motion_new;


def Extract_Boundary_Vertices(mesh_name, pyRefs):
    
    from .. import _meshutils_module
    
    #--- Extract boundary vertices from mesh
    
    PyVid_Out = [];
    PyCrd_Out = [];
    PyRef_Tab = [];
        
    _meshutils_module.py_Extract_Vertices_Ref (mesh_name, pyRefs, PyCrd_Out, PyVid_Out, PyRef_Tab);
        
    refmax=-1;
    for iRef in range(len(PyRef_Tab)/2):
        refmax = max(refmax,PyRef_Tab[2*iRef]);
    
    Bdr = [[] for x in range((refmax+1))];
    
    idv=0;
    for iRef in range(len(PyRef_Tab)/2):
                
        Nbv = PyRef_Tab[2*iRef+1];
        
        for i in range(Nbv):
            
            vid = PyVid_Out[idv];
            
            x = PyCrd_Out[3*idv+0];
            y = PyCrd_Out[3*idv+1];
            z = PyCrd_Out[3*idv+2];
            
            Bdr[PyRef_Tab[2*iRef]].append([vid,x,y,z])
            idv=idv+1;
    
    return Bdr;
    
def Get3Dto2DEquivArea(nozzle, xnew_tab):
    
    #from .. import nozzle as noz
    #
    #geometry = noz.geometry;
    #
    #majoraxisTmp = geometry.Bspline(nozzle.wall.majoraxis.coefs);
    #minoraxisTmp = geometry.Bspline(nozzle.wall.minoraxis.coefs);
    #centerTmp = geometry.Bspline(nozzle.wall.centerline.coefs);
    #
    #majoraxisVal = majoraxisTmp.radius(xnew_tab);
    #minoraxisVal = minoraxisTmp.radius(xnew_tab);
    #
    #fr1 = majoraxisTmp.radius
    #fr2 = minoraxisTmp.radius
    #fz  = centerTmp.radius
    #
    #params = np.zeros(100)
    #params[0] = 0.099; # z crd of the cut at the throat
    #params[1] = 0.122638;  # z crd of the flat exit (bottom of the shovel)
    #params[3] = 0.0;        # x_throat 
    #params[4] = 4.0;        # x_exit 
    
    ynew_tab = multif.HIGHF.MF_GetRadius (xnew_tab, nozzle)
    
    return ynew_tab;


def GenerateNozzleMesh_Deform (nozzle):
    
    from .. import _meshutils_module
    
    sys.stdout.write(" -- Generate nozzle mesh (deform)\n")
    
    h_tip = 0.012;
    
    motion_hdl = [];
    
    pathsrc = "%s/baseline_meshes_3d/" % (os.path.dirname(os.path.abspath(__file__)));
    
    #--- Extract boundary vertices from mesh
    
    #mesh_name    = "%sbaseline_%s_%s.su2" % (pathsrc, nozzle.method, nozzle.meshsize);
    mesh_name   = "%sbaseline_%s_%s.su2" % (pathsrc, nozzle.method.lower(), nozzle.cfd.mesh_size.lower());
    
    print mesh_name  
     
    if not os.path.exists(mesh_name):
        sys.stderr.write("  ## ERROR mesh generation deform: baseline mesh not found!\n \
        Expected: %s\n" % mesh_name);
        sys.exit(1);
    
    pyRefs = [1, 2, 3, 9];
    
    Bdr = Extract_Boundary_Vertices(mesh_name, pyRefs);
        
    #fil = open("mesh_motion.dat", "w");
    
    xwall  = nozzle.cfd.x_wall;
    ywall  = nozzle.cfd.y_wall;
    
    xwall_max = max(xwall);
    
    #for iref in range(len(Bdr)):
    #    Nbv = len(Bdr[iref]);
    #    
    #    if Nbv < 1 :
    #        continue;
    #        
    #    sys.stdout.write("- Ref %d : \n" % iref);
    #    
    #    for i in range(Nbv):
    #        sys.stdout.write("\t Ver %d : %lf %lf\n" % (Bdr[iref][i][0], Bdr[iref][i][1], Bdr[iref][i][2]));
    #    
    #sys.exit(1)
    
    #--- Write outside of nozzle's deformation file
    
    ref_outside = 3;
    
    xbas_max = max(float(l[1]) for l in Bdr[ref_outside])
    xbas_min = min(float(l[1]) for l in Bdr[ref_outside])
    
    #xx = [-0.67 , -0.65 , 0.1548, xwall[-1]];
    #yy = [0.7244, 0.4244, 0.41,   ywall[-1]+h_tip]
    
    xx = [-0.67 , -0.65 , 0.148,  xwall[-1]];
    yy = [0.7244, 0.7244, 0.7,  ywall[-1]+h_tip]
    
    #print "xwall %lf ywall end %lf h_tip %lf " % (xwall[-1],ywall[-1],h_tip)
    #sys.exit(1)
    
    tck = splrep(xx, yy, xb=xx[0], xe=xx[-1], k=2)
    
    f_outside = itp.interp1d(xx, yy, kind='cubic')
    
    ######## BEGIN PLOT
    ###x3 = np.linspace(xx[0], xx[-1], 200)
    ###y3 = splev(x3, tck)
    ###
    ###y32 = f_outside(x3)
    ###
    ###import matplotlib.pyplot as plt
    ###plt.plot(x3,y3, '-')
    ###plt.plot(x3,y32, '-')
    ###plt.plot(xx,yy, 'o')
    ###plt.show()
    ###
    ###sys.exit(1);
    ###
    ######## END PLOT
    
    Nbv = len(Bdr[ref_outside]);
    
    for i in range(Nbv):
        
        vid = Bdr[ref_outside][i][0];
        x = Bdr[ref_outside][i][1];
        
        #print "x %lf xbax in %lf %lf xwall in %lf %lf " % (x, xbas_min, xbas_max, xx[0], xwall_max)
        xnew = xx[0] + (x-xbas_min)/(xbas_max-xbas_min)*(xwall_max-xx[0]);
        #print "xnew %lf x %lf xbax in %lf %lf xwall in %lf %lf " % (xnew, x, xbas_min, xbas_max, xx[0], xwall_max)
        
        ynew = f_outside(xnew);
        
        motion_hdl.append([vid-1, xnew, ynew, Bdr[ref_outside][i][1], Bdr[ref_outside][i][2]]);
        
        #fil.write("%d %le %le\n" % (vid-1, xnew, ynew));
    
    #--- Project tip of nozzle
    
    ref_tip = 2;
    
    Nbv = len(Bdr[ref_tip]);
    
    ymax = max(float(l[2]) for l in Bdr[ref_tip])
    ymin = min(float(l[2]) for l in Bdr[ref_tip])
    
    for i in range(Nbv):
        
        vid = Bdr[ref_tip][i][0];
        x   = Bdr[ref_tip][i][1];
        y   = Bdr[ref_tip][i][2];
        
        #print "x %lf xbax in %lf %lf xwall in %lf %lf " % (x, xbas_min, xbas_max, xx[0], xwall_max)
        xnew = xwall[-1];
        #print "xnew %lf x %lf xbax in %lf %lf xwall in %lf %lf " % (xnew, x, xbas_min, xbas_max, xx[0], xwall_max)
        ynew = ywall[-1] + (y-ymin)/(ymax-ymin)*h_tip;
        
        #print "ref %d : ver %d : (%lf %lf) -> (%lf %lf)" % (ref_tip, vid, x, y, xnew, ynew)
        
        motion_hdl.append([vid-1, xnew, ynew, Bdr[ref_tip][i][1], Bdr[ref_tip][i][2]]);
        #fil.write("%d %le %le\n" % (vid-1, xnew, ynew));
        
    # --- Project inner wall
    
    ref_wall = 1;
    
    Nbv = len(Bdr[ref_wall]);
    
    xnew_tab = [];
    ynew_tab = [];
    dydx = [];
    
    xmin = min(float(l[1]) for l in Bdr[ref_wall])
    xmax = max(float(l[1]) for l in Bdr[ref_wall])
    
    for i in range(Nbv):        
        x = Bdr[ref_wall][i][1];
        xnew = xwall[0] + (x-xmin)/(xmax-xmin)*(xwall[-1]-xwall[0]);
        xnew_tab.append(xnew);
    
    
    if nozzle.param == "2D":
        _meshutils_module.py_BSplineGeo3LowF (nozzle.wall.knots, nozzle.wall.coefs, xnew_tab, ynew_tab, dydx);
    else :
        xnew_tab = np.array(xnew_tab);        
        ynew_tab = Get3Dto2DEquivArea(nozzle, xnew_tab);
        
        
    for i in range(Nbv):
        
        vid = Bdr[ref_wall][i][0];
        motion_hdl.append([vid-1, xnew_tab[i], ynew_tab[i], Bdr[ref_wall][i][1], Bdr[ref_wall][i][2]]);
        #fil.write("%d %le %le\n" % (vid-1, xnew_tab[i], ynew_tab[i]));
        
    # --- Project thrust marker
    
    # if nozzle.method == "EULER":
        
    #     ref_thrust = 9;
        
    #     Nbv = len(Bdr[ref_thrust]);
        
    #     x = Bdr[ref_thrust][0][1];
    #     xnew = xx[0] + (x-xbas_min)/(xbas_max-xbas_min)*(xwall_max-xx[0]);
        
    #     ynew_tab = [];
    #     dydx = [];
        
    #     if nozzle.param == "2D":
    #         _meshutils_module.py_BSplineGeo3LowF (nozzle.wall.knots, nozzle.wall.coefs, xnew_tab, ynew_tab, dydx);
    #     else :
    #         xnew_tab = np.array([xnew]);
    #         ynew_tab = Get3Dto2DEquivArea(nozzle, xnew_tab);
            
    #     ymax = ynew_tab[0];
        
    #     ymax_bas =  max(float(l[2]) for l in Bdr[ref_thrust])
        
    #     for i in range(Nbv):
            
    #         vid = Bdr[ref_thrust][i][0];
    #         y   = Bdr[ref_thrust][i][2];
            
    #         ynew = y/ymax_bas*ymax;
            
    #         motion_hdl.append([vid-1, xnew, ynew, Bdr[ref_thrust][i][1], Bdr[ref_thrust][i][2]]);
    #         #fil.write("%d %le %le\n" % (vid-1, xnew, ynew));
        
    
    #--- Write mesh motion file
    
    fil = open("mesh_motion.dat", "w");
    
    for i in range(len(motion_hdl)):
        fil.write("%d %le %le\n" % (motion_hdl[i][0], motion_hdl[i][1], motion_hdl[i][2]));
    
    fil.close();
    
    
    ###--- Debug : plot boundary
    ##
    ##import matplotlib.pyplot as plt
    ##
    ##x0 = []
    ##x1 = []
    ##y0 = []
    ##y1 = []
    ##
    ##for i in range(0,len(motion_hdl)):
    ##    x0.append(float(motion_hdl[i][1]))
    ##    y0.append(float(motion_hdl[i][2]))
    ##    
    ##    x1.append(float(motion_hdl[i][3]))
    ##    y1.append(float(motion_hdl[i][4]))
    ##    
    ##
    ##plt.plot(x1,y1, ".", label="Initial")
    ##plt.plot(x0,y0, ".", label="Moved")
    ##plt.legend()
    ##plt.show();
    
    
    
    # --- Setup config file
    
    config = SU2.io.Config();
    
    # -- Note : the values don't matter. All markers must be defined otherwise SU2 exits.
    config.MARKER_HEATFLUX = "( PhysicalLine6, 0.0, PhysicalLine7, 0.0 )"
    config.MARKER_INLET    = "( PhysicalLine1, 955.000000, 97585.000000, 1.0, 0.0, 0.0 )"
    config.MARKER_FAR      = "( ( PhysicalLine5, PhysicalLine4 ) )"
    config.MARKER_SYM      = "( ( PhysicalLine2 ) )"
    config.MARKER_OUTLET   = "( PhysicalLine3, 18754.000000)"
    
    config.MESH_FILENAME= "baseline_coarse.su2"
    config.DV_KIND= "SURFACE_FILE"
    config.DV_MARKER= "( PhysicalLine7 )"
    config.MOTION_FILENAME= "mesh_motion.dat"
    
    config.DEFORM_LINEAR_SOLVER  = "FGMRES"
    config.DEFORM_LINEAR_ITER    = "500"
    config.DEFORM_NONLINEAR_ITER = "5"
    config.DEFORM_CONSOLE_OUTPUT = "YES"
    config.DEFORM_TOL_FACTOR     = "1e-6"
    config.DEFORM_STIFFNESS_TYPE = "WALL_DISTANCE"
    
    config.HOLD_GRID_FIXED       = "NO"
    config.HOLD_GRID_FIXED_COORD = "(-1e6,-1e6,-1e6,1e6,1e6,1e6)"
    config.VISUALIZE_DEFORMATION = "YES"
    config.MARKER_MOVING         = "( PhysicalLine7 )"
    
    config.NUMBER_PART = 1;
    
    
    config.MESH_FILENAME          = mesh_name
    config.DV_KIND                = 'SURFACE_FILE'
    config.DV_MARKER              = '( ( PhysicalLine3, PhysicalLine2, PhysicalLine1, PhysicalLine9 ) )'
    config.MOTION_FILENAME        = 'mesh_motion.dat'
    config.DEFORM_LINEAR_SOLVER   = 'FGMRES'
    config.DEFORM_LINEAR_ITER     = 500
    config.DEFORM_NONLINEAR_ITER  = 10
    config.DEFORM_CONSOLE_OUTPUT  = 'YES'
    config.DEFORM_TOL_FACTOR      = 1e-06
    config.DEFORM_STIFFNESS_TYPE  = 'WALL_DISTANCE'
    config.HOLD_GRID_FIXED        = 'NO'
    config.HOLD_GRID_FIXED_COORD  = '(-1e6,-1e6,-1e6,1e6,1e6,1e6)'
    config.VISUALIZE_DEFORMATION  = 'YES'
    config.MARKER_MOVING          = '( PhysicalLine3, PhysicalLine2, PhysicalLine1, PhysicalLine9 )'
    config.NUMBER_PART            =  1
    config.MARKER_EULER           = '( ( PhysicalLine1, PhysicalLine2, PhysicalLine3 ) )'
    config.MARKER_INLET           = '( PhysicalLine8, 955.000000, 97585.000000, 1.0, 0.0, 0.0, PhysicalLine4,  228.016984, 22181.944264, 1.0, 0.0, 0.0 )'
    config.MARKER_FAR             = '( ( PhysicalLine5 ) )'
    config.MARKER_SYM             = '( ( PhysicalLine7 ) )'
    config.MARKER_OUTLET          = '( PhysicalLine6, 18754.000000)'
    config.MARKER_THRUST          = '( PhysicalLine9 )'
    config.MESH_OUT_FILENAME      = nozzle.cfd.mesh_name;
    
    config.SU2_RUN = nozzle.cfd.su2_run;
    
    
    # --- Run SU2
    
    info = SU2.run.DEF(config)

    
    
    

def CallGmsh (nozzle):
    import subprocess
    gmsh_executable = 'gmsh';
    try :
        cmd = [gmsh_executable, '-2', nozzle.tmpGeoNam, '-o', nozzle.tmpMshNam];
        print "%s" % cmd
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
    except:
        raise;
    
def MeshPrepro (nozzle):
    from .. import _meshutils_module
    
    out = _meshutils_module.py_MeshPrepro2D (nozzle.tmpMshNam, nozzle.cfd.mesh_name);
    
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
    
    xinlet = xwall[0];
    xoutlet = xwall[nx-1];
    
    CrdBox = [[0 for x in range(2)] for y in range(9)] 

    CrdBox[1][0] = xinlet;              CrdBox[1][1] = 0;
    CrdBox[2][0] = xoutlet;             CrdBox[2][1] = 0;
    CrdBox[3][0] = 6.6;                 CrdBox[3][1] = 0;
    CrdBox[4][0] = 6.6;                 CrdBox[4][1] = 4.5;
    CrdBox[5][0] = -0.67 + xinlet;      CrdBox[5][1] = 4.5;
    CrdBox[6][0] = -0.67 + xinlet;      CrdBox[6][1] = 0.7244;
    CrdBox[7][0] = 0.1548 + xinlet;     CrdBox[7][1] = 0.7244;
    CrdBox[8][0] = xoutlet;             CrdBox[8][1] = ywall[nx-1]+0.012;

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
    
    crdThrust = [0.0,0.0]
    
    vid = 11;
    for i in range(0,nx) :
        vid=vid+1;
        fil.write('Point(%d)  = {%lf, %lf, 0, %lf};\n'% (vid, xwall[nx-i-1], ywall[nx-i-1], sizWal));
        
        if ( math.fabs(x_thrust-xwall[nx-i-1]) < 1e-6 ):
            exit_vid = vid;
            crdThrust = [xwall[nx-i-1], ywall[nx-i-1]]
    if ( exit_vid == -1 ):
        sys.stderr.write(" ## ERROR : Coordinates for the thrust computation don't match.\n");
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
    #fil.write('Line(9)  = {6, 7};\n');
    #fil.write('Line(10) = {7, 8};\n');
    fil.write('BSpline(10) = {6, 7, 8};\n');
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
    
    savIdx = i;
    
    fil.write('Line(13) = {%d, 1};\n' % (i));
    
    # --- Plane surface
    
    fil.write('Line Loop(14) = {1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12,15, 13};\n');
    fil.write('Plane Surface(14) = {14};\n');
    
    hl1 = hl[0];
    hl2 = hl[1];
    hl3 = hl[2];
    hl4 = hl[3];
    hl5 = hl[4];
    
    fields = [[-100, 4.6, -10, CrdBox[6][1]+0.3,    hl3], \
              [-100, 4.6, -10, CrdBox[6][1],    hl1],\
              [-100, 4.6, -10, CrdBox[6][1],    hl5],\
              [-100, 0.3, -10, CrdBox[6][1],   hl4],\
              [-100, 100, -10, CrdBox[6][1]+0.3,    hl2],\
              [-100, xoutlet+0.04, -10, CrdBox[6][1], hl4]];
    
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
        fil.write('Field[%d].EdgesList = {10,11,12,15}; \n' % (NbrFld+2));
        fil.write('Field[%d].NodesList = {6,%d};      \n' % ((NbrFld+2), savIdx));
        fil.write('Field[%d].hfar = 1;                 \n' % (NbrFld+2));
        fil.write('Field[%d].hwall_n = %le;       \n' % ((NbrFld+2), ds));
        fil.write('Field[%d].hwall_t = 0.300000;       \n' % (NbrFld+2));
        fil.write('Field[%d].ratio = %le;              \n' % ((NbrFld+2),ratio));
        fil.write('Field[%d].thickness = %le;         \n' % ((NbrFld+2), thickness));
        fil.write('BoundaryLayer Field = %d;           \n' % (NbrFld+2));
        
        fil.write('Physical Line(1)  = {12, 15};                          \n');
        fil.write('Physical Line(2)  = {11};                          \n');
        fil.write('Physical Line(3)  = {10};                       \n');
        fil.write('Physical Line(4)  = {7, 8};                        \n');
        fil.write('Physical Line(5)  = {6};                           \n');
        fil.write('Physical Line(6)  = {3, 4, 5};                     \n');
        fil.write('Physical Line(7)  = {1, 2};                        \n');
        fil.write('Physical Line(8)  = {13};                           \n');
        
        fil.write('Physical Surface(21) = {14};                          \n');        
    
    else :
    
        fil.write('Delete {                                                  \n');
        fil.write('  Surface{14};                                            \n');
        fil.write('}                                                         \n');
        fil.write('Line(16) = {%d, 2};                                       \n'%exit_vid);
        fil.write('Line Loop(17) = {16, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12}; \n');
        fil.write('Plane Surface(18) = {17};                                 \n');
        fil.write('Line Loop(19) = {1, -16, 15, 13};                         \n');
        fil.write('Plane Surface(20) = {19};                                 \n');
        fil.write('Physical Surface(21) = {20, 18};                          \n');
        
        
        #fil.write('Physical Line(1)  = {1};                             \n');
        #fil.write('Physical Line(2)  = {2};                             \n');
        #fil.write('Physical Line(3)  = {3};                             \n');
        #fil.write('Physical Line(4)  = {4};                             \n');
        #fil.write('Physical Line(5)  = {5};                             \n');
        #fil.write('Physical Line(6)  = {6};                             \n');
        #fil.write('Physical Line(7)  = {7};                             \n');
        #fil.write('Physical Line(8)  = {8};                             \n');
        #fil.write('Physical Line(9)  = {9};                             \n');
        #fil.write('Physical Line(10) = {10};                            \n');
        #fil.write('Physical Line(11) = {11};                            \n');
        #fil.write('Physical Line(12) = {12,15};                         \n');
        #fil.write('Physical Line(13) = {13};                            \n');
        #fil.write('Physical Line(14) = {14};                            \n');
        
        
        fil.write('Physical Line(1)  = {12, 15};                          \n');
        fil.write('Physical Line(2)  = {11};                          \n');
        fil.write('Physical Line(3)  = {10};                       \n');
        fil.write('Physical Line(4)  = {7, 8};                        \n');
        fil.write('Physical Line(5)  = {6};                           \n');
        fil.write('Physical Line(6)  = {3, 4, 5};                     \n');
        fil.write('Physical Line(7)  = {1, 2};                        \n');
        fil.write('Physical Line(8)  = {13};                           \n');
        fil.write('Physical Line(9)  = {16};                           \n');
    
    
    fil.close();
    
    fil = open("%s.opt"%FilNam,'w');
    fil.write("Mesh.SaveElementTagType = 2;\n");
    fil.close();
    


def ExitGeoFile(exitNam, Mesh_options):
    
    # --- Options
    
    xwall  = Mesh_options.xwall;
    ywall  = Mesh_options.ywall;
    hl     = Mesh_options.hl;
    method = Mesh_options.method;
    
    ds        =  Mesh_options.ds;       
    ratio     =  Mesh_options.ratio;   
    thickness =  Mesh_options.thickness;
    
    x_thrust  = Mesh_options.x_thrust;
    
    nx=len(xwall);
    
    crdThrust=[0.0,0.0]
    for i in range(0,nx) :
        if ( math.fabs(x_thrust-xwall[nx-i-1]) < 1e-6 ):
            crdThrust = [xwall[nx-i-1], ywall[nx-i-1]]
    
    xinlet = xwall[0];
    xoutlet = xwall[nx-1];
    
    sizExit = 0.2*(xoutlet-crdThrust[0]);
    
    try:
        filExit = open(exitNam, 'w');
    except:
        sys.stderr.write("  ## ERROR : Could not open %s\n" % exitNam);
        sys.exit(0);

    sys.stdout.write("%s OPENED.\n" % exitNam);
    
    filExit.write('Point(1) = {%lf, %lf, 0, %lf};\n' % (crdThrust[0], crdThrust[1], sizExit));
    filExit.write('Point(2) = {%lf, %lf, 0, %lf};\n' % (crdThrust[0], 0.0, sizExit));
    
    x0 = 2*crdThrust[0]-xoutlet;
    
    filExit.write('Point(3) = {%lf, %lf, 0, %lf};\n' % (x0, crdThrust[1], sizExit));
    filExit.write('Point(4) = {%lf, %lf, 0, %lf};\n' % (x0, 0.0, sizExit));
    
    filExit.write('Point(5) = {%lf, %lf, 0, %lf};\n' % (xoutlet, crdThrust[1], sizExit));
    filExit.write('Point(6) = {%lf, %lf, 0, %lf};\n' % (xoutlet, 0.0, sizExit));
    
    filExit.write('Line(1)  = {1, 2};\n');
    filExit.write('Line(2)  = {3, 4};\n');
    filExit.write('Line(3)  = {5, 6};\n');
    filExit.write('Line(4)  = {3, 1};\n');
    filExit.write('Line(5)  = {1, 5};\n');
    filExit.write('Line(6)  = {4, 2};\n');
    filExit.write('Line(7)  = {2, 6};\n');
    
    filExit.write('Line Loop(1) = {1,-6,-2,4};\n');
    filExit.write('Plane Surface(1) = {1};\n');
    
    filExit.write('Line Loop(2) = {5,3,-7,-1};\n');
    filExit.write('Plane Surface(2) = {2};\n');
    
    filExit.write('Physical Line(1)  = {1};                             \n');
    filExit.write('Physical Line(2)  = {2};                             \n');
    filExit.write('Physical Line(3)  = {3};                             \n');
    filExit.write('Physical Line(4)  = {4};                             \n');
    filExit.write('Physical Line(5)  = {5};                             \n');
    filExit.write('Physical Line(6)  = {6};                             \n');
    filExit.write('Physical Line(7)  = {7};                             \n');
    
    filExit.write('Physical Surface(1) = {1};                          \n');
    filExit.write('Physical Surface(2) = {2};                          \n');
    
    filExit.close();
    
    filExit = open("%s.opt"%exitNam,'w');
    filExit.write("Mesh.SaveElementTagType = 2;\n");
    filExit.close();


def NozzleGeoFileRoundedEdges(FilNam, Mesh_options):
    
    # --- Options
    
    xwall  = Mesh_options.xwall;  # x-coordinates of points along inner nozzle wall
    ywall  = Mesh_options.ywall;  # r-coordinates of points along inner nozzle wall
    hl     = Mesh_options.hl;     # list of 5 characterisitic element lengths
    method = Mesh_options.method; # string; method of calculation
    
    ds        =  Mesh_options.ds; # 7e-06   
    ratio     =  Mesh_options.ratio; # 1.3
    thickness =  Mesh_options.thickness; # 0.02
    
    tet       = 0.012; # trailing edge thickness of nozzle
    dl        = 0.05;  # maximum length between which rounding is done
    
    x_thrust  = Mesh_options.x_thrust; # x-coordinate of vertical line where thrust is calculated
    
    # --- Domain definition
    
    nx = len(xwall);
    
    length = xwall[nx-1];
    
    CrdBox = [[0 for x in range(2)] for y in range(9)] 
    
    CrdBox[1][0] = 0;          CrdBox[1][1] = 0;
    CrdBox[2][0] = length;     CrdBox[2][1] = 0;
    #    CrdBox[3][0] = 1.5;        CrdBox[3][1] = 0;
    CrdBox[3][0] = 4.5;        CrdBox[3][1] = 0; 
    #    CrdBox[4][0] = 1.5;        CrdBox[4][1] = 2.5;
    CrdBox[4][0] = 4.5;        CrdBox[4][1] = 2.5;
    CrdBox[5][0] = -0.67;      CrdBox[5][1] = 2.5;
    CrdBox[6][0] = -0.67;      CrdBox[6][1] = 0.4244;
    CrdBox[7][0] = 0.1548;     CrdBox[7][1] = 0.4244;
    CrdBox[8][0] = length;     CrdBox[8][1] = ywall[nx-1]+tet;
    
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
    sizCur = 0.005;
        
    exit_vid = -1;
    
    # --- Write bounding points
    
    fil.write('Point(1) = {%lf, %lf, 0, %lf};\n' % (CrdBox[1][0], CrdBox[1][1], sizWal));
    fil.write('Point(2) = {%lf, %lf, 0, %lf};\n' % (x_thrust,     CrdBox[2][1], sizWal));
    fil.write('Point(3) = {%lf, %lf, 0, %lf};\n' % (CrdBox[3][0], CrdBox[3][1], sizSym));
    fil.write('Point(4) = {%lf, %lf, 0, %lf};\n' % (CrdBox[4][0], CrdBox[4][1], sizFar));
    fil.write('Point(5) = {%lf, %lf, 0, %lf};\n' % (CrdBox[5][0], CrdBox[5][1], sizFar));
    fil.write('Point(6) = {%lf, %lf, 0, %lf};\n' % (CrdBox[6][0], CrdBox[6][1], sizWal));
    
    # Begin first round on top surface of nozzle
    fil.write('Point(7) = {%lf, %lf, 0, %lf};\n' % (CrdBox[7][0]-dl, CrdBox[7][1], sizWal));
    fil.write('Point(8) = {%lf, %lf, 0, %lf};\n' % (CrdBox[7][0], CrdBox[7][1], sizWal));
    alpha = np.abs(math.atan2((CrdBox[8][1]-CrdBox[7][1]),(CrdBox[8][0]-CrdBox[7][0])));
    dx = dl*math.cos(alpha);
    dy = dl*math.sin(alpha);
    fil.write('Point(9) = {%lf, %lf, 0, %lf};\n' % (CrdBox[7][0]+dx, CrdBox[7][1]-dy, sizWal));    
    
    # Begin second round on top corner of nozzle
    fil.write('Point(10) = {%lf, %lf, 0, %lf};\n' % (CrdBox[8][0]-dx, CrdBox[8][1]+dy, sizWal));        
    fil.write('Point(11) = {%lf, %lf, 0, %lf};\n' % (CrdBox[8][0], CrdBox[8][1], sizCur));
    dy2 = tet/2; #min(dy,tet/2);
    print 'dy2: %f' % dy2
    fil.write('Point(12) = {%lf, %lf, 0, %lf};\n' % (CrdBox[8][0], CrdBox[8][1]-dy2, sizCur));
    
    # Begin third round on bottom corner of nozzle
    #fil.write('Point(13) = {%lf, %lf, 0, %lf};\n' % (CrdBox[8][0], CrdBox[8][1]-dy2, sizCur));
    fil.write('Point(14) = {%lf, %lf, 0, %lf};\n' % (length, CrdBox[8][1]-tet, sizCur));
    dx3 = min(dl,(length-x_thrust)/2);
    x3 = CrdBox[8][0]-dx3;
    y3 = np.interp(x3,xwall,ywall);
    fil.write('Point(15) = {%lf, %lf, 0, %lf};\n' % (x3, y3, sizCur));
    
    # Points to govern field specification
    fil.write('Point(16)  = {%lf, %lf, 0, %lf};\n'% (CrdBox[3][0], CrdBox[8][1], sizWal));
    fil.write('Point(17) = {%lf, %lf, 0, %lf};\n'% (CrdBox[3][0], CrdBox[7][1]+0.25*CrdBox[8][1], sizWal));
    fil.write('Point(18) = {%lf, %lf, 0, %lf};\n'% (CrdBox[6][0], CrdBox[7][1]+0.25*CrdBox[8][1], sizWal));
    
    # --- Add B-Spline control points
    
    vid = 18;
    for i in range(0,nx) :
    
        if ( xwall[nx-i-1] > x3 ):
            pass;
        else:
            fil.write('Point(%d)  = {%lf, %lf, 0, %lf};\n'% (vid, xwall[nx-i-1], ywall[nx-i-1], sizWal));
    
        vid=vid+1;        
        if ( math.fabs(x_thrust-xwall[nx-i-1]) < 1e-6 ):
            exit_vid = vid;
            print "x %lf exit_vid = %d" % (xwall[nx-i-1], exit_vid)
    
    if ( exit_vid == -1 ):
        print " ## ERROR : Coordinates for the thrust computation don't match.";
        sys.exit(1);
    
    # --- Add domain lines
    
    fil.write('Line(1)  = {1, 2};\n');
    fil.write('Line(2)  = {2, 3};\n');
    fil.write('Line(3)  = {3, 16};\n');
    fil.write('Line(4)  = {16, 17};\n');
    fil.write('Line(5)  = {17, 4};\n');
    fil.write('Line(6)  = {4, 5};\n');
    fil.write('Line(7)  = {5, 18};\n');
    fil.write('Line(8)  = {18, 6};\n');
    fil.write('Line(9)  = {6, 7};\n');
    fil.write('Line(10) = {9, 10};\n');
    #fil.write('Line(11) = {12, 13};\n');
    
    # --- Define B-splines for rounds
    fil.write('BSpline(12) = {7, 8, 9};\n');
    fil.write('BSpline(13) = {10, 11, 12};\n');
    #fil.write('BSpline(14) = {13, 14, 15};\n');
    fil.write('BSpline(14) = {12, 14, 15};\n');
    
    # --- Define B-Spline for inner wall shape
    
    # B-spline from exit to thrust integration line
    fil.write('BSpline(15) = { 15');
    eid=18+len([e for e in xwall if e > x3]); # some points are skipped
    for i in range(eid, exit_vid+1):
        fil.write(', %d' % i);
    fil.write('};\n');
    
    # B-spline from thrust integration line to inlet    
    fil.write('BSpline(16) = { %d'%(exit_vid));    
    for i in range(exit_vid+1, nx+18):
        fil.write(', %d' % i);
    fil.write('};\n');
    
    savIdx = i;
    
    fil.write('Line(17) = {%d, 1};\n' % (i));
    
    # --- Plane surface
    
    #fil.write('Line Loop(18) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 10, 13, 11, 14, 15, 16, 17};\n');
    fil.write('Line Loop(18) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 10, 13, 14, 15, 16, 17};\n');
    fil.write('Plane Surface(18) = {18};\n');
    
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
              [-100, length+0.04, -10, 0.5, hl4],\
              [length-tet/2, length+tet,CrdBox[8][1]-tet*1.5,CrdBox[8][1]+0.5*tet, hl4/3]];
    
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
        fil.write('Field[%d].EdgesList = {9,10,11,12,15}; \n' % (NbrFld+2));
        fil.write('Field[%d].NodesList = {6,%d};      \n' % ((NbrFld+2), savIdx));
        fil.write('Field[%d].hfar = 1;                 \n' % (NbrFld+2));
        fil.write('Field[%d].hwall_n = %le;       \n' % ((NbrFld+2), ds));
        fil.write('Field[%d].hwall_t = 0.300000;       \n' % (NbrFld+2));
        fil.write('Field[%d].ratio = %le;              \n' % ((NbrFld+2),ratio));
        fil.write('Field[%d].thickness = %le;         \n' % ((NbrFld+2), thickness));
        fil.write('BoundaryLayer Field = %d;           \n' % (NbrFld+2));
        
        
        #fil.write('Physical Line(1)  = {1};                             \n');
        #fil.write('Physical Line(2)  = {2};                             \n');
        #fil.write('Physical Line(3)  = {3};                             \n');
        #fil.write('Physical Line(4)  = {4};                             \n');
        #fil.write('Physical Line(5)  = {5};                             \n');
        #fil.write('Physical Line(6)  = {6};                             \n');
        #fil.write('Physical Line(7)  = {7};                             \n');
        #fil.write('Physical Line(8)  = {8};                             \n');
        #fil.write('Physical Line(9)  = {9};                             \n');
        #fil.write('Physical Line(10) = {10};                            \n');
        #fil.write('Physical Line(11) = {11};                            \n');
        #fil.write('Physical Line(12) = {12,15};                         \n');
        #fil.write('Physical Line(13) = {13};                            \n');
        #fil.write('Physical Line(14) = {14};                            \n');
        
        fil.write('Physical Line(1)  = {12, 15};                          \n');
        fil.write('Physical Line(2)  = {11};                          \n');
        fil.write('Physical Line(3)  = {9, 10};                       \n');
        fil.write('Physical Line(4)  = {7, 8};                        \n');
        fil.write('Physical Line(5)  = {6};                           \n');
        fil.write('Physical Line(6)  = {3, 4, 5};                     \n');
        fil.write('Physical Line(7)  = {1, 2};                        \n');
        fil.write('Physical Line(8)  = {13};                           \n');
        #fil.write('Physical Line(9)  = {19};                           \n');
        
        fil.write('Physical Surface(21) = {14};                          \n');
    
    
    else :
    
        fil.write('Delete {                                                  \n');
        fil.write('  Surface{18};                                            \n');
        fil.write('}                                                         \n');
        
        # Thrust integration line
        fil.write('Line(19) = {%d, 2};                                       \n'% exit_vid);
        
        # All flow exterior of nozzle (right of thrust integration line)
        #fil.write('Line Loop(20) = {19, 2, 3, 4, 5, 6, 7, 8, 9, 12, 10, 13, 11, 14, 15}; \n');
        fil.write('Line Loop(20) = {19, 2, 3, 4, 5, 6, 7, 8, 9, 12, 10, 13, 14, 15}; \n');
        fil.write('Plane Surface(21) = {20};                                 \n');
        
        # All flow internal to nozzle (left of thrust integration line)
        fil.write('Line Loop(22) = {1, -19, 16, 17};                         \n');
        fil.write('Plane Surface(23) = {22};                                 \n');
        
        # Make both surfaces above physical
        fil.write('Physical Surface(24) = {21, 23};                          \n');
        
        # Inner nozzle wall
        fil.write('Physical Line(1)  = {15, 16};                          \n');
        
        # Trailing edge of nozzle
        #fil.write('Physical Line(2)  = {13,11,14};                          \n');
        fil.write('Physical Line(2)  = {13,14};                          \n');
        
        # Exterior of nozzle
        fil.write('Physical Line(3)  = {9, 12, 10};                       \n');
        
        # Inlet of ambient flow
        fil.write('Physical Line(4)  = {7, 8};                        \n');
        
        # Ceiling
        fil.write('Physical Line(5)  = {6};                           \n');
        
        # Outlet of CFD domain
        fil.write('Physical Line(6)  = {3, 4, 5};                     \n');
        
        # Axis of symmetry
        fil.write('Physical Line(7)  = {1, 2};                        \n');
        
        # Inlet of nozzle
        fil.write('Physical Line(8)  = {17};                           \n');
        
        
        fil.write('Physical Line(9)  = {19};                           \n');
    
    
    fil.close();
    
    fil = open("%s.opt"%FilNam,'w');
    fil.write("Mesh.SaveElementTagType = 2;\n");
    fil.close();



