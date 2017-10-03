import os, time, sys, shutil, copy, math
import numpy as np
#import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
from scipy import interpolate as itp
import subprocess
from multif import _meshutils_module

def WriteGeo(FilNam, Ver, Spl, Lin, Loo, Phy, Siz, Bak, dim):
    
    NbrVer = len(Ver)-1;
    NbrSpl = len(Spl)-1;
    NbrLin = len(Lin)-1;
    NbrLoo = len(Loo)-1;
    NbrSiz = len(Siz);
    NbrPhy = len(Phy)-1;
    NbrBak = len(Bak)-1;
    
    try:
        fil = open(FilNam, 'w');
    except:
        sys.stderr.write("  ## ERROR : Could not open %s\n" % FilNam);
        sys.exit(0);
    
    sys.stdout.write("%s OPENED.\n" % FilNam);
    
    for i in range(0,NbrSiz):
        fil.write("cl%d=%lf;\n" % (i+1, Siz[i]));
    
    for i in range(1,NbrVer+1):
        fil.write("Point(%d) = {%lf,%lf,%lf,cl%d};\n" % (i, Ver[i][0], Ver[i][1], Ver[i][2], Ver[i][3]));
        
    for i in range(1,NbrSpl+1):
        nbv = len(Spl[i]);
        if ( nbv < 2 ) :
            sys.stdout.write("  ## ERROR : Invalid spline %d. Skip.\n" % i);
        
        if ( nbv == 2 ):
            kwd = "Line"
        else:
            kwd = "Spline"
        
        fil.write("%s(%d) = {%d" % (kwd,i,Spl[i][0]));
        for j in range(1,len(Spl[i])):
            fil.write(", %d" % Spl[i][j])
        fil.write("};\n");
    
    for i in range(1,NbrLoo+1):
        nbl = len(Loo[i]);
        
        if ( nbl < 3 ) :
            sys.stdout.write("  ## ERROR : Invalid line loop %d. Skip.\n" % i);
            
        if ( nbl <= 4 ):
            kwd = "Ruled"
        else:
            kwd = "Plane"
        
        fil.write( "Line Loop(%d) = {%d" % (i, Loo[i][0]));
        for j in range(1,nbl):
            fil.write(", %d" % Loo[i][j]);
        fil.write("};\n");    
        fil.write( "%s Surface(%d) = {%d};\n" % (kwd,i, i))
        
        
    # --------------------
    # --- 2D : physical lines + surf
    # --------------------
    
    if ( dim == 2 ):
        for i in range(1,NbrPhy+1):
            nbs = len(Phy[i]);
    
            if ( nbs < 1 ) :
                sys.stdout.write("  ## ERROR : Invalid physical group %d. Skip.\n" % i);
    
            fil.write( "Physical Line(%d) = {%d" % (i, Phy[i][0]));
            for j in range(1,nbs):
                fil.write(", %d" % Phy[i][j]);
            fil.write("};\n");
    
        if ( NbrLoo >= 1 ):    
            fil.write( "Physical Surface(1) = {1");
            for i in range(2,NbrLoo+1):
                fil.write(", %d" % i);
            fil.write("};\n");    
            
    
    
    # --------------------
    # --- 3D : physical surfaces + volume
    # --------------------
    
    if ( dim == 3 ):
        
        for i in range(1,NbrPhy+1):
            nbs = len(Phy[i]);
    
            if ( nbs < 1 ) :
                sys.stdout.write("  ## ERROR : Invalid physical group %d. Skip.\n" % i);
    
            fil.write( "Physical Surface(%d) = {%d" % (i, Phy[i][0]));
            for j in range(1,nbs):
                fil.write(", %d" % Phy[i][j]);
            fil.write("};\n");
        
        fil.write("Surface Loop(1) = { 1");
        for i in range(2,NbrLoo+1):
            fil.write(", %d" % i);
        fil.write("};\n");    
        fil.write("Volume(1) = {1};\n");
        fil.write("Physical Volume(1) = {1};\n");
    
    for i in range(1,NbrBak+1):
        fil.write( "Field[%d] = Box;       \n" % (i));
        fil.write( "Field[%d].VIn  = %lf;  \n" % (i, Bak[i][0]));
        fil.write( "Field[%d].VOut = %lf;  \n" % (i, Bak[i][1]));
        fil.write( "Field[%d].XMax = %lf;  \n" % (i, Bak[i][2]));
        fil.write( "Field[%d].XMin = %lf;  \n" % (i, Bak[i][3]));
        fil.write( "Field[%d].YMax = %lf;  \n" % (i, Bak[i][4]));
        fil.write( "Field[%d].YMin = %lf;  \n" % (i, Bak[i][5]));
        fil.write( "Field[%d].ZMax = %lf;  \n" % (i, Bak[i][6]));
        fil.write( "Field[%d].ZMin = %lf;  \n" % (i, Bak[i][7]));
        #fil.write( "Background Field = %d; \n" % (i));
    
    fil.write("Field[%d] = Min;              \n" % (NbrBak+1));
    fil.write("Field[%d].FieldsList = {%d" % (NbrBak+1, 1));
    for i in range(2,NbrBak+1):
        fil.write(", %d" % i);
    fil.write("};\n");    
    
    fil.write("Background Field = %d;        \n" % (NbrBak+1));
    
    fil.close();
    
    # --------------------
    # --- Write options
    # --------------------
    
    OptNam = "%s.opt" % FilNam;
    
    try:
        fil = open(OptNam, 'w');
    except:
        sys.stderr.write("  ## ERROR : Could not open %s\n" % OptNam);
        sys.exit(0);
        
    sys.stdout.write("%s OPENED.\n" % OptNam);
    
    fil.write("Mesh.SaveElementTagType = 2;\n");
    
    fil.close();    


# Function which approximates 2D axisymmetric nozzle with equivalent areas as
# 3D nonaxisymmetric nozzle. Requires the following inputs:
# x = vector of axial coordinates spanning length of nozzle
# fr1 = function which returns major axis radius given a value of x
# fr2 = function which returns minor axis radius given a value of x
# fz = function which returns z location of centerline given a value of x
# params = a list with the following components:
#          [z coord of cut at the throat, z coord of the flat exit, 0, 
#           x-coordinate of inlet, x-coordinate of exit]
def MF_GetRadius_old (x, fr1, fr2, fz, params):
	
	zcut0   = params[0]; # z crd of the cut at the throat
	zcut1   = params[1]; # z crd of the flat exit (bottom of the shovel)
	xthroat = params[3];
	xexit   = params[4];
	
	if ( isinstance(x, (list, tuple, np.ndarray)) ):
		rad = np.zeros(len(x));
		for i in range(0,len(x)):
			rad[i] = MF_GetRadius (x[i], fr1, fr2, fz, params);
		return rad;
	
	r1 = fr1(x);
	r2 = fr2(x);	
	
	z0  = fz(xthroat);
	r10 = fr1(xthroat);
	r20 = fr2(xthroat);
	# Theta is measured down from the vertical axis
	theta0 = math.acos((zcut0-z0)/r10);
	
	z1  = fz(xexit);
	r11 = fr1(xexit);
	r21 = fr2(xexit);
	# Theta is measured down from the vertical axis
	theta1 = math.acos((zcut1-z1)/r21);
	
	ftheta = itp.interp1d([xthroat,xexit],[theta0,theta1],kind='linear')
	theta = ftheta(x);

	alp = (xexit-x)/(xexit-xthroat);
	
	area = 0.25*r1*r2*math.pi;
	area = area + 0.5*r1*r2*(theta-np.pi/2 - 0.5*math.sin(2*theta-math.pi));
	areab = 0.5*r1*r2*(math.pi - theta + 0.5*math.sin(2*theta-math.pi));
	
	area = 2*(area + alp*areab);
	
	#return area; 
	
	rad = math.sqrt(area/math.pi);
	
	return rad;
	
def MF_GetRadius (x, nozzle):
    
    if ( isinstance(x, (list, tuple, np.ndarray)) ):
        rad = np.zeros(len(x));
        for i in range(0,len(x)):
        	rad[i] = MF_GetRadius (x[i],nozzle);
        return rad;
    
    from .. import nozzle as noz
    
    geometry = noz.geometry;
    
    majoraxisTmp = geometry.Bspline(nozzle.wall.majoraxis.coefs);
    minoraxisTmp = geometry.Bspline(nozzle.wall.minoraxis.coefs);
    
    fr1 = majoraxisTmp.radius
    fr2 = minoraxisTmp.radius
    
    theta_in  = 1.572865;
    theta_out = 1.855294;
    
    coefs_center = nozzle.wall.centerline.coefs;
    
    x_in   = coefs_center[0];
    x_out  = coefs_center[len(coefs_center)/2-1];
    
    #--- Get x, r1, r2, zcenter
    
    r1 = fr1(x);
    r2 = fr2(x);
    
    alp = (x-x_in)/(x_out-x_in);
    
    theta = alp*theta_out + (1.0-alp)*theta_in;
    
    #--- Compute area
    
    area_A = 0.5*theta*r1*r2; 
    area_B = 0.5*(r1*r2*math.fabs(math.cos(theta)*math.sin(theta))); # area of triangle (orig,P(theta),orig-zcut)
    area_C = 0.5*(math.pi-theta)*r1*r2 - area_B; # area below zcut
        
    area_tot = area_A+area_B + (1.0-alp) * area_C;
    
    #--- Verif area_C
    x0 = r1*math.sin(theta);
    y0 = r2*math.cos(theta);
    #area_C_int = 0.5*r1*r2*((1.0-2*y0)*x0 - 0.5*math.sin(2*x0));
    
    #print "area_C %lf area_C_int %lf " % (area_C, area_C_int)
    #sys.exit(1);
    
    area_tot = area_A+area_B + (1.0-alp) * area_C 
    
    rad = math.sqrt(2*area_tot/math.pi);
    
    return rad;    
#
# 
#def MF_DefineAxiSymCAD (FilNam, x_inp, fr1, fr2, fz, sizes, params):
#    
#    NbrLnk = 400;
#    
#    
#    from .. import nozzle as noz
#    
#    geometry = noz.geometry;
#    
#    majoraxisTmp = geometry.Bspline(nozzle.wall.majoraxis.coefs);
#    minoraxisTmp = geometry.Bspline(nozzle.wall.minoraxis.coefs);
#    
#    fr1 = majoraxisTmp.radius
#    fr2 = minoraxisTmp.radius
#    
#    coefs_center = nozzle.wall.centerline.coefs;
#    
#    x_in   = coefs_center[0];
#    x_out  = coefs_center[len(coefs_center)/2-1];
#    
#    xmax    = params[5];
#    ymax    = params[6];
#    zmax    = params[7];
#    
#    hl0 = sizes[3];  # Edge size in the 'small' plume area
#    hl1 = sizes[4];  # Edge size in the vicinity of the aircraft
#    
#    x = np.linspace(x_in, x_out, NbrLnk);
#    #y = MF_GetRadius (x, x_inp, fr1, fr2, fz, params);
#    y = MF_GetRadius (x,nozzle);
#    
#    
#    Ver = [0];
#    Spl = [0];
#    Lin = [0];
#    Loo = [0];
#    Phy = [0];   # Physical groups
#    Bak = [0];
#    
#    y0 = 0.8; 
#    y2 = 0.8*y0;
#    thickness = params[8];
#    
#    # Point sizes : 
#    Siz = [sizes[0], sizes[1], sizes[2]] # inside nozzle, aircraft surf, farfield
#    
#    OutVid = np.zeros(100);
#    OutLid = np.zeros(100);
#    
#    # --- Inner nozzle wall
#    
#    Spl.append([]);
#    for i in range(0,NbrLnk):
#        Ver.append([x[i], y[i], 0, 1]);
#        Spl[-1].append(len(Ver)-1);
#    
#    OutLid[1] = len(Spl)-1;
#    OutVid[8] = len(Ver)-1;
#    OutVid[9] = 1;
#    
#    # --- Box
#    
#    Ver.append([0,0,0,1]);
#    OutVid[1] = len(Ver)-1;
#    
#    Ver.append([xmax,0,0,3]);
#    OutVid[2] = len(Ver)-1;
#    
#    Ver.append([xmax,zmax,0,3]);
#    OutVid[3] = len(Ver)-1;
#    
#    Ver.append([0,zmax,0,3]);
#    OutVid[4] = len(Ver)-1;
#    
#    Ver.append([0,y0,0,2]);
#    OutVid[5] = len(Ver)-1;
#    
#    Ver.append([x_in,y2,0,2]);
#    OutVid[6] = len(Ver)-1;
#    
#    Ver.append([x[-1],y[-1]+thickness,0,2]);
#    OutVid[7] = len(Ver)-1;
#    
#    Ver.append([x[-1],0,0,1]);
#    OutVid[10] = len(Ver)-1;
#    
#    # --- Add lines
#    
#    Spl.append([OutVid[8],OutVid[7]]);
#    Spl.append([OutVid[7],OutVid[6],OutVid[5]]);
#    Spl.append([OutVid[5],OutVid[4]]);
#    Spl.append([OutVid[4],OutVid[3]]);
#    Spl.append([OutVid[3],OutVid[2]]);
#    Spl.append([OutVid[2],OutVid[10]]);
#    Spl.append([OutVid[10],OutVid[1]]);
#    Spl.append([OutVid[1],OutVid[9]]);
#    
#    # --- Add surface
#    
#    Loo.append([1, 2, 3, 4, 5, 6, 7, 8, 9]);
#    for i in range(1,7):
#        Phy.append([i]); # physical lines
#    Phy.append([7,8]);
#    Phy.append([9]);
#    
#    # --- Define background size fields
#    
#    # size_in size_out xmax xmin ymax ymin zmax zmin
#    hei = 0.45*(Ver[int(OutVid[5])][1]);
#    Bak.append([hl1, 1000, 100, -100, Ver[int(OutVid[5])][1]+hei, -1, 1, -1]);
#    
#    # size_in size_out xmax xmin ymax ymin zmax zmin
#    hei = 0.3*(Ver[int(OutVid[8])][1]);
#    Bak.append([hl0, 1000, 100, Ver[int(OutVid[8])][0], Ver[int(OutVid[8])][1]+hei, -1, 1, -1]);
#    
#    # --- Write .geo file
#    
#    WriteGeo(FilNam, Ver, Spl, Lin, Loo, Phy, Siz, Bak, 2)
#
def HF_GenerateExitMesh(nozzle, xloc=-1):
    
    geonam = "nozzle_exit.geo";
    HF_ExitPlaneCAD (geonam, nozzle, xloc);
    
    gmsh_executable = 'gmsh';
    try :
        cmd = [gmsh_executable, '-2', geonam, '-o', "nozzle_exit_hin.mesh"];
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
    except:
        raise;
    
    
def HF_ExitPlaneCAD (filename, nozzle, xloc):
    
    from .. import nozzle as noz
    
    geometry = noz.geometry;
    
    majoraxisTmp = geometry.Bspline(nozzle.wall.majoraxis.coefs);
    minoraxisTmp = geometry.Bspline(nozzle.wall.minoraxis.coefs);
    centerTmp    =  geometry.Bspline(nozzle.wall.centerline.coefs);    
    
    fr1 = majoraxisTmp.radius
    fr2 = minoraxisTmp.radius
    fz  = centerTmp.radius
    
    theta_in  = 1.572865;
    theta_out = 1.855294;
    
    coefs_center = nozzle.wall.centerline.coefs;
    
    x_in   = coefs_center[0];
    x_out  = coefs_center[len(coefs_center)/2-1];
    
    #--- Get x, r1, r2, zcenter
    
    if xloc > -0.1 :
        x = xloc;
    else :
        h_in = 0.06;
        x = x_out-h_in;
    
    r1 = fr1(x);
    r2 = fr2(x);
    zcen = fz(x);
    
    alp = (x-x_in)/(x_out-x_in);
    
    NbrLnk = 400;
    
    #--- 
    
    theta = alp*theta_out + (1.0-alp)*theta_in;
    
    Ver = [0];
    Spl = [0];
    Lin = [0];
    Loo = [0];
    Phy = [0];
    Bak = [0];
    dim = 3;
    
    Siz = [0.007, 0.008];
    
    # Upper part
    ttab = np.linspace(0.0,theta, NbrLnk);
    Spl.append([])
    for i in range(len(ttab)):
        Ver.append([x,r1*math.sin(ttab[i]), zcen+r2*math.cos(ttab[i]),1]);
        Spl[-1].append(len(Ver)-1);
    
    zcut = r2*math.cos(theta);
    
    # Lower part
    ttab = np.linspace(theta, math.pi, NbrLnk);
    Spl.append([])
    Spl[-1].append(len(Ver)-1);
    
    for i in range(1,len(ttab)):        
        z = zcen+alp*zcut+(1.0-alp)*(r2*math.cos(ttab[i]));
        Ver.append([x,r1*math.sin(ttab[i]),z,1]);
        Spl[-1].append(len(Ver)-1);
        
    # Line
    Spl.append([]);
    Spl[-1] = [len(Ver)-1,1];
    
    # Loop
    Loo.append([]);
    Loo[-1] = [1,2,3];
    
    # Physical
    Phy.append([1]);
        
    WriteGeo(filename, Ver, Spl, Lin, Loo, Phy, Siz, Bak, 3) ;
    
    

def HF_DefineCAD (FilNam, x_inp, fr1, fr2, fz, sizes, params):
    
    # --- Parameters
    
    NbrLnk = 400; # How many control points for each spline?
    
    xmax    = params[5];
    ymax    = params[6];
    zmax    = params[7];
    
    Box = np.zeros([3,2]);                # Computational domain size
    Box[0][0] = 0;      Box[0][1] = xmax; 
    Box[1][0] = 0;      Box[1][1] = ymax; 
    Box[2][0] = -zmax;    Box[2][1] = zmax; 
    
    hl0 = sizes[3];  # Edge size in the 'small' plume area
    hl1 = sizes[4];  # Edge size in the vicinity of the aircraft
    
    # Point sizes : 
    Siz = [sizes[0], sizes[1], sizes[2]] # inside nozzle, aircraft surf, farfield
    
    xthroat = params[3];
    xexit   = params[4];
    thickness = params[8];
    
    # -------------------------
    
    NbrSec = len(x_inp);
    
    Ver = [0];
    Spl = [0];
    Lin = [0];
    Loo = [0];
    Phy = [0];   # Physical groups
    Bak = [0];
    
    lid = np.zeros(4)
    
    SavIdx = np.zeros([NbrSec+10,4])
    SavShovel = np.zeros([NbrLnk+2,2]);
    
    SavSpl = np.zeros([NbrSec+10,7])
    
    OutVid = np.zeros(100);
    OutLid = np.zeros(100);
    
    sec_throat = 1;
    
    # --- Compute flat exit parameters
    
    # --- z crd that define the geometry
    zcut0 = params[0]; # z crd of the cut at the throat
    zcut1 = params[1]; # z crd of the flat exit (bottom of the shovel)
    zcut2 = params[2]; # z coordinate of the top of the shovel
    
    x2 = xexit;
    z2 = fz(x2);
    r12 = fr1(x2);
    r22 = fr2(x2);
    theta2 = math.acos((zcut2-z2)/r22);
    
    x0 = xthroat;
    z0 = fz(x0);
    r10 = fr1(x0);
    r20 = fr2(x0);
    theta0 = math.acos((zcut0-z0)/r10);
    
    x1 = xexit;
    z1 = fz(x1);
    r11 = fr1(x1);
    r21 = fr2(x1);
        
    theta1 = math.acos((zcut1-z1)/r21);
    ftheta = itp.interp1d([x0,x1],[theta0,theta1],kind='linear')
    dalp = 1./(NbrLnk+1);
    
    NbrVer = 0;
    NbrEll = 0;
    NbrSpl = 0;
    
    for i in range(0,NbrSec):    
        x  = x_inp[i];
        z  = fz(x);
        
        r1 = fr1(x);
        r2 = fr2(x);
        
        
        if ( math.fabs(x-xthroat) < 1e-12 ):
            sec_throat = i;
        
        theta = theta2;
        
        # --- Center point
        
        NbrVer = NbrVer+1;
        Ver.append([x,0,z+r2,1]);
        SavIdx[i,0] = NbrVer;
            
        ynew = r1*math.sin(theta)
        znew = z+r2*math.cos(theta)
        
        NbrVer = NbrVer+1;
        Ver.append([x, ynew, znew,1]);
        SavIdx[i,1] = NbrVer;
        
        if ( i == NbrSec-1 ) :
            SavShovel[0,0] = ynew;
            SavShovel[0,1] = znew;
        
        alp = dalp;
        for j in range(0, NbrLnk):
            thetaloc = alp*theta;
            NbrVer = NbrVer+1;
            Ver.append([x, r1*math.sin(thetaloc), z+r2*math.cos(thetaloc),1]);
            alp = alp+dalp;
        
        NbrSpl = NbrSpl+1;
        Spl.append([])
        Spl[-1].append(SavIdx[i,0]);
        for j in range(0, NbrLnk):
            Spl[-1].append(NbrVer-NbrLnk+j+1);
        Spl[-1].append(SavIdx[i,1]);
        
        SavSpl[i][0] = len(Spl)-1;
        
        if ( i > 0 ):
            x0 = x_inp[i-1];
            x1 = x_inp[i];
            
            # --- Add upper spline
            
            alp = dalp;
            for j in range(0, NbrLnk):
                x = alp*x1+(1.-alp)*x0;
                z = fz(x);
                r = fr2(x);
                alp = alp+dalp;
                NbrVer = NbrVer+1;
                Ver.append([x, 0, z+r,1])
            
            NbrSpl = NbrSpl+1;
            Spl.append([])
            Spl[-1].append(SavIdx[i-1,0]);
            for j in range(0, NbrLnk):
                Spl[-1].append(NbrVer-NbrLnk+j+1);
            Spl[-1].append( SavIdx[i,0]);
            
            SavSpl[i][3] = len(Spl)-1;
            
            # --- Add lower spline
            
            r0 = fr1(x_inp[i-1]);
            r1 = fr1(x_inp[i]);
            
            alp = dalp;
            for j in range(0, NbrLnk):
                x = alp*x1+(1.-alp)*x0;
                z = fz(x);
                r = fr1(x);
                theta = theta2;
                alp = alp+dalp;
                NbrVer = NbrVer+1;
                Ver.append([x, fr1(x)*math.sin(theta), z+fr2(x)*math.cos(theta),1]);
                            
            Spl.append([])
            Spl[-1].append(SavIdx[i-1,0]+1);
            for j in range(0, NbrLnk):
                Spl[-1].append(NbrVer-NbrLnk+j+1);
            Spl[-1].append(SavIdx[i,0]+1);
            
            SavSpl[i][4] = len(Spl)-1;
    
    # ----------------------------------------
    # --- Intermediate layer for shovel
    # ----------------------------------------
    
    # -> theta goes from theta2 to ftheta(x)
    
    for i in range(0,NbrSec):
        x  = x_inp[i];
        z  = fz(x);
        r1 = fr1(x);
        r2 = fr2(x);
        
        if ( x < xthroat-1e-12 ):
            theta = ftheta(2);
        else :
            theta = ftheta(x);
        
        if ( theta2 > theta ) :
            sys.write("ERROR : WRONG THETA VALUE\n");
            sys.exit(0);
        
        ynew = r1*math.sin(theta)
        znew = z+r2*math.cos(theta)
        
        NbrVer = NbrVer+1;
        Ver.append([x, ynew, znew,1])
        
        SavIdx[i,2] = NbrVer;
        
        SavShovel[NbrLnk+1,0] = ynew;
        SavShovel[NbrLnk+1,1] = znew;
        
        alp = dalp;
        for j in range(0, NbrLnk):
            thetaloc = alp*theta + (1.-alp)*theta2;
            ynew = r1*math.sin(thetaloc)
            znew = z+r2*math.cos(thetaloc)
            NbrVer = NbrVer+1;
            Ver.append([x, ynew, znew,1]);
            
            if ( i == NbrSec-1 ):
                SavShovel[j+1,0] = ynew;
                SavShovel[j+1,1] = znew;
            
            alp = alp+dalp;
      
        NbrSpl = NbrSpl+1;
        Spl.append([])
        Spl[-1].append(SavIdx[i,1]);
        for j in range(0, NbrLnk):
            Spl[-1].append(NbrVer-NbrLnk+j+1);
        Spl[-1].append(SavIdx[i,2]);
        
        SavSpl[i][1] = len(Spl)-1;
        
        if ( i > 0 ):
            x0 = x_inp[i-1];
            x1 = x_inp[i];
            r0 = fr1(x_inp[i-1]);
            r1 = fr1(x_inp[i]);
            
            alp = dalp;
            for j in range(0, NbrLnk):
                x = alp*x1+(1.-alp)*x0;
                z = fz(x);
                r = fr1(x);
                
                if ( x < xthroat-1e-12 ):
                    theta = ftheta(2);
                else :
                    theta = ftheta(x);
        
                ynew = fr1(x)*math.sin(theta);
                znew = z+fr2(x)*math.cos(theta);
    
                alp = alp+dalp;
                NbrVer = NbrVer+1;
                Ver.append([x, ynew, znew,1]);
    
                
            NbrSpl = NbrSpl+1;
            Spl.append([])
            Spl[-1].append(SavIdx[i-1,2]);
            for j in range(0, NbrLnk):
                Spl[-1].append(NbrVer-NbrLnk+j+1);
            Spl[-1].append(SavIdx[i,2]);
    
            SavSpl[i][5] = len(Spl)-1;
            
    # --------------------
    # --- Model bottom
    # --------------------
    
    dalp = 1./(NbrLnk+1);
    
    for i in range(0,NbrSec):
        x  = x_inp[i];
        z  = fz(x);        
        r1 = fr1(x);
        r2 = fr2(x);
        
        if ( x < 2-1e-12 ):
            theta = ftheta(2);
        else :
            theta = ftheta(x);
        
        dtheta  = (math.pi-theta) / (NbrLnk+1);    
        
        thetaloc = theta+dtheta;
        for j in range(0, NbrLnk+1):
            ynew = r1*math.sin(thetaloc);
            znew = z+r2*math.cos(thetaloc);
            
            if ( x > xthroat-1e-12 ):
                alp = (x-xthroat)/(xexit-xthroat);
                
                ynew0 = (1.-(j+1)*dalp)*r1*math.sin(theta);
                znew0 = z+r2*math.cos(theta);
                
                ynew = alp*(ynew0)+(1.-alp)*ynew;
                znew = alp*(znew0)+(1.-alp)*znew;
    
            NbrVer = NbrVer+1;
            Ver.append([x, ynew, znew,1]);
            
            thetaloc = thetaloc + dtheta;
        
        SavIdx[i,3] = NbrVer;
        
        NbrSpl = NbrSpl+1;
        Spl.append([])
        Spl[-1].append(SavIdx[i,2]);
        for j in range(0, NbrLnk):
            Spl[-1].append(NbrVer-NbrLnk+j);
        Spl[-1].append(SavIdx[i,3]);
        
        SavSpl[i][2] = len(Spl)-1;
        
        if ( i > 0 ) :
            x0 = x_inp[i-1];
            x1 = x_inp[i];
            
            r0 = fr1(x_inp[i-1]);
            r1 = fr1(x_inp[i]);
            
            alp = dalp;
            for j in range(0, NbrLnk):
                x = alp*x1+(1.-alp)*x0;
                z = fz(x);
                r = fr1(x);
                
                theta = math.pi
                ynew = fr1(x)*math.sin(theta);
                znew = z+fr2(x)*math.cos(theta);
                
                if ( x > xthroat-1e-12 ) :
                    alp0 = (x-xthroat)/(xexit-xthroat);
                    znew0 = z+fr2(x)*math.cos(ftheta(x));
                    znew = alp0*(znew0)+(1.-alp0)*znew;
                    
                alp = alp+dalp;
                NbrVer = NbrVer+1;
                Ver.append([x, ynew, znew,1]);
    
            NbrSpl = NbrSpl+1;
            Spl.append([])
            Spl[-1].append(SavIdx[i-1,3]);
            for j in range(0, NbrLnk):
                Spl[-1].append(NbrVer-NbrLnk+j+1);
            Spl[-1].append(SavIdx[i,3]);
            
            SavSpl[i][6] = len(Spl)-1;
            
    # --------------------
    # --- Add surfaces
    # --------------------
    
    # --- Before throat / top
    
    Phy.append([]);
    for i in range(1, sec_throat+1):
        Loo.append([SavSpl[i][0], -SavSpl[i][4], -SavSpl[i-1][0], SavSpl[i][3]]);
        Phy[-1].append(len(Loo)-1);
        Loo.append([SavSpl[i][4], SavSpl[i][1], -SavSpl[i][5], -SavSpl[i-1][1]]);
        Phy[-1].append(len(Loo)-1);
        
    # --- Before throat / bottom
    
    Phy.append([]);
    for i in range(1, sec_throat+1):
        Loo.append([SavSpl[i][5], SavSpl[i][2], -SavSpl[i][6], -SavSpl[i-1][2]]);
        Phy[-1].append(len(Loo)-1);
    
    # --- After throat / top
    
    Phy.append([]);
    for i in range(sec_throat+1, NbrSec):
        Loo.append([SavSpl[i][0], -SavSpl[i][4], -SavSpl[i-1][0], SavSpl[i][3]]);
        Phy[-1].append(len(Loo)-1);
        Loo.append([SavSpl[i][4], SavSpl[i][1], -SavSpl[i][5], -SavSpl[i-1][1]]);
        Phy[-1].append(len(Loo)-1);
        
    # --- After throat / bottom
    
    Phy.append([]);
    for i in range(sec_throat+1, NbrSec):
        Loo.append([SavSpl[i][5], SavSpl[i][2], -SavSpl[i][6], -SavSpl[i-1][2]]);
        Phy[-1].append(len(Loo)-1);
        
    # --------------------
    # --- Add shovel
    # --------------------
    
    x_shovel = 4.6;
    
    for j in range(0,NbrLnk+2):
        y = SavShovel[j,0]
        z = SavShovel[j,1]
        
        NbrVer = NbrVer+1;
        Ver.append([ x_shovel, y, z,1]);
    
    SavIdx[NbrSec][0] = NbrVer-NbrLnk-1;
    SavIdx[NbrSec][1] = NbrVer;
    
    NbrSpl = NbrSpl+1;
    Spl.append([])
    Spl[-1].append(NbrVer-NbrLnk-1);
    for j in range(0, NbrLnk):
        Spl[-1].append(NbrVer-NbrLnk-2+j+1);
    Spl[-1].append(NbrVer);
    
    SavSpl[NbrSec][0] = len(Spl)-1;
    
    # --- Tip
    
    NbrVer = NbrVer+1;
    xnew = x_shovel+0.2;
    ynew = 0;
    znew = zcut1;
    Ver.append([xnew, ynew, znew,1]);
    
    SavIdx[NbrSec][2] = NbrVer;
    
    NbrSpl = NbrSpl+1;
    Spl.append([ SavIdx[NbrSec-1][1], SavIdx[NbrSec][0]]);
    SavSpl[NbrSec][2] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([SavIdx[NbrSec-1][2], SavIdx[NbrSec][1]]);
    SavSpl[NbrSec][3] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([SavIdx[NbrSec][1], SavIdx[NbrSec][2]]);
    SavSpl[NbrSec][1] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([SavIdx[NbrSec-1][3], SavIdx[NbrSec][2]]);
    SavSpl[NbrSec][4] = len(Spl)-1;
    
    
    lid[0] = SavSpl[NbrSec-1][2];
    lid[1] = SavSpl[NbrSec][4];
    lid[2] = -SavSpl[NbrSec][1];
    lid[3] = -SavSpl[NbrSec][3];
    Loo.append([lid[0], lid[1], lid[2], lid[3]]);
        
    Phy.append([len(Loo)-1])
    
    lid[0] = SavSpl[NbrSec-1][1];
    lid[1] = SavSpl[NbrSec][3];
    lid[2] = -SavSpl[NbrSec][0];
    lid[3] = -SavSpl[NbrSec][2];
    Loo.append([lid[0], lid[1], lid[2], lid[3]]);
    
    Phy.append([len(Loo)-1])
    
    # --------------------
    # --- Add outside of shovel
    # --------------------
    
    x = x_shovel;
    z  = fz(x_inp[NbrSec-1]);
    
    r1 = fr1(x_inp[NbrSec-1])+ thickness;
    r2 = fr2(x_inp[NbrSec-1])+ thickness;
    
    theta = theta1;
    
    if ( theta2 > theta ) :
        sys.write("ERROR : WRONG THETA VALUE\n");
        sys.exit(0);
        
    NbrVer = NbrVer+1;
    ynew = r1*math.sin(theta2)
    znew = z+r2*math.cos(theta2)
    Ver.append([x, ynew, znew,1])
        
    OutVid[1] = len(Ver)-1;
    
    ynew = r1*math.sin(theta)
    znew = z+r2*math.cos(theta)
    NbrVer = NbrVer+1;
    Ver.append([x, ynew, znew,1])
    
    OutVid[2] = len(Ver)-1;
    
    alp = dalp;
    for j in range(0, NbrLnk):
        thetaloc = alp*theta + (1.-alp)*theta2;
        ynew = r1*math.sin(thetaloc)
        znew = z+r2*math.cos(thetaloc)
        NbrVer = NbrVer+1;
        Ver.append([x, ynew, znew,1])
        alp = alp+dalp;
    
    NbrSpl = NbrSpl+1;
    
    Spl.append([])
    Spl[-1].append(OutVid[1]);
    for j in range(0, NbrLnk):
        Spl[-1].append(NbrVer-NbrLnk+j+1);
    Spl[-1].append(OutVid[2]);
    
    OutLid[1] = len(Spl)-1;
    
    NbrVer = NbrVer+1;
    xnew = x_shovel+0.2;
    ynew = 0;
    idx = int(OutVid[2]);
    znew = Ver[idx][2];
    Ver.append([xnew, ynew, znew,1]);
    
    OutVid[3] = len(Ver)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[2], OutVid[3]]);
    
    OutLid[2] = len(Spl)-1;
    
    # --- Vertical line tip
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[3], SavIdx[NbrSec,2]]);
    
    OutLid[3] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[1], SavIdx[NbrSec,0]]);
    
    OutLid[4] = len(Spl)-1;
        
    # --------------------
    # --- Add outside 
    # --------------------
    
    x  = x_inp [NbrSec-1];
    z = fz(x);    
    r1 = fr1(x_inp[NbrSec-1])+ thickness;
    r2 = fr2(x_inp[NbrSec-1])+ thickness;
    
    theta = theta2;
    
    NbrVer = NbrVer+1;
    Ver.append([x, 0, z+r2,1])
    SavIdx[NbrSec+1,0] = NbrVer;
    OutVid[5] = len(Ver)-1;
        
    ynew = (r1)*math.sin(theta)
    znew = z+(r2)*math.cos(theta)
    NbrVer = NbrVer+1;
    Ver.append([x, ynew, znew,1]);
    
    OutVid[4] = len(Ver)-1;
    
    alp = dalp;
    for j in range(0, NbrLnk):
        thetaloc = alp*theta;
        ynew = (r1)*math.sin(thetaloc)
        znew = z+(r2)*math.cos(thetaloc)
        NbrVer = NbrVer+1;
        Ver.append([x, ynew, znew,1]);
        alp = alp+dalp;
    
    NbrSpl = NbrSpl+1;
    Spl.append([])
    Spl[-1].append(OutVid[5]);
    for j in range(0, NbrLnk):
        Spl[-1].append(NbrVer-NbrLnk+j+1);
    Spl[-1].append(OutVid[4]);
    
    OutLid[5] = len(Spl)-1;
        
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[1], OutVid[4]]);
    
    OutLid[6] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[4], SavIdx[NbrSec-1][1]]);
    
    OutLid[7] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[5], SavIdx[NbrSec-1][0]]);
    
    OutLid[8] = len(Spl)-1;
    
    # -----
        
    r1 = fr1(x_inp[NbrSec-1])+ thickness;
    r2 = fr2(x_inp[NbrSec-1])+ thickness;
    z = fz(x_inp[NbrSec-1]);
    
    theta = theta1;
    
    xnew = x_inp[NbrSec-1];
    ynew = r1*math.sin(theta)
    znew = z+r2*math.cos(theta)
        
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[6] = int(len(Ver)-1);
    
    alp = dalp;
    for j in range(0, NbrLnk):
        thetaloc = alp*theta + (1.-alp)*theta2;
        ynew = r1*math.sin(thetaloc)
        znew = z+r2*math.cos(thetaloc)
        NbrVer = NbrVer+1;
        Ver.append([xnew, ynew, znew,2])
        alp = alp+dalp;
    
    NbrSpl = NbrSpl+1;
    
    Spl.append([])
    Spl[-1].append(OutVid[4]);
    for j in range(0, NbrLnk):
        Spl[-1].append(NbrVer-NbrLnk+j+1);
    Spl[-1].append(OutVid[6]);
    
    OutLid[9] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[6], OutVid[2]]);
    
    OutLid[10] = len(Spl)-1;
    
    NbrVer = NbrVer+1;
    idx = int(OutVid[6]);
    znew = Ver[idx][2];
    Ver.append([xnew, 0, znew,2])
    OutVid[7] = len(Ver)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[6], OutVid[7]]);
    
    OutLid[11] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[7], OutVid[3]]);
    
    OutLid[12] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([SavIdx[NbrSec][1], OutVid[2]]);
    
    OutLid[13] = len(Spl)-1;
    
    #--- Add surfaces for the outside of the shovel
    
    Loo.append([OutLid[11], OutLid[12], -OutLid[2], -OutLid[10]]);
    Phy.append([len(Loo)-1])
    
    Loo.append([OutLid[2], OutLid[3], -SavSpl[NbrSec][1], OutLid[13]]);
    Phy.append([len(Loo)-1])
    
    Loo.append([-SavSpl[NbrSec][0], -OutLid[4], OutLid[1], -OutLid[13]]);
    Phy.append([len(Loo)-1])
    
    Loo.append([OutLid[10], -OutLid[1], OutLid[6], OutLid[9]]);
    Phy.append([len(Loo)-1])
    
    Loo.append([SavSpl[NbrSec][2], -OutLid[4], OutLid[6], OutLid[7]]);
    Phy.append([len(Loo)-1])
    
    Loo.append([OutLid[7], -SavSpl[NbrSec-1][0], -OutLid[8], OutLid[5]]);
    Phy.append([len(Loo)-1])
    
    # --- Outside aircraft / x=2 section
    
    x  = 2;
    z  = fz(xexit);
    r2 = 0.8
    r1 = fr1(xexit)*r2/fr2(xexit);
    
    xnew = x;
    ynew = 0;
    znew = z+r2;
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[8] = int(len(Ver)-1);
    
    xnew = x;
    ynew = r1*math.sin(theta1)
    znew = z+r2*math.cos(theta1)
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[9] = int(len(Ver)-1);
    
    xnew = x;
    ynew = r1*math.sin(theta2)
    znew = z+r2*math.cos(theta2)
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[13] = int(len(Ver)-1);
    
    xnew = x;
    ynew = 0
    znew = z-r2
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[15] = int(len(Ver)-1);
    
    # --- Outside aircraft / x=0 section
    
    x  = 0;
    z  = fz(xexit);
    r2 = 1.0;
    r1 = fr1(xexit)*r2/fr2(xexit);
    
    xnew = x;
    ynew = 0;
    znew = z+r2;
        
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[10] = int(len(Ver)-1);
    
    alp = dalp;
    for j in range(0, NbrLnk):
        thetaloc = alp*theta2;
        ynew = r1*math.sin(thetaloc)
        znew = z+r2*math.cos(thetaloc)
        NbrVer = NbrVer+1;
        Ver.append([xnew, ynew, znew,2])
        alp = alp+dalp;
    
    ynew = r1*math.sin(theta2)
    znew = z+r2*math.cos(theta2)  
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[11] = int(len(Ver)-1);
    
    NbrSpl = NbrSpl+1;
    
    Spl.append([])
    Spl[-1].append(OutVid[10]);
    for j in range(0, NbrLnk):
        Spl[-1].append(NbrVer-NbrLnk+j+1);
    Spl[-1].append(OutVid[11]);
    
    OutLid[14] = len(Spl)-1;
    
    # --- theta2 -> theta1
    
    alp = dalp;
    for j in range(0, NbrLnk):
        thetaloc = alp*theta1+(1.-alp)*theta2;
        ynew = r1*math.sin(thetaloc)
        znew = z+r2*math.cos(thetaloc)
        NbrVer = NbrVer+1;
        Ver.append([xnew, ynew, znew,2])
        alp = alp+dalp;
    
    
    ynew = r1*math.sin(theta1)
    znew = z+r2*math.cos(theta1)  
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[12] = int(len(Ver)-1);
    
    NbrSpl = NbrSpl+1;
    Spl.append([])
    Spl[-1].append(OutVid[11]);
    for j in range(0, NbrLnk):
        Spl[-1].append(NbrVer-NbrLnk+j+1);
    Spl[-1].append(OutVid[12]);
    
    OutLid[15] = len(Spl)-1;
    
    # --- theta1 -> pi
    
    alp = dalp;
    for j in range(0, NbrLnk):
        thetaloc = alp*math.pi+(1.-alp)*theta1;
        ynew = r1*math.sin(thetaloc)
        znew = z+r2*math.cos(thetaloc)
        NbrVer = NbrVer+1;
        Ver.append([xnew, ynew, znew,2])
        alp = alp+dalp;
    
    ynew = 0
    znew = z-r2 
    NbrVer = NbrVer+1;
    Ver.append([xnew, ynew, znew,2])
    OutVid[14] = int(len(Ver)-1);
    
    NbrSpl = NbrSpl+1;
    Spl.append([])
    Spl[-1].append(OutVid[12]);
    for j in range(0, NbrLnk):
        Spl[-1].append(NbrVer-NbrLnk+j+1);
    Spl[-1].append(OutVid[14]);
    
    OutLid[16] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[6], OutVid[9], OutVid[12]]);
    OutLid[17] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[4], OutVid[13], OutVid[11]]);
    OutLid[18] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[7], OutVid[15], OutVid[14]]);
    OutLid[19] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([OutVid[5], OutVid[8], OutVid[10]]);
    OutLid[20] = len(Spl)-1;
    
    NbrSpl = NbrSpl+1;
    Spl.append([SavIdx[0][0], SavIdx[0][3]]);
    OutLid[21] = len(Spl)-1;
    
    # --- Outside top
    Loo.append([OutLid[20], OutLid[14], -OutLid[18], -OutLid[5]]);
    Loo.append([OutLid[18], OutLid[15], -OutLid[17], -OutLid[9]]);
    
    NbrSur = len(Loo)-1
    Phy.append([NbrSur-1, NbrSur])
    
    # --- Outside bottom
    Loo.append([OutLid[17], OutLid[16], -OutLid[19], -OutLid[11]]);
    
    NbrSur = len(Loo)-1
    Phy.append([NbrSur])
    
    # --- Inlet
    Loo.append([SavSpl[0][0], SavSpl[0][1], SavSpl[0][2], -OutLid[21]]);
    
    NbrSur = len(Loo)-1
    Phy.append([NbrSur])
    
    # --------------------
    # --- Add box 
    # --------------------
    
    Ver.append([Box[0][0], Box[1][0], Box[2][0], 3])
    OutVid[15] = int(len(Ver)-1);
    
    Ver.append([Box[0][0], Box[1][1], Box[2][0], 3])
    OutVid[16] = int(len(Ver)-1);
    
    Ver.append([Box[0][0], Box[1][1], Box[2][1], 3])
    OutVid[17] = int(len(Ver)-1);
    
    Ver.append([Box[0][0], Box[1][0], Box[2][1], 3])
    OutVid[18] = int(len(Ver)-1);
    
    Ver.append([Box[0][1], Box[1][0], Box[2][0], 3])
    OutVid[19] = int(len(Ver)-1);
    
    Ver.append([Box[0][1], Box[1][1], Box[2][0], 3])
    OutVid[20] = int(len(Ver)-1);
    
    Ver.append([Box[0][1], Box[1][1], Box[2][1], 3])
    OutVid[21] = int(len(Ver)-1);
    
    Ver.append([Box[0][1], Box[1][0], Box[2][1], 3])
    OutVid[22] = int(len(Ver)-1);
    
    # --- Lines
    
    # - x=0
    Spl.append([OutVid[15], OutVid[16]]);
    OutLid[22] = len(Spl)-1;    
    
    Spl.append([OutVid[16], OutVid[17]]);
    OutLid[23] = len(Spl)-1;
    
    Spl.append([OutVid[17], OutVid[18]]);
    OutLid[24] = len(Spl)-1;
    
    Spl.append([OutVid[18], OutVid[10]]);
    OutLid[25] = len(Spl)-1;
    
    Spl.append([OutVid[14], OutVid[15]]);
    OutLid[26] = len(Spl)-1;
    
    # - x=end    
    Spl.append([OutVid[19], OutVid[20]]);
    OutLid[27] = len(Spl)-1;
    
    Spl.append([OutVid[20], OutVid[21]]);
    OutLid[28] = len(Spl)-1;
    
    Spl.append([OutVid[21], OutVid[22]]);
    OutLid[29] = len(Spl)-1;
    
    Spl.append([OutVid[22], OutVid[19]]);
    OutLid[30] = len(Spl)-1;
    
    # - 
    
    Spl.append([OutVid[20], OutVid[16]]);
    OutLid[31] = len(Spl)-1;
    
    Spl.append([OutVid[21], OutVid[17]]);
    OutLid[32] = len(Spl)-1;
    
    Spl.append([OutVid[22], OutVid[18]]);
    OutLid[33] = len(Spl)-1;
    
    Spl.append([OutVid[19], OutVid[15]]);
    OutLid[34] = len(Spl)-1;
    
    # --- Surfaces
    
    Loo.append([OutLid[22], OutLid[23], OutLid[24], OutLid[25], OutLid[14], OutLid[15], OutLid[16], OutLid[26]]);
    Phy.append([len(Loo)-1]);
    
    Loo.append([OutLid[31], OutLid[23], -OutLid[32], -OutLid[28]]);
    Phy.append([len(Loo)-1]);
    
    Loo.append([OutLid[32], OutLid[24], -OutLid[33], -OutLid[29]]);
    Phy.append([len(Loo)-1]);
    
    Spl.append([SavIdx[NbrSec-1][0], SavIdx[NbrSec-1][3]]);
    OutLid[35] = len(Spl)-1;
    
    # -- Symmetry plane
    
    Loo.append([OutLid[33], OutLid[25], -OutLid[20], OutLid[8], OutLid[35], SavSpl[NbrSec][4], -OutLid[3], -OutLid[12], OutLid[19], OutLid[26],  -OutLid[34], -OutLid[30]]);
    Phy.append([len(Loo)-1]);
    
    Loo.append([]);
    for i in range(1, NbrSec):
        Loo[-1].append(SavSpl[i][3]);
    Loo[-1].append(OutLid[35]);
    for i in range(1, NbrSec):
        Loo[-1].append(-SavSpl[NbrSec-i][6]);
    Loo[-1].append(-OutLid[21]);
    Phy.append([len(Loo)-1]);
    
    Loo.append([OutLid[27], OutLid[31], -OutLid[22], -OutLid[34]]);
    Phy.append([len(Loo)-1]);
    
    Loo.append([OutLid[27], OutLid[28], OutLid[29], OutLid[30]]);
    Phy.append([len(Loo)-1]);
    
    ## -- Exit plane (NOT physical)
    #
    #Loo.append([SavSpl[NbrSec-1][0],SavSpl[NbrSec-1][1],SavSpl[NbrSec-1][2],-OutLid[35]]);
    #SavExi = len(Loo)-1;
    
    # ---------------------------------
    # --- Define background sizefields
    # ---------------------------------
        
    # size_in size_out xmax xmin ymax ymin zmax zmin
    hei = 0.25*(Ver[int(OutVid[10])][2]-Ver[int(OutVid[14])][2]);
    Bak.append([hl1, 1000, 100, -100, Ver[int(OutVid[11])][1]+hei, -1, Ver[int(OutVid[10])][2]+hei, Ver[int(OutVid[14])][2]-hei]);
    
    # size_in size_out xmax xmin ymax ymin zmax zmin
    hei = 0.25*(Ver[int(OutVid[5])][2]-Ver[int(OutVid[7])][2]);
    Bak.append([hl0, 1000, 100, Ver[int(OutVid[5])][0], Ver[int(OutVid[4])][1]+hei, -1, Ver[int(OutVid[5])][2]+hei, Ver[int(OutVid[7])][2]-hei]);
    
    WriteGeo(FilNam, Ver, Spl, Lin, Loo, Phy, Siz, Bak, 3)
    
    # ---------------------------------
    # --- Write nozzle exit plane (for thrust computation)
    # ---------------------------------
    
    #savLoo = Loo[SavExi];
    Loo = [0];
    #Loo.append(savLoo);
    Loo.append([SavSpl[NbrSec-1][0],SavSpl[NbrSec-1][1],SavSpl[NbrSec-1][2],-OutLid[35]]);
    
    Phy = [0];
    Phy.append([1]);   # Physical groups
    
    ExiNam = "%s_exit.geo" %  os.path.splitext(os.path.basename(FilNam))[0]
    WriteGeo(ExiNam, Ver, Spl, Lin, Loo, Phy, Siz, Bak, 3)
    
def HF_GenerateMesh(nozzle):
    
    params = np.zeros(100)
    params[0] = -0.0739843; # z crd of the cut at the throat
    params[1] = -0.458624;  # z crd of the flat exit (bottom of the shovel)
    params[2] = -0.15;      # z coordinate of the top of the shovel
    params[3] = 2.0;        # x_throat 
    params[4] = 4.0;        # x_exit 
    #params[5] = 10;             # xmax
    params[5] = 7;              # xmax
    params[6] = 8;              # ymax
    params[7] = 10;             # zmax
    params[8] = 0.04;       # thickness of shovel + exit
    
    # --- Size parameters : 
    #         - inside nozzle
    #         - surface aircraft
    #         - farfield
    #         - (box) plume region
    #         - (box) aircraft vicinity
    #sizes = [0.02, 0.1, 1,  0.03, 0.08 ];
    #sizes = [0.035, 0.1, 1,  0.06, 0.1 ];
    sizes = [0.035, 0.2, 1,  0.09, 0.15 ];
    
    
    # --- Parse centerline bspline
    
    inpBsp = np.loadtxt('inputbspline.dat');
    
    k = 3;
    Nbc = len(inpBsp)/2;
    [t,c] = np.split(inpBsp,2)
    tck = [t,c,3];
    
    def fz(x):
        return splev(x, tck)
    
    # --- Get cross sections definition
    
    inpData = np.loadtxt('inputsection.dat');
    NbrSec = len(inpData);
    
    x  = inpData[:,0];
    r1 = inpData[:,1];
    r2 = inpData[:,2];
    
    tck1 = splrep(x, r1, xb=x[0], xe=x[-1], k=3, s=0);
    def fr1(x):
        return splev(x, tck1);
        
    tck2 = splrep(x, r2, xb=x[0], xe=x[-1], k=3, s=0);
    def fr2(x):
        return splev(x, tck2);
    
    cadNam     = "nozzle.geo";
    cadExitNam = "nozzle_exit.geo";
    
    HF_DefineCAD(cadNam,x,fr1,fr2, fz, sizes, params);
    
    gmsh_executable = 'gmsh';
    try :
        cmd = [gmsh_executable, '-3', cadNam, '-o', nozzle.cfd.mesh_name];
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
    except:
        raise;
        
        
    try :
        cmd = [gmsh_executable, '-2', cadExitNam, '-o', nozzle.cfd.exit_mesh_name];
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=None)
    except:
        raise;


def LoadBSplineParameters(FilNam):
    
    try:
        hdl = np.loadtxt(FilNam);
    except:
        sys.stderr.write('   ## ERROR: Unable to load baseline BSpline parameters.\n');
        sys.stderr.write('      %s NOT FOUND\n');
        sys.exit(1);
    
    NKnots_center = int(hdl[0]);
    NCoefs_center = int(hdl[1]);
    NKnots_r1     = int(hdl[2]);
    NCoefs_r1     = int(hdl[3]);
    NKnots_r2     = int(hdl[4]);
    NCoefs_r2     = int(hdl[5]);
    
    Knots_center = [];
    Coefs_center = [];
    Knots_r1     = [];
    Coefs_r1     = [];
    Knots_r2     = [];
    Coefs_r2     = [];
    
    cur = 6;
    for i in range(cur, cur+NKnots_center):
        Knots_center.append(hdl[i]);
    
    cur = i+1;
    
    for i in range(cur, cur+NCoefs_center):
        Coefs_center.append(hdl[i]);
        
    cur = i+1;
    
    for i in range(cur, cur+NKnots_r1):
        Knots_r1.append(hdl[i]);
        
    cur = i+1;
    
    for i in range(cur, cur+NCoefs_r1):
        Coefs_r1.append(hdl[i]);
        
    cur = i+1;
    
    for i in range(cur, cur+NKnots_r2):
        Knots_r2.append(hdl[i]);
        
    cur = i+1;
    
    for i in range(cur, cur+NCoefs_r2):
        Coefs_r2.append(hdl[i]);
        
        
    return Knots_center, Coefs_center, Knots_r1, Coefs_r1, Knots_r2, Coefs_r2;


def SwitchTriangleIdx(tri):
    
    vid = tri;
    
    if vid[1] < vid[0] :
        if vid[2] < vid[1] :
            tri[0] = vid[2];
            tri[1] = vid[0];
            tri[2] = vid[1];
        else:
            tri[0] = vid[1];
            tri[1] = vid[2];
            tri[2] = vid[0];
            
    elif vid[2] < vid[0] :
        tri[0] = vid[2];
        tri[1] = vid[0];
        tri[2] = vid[1];
        
    return tri;



def GetEdgeInwardNormal(crd0, crd1):
    
    vno = [0.0,0.0]
    
    vno[0] = crd0[1] - crd1[1];
    vno[1] = crd1[0] - crd0[0];
    
    return vno;

def GetTriangleArea3D(Coord):
    
    crd0 = Coord[0];
    crd1 = Coord[1];
    crd2 = Coord[2];
    
    x = [0.0,0.0,0.0];
    y = [0.0,0.0,0.0];
    z = [0.0,0.0,0.0];
    
    x[0] = 0.;
    y[0] = 0.;
    z[0] = 0.;
    
    x[1] = crd1[0] - crd0[0];
    y[1] = crd1[1] - crd0[1];
    z[1] = crd1[2] - crd0[2];
    
    x[2] = crd2[0] - crd0[0];
    y[2] = crd2[1] - crd0[1];
    z[2] = crd2[2] - crd0[2];
    
    #--- P0P1/\P0P2
    vnx =   y[1]*z[2] - y[2]*z[1];
    vny = - x[1]*z[2] + x[2]*z[1];
    vnz =   x[1]*y[2] - x[2]*y[1];
    
    #--- Area of the face : 2*Air = ||u1 vect u2||
    Area = 0.5*math.sqrt(math.fabs(vnx*vnx + vny*vny + vnz*vnz));
    
    return Area;

def GetTriangleArea2D(Coord):
    
    nDim = 3;
    
    a = [0.0,0.0,0.0];
    b = [0.0,0.0,0.0];
    
    vno = [[0.0,0.0,0.0],[0.0,0.0,0.0]]
    
    Coord_0 = Coord[0];
    Coord_1 = Coord[1];
    Coord_2 = Coord[2];
    
    vno[0] = GetEdgeInwardNormal(Coord_1, Coord_2);
    vno[1] = GetEdgeInwardNormal(Coord_2, Coord_0);
    
    print vno
    
    Area = math.fabs(vno[0][0]*vno[1][1]-vno[0][1]*vno[1][0]);
    
    return Area;
  


def Get2DMeshArea(MshNam):
    
    #--- Load mesh
    
    SolNam = ""
    
    Ver = [];
    Tri = [];
    Tet = [];
    Edg = [];
    Sol = [];
    
    _meshutils_module.py_ReadMesh(MshNam, SolNam, Ver, Tri, Tet, Edg, Sol);
    
    NbrTri = len(Tri)/4;
    Tri = np.reshape(Tri,(NbrTri, 4)).astype(int);
    
    NbrVer = len(Ver)/3;
    Ver = np.reshape(Ver,(NbrVer, 3));
    
    #--- Compute area
    
    NbrTri = len(Tri);
    
    Coord = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    
    AreaTot = 0.0;
    
    for iTri in range(0,NbrTri):
        
        for j in range(3):
            iVer = Tri[iTri][j]-1;
            Coord[j]  = Ver[iVer]; 

        area = GetTriangleArea3D(Coord);
        iVer = Tri[iTri][0]-1;
        
        AreaTot += area;
        
    return AreaTot;
 
  
def HF_GenerateMesh_Deform(nozzle):
    
    from .. import _meshutils_module
    
    RefUp   = [9];
    RefDown = [10];
    
    Knots_center =  nozzle.wall.centerline.knots;
    Coefs_center =  nozzle.wall.centerline.coefs;
    Knots_r1     =  nozzle.wall.majoraxis.knots;
    Coefs_r1     =  nozzle.wall.majoraxis.coefs;
    Knots_r2     =  nozzle.wall.minoraxis.knots;
    Coefs_r2     =  nozzle.wall.minoraxis.coefs;
    
    
    #pathsrc = "%s/baseline_meshes/" % (os.path.dirname(os.path.abspath(__file__)));
    #basNamGMF   = "%sbaseline_%s_%s.meshb" % (pathsrc, nozzle.method.lower(), nozzle.cfd.mesh_size.lower());
    #basNamSU2    = "baseline_%s_%s.su2" % (nozzle.method.lower(), nozzle.cfd.mesh_size.lower());
    #
    #_meshutils_module.py_ProjectNozzleWall3D(basNamGMF, RefUp, RefDown,
    #Knots_center, Coefs_center,
    #Knots_r1, Coefs_r1  ,
    #Knots_r2, Coefs_r2  ,
    # "mesh_motion.dat");
    
    
    pathsrc = "%s/baseline_meshes/" % (os.path.dirname(os.path.abspath(__file__)));
    basNamGMF   = "%sbaseline_1_%s_%s.meshb" % (pathsrc, nozzle.method.lower(), nozzle.cfd.mesh_size.lower());
    basNamDV = "%sbaseline_1_DV.dat" % pathsrc;
    basNamSU2    = "baseline_%s_%s.su2" % (nozzle.method.lower(), nozzle.cfd.mesh_size.lower());
    
    Knots_center_bas, Coefs_center_bas, \
    Knots_r1_bas, Coefs_r1_bas  , \
    Knots_r2_bas, Coefs_r2_bas  = LoadBSplineParameters(basNamDV);
        
    _meshutils_module.py_ProjectNozzleWall3D_DV(basNamGMF, RefUp, RefDown,
    Knots_center_bas, Coefs_center_bas,
    Knots_r1_bas, Coefs_r1_bas  ,
    Knots_r2_bas, Coefs_r2_bas  ,
    Knots_center, Coefs_center,
    Knots_r1, Coefs_r1  ,
    Knots_r2, Coefs_r2  ,
     "mesh_motion.dat");
     
    _meshutils_module.py_ConvertGMFToSU2(basNamGMF,"",basNamSU2);
    
    
    
    # --- Call SU2_DEF to deform baseline mesh
    # --- Setup config file
    
    from .. import SU2
    
    config = SU2.io.Config();
    
    # -- Note : the values don't matter. All markers must be defined otherwise SU2 exits.
    config.MARKER_EULER  = '( 7, 8, 9, 10, 12, 13, 14, 15, 16)'
    config.MARKER_INLET  = '(11, 601, 275000, 1.0, 0.0, 0.0 )'
    config.MARKER_FAR    = '( 1, 2, 3, 5, 6)'
    config.MARKER_SYM    = '( 4 )'
    
    config.MESH_FILENAME= basNamSU2;
    config.DV_KIND= "SURFACE_FILE"
    config.DV_MARKER= "( 9, 10 )"
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
    config.MARKER_MOVING         = "( 9, 10 )"
    
    config.MESH_OUT_FILENAME      = nozzle.cfd.mesh_name;
    
    config.NUMBER_PART = nozzle.partitions;
    
    # --- Run SU2
    
    info = SU2.run.DEF(config)


