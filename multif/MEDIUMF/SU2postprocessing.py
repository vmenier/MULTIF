# -*- coding: utf-8 -*-

import os, sys
from .. import SU2
from runSU2 import CheckSU2Convergence

import numpy as np
from scipy.interpolate import interp1d 

from .. import _meshutils_module
from multif import _mshint_module

def PostProcess ( nozzle, output ):
    
    # --- Check residual convergence
    
    from multif import MEDIUMF
    history, FinRes, ResDif = MEDIUMF.checkResidual();
    
    sys.stdout.write("Fluid calculation convergence: Final density residual = %.2lf, Residual reduction = %.2lf\n" % (FinRes, ResDif));
    
    # --- Interpolate and extract all the solution fields along the nozzle exit
    # SolExtract : Solution (x, y, sol1, sol2, etc.)
    # Size : [NbrVer, SolSiz]
    # Header : Name of the solution fields (Conservative_1, Mach, etc.)    
    #SolExtract, Size, Header  = ExtractSolutionAtExit(nozzle);

    # Extract solution at wall
    SolExtract_w, Size_w, idHeader_w  = ExtractSolutionAtWall(nozzle);
    iPres = idHeader_w['Pressure'];
    iTemp = idHeader_w['Temperature'];
    Pres = SolExtract_w[:,iPres];
    Temp = SolExtract_w[:,iTemp];
        
    # --- Assign responses
    if 'THRUST' in nozzle.responses:
                
        if nozzle.method == "EULER":
            SolExtract, Size, Header  = ExtractSolutionAtExit(nozzle);
            #SolExtract, Size, Header  = ExtractExitRANS (nozzle);
            nozzle.responses['THRUST'] = ComputeThrust ( nozzle, SolExtract, Size, Header );
            #nozzle.responses['THRUST'] = Get_Thrust_File(nozzle);
        else :
            
            SolExtract, Size, Header  = ExtractExitRANS ("exit.mesh", "nozzle.su2", "nozzle.dat");
            nozzle.responses['THRUST'] = ComputeThrust ( nozzle, SolExtract, Size, Header );
            #nozzle.responses['THRUST'] = ComputeThrust_RANS (nozzle);
        
    if 'WALL_PRES_AVG' in nozzle.responses or 'WALL_TEMP_AVG' in nozzle.responses:
        
        AreaTot, PresAvg, TempAvg = MF_Integrate_Sol_Wall(nozzle);
        
        if 'WALL_PRES_AVG' in nozzle.responses:
            nozzle.responses['WALL_PRES_AVG'] = PresAvg;
        
        if 'WALL_TEMP_AVG' in nozzle.responses:
            nozzle.responses['WALL_TEMP_AVG'] = TempAvg;
            
    if 'WALL_TEMPERATURE' in nozzle.responses:
        sys.stderr.write(' ## ERROR : WALL_TEMPERATURE not currently available from SU2\n\n');
        sys.exit(1);
        #nozzle.responses['WALL_TEMPERATURE'] = 0;
        
    if 'WALL_PRESSURE' in nozzle.responses:

        func = interp1d(SolExtract_w[:,0],  Pres, kind='linear');        
        nozzle.responses['WALL_PRESSURE'] = np.squeeze(func(nozzle.outputLocations['WALL_PRESSURE']));

        # --- CHECK INTERPOLATION :
        #import matplotlib.pyplot as plt
        #plt.plot(SolExtract[:,0], Pres, "-", label="EXTRACT")
        #plt.plot(nozzle.outputLocations['WALL_PRESSURE'], nozzle.responses['WALL_PRESSURE'], "-", label="OUTPUT")
        #plt.legend()
        #plt.show();
        #sys.exit(1);
        
    if 'PRESSURE' in nozzle.responses:
            
        x = nozzle.outputLocations['PRESSURE'][:,0];
        y = nozzle.outputLocations['PRESSURE'][:,1];
        
        nozzle.responses['PRESSURE'] = np.squeeze(ExtractSolutionAtXY (x, y, ["Pressure"]));
        
    if 'VELOCITY' in nozzle.responses:

        x = nozzle.outputLocations['VELOCITY'][:,0];
        y = nozzle.outputLocations['VELOCITY'][:,1];
        #cons = ExtractSolutionAtXY (x, y, ["Conservative_1","Conservative_2","Conservative_3"]);
        cons = ExtractSolutionAtXY (x, y, ["Density","X-Momentum","Y-Momentum"]);
        
        
        nozzle.responses['VELOCITY'] = [[],[],[]]
        for i in range(len(cons)):
            nozzle.responses['VELOCITY'][0].append(cons[i][1]/cons[i][0]); 
            nozzle.responses['VELOCITY'][1].append(cons[i][2]/cons[i][0]); 
            nozzle.responses['VELOCITY'][2].append(0.0); 

    if 'SU2_RESIDUAL' in nozzle.responses:

        nozzle.responses['SU2_RESIDUAL'] = FinRes

    if output == 'verbose':
        sys.stdout.write('SU2 responses obtained\n');
    
    return 0;  
        

def ExtractSolutionAtExit ( nozzle ):

    mesh_name    = nozzle.cfd.mesh_name;

    restart_name = nozzle.cfd.restart_name;
    
    #pyResult = [];
    #pyInfo   = [];
    #pyHeader = [];
    #
    ##pyBox = [nozzle.length,nozzle.length,0,nozzle.height+1e-20];
    #pyBox = [nozzle.cfd.x_thrust,nozzle.cfd.x_thrust,-1e-20,nozzle.cfd.y_thrust+1e-20];
    #
    #_meshutils_module.py_ExtractAlongLine (mesh_name, restart_name, pyBox, pyResult, pyInfo, pyHeader);
    
    # --- Extract CFD solution at the inner wall    
    
    pyResult = [];
    pyInfo   = [];
    pyHeader = [];
    
    pyRef = [9];
    
    _meshutils_module.py_ExtractAtRef (mesh_name, restart_name, pyRef, pyResult, pyInfo, pyHeader);
    
    NbrRes = pyInfo[0];
    ResSiz = pyInfo[1];
    
    Result = np.asarray(pyResult);
    
    OutResult = np.reshape(Result,(NbrRes, ResSiz));
    Out_sort = OutResult[OutResult[:,1].argsort()]
    
    idHeader = dict();
    for iFld in range(0,len(pyHeader)):
        idHeader[pyHeader[iFld]] = iFld+3;

    return OutResult, pyInfo, idHeader;
    
    
    #return Out_sort, pyInfo, pyHeader;


def Get_Thrust_File(nozzle):
    
    thrust_filename = "thrust_nodef.dat"
    
    if not os.path.isfile(thrust_filename) :
        sys.stderr.write("  ## ERROR Get thrust : file %s not found.\n \
        Are you using the right version of SU2?\n" % thrust_filename);
        return -1;
        
    thrust = np.loadtxt(thrust_filename);
    return float(thrust);



def ComputeThrust_2 ( options, SolExtract, Size, Header )    :

    # T = 2PI * Int_{0}^{R} (rho U ( U - U0) + P - Po ) r dr

    NbrVer = Size[0];
    SolSiz = Size[1];

    print SolSiz
    print Header

    if len(Header) != SolSiz-2 and len(Header) != SolSiz-3:
        sys.stderr.write("  ## ERROR : ComputeThrust : Inconsistent solution header (%d != %d).\n" % (len(Header),SolSiz-2));
        sys.exit(0);

    ## --- Get solution field indices

    # Deprecated keywords : 
    #iMach  = Header['Mach'];
    #iTem   = Header['Temperature'];
    #iCons1 = Header['Conservative_1'];
    #iCons2 = Header['Conservative_2'];
    #iCons3 = Header['Conservative_3'];
    #iCons4 = Header['Conservative_4'];
    #iPres  = Header['Pressure'];
    
    # New keywords (> SU2 Raven 5.0)
    iMach  = Header['Mach'];
    iTem   = Header['Temperature'];
    iCons1 = Header['Density'];
    iCons2 = Header['X-Momentum'];
    iCons3 = Header['Y-Momentum'];
    iCons4 = Header['Energy'];
    iPres  = Header['Pressure'];
    
    # --- Compute thrust    

    Thrust = 0;

    #freestream.P = atm.P; % Pa, atmospheric pressure
    #freestream.T = atm.T; % K, atmospheric temperature
    #freestream.M = mach;
    #freestream.U = freestream.M*sqrt(fluid.gam*fluid.R*freestream.T);

    #P0  = nozzle.environment.P;
    #M0  = nozzle.mission.mach;
    #Gam = nozzle.fluid.gam;
    #Rs  = nozzle.fluid.R;
    #T0  = nozzle.environment.T;
    #U0  = M0*np.sqrt(Gam*Rs*T0);
    
    P0  = options["P0"];
    M0  = options["M0"];
    Gam = options["Gam"];
    Rs  = options["Rs"];
    T0  = options["T0"];
    U0  = options["U0"];
    


    #for iVer in range(1, NbrVer) :
    #    

    #    y    = float(SolExtract[iVer][1]);

    #    
    #    #if y > nozzle.height-1e-6:
    #    #    print "REMOVE POINT %d" % iVer
    #    #    continue;
    #    
    #    rho  = 0.5*(SolExtract[iVer][2+iCons1] +  SolExtract[iVer-1][2+iCons1]);
    #    rhoU = 0.5*(SolExtract[iVer][2+iCons2] +  SolExtract[iVer-1][2+iCons2]);
    #    Pres = 0.5*(SolExtract[iVer][2+iPres]  +  SolExtract[iVer-1][2+iPres] );
    #    Mach = 0.5*(SolExtract[iVer][2+iMach]  +  SolExtract[iVer-1][2+iMach] );
    #    Temp = 0.5*(SolExtract[iVer][2+iTem]   +  SolExtract[iVer-1][2+iTem]  );
    #    
    #    
    #    U = rhoU/rho;
    #    
    #    dy = y - SolExtract[iVer-1][1];
    #    
    #    #print "%lf %lf %lf %lf %lf %lf" % (y, rho, rhoU, Pres, Mach, Temp);
    #            
    #    Thrust = Thrust + dy*(rhoU*(U-U0)+Pres-P0);
    #
    #print "THRUST = %lf" % Thrust

    #y = SolExtract[iVer][];

    NbrVer = len(SolExtract);

    y   = np.zeros(NbrVer);
    sol = np.zeros([NbrVer,5]);

    for i in range(0,NbrVer):
        y[i] = float(SolExtract[i][1]);

        sol[i][0] = float(SolExtract[i][iCons1]);
        sol[i][1] = float(SolExtract[i][iCons2]);
        sol[i][2] = float(SolExtract[i][iPres]);


    fsol = [];
    for j in range(0,3):
        fsol.append(interp1d(y,sol[:,j], kind='linear'));

    nbv = 4000;
    ynew = np.linspace(0,y[-1],nbv);

    tabrho  = fsol[0](ynew);
    tabrhoU = fsol[1](ynew);
    tabPres = fsol[2](ynew);

    fil = open('exit.dat','w')

    for i in range(1, nbv) :

        rho  = 0.5*(tabrho[i-1]+tabrho[i]);
        rhoU = 0.5*(tabrhoU[i-1]+tabrhoU[i]);
        Pres = 0.5*(tabPres[i-1]+tabPres[i]);

        U = rhoU/rho;

        dy = ynew[i]-ynew[i-1];

        fil.write("%lf %lf %lf %lf\n" % (ynew[i], rho, rhoU, Pres));

        Thrust = Thrust + dy*(rhoU*(U-U0)+Pres-P0);    

    fil.close();

    return Thrust;


def MF_Integrate_Sol_Wall(nozzle):
    
    mesh_name    = nozzle.cfd.mesh_name;
    restart_name = nozzle.cfd.restart_name;

    pyResult = [];
    pyInfo   = [];
    pyHeader = [];
    
    pyRef = [1];
    
    _meshutils_module.py_ExtractAtRef (mesh_name, restart_name, pyRef, pyResult, pyInfo, pyHeader);
    
    NbrRes = pyInfo[0];
    ResSiz = pyInfo[1];
    
    Result = np.asarray(pyResult);
    
    OutResult = np.reshape(Result,(NbrRes, ResSiz));
    Out_sort = OutResult[OutResult[:,0].argsort()]
    
    idHeader = dict();
    for iFld in range(0,len(pyHeader)):
        idHeader[pyHeader[iFld]] = iFld+3;
    
    iTemp = idHeader['Temperature'];
    iPres = idHeader['Pressure']; 
    
    pres_avg = 0.0;
    temp_avg = 0.0;
    
    LenTot = 0.0;
    for i in range(len(Out_sort)-1):
        
        nrm = 0.0;
        for j in range(2):
            nrm += (Out_sort[i+1][j]-Out_sort[i][j])*(Out_sort[i+1][j]-Out_sort[i][j]);
        nrm = np.sqrt(nrm);
        
        LenTot += nrm;
        
        pres = 0.5*(Out_sort[i][iPres]+Out_sort[i+1][iPres]);
        temp = 0.5*(Out_sort[i][iTemp]+Out_sort[i+1][iTemp]);
        
        pres_avg += nrm*pres;
        temp_avg += nrm*temp;
    
    pres_avg /= LenTot;
    temp_avg /= LenTot;
    
    
    return LenTot, pres_avg, temp_avg;

    

def ComputeThrust ( nozzle, SolExtract, Size, Header )    :
    
    # T = 2PI * Int_{0}^{R} (rho U ( U - U0) + P - Po ) r dr
    
    NbrVer = Size[0];
    SolSiz = Size[1];
    
    print SolSiz
    print Header
    
    if len(Header) != SolSiz-2 and len(Header) != SolSiz-3:
        sys.stderr.write("  ## ERROR : ComputeThrust : Inconsistent solution header (%d != %d).\n" % (len(Header),SolSiz-2));
        sys.exit(0);
    
    ## --- Get solution field indices
    
    # Deprecated keywords:
    #iMach  = Header['Mach'];
    #iTem   = Header['Temperature'];
    #iCons1 = Header['Conservative_1'];
    #iCons2 = Header['Conservative_2'];
    #iCons3 = Header['Conservative_3'];
    #iCons4 = Header['Conservative_4'];
    #iPres  = Header['Pressure'];
    
    # New keywords (> SU2 Raven 5.0)
    iMach  = Header['Mach'];
    iTem   = Header['Temperature'];
    iCons1 = Header['Density'];
    iCons2 = Header['X-Momentum'];
    iCons3 = Header['Y-Momentum'];
    iCons4 = Header['Energy'];
    iPres  = Header['Pressure'];
    
    
    # --- Compute thrust    
        
    Thrust = 0;
    
    #freestream.P = atm.P; % Pa, atmospheric pressure
    #freestream.T = atm.T; % K, atmospheric temperature
    #freestream.M = mach;
    #freestream.U = freestream.M*sqrt(fluid.gam*fluid.R*freestream.T);
    
    P0  = nozzle.environment.P;
    M0  = nozzle.mission.mach;
    Gam = nozzle.fluid.gam;
    Rs  = nozzle.fluid.R;
    T0  = nozzle.environment.T;
    U0  = M0*np.sqrt(Gam*Rs*T0);
        
    NbrVer = len(SolExtract);
    
    y   = np.zeros(NbrVer);
    sol = np.zeros([NbrVer,5]);
    
    tabrho  = [];
    tabrhoU = [];
    tabPres = [];
    
    for i in range(0,NbrVer):
        y[i] = float(SolExtract[i][1]);
        
        tabrho.append(float(SolExtract[i][iCons1]));
        tabrhoU.append(float(SolExtract[i][iCons2]));
        tabPres.append(float(SolExtract[i][iPres]));
        
    
    
    #fsol = [];
    #for j in range(0,3):
    #    fsol.append(interp1d(y,sol[:,j], kind='linear'));
    #
    #nbv = 4000;
    #ynew = np.linspace(0,y[-1],nbv);
        
    #tabrho  = fsol[0](ynew);
    #tabrhoU = fsol[1](ynew);
    #tabPres = fsol[2](ynew);

    ynew = y;
    nbv = len(ynew);
    
    fil = open('exit.dat','w')
    
    for i in range(1, nbv) :
        
        rho  = 0.5*(tabrho[i-1]+tabrho[i]);
        rhoU = 0.5*(tabrhoU[i-1]+tabrhoU[i]);
        Pres = 0.5*(tabPres[i-1]+tabPres[i]);
                
        U = rhoU/rho;
        
        dy = ynew[i]-ynew[i-1];
        
        fil.write("%lf %lf %lf %lf\n" % (ynew[i], rho, rhoU, Pres));
                    
        Thrust = Thrust + dy*(rhoU*(U-U0)+Pres-P0);    
    
    fil.close();
    
    return Thrust;


def ExtractExitRANS (exit_name, mesh_name, sol_name)    :
    
    exitNam = "exit.mesh";
    
    # --- Interpolate solution
    
    info = [];
    Crd  = [];
    Tri  = [];
    Tet  = [];
    Sol  = [];
    Header = [];
    
    out = _mshint_module.py_Interpolation (exitNam,mesh_name, sol_name,\
        info, Crd, Tri, Tet, Sol, Header);
    
    # --- Extract CFD solution at the inner wall    
    
    pyResult = [];
    pyInfo   = [];
    pyHeader = [];
    
    pyRef = [1];
    
    _meshutils_module.py_ExtractAtRef ("exit.mesh", "exit.solb", pyRef, pyResult, pyInfo, pyHeader);
    
    NbrRes = pyInfo[0];
    ResSiz = pyInfo[1];
        
    Result = np.asarray(pyResult);
    
    OutResult = np.reshape(Result,(NbrRes, ResSiz));
        
    #print "%d results" % NbrRes;
    #
    #for i in range(0,10):
    #    print "%d : (%lf,%lf) : rho = %lf" % (i, OutResult[i][0], OutResult[i][1], OutResult[i][3]);
    
    # --- Get solution field indices
    
    iMach  = -1;
    iTem   = -1;
    iCons1 = -1;
    iCons2 = -1;
    iCons3 = -1;
    iCons4 = -1;
    iPres  = -1;
    
    idHeader = dict();
    for iFld in range(0,len(Header)):
        idHeader[Header[iFld]] = iFld+3;
    
    return OutResult, pyInfo, idHeader;
    
    
def ExtractSolutionAtWall (nozzle):
    
    # --- Extract CFD solution at the inner wall    

    mesh_name    = nozzle.cfd.mesh_name;

    restart_name = nozzle.cfd.restart_name;
    
    pyResult = [];
    pyInfo   = [];
    pyHeader = [];
    
    pyRef = [1];
    
    _meshutils_module.py_ExtractAtRef (mesh_name, restart_name, pyRef, pyResult, pyInfo, pyHeader);
    
    NbrRes = pyInfo[0];
    ResSiz = pyInfo[1];
        
    Result = np.asarray(pyResult);
    
    OutResult = np.reshape(Result,(NbrRes, ResSiz));
        
    #print "%d results" % NbrRes;
    #
    #for i in range(0,10):
    #    print "%d : (%lf,%lf) : rho = %lf" % (i, OutResult[i][0], OutResult[i][1], OutResult[i][2]);
        
    idHeader = dict();
    for iFld in range(0,len(pyHeader)):
        idHeader[pyHeader[iFld]] = iFld+3;

    return OutResult, pyInfo, idHeader;



def WriteGMFMesh2D(MshNam, Ver, Tri):
    
    f = open(MshNam, 'wb');
    
    NbrVer = len(Ver);
    NbrTri = len(Tri);
    
    f.write("MeshVersionFormatted\n2\nDimension\n2\n\n");


    #--- Write vertices
    f.write("Vertices\n%d\n" % NbrVer);
    for i in range(0,NbrVer):
        f.write("%lf %lf %d\n" % (Ver[i][0],Ver[i][1],0));
        
    #--- Write triangles
    f.write("\n\nTriangles\n%d\n" % NbrTri);
    for i in range(0,NbrTri):
        f.write("%d %d %d 1\n" % (Tri[i][0],Tri[i][1], Tri[i][2]));
        
    f.write("\nEnd\n");
    
    if f :
        f.close();


def ExtractSolutionAtXY (x, y, tagField):
    
    Ver = [];
    Tri = [];
    
    # --- Create structured mesh
    
    x, id_x = np.unique(x.round(decimals=4),return_inverse=True);
    y, id_y = np.unique(y.round(decimals=4),return_inverse=True);
    
    Ni = len(x);
    Nj = len(y);
    
    if Ni == 1:
        x = np.append(x, x[0]+0.1);
        Ni=Ni+1;
        
    if Nj == 1:
        y = np.append(y, y[0]+0.1);
        Nj=Nj+1;
        
    NbrVer = Ni*Nj;
    
    for i in range(0,Ni):
        for j in range(0,Nj):
            Ver.append([x[i],y[j],0.0]);
    
    NbrTri = 0;
    
    for i in range(2,Ni+1):
        for j in range(2,Nj+1):
            ind = (i-2)*Nj+j-1;
    
            Tri.append([ind, ind+Nj, ind+Nj+1]);
            Tri.append([ind, ind+1+Nj, ind+1,1]);

    WriteGMFMesh2D('nozzle_extract.mesh', Ver, Tri);
    
    info = [];
    Crd  = [];
    Tri  = [];
    Tet  = [];
    Sol  = [];
    Header = [];
    
    out = _mshint_module.py_Interpolation ("nozzle_extract.mesh", "nozzle.su2", "nozzle.dat",\
        info, Crd, Tri, Tet, Sol, Header);
    
    dim    = info[0];
    NbrVer = info[1]; 
    NbrTri = info[2];
    NbrTet = info[3];
    SolSiz = info[4];
    
    NbrFld = len(tagField);
    iFldTab = [];
        
    for iTab in range(NbrFld):
        
        tag = tagField[iTab];
        iFld = -1        
        
        for i in range(0,len(Header)):
            if ( Header[i] == tag ):
                iFld = i;
                break;
        
        if iFld == -1 :
            sys.stderr.write("  ## ERROR Extraction solution : required field not found in solution file.\n");
            sys.exit(1);
        
        iFldTab.append(iFld);
        
    OutSol = [];
    
    for i in range(len(id_x)):
        
        ii = id_x[i]+1;
        jj = id_y[i]+1;
        
        iVer = (ii-1)*Nj+jj;
        
        idx  = (iVer-1)*dim;
        idxs = (iVer-1)*SolSiz;
        
        
        solTab = [];
        
        for iTab in range(NbrFld):
            iFld = iFldTab[iTab];
            solTab.append(Sol[idxs+iFld+1]);
        
        OutSol.append(solTab);
        
    return OutSol;

    
def VerificationPostPro(nozzle):
    
    # --- Load validation data
    
    pathsrc = "%s/verification_data/" % (os.path.dirname(os.path.abspath(__file__)));
    
    #--- Extract boundary vertices from mesh
    
    wall_name    = "%swall_%s_%s.dat" % (pathsrc, nozzle.method, nozzle.meshsize);
    
    wall_hdl = np.loadtxt(wall_name);
    
    # --- Extract CFD solution at the inner wall    
    
    mesh_name    = nozzle.mesh_name;

    restart_name = nozzle.restart_name;
    
    pyResult = [];
    pyInfo   = [];
    pyHeader = [];
    
    pyRef = [1];
    
    _meshutils_module.py_ExtractAtRef (mesh_name, restart_name, pyRef, pyResult, pyInfo, pyHeader);
    
    NbrRes = pyInfo[0];
    ResSiz = pyInfo[1];
        
    Result = np.asarray(pyResult);
    
    OutResult = np.reshape(Result,(NbrRes, ResSiz));
    
    # --- Get solution field indices
        
    idHeader = dict();
    for iFld in range(0,len(pyHeader)):
        idHeader[pyHeader[iFld]] = iFld+3;
    
    Pres = OutResult[:,idHeader['Pressure']];
    Temp = OutResult[:,idHeader['Temperature']];
    
    fil = open("extraction_wall_valid.dat",'w');
    
    for i in range(len(wall_hdl)):
        fil.write("%lf %lf %lf %lf\n" % (wall_hdl[i,0], wall_hdl[i,1], wall_hdl[i,2], wall_hdl[i,3]));
    fil.close();
        
    fil = open("extraction_wall.dat",'w');
    
    for i in range(len(Pres)):
        fil.write("%lf %lf %lf %lf\n" % (OutResult[i,0], OutResult[i,1], Pres[i], Temp[i]));
    fil.close();
    
    
    
    
    
