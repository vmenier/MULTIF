from multif import _mshint_module
from multif import _meshutils_module
import numpy as np
import sys, math
import os
import multif

from .. import SU2

def PostProcess ( nozzle, output ):
    
    # --- Check residual convergence
    
    from multif import MEDIUMF
    history, FinRes, ResDif = MEDIUMF.checkResidual();
    sys.stdout.write("Fluid calculation convergence: Final density residual = %.2lf, Residual reduction = %.2lf\n" % (FinRes, ResDif));
    
    #IniRes, FinRes = CheckConvergence(nozzle);
    #ResDif = FinRes - IniRes;
    #sys.stdout.write("Initial res = %le, Final res = %lf, Diff = %lf\n" % (IniRes, FinRes, ResDif));
    
    # --- Interpolate and extract all the solution fields along the nozzle exit
    # SolExtract : Solution (x, y, sol1, sol2, etc.)
    # Size : [NbrVer, SolSiz]
    # Header : Name of the solution fields (Conservative_1, Mach, etc.)    
    #SolExtract, Size, Header  = ExtractSolutionAtExit(nozzle);

    # Extract solution at wall
# XXX    SolExtract_w, Size_w, idHeader_w  = ExtractSolutionAtWall(nozzle);
# XXX    iPres = idHeader_w['Pressure'];
# XXX    iTemp = idHeader_w['Temperature'];
# XXX    Pres = SolExtract_w[:,iPres];
# XXX    Temp = SolExtract_w[:,iTemp];
        
    # --- Assign responses
    if 'THRUST' in nozzle.responses:
        # XXX Ensure thrust calculation is correct     
        
        nozzle.responses['THRUST'] = HF_Compute_Thrust_Wrap (nozzle);
        
#        SolExtract, Size, Header  = ExtractSolutionAtExit(nozzle);
#        nozzle.responses['THRUST'] = ComputeThrust ( nozzle, SolExtract, Size, Header );
        #nozzle.responses['THRUST'] = Get_Thrust_File(nozzle);
        
    if 'WALL_TEMPERATURE' in nozzle.responses:
        sys.stderr.write(' ## ERROR : WALL_TEMPERATURE not currently available from SU2\n\n');
        sys.exit(1);
        #nozzle.responses['WALL_TEMPERATURE'] = 0;
    
    if 'WALL_PRES_AVG' in nozzle.responses or 'WALL_TEMP_AVG' in nozzle.responses:
        
        AreaTot, PresAvg, TempAvg = HF_Integrate_Sol_Wall(nozzle);
        
        if 'WALL_PRES_AVG' in nozzle.responses:
            nozzle.responses['WALL_PRES_AVG'] = PresAvg;
        
        if 'WALL_TEMP_AVG' in nozzle.responses:
            nozzle.responses['WALL_TEMP_AVG'] = TempAvg;
    
    if 'WALL_PRESSURE' in nozzle.responses:

# XXX        func = interp1d(SolExtract_w[:,0],  Pres, kind='linear');        
# XXX        nozzle.responses['WALL_PRESSURE'] = np.squeeze(func(nozzle.outputLocations['WALL_PRESSURE']));
        nozzle.responses['WALL_PRESSURE'] = [-1]*len(nozzle.outputLocations['WALL_PRESSURE']) # temporary
  
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
        
# XXX        nozzle.responses['PRESSURE'] = np.squeeze(ExtractSolutionAtXY (x, y, ["Pressure"]));
        nozzle.responses['PRESSURE'] = [-1]*len(x) # temporary
        
    if 'VELOCITY' in nozzle.responses:

        x = nozzle.outputLocations['VELOCITY'][:,0];
        y = nozzle.outputLocations['VELOCITY'][:,1];
# XXX        cons = ExtractSolutionAtXY (x, y, ["Conservative_1","Conservative_2","Conservative_3"]);
        # cons = [[1.,0.,0.]] # temporary
                    
        nozzle.responses['VELOCITY'] = [[],[],[]]
        # for i in range(len(cons)):
        #     nozzle.responses['VELOCITY'][0].append(cons[i][1]/cons[i][0]); 
        #     nozzle.responses['VELOCITY'][1].append(cons[i][2]/cons[i][0]); 
        #     nozzle.responses['VELOCITY'][2].append(0.0); 

        # Return -1 showing lack of data
        for i in range(len(x)):
            nozzle.responses['VELOCITY'][0].append(-1.);
            nozzle.responses['VELOCITY'][1].append(-1.);
            nozzle.responses['VELOCITY'][2].append(-1.);

    if output == 'verbose':
        sys.stdout.write('SU2 responses obtained\n');
    
    return 0;  
        
 
def CheckConvergence ( nozzle ) :    
    
    plot_format      = nozzle.cfd.output_format;
    plot_extension   = SU2.io.get_extension(plot_format)
    history_filename = nozzle.cfd.conv_filename + plot_extension
    #special_cases    = SU2.io.get_specialCases(config)
    
    history      = SU2.io.read_history( history_filename )
    
    plot = SU2.io.read_plot(history_filename);
    
    RhoRes = history['Res_Flow[0]'];
    NbrIte = len(RhoRes);
    
    IniRes = RhoRes[0];
    FinRes = RhoRes[NbrIte-1];
    
    #print "Initial res = %le, Final res = %lf, DIFF = %lf\n" % (IniRes, FinRes, ResDif);
    return IniRes, FinRes;

def HF_Compute_Thrust_Wrap (nozzle):

    
    try:
        #multif.HIGHF.HF_GenerateExitMesh(nozzle);
        multif.HIGHF.HF_GenerateExitMesh_Direct(nozzle);
    except ImportError as e: 
        print 'Error '
        print e
        print
        raise
    
    P0  = nozzle.environment.P;
    M0  = nozzle.mission.mach;
    Gam = nozzle.fluid.gam;
    Rs  = nozzle.fluid.R;
    T0  = nozzle.environment.T;
    U0  = M0*np.sqrt(Gam*Rs*T0);
    
    options = { 'mesh_name' : nozzle.cfd.mesh_name,   \
                'restart_name' : nozzle.cfd.restart_name,   \
                'P0'          : P0 ,  \
                'M0'          : M0 ,  \
                'Gam'          : Gam,  \
                'Rs'          : Rs ,  \
                'T0'          : T0 ,  \
                'U0'          : U0   \
               }
    
    return HF_Compute_Thrust (options);
    


def HF_InterpolateToStructural (options):
    
    info = [];
    Crd  = [];
    Tri  = [];
    Tet  = [];
    Sol  = [];
    Header = [];
    
    structNam = "./visu.meshb"
    
    structNamLoc = "%s/structural.meshb" % (os.getcwd());
    
    try :
        os.symlink(structNam, structNamLoc);
    except:
        # ---
        sys.stdout.write("%s already exists\n" % exitNamLoc);
    
    out = _mshint_module.py_Interpolation (structNamLoc, options["mesh_name"], options["restart_name"],\
    info, Crd, Tri, Tet, Sol, Header);
    
    dim    = info[0];
    NbrVer = info[1]; 
    NbrTri = info[2];
    NbrTet = info[3];
    SolSiz = info[4];
    
    return Crd, Tri, Sol, SolSiz
    

def HF_Integrate_Sol_Wall(nozzle):
    
    MshNam = "nozzle.su2";
    SolNam = "nozzle.dat"
    
    #--- Extract surface patches from fluid mesh
    
    pyVer = [];
    pyTri = [];
    pySol = [];
    Ref = [9,10];
    
    _meshutils_module.py_ExtractSurfacePatches (MshNam, SolNam, pyVer, pyTri, pySol, Ref);
    
    pyVer = map(float, pyVer);
    pyTri = map(int, pyTri);
    pySol = map(float, pySol);
    
    NbrVer = len(pyVer)/3;
    Ver = np.array(pyVer).reshape(NbrVer,3).tolist();
    
    NbrTri = len(pyTri)/4;
    Tri = np.array(pyTri).reshape(NbrTri,4).tolist();
    
    SolSiz = len(pySol)/NbrVer;
    Sol = np.array(pySol).reshape(NbrVer,SolSiz).tolist();
    
    
    iPres = 5;
    iTemp = 6;
    
    if nozzle.method == 'RANS':
        iPres += 2;
        iTemp += 2;
    
    PresAvg = 0.0;
    TempAvg = 0.0;
    AreaTot = 0.0;
    
    Coord = [0,0,0]
    for iTri in range(NbrTri):
        
        pres = 0.0;
        temp = 0.0;
        
        for j in range(3):
            iVer = Tri[iTri][j]-1;
            Coord[j]  = Ver[iVer]; 
            pres += Sol[iVer][iPres];
            temp += Sol[iVer][iTemp];

        pres /= 3;
        temp /= 3;
        
        area = multif.HIGHF.hf_meshgeneration.GetTriangleArea3D(Coord);
        AreaTot += area;
        
        PresAvg += area*pres;
        TempAvg += area*temp;
    
    PresAvg /= AreaTot;
    TempAvg /= AreaTot;
    return AreaTot, PresAvg, TempAvg;


#
#def HF_Compute_Thrust (options):
#    
#    info = [];
#    Crd  = [];
#    Tri  = [];
#    Tet  = [];
#    Sol  = [];
#    Header = [];
#    
#    exitNam = "%s/baseline_meshes/nozzle_exit.mesh" % (os.path.dirname(os.path.abspath(__file__)));
#    exitNamLoc = "./nozzle_exit_hin.mesh";
#    
#    try :
#        os.symlink(exitNam, exitNamLoc);
#    except:
#        # ---
#        sys.stdout.write("%s already exists\n" % exitNamLoc);
#    
#    out = _mshint_module.py_Interpolation (exitNamLoc, options["mesh_name"], options["restart_name"],\
#    info, Crd, Tri, Tet, Sol, Header);
#    
#    dim    = info[0];
#    NbrVer = info[1]; 
#    NbrTri = info[2];
#    NbrTet = info[3];
#    SolSiz = info[4];
#    
#    #for i in range(40):
#    #    print "SOL %d = %.5le" % (i, Sol[i])
#    #
#    #for i in range(1,4):
#    #    idx = (i-1)*dim;
#    #    idxs = (i-1)*SolSiz;
#    #    print "Ver %d : %lf %lf %lf / idxs %d  Sol = %le %le ..." % (i, Crd[idx], Crd[idx+1], Crd[idx+2], idxs, Sol[idxs+0], Sol[idxs+1])
#    
#    # --- Get solution field indices
#    
#    idHeader = dict();
#    for iFld in range(0,len(Header)):
#        idHeader[Header[iFld]] = iFld+1;
#    
#    #iMach  = idHeader['Mach'];
#    #iTem   = idHeader['Temperature'];
#    #iCons1 = idHeader['Conservative_1'];
#    #iCons2 = idHeader['Conservative_2'];
#    #iCons3 = idHeader['Conservative_3'];
#    #iCons4 = idHeader['Conservative_4'];
#    #iPres  = idHeader['Pressure'];
#    
#    # New keywords (> SU2 Raven 5.0)
#    
#    iMach  = idHeader['Mach'];
#    iTem   = idHeader['Temperature'];
#    iCons1 = idHeader['Density'];
#    iCons2 = idHeader['X-Momentum'];
#    iCons3 = idHeader['Y-Momentum'];
#    iCons4 = idHeader['Energy'];
#    iPres  = idHeader['Pressure'];
#      
#    # --- Compute thrust    
#  
#    Thrust = 0.0;
#  
#    #P0  = nozzle.environment.P;
#    #M0  = nozzle.mission.mach;
#    #Gam = nozzle.fluid.gam;
#    #Rs  = nozzle.fluid.R;
#    #T0  = nozzle.environment.T;
#    #U0  = M0*np.sqrt(Gam*Rs*T0);
#    
#    P0  = options["P0"];
#    M0  = options["M0"];
#    Gam = options["Gam"];
#    Rs  = options["Rs"];
#    T0  = options["T0"];
#    U0  = options["U0"];
#    
#    v = np.zeros([3,3]);
#    a = np.zeros(3);
#    b = np.zeros(3);
#    
#    areatot = 0;
#    
#    for iTri in range(1, NbrTri) :
#        idt = 3*(iTri-1);
#        
#        rho  = 0.0;
#        rhoU = 0.0;
#        Pres = 0.0;
#        Mach = 0.0;
#        Temp = 0.0;
#        
#        for j in range(0,3):
#            iVer = int(Tri[idt+j]);
#            
#            idv = 3*(iVer-1);
#            ids = SolSiz*(iVer-1);
#            
#            for d in range(0,3):
#                v[j][d] = Crd[idv+d];
#            
#            
#            #if iTri == 1:
#            #    print "rho %d (%d) = %.3le (ids %d)" % (j,iVer,Sol[ids+iCons1], ids)
#                
#            rho  = rho  + Sol[ids+iCons1];
#            rhoU = rhoU + Sol[ids+iCons2];
#            Pres = Pres + Sol[ids+iPres];
#            Mach = Mach + Sol[ids+iMach];
#            Temp = Temp + Sol[ids+iTem];    
#        
#        us3 = 1./3.;
#        rho  = us3*rho; 
#        rhoU = us3*rhoU; 
#        Pres = us3*Pres; 
#        Mach = us3*Mach; 
#        Temp = us3*Temp; 
#        
#        # --- Compute triangle area
#        
#        for d in range(0,3):
#            a[d] = v[1][d] - v[0][d];
#            b[d] = v[2][d] - v[0][d];
#        
#        area = 0.;
#        
#        area = area + (a[1]*b[2]-a[2]*b[1])*(a[1]*b[2]-a[2]*b[1]);
#        area = area + (a[2]*b[0]-a[0]*b[2])*(a[2]*b[0]-a[0]*b[2]);
#        area = area + (a[0]*b[1]-a[1]*b[0])*(a[0]*b[1]-a[1]*b[0]);
#        
#        area = 0.5*np.sqrt(area);
#        
#        areatot += area;
#                
#        # --- Compute thrust
#                    
#        if rho < 1e-30:
#            sys.stderr.write("## ERROR HF_Compute_Thrust: rho == 0 for triangle %d\n" % iTri);
#            sys.exit(1);
#        
#        U = rhoU/rho;
#        
#        if iTri < 10:
#            print "Tri %d Pres %lf rho %lf U %lf rhoU %lf U0 %lf P0 %lf area %lf" % (iTri, Pres, rho, U, rhoU, U0, P0, area)
#        
#        try :
#            Thrust += area*(rhoU*(U-U0)+Pres-P0);
#        except :
#            Thrust = -1;
#        #Thrust = Thrust + area*(rhoU*(U*U0-U0)+P0*Pres-P0);
#    
#    Thrust = 2*Thrust; # symmetry
#    
#    print "AREA TOT %lf " % areatot;
#    
#    return Thrust;
#    
#

def HF_Compute_Thrust (options):
        
    MshNam = "nozzle.su2"
    SolNam = "nozzle.dat"
    
    #--- Extract surface patches from fluid mesh
    
    pyVer = [];
    pyTri = [];
    pySol = [];
    Ref = [19];
    
    _meshutils_module.py_ExtractSurfacePatches (MshNam, SolNam, pyVer, pyTri, pySol, Ref);
    
    pyVer = map(float, pyVer);
    pyTri = map(int, pyTri);
    pySol = map(float, pySol);
    
    NbrVer = len(pyVer)/3;
    Ver = np.array(pyVer).reshape(NbrVer,3).tolist();
    
    NbrTri = len(pyTri)/4;
    Tri = np.array(pyTri).reshape(NbrTri,4).tolist();
    
    SolSiz = len(pySol)/NbrVer;
    Sol = np.array(pySol).reshape(NbrVer,SolSiz).tolist();
        
    iCons1 = 0; #idHeader['Density'];
    iCons2 = 1; #idHeader['X-Momentum'];
    iCons3 = 2; #idHeader['Y-Momentum'];
    iCons4 = 3; #idHeader['Energy'];
    iPres  = 5; #idHeader['Pressure'];
    
    #for iVer in range(10):
    #    
    #    dens = Sol[iVer][iCons1];
    #    
    #    velx = Sol[iVer][iCons2]/dens;
    #    vely = Sol[iVer][iCons3]/dens;
    #    velz = Sol[iVer][iCons4]/dens;
    #    
    #    vel  = math.sqrt(velx*velx+vely*vely+velz*velz);
    #    pres = Sol[iVer][iPres];
    #    
    #    print "%d %lf %lf %lf  %lf  %lf  %lf" % (iVer, Ver[iVer][0], Ver[iVer][1], Ver[iVer][2], dens, vel, pres);
    
    #--- Compute Thrust
    
    pres_inf  = options["P0"];
    vel_inf   = options["U0"];
    
    Thrust = 0.0;
    
    v = np.zeros([3,3]);
    a = np.zeros(3);
    b = np.zeros(3);
    
    areatot = 0;
    
    for iTri in range(NbrTri) :
            
        #--- Compute triangle area
        
        for j in range(0,3):
            iVer = int(Tri[iTri][j])-1;
            
            for d in range(0,3):
                v[j][d] = float(Ver[iVer][d]);
        
        for d in range(0,3):
            a[d] = v[1][d] - v[0][d];
            b[d] = v[2][d] - v[0][d];
        
        area = 0.;
        area = area + (a[1]*b[2]-a[2]*b[1])*(a[1]*b[2]-a[2]*b[1]);
        area = area + (a[2]*b[0]-a[0]*b[2])*(a[2]*b[0]-a[0]*b[2]);
        area = area + (a[0]*b[1]-a[1]*b[0])*(a[0]*b[1]-a[1]*b[0]);
        area = 0.5*np.sqrt(area);
        
        areatot += area;
    
        #--- Increment thrust value
        
        for j in range(0,3):
            
            iVer = int(Tri[iTri][j])-1;
            
            dens = Sol[iVer][iCons1];
            
            velx = Sol[iVer][iCons2]/dens;
            vely = Sol[iVer][iCons3]/dens;
            velz = Sol[iVer][iCons4]/dens;
            
            vel  = math.sqrt(velx*velx+vely*vely+velz*velz);
            pres = Sol[iVer][iPres];
            
            Thrust += area/3.0*(dens*vel*(vel-vel_inf)+pres-pres_inf);
    
    Thrust = 2*Thrust; # symmetry
    
    print "THRUST %lf" % Thrust
    print "AREA TOT %lf " % areatot;
    
    return Thrust;
    
