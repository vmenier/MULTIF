import os, time, sys, shutil, copy, math
import numpy as np
#import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
from scipy import interpolate as itp
import subprocess
from multif import _meshutils_module
from multif import _mshint_module


# --- hf_SurfaceFluidMeshFull:
#   Extracts nozzle wall surface mesh
#   Doubles the mesh (not just one symmetric part)
#   Writes the full surface mesh
def hf_SurfaceFluidMeshFull (MshNam, SolNam, BasNamOut):
    
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
    
    #--- Double the mesh for interpolation
    
    tag = np.zeros(NbrVer);
    
    for i in range(NbrVer):
        
        y = Ver[i][1];
        
        if y > 1e-12:
            # Create a new vertex
            Ver.append([Ver[i][0],-Ver[i][1],Ver[i][2]])
            Sol.append(Sol[i]);
            tag[i] = len(Ver)-1;
        else:
            tag[i] = i;
        
    for i in range(NbrTri):
        Tri.append([tag[Tri[i][0]-1]+1,tag[Tri[i][1]-1]+1,tag[Tri[i][2]-1]+1,Tri[i][3]]);
    
    #--- Write out mesh
        
    Ver = np.array(Ver).reshape(3*len(Ver)).tolist();
    Tri = np.array(Tri).reshape(4*len(Tri)).tolist();
    Sol = np.array(Sol).reshape(SolSiz*len(Sol)).tolist();
    
    for i in range(len(Tri)-10,len(Tri)):
        print Tri[i]
    
    Ver = map(float, Ver);
    Tri = map(int,   Tri);
    Sol = map(float, Sol);
    Tet = [];
    Edg = [];
    
    _meshutils_module.py_WriteMesh("%s"%BasNamOut, Ver, Tri, Tet, Edg, Sol);
    

def hf_FluidStructureInterpolation(MshNam_str, MshNam_cfd, SolNam_cfd):
    
    refTab = [85, 79, 55, 67, 61, 49, 37, 31, 43];
    
    #MshNam_cfd = "nozzle.su2";
    #SolNam_cfd = "nozzle.dat";
    BasNamOut  = "nozzle_cfd_full"
    
    #--- Extract nozzle wall surface mesh + solution
    #    and double the mesh along the y direction (i.e. no symmetry)
    
    hf_SurfaceFluidMeshFull (MshNam_cfd, SolNam_cfd, BasNamOut);
    
    #--- Interpolate solution
    
    info   = [];
    Crd    = [];
    Tri    = [];
    Tet    = [];
    Sol    = [];
    Header = [];
    
    out = _mshint_module.py_Interpolation (MshNam_str, "%s.meshb"%BasNamOut, "%s.solb"%BasNamOut,\
        info, Crd, Tri, Tet, Sol, Header);
        
    #--- Get Pres and Temp indices
    
    header=np.genfromtxt(SolNam_cfd,dtype='str', max_rows=1)
    idHeader = dict();
    for i in range(len(header)):
        idHeader[header[i].strip('\"')] = i;
    
    iPres = idHeader['Pressure']-3;
    iTemp = idHeader['Temperature']-3;
        
    #--- Return Pres and Temp
    
    dim    = info[0];
    NbrVer = info[1]; 
    NbrTri = info[2];
    NbrTet = info[3];
    SolSiz = info[4];
    
    Ver = np.array(Crd).reshape(NbrVer,3).tolist();
    Tri = np.array(Tri).reshape(NbrTri,3).tolist();
    Sol = np.array(Sol).reshape(NbrVer+1,SolSiz).tolist();
    
    Pres = np.zeros(NbrVer+1);
    Temp = np.zeros(NbrVer+1);
    
    for i in range(1,NbrVer):
        Pres[i] = Sol[i-1][iPres];
        Temp[i] = Sol[i-1][iTemp];
    
    
    return Ver, Tri, Pres, Temp;
    