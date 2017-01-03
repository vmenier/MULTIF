# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import multif
from .. import SU2

from .. import _meshutils_module
import ctypes
import numpy as np

from .. import nozzle as nozzlemod


def PostProcessing (nozzle):
    
    # --- Check residual convergence
    
    IniRes, FinRes = CheckConvergence(nozzle);
    ResDif = FinRes - IniRes;
    sys.stdout.write("Initial res = %le, Final res = %lf, Diff = %lf\n" % (IniRes, FinRes, ResDif));
    
    # --- Interpolate and extract all the solution fields along the nozzle exit
    # SolExtract : Solution (x, y, sol1, sol2, etc.)
    # Size : [NbrVer, SolSiz]
    # Header : Name of the solution fields (Conservative_1, Mach, etc.)
    
    SolExtract, Size, Header  = ExtractSolutionAtExit(nozzle);
    
    nozzle.thrust = -1;
    nozzle.volume = -1;
    
    if nozzle.GetOutput['THRUST'] == 1:            
        nozzle.thrust = ComputeThrust ( nozzle, SolExtract, Size, Header );
        
    if nozzle.GetOutput['VOLUME'] == 1 or nozzle.GetOutput['MASS'] == 1 or nozzle.GetOutput['MASS_WALL_ONLY'] == 1:
         volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle);
         nozzle.volume = np.sum(volume);
         nozzle.mass = np.sum(mass);
         n_layers = len(nozzle.wall.layer);
         nozzle.mass_wall_only = np.sum(mass[:n_layers]);
         

def CheckConvergence ( nozzle ) :
    
    # filenames
    plot_format      = nozzle.OUTPUT_FORMAT;
    plot_extension   = SU2.io.get_extension(plot_format)
    history_filename = nozzle.CONV_FILENAME + plot_extension
    #special_cases    = SU2.io.get_specialCases(config)
    
    history      = SU2.io.read_history( history_filename )
    
    plot = SU2.io.read_plot(history_filename);
    
    RhoRes = history['Res_Flow[0]'];
    NbrIte = len(RhoRes);
    
    IniRes = RhoRes[0];
    FinRes = RhoRes[NbrIte-1];
        
    #print "Initial res = %le, Final res = %lf, DIFF = %lf\n" % (IniRes, FinRes, ResDif);
    return IniRes, FinRes;
    
    
def ExtractSolutionAtExit ( nozzle ):

	mesh_name    = nozzle.mesh_name;
	restart_name = nozzle.restart_name;
	
	pyResult = [];
	pyInfo   = [];
	pyHeader = [];
	
	pyBox = [nozzle.length,nozzle.length,0,nozzle.height+1e-20];
	
	_meshutils_module.py_ExtractAlongLine (mesh_name, restart_name, pyBox, pyResult, pyInfo, pyHeader);
	
	NbrRes = pyInfo[0];
	ResSiz = pyInfo[1];
	
	Result = np.asarray(pyResult);
	
	OutResult = np.reshape(Result,(NbrRes, ResSiz));
	
	Out_sort = OutResult[OutResult[:,1].argsort()]
	
	return Out_sort, pyInfo, pyHeader;
	

def ComputeThrust ( nozzle, SolExtract, Size, Header )	:
	
	# T = 2PI * Int_{0}^{R} (rho U ( U - U0) + P - Po ) r dr
	
	NbrVer = Size[0];
	SolSiz = Size[1];
	
	if len(Header) != SolSiz-2:
		sys.stderr.write("  ## ERROR : ComputeThrust : Inconsistent solution header.\n");
		sys.exit(0);
	
	# --- Get solution field indices
	
	iMach  = -1;
	iTem   = -1;
	iCons1 = -1;
	iCons2 = -1;
	iCons3 = -1;
	iCons4 = -1;
	iPres  = -1;
	
	for iFld in range(0,len(Header)):
		if Header[iFld] == 'Mach':
			iMach = iFld;
		elif Header[iFld] == 'Temperature':
			iTem = iFld;
		elif Header[iFld] == 'Conservative_1':
			iCons1 = iFld;
		elif Header[iFld] == 'Conservative_2':
			iCons2 = iFld;
		elif Header[iFld] == 'Conservative_3':
			iCons3 = iFld;
		elif Header[iFld] == 'Conservative_4':
			iCons4 = iFld;
		elif Header[iFld] == 'Pressure':
			iPres = iFld;
	
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
	
	for iVer in range(1, NbrVer) :
		
		y    = SolExtract[iVer][1];
		rho  = SolExtract[iVer][2+iCons1];
		rhoU = SolExtract[iVer][2+iCons2];
		Pres = SolExtract[iVer][2+iPres];
		Mach = SolExtract[iVer][2+iMach];
		Temp = SolExtract[iVer][2+iTem];
		
		
		U = rhoU/rho;
		
		dy = y - SolExtract[iVer-1][1];
		
		print "%lf %lf %lf %lf %lf %lf" % (y, rho, rhoU, Pres, Mach, Temp);
				
		Thrust = Thrust + dy*(rhoU*(U-U0)+Pres-P0);
		
	return Thrust;
	
def ExtractSolutionAtWall (nozzle):
    
    # --- Extract CFD solution at the inner wall
    
    mesh_name    = nozzle.mesh_name;
    restart_name = nozzle.restart_name;
    
    pyResult = [];
    pyInfo   = [];
    pyHeader = [];
    
    pyRef = [7];
    
    _meshutils_module.py_ExtractAtRef (mesh_name, restart_name, pyRef, pyResult, pyInfo, pyHeader);
    
    NbrRes = pyInfo[0];
    ResSiz = pyInfo[1];
        
    Result = np.asarray(pyResult);
    
    OutResult = np.reshape(Result,(NbrRes, ResSiz));
        
    #print "%d results" % NbrRes;
    #
    #for i in range(0,10):
    #    print "%d : (%lf,%lf) : rho = %lf" % (i, OutResult[i][0], OutResult[i][1], OutResult[i][2]);
    
    # --- Get solution field indices
    
    iMach  = -1;
    iTem   = -1;
    iCons1 = -1;
    iCons2 = -1;
    iCons3 = -1;
    iCons4 = -1;
    iPres  = -1;
        
    idHeader = dict();
    for iFld in range(0,len(pyHeader)):
        idHeader[pyHeader[iFld]] = iFld+2;
    
    
    return OutResult, pyInfo, idHeader;
    
    
    
    
