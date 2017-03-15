# -*- coding: utf-8 -*-


"""

R. Fenrich & V. Menier, July 2016

""" 

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
from .. import SU2
import multif

import material
import component
import inlet
import environment
import fluid
import mission
import tolerance
import geometry

from .. import _meshutils_module
import ctypes
import numpy as np

from parserDV import *

class Nozzle:
    def __init__(self):
        pass
    
    def AssignListDV(self,config,dvList,key,NbrDV,sizeDV):
        if key in config:
    
            dv = config[key].strip('()');
            dv = dv.split(",");
            dv_size = len(dv);
            
            if  dv_size != sizeDV:
                sys.stderr.write('  ## ERROR : Inconsistent number of '      \
                  'design variables found in %s: %d instead of %d\n\n' %     \
                  (key,dv_size,sizeDV));
                sys.exit(0);
            
            tag = np.zeros(NbrDV);
            for iw in range(0,dv_size):
    
                val = int(dv[iw]);
    
                if ( val > NbrDV or val < 0 ):
                    sys.stderr.write('  ## ERROR : Inconsistent design '     \
                      'variables provided in %s (idx=%d). Check DV_LIST '    \
                      'and %s_DV keyword specifications.\n\n' % (key,val));
                    sys.exit(0);
    
                dvList.append(val-1);
    
                if val-1 >= 0 :
                    tag[val-1] = 1;
                
            for iw in range(0,NbrDV):
                if tag[iw] != 1 :
                    sys.stderr.write('  ## ERROR : Design variable %d '      \
                      'in %s is not defined.\n\n' % (iw+1,key));
                    sys.exit(0);
        
        else :
            sys.stderr.write('\n ## ERROR : Expected %s in config file.'     \
              '\n\n' % key);
            sys.exit(0);    

        
    def SetupDV (self, config, output='verbose'):
        nozzle = self;
        
        NbrDVTot = 0;                
        if 'DV_LIST' in config:

            # Build a list of all possible expected keys. The list contains
            # other lists which represent possible more specific keys that can
            # be specified.
            dv_keys = list();
            dv_n = list(); # record allowable number of dv for each key
            dv_keys.append(['WALL']);
            dv_n.append([len(nozzle.wall.coefs)]);
            #dv_keys.append(['WALL_CENTERLINE']);
            #dv_n.append([len(nozzle.wall.centerline.coefs)]);
            if hasattr(nozzle.wall,'temperature'):
                if nozzle.wall.temperature.param == 'PIECEWISE_LINEAR':
                    ltemp = ['WALL_TEMP', 'WALL_TEMP_LOCATIONS',
                             'WALL_TEMP_VALUES'];
                    dv_keys.append(ltemp);
                    dv_n.append([2*len(nozzle.wall.temperature.thicknessNodes),
                                 len(nozzle.wall.temperature.thicknessNodes),
                                 len(nozzle.wall.temperature.thicknessNodes)]);
                elif nozzle.wall.temperature.param == 'KLE':
                    dv_keys.append('WALL_TEMP_COEFS');
                    dv_n.append(len(nozzle.wall.tempearture.coefs));
            for i in range(len(nozzle.wall.layer)): # assume piecewise linear
                if nozzle.wall.layer[i].param == 'PIECEWISE_LINEAR':
                    ltemp = [nozzle.wall.layer[i].name, 
                             nozzle.wall.layer[i].name + '_THICKNESS_LOCATIONS', 
                             nozzle.wall.layer[i].name + '_THICKNESS_VALUES'];
                    dv_keys.append(ltemp);
                    dv_n.append([2*len(nozzle.wall.layer[i].thicknessNodes),
                                 len(nozzle.wall.layer[i].thicknessNodes),
                                 len(nozzle.wall.layer[i].thicknessNodes)]);
                elif nozzle.wall.layer[i].param == 'CONSTANT':
                    ltemp = [nozzle.wall.layer[i].name + '_THICKNESS'];
                    dv_keys.append(ltemp);
                    dv_n.append([1]);
            dv_keys.append(['BAFFLES','BAFFLES_LOCATION','BAFFLES_THICKNESS',
                            'BAFFLES_HEIGHT']);
            dv_n.append([3*nozzle.baffles.n,nozzle.baffles.n,
                        nozzle.baffles.n,nozzle.baffles.n]);
            dv_keys.append(['STRINGERS','STRINGERS_BREAK_LOCATIONS',
                            'STRINGERS_HEIGHT_VALUES',
                            'STRINGERS_THICKNESS_VALUES']);
            dv_n.append([3*len(nozzle.stringers.thicknessNodes),
                         len(nozzle.stringers.thicknessNodes),
                         len(nozzle.stringers.thicknessNodes),
                         len(nozzle.stringers.thicknessNodes)]);
            for k in nozzle.materials:
                ltemp = [k, k + '_DENSITY', k + '_ELASTIC_MODULUS',
                         k + '_SHEAR_MODULUS', k + '_POISSON_RATIO',
                         k + '_MUTUAL_INFLUENCE_COEFS',
                         k + '_THERMAL_CONDUCTIVITY', 
                         k + '_THERMAL_EXPANSION_COEF',
                         k + '_PRINCIPLE_FAILURE_STRAIN',
                         k + '_LOCAL_FAILURE_STRAIN',
                         k + '_YIELD_STRESS',
                         k + '_MAX_SERVICE_TEMPERATURE']
                dv_keys.append(ltemp);
                if nozzle.materials[k].type == 'ISOTROPIC':
                    dv_n.append([[2,5],1,1,0,1,0,1,1,[1,2],[1,2],1,1]);
                else: # ANISOTROPIC_SHELL
                    dv_n.append([12,1,2,1,1,2,[2,3],3,5,5,[1,2],1]);              
            dv_keys.append(['INLET_PSTAG']);
            dv_n.append([1]);
            dv_keys.append(['INLET_TSTAG']);
            dv_n.append([1]);
            dv_keys.append(['ATM_PRES']);
            dv_n.append([1]);
            dv_keys.append(['ATM_TEMP']);
            dv_n.append([1]);
            dv_keys.append(['HEAT_XFER_COEF_TO_ENV']);
            dv_n.append([1]);
            
            # Extract listed design variables
            hdl = config['DV_LIST'].strip('()');
            hdl = [x.strip() for x in hdl.split(",")];

            dv_keys_size = len(hdl)/2;
            
            # First check that design variables have not been overspecified
            for k in dv_keys:
                if len(k) > 1:
                    if k[0] in hdl:
                        for i in range(1,len(k)):
                            if k[i] in hdl:
                                sys.stderr.write('\n ## ERROR: %s and %s '  \
                                'cannot both be specified in DV_LIST\n\n' % \
                                (k[0],k[i]));
                                sys.exit(0);
            
            nozzle.DV_Tags = []; # tags: e.g. wall, tstag etc.
            nozzle.DV_Head = []; # correspondance btw DV_Tags and DV_List
            
            for i in range(0,2*dv_keys_size,2):
                key = hdl[i];
                NbrDV = int(hdl[i+1]);

                # Check to make sure key is acceptable
                check = 0;
                strMsg = '';
                for j in dv_keys:
                    for k in j:
                        strMsg = strMsg + ',' + k;
                        if key == k:
                            check = 1;
                if check != 1:
                    sys.stderr.write('\n ## ERROR : Unknown design variable ' \
                      'key : %s\n\n' % key);
                    sys.stderr.write('            Expected = %s\n\n' % strMsg);
                    sys.exit(0);
                
                # Append important information
                nozzle.DV_Tags.append(key);
                nozzle.DV_Head.append(NbrDVTot);   
                 
                NbrDVTot = NbrDVTot + NbrDV;                    

                if key == 'WALL':                    
                    nozzle.wall.dv = [];
                    sizeTemp = nozzle.wall.coefs_size; # number of DV for checking
                    self.AssignListDV(config,nozzle.wall.dv,'WALL_COEFS_DV',
                                      NbrDV,sizeTemp);
                    continue;
                    
                if key == 'WALL_TEMP':
                    nozzle.wall.temperature.dv = [];
                    sizeTemp = 2*len(nozzle.wall.temperature.thicknessNodes);
                    self.AssignListDV(config,nozzle.wall.temperature.dv,
                                      'WALL_TEMP_DV',NbrDV,sizeTemp);
                    continue;
                
                # Check all layers with non-specific names, e.g. LAYER1, etc.
                check = 0;
                for j in range(len(nozzle.wall.layer)):
                    if key == nozzle.wall.layer[j].name:
                        nozzle.wall.layer[j].dv = [];
                        # number of DV for checking; assumes 2xN array
                        sizeTemp = 2*len(nozzle.wall.layer[j].thicknessNodes);
                        self.AssignListDV(config,nozzle.wall.layer[j].dv,
                                          nozzle.wall.layer[j].formalName+'_DV',
                                          NbrDV,sizeTemp);
                        check = 1;
                if check == 1:
                    continue;
                    
                if key == 'BAFFLES':                    
                    nozzle.baffles.dv = [];
                    sizeTemp = nozzle.baffles.n*3; # number of DV for checking
                    self.AssignListDV(config,nozzle.baffles.dv,'BAFFLES_DV',
                                      NbrDV,sizeTemp);
                    continue;
                    
                if key == 'STRINGERS':
                    nozzle.stringers.dv = [];
                    sizeTemp = len(nozzle.stringers.thicknessNodes)*3; # number of DV for checking
                    self.AssignListDV(config,nozzle.stringers.dv,'STRINGERS_DV',
                                      NbrDV,sizeTemp);
                    continue;
                    
                # Check all design variables for the right quantity
                for j in range(len(dv_keys)):
                    for k in range(len(dv_keys[j])):
                        if key == dv_keys[j][k]:
                            if isinstance(dv_n[j][k], int):
                                if NbrDV != dv_n[j][k]:
                                    sys.stderr.write('\n ## ERROR : Inconsistent' \
                                      ' number of DV for %s definition (%d '      \
                                      'provided instead of %d).\n\n' %            \
                                      (key,NbrDV,dv_n[j][k]));
                                    sys.exit(0);
                            else: # if list is provided
                                check = 0;
                                for m in dv_n[j][k]:
                                    if NbrDV == m:
                                        check = 1;
                                if check == 0:
                                    sys.stderr.write('\n ## ERROR : Inconsistent' \
                                      ' number of DV for %s definition (%d '      \
                                      'provided instead of %d).\n\n' %            \
                                      (key,NbrDV,m));
                                    sys.exit(0);                                        
                
            # for i in keys
        
        else :
            sys.stdout.write('\n  -- Info : No design variable set was '      \
              'defined. Running the baseline parameters.\n\n');
        
        nozzle.NbrDVTot = NbrDVTot;
        
        if NbrDVTot > 0 :
            nozzle.DV_Head.append(NbrDVTot);
            
        if output == 'verbose':
            sys.stdout.write('Setup Design Variables complete\n');            
            

    def SetupFidelityLevels (self, config, flevel, output='verbose'):
        
        nozzle = self;
        
        fidelity_tags = config['FIDELITY_LEVELS_TAGS'].strip('()');
        fidelity_tags = fidelity_tags.split(",");

        nozzle.MaxCFL = 30.0;

        NbrFidLev = len(fidelity_tags);

        analysisType = -1;

        if NbrFidLev < 1:
            sys.stdout.write("  ## ERROR : No fidelity level was defined.\n");
            sys.exit(0);

        if flevel >= NbrFidLev:
          sys.stderr.write("\n  ## ERROR : Fidelity level %d not defined.\n" % flevel);
          sys.exit(0)					

        if output == 'verbose':
          sys.stdout.write('\n%d fidelity level(s) defined. Summary :\n'      \
            % (NbrFidLev));

          sys.stdout.write('-' * 90);
          sys.stdout.write('\n%s | %s | %s\n' % ("Level #".ljust(10),         \
            "Tag".ljust(10),"Description.".ljust(70)));
          sys.stdout.write('-' * 90);
          sys.stdout.write('\n');

        for i in range(NbrFidLev):
            tag = fidelity_tags[i];
            kwd = "DEF_%s" % tag;

            if kwd not in config:
                sys.stderr.write("\n  ## ERROR : The fidelity level tagged '  \
                  '%s is not defined.\n\n" % kwd);
                sys.exit(0);

            cfgLvl = config[kwd].strip('()');
            cfgLvl = cfgLvl.split(",");

            method = cfgLvl[0];

            description = "";

            if i == flevel:
                nozzle.method = method;

            if method == 'NONIDEALNOZZLE':

                tol = float(cfgLvl[1]);
                if tol < 1e-16:
                    sys.stderr.write("\n ## ERROR : Wrong tolerance for "     \
                     "fidelity level %d (tagged %s)\n\n" % (i,tag));
                    sys.exit(0);
                description += "Quasi-1D, ODE solver relative and absolute "  \
                    "tolerance set to %le" % (tol);
                        
                # Set thermostructural parameters if necessary
                if len(cfgLvl) == 5:
                    if cfgLvl[3] == 'LINEAR':
                        description += ", linear structural analysis"
                    elif cfgLvl[3] == 'NONLINEAR':
                        description += ", nonlinear structural analysis"
                    else:
                        sys.stderr.write('\n ## ERROR: Only LINEAR or '   \
                          'NONLINEAR can be specified for structural '    \
                          'analysis type. %s specified instead in.\n\n'   \
                          % cfgLvl[3]);
                        sys.exit(0); 
                    description += ', thermostructural fidelity level %s' % cfgLvl[4];
                        
                if i == flevel:                    
                    # Set tolerances
                    nozzle.tolerance = tolerance.Tolerance();
                    nozzle.tolerance.setRelTol(tol);
                    nozzle.tolerance.setAbsTol(tol);

                    # Set analysis type
                    try:
                        analysisType = cfgLvl[2];
                    except:
                        sys.stderr.write('\n ## WARNING : Analysis type '     \
                          'could not be determined. THERMOSTRUCTURAL, '       \
                          'THERMAL, or STRUCTURAL keyword must be provided '  \
                          'in the model definition of fidelity level %d.\n\n' \
                          % flevel);
                        sys.exit(0);
                        
                    # Set thermostructural parameters if necessary
                    if len(cfgLvl) == 5:
                        if cfgLvl[3] == 'LINEAR':
                            nozzle.linearStructuralAnalysis = 1;
                        elif cfgLvl[3] == 'NONLINEAR':
                            nozzle.linearStructuralAnalysis = 0;
                            
                        nozzle.thermostructuralFidelityLevel = float(cfgLvl[4]);
                        if( nozzle.thermostructuralFidelityLevel < 0 or 
                        nozzle.thermostructuralFidelityLevel > 1):
                            sys.stderr.write('\n ## ERROR: thermostructural ' \
                              'fidelity level must range from 0 (low) to 1 '\
                              '(high). %f provided instead.\n\n' % 
                              nozzle.thermostructuralFidelityLevel);
                            sys.exit(0);
                    else:
                        nozzle.linearStructuralAnalysis = 1;
                        nozzle.thermostructuralFidelityLevel = 0.5;
                        
            elif method == 'RANS' or method == 'EULER':
                dim = cfgLvl[1];
                if dim != '2D' and dim != '3D' :
                    sys.stderr.write("\n ## ERROR : Wrong dimension for "     \
                      "fidelity level %d (tagged %s) : only 2D or 3D "        \
                      "simulations\n\n" % (i,tag));
                    sys.exit(0);
                
                nozzle.Dim = dim;
                
                meshsize = cfgLvl[2];    
                if( meshsize != 'COARSE' 
                and meshsize != 'MEDIUM' 
                and meshsize != 'FINE' ):
                    sys.stderr.write("\n ## ERROR : Wrong mesh level for "    \
                      "fidelity level %d (tagged %s) : must be set to "       \
                      "either COARSE, MEDIUM or FINE" % (i,tag));
                    sys.exit(0);
                description += "%s %s CFD, %s mesh" \
                  % (dim, method, meshsize);
                  
                # --- Setup convergence parameter
                
                if 'SU2_CONVERGENCE_ORDER' in config:
                    nozzle.su2_convergence_order = int(config['SU2_CONVERGENCE_ORDER']);
                else:
                    nozzle.su2_convergence_order = 3;
                description += ", relative convergence order %i" % nozzle.su2_convergence_order;
                

                if 'SU2_OUTPUT_FORMAT' in config:
                    nozzle.OUTPUT_FORMAT = config['SU2_OUTPUT_FORMAT'];
                else:
                    nozzle.OUTPUT_FORMAT = 'TECPLOT';

                # --- Setup max iterations for SU2
                if 'SU2_MAX_ITERATIONS' in config:
                    nozzle.su2_max_iterations = int(config['SU2_MAX_ITERATIONS']);
                elif nozzle.method == 'EULER':
                    nozzle.su2_max_iterations = 300;
                else:
                    nozzle.su2_max_iterations = 2000;                  
                description += ", max iterations %i" % nozzle.su2_max_iterations;
                  
                # Set thermostructural parameters if necessary
                if len(cfgLvl) == 6:
                    if cfgLvl[4] == 'LINEAR':
                        description += ", linear structural analysis"
                    elif cfgLvl[4] == 'NONLINEAR':
                        description += ", nonlinear structural analysis"
                    else:
                        sys.stderr.write('\n ## ERROR: Only LINEAR or '   \
                          'NONLINEAR can be specified for structural '    \
                          'analysis type. %s specified instead in.\n\n'   \
                          % cfgLvl[4]);
                        sys.exit(0);                    
                    description += ', thermostructural fidelity level %s' % cfgLvl[5];
                    
                if i == flevel:
                    # Set mesh level
                    nozzle.meshsize = meshsize;
                    
                    nozzle.bl_ds        = 0.000007;
                    nozzle.bl_ratio     = 1.3; 
                    nozzle.bl_thickness = 0.02;
					
                    if ( method == "EULER" ):
                    	scaleMesh = 1.0;
                    	if meshsize == 'COARSE':
                    	    scaleMesh = 2.5;
                    	elif meshsize == 'MEDIUM':
                    	    scaleMesh = 1.1;
                    	elif meshsize == 'FINE':
                    	    scaleMesh = 0.5;
                    else : 
                    	scaleMesh = 1.0;
                    	if meshsize == 'COARSE':
                    	    scaleMesh = 2.5;
                    	elif meshsize == 'MEDIUM':
                    	    scaleMesh = 1.5;
                    	elif meshsize == 'FINE':
                    	    scaleMesh = 0.9;

                    nozzle.meshhl = scaleMesh*np.asarray([0.1, 0.07, 0.06, 0.006, 0.0108]);

                    # Set analysis type
                    try:
                        analysisType = cfgLvl[3];
                    except:
                        sys.stderr.write('\n ## WARNING : Analysis type '     \
                          'could not be determined. AEROTHERMOSTRUCTURAL, '   \
                          'AEROTHERMAL, AEROSTRUCTURAL, or AERO keyword must' \
                          ' be provided in the model definition of '          \
                          'fidelity level %d.\n\n' % flevel);
                        sys.exit(0);
                        
                    # Set thermostructural parameters if necessary
                    if len(cfgLvl) == 6:
                        if cfgLvl[4] == 'LINEAR':
                            nozzle.linearStructuralAnalysis = 1;
                        elif cfgLvl[4] == 'NONLINEAR':
                            nozzle.linearStructuralAnalysis = 0;      
                            
                        nozzle.thermostructuralFidelityLevel = float(cfgLvl[5]);
                        if( nozzle.thermostructuralFidelityLevel < 0 or 
                        nozzle.thermostructuralFidelityLevel > 1):
                            sys.stderr.write('\n ## ERROR: thermostructural ' \
                              'fidelity level must range from 0 (low) to 1 '\
                              '(high). %f provided instead.\n\n' % 
                              nozzle.thermostructuralFidelityLevel);
                            sys.exit(0);
                    else:
                        nozzle.linearStructuralAnalysis = 1;
                        nozzle.thermostructuralFidelityLevel = 0.5;
                            
            else :
                sys.stderr.write("\n ## ERROR : Unknown governing method "    \
                  "(%s) for fidelity level %s.\n\n" % (method, tag));
                sys.stderr.write("  Note: it must be either NONIDEALNOZZLE,"  \
                  "EULER, or RANS\n");
                sys.exit(0);

            if output == 'verbose':
                sys.stdout.write("   %s | %s | %s \n" %                       \
                  ( ("%d" % i).ljust(7), tag.ljust(10),                       \
                  textwrap.fill(description, 70,                              \
                  subsequent_indent="".ljust(26))) );
            elif output == 'quiet':
                pass
            else:
                raise ValueError('keyword argument output can only be set '   \
                  'to "verbose" or "quiet" mode')
                  
        # Setup thermal and structural analysis
        if analysisType == 'AEROTHERMOSTRUCTURAL':
            nozzle.thermalFlag = 1;
            nozzle.structuralFlag = 1;
        elif analysisType == 'AEROTHERMAL':
            nozzle.thermalFlag = 1;
            nozzle.structuralFlag = 0;
        elif analysisType == 'AEROSTRUCTURAL':
            nozzle.thermalFlag = 0;
            nozzle.structuralFlag = 1;
        elif analysisType == 'AERO':
            nozzle.thermalFlag = 0;
            nozzle.structuralFlag = 0;
        else:
            sys.stderr.write('\n ## ERROR: AEROTHERMOSTRUCTURAL, '        \
              'AEROTHERMAL, AEROSTRUCTURAL, or AERO must be '        \
              'provided as a keyword for analyis '    \
              'type. %s provided instead.\n\n' % analysisType);
            sys.exit(0);                  

        if output == 'verbose':
            sys.stdout.write('-' * 90);
            sys.stdout.write('\n\n');
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to '    \
              '"verbose" or "quiet" mode')

        if flevel >= NbrFidLev :
            sys.stderr.write("\n ## ERROR : Level %d not defined !! "         \
              "\n\n" % flevel);
            sys.exit(0);

        if output == 'verbose':
            sys.stdout.write('  -- Info : Fidelity level to be run : '        \
              '%d\n' % flevel);
            sys.stdout.write('Analysis type is %s.\n\n' % analysisType);  
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to '    \
              '"verbose" or "quiet" mode')
              
        if output == 'verbose':
            sys.stdout.write('Setup Fidelity Levels complete\n');              
        
    def SetupMission(self, config, output='verbose'):
    
        nozzle = self;
        
        if( ('MISSION' in config) and (('INLET_PSTAG' in config) 
        or ('ALTITUDE' in config) or ('INLET_TSTAG' in config) 
        or ('MACH' in config)) ):
      
            sys.stderr.write('\n ## ERROR : MISSION cannot be specified '     \
              'in conjunction with INLET_PSTAG, INLET_TSTAG, ALTITUDE, or '   \
              'MACH. Either MISSION or all four following quantities must '   \
              'be specified: INLET_PSTAG, INLET_TSTAG, ALTITUDE, MACH.\n\n');
            sys.exit(0);
      
        elif 'MISSION' in config:
        
          mission_id = int(config["MISSION"]);

          if(mission_id == 0): # standard top-of-climb, 40,000 ft case
              altitude = 40000. # ft
              mach     = 0.511
              inletTs  = 955.0 # K
              inletPs  = 97585. # Pa
          else:
              sys.stdout.write('\n ## WARNING : Previously available '        \
                'missions (1 through 5) are available for backwards '         \
                'compatability, but you should be using MISSION= 0 which '    \
                'corresponds to the conditions at max climb rate.\n\n');
              if(mission_id == 1): # static sea-level thrust case
                  altitude = 0.
                  mach     = 0.01
                  inletTs  = 888.3658
                  inletPs  = 3.0550e5
              elif(mission_id == 2): # intermediate case
                  altitude = 15000.
                  mach     = 0.5
                  inletTs  = 942.9857
                  inletPs  = 2.3227e5
              elif(mission_id == 3): # high speed, high altitude case
                  altitude = 35000.
                  mach     = 0.9
                  inletTs  = 1021.5
                  inletPs  = 1.44925e5
              elif(mission_id == 4): # case with shock in nozzle
                  altitude = 0.
                  mach     = 0.01
                  inletTs  = 900.
                  inletPs  = 1.3e5
              elif(mission_id == 5): # subsonic flow
                  altitude = 0.
                  mach     = 0.01
                  inletTs  = 900.
                  inletPs  = 1.1e5	
              else : 
                  sys.stderr.write('\n ## ERROR : UNKNOWN MISSION ID %d !! '     \
                                   '\n\n' % mission);
                  sys.exit(0);
 
        else: # user must specify altitude, mach, inletTs, and inletPs
          
          if( ('INLET_PSTAG' in config) and ('INLET_TSTAG' in config) 
          and ('ALTITUDE' in config) and ('MACH' in config)):
          
            mission_id = -1; # to denote custom mission
            altitude = float(config['ALTITUDE']);
            mach = float(config['MACH']);
            inletTs = float(config['INLET_TSTAG']);
            inletPs = float(config['INLET_PSTAG']);
          
          else:
          
            sys.stderr.write('\n ## ERROR : If MISSION is not specified '    \
                             'then INLET_PSTAG, INLET_TSTAG, ALTITUDE, and ' \
                             'MACH must all be specified\n\n');
            sys.exit(0);
        
        nozzle.mission = mission.Mission(mission_id);  
        nozzle.mission.setMach(mach);
        nozzle.inlet   = inlet.Inlet(inletPs,inletTs);

        # --- Setup heat transfer coeff. from external nozzle wall to env.

        if 'HEAT_XFER_COEF_TO_ENV' in config:
            hInf = float(config['HEAT_XFER_COEF_TO_ENV'].strip('()'));
        else:
            hInf = 12.62; # W/m^2/K
        if mission_id == 0: # set pressure/temp manually so exact pressure and temp are known
            nozzle.environment = environment.Environment(altitude,hInf);
            nozzle.environment.setPressure(18754.);
            nozzle.environment.setTemperature(216.7);
        else:
            nozzle.environment = environment.Environment(altitude,hInf);

        # --- Setup fluid

        heatRatio = 1.4;
        gasCst    = 287.06;
        nozzle.fluid = fluid.Fluid(heatRatio, gasCst);
            
        if output == 'verbose':
            sys.stdout.write('Setup Mission complete\n');
    
    
    def SetupBSplineCoefs(self, config, output='verbose'):
        # Assume B-spline is a 3rd-degree B-spline. Thus, given the coefs, the
        # knots can be calculated assuming evenly-spaced knots with a B-spline
        # that terminates at both ends and not earlier.
        
        nozzle = self;
        
        # Set up nozzle wall
        nozzle.wall = component.AxisymmetricWall();
        
        #wall_keys = ('WALL','WALL_KNOTS','WALL_COEFS');
        wall_keys = ('WALL','WALL_COEFS');

        if all (key in config for key in wall_keys):

            # --- Get coefs

            hdl = config['WALL_COEFS'].strip('()');
            hdl = hdl.split(",");

            coefs_size = len(hdl);

            coefs = [];

            for i in range(0,coefs_size):
                coefs.append(float(hdl[i]));

            # --- Get knots

#            hdl = config['GEOM_WALL_KNOTS'].strip('()');
#            hdl = hdl.split(",");
#
#            knots_size = len(hdl);
#
#            knots = [];
#
#            for i in range(0,knots_size):
#                knots.append(float(hdl[i]));

            #print "%d coefs, %d knots\n" % (coefs_size, knots_size);

        else:
            string = '';
            for key in wall_keys :
                string = "%s %s " % (string,key);
            sys.stderr.write('\n ## ERROR : NO INNER WALL DEFINITION IS '    \
                             'PROVIDED.\n\n');
            sys.stderr.write('             Expected = %s\n\n' % string);
            sys.exit(0);
            
        # Calculate knots using assumptions about 3rd degree B-spline
        n = len(coefs)/2-3;
        knots =  [0,0,0,0] + range(1,n) + [n,n,n,n];
        
        nozzle.wall.coefs =  coefs;
        nozzle.wall.knots = [float(i) for i in knots];
        nozzle.wall.coefs_size = coefs_size;
        
        if output == 'verbose':
            sys.stdout.write('Setup B-Spline Coefs complete\n');        


    def SetupWallTemp(self, config, output='verbose'): 
        
        nozzle = self;

        nozzle.wall_temp = 0;
        
        if 'WALL_TEMP' not in config:
            return 0;
        else:
            if output == 'verbose':
                sys.stdout.write('Wall temp will be fixed\n');
            nozzle.wall_temp = 1;
                
        definition = config['WALL_TEMP'].strip('()');
        if definition == 'KLE':
            
            sys.stderr.write('\n ## ERROR: KLE definition for wall temp not ' \
              'implemented yet.\n\n');
            sys.exit(0);
            
            #nozzle.wall.temperature = component.Distribution('WALL_TEMP');
            #nozzle.wall.temperature.param = 'KLE';
            #nozzle.wall.temperature.coefs = config['WALL_TEMP_COEFS'].strip('()');
            
        elif definition == 'PIECEWISE_LINEAR':
                        
            try:
                thicknessNodes = self.ParseThickness(config,'WALL_TEMP',loc='_LOCATIONS',val='_VALUES');
            except:
                try:
                    xCoordFile = config['WALL_TEMP_LOCATIONS'].strip('()');
                    yCoordFile = config['WALL_TEMP_VALUES'].strip('()');
                    xCoord = np.loadtxt(xCoordFile);
                    yCoord = np.loadtxt(yCoordFile);
                    thicknessNodes = [[0 for i in range(2)] for j in range(yCoord.size)];
                    
                    for i in range(xCoord.size):
                        thicknessNodes[i][0] = xCoord[i];
                        thicknessNodes[i][1] = yCoord[i];                    
                except:                        
                    sys.stderr.write('\n ## ERROR : Piecewise linear definition ' \
                         'could not be parsed for WALL_TEMP.\n\n');
                    sys.exit(0);
                
            nozzle.wall.temperature = component.Distribution('WALL_TEMP');
            nozzle.wall.temperature.param = 'PIECEWISE_LINEAR';
            nozzle.wall.temperature.thicknessNodes = thicknessNodes;                
                
        else:
            sys.stderr.write('\n ## ERROR: %s wall temp definition not '      \
              'accepted. Only KLE and PIECEWISE_LINEAR are implemented.\n\n');
            sys.exit(0);
            
        if output == 'verbose':
            sys.stdout.write('Setup Wall Temp complete\n'); 
            
    
    def ParseThickness(self, config, key, loc='_THICKNESS_LOCATIONS', val = '_THICKNESS_VALUES'):
        
        loc_name = key + loc;
        val_name = key + val;
        
        wall_keys = (loc_name, val_name);
        
        if all (key in config for key in wall_keys):
        
            hdl = config[loc_name].strip('()');
            hdl = hdl.split(",");
  
            size_loc = len(hdl);
          
            wall_thickness = [[0 for i in range(2)] for j in range(size_loc)];
                        
            for i in range(0,size_loc):
                wall_thickness[i][0] = float(hdl[i]);
                if wall_thickness[i][0] < 0.0 or wall_thickness[i][0] > 1.0:
                    sys.stderr.write('\n ## ERROR : Invalid wall thickness ' \
                                     'definition (%s).\n' % key);
                    sys.stderr.write('              All values must be '     \
                                     'between 0 and 1.\n\n');
                    sys.exit(0);
            
            hdl = config[val_name].strip('()');
            hdl = hdl.split(",");
            size_val = len(hdl);
            
            if size_val != size_loc :
                sys.stderr.write('\n ## ERROR : Inconsistent wall thickness' \
                                 'definition (%s).\n' % key);
                sys.stderr.write('              Same number of locations '   \
                                 'and values required.\n\n');
                sys.exit(0);
            
            for i in range(0,size_val):
                wall_thickness[i][1] = float(hdl[i]);
            
            if wall_thickness[0][0] != 0.0 or wall_thickness[size_val-1][0] != 1.0:
                sys.stderr.write('\n ## ERROR : Invalid wall thickness '     \
                                 'definition (%s).\n' % key);
                sys.stderr.write('              First and last loc values '  \
                                 'must be 0 and 1 resp.\n\n');
                sys.exit(0);
                
            return wall_thickness;
        
        else:
            raise; 

    
    # Setup and assign thickness distribution and material to a single layer
    def SetupLayer(self,config,layer,matStructure):
        
        # Determine name of layer
        name = layer.name;        
        if name != 'THERMAL_LAYER' and name != 'LOAD_LAYER':
            sys.stderror.write('\n ## ERROR : Only THERMAL_LAYER and '       \
                               'LOAD_LAYER are accepted layer names\n\n');
            sys.exit(0);
        
        # Setup thickness keys
        thickness_keys = ['_THICKNESS_LOCATIONS','_THICKNESS_VALUES'];
        for i in range(len(thickness_keys)):
            thickness_keys[i] = name + thickness_keys[i];
        thickness_keys = tuple(thickness_keys)
            
        # Setup material keys
        material_keys = ['_DENSITY','_ELASTIC_MODULUS','_POISSON_RATIO',
                         '_THERMAL_CONDUCTIVITY','_THERMAL_EXPANSION_COEF',
                         '_YIELD_STRESS','_PRINCIPLE_FAILURE_STRAIN',
                         '_LOCAL_FAILURE_STRAIN','_MAX_SERVICE_TEMPERATURE'];  
        for i in range(len(material_keys)):
            material_keys[i] = name + material_keys[i];
        material_keys = tuple(material_keys)
            
        # Update thickness distribution of layer
        if all (k in config for k in thickness_keys):
            try:
                layer.thicknessNodes = self.ParseThickness(config,name);
            except:
                sys.stderr.write('\n ## ERROR : Thickness definition '       \
                     'could not be parsed for %s.\n\n' % name);
                sys.exit(0);
        else:
            if name == 'LOAD_LAYER':
                layer.thicknessNodes = [[0.0,0.016], [1.0, 0.016]];
            else: # THERMAL_LAYER
                layer.thicknessNodes = [[0.0,0.01], [1.0, 0.01]];

        # Run through keys first and see if material is anisotropic
        isotropicFlag = 1 # assume isotropic
        for key in material_keys:
            var = config[key].strip('()');
            var = var.split(',');
            if len(var) == 2:
                isotropicFlag = 0; 
                break;
            elif len(var) == 1:
                pass;
            else:
                sys.stderr.write('\n ## ERROR: Key %s could not be parsed '  \
                'for %s layer' % (key,name));
                sys.exit(0);

        # Define material
        if isotropicFlag: # isotropic
            layer.material = material.Material('isotropic',matStructure);
        else:
            layer.material = material.Material('anisotropic',matStructure);
        
        # Now assign material properties to material
        for key in material_keys:
            
            var = config[key].strip('()');
            var = var.split(',');
            if len(var) >= 2:
                var = [float(i) for i in var]
            elif len(var) == 1:
                var = float(var[0]);
            else:
                sys.stderr.write('\n ## ERROR: Key %s could not be parsed '  \
                'for %s layer' % (key,name));
                sys.exit(0);
              
            if key == name + '_DENSITY':
                layer.material.setDensity(var);
            elif key == name + '_ELASTIC_MODULUS':
                layer.material.setElasticModulus(var);
            elif key == name + '_POISSON_RATIO':
                layer.material.setPoissonRatio(var);
            elif key == name + '_THERMAL_CONDUCTIVITY':
                layer.material.setThermalConductivity(var);
            elif key == name + '_THERMAL_EXPANSION_COEF':
                layer.material.setThermalExpansionCoef(var);
            elif ( key == name + '_YIELD_STRESS' or
                   key == name + '_PRINCIPLE_FAILURE_STRAIN' or
                   key == name + '_LOCAL_FAILURE_STRAIN'):
                layer.material.setFailureType(var); 
            elif key == name + '_MAX_SERVICE_TEMPERATURE':
                layer.material.setMaxServiceTemperature(var);
            else:
                sys.stderr.write('\n ## ERROR: Key %s could not implemented' \
                ' for %s layer' % (key,name));
                sys.exit(0);
    
    def SetupLayerThickness(self, config, layer):
        
        if layer.param == 'PIECEWISE_LINEAR':            
            # Setup thickness keys
            t_keys = ['_THICKNESS_LOCATIONS','_THICKNESS_VALUES'];
            for i in range(len(t_keys)):
                t_keys[i] = layer.formalName + t_keys[i];
            t_keys = tuple(t_keys)

            # Assign thicknesses
            if all (k in config for k in t_keys):
                try:
                    layer.thicknessNodes = self.ParseThickness(config,layer.formalName);
                except:
                    sys.stderr.write('\n ## ERROR : Thickness definition '   \
                         'could not be parsed for %s.\n\n' % layer.name);
                    sys.exit(0);
            else: # assume smart values for constant thickness layer
                if layer.name == 'THERMAL_LAYER':
                    layer.thicknessNodes = [[0.0,0.01], [1.0, 0.01]];
                elif layer.name == 'LOAD_LAYER_INSIDE':
                    layer.thicknessNodes = [[0.0,0.004], [1.0, 0.004]];
                elif layer.name == 'LOAD_LAYER_MIDDLE':
                    layer.thicknessNodes = [[0.0,0.008], [1.0, 0.008]];   
                elif layer.name == 'LOAD_LAYER_OUTSIDE':
                    layer.thicknessNodes = [[0.0,0.004], [1.0, 0.004]];                      
                else: # just assume a layer 1 cm thick
                    layer.thicknessNodes = [[0.0,0.01], [1.0, 0.01]];  
        
        elif layer.param == 'CONSTANT': # implement as piecewise linear
            t_value = float(config[layer.formalName + '_THICKNESS']);
            layer.thicknessNodes = [[0.0,t_value], [1.0, t_value]];
        else:
            sys.stderr.write('\n ## ERROR: Parameterization %s '     \
            'is not valid for %s. Only PIECEWISE_LINEAR or '    \
            'CONSTANT is accepted\n' % (layer.param, layer.name));
            sys.exit(0);
            
    def SetupMaterials(self, config,output='verbose'):
        
        nozzle = self;
        
        nozzle.materials = {};
 
        # Discover all materials and assign data
        i = 1;
        while i:
            key = 'MATERIAL' + str(i);
            if key in config:
                # Parse information related to layer
                info = config[key].strip('()');
                info = [x.strip() for x in info.split(',')];
                if len(info) == 2:
                    name = info[0];
                    type_prop = info[1]; # ISOTROPIC or ANISOTROPIC
                    
                    if type_prop == 'ISOTROPIC' or type_prop == 'ANISOTROPIC_SHELL':
                        
                        nozzle.materials[name] = material.Material(name,type_prop,'single');
                        
                    elif type_prop == 'FIXED_RATIO_PANEL':
                        nozzle.materials[name] = material.Material(name,type_prop,'panel');
                        keyLayers = 'MATERIAL' + str(i) + '_LAYERS';
                        keyRatios = 'MATERIAL' + str(i) + '_THICKNESS_RATIOS';
                        if keyLayers in config:
                            layerNames = config[keyLayers].strip('()');
                            layerNames = [x.strip() for x in layerNames.split(',')];
                        else:
                            sys.stderr.write('\n ## ERROR: %s not found in '  \
                              'config file. Layer names must be defined.\n\n' \
                              % keyLayers);
                            sys.exit(0);
                        if keyRatios in config:
                            layerRatios = config[keyRatios].strip('()');
                            layerRatios = [float(x.strip()) for x in layerRatios.split(',')];
                        else:
                            sys.stderr.write('\n ## ERROR: %s not found in '  \
                              'config file. Ratios of layer thickness to '    \
                              'total panel thickness must be defined.\n\n'    \
                              % keyRatios);
                            sys.exit(0);                            
                        
                        # Set panel layers and associate each layer with a pre-defined
                        # material
                        matAddress = list();
                        for j in range(len(layerNames)):
                            matAddress.append(nozzle.materials[layerNames[j]]);
                        nozzle.materials[name].setPanelLayers(layerNames,layerRatios,matAddress);   
                                         
                    else:
                        sys.stderr.write('\n ## ERROR : Only ISOTROPIC, '     \
                          'ANISOTROPIC_SHELL, or FIXED_RATIO_PANEL materials are '  \
                          'implemented. Received %s instead.\n\n' % type_prop);
                        sys.exit(0);
                    
                    # Setup material keys
                    m_keys = ['_DENSITY','_ELASTIC_MODULUS','_SHEAR_MODULUS', 
                              '_POISSON_RATIO','_MUTUAL_INFLUENCE_COEFS',
                              '_THERMAL_CONDUCTIVITY',
                              '_THERMAL_EXPANSION_COEF',
                              '_PRINCIPLE_FAILURE_STRAIN',
                              '_LOCAL_FAILURE_STRAIN',
                              '_YIELD_STRESS',
                              '_MAX_SERVICE_TEMPERATURE'];
                    for j in range(len(m_keys)):
                        m_keys[j] = key + m_keys[j];
                    m_keys = tuple(m_keys)
                    
                    # Assign material properties
                    for k in m_keys:
                        
                        if k not in config:
                            if output == 'verbose':
                                sys.stdout.write('%s not found in config ' \
                                'file. Skipping assignment...\n' % k);
                            continue;
                        
                        var = config[k].strip('()');
                        var = var.split(',');
                        if len(var) > 1:
                            var = [float(x) for x in var]
                        elif len(var) == 1:
                            var = float(var[0]);
                        else:
                            sys.stderr.write('\n ## ERROR: Key %s could not be parsed '  \
                            'for %s material\n\n' % (k,name));
                            sys.exit(0);
                        
                        if k == key + '_DENSITY':
                            nozzle.materials[name].setDensity(var);
                        elif k == key + '_ELASTIC_MODULUS':
                            nozzle.materials[name].setElasticModulus(var);
                        elif k == key + '_SHEAR_MODULUS':
                            nozzle.materials[name].setShearModulus(var);
                        elif k == key + '_POISSON_RATIO':
                            nozzle.materials[name].setPoissonRatio(var);
                        elif k == key + '_MUTUAL_INFLUENCE_COEFS':
                            nozzle.materials[name].setMutualInfluenceCoefs(var);
                        elif k == key + '_THERMAL_CONDUCTIVITY':
                            nozzle.materials[name].setThermalConductivity(var);
                        elif k == key + '_THERMAL_EXPANSION_COEF':
                            nozzle.materials[name].setThermalExpansionCoef(var);
                        elif k == key + '_PRINCIPLE_FAILURE_STRAIN':
                            nozzle.materials[name].setFailureType('PRINCIPLE_FAILURE_STRAIN',var);
                        elif k == key + '_LOCAL_FAILURE_STRAIN':
                            nozzle.materials[name].setFailureType('LOCAL_FAILURE_STRAIN',var);
                        elif k == key + '_YIELD_STRESS':
                            nozzle.materials[name].setFailureType('VON_MISES',var);
                        elif k == key + '_MAX_SERVICE_TEMPERATURE':
                            nozzle.materials[name].setMaxServiceTemperature(var);
                        else:
                            sys.stderr.write('\n ## ERROR: Key %s could not implemented' \
                            ' for %s material\n\n' % (k,name));
                            sys.exit(0);
                    
                else: # concise material format
                    sys.stderr.write('\n ## ERROR: Concise material format ' \
                    'not implemented yet\n\n');
                    sys.exit(0);
                
            else:
                break;
            i += 1
        
        if output == 'verbose':
            sys.stdout.write('%d materials processed\n' % (i-1));
            sys.stdout.write('Setup Materials complete\n');

    
    def SetupWallLayers(self, config, output='verbose'):
  
        nozzle = self;
        
        nozzle.wall.layer = list();
        
        # Discover all layers and assign data
        i = 1;
        while i:
            key = 'LAYER' + str(i);
            if key in config:
                # Parse information related to layer
                info = config[key].strip('()');
                name, param, material = [x.strip() for x in info.split(',')];
                nozzle.wall.layer.append(component.AxisymmetricWall(name));
                nozzle.wall.layer[i-1].param = param;
                nozzle.wall.layer[i-1].formalName = key;
                
                # Assign thickness distribution
                self.SetupLayerThickness(config,nozzle.wall.layer[i-1]);
                # Assign material
                nozzle.wall.layer[i-1].material=nozzle.materials[material];
                
            else:
                break;
            i += 1
        
        #nozzle.wall.load_layer    = nozzle.wall.layer[1];
        #nozzle.wall.thermal_layer = nozzle.wall.layer[0];

        if output == 'verbose':
            sys.stdout.write('%d layers processed\n' % (i-1));
            sys.stdout.write('Setup Wall Layers complete\n');        
        #print nozzle.wall.layer[1].__dict__

        
    def SetupBaffles(self, config, output='verbose'):
        
        nozzle = self;
        
        info = config['BAFFLES'].strip('()');
        n, material = [x.strip() for x in info.split(',')];
        nozzle.baffles = component.Baffles(n);
        nozzle.baffles.material = nozzle.materials[material];
        
        info = config['BAFFLES_LOCATION'].strip('()');
        ltemp = [x.strip() for x in info.split(',')];
        nozzle.baffles.location = [float(x) for x in ltemp];
        
        info = config['BAFFLES_HEIGHT'].strip('()');
        ltemp = [x.strip() for x in info.split(',')];        
        nozzle.baffles.height = [float(x) for x in ltemp];
        
        info = config['BAFFLES_THICKNESS'].strip('()');
        ltemp = [x.strip() for x in info.split(',')];          
        nozzle.baffles.thickness = [float(x) for x in ltemp];  
        
        if output == 'verbose':
            sys.stdout.write('Setup Baffles complete\n');        

        
    def SetupStringers(self, config, output='verbose'):
        
        nozzle = self;
        
        info = config['STRINGERS'].strip('()');
        n, material = [x.strip() for x in info.split(',')];
        nozzle.stringers = component.Stringers(n);
        nozzle.stringers.material = nozzle.materials[material];
        
        if ( 'STRINGERS_BREAK_LOCATIONS' not in config or 
             config['STRINGERS_BREAK_LOCATIONS'] == 'BAFFLES_LOCATION' ):
                 if output == 'verbose':
                     sys.stdout.write('Stringer break locations will be set' \
                       ' from baffle locations\n');
                 key = '';
                 k_loc = 'BAFFLES_LOCATION';
                 if ( 'STRINGERS_HEIGHT_VALUES' not in config or
                     config['STRINGERS_HEIGHT_VALUES'] == 'BAFFLES_HEIGHT' ):
                         if output == 'verbose':
                             sys.stdout.write('Stringer heights will be set'  \
                               ' from baffle heights\n');
                         h_val = 'BAFFLES_HEIGHT';
                 else:
                     h_val = 'STRINGERS_HEIGHT_VALUES';
                 t_val = 'STRINGERS_THICKNESS_VALUES';
        else:
            if ( 'STRINGERS_HEIGHT_VALUES' not in config or
                config['STRINGERS_HEIGHT_VALUES'] == 'BAFFLES_HEIGHT' ):
                    sys.stderr.write('\n ## ERROR: Stringer height values '   \
                      'cannot be set to baffle height values if stringers '   \
                      'break locations are not set to baffles location.\n\n')
                    sys.exit(0);
            key = 'STRINGERS';
            k_loc = '_BREAK_LOCATIONS';
            t_val = '_THICKNESS_VALUES';
            h_val = '_HEIGHT_VALUES';
        
        # Assign thickness distribution
        try:
            nozzle.stringers.thicknessNodes = self.ParseThickness(config,key,loc=k_loc,val=t_val);
        except:
            sys.stderr.write('\n ## ERROR : Thickness definition '   \
                 'could not be parsed for STRINGERS.\n\n');
            sys.exit(0);
        
        # Assign height distribution
        try:
            nozzle.stringers.heightNodes = self.ParseThickness(config,key,loc=k_loc,val=h_val);
        except:
            sys.stderr.write('\n ## ERROR : Height definition '   \
                 'could not be parsed for STRINGERS.\n\n');
            sys.exit(0);             

        #info = config['STRINGERS_HEIGHT'].strip('()');
        #ltemp = [x.strip() for x in info.split(',')];        
        #nozzle.stringers.height = [float(x) for x in ltemp];
        
        #info = config['STRINGERS_THICKNESS'].strip('()');
        #ltemp = [x.strip() for x in info.split(',')];          
        #nozzle.stringers.thickness = [float(x) for x in ltemp];
        
        if output == 'verbose':
            sys.stdout.write('Setup Stringers complete\n');        
        

    def SetupWall (self, config, output='verbose'):
        
        nozzle = self;
        
        # --- SHAPE OF INNER WALL (B-SPLINE)
        
        coefs = nozzle.wall.coefs;
        knots = nozzle.wall.knots;
        
        coefs_size = nozzle.wall.coefs_size;
        
        nozzle.height = coefs[coefs_size-1];
        nozzle.length = coefs[coefs_size/2-1];
        
        nozzle.xthrust = -1; # --- x crd for the thrust integration
        nozzle.ythrust = -1;
        
        if nozzle.method == 'RANS' or nozzle.method == 'EULER':
            x    = [];
            y    = [];
        
            nx = 4000; # use 4000 points to interpolate inner wall shape
        
            _meshutils_module.py_BSplineGeo3 (knots, coefs, x, y, nx);
        
            nozzle.xwall = x;
            nozzle.ywall = y;
            
            dx_exit = max(1.3*nozzle.meshhl[3], 0.001);
            for i in range(0,nx) :
            	if (  x[nx-i-1] < x[-1]-dx_exit  ):
            		nozzle.x_thrust = x[nx-i-1];
            		nozzle.y_thrust = y[nx-i-1];
            		break;
            
        coefsnp = np.empty([2, coefs_size/2]);
    
        for i in range(0, coefs_size/2):
            coefsnp[0][i] = coefs[i];
            coefsnp[1][i] = coefs[i+coefs_size/2];        
        
        nozzle.wall.geometry = geometry.Bspline(coefsnp);
        
        # --- Setup inner nozzle wall temperature if necessary
        if hasattr(nozzle.wall,'temperature'):
            if nozzle.wall.temperature.param == 'PIECEWISE_LINEAR':
                tsize = len(nozzle.wall.temperature.thicknessNodes);
                thicknessNodeArray = np.zeros(shape=(2,tsize));
                for j in range(tsize):
                    thicknessNodeArray[0][j] = nozzle.wall.temperature.thicknessNodes[j][0]*nozzle.length;
                    thicknessNodeArray[1][j] = nozzle.wall.temperature.thicknessNodes[j][1];
                nozzle.wall.temperature.geometry = geometry.PiecewiseLinear(thicknessNodeArray);                
            elif nozzle.wall.temperature.param == 'KLE':
                sys.stderr.write('\n ## ERROR: KLE not implemented yet for '  \
                  'inner nozzle wall temperature definition.\n\n');
                sys.exit(0);
        
        # --- Setup thickness of each layer
        
        for i in range(len(nozzle.wall.layer)):
            
            tsize = len(nozzle.wall.layer[i].thicknessNodes);
            thicknessNodeArray = np.zeros(shape=(2,tsize))        
            for j in range(tsize):
                thicknessNodeArray[0][j] = nozzle.wall.layer[i].thicknessNodes[j][0]*nozzle.length;
                thicknessNodeArray[1][j] = nozzle.wall.layer[i].thicknessNodes[j][1];
            nozzle.wall.layer[i].thickness = geometry.PiecewiseLinear(thicknessNodeArray);

        # --- Setup exterior wall
        nozzle.exterior = component.AxisymmetricWall('exterior');
        
        # Determine *approx* height of nozzle exit at outside of outermost layer
        exitHeight = nozzle.wall.geometry.radius(nozzle.length);
        for i in range(len(nozzle.wall.layer)):
            exitHeight += nozzle.wall.layer[i].thickness.radius(nozzle.length);
        shapeDefinition = np.array([[0., 0.1548, nozzle.length],
                                    [0.4244, 0.4244, exitHeight + 0.012]]);
        nozzle.exterior.geometry = geometry.PiecewiseLinear(shapeDefinition);      
          
        # --- Setup thickness distribution for stringers
        tsize = len(nozzle.stringers.thicknessNodes);
        thicknessNodeArray = np.zeros(shape=(2,tsize))
        for j in range(tsize):
            thicknessNodeArray[0][j] = nozzle.stringers.thicknessNodes[j][0]*nozzle.length;
            thicknessNodeArray[1][j] = nozzle.stringers.thicknessNodes[j][1];    
        nozzle.stringers.thickness = geometry.PiecewiseLinear(thicknessNodeArray);
        
        # --- Setup height distribution for stringers
        tsize = len(nozzle.stringers.heightNodes);
        thicknessNodeArray = np.zeros(shape=(2,tsize))
        for j in range(tsize):
            thicknessNodeArray[0][j] = nozzle.stringers.heightNodes[j][0]*nozzle.length;
            thicknessNodeArray[1][j] = nozzle.stringers.heightNodes[j][1];    
        nozzle.stringers.height = geometry.PiecewiseLinear(thicknessNodeArray);        
            
        # --- Rescale x-coordinates of baffles
        nozzle.baffles.location = [q*nozzle.length for q in nozzle.baffles.location];

        # --- Useful x vs. r data for updating baffle and stringer heights
        n = 10000 # 1e4
        x = np.linspace(0,nozzle.length,n)
        rList = geometry.layerCoordinatesInGlobalFrame(nozzle,x);
        
        # --- Setup external geometry
        if nozzle.stringers.n > 0:
            innerBaffleRadius = np.interp(nozzle.length,x,rList[-2]);
        else:
            innerBaffleRadius = np.interp(nozzle.length,x,rList[-1]);   
        shapeDefinition = np.array([[0., 0.1548, nozzle.length],
                                    [0.4244, 0.4244, innerBaffleRadius + 0.012]]);
        nozzle.exterior.geometry = geometry.PiecewiseLinear(shapeDefinition);  
        
        # --- Check for intersection of nozzle wall with external geometry
        wallDiff = rList[-2] - nozzle.exterior.geometry.radius(x);
        if max(wallDiff) > 0:
            sys.stdout.write('\nWARNING! Nozzle wall intersects with exterior geometry\n');
            wallBend = 0.1548;
            for i in range(len(nozzle.wall.layer)):
                wallBend = max([wallBend, nozzle.wall.layer[i].thickness.nodes[0,-2]]);
            shapeDefinition = np.array([[0., wallBend, nozzle.length],
                                        [0.4244, 0.4244, innerBaffleRadius + 0.012]]);
            nozzle.exterior.geometry = geometry.PiecewiseLinear(shapeDefinition); 
            sys.stdout.write('WARNING! Altered exterior nozzle geometry to fit nozzle wall shape\n\n');
       
        # --- Update height of baffles to coincide with exterior wall shape
        for i in range(len(nozzle.baffles.location)):
            loc = nozzle.baffles.location[i];
            # inner baffle radius is outside radius of outermost layer
            if nozzle.stringers.n > 0:
                innerBaffleRadius = np.interp(loc,x,rList[-2]);
            else:
                innerBaffleRadius = np.interp(loc,x,rList[-1]);
            outerBaffleRadius = nozzle.exterior.geometry.radius(loc);
            nozzle.baffles.height[i] = outerBaffleRadius - innerBaffleRadius;
            if output == 'verbose':
                sys.stdout.write('Baffle %i height resized to %f\n' % 
                  (i+1,nozzle.baffles.height[i]));
            
        # --- If stringers are dependent on baffles, re-update stringers
        if ( 'STRINGERS_HEIGHT_VALUES' not in config or
            config['STRINGERS_HEIGHT_VALUES'] == 'BAFFLES_HEIGHT' ):  
                for i in range(len(nozzle.stringers.heightNodes)):
                    nozzle.stringers.heightNodes[i][1] = nozzle.baffles.height[i];
       
        # --- Setup height distribution for stringers
        tsize = len(nozzle.stringers.heightNodes);
        thicknessNodeArray = np.zeros(shape=(2,tsize))
        for j in range(tsize):
            thicknessNodeArray[0][j] = nozzle.stringers.heightNodes[j][0]*nozzle.length;
            thicknessNodeArray[1][j] = nozzle.stringers.heightNodes[j][1];    
        nozzle.stringers.height = geometry.PiecewiseLinear(thicknessNodeArray);        
        
        # --- Check that stringer heights are not above baffle heights
        for i in range(len(nozzle.baffles.location)):
            baffleLocation = nozzle.baffles.location[i];
            baffleHeight = nozzle.baffles.height[i];
            stringerHeight = nozzle.stringers.height.radius(baffleLocation);
            if stringerHeight > baffleHeight:
                sys.stderr.write('\n ## ERROR: Stringers must have heights'  \
                  ' that remain below the baffle height.\n\n');
                sys.exit(0);

        if output == 'verbose':
            sys.stdout.write('Setup Wall complete\n');


    def ParseDV (self, config, output='verbose'):
        
        nozzle = self;
        
        if 'INPUT_DV_FORMAT' in config:
            inputDVformat = config['INPUT_DV_FORMAT'];
        else :
            sys.stderr.write('\n ## ERROR : Input DV file format not '        \
              'specified. (INPUT_DV_FORMAT expected: PLAIN or DAKOTA)\n\n');
            sys.exit(0);        
                
        if 'INPUT_DV_NAME' in config:
            filename = config['INPUT_DV_NAME'];
        else :
            sys.stderr.write('\n ## ERROR : Input DV file name not '          \
              'specified. (INPUT_DV_NAME expected)\n\n');
            sys.exit(0);
                
        if inputDVformat == 'PLAIN':
            DV_List, OutputCode, Derivatives_DV = ParseDesignVariables_Plain(filename);    
            NbrDV = len(DV_List);                    
        elif inputDVformat == 'DAKOTA' :
            DV_List, OutputCode, Derivatives_DV = ParseDesignVariables_Dakota(filename);    
            NbrDV = len(DV_List);
        else:
            sys.stderr.write('\n ## ERROR : Unknown DV input file format '    \
              '%s\n\n' % inputDVformat);
            sys.exit(0);
        
        if NbrDV != nozzle.NbrDVTot : 
            sys.stderr.write('\n ## Error : Inconsistent number of design '   \
              'variables are given in %s\n\n' % filename);
            sys.stderr.write('             %d given, %d expected\n' %         \
              (NbrDV, nozzle.NbrDVTot ));
            sys.exit(0);
    
        nozzle.DV_List = DV_List;
        
        if output == 'verbose':
            sys.stdout.write('Setup Parse Design Variables complete\n');        
    
    
    def UpdateDV(self, config, output='verbose'):
        
        nozzle = self;
        
        NbrTags = len(nozzle.DV_Tags);
        
        prt_name = [];
        prt_basval = [];
        prt_newval = [];
        
        NbrChanged = 0; # Count the total number of changed parameters
                        # Note: different from the number of DV, 
                        #   because one DV might correspond to more BSP coefs
        
        for iTag in range(NbrTags):
            Tag = nozzle.DV_Tags[iTag];
            NbrDV = nozzle.DV_Head[iTag+1] - nozzle.DV_Head[iTag];
            
            if Tag == 'WALL':
                for iCoef in range(len(nozzle.wall.dv)):
                    id_dv = nozzle.DV_Head[iTag] + nozzle.wall.dv[iCoef];                    
                    # --- Update coef iCoef if required
                    if id_dv >= nozzle.DV_Head[iTag]:
                        prt_name.append('Bspline coef #%d' % (iCoef+1));
                        prt_basval.append('%.4lf'% nozzle.wall.coefs[iCoef]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.wall.coefs[iCoef] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                continue;
            
            # Update general piecewise linear wall temp definition
            if Tag == 'WALL_TEMP':
                lsize = len(nozzle.wall.temperature.thicknessNodes);
                brk = np.max(nozzle.wall.temperature.dv[:lsize])+1;
                for iCoord in range(len(nozzle.wall.temperature.dv)):
                    id_dv = nozzle.DV_Head[iTag] + nozzle.wall.temperature.dv[iCoord];
                    # Update coordinate in thickness array if required
                    if id_dv < nozzle.DV_Head[iTag]:
                        pass
                    elif id_dv < nozzle.DV_Head[iTag]+brk:
                        prt_name.append('wall temp location #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodes[iCoord][0]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.wall.temperature.thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                    else: # id_dv > nozzle.DV_Head[iTag]+brk
                        prt_name.append('wall temp value #%d' % (iCoord+1-lsize));
                        prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodes[iCoord-lsize][1]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.wall.temperature.thicknessNodes[iCoord-lsize][1] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                continue;
                
            # Update specific thickness or values for piecewise linear wall temp def
            check = 0;
            if Tag == 'WALL_TEMP_LOCATIONS':
                for iCoord in range(len(nozzle.wall.temperature.thicknessNodes)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('wall temp location #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodes[iCoord][0]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.wall.temperature.thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                check = 1;
            elif Tag == 'WALL_TEMP_VALUES':
                for iCoord in range(len(nozzle.wall.temperature.thicknessNodes)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                    prt_name.append('wall temp value #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodes[iCoord][1]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.wall.temperature.thicknessNodes[iCoord][1] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                check = 1;
            if check == 1:
                continue;                    

            # Update all layers with non-specific names, e.g. LAYER1, etc.
            check = 0;
            for j in range(len(nozzle.wall.layer)):
                if Tag == nozzle.wall.layer[j].name:
                    lsize = len(nozzle.wall.layer[j].thicknessNodes);
                    brk = np.max(nozzle.wall.layer[j].dv[:lsize])+1;
                    for iCoord in range(len(nozzle.wall.layer[j].dv)):
                        id_dv = nozzle.DV_Head[iTag] + nozzle.wall.layer[j].dv[iCoord];
                        # Update coordinate in thickness array if required
                        if id_dv < nozzle.DV_Head[iTag]:
                            pass
                        elif id_dv < nozzle.DV_Head[iTag]+brk:
                            prt_name.append('%s thickness location #%d' % (nozzle.wall.layer[j].name,iCoord+1));
                            prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodes[iCoord][0]);
                            prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                            nozzle.wall.layer[j].thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                            NbrChanged = NbrChanged+1;
                        else: # id_dv > nozzle.DV_Head[iTag]+brk
                            prt_name.append('%s thickness value #%d' % (nozzle.wall.layer[j].name,iCoord+1-lsize));
                            prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodes[iCoord-lsize][1]);
                            prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                            nozzle.wall.layer[j].thicknessNodes[iCoord-lsize][1] = nozzle.DV_List[id_dv];
                            NbrChanged = NbrChanged+1;
                    check = 1;
            if check == 1:
                continue;
                
            # Update all piecewise-linear layers with specific names, e.g. LAYER1_THICKNESS_VALUES
            check = 0;
            for j in range(len(nozzle.wall.layer)):
                if Tag == nozzle.wall.layer[j].name + '_THICKNESS_LOCATIONS':
                    for iCoord in range(len(nozzle.wall.layer[j].thicknessNodes)):
                        id_dv = nozzle.DV_Head[iTag] + iCoord;
                        prt_name.append('%s thickness location #%d' % (nozzle.wall.layer[j].name,iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodes[iCoord][0]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.wall.layer[j].thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                    check = 1;
                elif Tag == nozzle.wall.layer[j].name + '_THICKNESS_VALUES':
                    for iCoord in range(len(nozzle.wall.layer[j].thicknessNodes)):
                        id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                        prt_name.append('%s thickness value #%d' % (nozzle.wall.layer[j].name,iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodes[iCoord][1]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.wall.layer[j].thicknessNodes[iCoord][1] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                    check = 1;
            if check == 1:
                continue;    
                
            # Update all constant layers with specific names, e.g. LAYER1_THICKNESS
            check = 0;
            for j in range(len(nozzle.wall.layer)):
                if Tag == nozzle.wall.layer[j].name + '_THICKNESS':
                    id_dv = nozzle.DV_Head[iTag];
                    prt_name.append('%s thickness' % nozzle.wall.layer[j].name);
                    prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodes[0][1]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.wall.layer[j].thicknessNodes[0][1] = nozzle.DV_List[id_dv];
                    nozzle.wall.layer[j].thicknessNodes[1][1] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                    check = 1;
            if check == 1:
                continue;

            if Tag == 'BAFFLES':
                lsize = len(nozzle.baffles.dv)/3;
                brk1 = np.max(nozzle.baffles.dv[:lsize])+1;
                brk2 = np.max(nozzle.baffles.dv[lsize:2*lsize])+1;
                for iCoord in range(len(nozzle.baffles.dv)):
                    id_dv = nozzle.DV_Head[iTag] + nozzle.baffles.dv[iCoord];                    
                    # --- Update coordinate in baffle definition if required
                    if id_dv < nozzle.DV_Head[iTag]:
                        pass;
                    elif id_dv < nozzle.DV_Head[iTag] + brk1:
                        prt_name.append('baffle location #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.baffles.location[iCoord]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.baffles.location[iCoord] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                        # If stringer location dependent on baffle location
                        if ('STRINGERS_BREAK_LOCATIONS' not in config or 
                            config['STRINGERS_BREAK_LOCATIONS'] == 'BAFFLES_LOCATION'):
                                prt_name.append('stringer location #%d' % (iCoord+1));
                                prt_basval.append('%.4lf'% nozzle.baffles.location[iCoord]);
                                prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                                nozzle.stringers.thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                                nozzle.stringers.heightNodes[iCoord][0] = nozzle.DV_List[id_dv];                        
                                NbrChanged = NbrChanged+1;
                    elif id_dv < nozzle.DV_Head[iTag] + brk2:
                        prt_name.append('baffle thickness #%d' % (iCoord+1-lsize));
                        prt_basval.append('%.4lf'% nozzle.baffles.thickness[iCoord-lsize]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.baffles.thickness[iCoord-lsize] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;                        
                    else:
                        prt_name.append('baffle height #%d' % (iCoord+1-2*lsize));
                        prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord-2*lsize]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.baffles.height[iCoord-2*lsize] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;   
                        # If stringer height depends on baffle height
                        if ('STRINGERS_HEIGHT_VALUES' not in config or 
                            config['STRINGERS_HEIGHT_VALUES'] == 'BAFFLES_HEIGHT'):
                                prt_name.append('stringer height #%d' % (iCoord+1-2*lsize));
                                prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord-2*lsize]);
                                prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                                nozzle.stringers.heightNodes[iCoord-2*lsize][1] = nozzle.DV_List[id_dv];                     
                                NbrChanged = NbrChanged+1;                        
                continue;
                
            if Tag == 'BAFFLES_LOCATION':
                for iCoord in range(len(nozzle.baffles.location)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('baffle location #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.baffles.location[iCoord]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.baffles.location[iCoord] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                    # If stringer location dependent on baffle location
                    if ('STRINGERS_BREAK_LOCATIONS' not in config or 
                        config['STRINGERS_BREAK_LOCATIONS'] == 'BAFFLES_LOCATION'):
                            prt_name.append('stringer location #%d' % (iCoord+1));
                            prt_basval.append('%.4lf'% nozzle.baffles.location[iCoord]);
                            prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                            nozzle.stringers.thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                            nozzle.stringers.heightNodes[iCoord][0] = nozzle.DV_List[id_dv];                        
                            NbrChanged = NbrChanged+1;                    
                continue;
               
            if Tag == 'BAFFLES_HEIGHT':
                for iCoord in range(len(nozzle.baffles.height)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('baffle height #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.baffles.height[iCoord] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                    # If stringer height depends on baffle height
                    if ('STRINGERS_HEIGHT_VALUES' not in config or 
                        config['STRINGERS_HEIGHT_VALUES'] == 'BAFFLES_HEIGHT'):
                            prt_name.append('stringer height #%d' % (iCoord+1));
                            prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord]);
                            prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                            nozzle.stringers.heightNodes[iCoord][1] = nozzle.DV_List[id_dv];                     
                            NbrChanged = NbrChanged+1;
                continue;
            
            if Tag == 'BAFFLES_THICKNESS':
                for iCoord in range(len(nozzle.baffles.thickness)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('baffle thickness #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.baffles.thickness[iCoord]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.baffles.thickness[iCoord] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                continue;
                
            if Tag == 'STRINGERS':
                lsize = len(nozzle.stringers.dv)/3;
                brk1 = np.max(nozzle.stringers.dv[:lsize])+1;
                brk2 = np.max(nozzle.stringers.dv[lsize:2*lsize])+1;
                for iCoord in range(len(nozzle.stringers.dv)):
                    id_dv = nozzle.DV_Head[iTag] + nozzle.stringers.dv[iCoord];  
                    # Update coordinate in thickness array if required
                    if id_dv < nozzle.DV_Head[iTag]:
                        pass
                    elif id_dv < nozzle.DV_Head[iTag]+brk1:
                        prt_name.append('stringer break location #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodes[iCoord][0]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.stringers.thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                        nozzle.stringers.heightNodes[iCoord][0] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                    elif id_dv < nozzle.DV_Head[iTag] + brk2:
                        prt_name.append('stringer height value #%d' % (iCoord+1-lsize));
                        prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodes[iCoord-lsize][1]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.stringers.heightNodes[iCoord-lsize][1] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                    else: # id_dv > nozzle.DV_Head[iTag]+brk
                        prt_name.append('stringer thickness value #%d' % (iCoord+1-2*lsize));
                        prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodes[iCoord-2*lsize][1]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.stringers.thicknessNodes[iCoord-2*lsize][1] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                continue;    
                
            # Update all layers with specific names, e.g. LAYER1_THICKNESS_VALUES
            if Tag == 'STRINGERS_BREAK_LOCATIONS':
                for iCoord in range(len(nozzle.stringers.thicknessNodes)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('stringer break location #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodes[iCoord][0]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.stringers.thicknessNodes[iCoord][0] = nozzle.DV_List[id_dv];
                    nozzle.stringers.heightNodes[iCoord][0] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                    continue;
            
            if Tag == 'STRINGERS_HEIGHT_VALUES':
                for iCoord in range(len(nozzle.stringers.heightNodes)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                    prt_name.append('stringer height value #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.stringers.heightNodes[iCoord][1]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.stringers.heightNodes[iCoord][1] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                    continue;
                    
            if Tag == 'STRINGERS_THICKNESS_VALUES':
                for iCoord in range(len(nozzle.stringers.thicknessNodes)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                    prt_name.append('stringer thickness value #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodes[iCoord][1]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    nozzle.stringers.thicknessNodes[iCoord][1] = nozzle.DV_List[id_dv];
                    NbrChanged = NbrChanged+1;
                    continue;                    
                        
            # Update all materials with non-specific names, e.g. MATERIAL1, etc.
            check = 0;
            for k in nozzle.materials:
                id_dv = nozzle.DV_Head[iTag];
                
                # Skip fixed_ratio_panel materials which use predefined materials
                #if nozzle.materials[k].type == 'FIXED_RATIO_PANEL':
                #    continue;
                
                if Tag == k: # Update material with non-specific names                
                    if NbrDV == 2: # update density and thermal conductivity
                    
                        prt_name.append('%s density' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getDensity());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setDensity(nozzle.DV_List[id_dv]);

                        prt_name.append('%s thermal conductivity' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalConductivity());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                        nozzle.materials[k].setThermalConductivity(nozzle.DV_List[id_dv+1]);
                        
                        NbrChanged = NbrChanged+2;
                        
                    elif NbrDV == 7:
                        
                        prt_name.append('%s density' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getDensity());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setDensity(nozzle.DV_List[id_dv]);
                        
                        prt_name.append('%s elastic modulus' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getElasticModulus());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                        nozzle.materials[k].setElasticModulus(nozzle.DV_List[id_dv+1]);

                        prt_name.append('%s poisson ratio' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getPoissonRatio());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+2]);
                        nozzle.materials[k].setPoissonRatio(nozzle.DV_List[id_dv+2]);

                        prt_name.append('%s thermal conductivity' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalConductivity());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+3]);
                        nozzle.materials[k].setThermalConductivity(nozzle.DV_List[id_dv+3]);

                        prt_name.append('%s thermal expansion coef' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalExpansionCoef());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+4]);
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.DV_List[id_dv+4]);
                        
                        prt_name.append('%s failure limit' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getFailureLimit());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+5]);
                        nozzle.materials[k].setFailureType(nozzle.materials[k].failureType,nozzle.DV_List[id_dv+5]);

                        prt_name.append('%s max service temperature' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].Tmax);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+6]);
                        nozzle.materials[k].setMaxServiceTemperature(nozzle.DV_List[id_dv+6]);                        
                        
                        NbrChanged = NbrChanged+7;
                        
                    elif NbrDV == 15: # XXX This needs to be updated
                        
                        prt_name.append('%s density' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getDensity());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setDensity(nozzle.DV_List[id_dv]);

                        ltemp = nozzle.materials[k].getElasticModulus();
                        prt_name.append('%s elastic modulus E1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                        prt_name.append('%s elastic modulus E2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+2]);                        
                        nozzle.materials[k].setElasticModulus(nozzle.DV_List[id_dv+1:id_dv+3]);
                        
                        prt_name.append('%s shear modulus' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getShearModulus());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+3]);
                        nozzle.materials[k].setShearModulus(nozzle.DV_List[id_dv+3]);
                        
                        prt_name.append('%s poisson ratio' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getPoissonRatio());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+4]);
                        nozzle.materials[k].setPoissonRatio(nozzle.DV_List[id_dv+4]);
                        
                        ltemp = nozzle.materials[k].getMutualInfluenceCoefs();
                        prt_name.append('%s mutual influence coef u1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+5]);
                        prt_name.append('%s mutual influence coef u2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+6]);                        
                        nozzle.materials[k].setMutualInfluenceCoefs(nozzle.DV_List[id_dv+5:id_dv+7]);       
                        
                        ltemp = nozzle.materials[k].getThermalConductivity();
                        prt_name.append('%s thermal conductivity k1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+7]);
                        prt_name.append('%s thermal conductivity k2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+8]);                        
                        nozzle.materials[k].setThermalConductivity(nozzle.DV_List[id_dv+7:id_dv+9]);  

                        ltemp = nozzle.materials[k].getThermalExpansionCoef();
                        prt_name.append('%s thermal expansion coef a1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+9]);
                        prt_name.append('%s thermal expansion coef a2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+10]);           
                        prt_name.append('%s thermal expansion coef a12' % k);
                        prt_basval.append('%.2le' % ltemp[2]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+11]);                         
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.DV_List[id_dv+9:id_dv+12]);
                        
                        NbrChanged = NbrChanged+12;

                    else:
                        raise RuntimeError('%d design variables not accepted' \
                              ' for assignment for %s material' % (NbrDV,k));
                    
                    check = 1;
                              
                elif Tag == k + '_DENSITY':
                    prt_name.append('%s density' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].getDensity());
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    nozzle.materials[k].setDensity(nozzle.DV_List[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;
                elif Tag == k + '_ELASTIC_MODULUS':
                    if NbrDV == 1:
                        prt_name.append('%s elastic modulus' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getElasticModulus());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setElasticModulus(nozzle.DV_List[id_dv]);
                        NbrChanged = NbrChanged + 1;
                    elif NbrDV == 2:
                        ltemp = nozzle.materials[k].getElasticModulus();
                        prt_name.append('%s elastic modulus E1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s elastic modulus E2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);                        
                        nozzle.materials[k].setElasticModulus(nozzle.DV_List[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                    check = 1;                        
                elif Tag == k + '_SHEAR_MODULUS':
                    prt_name.append('%s shear modulus' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].getShearModulus());
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    nozzle.materials[k].setShearModulus(nozzle.DV_List[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;
                elif Tag == k + '_POISSON_RATIO':
                    prt_name.append('%s poisson ratio' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].getPoissonRatio());
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    nozzle.materials[k].setPoissonRatio(nozzle.DV_List[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;
                elif Tag == k + '_MUTUAL_INFLUENCE_COEFS':
                    ltemp = nozzle.materials[k].getMutualInfluenceCoefs();
                    prt_name.append('%s mutual influence coef u1' % k);
                    prt_basval.append('%.2le' % ltemp[0]);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    prt_name.append('%s mutual influence coef u2' % k);
                    prt_basval.append('%.2le' % ltemp[1]);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]); 
                    if NbrDV == 1:
                        nozzle.materials[k].setMutualInfluenceCoefs(nozzle.DV_List[id_dv]);  
                    elif NbrDV == 2:                                               
                        nozzle.materials[k].setMutualInfluenceCoefs(nozzle.DV_List[id_dv:id_dv+2]);  
                    NbrChanged = NbrChanged + 2;
                    check = 1;
                elif Tag == k + '_THERMAL_CONDUCTIVITY':
                    if NbrDV == 1:
                        prt_name.append('%s thermal conductivity' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalConductivity());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setThermalConductivity(nozzle.DV_List[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        ltemp = nozzle.materials[k].getThermalConductivity();
                        prt_name.append('%s thermal conductivity k1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s thermal conductivity k2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);                        
                        nozzle.materials[k].setThermalConductivity(nozzle.DV_List[id_dv:id_dv+2]);    
                        NbrChanged = NbrChanged + 2;
                        check = 1;
                    elif NbrDV == 3:
                        ltemp = nozzle.materials[k].getThermalConductivity();
                        prt_name.append('%s thermal conductivity k1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s thermal conductivity k2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);          
                        prt_name.append('%s thermal conductivity k3' % k);
                        prt_basval.append('%.2le' % ltemp[2]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+2]);                                                
                        nozzle.materials[k].setThermalConductivity(nozzle.DV_List[id_dv:id_dv+3]);    
                        NbrChanged = NbrChanged + 3;
                        check = 1;                       
                elif Tag == k + '_THERMAL_EXPANSION_COEF':
                    if NbrDV == 1:
                        prt_name.append('%s thermal expansion coef' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalExpansionCoef());
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.DV_List[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 3:
                        ltemp = nozzle.materials[k].getThermalExpansionCoef();
                        prt_name.append('%s thermal expansion coef a1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s thermal expansion coef a2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);           
                        prt_name.append('%s thermal expansion coef a12' % k);
                        prt_basval.append('%.2le' % ltemp[2]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+2]);                         
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.DV_List[id_dv:id_dv+3]);
                        NbrChanged = NbrChanged + 3;
                        check = 1;
                elif Tag == k + '_PRINCIPLE_FAILURE_STRAIN':
                    if NbrDV == 1:
                        prt_name.append('%s principal failure strain' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setFailureType('PRINCIPLE_FAILURE_STRAIN',nozzle.DV_List[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        prt_name.append('%s principal failure strain (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s principal failure strain (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);                        
                        nozzle.materials[k].setFailureType('PRINCIPLE_FAILURE_STRAIN',nozzle.DV_List[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                        check = 1;                        
                    elif NbrDV == 5:
                        prt_name.append('%s principal failure strain e1 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s principal failure strain e1 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);   
                        prt_name.append('%s principal failure strain e2 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[2]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+2]);
                        prt_name.append('%s principal failure strain e2 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[3]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+3]);   
                        prt_name.append('%s principal failure strain e12' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[4]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+4]);                          
                        nozzle.materials[k].setFailureType('PRINCIPLE_FAILURE_STRAIN',nozzle.DV_List[id_dv:id_dv+5]);
                        NbrChanged = NbrChanged + 5;
                        check = 1; 
                elif Tag == k + '_LOCAL_FAILURE_STRAIN':
                    if NbrDV == 1:
                        prt_name.append('%s local failure strain' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setFailureType('LOCAL_FAILURE_STRAIN',nozzle.DV_List[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        prt_name.append('%s local failure strain (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s local failure strain (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);                        
                        nozzle.materials[k].setFailureType('LOCAL_FAILURE_STRAIN',nozzle.DV_List[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                        check = 1;                        
                    elif NbrDV == 5:
                        prt_name.append('%s local failure strain e1 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s local failure strain e1 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);   
                        prt_name.append('%s local failure strain e2 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[2]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+2]);
                        prt_name.append('%s local failure strain e2 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[3]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+3]);   
                        prt_name.append('%s local failure strain e12' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[4]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+4]);                          
                        nozzle.materials[k].setFailureType('LOCAL_FAILURE_STRAIN',nozzle.DV_List[id_dv:id_dv+5]);
                        NbrChanged = NbrChanged + 5;
                        check = 1; 
                elif Tag == k + '_YIELD_STRESS':
                    if NbrDV == 1:
                        prt_name.append('%s yield stress' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].yieldStress);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        nozzle.materials[k].setFailureType('VON_MISES',nozzle.DV_List[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        prt_name.append('%s yield stress direction 1' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].yieldStress[0]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                        prt_name.append('%s yield stress direction 2' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].yieldStress[1]);
                        prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);                        
                        nozzle.materials[k].setFailureType('VON_MISES',nozzle.DV_List[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                        check = 1;                
                elif Tag == k + '_MAX_SERVICE_TEMPERATURE':
                    prt_name.append('%s max service temp' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].Tmax);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    nozzle.materials[k].setMaxServiceTemperature(nozzle.DV_List[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;            
            if check == 1:
                continue;

            if Tag == 'INLET_PSTAG':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Inlet stagnation pressure');
                prt_basval.append('%.2lf'% nozzle.inlet.Pstag);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);                
                nozzle.inlet.setPstag(nozzle.DV_List[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;
                
            if Tag == 'INLET_TSTAG':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Inlet stagnation temperature');
                prt_basval.append('%.2lf'% nozzle.inlet.Tstag);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);                
                nozzle.inlet.setTstag(nozzle.DV_List[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;

            if Tag == 'ATM_PRES':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Atmospheric pressure');
                prt_basval.append('%.2lf'% nozzle.environment.P);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);                
                nozzle.environment.setPressure(nozzle.DV_List[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;
                
            if Tag == 'ATM_TEMP':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Atmospheric temperature');
                prt_basval.append('%.2lf'% nozzle.environment.T);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);                
                nozzle.environment.setTemperature(nozzle.DV_List[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;

            if Tag == 'HEAT_XFER_COEF_TO_ENV':
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Heat xfer coef. to environment');
                prt_basval.append('%.2lf'% nozzle.environment.hInf);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);                
                nozzle.environment.setHeatTransferCoefficient(nozzle.DV_List[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;
        
        # --- Print summary
        if output == 'verbose':
          sys.stdout.write('\n%d parameter(s) updated according to %d design variable(s). Summary:\n' % (NbrChanged, nozzle.NbrDVTot));
      
          sys.stdout.write('-' * 79);
          sys.stdout.write('\n%s | %s | %s\n' % ("DV name".ljust(45), "Baseline value".ljust(20),"Updated value".ljust(20)));
          sys.stdout.write('-' * 79);
          sys.stdout.write('\n');
          for i in range(0,len(prt_name)):
              sys.stdout.write('%s | %s | %s\n' % (prt_name[i].ljust(45), prt_basval[i].ljust(20),prt_newval[i].ljust(20)));
          sys.stdout.write('-' * 79);    
          sys.stdout.write('\n\n');
        elif output == 'quiet':
          pass
        else:
          raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')
          
        if output == 'verbose':
            sys.stdout.write('Setup Update Design Variables complete\n');          
        
    def SetupOutputFunctions (self, config, output='verbose'):
        
        nozzle = self;
        
        nozzle.Output_Tags = [];
        
        nozzle.GetOutput = dict();
        
        nozzle.OutputLocations = dict(); # for field variable outputs
        
        if 'OUTPUT_NAME' in config:
            nozzle.Output_Name = config['OUTPUT_NAME'];
        else :
            sys.stderr.write("\n ## ERROR : Output function file name not "   \
              "specified in config file. (OUTPUT_NAME expected)\n\n");
            sys.exit(0);
            
        if 'OUTPUT_FORMAT' in config:
            if config['OUTPUT_FORMAT'] == 'PLAIN':
                nozzle.outputFormat = 'PLAIN';
            elif config['OUTPUT_FORMAT'] == 'DAKOTA':
                nozzle.outputFormat = 'DAKOTA';
            else:
                sys.stderr.write("\n ## ERROR : Output function file format " \
                  "can only be PLAIN or DAKOTA\n\n");
                sys.exit(0);                    
        else :
            nozzle.outputFormat = 'PLAIN';
            
        # Determine which components have stresses and/or temperatures calculated
        nozzle.stressComponentList = list();
        nozzle.tempComponentList = list();
        # Add elements to stress component list
        for i in range(len(nozzle.wall.layer)):
            nozzle.stressComponentList.append(nozzle.wall.layer[i].name);
        nozzle.stressComponentList.append('STRINGERS');
        for i in range(nozzle.baffles.n):
            nozzle.stressComponentList.append('BAFFLE' + str(i+1));
        # Add elements to temperature component list
        for i in range(len(nozzle.wall.layer)):
            nozzle.tempComponentList.append(nozzle.wall.layer[i].name);
        
        # --- Initialize outputs
        nozzle.mass = -1;
        nozzle.mass_wall_only = -1;
        nozzle.volume = -1;
        nozzle.thrust = -1;
        nozzle.ks_total_stress = list();
        nozzle.pn_total_stress = list();
        nozzle.max_total_stress = list();
        nozzle.max_thermal_stress = list();
        nozzle.max_mechanical_stress = list();
        nozzle.max_failure_criteria = list();
        nozzle.ks_failure_criteria = list();
        nozzle.pn_failure_criteria = list();
        for i in range(len(nozzle.stressComponentList)):
            nozzle.ks_total_stress.append(-1);
            nozzle.pn_total_stress.append(-1);
            nozzle.max_total_stress.append(-1);
            nozzle.max_thermal_stress.append(-1);
            nozzle.max_mechanical_stress.append(-1);
            nozzle.max_failure_criteria.append(-1);
            nozzle.ks_failure_criteria.append(-1);
            nozzle.pn_failure_criteria.append(-1);
        nozzle.ks_temperature = list();
        nozzle.pn_temperature = list();
        nozzle.max_temperature = list();
        nozzle.ks_temp_ratio = list();
        nozzle.pn_temp_ratio = list();
        nozzle.max_temp_ratio = list();
        for i in range(len(nozzle.tempComponentList)):
            nozzle.ks_temperature.append(-1);
            nozzle.pn_temperature.append(-1);
            nozzle.max_temperature.append(-1); 
            nozzle.ks_temp_ratio.append(-1);
            nozzle.pn_temp_ratio.append(-1);
            nozzle.max_temp_ratio.append(-1);
        nozzle.wall_temperature = -1;
        nozzle.wall_pressure = -1;
        nozzle.pressure = -1;
        nozzle.velocity = -1;
        
        if 'OUTPUT_FUNCTIONS' in config:
            
            # Possible output functions
            dv_scalar = ['MASS', 'MASS_WALL_ONLY', 'VOLUME', 'THRUST'];
            dv_stress = ['KS_TOTAL_STRESS', 'PN_TOTAL_STRESS', 
                         'MAX_TOTAL_STRESS', 'MAX_MECHANICAL_STRESS',
                         'MAX_THERMAL_STRESS'];
            dv_failure = ['MAX_FAILURE_CRITERIA', 'KS_FAILURE_CRITERIA',
                          'PN_FAILURE_CRITERIA'];
            dv_temp = ['KS_TEMPERATURE','PN_TEMPERATURE','MAX_TEMPERATURE'];
            dv_tempRatio = ['KS_TEMP_RATIO','PN_TEMP_RATIO','MAX_TEMP_RATIO'];
            dv_field = ['WALL_TEMPERATURE','WALL_PRESSURE','PRESSURE','VELOCITY'];
            dv_keys = dv_scalar + dv_stress + dv_failure + dv_temp + \
                      dv_tempRatio + dv_field;
            
            # Stress and temperature outputs can have prefixes appended to them
            stress_prefix = [];
            temp_prefix = [];
            for i in range(len(nozzle.wall.layer)):
                stress_prefix.append(nozzle.wall.layer[i].name);
                temp_prefix.append(nozzle.wall.layer[i].name);
            if nozzle.stringers.n > 0:
                stress_prefix.append('STRINGERS');
            for i in range(nozzle.baffles.n):
                stress_prefix.append('BAFFLE' + str(i+1));
                
            for key in dv_keys:
                if key in dv_scalar:
                    nozzle.GetOutput[key] = 0;
                elif key in dv_stress:
                    nozzle.GetOutput[key] = [0]*len(stress_prefix);
                elif key in dv_failure:
                    nozzle.GetOutput[key] = [0]*len(stress_prefix);
                elif key in dv_temp:
                    nozzle.GetOutput[key] = [0]*len(temp_prefix);
                elif key in dv_tempRatio:
                    nozzle.GetOutput[key] = [0]*len(temp_prefix);
                elif key in dv_field:
                    nozzle.GetOutput[key] = 0;

            hdl = config['OUTPUT_FUNCTIONS'].strip('()');
            hdl = hdl.split(",");
            
            for i in range(len(hdl)):
                
                key = hdl[i].strip();
                        
                if( key in dv_scalar ):
                    nozzle.GetOutput[key] = 1;
                    nozzle.Output_Tags.append(key);                                  
                elif( key == 'MAX_THERMAL_STRESS' ):
                    if nozzle.structuralFlag == 1 and nozzle.thermalFlag == 1:
                        nozzle.GetOutput[key] = [1]*len(stress_prefix);
                        nozzle.Output_Tags.append(key);
                    else:
                        sys.stderr.write('\n ## ERROR : %s cannot be '        \
                          'returned if THERMAL_ANALYSIS= 0 or '               \
                          'STRUCTURAL_ANALYSIS= 0.\n\n' % key);
                        sys.exit(0);  
                elif( key in dv_stress ):
                    if nozzle.structuralFlag == 1:
                        nozzle.GetOutput[key] = [1]*len(stress_prefix);
                        nozzle.Output_Tags.append(key);
                    else:
                        sys.stderr.write('\n ## ERROR : %s cannot be '        \
                          'returned if STRUCTURAL_ANALYSIS= 0.\n\n' % key);
                        sys.exit(0);        
                elif( key in dv_failure ):
                    if nozzle.structuralFlag == 1:
                        nozzle.GetOutput[key] = [1]*len(stress_prefix);
                        nozzle.Output_Tags.append(key);
                    else:
                        sys.stderr.write('\n ## ERROR : %s cannot be '        \
                          'returned if STRUCTURAL_ANALYSIS= 0.\n\n' % key);
                        sys.exit(0);                          
                elif( key in dv_temp):
                    if nozzle.thermalFlag == 1:
                        nozzle.GetOutput[key] = [1]*len(temp_prefix);
                        nozzle.Output_Tags.append(key);
                    else:
                        sys.stderr.write('\n ## ERROR : %s cannot be '        \
                          'returned if THERMAL_ANALYSIS= 0\n\n' % key);  
                elif( key in dv_tempRatio):
                    if nozzle.thermalFlag == 1:
                        nozzle.GetOutput[key] = [1]*len(temp_prefix);
                        nozzle.Output_Tags.append(key);
                    else:
                        sys.stderr.write('\n ## ERROR : %s cannot be '        \
                          'returned if THERMAL_ANALYSIS= 0\n\n' % key);   
                elif( key in dv_field ):
                    nozzle.GetOutput[key] = 1;
                    nozzle.Output_Tags.append(key);
                    # Also extract necessary information for field
                    if( key+'_LOCATIONS' in config ):
                        loc = config[key+'_LOCATIONS'].strip('()')
                        loc = loc.split(';')
                        try: # list of numbers given
                            float(loc[0].split(',')[0])
                            tmp = 1
                        except: # filename given
                            tmp = 0
                        if( tmp ): # list of numbers given
                            locList = []
                            for item in loc:
                                item2 = item.split(',')
                                locList.append([float(e) for e in item2])
                            nozzle.OutputLocations[key] = np.squeeze(np.array(locList))                        
                        else: # filename given
                            nozzle.OutputLocations[key] = np.loadtxt(loc[0])                               
                    else:
                        sys.stderr.write('\n ## ERROR : key %s_LOCATIONS ' \
                          'not found in config file to specify location of ' \
                          'requested output responses\n\n')
                else: # cycle through all possible more specific names for stress and temp
                    # check stress first         
                    assigned = 0;
                    for prefix in stress_prefix:
                        for k in dv_stress:
                            if key == prefix + '_' + k:
                                nozzle.GetOutput[k][stress_prefix.index(prefix)] = 1;
                                nozzle.Output_Tags.append(key);
                                assigned = 1;
                                continue;
                    for prefix in stress_prefix:
                        for k in dv_failure:
                            if key == prefix + '_' + k:
                                nozzle.GetOutput[k][stress_prefix.index(prefix)] = 1;
                                nozzle.Output_Tags.append(key);
                                assigned = 1;
                                continue;                                
                    for prefix in temp_prefix:
                        for k in dv_temp:
                            if key == prefix + '_' + k:
                                nozzle.GetOutput[k][temp_prefix.index(prefix)] = 1;
                                nozzle.Output_Tags.append(key);
                                assigned = 1;
                                continue;
                    for prefix in temp_prefix:
                        for k in dv_tempRatio:
                            if key == prefix + '_' + k:
                                nozzle.GetOutput[k][temp_prefix.index(prefix)] = 1;
                                nozzle.Output_Tags.append(key);
                                assigned = 1;
                                continue;
                    if assigned == 0:
                        string = '';
                        for k in dv_keys :
                            string = "%s %s " % (string,k);
                        sys.stderr.write('\n ## ERROR : Unknown output '          \
                          'function name : %s\n' % key);
                        sys.stderr.write('             Expected = %s\n\n' % string);
                        sys.exit(0);
                
            # for i in keys
        
        if len(nozzle.Output_Tags) == 0 :
            sys.stderr.write("\n  ## Error : No output function was given.\n\n");
            sys.exit(0);
        
        if output == 'verbose':
            sys.stdout.write('Setup Output Functions complete\n');


    def WriteOutputFunctions_Plain (self, output='verbose'):
      
        nozzle = self;

        filename = nozzle.Output_Name;
        
        if output == 'verbose':
            sys.stdout.write('\n');
            string = " Post-processing ";
            nch = (60-len(string))/2;
            sys.stdout.write('-' * nch);
            sys.stdout.write(string);
            sys.stdout.write('-' * nch);
            sys.stdout.write('\n\n');
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to '    \
              '"verbose" or "quiet" mode')
        
        try:
            fil = open(filename, 'w');
        except:
            sys.stderr.write("  ## ERROR : Could not open output file %s\n" % filename);
            sys.exit(0);
  
        if output == 'verbose':
            sys.stdout.write('  -- Info : Output functions file : %s\n' % filename);
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to '    \
              '"verbose" or "quiet" mode')
              
        # Stress and temperature outputs can have prefixes appended to them
        stress_prefix = [];
        temp_prefix = [];
        for i in range(len(nozzle.wall.layer)):
            stress_prefix.append(nozzle.wall.layer[i].name);
            temp_prefix.append(nozzle.wall.layer[i].name);
        if nozzle.stringers.n > 0:
            stress_prefix.append('STRINGERS');
        for i in range(nozzle.baffles.n):
            stress_prefix.append('BAFFLE' + str(i+1)); 

        # For printing to screen
        prt_item = list();
        prt_comp = list();
        prt_val = list();             

        for i in range(0, len(nozzle.Output_Tags)):
          
            tag = nozzle.Output_Tags[i];

            tag = nozzle.Output_Tags[i];
    
            if tag == 'MASS':
                fil.write('%0.16f\n' % nozzle.mass);
                prt_item.append('mass');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.mass); 
                    
            elif tag == 'MASS_WALL_ONLY':
                fil.write('%0.16f\n' % nozzle.mass_wall_only);
                prt_item.append('wall mass');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.mass_wall_only);                  
                    
            elif tag == 'VOLUME':
                fil.write('%0.16f\n' % nozzle.volume);
                prt_item.append('volume');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.volume); 
                    
            elif tag == 'THRUST':
                fil.write('%0.16f\n' % nozzle.thrust);
                prt_item.append('thrust');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.thrust); 

            elif tag == 'WALL_TEMPERATURE':
                for i in range(nozzle.wall_temperature.size):
                    fil.write('%0.16f\n' % nozzle.wall_temperature[i]);
                    prt_item.append('wall temp loc %i' % i);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.wall_temperature[i]);
                    
            elif tag == 'WALL_PRESSURE':
                for i in range(nozzle.wall_pressure.size):
                    fil.write('%0.16f\n' % nozzle.wall_pressure[i]);
                    prt_item.append('wall pressure loc %i' % i);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.wall_pressure[i]);
                    
            elif tag == 'PRESSURE':
                for i in range(nozzle.pressure.size):
                    fil.write('%0.16f\n' % nozzle.pressure[i]);
                    prt_item.append('pressure loc %i' % i);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.pressure[i]);  
                    
            elif tag == 'VELOCITY':
                nr, nc = nozzle.velocity.shape
                for i in range(nr):
                    fil.write('%0.16f\n' % nozzle.velocity[i,0]);
                    fil.write('%0.16f\n' % nozzle.velocity[i,1]);
                    fil.write('%0.16f\n' % nozzle.velocity[i,2]);
                    prt_item.append('velocity (x-dir) loc %i' % i);
                    prt_item.append('velocity (y-dir) loc %i' % i);
                    prt_item.append('velocity (z-dir) loc %i' % i);
                    prt_comp.append('');
                    prt_comp.append('');
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.velocity[i,0]);    
                    prt_val.append('%0.16f' % nozzle.velocity[i,1]);    
                    prt_val.append('%0.16f' % nozzle.velocity[i,2]);    
                                
            elif tag == 'KS_TOTAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.ks_total_stress[j]);
                    prt_item.append('ks total stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_total_stress[j]);
                    
            elif tag == 'PN_TOTAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.pn_total_stress[j]);
                    prt_item.append('pn total stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_total_stress[j]);
                    
            elif tag == 'MAX_TOTAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.max_total_stress[j]);
                    prt_item.append('max total stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_total_stress[j]);
                    
            elif tag == 'MAX_THERMAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.max_thermal_stress[j]);
                    prt_item.append('max thermal stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_thermal_stress[j]);
                    
            elif tag == 'MAX_MECHANICAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.max_mechanical_stress[j]);
                    prt_item.append('max mechanical stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_mechanical_stress[j]);
                    
            elif tag == 'KS_FAILURE_CRITERIA':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.ks_failure_criteria[j]);
                    prt_item.append('ks failure criteria');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_failure_criteria[j]);
                    
            elif tag == 'PN_FAILURE_CRITERIA':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.pn_failure_criteria[j]);
                    prt_item.append('pn failure criteria');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_failure_criteria[j]);
                    
            elif tag == 'MAX_FAILURE_CRITERIA':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f\n' % nozzle.max_failure_criteria[j]);
                    prt_item.append('max failure criteria');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_failure_criteria[j]);
                    
            elif tag == 'KS_TEMPERATURE':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f\n' % nozzle.ks_temperature[j]);
                    prt_item.append('ks temperature');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_temperature[j]);
                    
            elif tag == 'PN_TEMPERATURE':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f\n' % nozzle.pn_temperature[j]);
                    prt_item.append('pn temperature');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_temperature[j]);
                    
            elif tag == 'MAX_TEMPERATURE':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f\n' % nozzle.max_temperature[j]);
                    prt_item.append('max temperature');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_temperature[j]);  
                    
            elif tag == 'KS_TEMP_RATIO':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f\n' % nozzle.ks_temp_ratio[j]);
                    prt_item.append('ks temp ratio');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_temp_ratio[j]);                    

            elif tag == 'PN_TEMP_RATIO':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f\n' % nozzle.pn_temp_ratio[j]);
                    prt_item.append('pn temp ratio');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_temp_ratio[j]);
                    
            elif tag == 'MAX_TEMP_RATIO':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f\n' % nozzle.max_temp_ratio[j]);
                    prt_item.append('max temp ratio');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_temp_ratio[j]);
                    
            else:
                # Cycle through all possible specific stress outputs
                for prefix in stress_prefix:
                    j = stress_prefix.index(prefix);
                    if tag == prefix + '_KS_TOTAL_STRESS':
                        fil.write('%0.16f\n' % nozzle.ks_total_stress[j]);
                        prt_item.append('ks total stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_total_stress[j]);
                    elif tag == prefix + '_PN_TOTAL_STRESS':
                        fil.write('%0.16f\n' % nozzle.pn_total_stress[j]);
                        prt_item.append('pn total stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_total_stress[j]);
                    elif tag == prefix + '_MAX_TOTAL_STRESS':
                        fil.write('%0.16f\n' % nozzle.max_total_stress[j]);
                        prt_item.append('max total stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_total_stress[j]);
                    elif tag == prefix + '_MAX_THERMAL_STRESS':
                        fil.write('%0.16f\n' % nozzle.max_thermal_stress[j]);
                        prt_item.append('max thermal stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_thermal_stress[j]);
                # Cycle through all possible specific failure criteria outputs
                for prefix in stress_prefix:
                    j = stress_prefix.index(prefix);
                    if tag == prefix + '_KS_FAILURE_CRITERIA':
                        fil.write('%0.16f\n' % nozzle.ks_failure_criteria[j]);
                        prt_item.append('ks failure criteria');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_failure_criteria[j]);
                    elif tag == prefix + '_PN_FAILURE_CRITERIA':
                        fil.write('%0.16f\n' % nozzle.pn_failure_criteria[j]);
                        prt_item.append('pn failure criteria');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_failure_criteria[j]);
                    elif tag == prefix + '_MAX_FAILURE_CRITERIA':
                        fil.write('%0.16f\n' % nozzle.max_failure_criteria[j]);
                        prt_item.append('max failure criteria');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_failure_criteria[j]);
                # Cycle through all possible specific temperature outputs
                for prefix in temp_prefix:
                    j = temp_prefix.index(prefix);
                    if tag == prefix + '_KS_TEMPERATURE':
                        fil.write('%0.16f\n' % nozzle.ks_temperature[j]);
                        prt_item.append('ks temperature');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_temperature[j]);
                    elif tag == prefix + '_PN_TEMPERATURE':
                        fil.write('%0.16f\n' % nozzle.pn_temperature[j]);
                        prt_item.append('pn temperature');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_temperature[j]);
                    elif tag == prefix + '_MAX_TEMPERATURE':
                        fil.write('%0.16f\n' % nozzle.max_temperature[j]);
                        prt_item.append('max temperature');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_temperature[j]);
                # Cycle through all possible specific temperature ratio outputs
                for prefix in temp_prefix:
                    j = temp_prefix.index(prefix);
                    if tag == prefix + '_KS_TEMP_RATIO':
                        fil.write('%0.16f\n' % nozzle.ks_temp_ratio[j]);
                        prt_item.append('ks temp ratio');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_temp_ratio[j]);
                    elif tag == prefix + '_PN_TEMP_RATIO':
                        fil.write('%0.16f\n' % nozzle.pn_temp_ratio[j]);
                        prt_item.append('pn temp ratio');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_temp_ratio[j]);
                    elif tag == prefix + '_MAX_TEMP_RATIO':
                        fil.write('%0.16f\n' % nozzle.max_temp_ratio[j]);
                        prt_item.append('max temp ratio');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_temp_ratio[j]);
      
        if output == 'verbose':
            sys.stdout.write('\n');
        fil.close();
        
        # --- Print summary
        if output == 'verbose':
            sys.stdout.write('-' * 65);
            sys.stdout.write('\n%s | %s | %s\n' % ("Item".ljust(25), "Component".ljust(20),"Value".ljust(20)));
            sys.stdout.write('-' * 65);
            sys.stdout.write('\n');
            for i in range(0,len(prt_item)):
                sys.stdout.write('%s | %s | %s\n' % (prt_item[i].ljust(25), prt_comp[i].ljust(20),prt_val[i].ljust(20)));
            sys.stdout.write('-' * 65);    
            sys.stdout.write('\n\n');            
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')        

        
    def WriteOutputFunctions_Dakota (self,output='verbose'):
    
        nozzle = self;

        filename = nozzle.Output_Name;
        
        if output == 'verbose':
            sys.stdout.write('\n');
            string = " Post-processing ";
            nch = (60-len(string))/2;
            sys.stdout.write('-' * nch);
            sys.stdout.write(string);
            sys.stdout.write('-' * nch);
            sys.stdout.write('\n\n');
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to '    \
              '"verbose" or "quiet" mode')
        
        try:
            fil = open(filename, 'w');
        except:
            sys.stderr.write("  ## ERROR : Could not open output file %s\n" % filename);
            sys.exit(0);
  
        if output == 'verbose':
            sys.stdout.write('  -- Info : Output functions file : %s\n' % filename);
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to '    \
              '"verbose" or "quiet" mode')
              
        # Stress and temperature outputs can have prefixes appended to them
        stress_prefix = [];
        temp_prefix = [];
        for i in range(len(nozzle.wall.layer)):
            stress_prefix.append(nozzle.wall.layer[i].name);
            temp_prefix.append(nozzle.wall.layer[i].name);
        if nozzle.stringers.n > 0:
            stress_prefix.append('STRINGERS');
        for i in range(nozzle.baffles.n):
            stress_prefix.append('BAFFLE' + str(i+1)); 

        # For printing to screen
        prt_item = list();
        prt_comp = list();
        prt_val = list();             

        for i in range(0, len(nozzle.Output_Tags)):
          
            tag = nozzle.Output_Tags[i];

            tag = nozzle.Output_Tags[i];
    
            if tag == 'MASS':
                fil.write('%0.16f mass\n' % nozzle.mass);
                prt_item.append('mass');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.mass); 
                    
            elif tag == 'MASS_WALL_ONLY':
                fil.write('%0.16f mass_wall_only\n' % nozzle.mass_wall_only);
                prt_item.append('wall mass');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.mass_wall_only);                  
                    
            elif tag == 'VOLUME':
                fil.write('%0.16f volume\n' % nozzle.volume);
                prt_item.append('volume');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.volume); 
                    
            elif tag == 'THRUST':
                fil.write('%0.16f thrust\n' % nozzle.thrust);
                prt_item.append('thrust');
                prt_comp.append('');
                prt_val.append('%0.16f' % nozzle.thrust); 
                
            elif tag == 'WALL_TEMPERATURE':
                for i in range(nozzle.wall_temperature.size):
                    fil.write('%0.16f wall_temp_%i\n' % (nozzle.wall_temperature[i],i));
                    prt_item.append('wall temp loc %i' % i);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.wall_temperature[i]);
                    
            elif tag == 'WALL_PRESSURE':
                for i in range(nozzle.wall_pressure.size):
                    fil.write('%0.16f wall_pressure_%i\n' % (nozzle.wall_pressure[i],i));
                    prt_item.append('wall pressure loc %i' % i);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.wall_pressure[i]);
                    
            elif tag == 'PRESSURE':
                for i in range(nozzle.pressure.size):
                    fil.write('%0.16f pressure_%i\n' % (nozzle.pressure[i],i));
                    prt_item.append('pressure loc %i' % i);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.pressure[i]);  
                    
            elif tag == 'VELOCITY':
                nr, nc = nozzle.velocity.shape
                for i in range(nr):
                    fil.write('%0.16f velocity_x_%i\n' % (nozzle.velocity[i,0],i));
                    fil.write('%0.16f velocity_y_%i\n' % (nozzle.velocity[i,1],i));
                    fil.write('%0.16f velocity_z_%i\n' % (nozzle.velocity[i,2],i));
                    prt_item.append('velocity (x-dir) loc %i' % i);
                    prt_item.append('velocity (y-dir) loc %i' % i);
                    prt_item.append('velocity (z-dir) loc %i' % i);
                    prt_comp.append('');
                    prt_comp.append('');
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.velocity[i,0]);    
                    prt_val.append('%0.16f' % nozzle.velocity[i,1]);    
                    prt_val.append('%0.16f' % nozzle.velocity[i,2]);                 
                                
            elif tag == 'KS_TOTAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f ks_total_stress_%s\n' % (nozzle.ks_total_stress[j],nozzle.stressComponentList[j]));
                    prt_item.append('ks total stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_total_stress[j]);
                    
            elif tag == 'PN_TOTAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f pn_total_stress_%s\n' % (nozzle.pn_total_stress[j],nozzle.stressComponentList[j]));
                    prt_item.append('pn total stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_total_stress[j]);
                    
            elif tag == 'MAX_TOTAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f max_total_stress_%s\n' % (nozzle.max_total_stress[j],nozzle.stressComponentList[j]));
                    prt_item.append('max total stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_total_stress[j]);
                    
            elif tag == 'MAX_THERMAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f max_thermal_stress_%s\n' % (nozzle.max_thermal_stress[j],nozzle.stressComponentList[j]));
                    prt_item.append('max thermal stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_thermal_stress[j]);
                    
            elif tag == 'MAX_MECHANICAL_STRESS':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f max_mechanical_stress_%s\n' % (nozzle.max_mechanical_stress[j],nozzle.stressComponentList[j]));
                    prt_item.append('max mechanical stress');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_mechanical_stress[j]);
                    
            elif tag == 'KS_FAILURE_CRITERIA':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f ks_failure_criteria_%s\n' % (nozzle.ks_failure_criteria[j],nozzle.stressComponentList[j]));
                    prt_item.append('ks failure criteria');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_failure_criteria[j]);
                    
            elif tag == 'PN_FAILURE_CRITERIA':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f pn_failure_criteria_%s\n' % (nozzle.pn_failure_criteria[j],nozzle.stressComponentList[j]));
                    prt_item.append('pn failure criteria');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_failure_criteria[j]);
                    
            elif tag == 'MAX_FAILURE_CRITERIA':
                for j in range(len(nozzle.stressComponentList)):
                    fil.write('%0.16f max_failure_criteria_%s\n' % (nozzle.max_failure_criteria[j],nozzle.stressComponentList[j]));
                    prt_item.append('max failure criteria');
                    prt_comp.append('%s' % nozzle.stressComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_failure_criteria[j]);
                    
            elif tag == 'KS_TEMPERATURE':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f ks_temperature_%s\n' % (nozzle.ks_temperature[j],nozzle.tempComponentList[j]));
                    prt_item.append('ks temperature');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_temperature[j]);
                    
            elif tag == 'PN_TEMPERATURE':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f pn_temperature_%s\n' % (nozzle.pn_temperature[j],nozzle.tempComponentList[j]));
                    prt_item.append('pn temperature');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_temperature[j]);
                    
            elif tag == 'MAX_TEMPERATURE':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f max_temperature_%s\n' % (nozzle.max_temperature[j],nozzle.tempComponentList[j]));
                    prt_item.append('max temperature');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_temperature[j]);  
                    
            elif tag == 'KS_TEMP_RATIO':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f ks_temp_ratio_%s\n' % (nozzle.ks_temp_ratio[j],nozzle.tempComponentList[j]));
                    prt_item.append('ks temp ratio');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.ks_temp_ratio[j]);                    

            elif tag == 'PN_TEMP_RATIO':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f pn_temp_ratio_%s\n' % (nozzle.pn_temp_ratio[j],nozzle.tempComponentList[j]));
                    prt_item.append('pn temp ratio');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.pn_temp_ratio[j]);
                    
            elif tag == 'MAX_TEMP_RATIO':
                for j in range(len(nozzle.tempComponentList)):
                    fil.write('%0.16f max_temp_ratio_%s\n' % (nozzle.max_temp_ratio[j],nozzle.tempComponentList[j]));
                    prt_item.append('max temp ratio');
                    prt_comp.append('%s' % nozzle.tempComponentList[j]);
                    prt_val.append('%0.16f' % nozzle.max_temp_ratio[j]);
                    
            else:
                # Cycle through all possible specific stress outputs
                for prefix in stress_prefix:
                    j = stress_prefix.index(prefix);
                    if tag == prefix + '_KS_TOTAL_STRESS':
                        fil.write('%0.16f ks_total_stress_%s\n' % (nozzle.ks_total_stress[j],prefix));
                        prt_item.append('ks total stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_total_stress[j]);
                    elif tag == prefix + '_PN_TOTAL_STRESS':
                        fil.write('%0.16f pn_total_stress_%s\n' % (nozzle.pn_total_stress[j],prefix));
                        prt_item.append('pn total stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_total_stress[j]);
                    elif tag == prefix + '_MAX_TOTAL_STRESS':
                        fil.write('%0.16f max_total_stress_%s\n' % (nozzle.max_total_stress[j],prefix));
                        prt_item.append('max total stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_total_stress[j]);
                    elif tag == prefix + '_MAX_THERMAL_STRESS':
                        fil.write('%0.16f max_thermal_stress_%s\n' % (nozzle.max_thermal_stress[j],prefix));
                        prt_item.append('max thermal stress');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_thermal_stress[j]);
                # Cycle through all possible specific failure criteria outputs
                for prefix in stress_prefix:
                    j = stress_prefix.index(prefix);
                    if tag == prefix + '_KS_FAILURE_CRITERIA':
                        fil.write('%0.16f ks_failure_criteria_%s\n' % (nozzle.ks_failure_criteria[j],prefix));
                        prt_item.append('ks failure criteria');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_failure_criteria[j]);
                    elif tag == prefix + '_PN_FAILURE_CRITERIA':
                        fil.write('%0.16f pn_failure_criteria_%s\n' % (nozzle.pn_failure_criteria[j],prefix));
                        prt_item.append('pn failure criteria');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_failure_criteria[j]);
                    elif tag == prefix + '_MAX_FAILURE_CRITERIA':
                        fil.write('%0.16f max_failure_criteria_%s\n' % (nozzle.max_failure_criteria[j],prefix));
                        prt_item.append('max failure criteria');
                        prt_comp.append('%s' % stress_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_failure_criteria[j]);
                # Cycle through all possible specific temperature outputs
                for prefix in temp_prefix:
                    j = temp_prefix.index(prefix);
                    if tag == prefix + '_KS_TEMPERATURE':
                        fil.write('%0.16f ks_temperature_%s\n' % (nozzle.ks_temperature[j],prefix));
                        prt_item.append('ks temperature');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_temperature[j]);
                    elif tag == prefix + '_PN_TEMPERATURE':
                        fil.write('%0.16f pn_temperature_%s\n' % (nozzle.pn_temperature[j],prefix));
                        prt_item.append('pn temperature');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_temperature[j]);
                    elif tag == prefix + '_MAX_TEMPERATURE':
                        fil.write('%0.16f max_temperature_%s\n' % (nozzle.max_temperature[j],prefix));
                        prt_item.append('max temperature');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_temperature[j]);
                # Cycle through all possible specific temperature ratio outputs
                for prefix in temp_prefix:
                    j = temp_prefix.index(prefix);
                    if tag == prefix + '_KS_TEMP_RATIO':
                        fil.write('%0.16f ks_temp_ratio_%s\n' % (nozzle.ks_temp_ratio[j],prefix));
                        prt_item.append('ks temp ratio');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.ks_temp_ratio[j]);
                    elif tag == prefix + '_PN_TEMP_RATIO':
                        fil.write('%0.16f pn_temp_ratio_%s\n' % (nozzle.pn_temp_ratio[j],prefix));
                        prt_item.append('pn temp ratio');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.pn_temp_ratio[j]);
                    elif tag == prefix + '_MAX_TEMP_RATIO':
                        fil.write('%0.16f max_temp_ratio_%s\n' % (nozzle.max_temp_ratio[j],prefix));
                        prt_item.append('max temp ratio');
                        prt_comp.append('%s' % temp_prefix[j]);
                        prt_val.append('%0.16f' % nozzle.max_temp_ratio[j]);
      
        if output == 'verbose':
            sys.stdout.write('\n');
        fil.close();
        
        # --- Print summary
        if output == 'verbose':
            sys.stdout.write('-' * 65);
            sys.stdout.write('\n%s | %s | %s\n' % ("Item".ljust(25), "Component".ljust(20),"Value".ljust(20)));
            sys.stdout.write('-' * 65);
            sys.stdout.write('\n');
            for i in range(0,len(prt_item)):
                sys.stdout.write('%s | %s | %s\n' % (prt_item[i].ljust(25), prt_comp[i].ljust(20),prt_val[i].ljust(20)));
            sys.stdout.write('-' * 65);    
            sys.stdout.write('\n\n');            
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')        

        
    def Draw (self, output='verbose'):
        
        nozzle = self;
        
        sys.stdout.write("  -- Output a vectorized picture of the nozzle and material thicknesses.\n");
        
        FilNam = "nozzle.svg"
        
        wid = 1.3*nozzle.length;
        hei = 1.3*nozzle.length;
        
        SVGDim = [500, 500];
        margin = 50;
        
        xtab       = [];
        ytab       = [];

        
        nx = 100;
        _meshutils_module.py_BSplineGeo3 (nozzle.knots, nozzle.coefs, xtab, ytab, nx);
        
        # --- Get nozzle shape

        
        # --- Define scaling

        ymin = min(ytab);
        ymax = max(ytab);
        
        
        Box = [[0, nozzle.length],[ymin, ymax]];
                
        wid = Box[0][1] - Box[0][0];
        hei = Box[1][1] - Box[1][0];
        
        if ( wid > hei ) :
            sca = SVGDim[0] / wid;
        
        else :
            sca = SVGDim[1] / hei;
        
        wid =  wid*sca + 2*margin;
        hei =  hei*sca + 2*margin;
        #print Box
        
        #print "WID = %lf, HEI = %lf" % (wid, hei);
        
        for i in range(0,len(xtab)):
            xtab[i] = margin+sca*xtab[i];
            ytab[i] = -margin+ hei - (sca*ytab[i] - ymin * sca);
            
        
        
        #--- Write file
        
        try:
            fil = open(FilNam, 'w');
        except:
            sys.stderr.write("  ## ERROR : Could not open %s\n" % FilNam);
            return;
            
        fil.write("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
        fil.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
        fil.write("<svg version=\"1.1\" id=\"Layer_1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" x=\"%lfpx\" y=\"%lfpx\" width=\"%lfpx\" height=\"%lfpx\"  viewBox=\"0 0 %lf %lf\" enable-background=\"new 0 0 %lf %lf\" xml:space=\"preserve\">" % (0, 0, wid, hei, wid, hei, wid, hei));
        
        #fprintf(OutFil,"<g id=\"BdryEdges\" fill=\"blue\">");
        
        
        for i in range(1,len(xtab)):
            x = sca*xtab[i];
            y = sca*ytab[i];
            hl = sca*nozzle.wall.lower_thickness.radius(x);
            hu = sca*nozzle.wall.upper_thickness.radius(x);
            #print "x = %lf y = %lf lower thickness = %lf upper thickness = %lf " % (x, y, hl, hu);    
            fil.write("<g id=\"BdryEdges\" fill=\"blue\">");
            fil.write("<line fill=\"none\" stroke=\"%s\" stroke-miterlimit=\"1\" x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\"/>" \
            % ("black", xtab[i-1],  ytab[i-1],  xtab[i],  ytab[i]));
            fil.write("</g>\n");
            
            fil.write("<polygon id=\"tri\" points=\"%lf,%lf %lf,%lf %lf,%lf %lf,%lf\" style=\"fill:#717D7E;stroke:%s;stroke-width:0;fill-rule:nonzero;\" />" % \
              (xtab[i-1],  ytab[i-1],  xtab[i],  ytab[i],  xtab[i],  ytab[i]-hl, xtab[i-1],  ytab[i-1]-hl, "black"));
        
            fil.write("<polygon id=\"tri\" points=\"%lf,%lf %lf,%lf %lf,%lf %lf,%lf\" style=\"fill:#2E4053;stroke:%s;stroke-width:0;fill-rule:nonzero;\" />" % \
              (xtab[i-1],  ytab[i-1]-hl,  xtab[i],  ytab[i]-hl,  xtab[i],  ytab[i]-hl-hu, xtab[i-1],  ytab[i-1]-hl-hu, "black"));
        
        fil.write("</svg>");
        
        sys.stdout.write("  -- Info : nozzle.svg OPENED.\n\n");

def NozzleSetup( config_name, flevel, output='verbose', partitions=1 ):
    import tempfile
        
    if not os.path.isfile(config_name) :
        sys.stderr.write("  ## ERROR : could not find configuration file %s\n\ns" % config_name);
        sys.exit(0);
    
    nozzle = multif.nozzle.nozzle.Nozzle()
    
    config = SU2.io.Config(config_name)
    
    if output == 'verbose':
        print config;

    nozzle.mesh_name    =  'nozzle.su2';
    nozzle.restart_name =  'nozzle.dat';
    nozzle.CONV_FILENAME = 'history'

    nozzle.partitions = partitions;
    
    if 'TEMP_RUN_DIR' in config:
        if config['TEMP_RUN_DIR'] == 'YES':
            nozzle.runDir       =  tempfile.mkdtemp();    
        else:
            nozzle.runDir = '';
    else:
        nozzle.runDir = '';
    
    # --- Path to SU2 exe
    
    if 'SU2_RUN' in config:
        nozzle.SU2_RUN = config['SU2_RUN'];
    else:
        nozzle.SU2_RUN = os.environ['SU2_RUN'];
        
    if 'SU2_OUTPUT_FORMAT' in config:
        nozzle.su2_output_format = config['SU2_OUTPUT_FORMAT'];
    else:
        nozzle.su2_output_format = 'TECPLOT';
    
    # --- Parse fidelity levels
    
    nozzle.SetupFidelityLevels(config, flevel, output);
    
    # --- Set flight regime + fluid
    
    nozzle.SetupMission(config,output);
    
    # --- Setup inner wall & parameterization (B-spline)
    
    nozzle.SetupBSplineCoefs(config,output);
    
    # --- Setup inner wall temperature definition (if necessary)
    
    nozzle.SetupWallTemp(config,output);
    
    # --- Setup materials
    
    nozzle.SetupMaterials(config,output);
    
    # --- Setup wall layer thickness(es) and material(s)
    
    nozzle.SetupWallLayers(config,output);
    
    # --- Setup baffles
    
    nozzle.SetupBaffles(config,output);
    
    # --- Setup stringers
    
    nozzle.SetupStringers(config,output);
    
    # --- Setup DV definition
        
    nozzle.SetupDV(config,output);
    
    # --- If input DV are provided, parse them and update nozzle
    
    if nozzle.NbrDVTot > 0 :
        
        # Parse DV from input DV file (plain or dakota format)
        nozzle.ParseDV(config,output);
        
        # Update DV using values provided in input DV file
        nozzle.UpdateDV(config,output);
    
    # --- Computer inner wall's B-spline and thermal and load layer thicknesses
    #     B-spline coefs, and thickness node arrays may have been updated by
    #     the design variables input file; update exterior geometry & baffles
    
    nozzle.SetupWall(config,output);
    
    # --- Get output functions to be returned
    
    nozzle.SetupOutputFunctions(config,output);
        
    #sys.exit(1);
    return nozzle;    
