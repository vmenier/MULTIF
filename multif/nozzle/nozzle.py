# -*- coding: utf-8 -*-


"""

R. Fenrich & V. Menier, July 2016

""" 

import os, sys, tempfile
import textwrap
import copy

import numpy as np

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

from parserDV import *

class CFD:
    def __init__(self):
        pass

class Nozzle:
    def __init__(self):
        pass

     
    def AssignListDV(self,config,dvList,key,NbrDV,sizeDV):
        
        if isinstance(key,tuple): # multiple keys provided
            dv = [];
            for i in range(len(key)):
                
                if key[i] not in config:
                    sys.stderr.write('\n ## ERROR: Expected %s in config ' \
                      'file.\n\n' % key[i]);
                    sys.exit(0);
                    
                tmp = config[key[i]].strip('()');
                tmp = tmp.split(',');
                dv = dv + tmp;
                
            dv_size = len(dv);
        
        else:
            if key not in config:
                sys.stderr.write('\n ## ERROR: Expected %s in config ' \
                  'file.\n\n' % key);
                sys.exit(0);            
            dv = config[key].strip('()');
            dv = dv.split(",");
            dv_size = len(dv);         
        
        # Check that user-stated number of design variables in DV_LIST matches
        # that from component definition (usually a COMPONENT_DV keyword)
        if  dv_size != sizeDV:
			sys.stderr.write('  ## ERROR : Inconsistent number of '      \
			  'design variables found in %s: %d instead of %d\n\n' %     \
			  (key,dv_size,sizeDV));
			print dv
			sys.exit(0);
        
        tag = np.zeros(NbrDV);
        for iw in range(0,dv_size):

            val = int(dv[iw]);

            if ( val > NbrDV or val < 0 ):
                sys.stderr.write('  ## ERROR : Inconsistent design '     \
                  'variables provided in %s (idx=%d). Check DV_LIST '    \
                  'and %s_DV keyword specifications.\n\n' % (key,val,key));
                sys.stderr.write('idDV %d NbrDV %d\n' % (val, NbrDV));
                sys.exit(0);

            dvList.append(val-1);

            if val-1 >= 0 :
                tag[val-1] = 1;
        
        # Check for missing design variables
        for iw in range(0,NbrDV):
            if tag[iw] != 1 :
                sys.stderr.write('  ## ERROR : Design variable %d '      \
                  'in %s is not defined.\n\n' % (iw+1,key));
                sys.exit(0);  

        
    def SetupDV (self, config, output='verbose'):
        nozzle = self;
        
        nozzle.DV_Tags = []; # tags: e.g. wall, tstag etc.
        nozzle.DV_Head = []; # correspondence b/w DV_Tags and dvList
        NbrDVTot = 0; # keep track of number of DV            
             
        if 'DV_LIST' in config:

            # Build a list of all possible expected keys. The list contains
            # other lists which represent possible more specific keys that can
            # be specified.
            dv_keys = list();
            dv_n = list(); # record allowable number of dv for each key
            
            # Wall (inner wall shape)
            if nozzle.param == '3D':
                dv_keys.append(['WALL']);
                dv_n.append([nozzle.wall.centerline.coefs_size + 
                             nozzle.wall.majoraxis.coefs_size + 
                             nozzle.wall.minoraxis.coefs_size]);
                dv_keys.append(['WALL_SHOVEL_HEIGHT']);
                dv_n.append([1]);
                dv_keys.append(['WALL_SHOVEL_START_ANGLE']);
                dv_n.append([1]);
            else: # assume 2D parameterization
                dv_keys.append(['WALL']);
                dv_n.append([len(nozzle.wall.coefs)]);
            
            # Wall temperature
            if nozzle.param == '3D':
                if hasattr(nozzle.wall,'temperature'):
                    sys.stdout.write('\nWARNING: WALL_TEMP keys are being ' \
                      'ignored in DV_LIST since 3D parameterization is used.\n\n');
            else: # assume 2D parameterization
                if hasattr(nozzle.wall,'temperature'):
                    if nozzle.wall.temperature.param == 'PIECEWISE_LINEAR':
                        ltemp = ['WALL_TEMP', 'WALL_TEMP_LOCATIONS',
                                 'WALL_TEMP_VALUES'];
                        dv_keys.append(ltemp);
                        dv_n.append([2*nozzle.wall.temperature.nBreaks,
                                     nozzle.wall.temperature.nBreaks,
                                     nozzle.wall.temperature.nBreaks]);
                    elif nozzle.wall.temperature.param == 'KLE':
                        dv_keys.append('WALL_TEMP_COEFS');
                        dv_n.append(len(nozzle.wall.temperature.coefs));
            
            # Wall layers
            if nozzle.param == '3D':
                for i in range(len(nozzle.wall.layer)): # assume piecewise linear
                    if nozzle.wall.layer[i].param == 'PIECEWISE_BILINEAR':
                        ltemp = [nozzle.wall.layer[i].name, 
                                 nozzle.wall.layer[i].name + '_THICKNESS_LOCATIONS',
                                 nozzle.wall.layer[i].name + '_THICKNESS_ANGLES', 
                                 nozzle.wall.layer[i].name + '_THICKNESS_VALUES'];
                        dv_keys.append(ltemp);
                        dv_n.append([nozzle.wall.layer[i].nAxialBreaks+
                                     nozzle.wall.layer[i].nAngularBreaks+
                                     nozzle.wall.layer[i].nBreaks,
                                     nozzle.wall.layer[i].nAxialBreaks,
                                     nozzle.wall.layer[i].nAngularBreaks,
                                     nozzle.wall.layer[i].nBreaks]);
                    elif nozzle.wall.layer[i].param == 'PIECEWISE_LINEAR':
                        sys.stderr.write('\n## ERROR: 3D parameterization ' \
                          'does not implement a PIECEWISE_LINEAR param., ' \
                          'only a PIECEWISE_BILINEAR param.\n\n');
                        sys.exit(0);
                    elif nozzle.wall.layer[i].param == 'CONSTANT':
                        ltemp = [nozzle.wall.layer[i].name + '_THICKNESS'];
                        dv_keys.append(ltemp);
                        dv_n.append([1]);
            else: # assume 2D parameterization
                for i in range(len(nozzle.wall.layer)): # assume piecewise linear
                    if nozzle.wall.layer[i].param == 'PIECEWISE_LINEAR':
                        ltemp = [nozzle.wall.layer[i].name, 
                                 nozzle.wall.layer[i].name + '_THICKNESS_LOCATIONS', 
                                 nozzle.wall.layer[i].name + '_THICKNESS_VALUES'];
                        dv_keys.append(ltemp);
                        dv_n.append([2*nozzle.wall.layer[i].nBreaks,
                                     nozzle.wall.layer[i].nBreaks,
                                     nozzle.wall.layer[i].nBreaks]);
                    elif nozzle.wall.layer[i].param == 'CONSTANT':
                        ltemp = [nozzle.wall.layer[i].name + '_THICKNESS'];
                        dv_keys.append(ltemp);
                        dv_n.append([1]);                            
            
            # Baffles
            dv_keys.append(['BAFFLES','BAFFLES_LOCATION','BAFFLES_THICKNESS',
                            'BAFFLES_HEIGHT','BAFFLES_HALF_WIDTH']);
            dv_n.append([3*nozzle.baffles.n+1,nozzle.baffles.n,
                        nozzle.baffles.n,nozzle.baffles.n,1]);
                        
            # Stringers
            if nozzle.param == '3D':
                dv_keys.append(['STRINGERS','STRINGERS_BREAK_LOCATIONS',
                                'STRINGERS_ANGLES','STRINGERS_HEIGHT_VALUES',
                                'STRINGERS_THICKNESS_VALUES']);
                dv_n.append([nozzle.stringers.nAxialBreaks+
                             nozzle.stringers.n+
                             nozzle.stringers.nThicknessValues,
                             nozzle.stringers.nAxialBreaks,
                             nozzle.stringers.n,
                             nozzle.stringers.nThicknessValues,
                             nozzle.stringers.nThicknessValues]);                
            else: # assume 2D parameterization
                dv_keys.append(['STRINGERS','STRINGERS_BREAK_LOCATIONS',
                                'STRINGERS_HEIGHT_VALUES',
                                'STRINGERS_THICKNESS_VALUES']);
                dv_n.append([3*nozzle.stringers.nBreaks,
                             nozzle.stringers.nBreaks,
                             nozzle.stringers.nBreaks,
                             nozzle.stringers.nBreaks]);
                         
            # Materials
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
                    
            # Inlet, atmosphere, and heat transfer properties
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
            
            # First check that design variables have not been overspecified,
            # i.e. BAFFLES and BAFFLES_LOCATION or BAFFLES_THICKNESS have not
            # both been specified in DV_LIST.
            for k in dv_keys:
                if len(k) > 1:
                    if k[0] in hdl:
                        for i in range(1,len(k)):
                            if k[i] in hdl:
                                sys.stderr.write('\n ## ERROR: %s and %s '  \
                                'cannot both be specified in DV_LIST\n\n' % \
                                (k[0],k[i]));
                                sys.exit(0);
            
            # Loop through all user-specified design variables and assign
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

                # Search for design variables whose definitions differ b/w
                # 2D and 3D parameterizations
                if nozzle.param == '3D':
                    
                    if key == 'WALL':
                        nozzle.wall.dv = [];
                        sizeTemp = nozzle.wall.centerline.coefs_size + \
                                   nozzle.wall.majoraxis.coefs_size + \
                                   nozzle.wall.minoraxis.coefs_size;
                        nozzle.AssignListDV(config,nozzle.wall.dv, 
                          ('WALL_COEFS1_DV','WALL_COEFS2_DV','WALL_COEFS3_DV'),
                          NbrDV,sizeTemp);
                        continue;
                        
                    # Check all layers with non-specific names, e.g. LAYER1, etc.
                    check = 0;
                    for j in range(len(nozzle.wall.layer)):
                        if key == nozzle.wall.layer[j].name:
                            nozzle.wall.layer[j].dv = [];
                            # total number of possible DV for checking, assumes
                            # piecewise bilinear definition
                            sizeTemp = nozzle.wall.layer[j].nAxialBreaks + \
                                       nozzle.wall.layer[j].nAngularBreaks + \
                                       nozzle.wall.layer[j].nBreaks;
                            nozzle.AssignListDV(config,nozzle.wall.layer[j].dv,
                                              nozzle.wall.layer[j].formalName+'_DV',
                                              NbrDV,sizeTemp);
                            check = 1;
                    if check == 1:
                        continue;

                    if key == 'STRINGERS':
                        nozzle.stringers.dv = [];
                        # total number of possible DV for checking
                        sizeTemp = nozzle.stringers.nAxialBreaks + \
                                   nozzle.stringers.n + \
                                   nozzle.stringers.nThicknessValues;
                        nozzle.AssignListDV(config,nozzle.stringers.dv,'STRINGERS_DV',
                                          NbrDV,sizeTemp);
                        continue;
                
                # Setup design variables according to 2D parameterization
                else:
                    if key == 'WALL':
                        nozzle.wall.dv = [];
                        sizeTemp = nozzle.wall.coefs_size; # number of DV for checking
                        nozzle.AssignListDV(config,nozzle.wall.dv,'WALL_COEFS_DV',
                                          NbrDV,sizeTemp);
                        continue;

                    # Check all layers with non-specific names, e.g. LAYER1, etc.
                    check = 0;
                    for j in range(len(nozzle.wall.layer)):
                        if key == nozzle.wall.layer[j].name:
                            nozzle.wall.layer[j].dv = [];
                            # number of DV for checking; assumes 2xN array
                            sizeTemp = 2*nozzle.wall.layer[j].nBreaks;
                            nozzle.AssignListDV(config,nozzle.wall.layer[j].dv,
                                              nozzle.wall.layer[j].formalName+'_DV',
                                              NbrDV,sizeTemp);
                            check = 1;
                    if check == 1:
                        continue;

                    if key == 'STRINGERS':
                        nozzle.stringers.dv = [];
                        sizeTemp = nozzle.stringers.nBreaks*3; # number of DV for checking
                        self.AssignListDV(config,nozzle.stringers.dv,'STRINGERS_DV',
                                          NbrDV,sizeTemp);
                        continue;
                        
                # Design variable definitions shared between 2D and 3D param.
                        
                if key == 'WALL_TEMP':
                    nozzle.wall.temperature.dv = [];
                    sizeTemp = 2*nozzle.wall.temperature.nBreaks;
                    self.AssignListDV(config,nozzle.wall.temperature.dv,
                                      'WALL_TEMP_DV',NbrDV,sizeTemp);
                    continue;
                                            
                if key == 'BAFFLES':                    
                    nozzle.baffles.dv = [];
                    sizeTemp = nozzle.baffles.n*3+1; # number of DV for checking
                    self.AssignListDV(config,nozzle.baffles.dv,'BAFFLES_DV',
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
        
        if output == 'verbose':
            print config;

        # Setup parameterization
        if 'PARAMETERIZATION' in config:
            if config['PARAMETERIZATION'] == '2D' or config['PARAMETERIZATION'] == '3D':
                nozzle.param = config['PARAMETERIZATION'];
            else:
                sys.stderr.write('\n ## ERROR: PARAMETERIZATION=2D or 3D must' \
                  ' be specified. %s found instead.\n\n' % config['PARAMETERIZATION']);
                sys.exit(0);
        else:
            sys.stderr.write('\n ## ERROR: Keyword PARAMETERIZATION not found ' \
              'in the config file. Please specify PARAMETERIZATION= 2D or ' \
              '3D.\n\n');
            sys.exit(0);

        # Determine fidelity levels        
        fidelity_tags = config['FIDELITY_LEVELS_TAGS'].strip('()');
        fidelity_tags = fidelity_tags.split(",");

        NbrFidLev = len(fidelity_tags);

        if NbrFidLev < 1:
            sys.stdout.write("\n ## ERROR: No fidelity level was defined.\n");
            sys.exit(0);

        if flevel >= NbrFidLev:
            sys.stderr.write("\n ## ERROR: Fidelity level %d not defined.\n" % flevel);
            sys.exit(0);			

        if output == 'verbose':
            sys.stdout.write('\n%d fidelity level(s) defined using %s ' \
              'parameterization. Summary :\n' % (NbrFidLev, nozzle.param));

            sys.stdout.write('-' * 90);
            sys.stdout.write('\n%s | %s | %s\n' % ("Level #".ljust(10),         \
                "Tag".ljust(10),"Description.".ljust(70)));
            sys.stdout.write('-' * 90);
            sys.stdout.write('\n');

        # Cycle through each fidelity level and compile information. Assign 
        # information to nozzle if fidelity level is the one requested to be run.
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

            # Assign nozzle analysis method for fidelity level i
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
                    
                    # Set nozzle analysis dimension
                    nozzle.dim = '1D';
                    
                    # Set tolerances
                    nozzle.tolerance = tolerance.Tolerance();
                    nozzle.tolerance.setRelTol(tol);
                    nozzle.tolerance.setAbsTol(tol);
                    nozzle.tolerance.exitTempPercentError = 10*tol;
                    nozzle.solverApparentThroatLocation = tol;

                    # Set analysis type
                    try:
                        analysisType = cfgLvl[2];
                    except:
                        sys.stderr.write('\n ## WARNING : Analysis type '     \
                          'could not be determined. AEROTHERMOSTRUCTURAL, '       \
                          'AEROTHERMAL, or AERO keyword must be provided '  \
                          'in the model definition of fidelity level %d.\n\n' \
                          % flevel);
                        sys.exit(0);
                        
                    # Set thermostructural parameters if necessary
                    if len(cfgLvl) == 5:
                        if cfgLvl[3] == 'LINEAR':
                            nozzle.linearStructuralAnalysisFlag = 1;
                        elif cfgLvl[3] == 'NONLINEAR':
                            nozzle.linearStructuralAnalysisFlag = 0;
                            
                        nozzle.thermostructuralFidelityLevel = float(cfgLvl[4]);
                        if( nozzle.thermostructuralFidelityLevel < 0 or 
                        nozzle.thermostructuralFidelityLevel > 1):
                            sys.stderr.write('\n ## ERROR: thermostructural ' \
                              'fidelity level must range from 0 (low) to 1 '\
                              '(high). %f provided instead.\n\n' % 
                              nozzle.thermostructuralFidelityLevel);
                            sys.exit(0);
                    else:
                        nozzle.linearStructuralAnalysisFlag = 1;
                        nozzle.thermostructuralFidelityLevel = 0.5;
                        
            elif method == 'RANS' or method == 'EULER':

                nozzle.cfd.max_cfl = 30.0;

                dim = cfgLvl[1];
                if dim != '2D' and dim != '3D':
                    sys.stderr.write("\n ## ERROR : Wrong dimension for "     \
                      "fidelity level %d (tagged %s) : only 2D or 3D "        \
                      "simulations\n\n" % (i,tag));
                    sys.exit(0);
                
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
                    nozzle.cfd.su2_convergence_order = int(config['SU2_CONVERGENCE_ORDER']);
                else:
                    nozzle.cfd.su2_convergence_order = 6;
                description += ", relative convergence order %i" % nozzle.cfd.su2_convergence_order;
                
                if 'SU2_OUTPUT_FORMAT' in config:
                    nozzle.cfd.output_format = config['SU2_OUTPUT_FORMAT'];
                else:
                    nozzle.cfd.output_format = 'TECPLOT';

                # --- Setup max iterations for SU2
                if 'SU2_MAX_ITERATIONS' in config:
                    nozzle.cfd.su2_max_iterations = int(config['SU2_MAX_ITERATIONS']);
                elif nozzle.dim == '2D':
                    if nozzle.method == 'EULER':
                        nozzle.cfd.su2_max_iterations = 600;
                    else:
                        nozzle.cfd.su2_max_iterations = 1000;   
                else:
                    if nozzle.method == 'EULER':
                        nozzle.cfd.su2_max_iterations = 1200;
                    else:
                        nozzle.cfd.su2_max_iterations = 5000;
                description += ", max iterations %i" % nozzle.cfd.su2_max_iterations;
                  
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
                    
                    # Set nozzle analysis dimension
                    nozzle.dim = dim;

                    # Set mesh level
                    nozzle.cfd.mesh_size = meshsize;
                    nozzle.cfd.bl_ds        = 0.000007;
                    nozzle.cfd.bl_ratio     = 1.3; 
                    nozzle.cfd.bl_thickness = 0.02;
                    nozzle.cfd.bl_yplus     = 1.0;
					
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
                    	    nozzle.cfd.bl_yplus = 3.0;
                    	elif meshsize == 'MEDIUM':
                    	    scaleMesh = 1.5;
                    	    nozzle.cfd.bl_yplus = 2.0;
                    	elif meshsize == 'FINE':
                    	    scaleMesh = 0.9;
                    	    nozzle.cfd.bl_yplus = 1.0;	

                    nozzle.cfd.meshhl = scaleMesh*np.asarray([0.1, 0.07, 0.06, 0.006, 0.0108]);

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
                            nozzle.linearStructuralAnalysisFlag = 1;
                        elif cfgLvl[4] == 'NONLINEAR':
                            nozzle.linearStructuralAnalysisFlag = 0;      
                            
                        nozzle.thermostructuralFidelityLevel = float(cfgLvl[5]);
                        if( nozzle.thermostructuralFidelityLevel < 0 or 
                        nozzle.thermostructuralFidelityLevel > 1):
                            sys.stderr.write('\n ## ERROR: thermostructural ' \
                              'fidelity level must range from 0 (low) to 1 '\
                              '(high). %f provided instead.\n\n' % 
                              nozzle.thermostructuralFidelityLevel);
                            sys.exit(0);
                    else:
                        nozzle.linearStructuralAnalysisFlag = 1;
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
                  
        # Setup thermal and structural analysis
        if analysisType == 'AEROTHERMOSTRUCTURAL':
            nozzle.aeroFlag = 1;
            nozzle.thermalFlag = 1;
            nozzle.structuralFlag = 1;
        elif analysisType == 'AEROTHERMAL':
            nozzle.aeroFlag = 1;
            nozzle.thermalFlag = 1;
            nozzle.structuralFlag = 0;
        elif analysisType == 'AEROSTRUCTURAL':
            nozzle.aeroFlag = 1;
            nozzle.thermalFlag = 0;
            nozzle.structuralFlag = 1;
        elif analysisType == 'THERMOSTRUCTURAL':
            nozzle.aeroFlag = 0;
            nozzle.thermalFlag = 1;
            nozzle.structuralFlag = 1;
        elif analysisType == 'AERO':
            nozzle.aeroFlag = 1;
            nozzle.thermalFlag = 0;
            nozzle.structuralFlag = 0;
        else:
            sys.stderr.write('\n ## ERROR: AEROTHERMOSTRUCTURAL, '        \
              'THERMOSTRUCTURAL, AEROTHERMAL, AEROSTRUCTURAL, or AERO ' \
              'must be provided as a keyword for analyis '    \
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
    
    
    def SetupInnerWall(self, config, output='verbose'):
        # Assume B-spline is a 3rd-degree B-spline. Thus, given the coefs, the
        # knots can be calculated assuming evenly-spaced knots with a B-spline
        # that terminates at both ends and not earlier.
        
        nozzle = self;
        
        # Set up nozzle wall
        if nozzle.param == '3D':
            
            nozzle.wall = component.NonaxisymmetricWall();
            wall_keys = ('WALL','WALL_COEFS1','WALL_COEFS2','WALL_COEFS3')
            
            if all (key in config for key in wall_keys):
                
                # Check that a B-spline definition will be provided
                hdl = config['WALL'].strip('()');
                hdl = hdl.split(',');
                if not all (item == 'BSPLINE' for item in hdl):
                    sys.stderr.write('\n ## ERROR: 3D parameterization ' \
                      'currently only accepts a BSPLINE parameterization for ' \
                      'nozzle wall definition.\n\n');
                    sys.exit(0);
           
                # Obtain B-spline coefs for nozzle centerline
                hdl = config['WALL_COEFS1'].strip('()');
                hdl = hdl.split(',');
                nozzle.wall.centerline = component.Spline('CENTERLINE');
                nozzle.wall.centerline.coefs = [float(k) for k in hdl];
                n = len(hdl)/2-3; # calc knots w/ 3rd-deg B-spline assumptions
                knots = [0,0,0,0] + range(1,n) + [n,n,n,n];
                nozzle.wall.centerline.knots = [float(k) for k in knots];
                nozzle.wall.centerline.coefs_size = len(hdl);
                
                # Obtain B-spline coefs for nozzle elliptical cross section 
                # major axis
                hdl = config['WALL_COEFS2'].strip('()');
                hdl = hdl.split(',');
                nozzle.wall.majoraxis = component.Spline('MAJORAXIS');
                nozzle.wall.majoraxis.coefs = [float(k) for k in hdl];
                n = len(hdl)/2-3; # calc knots w/ 3rd-deg B-spline assumptions
                knots = [0,0,0,0] + range(1,n) + [n,n,n,n];
                nozzle.wall.majoraxis.knots = [float(k) for k in knots];
                nozzle.wall.majoraxis.coefs_size = len(hdl);
                
                # Obtain B-spline coefs for nozzle elliptical cross section 
                # minor axis
                hdl = config['WALL_COEFS3'].strip('()');
                hdl = hdl.split(',');
                nozzle.wall.minoraxis = component.Spline('MINORAXIS');
                nozzle.wall.minoraxis.coefs = [float(k) for k in hdl];
                n = len(hdl)/2-3; # calc knots w/ 3rd-deg B-spline assumptions
                knots = [0,0,0,0] + range(1,n) + [n,n,n,n];
                nozzle.wall.minoraxis.knots = [float(k) for k in knots];
                nozzle.wall.minoraxis.coefs_size = len(hdl);

            else:
                string = '';
                for key in wall_keys :
                    string = "%s %s " % (string,key);
                sys.stderr.write('\n ## ERROR: No inner wall definition is ' \
                                 'provided.\n\n');
                sys.stderr.write('             Expected = %s\n\n' % string);
                sys.exit(0);
                
            # Setup parameters governing shovel exit
            nozzle.wall.shovel_height = float(config['WALL_SHOVEL_HEIGHT'].strip('()'));
            nozzle.wall.shovel_start_angle = float(config['WALL_SHOVEL_START_ANGLE'].strip('()'));
                
        else: # assume nozzle.param = 2D
        
            nozzle.wall = component.AxisymmetricWall();
            wall_keys = ('WALL','WALL_COEFS');
            
            if all (key in config for key in wall_keys):

                hdl = config['WALL_COEFS'].strip('()');
                hdl = hdl.split(",");
                nozzle.wall.coefs = [float(k) for k in hdl];
                n = len(hdl)/2-3; # calc knots w/ 3rd-deg B-spline assumptions
                knots = [0,0,0,0] + range(1,n) + [n,n,n,n];
                nozzle.wall.knots = [float(k) for k in knots];
                nozzle.wall.coefs_size = len(hdl);

            else:
                string = '';
                for key in wall_keys :
                    string = "%s %s " % (string,key);
                sys.stderr.write('\n ## ERROR: No inner wall definition is ' \
                                 'provided.\n\n');
                sys.stderr.write('             Expected = %s\n\n' % string);
                sys.exit(0);
        
        if output == 'verbose':
            sys.stdout.write('Setup B-Spline Coefs complete\n');        


    def SetupWallTemp(self, config, output='verbose'): 
        
        nozzle = self;

        nozzle.wallTempFlag = 0;
        
        if 'WALL_TEMP' not in config:
            return 0;
        else:
            if output == 'verbose':
                sys.stdout.write('Wall temp will be fixed\n');
            nozzle.wallTempFlag = 1;
                
        definition = config['WALL_TEMP'].strip('()');
        if definition == 'KLE':
            
            sys.stderr.write('\n ## ERROR: KLE definition for wall temp not ' \
              'implemented yet.\n\n');
            sys.exit(0);
            
            #nozzle.wall.temperature = component.Distribution('WALL_TEMP');
            #nozzle.wall.temperature.param = 'KLE';
            #nozzle.wall.temperature.coefs = config['WALL_TEMP_COEFS'].strip('()');
            
        elif definition == 'PIECEWISE_LINEAR':
                        
            loc = config['WALL_TEMP_LOCATIONS'].strip('()');
            val = config['WALL_TEMP_VALUES'].strip('()');
            if ',' not in loc:
                xCoord = np.loadtxt(loc);
                yCoord = np.loadtxt(val);
            else:
                xCoord = np.array([float(a) for a in loc.split(',')]);
                yCoord = np.array([float(a) for a in val.split(',')]);
            thicknessNodesNonDim = np.transpose(np.array([xCoord,yCoord]));
                
            nozzle.wall.temperature = component.Distribution('WALL_TEMP');
            nozzle.wall.temperature.param = 'PIECEWISE_LINEAR';
            nozzle.wall.temperature.thicknessNodesNonDim = thicknessNodesNonDim; 
            nozzle.wall.temperature.nBreaks = xCoord.size;                            
                
        else:
            sys.stderr.write('\n ## ERROR: %s wall temp definition not '      \
              'accepted. Only KLE and PIECEWISE_LINEAR are implemented.\n\n');
            sys.exit(0);
            
        if output == 'verbose':
            sys.stdout.write('Setup Wall Temp complete\n'); 
            

    def ParseThicknessBilinear(self, config, key, loc='_THICKNESS_LOCATIONS', ang='_THICKNESS_ANGLES', val='_THICKNESS_VALUES'):
        
        # Keys in config file governing bilinear thickness distribution
        loc_name = key + loc;
        ang_name = key + ang;
        val_name = key + val;

        # Obtain axial coordinate (normalized to nozzle length) of bilinear break locations
        x = config[loc_name].strip('()');
        x = x.split(',');
        x = [float(i) for i in x];
        
        # Obtain angular coordinate (in degrees) of bilinear break locations
        theta = config[ang_name].strip('()');
        theta = theta.split(',');
        theta = [float(i) for i in theta];
        
        # Obtain thickness values at each break location
        t = config[val_name].strip('()');
        t = t.split(',');
        t = [float(i) for i in t];
        
        # Check that there are enough thickness values for each break node
        if len(t) != len(x)*len(theta):
            sys.stderr.write('\n## ERROR: Number of thickness values does ' \
              'not match number of thickness break locations for %s.\n\n' % \
              key);
            sys.exit(0);
            
        # Check that axial coordinate is normalized
        if sum([1 for k in x if k < 0.] + [1 for k in x if k > 1.]) > 0:
            sys.stderr.write('\n## ERROR: %s should contain axial coord. ' \
              'values normalized to the length of the nozzle.\n\n' % loc_name);
            sys.exit(0);
        
        # Check that angles are not too large or small
        if sum([1 for k in x if k < 0.] + [1 for k in x if k > 360.]) > 0.:
            sys.stderr.write('\n## ERROR: %s should contain angular location ' \
              'of bilinear node breaks in degrees, specified between 0 and ' \
              '360 with 0 corresponding to y axis.\n\n' % ang_name);
            sys.exit(0);
            
        # Check that first and last normalized x-coordinate is 0 and 1
        if x[0] != 0 and x[-1] != 1.0:
            sys.stderr.write('\n## ERROR: %s normalized axial coord. should' \
              ' start at 0 and end at 1.\n\n' % key);
            sys.exit(0);
            
        # Assemble thickness definition into a Numpy array
        wallThickness = np.zeros((len(t),3));
        for i in range(len(t)):
            wallThickness[i,0] = x[np.mod(i,len(x))]; # axial coordinate at node
            wallThickness[i,1] = theta[int(np.floor(i/len(x)))]; # angular coord. at node
            wallThickness[i,2] = t[i]; # thickness at node
            
        # Number of breaks in axial direction and angular direction
        nAxialBreaks = len(x);
        nAngularBreaks = len(theta);
        nBreaks = len(t);
        
        return wallThickness, nAxialBreaks, nAngularBreaks, nBreaks;
            
        
    def ParseThicknessLinear(self, config, key, loc='_THICKNESS_LOCATIONS', val = '_THICKNESS_VALUES'):
        
        # Keys in config file giving piecewise linear thickness distribution
        loc_name = key + loc;
        val_name = key + val;

        # Obtain axial coordinate (normalized to nozzle length) of bilinear break locations
        x = config[loc_name].strip('()');
        x = x.split(',');
        x = [float(i) for i in x];

        # Obtain thickness values at each break location
        t = config[val_name].strip('()');
        t = t.split(',');
        t = [float(i) for i in t];

        # Check that there are enough thickness values for each break node
        if len(t) != len(x):
            sys.stderr.write('\n## ERROR: Number of thickness values does ' \
              'not match number of thickness break locations for %s.\n\n' % \
              key);
            sys.exit(0);
            
        # Check that axial coordinate is normalized
        if sum([1 for k in x if k < 0.] + [1 for k in x if k > 1.]) > 0:
            sys.stderr.write('\n## ERROR: %s should contain axial coord. ' \
              'values normalized to the length of the nozzle.\n\n' % loc_name);
            sys.exit(0);

        # Check that first and last normalized x-coordinate is 0 and 1
        if x[0] != 0 and x[-1] != 1.0:
            sys.stderr.write('\n## ERROR: %s normalized axial coord. should' \
              ' start at 0 and end at 1.\n\n' % key);
            sys.exit(0);
      
        # Assemble thickness definition into Numpy array
        wallThickness = np.zeros((len(t),2));
        for i in range(len(t)):
            wallThickness[i,0] = x[np.mod(i,len(x))]; # axial coordinate at node
            wallThickness[i,1] = t[i]; # thickness at node

        # Number of breaks in piecewise linear function 
        nBreaks = len(x);
        
#        wallThickness = [[0 for i in range(2)] for j in range(nBreaks)];
#                    
#        for i in range(0,size_loc):
#            wallThickness[i][0] = float(hdl[i]);
#        
#        hdl = config[val_name].strip('()');
#        hdl = hdl.split(",");
#        size_val = len(hdl);
#        
#        for i in range(0,size_val):
#            wall_thickness[i][1] = float(hdl[i]);
            
        return wallThickness, nBreaks;

    
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
                layer.thicknessNodesNonDim = self.ParseThicknessLinear(config,name);
            except:
                sys.stderr.write('\n ## ERROR : Thickness definition '       \
                     'could not be parsed for %s.\n\n' % name);
                sys.exit(0);
        else:
            if name == 'LOAD_LAYER':
                layer.thicknessNodesNonDim = [[0.0,0.016], [1.0, 0.016]];
            else: # THERMAL_LAYER
                layer.thicknessNodesNonDim = [[0.0,0.01], [1.0, 0.01]];

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
        
        if layer.param == 'PIECEWISE_BILINEAR':
            # Setup thickness keys
            t_keys = ['_THICKNESS_LOCATIONS','_THICKNESS_ANGLES','_THICKNESS_VALUES'];
            for i in range(len(t_keys)):
                t_keys[i] = layer.formalName + t_keys[i];
            t_keys = tuple(t_keys);

            # Assign thicknesses          
            if all (k in config for k in t_keys):
                t, n, m, p = self.ParseThicknessBilinear(config,layer.formalName);
                layer.thicknessNodesNonDim = t; # Numpy array (n x 3) with break locations & thickness values
                layer.nAxialBreaks = n; # number of breaks in axial direction
                layer.nAngularBreaks = m; # number of breaks in a plane perp. to axis at axial break
                layer.nBreaks = p; # total number of breaks
            else:
                sys.stderr.write('\n ## ERROR: All the following keys should ' \
                  'be found in the config file: %s.\n\n' % t_keys);
                sys.exit(0);
        
        elif layer.param == 'PIECEWISE_LINEAR':
            # Setup thickness keys
            t_keys = ['_THICKNESS_LOCATIONS','_THICKNESS_VALUES'];
            for i in range(len(t_keys)):
                t_keys[i] = layer.formalName + t_keys[i];
            t_keys = tuple(t_keys);

            # Assign thicknesses
            if all (k in config for k in t_keys):
                t, m = self.ParseThicknessLinear(config,layer.formalName);
                layer.thicknessNodesNonDim = t;
                layer.nBreaks = m;                
            else:
                sys.stderr.write('\n ## ERROR: All the following keys should ' \
                  'be found in the config file: %s.\n\n' % t_keys);
                sys.exit(0); 
        
        elif layer.param == 'CONSTANT': # implement as piecewise linear
            t_value = float(config[layer.formalName + '_THICKNESS']);
            layer.thicknessNodesNonDim = np.array([[0.,t_value],[1.,t_value]]);
            layer.nBreaks = 2;
            layer.thicknessValue = t_value;
        else:
            sys.stderr.write('\n ## ERROR: Parameterization %s '     \
            'is not valid for %s. Only PIECEWISE_BILINEAR (for 3D param ' \
            'or PIECEWISE_LINEAR or '    \
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
                
                if nozzle.param == '3D':
                    
                    nozzle.wall.layer.append(component.NonaxisymmetricWall(name));
                    nozzle.wall.layer[i-1].param = param;
                    nozzle.wall.layer[i-1].formalName = key;
                    
                    # Assign thickness distribution
                    nozzle.SetupLayerThickness(config,nozzle.wall.layer[i-1]);
                    # Assign material
                    nozzle.wall.layer[i-1].material=nozzle.materials[material];
                    
                else: # assume parameterization is 2D
                
                    nozzle.wall.layer.append(component.AxisymmetricWall(name));
                    nozzle.wall.layer[i-1].param = param;
                    nozzle.wall.layer[i-1].formalName = key;
                    
                    # Assign thickness distribution
                    nozzle.SetupLayerThickness(config,nozzle.wall.layer[i-1]);
                    # Assign material
                    nozzle.wall.layer[i-1].material=nozzle.materials[material];
                
            else:
                break;
            i += 1

        if output == 'verbose':
            sys.stdout.write('%d layers processed\n' % (i-1));
            sys.stdout.write('Setup Wall Layers complete\n');

        
    def SetupBaffles(self, config, output='verbose'):
        
        nozzle = self;
        
        info = config['BAFFLES'].strip('()');
        n, material = [x.strip() for x in info.split(',')];
        nozzle.baffles = component.Baffles(n);
        nozzle.baffles.material = nozzle.materials[material];
        
        # Non-dimensional baffles location (normalized w.r.t. nozzle length)
        info = config['BAFFLES_LOCATION'].strip('()');
        ltemp = [x.strip() for x in info.split(',')];
        nozzle.baffles.locationNonDim = [float(x) for x in ltemp];
        
        info = config['BAFFLES_HEIGHT'].strip('()');
        nozzle.baffles.heightSetToExterior = 0;
        if( info == 'EXTERIOR' ):
            nozzle.baffles.height = [];
            nozzle.baffles.heightSetToExterior = 1;
        else:
            ltemp = [x.strip() for x in info.split(',')];        
            nozzle.baffles.height = [float(x) for x in ltemp];
        
        info = config['BAFFLES_THICKNESS'].strip('()');
        ltemp = [x.strip() for x in info.split(',')];          
        nozzle.baffles.thickness = [float(x) for x in ltemp];

        # For 3D param, the dimensional distance the baffle extends from the 
        # centerline in the transverse direction
        if( 'BAFFLES_HALF_WIDTH' in config ):
            info = config['BAFFLES_HALF_WIDTH'].strip('()');
            nozzle.baffles.halfWidth = float(info);
        
        if output == 'verbose':
            sys.stdout.write('Setup Baffles complete\n');        

        
    def SetupStringers(self, config, output='verbose'):
        
        nozzle = self;
        
        info = config['STRINGERS'].strip('()');
        n, material = [x.strip() for x in info.split(',')];
        nozzle.stringers = component.Stringers(n);
        nozzle.stringers.material = nozzle.materials[material];
        nozzle.stringers.absoluteRadialCoord = 0; # 0 if stringer height is measured relative to outside layer, 1 if measured relative to nozzle axis
        nozzle.stringers.locationDefinition = 'INDEPENDENT';
        nozzle.stringers.heightDefinition = 'INDEPENDENT';
        
        if ( 'STRINGERS_BREAK_LOCATIONS' not in config or 
             config['STRINGERS_BREAK_LOCATIONS'] == 'BAFFLES_LOCATION' ):
             
                 nozzle.stringers.locationDefinition = 'BAFFLES_LOCATION';
             
                 if output == 'verbose':
                     sys.stdout.write('Stringer break locations will be set' \
                       ' from baffle locations\n');
                 key = '';
                 k_loc = 'BAFFLES_LOCATION';
                 k_ang = 'STRINGERS_ANGLES';
                 if ( 'STRINGERS_HEIGHT_VALUES' not in config ):
                     if output == 'verbose':
                         sys.stdout.write('Stringer heights will be set'  \
                               ' to exterior geometry\n');
                     h_val = 'EXTERIOR';     
                     nozzle.stringers.absoluteRadialCoord = 1;
                 elif ( config['STRINGERS_HEIGHT_VALUES'] == 'BAFFLES_HEIGHT' ):
                     if output == 'verbose':
                         sys.stdout.write('Stringer heights will be set'  \
                               ' from baffle heights\n');
                     h_val = 'BAFFLES_HEIGHT';             
                 elif ( config['STRINGERS_HEIGHT_VALUES'] == 'EXTERIOR' ):
                     if output == 'verbose':
                         sys.stdout.write('Stringer heights will be set'  \
                               ' to exterior geometry\n');
                     h_val = 'EXTERIOR';
                     nozzle.stringers.absoluteRadialCoord = 1;
                 else:
                     h_val = 'STRINGERS_HEIGHT_VALUES';
                 t_val = 'STRINGERS_THICKNESS_VALUES';
                 
                 nozzle.stringers.heightDefinition = h_val;
                 
        else:
            if ( 'STRINGERS_HEIGHT_VALUES' not in config or
                config['STRINGERS_HEIGHT_VALUES'] == 'BAFFLES_HEIGHT'):
                    sys.stderr.write('\n ## ERROR: Stringer height values '   \
                      'cannot be set to baffle height values if stringers '   \
                      'break locations are not set to baffles location.\n\n')
                    sys.exit(0);

            key = 'STRINGERS';
            k_loc = '_BREAK_LOCATIONS';
            k_ang = '_ANGLES';
            t_val = '_THICKNESS_VALUES';
            h_val = '_HEIGHT_VALUES';
            
            if ( config['STRINGERS_HEIGHT_VALUES'] == 'EXTERIOR' ):
                h_val = 'EXTERIOR';
                nozzle.stringers.heightDefinition = h_val;
                nozzle.stringers.absoluteRadialCoord = 1;
        
        # Assign thickness distribution
        if nozzle.param == '3D':
            t, m, n, p = self.ParseThicknessBilinear(
              config,key,loc=k_loc,ang=k_ang,val=t_val);
            if( nozzle.stringers.n != n ):
                sys.stderr.write('\n ## ERROR: %i stringers defined for 3D ' \
                  'parameterization, but only %i stringer angles provided.' \
                  % (nozzle.stringers.n,n));
                sys.exit(0);
            nozzle.stringers.thicknessNodesNonDim = t;
            nozzle.stringers.nAxialBreaks = m;
            nozzle.stringers.nAngularBreaks = n;
            nozzle.stringers.nThicknessValues = p;
        else: # assume 2D parameterization            
            t, m = self.ParseThicknessLinear(config,key,loc=k_loc,val=t_val);
            nozzle.stringers.thicknessNodesNonDim = t;
            nozzle.stringers.nBreaks = m;
        
        # Assign height distribution
        if config['STRINGERS_HEIGHT_VALUES'] == 'EXTERIOR':
            if nozzle.param == '3D':
                nozzle.stringers.heightNodesNonDim = np.array([[0,0,0],[0,0,0]],dtype='float64'); # temporary assignment
            else:
                nozzle.stringers.heightNodesNonDim = np.array([[0,0],[0,0]],dtype='float64');
        else:
            if nozzle.param == '3D':
                t, m, n, p = self.ParseThicknessBilinear(
                  config,key,loc=k_loc,ang=k_ang,val=h_val);
                nozzle.stringers.heightNodesNonDim = t;
            else: # assume 2D parameterization
                nozzle.stringers.heightNodesNonDim, m = self.ParseThicknessLinear(
                  config,key,loc=k_loc,val=h_val);    
        
        if output == 'verbose':
            sys.stdout.write('Setup Stringers complete\n'); 
        

    def SetupWall (self, output='verbose'):
        
        nozzle = self;
        
        # ====================================================================
        # Setup shape of inner wall (B-spline)
        # ====================================================================
        
        if nozzle.param == '3D':

            nozzle.xinlet = nozzle.wall.centerline.coefs[0];
            nozzle.xoutlet = nozzle.wall.centerline.coefs[nozzle.wall.centerline.coefs_size/2-1];
            nozzle.zinlet = nozzle.wall.centerline.coefs[nozzle.wall.centerline.coefs_size/2];
            nozzle.length = nozzle.wall.centerline.coefs[nozzle.wall.centerline.coefs_size/2-1] - \
                            nozzle.wall.centerline.coefs[0];

            # The following parameters are currently hard-coded but could be obtained from 
            # nozzle.wall.shovel_start_angle and nozzle.wall.shovel_height
            nozzle.wall.shovel_start_angle = 1.572865;
            nozzle.wall.shovel_height = 0.122638;
            nozzle.wall.shovel_end_angle = np.arccos((nozzle.wall.shovel_height - \
                nozzle.wall.centerline.coefs[-1])/nozzle.wall.minoraxis.coefs[-1])

            if nozzle.dim == '3D':
                
                nozzle.wall.centerline.geometry = geometry.Bspline(nozzle.wall.centerline.coefs);
                nozzle.wall.majoraxis.geometry  = geometry.Bspline(nozzle.wall.majoraxis.coefs);
                nozzle.wall.minoraxis.geometry  = geometry.Bspline(nozzle.wall.minoraxis.coefs);            
                
                # Build equivalent nozzle shape based on equivalent area
                majoraxisTmp = geometry.Bspline(nozzle.wall.majoraxis.coefs);
                minoraxisTmp = geometry.Bspline(nozzle.wall.minoraxis.coefs);
                centerTmp = geometry.Bspline(nozzle.wall.centerline.coefs);
                nx = 2000; # Check accuracy and effect of this interpolation
                x = np.linspace(nozzle.xinlet,nozzle.xoutlet,num=nx);               
                
                fr1 = majoraxisTmp.radius
                fr2 = minoraxisTmp.radius
                fz = centerTmp.radius
                
                xi = nozzle.wall.centerline.coefs[0];
                inletRadius = np.sqrt(majoraxisTmp.radius(xi)*minoraxisTmp.radius(xi));
                #params = np.zeros(5);
                #params[0] = inletRadius*np.sin(nozzle.wall.shovel_start_angle*np.pi/180.);
                #params[1] = nozzle.wall.centerline.coefs[-1] + nozzle.wall.shovel_height;
                ##params[2]                
                #params[3] = nozzle.xinlet;
                #params[4] = nozzle.xoutlet;
                #equivRadius = multif.HIGHF.MF_GetRadius (x, fr1, fr2, fz, params)
                
                equivRadius = multif.HIGHF.MF_GetRadius (x, nozzle);
                
                shape2d = np.transpose(np.array([x,equivRadius]))
                nozzle.wall.geometry = geometry.PiecewiseLinear(shape2d);	
                
            # Map 3D parameterization to 2D definition using equivalent area
            else:

                # No shovel is present in axisymmetric geometry
                nozzle.wall.shovel_start_angle = np.pi;
                nozzle.wall.shovel_end_angle = np.pi;
                
                # Build equivalent nozzle shape based on equivalent area
                majoraxisTmp = geometry.Bspline(nozzle.wall.majoraxis.coefs);
                minoraxisTmp = geometry.Bspline(nozzle.wall.minoraxis.coefs);
                centerTmp = geometry.Bspline(nozzle.wall.centerline.coefs);
                nx = 2000; # Check accuracy and effect of this interpolation
                x = np.linspace(nozzle.xinlet,nozzle.xoutlet,num=nx);                
                
                #fr1 = majoraxisTmp.radius
                #fr2 = minoraxisTmp.radius
                #fz = centerTmp.radius             

                xi = nozzle.wall.centerline.coefs[0];
                inletRadius = np.sqrt(majoraxisTmp.radius(xi)*minoraxisTmp.radius(xi));
                #params = np.zeros(5);
                #params[0] = inletRadius*np.sin(nozzle.wall.shovel_start_angle*np.pi/180.);
                #params[1] = nozzle.wall.centerline.coefs[-1] + nozzle.wall.shovel_height;
                ##params[2]                
                #params[3] = nozzle.xinlet;
                #params[4] = nozzle.xoutlet;
                #equivRadius = multif.HIGHF.MF_GetRadius (x, fr1, fr2, fz, params);
                
                equivRadius = multif.MEDIUMF.Get3Dto2DEquivArea(nozzle,x);
                
                shape2d = np.transpose(np.array([x,equivRadius]));
                nozzle.wall.geometry = geometry.PiecewiseLinear(shape2d);
                
                # The following info is required by SU2 for 2D geometries
                if nozzle.dim == '2D':
                	
                    nx = 4000; # use 4000 points to interpolate inner wall shape
                    x = np.linspace(nozzle.xinlet,nozzle.xoutlet,num=nx);
                    y = multif.MEDIUMF.Get3Dto2DEquivArea(nozzle, x);
                    
                    nozzle.cfd.x_wall = x;
                    nozzle.cfd.y_wall = y;

                    dx_exit = max(1.3*nozzle.cfd.meshhl[3], 0.001);
                    for i in range(0,nx) :
                    	if (  x[nx-i-1] < x[-1]-dx_exit  ):
                    		nozzle.cfd.x_thrust = x[nx-i-1];
                    		nozzle.cfd.y_thrust = y[nx-i-1];
                    		break;
                      
        else: # assume 2D parameterization

            # No shovel is present in axisymmetric geometry
            nozzle.wall.shovel_start_angle = np.pi;
            nozzle.wall.shovel_end_angle = np.pi;

            nozzle.xinlet = nozzle.wall.coefs[0];
            nozzle.xoutlet = nozzle.wall.coefs[nozzle.wall.coefs_size/2-1];
            nozzle.zinlet = 0.;
            nozzle.length = nozzle.wall.coefs[nozzle.wall.coefs_size/2-1] - \
                            nozzle.wall.coefs[0];
        	
            # Upscale 2D param to 3D, this may be useful for comparisons
            if nozzle.dim == '3D':

                xi = nozzle.wall.coefs[0]; # inlet x-coord
                xe = nozzle.wall.coefs[nozzle.wall.coefs_size/2-1]; # exit x-coord   

                centerCoefs = np.array([xi, xi, xe, xe, 0., 0., 0., 0.]);

                # Generate axisymmetric shape, which is used later
                nozzle.wall.geometry = geometry.Bspline(nozzle.wall.coefs);

                # Create centerline
                nozzle.wall.centerline = component.Spline('CENTERLINE');
                nozzle.wall.centerline.coefs = list(centerCoefs);
                n = len(centerCoefs)/2-3;
                knots = [0,0,0,0] + range(1,n) + [n,n,n,n];
                nozzle.wall.centerline.knots = [float(k) for k in knots];
                nozzle.wall.centerline.coefs_size = len(centerCoefs);
                nozzle.wall.centerline.geometry = geometry.Bspline(centerCoefs);

                # Create major axis
                nozzle.wall.majoraxis = component.Spline('MAJORAXIS');
                nozzle.wall.majoraxis.coefs = nozzle.wall.coefs;
                nozzle.wall.majoraxis.knots = nozzle.wall.knots;
                nozzle.wall.majoraxis.coefs_size = nozzle.wall.coefs_size;              
                nozzle.wall.majoraxis.geometry = geometry.Bspline(nozzle.wall.coefs);
                
                # Create minor axis
                nozzle.wall.minoraxis = nozzle.wall.majoraxis;
                                
            else:
                                            
                nozzle.wall.geometry = geometry.Bspline(nozzle.wall.coefs);
                         
                if nozzle.method == 'RANS' or nozzle.method == 'EULER':  
                    x = [];
                    y = [];
                    nx = 4000; # use 4000 points to interpolate inner wall shape            
                    _meshutils_module.py_BSplineGeo3 (nozzle.wall.knots, \
                                                      nozzle.wall.coefs, x, y, nx);
                
                    nozzle.cfd.x_wall = x;
                    nozzle.cfd.y_wall = y;
                    
                    dx_exit = max(1.3*nozzle.cfd.meshhl[3], 0.001);   
                    for i in range(0,nx) :
                    	if (  x[nx-i-1] < x[-1]-dx_exit  ):
                    		nozzle.cfd.x_thrust = x[nx-i-1];
                    		nozzle.cfd.y_thrust = y[nx-i-1];
                    		break;

        # ====================================================================        
        # Setup inner wall temperature if necessary
        # ====================================================================  
      
        if hasattr(nozzle.wall,'temperature'):
            if nozzle.param == '3D':
                sys.stderr.write('\nWARNING: Wall temp control not setup ' \
                ' for 3D parameterization.\n\n');
                sys.exit(0);
            else: # assume 2D parameterization           
                if nozzle.wall.temperature.param == 'PIECEWISE_LINEAR':
                    nozzle.wall.temperature.thicknessNodes = \
                        copy.copy(nozzle.wall.temperature.thicknessNodesNonDim);
                    nozzle.wall.temperature.thicknessNodes[:,0] = \
                        nozzle.wall.temperature.thicknessNodes[:,0]* \
                        nozzle.length + nozzle.xinlet;
                    nozzle.wall.temperature.geometry = \
                      geometry.PiecewiseLinear(nozzle.wall.temperature.thicknessNodes);                
                elif nozzle.wall.temperature.param == 'KLE':
                    sys.stderr.write('\n ## ERROR: KLE not implemented yet for '  \
                      'inner nozzle wall temperature definition.\n\n');
                    sys.exit(0);

        # ====================================================================        
        # Setup thickness of each layer
        # ====================================================================
        
        for i in range(len(nozzle.wall.layer)):
            
            # Dimensionalize thickness node x-coordinate
            nozzle.wall.layer[i].thicknessNodes = \
                copy.copy(nozzle.wall.layer[i].thicknessNodesNonDim);
            nozzle.wall.layer[i].thicknessNodes[:,0] = \
                nozzle.wall.layer[i].thicknessNodes[:,0]*nozzle.length + nozzle.xinlet;
            
            if nozzle.param == '3D':
                
                if nozzle.dim == '3D':    
                    
                    if( nozzle.wall.layer[i].param == 'PIECEWISE_BILINEAR' ):
                        
                        nozzle.wall.layer[i].thickness = geometry.PiecewiseBilinear( \
                            nozzle.wall.layer[i].nAxialBreaks, \
                            nozzle.wall.layer[i].nAngularBreaks, \
                            nozzle.wall.layer[i].thicknessNodes);
                            
                    elif( nozzle.wall.layer[i].param == 'CONSTANT' ):
                        
                        nozzle.wall.layer[i].thickness = geometry.PiecewiseLinear( \
                            nozzle.wall.layer[i].thicknessNodes);                       
                
                # Map 3D parameterization to 2D definition using average thickness
                else:
                    
                    if( nozzle.wall.layer[i].param == 'PIECEWISE_BILINEAR' ):
                        
                        nx = nozzle.wall.layer[i].nAxialBreaks;
                    
                        # Find mean thickness at each x-coordinate station
                        xStations = nozzle.wall.layer[i].thicknessNodes[0:nx,0];
                        yStations = list(nozzle.wall.layer[i].thicknessNodes[:,1])[::nx];
                        yStations.append(360.);
                        yStations = np.array(yStations);
                        dy = yStations[1:]-yStations[0:-1];
                        tMean = np.zeros((nx,));
                        for j in range(nx):
                            tTemp = list(nozzle.wall.layer[i].thicknessNodes[:,2])[j::nx];
                            tTemp.append(tTemp[0]);
                            tTemp = np.array(tTemp);
                            tMean[j] = np.sum((tTemp[1:]+tTemp[0:-1])/2.*dy)/360.;
                        thicknessNodes = np.transpose(np.array([xStations,tMean]));
    
                        nozzle.wall.layer[i].thickness = geometry.PiecewiseLinear( \
                            thicknessNodes);
                        
                    elif( nozzle.wall.layer[i].param == 'CONSTANT' ):
                        
                        nozzle.wall.layer[i].thickness = geometry.PiecewiseLinear( \
                            nozzle.wall.layer[i].thicknessNodes);                        
                
            else: # assume 2D parameterization
            
                # upscale 2D param to 3D, this may be useful for comparisons
                if nozzle.dim == '3D':
                    
                    nx = nozzle.wall.layer[i].thicknessNodes.size/2;
                    nAngles = [0., 90., 180., 270.];
                    thicknessNodes = np.zeros((nx*4,3));
                    for j in range(len(nAngles)):
                        for k in range(nx):
                            thicknessNodes[k+j*nx,0] = nozzle.wall.layer[i].thicknessNodes[k,0];
                            thicknessNodes[k+j*nx,1] = nAngles[j];
                            thicknessNodes[k+j*nx,2] = nozzle.wall.layer[i].thicknessNodes[k,1];
                            
                    nozzle.wall.layer[i].thickness = geometry.PiecewiseBilinear( \
                        nx, 4, thicknessNodes);
                        
                else:
                    nozzle.wall.layer[i].thickness = geometry.PiecewiseLinear( \
                        nozzle.wall.layer[i].thicknessNodes);

        # ====================================================================
        # Setup nozzle exterior
        # ====================================================================

        nozzle.exterior = component.NonaxisymmetricWall('exterior');
        nozzle.exterior.geometry = {};   

        if nozzle.dim == '3D' and nozzle.param == '3D':
            zoutlettop = nozzle.wall.centerline.geometry.radius( 
                nozzle.wall.centerline.geometry.xend) + \
                nozzle.wall.minoraxis.geometry.radius(
                nozzle.wall.minoraxis.geometry.xend);
            zoutletbottom = nozzle.wall.shovel_height;
        elif nozzle.dim == '3D' and nozzle.param == '2D': # no shovel
            xoutlet = nozzle.wall.centerline.geometry.xend;
            zoutlettop = nozzle.wall.centerline.geometry.radius(xoutlet) + \
                nozzle.wall.minoraxis.geometry.radius(xoutlet) + \
                sum([nozzle.wall.layer[i].thickness.height(xoutlet,90.)
                    for i in range(len(nozzle.wall.layer))]);
            zoutletbottom = nozzle.wall.centerline.geometry.radius(xoutlet) - \
                nozzle.wall.minoraxis.geometry.radius(xoutlet) - \
                sum([nozzle.wall.layer[i].thickness.height(xoutlet,270.)
                    for i in range(len(nozzle.wall.layer))]);
        elif nozzle.dim != '3D': # regardless of param, there is no shovel
            xoutlet = nozzle.wall.geometry.xend;
            zoutlettop = nozzle.wall.geometry.radius(xoutlet) + \
                sum(nozzle.wall.layer[i].thickness.radius(xoutlet) 
                    for i in range(len(nozzle.wall.layer)));
            zoutletbottom = -zoutlettop;
        
        # Approximate top surface of internal aircraft cavity which nozzle
        # lies within
        nozzle.exterior.geometry['top'] = geometry.EllipticalExterior('top',
            nozzle.xoutlet, zoutlettop=zoutlettop, zoutletbottom=zoutletbottom);

        # Approximate bottom surface of internal aircraft cavity which 
        # nozzle lies within
        nozzle.exterior.geometry['bottom'] = geometry.EllipticalExterior(
            'bottom', nozzle.xoutlet, zoutlettop=zoutlettop, 
            zoutletbottom=zoutletbottom);

        # # Test nozzle exterior geometry
        # import matplotlib.pyplot as plt
        # from mpl_toolkits.mplot3d import Axes3D 
        # fig = plt.figure()
        # ax = plt.axes(projection='3d')
        # xsection = np.linspace(nozzle.xinlet,nozzle.xoutlet,10);
        # for i in range(len(xsection)):
        #     # Plot top of aircraft cavity
        #     xloc = nozzle.length*float(i)/float(len(xsection)-1) + nozzle.xinlet;
        #     antmp = np.linspace(0.,np.pi,100);
        #     xtmp = xloc*np.ones((len(antmp),1));
        #     ytmp = np.zeros((len(antmp),1));
        #     ztmp = np.zeros((len(antmp),1));
        #     for j in range(len(antmp)):
        #         tmp1, tmp2 = nozzle.exterior.geometry['top'].coord(xloc,antmp[j]);
        #         ytmp[j] = tmp1;
        #         ztmp[j] = tmp2;
        #     ax.scatter(xtmp,ytmp,ztmp);

        #     # Plot bottom of aircraft cavity
        #     xloc = nozzle.length*float(i)/float(len(xsection)-1) + nozzle.xinlet;
        #     antmp = np.linspace(np.pi,2*np.pi,100);
        #     xtmp = xloc*np.ones((len(antmp),1));
        #     ytmp = np.zeros((len(antmp),1));
        #     ztmp = np.zeros((len(antmp),1));
        #     for j in range(len(antmp)):
        #         tmp1, tmp2 = nozzle.exterior.geometry['bottom'].coord(xloc,antmp[j]);
        #         ytmp[j] = tmp1;
        #         ztmp[j] = tmp2;
        #     ax.scatter(xtmp,ytmp,ztmp);

        # plt.show();

        # ====================================================================
        # Setup baffles
        # ====================================================================

        # Dimensionalize baffle x-location
        nozzle.baffles.location = [q*nozzle.length + nozzle.xinlet for q in \
                                   nozzle.baffles.locationNonDim];

        if( nozzle.baffles.heightSetToExterior ):

            # There is nothing to be done. Mass and volume will be approximated
            # correctly and baffle shape only affects structural FEA. This
            # is taken into account correctly by AERO-S.
            pass;

        else:

            raise NotImplementedError("Baffle height is set fixed according " + \
                "to shape of internally-defined external geometry.");

        # ====================================================================
        # Setup stringers
        # ====================================================================

        # Dimensionalize stringers thickness x-coordinate
        nozzle.stringers.thicknessNodes = \
            copy.copy(nozzle.stringers.thicknessNodesNonDim);
        nozzle.stringers.thicknessNodes[:,0] = \
            nozzle.stringers.thicknessNodes[:,0]*nozzle.length + nozzle.xinlet;
        nozzle.stringers.heightNodes = \
            copy.copy(nozzle.stringers.heightNodesNonDim);
        nozzle.stringers.heightNodes[:,0] = \
            nozzle.stringers.heightNodes[:,0]*nozzle.length + nozzle.xinlet;
            
        if nozzle.param == '3D':
            
            if nozzle.dim == '3D':
            
                # Setup thickness distribution
                nozzle.stringers.thickness = [];
                ns = nozzle.stringers.n;
                na = nozzle.stringers.nAxialBreaks;
                for i in range(ns):
                    # extract 0th and 2nd column of data
                    nozzle.stringers.thickness.append(geometry.PiecewiseLinear( \
                           nozzle.stringers.thicknessNodes[i*na:i*na+na,0::2]));
                    nozzle.stringers.thickness[-1].angle = nozzle.stringers.thicknessNodes[i*na,1];
                
                # Setup height distribution
                if( nozzle.stringers.heightDefinition == 'EXTERIOR' or \
                    nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT' ):
                    # There is no height to assign as it is defined by the 
                    # exterior. Volume, mass, and structural FEA are taken 
                    # care of accordingly.
                    if nozzle.stringers.n != 2:
                        raise NotImplementedError("Stringers definition with " + \
                        "EXTERIOR or BAFFLES_HEIGHT height definition can only have 2 " + \
                        "stringers located on the top and bottom of nozzle.")                                
                    pass;
                else:
                    nozzle.stringers.height = [];
                    for i in range(ns):
                        # extract 0th and 2nd column of data
                        nozzle.stringers.height.append(geometry.PiecewiseLinear( \
                               nozzle.stringers.heightNodes[i*na:i*na+na,0::2]));    
                        nozzle.stringers.height[-1].angle = nozzle.stringers.heightNodes[i*na,1];
                           
            # Map 3D parameterization to 2D definition using average thickness              
            else:
                
                # Setup thickness distribution                
                nx = nozzle.stringers.nAxialBreaks;               
 
                # Setup height distribution
                if( nozzle.stringers.heightDefinition == 'EXTERIOR' or \
                    nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT' ):

                    if nozzle.stringers.n != 2:
                        raise NotImplementedError("Stringers definition with " + \
                        "EXTERIOR or BAFFLES_HEIGHT height definition can only have 2 " + \
                        "stringers located on the top and bottom of nozzle.")

                    nozzle.stringers.thickness = [];
                    ns = nozzle.stringers.n;
                    na = nozzle.stringers.nAxialBreaks;
                    for i in range(ns):
                        # extract 0th and 2nd column of data
                        nozzle.stringers.thickness.append(geometry.PiecewiseLinear( \
                            nozzle.stringers.thicknessNodes[i*na:i*na+na,0::2]));
                        nozzle.stringers.thickness[-1].angle = nozzle.stringers.thicknessNodes[i*na,1];

                    # There is no height to assign as it is defined by the 
                    # exterior. Volume, mass, and structural FEA are taken 
                    # care of accordingly.
                    
                    # hMean = np.zeros((nx,));
                    # xTmp = np.linspace(nozzle.stringers.thickness.xstart,nozzle.stringers.thickness.xend,100)
                    # hTop = nozzle.exterior.geometry['top'].coord(xTmp,np.pi/2)[1]
                    # hBot = nozzle.exterior.geometry['bottom'].coord(xTmp,-np.pi/2)[1]
                    # hMean = (np.abs(hTop)+np.abs(hBot))/2
                    # heightNodes = np.transpose(np.vstack((xTmp,hMean)))
                    # nozzle.stringers.height = geometry.PiecewiseLinear(heightNodes);

                else:

                    # Find mean thickness at each x-coordinate station
                    xStations = nozzle.stringers.thicknessNodes[0:nx,0];
                    
                    tMean = np.zeros((nx,));
                    for j in range(nx):
                        tMean[j] = np.mean(list(nozzle.stringers.thicknessNodes[:,2])[j::nx]);
                    thicknessNodes = np.transpose(np.array([xStations,tMean]));
                    
                    nozzle.stringers.thickness = geometry.PiecewiseLinear(thicknessNodes);

                    hMean = np.zeros((nx,));
                    for j in range(nx):
                        hMean[j] = np.mean(list(nozzle.stringers.heightNodes[:,2])[j::nx]);
                    heightNodes = np.transpose(np.array([xStations,hMean]));
                
                    nozzle.stringers.height = geometry.PiecewiseLinear(heightNodes);
        
        else:
            
            # upscale 2D param to 3D
            if nozzle.dim == '3D':
                
                # Setup thickness distribution
                nozzle.stringers.thickness = [];
                ns = nozzle.stringers.n;
                for i in range(ns):
                    nozzle.stringers.thickness.append(geometry.PiecewiseLinear( \
                           nozzle.stringers.thicknessNodes));
                
                # Setup height distribution
                if( nozzle.stringers.heightDefinition == 'EXTERIOR' or \
                    nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT' ):
                    # There is no height to assign as it is defined by the 
                    # exterior. Volume, mass, and structural FEA are taken 
                    # care of accordingly.
                    if nozzle.stringers.n != 2:
                        raise NotImplementedError("Stringers definition with " + \
                        "EXTERIOR or BAFFLES_HEIGHT height definition can only have 2 " + \
                        "stringers located on the top and bottom of nozzle.")                                
                    pass;
                else:                                                    
                    nozzle.stringers.height = [];
                    for i in range(ns):
                        nozzle.stringers.height.append(geometry.PiecewiseLinear( \
                            nozzle.stringers.heightNodes));
                
            else:
                nozzle.stringers.thickness = geometry.PiecewiseLinear( \
                       nozzle.stringers.thicknessNodes);

                # Setup height distribution
                if( nozzle.stringers.heightDefinition == 'EXTERIOR' or \
                    nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT' ):
                    # There is no height to assign as it is defined by the 
                    # exterior. Volume, mass, and structural FEA are taken 
                    # care of accordingly.
                    if nozzle.stringers.n != 2:
                        raise NotImplementedError("Stringers definition with " + \
                        "EXTERIOR or BAFFLES_HEIGHT height definition can only have 2 " + \
                        "stringers located on the top and bottom of nozzle.")                                
                    pass;
                else:
                    nozzle.stringers.height = geometry.PiecewiseLinear( \
                        nozzle.stringers.heightNodes);
            
        # ====================================================================
        # Previous code for axisymmetric stringers/baffles setup
        # ====================================================================

        # nozzle.exterior = component.AxisymmetricWall('exterior');
        
        # # Useful x vs. r data for updating baffle and stringer heights
        # n = 10000
        # x = np.linspace(nozzle.xinlet,nozzle.xoutlet,n)
        # rList = geometry.layerCoordinatesInGlobalFrame(nozzle,x);
    
        # # Setup external geometry
        # if nozzle.stringers.n > 0:
        #     innerBaffleRadius = np.interp(nozzle.xoutlet,x,rList[-2]);
        # else:
        #     innerBaffleRadius = np.interp(nozzle.xoutlet,x,rList[-1]);   
        # #shapeDefinition = np.transpose(np.array([[nozzle.xinlet, nozzle.xinlet+0.1548, nozzle.xoutlet],
        # #                  [0.4244, 0.4244, innerBaffleRadius + 0.012]]));
        # shapeDefinition = np.transpose(np.array([[nozzle.xinlet, nozzle.xinlet+0.1548, nozzle.xoutlet],
        #                   [0.8244, 0.8244, innerBaffleRadius + 0.012]]));
        # nozzle.exterior.geometry = geometry.PiecewiseLinear(shapeDefinition);

        # # Check for intersection of nozzle wall with external geometry
        # if nozzle.stringers.n > 0:
        #     wallDiff = rList[-2] - nozzle.exterior.geometry.radius(x);
        # else:
        #     wallDiff = rList[-1] - nozzle.exterior.geometry.radius(x);
        # if max(wallDiff) > 0:
        #     print max(wallDiff)
        #     print rList[-2]
        #     print nozzle.exterior.geometry.radius(x);
        #     sys.stderr.write('\n ## ERROR: External nozzle wall ' \
        #       'intersects exterior geometry. Terminating.\n\n');
        #     sys.exit(0);
    
        # # Update height of baffles to coincide with exterior wall shape
        # if( nozzle.baffles.heightSetToExterior ):

        #     nozzle.baffles.height = [];
            
        #     for i in range(len(nozzle.baffles.location)):
        #         loc = nozzle.baffles.location[i];
        #         # inner baffle radius is outside radius of outermost layer
        #         if nozzle.stringers.n > 0:
        #             innerBaffleRadius = np.interp(loc,x,rList[-2]);
        #         else:
        #             innerBaffleRadius = np.interp(loc,x,rList[-1]);
        #         outerBaffleRadius = nozzle.exterior.geometry.radius(loc);
        #         nozzle.baffles.height.append(outerBaffleRadius - innerBaffleRadius);
        #         if output == 'verbose':
        #             sys.stdout.write('Baffle %i height resized to %f\n' % 
        #               (i+1,nozzle.baffles.height[i]));

        # # If stringers are dependent on baffles, re-update stringers
        # if( nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT' ):
        #     nozzle.stringers.height = nozzle.exterior.geometry;
                                    
        # # Update height distribution for stringers
        # if( nozzle.stringers.heightDefinition == 'EXTERIOR' ):
        #         nozzle.stringers.height = nozzle.exterior.geometry

        if output == 'verbose':
            sys.stdout.write('Setup Wall complete\n');
            
            
    def InitializeOutput(self, qoi, value=-1, gradients=None):
    
        nozzle = self;
        
        # --- Define possible QOI names

        # Requiring only nozzle shape info
        shape_scalar = ['MASS', 'MASS_WALL_ONLY', 'VOLUME'];
        
        # Requiring aero analysis
        aero_scalar = ['THRUST', 'WALL_PRES_AVG', 'WALL_TEMP_AVG', 
                       'SU2_RESIDUAL'];
        aero_field = ['WALL_PRESSURE','PRESSURE','VELOCITY'];
        
        # Requiring structural analysis
        struct_scalar = ['KS_TOTAL_STRESS', 'PN_TOTAL_STRESS', 
                     'MAX_TOTAL_STRESS', 'MAX_MECHANICAL_STRESS',
                     'MAX_THERMAL_STRESS', 'MAX_FAILURE_CRITERIA', 
                     'KS_FAILURE_CRITERIA', 'PN_FAILURE_CRITERIA',
                     'FAILURE_CRITERIA'];
        
        # Requiring thermal analysis
        therm_scalar = ['KS_TEMPERATURE','PN_TEMPERATURE','MAX_TEMPERATURE',
                     'KS_TEMP_RATIO','PN_TEMP_RATIO','MAX_TEMP_RATIO',
                     'TEMP_RATIO'];
        therm_field = ['WALL_TEMPERATURE'];

        # Prefixes which can be appended to stress and temperature outputs
        prefix = [];
        for i in range(len(nozzle.wall.layer)):
            prefix.append(nozzle.wall.layer[i].name);
        if nozzle.stringers.n > 0:
            prefix.append('STRINGERS');
        for i in range(nozzle.baffles.n):
            prefix.append('BAFFLE' + str(i+1));
        nozzle.prefixLabels = prefix;
        
        # --- Check that QOI name is acceptable
        generalResponseNames = shape_scalar + aero_scalar + aero_field + \
            struct_scalar + therm_scalar + therm_field;
        specificResponseNames = []
        for p in prefix:
            for item in struct_scalar:
                specificResponseNames.append(p + '_' + item);
        for p in prefix:
            for item in therm_scalar:
                specificResponseNames.append(p + '_' + item);
				
        # Initialize responses and gradients
        if isinstance(qoi,list):
            for item in qoi:
                if (item in generalResponseNames) or (item in specificResponseNames):
                    nozzle.responses[item] = value;
                    nozzle.gradients[item] = gradients;
                else: 
                    sys.stderr.write('  ## ERROR : %s not accepted as output QOI\n\n' % item);
                    sys.exit(1);
        else:
            if (qoi in generalResponseNames) or (qoi in specificResponseNames):
                nozzle.responses[qoi] = value;
                nozzle.gradients[qoi] = gradients;
            else: 
                sys.stderr.write('  ## ERROR : %s not accepted as output QOI\n\n' % qoi);
                sys.exit(1);
        
        #print nozzle.gradients
        #sys.exit(1)
        return;
      
        
    def SetOutput(self, qoi, value=None, gradients=None):
    
        nozzle = self;
        
        if isinstance(qoi,list):
            for i in range(len(qoi)):
            
                if value(i) is not None:
                    if qoi(i) in nozzle.responses:
                        nozzle.responses[qoi[i]] = value[i];
                    else:
                        sys.stderror.write('  ## ERROR : %s not initialized for output QOI response\n\n' % qoi[i]);
                        sys.exit(1);
                
                if gradients[i] is not None:
                    if qoi(i) in nozzle.gradients:
                        nozzle.gradients[qoi[i]] = gradients[i];
                    else:
                        sys.stderror.write('  ## ERROR : %s not initialized for output QOI gradients\n\n' % qoi[i]);
                        sys.exit(1);
        else:

            if value is not None:
                if qoi in nozzle.responses:
                    nozzle.responses[qoi] = value;
                else:
                    sys.stderror.write('  ## ERROR : %s not initialized for output QOI response\n\n' % qoi);
                    sys.exit(1);
            
            if gradients is not None:
                if qoi in nozzle.gradients:
                    nozzle.gradients[qoi] = gradients;
                else:
                    sys.stderror.write('  ## ERROR : %s not initialized for output QOI gradients\n\n' % qoi);
                    sys.exit(1);                    
        return;

        
    def GetOutput(self, qoi):
    
        nozzle = self;
        
        if isinstance(qoi,list):
            responseList = [];
            gradientList = [];
            for i in range(len(qoi)):            
                responseList.append(nozzle.responses[qoi[i]]);
                gradientList.append(nozzle.gradients[qoi[i]]);
            
            return responseList, gradientList;
            
        else:
        
            if nozzle.responses[qoi] is None and nozzle.gradients[qoi] is None:
            
                return None;
                
            elif nozzle.responses[qoi] is None:
            
                return nozzle.gradients;
                
            elif nozzle.gradients[qoi] is None:
            
                return nozzle.responses[qoi];
                
            else:

                return nozzle.responses[qoi], nozzle.gradients[qoi];             
            

    def ParseDV(self, config, output='verbose'):
        
        nozzle = self;
        
        if 'INPUT_DV_FORMAT' in config:
            nozzle.inputDVformat = config['INPUT_DV_FORMAT'];
        else :
            sys.stderr.write('\n ## ERROR : Input DV file format not '        \
              'specified. (INPUT_DV_FORMAT expected: PLAIN or DAKOTA)\n\n');
            sys.exit(0);        
                
        if 'INPUT_DV_NAME' in config:
            nozzle.inputDVfilename = config['INPUT_DV_NAME'];
        else :
            sys.stderr.write('\n ## ERROR : Input DV file name not '          \
              'specified. (INPUT_DV_NAME expected)\n\n');
            sys.exit(0);

        #if nozzle.inputDVformat == 'PLAIN' or nozzle.inputDVformat == 'DAKOTA' or :
        nozzle.ParseDesignVariables(nozzle);
		
        for i in range(0, len(nozzle.outputTags)):     
            
            tag  = nozzle.outputTags[i];
            
            # 6 Get Hessian, gradient, and value Get Hessian and gradient
            # 5 Get Hessian and value
            # 4 Get Hessian
            # 3 Get gradient and value
            # 2 Get gradient
            # 1 Get value
            # 0 No data required, function is inactive
            code = nozzle.outputCode[i];
            
            if code == 1: # Output value only
                nozzle.SetOutput(tag,value=-1,gradients=None);
            elif code == 2: # Output gradients only
                nozzle.SetOutput(tag,value=None,gradients=[]);
            elif code == 3: # Output value and gradients
                nozzle.SetOutput(tag,value=-1,gradients=[]);

            else:
                # code = 0: no data required, function is inactive
                # code = 4: get Hessian
                # code = 5: get Hessian and value
                # code = 6: get Hessian, gradient, and value
                sys.stderr.write('\n ## ERROR : code %i in DV input file not available\n' % code);
                sys.exit(1);
        
        #if output == 'verbose':
        #    sys.stdout.write('Setup Parse Design Variables complete\n');    
        #    sys.stdout.write('Design variable list:\n');
        #    print nozzle.dvList
        #    sys.stdout.write('Output codes:\n');
        #    print nozzle.outputCode
        #    sys.stdout.write('Derivative design variables:\n');
        #    print nozzle.derivativesDV       
    
    
    def ParseDesignVariables(self, output='verbose'):
    	
    	nozzle = self;
    	
    	inputDVformat = nozzle.inputDVformat;
    	filename = nozzle.inputDVfilename;
    	
    	if inputDVformat == 'PLAIN':
    		dvList, outputCode, derivativesDV = ParseDesignVariables_Plain(filename); 
    		# Must assign outputCode separately here
    		if nozzle.gradientsFlag == 1:
    		    for i in range(len(nozzle.outputTags)):
    		        outputCode.append(3); # output function value and gradient
    		else:
    		    for i in range(len(nozzle.outputTags)):
    		        outputCode.append(1); # output function value only   
    		NbrDV = len(dvList);
    	elif inputDVformat == 'DAKOTA' :
    		dvList, outputCode, derivativesDV = ParseDesignVariables_Dakota(filename);    
    		NbrDV = len(dvList);
    	elif inputDVformat == 'SAMPLES' :
    		dvList = [];
    		outputCode = [];
    		derivativesDV = [];
    		
    		for i in range(nozzle.NbrDVTot):
				#id_dv = nozzle.DV_Head[iTag] + nozzle.wall.dv[i];
    			dvList.append(-1.0);
    		NbrDV = nozzle.NbrDVTot;
    		
    		for i in range(len(nozzle.outputTags)):
    		    outputCode.append(1); # output function value
    	
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
        
        nozzle.dvList = dvList;
        nozzle.outputCode = outputCode;
        nozzle.derivativesDV = derivativesDV;
    
    
    def UpdateDV(self, output='verbose'):
        
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
                
                if nozzle.param == '3D':

                    # Loop through all possible nozzle wall design variables
                    for i in range(len(nozzle.wall.dv)):
                        # Find placement of design variable in input list
                        id_dv = nozzle.DV_Head[iTag] + nozzle.wall.dv[i];
                        # Update B-spline coefficient if required
                        if id_dv >= nozzle.DV_Head[iTag]:
                            # Centerline coef needs updating (WALL_COEFS1)
                            if i < nozzle.wall.centerline.coefs_size:
                                prt_name.append('Centerline Bspline coef #%d' % (i+1));
                                prt_basval.append('%.4lf'%nozzle.wall.centerline.coefs[i]);
                                prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                                nozzle.wall.centerline.coefs[i] = nozzle.dvList[id_dv];
                                NbrChanged = NbrChanged+1;
                            # Major axis coef needs updating (WALL_COEFS2)
                            elif i < nozzle.wall.centerline.coefs_size + \
                                     nozzle.wall.majoraxis.coefs_size:
                                iLocal = i - nozzle.wall.centerline.coefs_size;
                                prt_name.append('Major axis Bspline coef #%d' % (iLocal+1));
                                prt_basval.append('%.4lf'%nozzle.wall.majoraxis.coefs[iLocal]);
                                prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                                nozzle.wall.majoraxis.coefs[iLocal] = nozzle.dvList[id_dv];
                                NbrChanged = NbrChanged+1;
                            else: # Minor axis coef needs updating (WALL_COEFS3)
                                iLocal = i - nozzle.wall.centerline.coefs_size - nozzle.wall.majoraxis.coefs_size;
                                prt_name.append('Minor axis Bspline coef #%d' % (iLocal+1));
                                prt_basval.append('%.4lf'%nozzle.wall.minoraxis.coefs[iLocal]);
                                prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                                nozzle.wall.minoraxis.coefs[iLocal] = nozzle.dvList[id_dv];
                                NbrChanged = NbrChanged+1;
                    continue;
                
                else: # assume 2D parameterization
                    # Loop through all possible nozzle wall design variables
                    for iCoef in range(len(nozzle.wall.dv)):
                        # Find placement of design variable in input list
                        id_dv = nozzle.DV_Head[iTag] + nozzle.wall.dv[iCoef];                    
                        # Update B-spline coefficient if required
                        if id_dv >= nozzle.DV_Head[iTag]:
                            prt_name.append('Bspline coef #%d' % (iCoef+1));
                            prt_basval.append('%.4lf'% nozzle.wall.coefs[iCoef]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.wall.coefs[iCoef] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;
                    continue;
                    
            if Tag == 'WALL_SHOVEL_HEIGHT':
                id_dv = nozzle.DV_Head[iTag];
                prt_name.append('Wall exit shovel height');
                prt_basval.append('%.4lf'%nozzle.wall.shovel_height);
                prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                nozzle.wall.shovel_height = nozzle.dvList[id_dv];
                NbrChanged = NbrChanged+1;
                continue;
            
            if Tag == 'WALL_SHOVEL_START_ANGLE':
                id_dv = nozzle.DV_Head[iTag];
                prt_name.append('Wall inlet shovel start angle');
                prt_basval.append('%.4lf'%nozzle.wall.shovel_start_angle);
                prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                nozzle.wall.shovel_start_angle = nozzle.dvList[id_dv];
                NbrChanged = NbrChanged+1;
                continue;                
            
            # Update general piecewise linear wall temp definition
            if Tag == 'WALL_TEMP':
                lsize = nozzle.wall.temperature.nBreaks;
                brk = np.max(nozzle.wall.temperature.dv[:lsize])+1;
                for iCoord in range(len(nozzle.wall.temperature.dv)):
                    id_dv = nozzle.DV_Head[iTag] + nozzle.wall.temperature.dv[iCoord];
                    # Update coordinate in thickness array if required
                    if id_dv < nozzle.DV_Head[iTag]:
                        pass
                    elif id_dv < nozzle.DV_Head[iTag]+brk:
                        prt_name.append('wall temp location #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodesNonDim[iCoord,0]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.wall.temperature.thicknessNodesNonDim[iCoord,0] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;
                    else: # id_dv > nozzle.DV_Head[iTag]+brk
                        prt_name.append('wall temp value #%d' % (iCoord+1-lsize));
                        prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodesNonDim[iCoord-lsize,1]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.wall.temperature.thicknessNodesNonDim[iCoord-lsize,1] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;
                continue;
                
            # Update specific thickness or values for piecewise linear wall temp def
            check = 0;
            if Tag == 'WALL_TEMP_LOCATIONS':
                for iCoord in range(nozzle.wall.temperature.nBreaks):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('wall temp location #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodesNonDim[iCoord,0]);
                    prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                    nozzle.wall.temperature.thicknessNodesNonDim[iCoord,0] = nozzle.dvList[id_dv];
                    NbrChanged = NbrChanged+1;
                check = 1;
            elif Tag == 'WALL_TEMP_VALUES':
                for iCoord in range(nozzle.wall.temperature.nBreaks):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                    prt_name.append('wall temp value #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.wall.temperature.thicknessNodesNonDim[iCoord,1]);
                    prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                    nozzle.wall.temperature.thicknessNodesNonDim[iCoord,1] = nozzle.dvList[id_dv];
                    NbrChanged = NbrChanged+1;
                check = 1;
            if check == 1:
                continue;                    

            # Update all layers with non-specific names, e.g. LAYER1, etc.
            check = 0;
            for j in range(len(nozzle.wall.layer)):
                if Tag == nozzle.wall.layer[j].name:
                    
                    if nozzle.param == '3D':
                        # Loop through all possible nozzle wall layer design vars
                        for i in range(len(nozzle.wall.layer[j].dv)):
                            # Find placement of design variable in input list
                            id_dv = nozzle.DV_Head[iTag] + nozzle.wall.layer[j].dv[i];
                            # Update design variable value if required
                            if id_dv >= nozzle.DV_Head[iTag]:
                                # Update axial break (THICKNESS_LOCATIONS)
                                if i < nozzle.wall.layer[j].nAxialBreaks:
                                    prt_name.append('%s thickness location #%d' % (nozzle.wall.layer[j].name,i+1));
                                    prt_basval.append('%.4lf'%nozzle.wall.layer[j].thicknessNodesNonDim[i,0]);
                                    prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                                    # Update all axial break locations
                                    for k in range(i,nozzle.wall.layer[j].nBreaks,nozzle.wall.layer[j].nAxialBreaks):
                                        nozzle.wall.layer[j].thicknessNodesNonDim[k,0] = nozzle.dvList[id_dv];
                                    NbrChanged = NbrChanged+1;
                                # Update angular break (THICKNESS_ANGLES)
                                elif i < nozzle.wall.layer[j].nAxialBreaks + \
                                         nozzle.wall.layer[j].nAngularBreaks:
                                    iLocal = i - nozzle.wall.layer[j].nAxialBreaks;
                                    prt_name.append('%s thickness angle #%d' % (nozzle.wall.layer[j].name,iLocal+1));
                                    prt_basval.append('%.4lf'%nozzle.wall.layer[j].thicknessNodesNonDim[iLocal,1]);
                                    prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                                    # Update all angular break locations
                                    angBreakLoc = iLocal*nozzle.wall.layer[j].nAxialBreaks;
                                    for k in range(angBreakLoc,angBreakLoc+nozzle.wall.layer[j].nAngularBreaks):
                                        nozzle.wall.layer[j].thicknessNodesNonDim[k,1] = nozzle.dvList[id_dv];
                                    NbrChanged = NbrChanged+1;                                    
                                # Update thickness value (THICKNESS_VALUES)
                                else:
                                    iLocal = i - nozzle.wall.layer[j].nAxialBreaks - \
                                      nozzle.wall.layer[j].nAngularBreaks;
                                    prt_name.append('%s thickness value #%d' % (nozzle.wall.layer[j].name,iLocal+1));
                                    prt_basval.append('%.4lf'%nozzle.wall.layer[j].thicknessNodesNonDim[iLocal,2]);
                                    prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                                    # Update thickness value
                                    nozzle.wall.layer[j].thicknessNodesNonDim[iLocal,2] = nozzle.dvList[id_dv];
                                    NbrChanged = NbrChanged+1;   
                        check = 1;

                    else: # assume 2D parameterization
                        
                        lsize = nozzle.wall.layer[j].nBreaks;
                        brk = np.max(nozzle.wall.layer[j].dv[:lsize])+1;
                        for iCoord in range(len(nozzle.wall.layer[j].dv)):
                            id_dv = nozzle.DV_Head[iTag] + nozzle.wall.layer[j].dv[iCoord];
                            # Update coordinate in thickness array if required
                            if id_dv < nozzle.DV_Head[iTag]:
                                pass
                            elif id_dv < nozzle.DV_Head[iTag]+brk:
                                prt_name.append('%s thickness location #%d' % (nozzle.wall.layer[j].name,iCoord+1));
                                prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodesNonDim[iCoord,0]);
                                prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                                nozzle.wall.layer[j].thicknessNodesNonDim[iCoord,0] = nozzle.dvList[id_dv];
                                NbrChanged = NbrChanged+1;
                            else: # id_dv > nozzle.DV_Head[iTag]+brk
                                prt_name.append('%s thickness value #%d' % (nozzle.wall.layer[j].name,iCoord+1-lsize));
                                prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodesNonDim[iCoord-lsize,1]);
                                prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                                nozzle.wall.layer[j].thicknessNodesNonDim[iCoord-lsize,1] = nozzle.dvList[id_dv];
                                NbrChanged = NbrChanged+1;
                        check = 1;
                        
            if check == 1:
                continue;
                
            # Update all piecewise-(bi)linear layers with specific names, e.g. LAYER1_THICKNESS_VALUES
            check = 0;
            for j in range(len(nozzle.wall.layer)):
                
                if nozzle.param == '3D':
                    
                    if Tag == nozzle.wall.layer[j].name + '_THICKNESS_LOCATIONS':
                        
                        # Update axial break (THICKNESS_LOCATIONS)
                        for i in range(nozzle.wall.layer[j].nAxialBreaks):
                            id_dv = nozzle.DV_Head[iTag] + i;
                            prt_name.append('%s thickness location #%d' % (nozzle.wall.layer[j].name,i+1));
                            prt_basval.append('%.4lf'%nozzle.wall.layer[j].thicknessNodesNonDim[i,0]);
                            prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                            # Update all axial break locations
                            for k in range(i,nozzle.wall.layer[j].nBreaks,nozzle.wall.layer[j].nAxialBreaks):
                                nozzle.wall.layer[j].thicknessNodesNonDim[k,0] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;                        
                        check = 1;
                        
                    elif Tag == nozzle.wall.layer[j].name + '_THICKNESS_ANGLES':
                        
                        # Update angular break (THICKNESS_ANGLES)
                        for i in range(nozzle.wall.layer[j].nAngularBreaks):
                            id_dv = nozzle.DV_Head[iTag] + i;
                            prt_name.append('%s thickness angle #%d' % (nozzle.wall.layer[j].name,i));
                            prt_basval.append('%.4lf'%nozzle.wall.layer[j].thicknessNodesNonDim[i,1]);
                            prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                            # Update all angular break locations
                            angBreakLoc = i*nozzle.wall.layer[j].nAxialBreaks;
                            for k in range(angBreakLoc,angBreakLoc+nozzle.wall.layer[j].nAngularBreaks):
                                nozzle.wall.layer[j].thicknessNodesNonDim[k,1] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;                                    
                        check = 1;
                        
                    elif Tag == nozzle.wall.layer[j].name + '_THICKNESS_VALUES':
                        
                        # Update thickness values (THICKNESS_VALUES)
                        for i in range(nozzle.wall.layer[j].nBreaks):
                            id_dv = nozzle.DV_Head[iTag] + i;
                            prt_name.append('%s thickness value #%d' % (nozzle.wall.layer[j].name,i+1));
                            prt_basval.append('%.4lf'%nozzle.wall.layer[j].thicknessNodesNonDim[i,2]);
                            prt_newval.append('%.4lf'%nozzle.dvList[id_dv]);
                            # Update thickness value
                            nozzle.wall.layer[j].thicknessNodesNonDim[i,2] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;   
                        check = 1;
                
                else: # assume 2D parameterization
                
                    if Tag == nozzle.wall.layer[j].name + '_THICKNESS_LOCATIONS':
                        for iCoord in range(nozzle.wall.layer[j].nBreaks):
                            id_dv = nozzle.DV_Head[iTag] + iCoord;
                            prt_name.append('%s thickness location #%d' % (nozzle.wall.layer[j].name,iCoord+1));
                            prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodesNonDim[iCoord,0]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.wall.layer[j].thicknessNodesNonDim[iCoord,0] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;
                        check = 1;
                    elif Tag == nozzle.wall.layer[j].name + '_THICKNESS_VALUES':
                        for iCoord in range(nozzle.wall.layer[j].nBreaks):
                            id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                            prt_name.append('%s thickness value #%d' % (nozzle.wall.layer[j].name,iCoord+1));
                            prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodesNonDim[iCoord,1]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.wall.layer[j].thicknessNodesNonDim[iCoord,1] = nozzle.dvList[id_dv];
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
                    prt_basval.append('%.4lf'% nozzle.wall.layer[j].thicknessNodesNonDim[0,1]);
                    prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                    nozzle.wall.layer[j].thicknessNodesNonDim[0,1] = nozzle.dvList[id_dv];
                    nozzle.wall.layer[j].thicknessNodesNonDim[1,1] = nozzle.dvList[id_dv];
                    nozzle.wall.layer[j].thicknessValue = nozzle.dvList[id_dv];
                    NbrChanged = NbrChanged+1;
                    check = 1;
            if check == 1:
                continue;

            if Tag == 'BAFFLES':
                lsize = (len(nozzle.baffles.dv)-1)/3;
                brk1 = np.max(nozzle.baffles.dv[:lsize])+1;
                brk2 = np.max(nozzle.baffles.dv[lsize:2*lsize])+1;
                brk3 = np.max(nozzle.baffles.dv[2*lsize:3*lsize])+1;
                for iCoord in range(len(nozzle.baffles.dv)):
                    id_dv = nozzle.DV_Head[iTag] + nozzle.baffles.dv[iCoord];                    
                    # --- Update coordinate in baffle definition if required
                    if id_dv < nozzle.DV_Head[iTag]:
                        pass;
                    elif id_dv < nozzle.DV_Head[iTag] + brk1:
                        prt_name.append('baffle location #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.baffles.locationNonDim[iCoord]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.baffles.locationNonDim[iCoord] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;
                        # If stringer location dependent on baffle location
                        if (nozzle.stringers.locationDefinition == 'BAFFLES_LOCATION'):
                                prt_name.append('stringer location #%d' % (iCoord+1));
                                prt_basval.append('%.4lf'% nozzle.baffles.locationNonDim[iCoord]);
                                prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                                nozzle.stringers.thicknessNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];
                                if (nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT'):
                                    nozzle.stringers.heightNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];                        
                                NbrChanged = NbrChanged+1;
                    elif id_dv < nozzle.DV_Head[iTag] + brk2:
                        prt_name.append('baffle thickness #%d' % (iCoord+1-lsize));
                        prt_basval.append('%.4lf'% nozzle.baffles.thickness[iCoord-lsize]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.baffles.thickness[iCoord-lsize] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;                        
                    elif id_dv < nozzle.DV_Head[iTag] + brk3:
                        prt_name.append('baffle height #%d' % (iCoord+1-2*lsize));
                        prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord-2*lsize]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.baffles.height[iCoord-2*lsize] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;   
                        # If stringer height depends on baffle height
                        if (nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT'):
                                prt_name.append('stringer height #%d' % (iCoord+1-2*lsize));
                                prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord-2*lsize]);
                                prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                                nozzle.stringers.heightNodesNonDim[iCoord-2*lsize][1] = nozzle.dvList[id_dv];                     
                                NbrChanged = NbrChanged+1;     
                    else:
                        prt_name.append('baffle half width');
                        prt_basval.append('%.4lf' % nozzle.baffles.halfWidth);
                        prt_newval.append('%.4lf' % nozzle.dvList[id_dv]);
                        nozzle.baffles.halfWidth = nozzle.dvList[id_dv];  
                        NbrChanged += 1;                 
                continue;
                
            if Tag == 'BAFFLES_LOCATION':
                for iCoord in range(len(nozzle.baffles.locationNonDim)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('baffle location #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.baffles.locationNonDim[iCoord]);
                    prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                    nozzle.baffles.locationNonDim[iCoord] = nozzle.dvList[id_dv];
                    NbrChanged = NbrChanged+1;
                    # If stringer location dependent on baffle location
                    if (nozzle.stringers.locationDefinition == 'BAFFLES_LOCATION'):
                            prt_name.append('stringer location #%d' % (iCoord+1));
                            prt_basval.append('%.4lf'% nozzle.baffles.locationNonDim[iCoord]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.stringers.thicknessNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];
                            nozzle.stringers.heightNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];                        
                            NbrChanged = NbrChanged+1;                    
                continue;
               
            if Tag == 'BAFFLES_HEIGHT':
                for iCoord in range(len(nozzle.baffles.height)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('baffle height #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord]);
                    prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                    nozzle.baffles.height[iCoord] = nozzle.dvList[id_dv];
                    NbrChanged = NbrChanged+1;
                    # If stringer height depends on baffle height
                    if (nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT'):
                            prt_name.append('stringer height #%d' % (iCoord+1));
                            prt_basval.append('%.4lf'% nozzle.baffles.height[iCoord]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.stringers.heightNodesNonDim[iCoord][1] = nozzle.dvList[id_dv];                     
                            NbrChanged = NbrChanged+1;
                continue;
            
            if Tag == 'BAFFLES_THICKNESS':
                for iCoord in range(len(nozzle.baffles.thickness)):
                    id_dv = nozzle.DV_Head[iTag] + iCoord;
                    prt_name.append('baffle thickness #%d' % (iCoord+1));
                    prt_basval.append('%.4lf'% nozzle.baffles.thickness[iCoord]);
                    prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                    nozzle.baffles.thickness[iCoord] = nozzle.dvList[id_dv];
                    NbrChanged = NbrChanged+1;
                continue;

            if Tag == 'BAFFLES_HALF_WIDTH':
                id_dv = nozzle.DV_Head[iTag];
                prt_name.append('baffle half width');
                prt_basval.append('%.4lf' % nozzle.baffles.halfWidth);
                prt_newval.append('%.4lf' % nozzle.dvList[id_dv]);
                nozzle.baffles.halfWidth = nozzle.dvList[id_dv];  
                NbrChanged += 1;  
                
            if Tag == 'STRINGERS':
                
                if nozzle.param == '3D':
                    
                    # Loop through all possible stringer design variables
                    for i in range(len(nozzle.stringers.dv)):
                        # Obtain design variable index in input file
                        id_dv = nozzle.DV_Head[iTag] + nozzle.stringers.dv[i];
                        # Update stringer parameter if required
                        if id_dv >= nozzle.DV_Head[iTag]:
                            # Update stringer break locations
                            if i < nozzle.stringers.nAxialBreaks:
                                prt_name.append('stringer break location #%d' % (i+1));
                                prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[i,0]);
                                prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                                jtmp = i;
                                for j in range(nozzle.stringers.n):
                                    nozzle.stringers.thicknessNodesNonDim[jtmp,0] = nozzle.dvList[id_dv];
                                    jtmp += nozzle.stringers.nAxialBreaks;
                                    NbrChanged = NbrChanged+1;
                            # Update stringer angles
                            elif i < nozzle.stringers.nAxialBreaks + \
                              nozzle.stringers.n:
                                jtmp = (i-nozzle.stringers.nAxialBreaks)*nozzle.stringers.nAxialBreaks;
                                prt_name.append('stringer angle #%d' % (i+1-nozzle.stringers.nAxialBreaks));
                                prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[jtmp,1]);
                                prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                                jtmp = 0;
                                for j in range(nozzle.stringers.n):
                                    nozzle.stringers.thicknessNodesNonDim[jtmp:jtmp+nozzle.stringers.nAxialBreaks,1] = nozzle.dvList[id_dv]*np.ones((nozzle.stringers.nAxialBreaks,));
                                    jtmp += nozzle.stringers.nAxialBreaks;
                                    NbrChanged = NbrChanged+1;
                            # Update stringer heights
                            #elif i < :
                            #    pass;             
                            # Update stringer thicknesses
                            else:
                                jtmp = i - nozzle.stringers.nAxialBreaks - nozzle.stringers.n;
                                prt_name.append('stringer thickness #%d' % (jtmp+1));
                                prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[jtmp,2]);
                                prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                                nozzle.stringers.thicknessNodesNonDim[jtmp,2] = nozzle.dvList[id_dv];
                                NbrChanged += 1;                     
                    
                else: # assume 2D parameterization
                    
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
                            prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[iCoord][0]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.stringers.thicknessNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];
                            nozzle.stringers.heightNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;
                        elif id_dv < nozzle.DV_Head[iTag] + brk2:
                            prt_name.append('stringer height value #%d' % (iCoord+1-lsize));
                            prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[iCoord-lsize][1]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.stringers.heightNodesNonDim[iCoord-lsize][1] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;
                        else: # id_dv > nozzle.DV_Head[iTag]+brk
                            prt_name.append('stringer thickness value #%d' % (iCoord+1-2*lsize));
                            prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[iCoord-2*lsize][1]);
                            prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                            nozzle.stringers.thicknessNodesNonDim[iCoord-2*lsize][1] = nozzle.dvList[id_dv];
                            NbrChanged = NbrChanged+1;
                continue;    
                
            # Update all layers with specific names, e.g. LAYER1_THICKNESS_VALUES
            if Tag == 'STRINGERS_BREAK_LOCATIONS':
                
                if nozzle.param == '3D':
                    sys.stderr.write('\n## ERROR: Use of STRINGERS_BREAK_LOCATIONS ' \
                      'as a design variable when using the 3D parameterization' \
                      ' is currently not implemented.\n\n');
                    sys.exit(0);
                else: # assume 2D parameterization
                    for iCoord in range(nozzle.stringers.nBreaks):
                        id_dv = nozzle.DV_Head[iTag] + iCoord;
                        prt_name.append('stringer break location #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[iCoord][0]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.stringers.thicknessNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];
                        nozzle.stringers.heightNodesNonDim[iCoord][0] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;
                    continue;
            
            if Tag == 'STRINGERS_HEIGHT_VALUES':
                
                if nozzle.param == '3D':
                    sys.stderr.write('\n## ERROR: Use of STRINGERS_HEIGHT_VALUES ' \
                      'as a design variable when using the 3D parameterization' \
                      ' is currently not implemented.\n\n');
                    sys.exit(0);
                else: # assume 2D parameterization
                    for iCoord in range(len(nozzle.stringers.heightNodesNonDim)):
                        id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                        prt_name.append('stringer height value #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.stringers.heightNodesNonDim[iCoord][1]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.stringers.heightNodesNonDim[iCoord][1] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;
                    continue;
                    
            if Tag == 'STRINGERS_THICKNESS_VALUES':

                if nozzle.param == '3D':                    
                    for i in range(nozzle.stringers.nThicknessValues):
                        id_dv = nozzle.DV_Head[iTag] + i;
                        prt_name.append('stringer thickness value #%d' % (i+1));
                        prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[i,2]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.stringers.thicknessNodesNonDim[i,2] = nozzle.dvList[id_dv];
                        NbrChanged = NbrChanged+1;
                    continue;
                        
                else: # assume 2D parameterization
                    for iCoord in range(len(nozzle.stringers.thicknessNodesNonDim)):
                        id_dv = nozzle.DV_Head[iTag] + iCoord;                  
                        prt_name.append('stringer thickness value #%d' % (iCoord+1));
                        prt_basval.append('%.4lf'% nozzle.stringers.thicknessNodesNonDim[iCoord,1]);
                        prt_newval.append('%.4lf'% nozzle.dvList[id_dv]);
                        nozzle.stringers.thicknessNodesNonDim[iCoord,1] = nozzle.dvList[id_dv];
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
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setDensity(nozzle.dvList[id_dv]);

                        prt_name.append('%s thermal conductivity' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalConductivity());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);
                        nozzle.materials[k].setThermalConductivity(nozzle.dvList[id_dv+1]);
                        
                        NbrChanged = NbrChanged+2;
                        
                    elif NbrDV == 7:
                        
                        prt_name.append('%s density' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getDensity());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setDensity(nozzle.dvList[id_dv]);
                        
                        prt_name.append('%s elastic modulus' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getElasticModulus());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);
                        nozzle.materials[k].setElasticModulus(nozzle.dvList[id_dv+1]);

                        prt_name.append('%s poisson ratio' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getPoissonRatio());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+2]);
                        nozzle.materials[k].setPoissonRatio(nozzle.dvList[id_dv+2]);

                        prt_name.append('%s thermal conductivity' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalConductivity());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+3]);
                        nozzle.materials[k].setThermalConductivity(nozzle.dvList[id_dv+3]);

                        prt_name.append('%s thermal expansion coef' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalExpansionCoef());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+4]);
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.dvList[id_dv+4]);
                        
                        prt_name.append('%s failure limit' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getFailureLimit());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+5]);
                        nozzle.materials[k].setFailureType(nozzle.materials[k].failureType,nozzle.dvList[id_dv+5]);

                        prt_name.append('%s max service temperature' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].Tmax);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+6]);
                        nozzle.materials[k].setMaxServiceTemperature(nozzle.dvList[id_dv+6]);                        
                        
                        NbrChanged = NbrChanged+7;
                        
                    elif NbrDV == 15: # XXX This needs to be updated
                        
                        prt_name.append('%s density' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getDensity());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setDensity(nozzle.dvList[id_dv]);

                        ltemp = nozzle.materials[k].getElasticModulus();
                        prt_name.append('%s elastic modulus E1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);
                        prt_name.append('%s elastic modulus E2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+2]);                        
                        nozzle.materials[k].setElasticModulus(nozzle.dvList[id_dv+1:id_dv+3]);
                        
                        prt_name.append('%s shear modulus' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getShearModulus());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+3]);
                        nozzle.materials[k].setShearModulus(nozzle.dvList[id_dv+3]);
                        
                        prt_name.append('%s poisson ratio' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getPoissonRatio());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+4]);
                        nozzle.materials[k].setPoissonRatio(nozzle.dvList[id_dv+4]);
                        
                        ltemp = nozzle.materials[k].getMutualInfluenceCoefs();
                        prt_name.append('%s mutual influence coef u1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+5]);
                        prt_name.append('%s mutual influence coef u2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+6]);                        
                        nozzle.materials[k].setMutualInfluenceCoefs(nozzle.dvList[id_dv+5:id_dv+7]);       
                        
                        ltemp = nozzle.materials[k].getThermalConductivity();
                        prt_name.append('%s thermal conductivity k1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+7]);
                        prt_name.append('%s thermal conductivity k2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+8]);                        
                        nozzle.materials[k].setThermalConductivity(nozzle.dvList[id_dv+7:id_dv+9]);  

                        ltemp = nozzle.materials[k].getThermalExpansionCoef();
                        prt_name.append('%s thermal expansion coef a1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+9]);
                        prt_name.append('%s thermal expansion coef a2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+10]);           
                        prt_name.append('%s thermal expansion coef a12' % k);
                        prt_basval.append('%.2le' % ltemp[2]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+11]);                         
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.dvList[id_dv+9:id_dv+12]);
                        
                        NbrChanged = NbrChanged+12;

                    else:
                        raise RuntimeError('%d design variables not accepted' \
                              ' for assignment for %s material' % (NbrDV,k));
                    
                    check = 1;
                              
                elif Tag == k + '_DENSITY':
                    prt_name.append('%s density' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].getDensity());
                    prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                    nozzle.materials[k].setDensity(nozzle.dvList[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;
                elif Tag == k + '_ELASTIC_MODULUS':
                    if NbrDV == 1:
                        prt_name.append('%s elastic modulus' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getElasticModulus());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setElasticModulus(nozzle.dvList[id_dv]);
                        NbrChanged = NbrChanged + 1;
                    elif NbrDV == 2:
                        ltemp = nozzle.materials[k].getElasticModulus();
                        prt_name.append('%s elastic modulus E1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s elastic modulus E2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);                        
                        nozzle.materials[k].setElasticModulus(nozzle.dvList[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                    check = 1;                        
                elif Tag == k + '_SHEAR_MODULUS':
                    prt_name.append('%s shear modulus' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].getShearModulus());
                    prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                    nozzle.materials[k].setShearModulus(nozzle.dvList[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;
                elif Tag == k + '_POISSON_RATIO':
                    prt_name.append('%s poisson ratio' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].getPoissonRatio());
                    prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                    nozzle.materials[k].setPoissonRatio(nozzle.dvList[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;
                elif Tag == k + '_MUTUAL_INFLUENCE_COEFS':
                    ltemp = nozzle.materials[k].getMutualInfluenceCoefs();
                    prt_name.append('%s mutual influence coef u1' % k);
                    prt_basval.append('%.2le' % ltemp[0]);
                    prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                    prt_name.append('%s mutual influence coef u2' % k);
                    prt_basval.append('%.2le' % ltemp[1]);
                    prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]); 
                    if NbrDV == 1:
                        nozzle.materials[k].setMutualInfluenceCoefs(nozzle.dvList[id_dv]);  
                    elif NbrDV == 2:                                               
                        nozzle.materials[k].setMutualInfluenceCoefs(nozzle.dvList[id_dv:id_dv+2]);  
                    NbrChanged = NbrChanged + 2;
                    check = 1;
                elif Tag == k + '_THERMAL_CONDUCTIVITY':
                    if NbrDV == 1:
                        prt_name.append('%s thermal conductivity' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalConductivity());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setThermalConductivity(nozzle.dvList[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        ltemp = nozzle.materials[k].getThermalConductivity();
                        prt_name.append('%s thermal conductivity k1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s thermal conductivity k2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);                        
                        nozzle.materials[k].setThermalConductivity(nozzle.dvList[id_dv:id_dv+2]);    
                        NbrChanged = NbrChanged + 2;
                        check = 1;
                    elif NbrDV == 3:
                        ltemp = nozzle.materials[k].getThermalConductivity();
                        prt_name.append('%s thermal conductivity k1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s thermal conductivity k2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);          
                        prt_name.append('%s thermal conductivity k3' % k);
                        prt_basval.append('%.2le' % ltemp[2]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+2]);                                                
                        nozzle.materials[k].setThermalConductivity(nozzle.dvList[id_dv:id_dv+3]);    
                        NbrChanged = NbrChanged + 3;
                        check = 1;                       
                elif Tag == k + '_THERMAL_EXPANSION_COEF':
                    if NbrDV == 1:
                        prt_name.append('%s thermal expansion coef' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].getThermalExpansionCoef());
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.dvList[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 3:
                        ltemp = nozzle.materials[k].getThermalExpansionCoef();
                        prt_name.append('%s thermal expansion coef a1' % k);
                        prt_basval.append('%.2le' % ltemp[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s thermal expansion coef a2' % k);
                        prt_basval.append('%.2le' % ltemp[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);           
                        prt_name.append('%s thermal expansion coef a12' % k);
                        prt_basval.append('%.2le' % ltemp[2]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+2]);                         
                        nozzle.materials[k].setThermalExpansionCoef(nozzle.dvList[id_dv:id_dv+3]);
                        NbrChanged = NbrChanged + 3;
                        check = 1;
                elif Tag == k + '_PRINCIPLE_FAILURE_STRAIN':
                    if NbrDV == 1:
                        prt_name.append('%s principal failure strain' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setFailureType('PRINCIPLE_FAILURE_STRAIN',nozzle.dvList[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        prt_name.append('%s principal failure strain (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s principal failure strain (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);                        
                        nozzle.materials[k].setFailureType('PRINCIPLE_FAILURE_STRAIN',nozzle.dvList[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                        check = 1;                        
                    elif NbrDV == 5:
                        prt_name.append('%s principal failure strain e1 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s principal failure strain e1 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);   
                        prt_name.append('%s principal failure strain e2 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[2]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+2]);
                        prt_name.append('%s principal failure strain e2 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[3]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+3]);   
                        prt_name.append('%s principal failure strain e12' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[4]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+4]);                          
                        nozzle.materials[k].setFailureType('PRINCIPLE_FAILURE_STRAIN',nozzle.dvList[id_dv:id_dv+5]);
                        NbrChanged = NbrChanged + 5;
                        check = 1; 
                elif Tag == k + '_LOCAL_FAILURE_STRAIN':
                    if NbrDV == 1:
                        prt_name.append('%s local failure strain' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setFailureType('LOCAL_FAILURE_STRAIN',nozzle.dvList[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        prt_name.append('%s local failure strain (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s local failure strain (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);                        
                        nozzle.materials[k].setFailureType('LOCAL_FAILURE_STRAIN',nozzle.dvList[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                        check = 1;                        
                    elif NbrDV == 5:
                        prt_name.append('%s local failure strain e1 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s local failure strain e1 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);   
                        prt_name.append('%s local failure strain e2 (tension)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[2]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+2]);
                        prt_name.append('%s local failure strain e2 (compression)' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[3]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+3]);   
                        prt_name.append('%s local failure strain e12' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].failureStrain[4]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+4]);                          
                        nozzle.materials[k].setFailureType('LOCAL_FAILURE_STRAIN',nozzle.dvList[id_dv:id_dv+5]);
                        NbrChanged = NbrChanged + 5;
                        check = 1; 
                elif Tag == k + '_YIELD_STRESS':
                    if NbrDV == 1:
                        prt_name.append('%s yield stress' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].yieldStress);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        nozzle.materials[k].setFailureType('VON_MISES',nozzle.dvList[id_dv]);
                        NbrChanged = NbrChanged + 1;
                        check = 1;
                    elif NbrDV == 2:
                        prt_name.append('%s yield stress direction 1' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].yieldStress[0]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                        prt_name.append('%s yield stress direction 2' % k);
                        prt_basval.append('%.2le' % nozzle.materials[k].yieldStress[1]);
                        prt_newval.append('%.2le'% nozzle.dvList[id_dv+1]);                        
                        nozzle.materials[k].setFailureType('VON_MISES',nozzle.dvList[id_dv:id_dv+2]);
                        NbrChanged = NbrChanged + 2;
                        check = 1;                
                elif Tag == k + '_MAX_SERVICE_TEMPERATURE':
                    prt_name.append('%s max service temp' % k);
                    prt_basval.append('%.2le' % nozzle.materials[k].Tmax);
                    prt_newval.append('%.2le'% nozzle.dvList[id_dv]);
                    nozzle.materials[k].setMaxServiceTemperature(nozzle.dvList[id_dv]);
                    NbrChanged = NbrChanged + 1;
                    check = 1;            
            if check == 1:
                continue;

            if Tag == 'INLET_PSTAG':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Inlet stagnation pressure');
                prt_basval.append('%.2lf'% nozzle.inlet.Pstag);
                prt_newval.append('%.2lf'% nozzle.dvList[id_dv]);                
                nozzle.inlet.setPstag(nozzle.dvList[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;
                
            if Tag == 'INLET_TSTAG':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Inlet stagnation temperature');
                prt_basval.append('%.2lf'% nozzle.inlet.Tstag);
                prt_newval.append('%.2lf'% nozzle.dvList[id_dv]);                
                nozzle.inlet.setTstag(nozzle.dvList[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;

            if Tag == 'ATM_PRES':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Atmospheric pressure');
                prt_basval.append('%.2lf'% nozzle.environment.P);
                prt_newval.append('%.2lf'% nozzle.dvList[id_dv]);                
                nozzle.environment.setPressure(nozzle.dvList[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;
                
            if Tag == 'ATM_TEMP':                
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Atmospheric temperature');
                prt_basval.append('%.2lf'% nozzle.environment.T);
                prt_newval.append('%.2lf'% nozzle.dvList[id_dv]);                
                nozzle.environment.setTemperature(nozzle.dvList[id_dv]);
                NbrChanged = NbrChanged + 1;
                continue;

            if Tag == 'HEAT_XFER_COEF_TO_ENV':
                id_dv = nozzle.DV_Head[iTag];                
                prt_name.append('Heat xfer coef. to environment');
                prt_basval.append('%.2lf'% nozzle.environment.hInf);
                prt_newval.append('%.2lf'% nozzle.dvList[id_dv]);                
                nozzle.environment.setHeatTransferCoefficient(nozzle.dvList[id_dv]);
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

            
    def SetupResponsesAndGradients (self, config, output='verbose'):
        
        nozzle = self;
        
        nozzle.outputTags = []; # record output tags in config file
        
        nozzle.outputLocations = dict(); # for field variable outputs
        
        if 'OUTPUT_NAME' in config:
            nozzle.outputFile = config['OUTPUT_NAME'];
        else :
            sys.stderr.write("\n ## ERROR : Output function file name not "   \
              "specified in config file. (outputName expected)\n\n");
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
        
        outputTags = []
        if 'OUTPUT_FUNCTIONS' in config:
            
            # Obtain output functions of interest
            hdl = config['OUTPUT_FUNCTIONS'].strip('()');
            hdl = hdl.split(",");
            
            for i in range(len(hdl)):
                outputTags.append(hdl[i].strip());
            nozzle.outputTags = outputTags;
            
	        # Initialize gradients computation flags	
            nozzle.gradientsFlag = 0;
            nozzle.gradientsMethod = 'FINITE_DIFF';
	        
            # Assign information for gradient calculation if necessary
            if ('OUTPUT_GRADIENTS' in config) and (config['OUTPUT_GRADIENTS'] == 'YES'):

                nozzle.gradientsFlag = 1;
                if 'OUTPUT_GRADIENTS_FILENAME' in config:
                    nozzle.gradientsFile = config['OUTPUT_GRADIENTS_FILENAME'];	
                else:
                    nozzle.gradientsFile = 'gradients.dat';		    

                # Set gradient computation method & required information
                if 'GRADIENTS_COMPUTATION_METHOD' in config:

                    if config['GRADIENTS_COMPUTATION_METHOD'] == 'FINITE_DIFF': 
                        nozzle.gradientsMethod = 'FINITE_DIFF';

                        if 'FD_STEP_SIZE' in config: # absolute forward finite difference step size
                            nozzle.fd_step_size = config['FD_STEP_SIZE'].strip('()');
                            if ',' in nozzle.fd_step_size: # Different step size for all vars
                                nozzle.fd_step_size = nozzle.fd_step_size.split(',');
                                nozzle.fd_step_size = [float(i) for i in nozzle.fd_step_size];
                            else: # Single step size for all variables
                                nozzle.fd_step_size = float(nozzle.fd_step_size);
                        else:
                            sys.stderr.write("  ## ERROR : FD_STEP_SIZE not specified in config file.\n\n");
                            sys.exit(1);

                    elif config['GRADIENTS_COMPUTATION_METHOD'] == 'ADJOINT':
                        nozzle.gradientsMethod = 'ADJOINT';
                        
                        # Some derivatives may need finite differencing
                        if 'FD_STEP_SIZE' in config: # absolute forward finite difference step size
                            nozzle.fd_step_size = config['FD_STEP_SIZE'].strip('()');
                            if ',' in nozzle.fd_step_size: # Different step size for all vars
                                nozzle.fd_step_size = nozzle.fd_step_size.split(',');
                                nozzle.fd_step_size = [float(i) for i in nozzle.fd_step_size];
                            else: # Single step size for all variables
                                nozzle.fd_step_size = float(nozzle.fd_step_size); 
                        else:
                            sys.stdout.write('Even though ADJOINT gradients have been specified, '
                              'specifying FD_STEP_SIZE is a good idea too.');

                    else:
                        sys.stderr.write("  ## ERROR : Invalid entry for GRADIENTS_COMPUTATION_METHOD option (expected: ADJOINT or FINITE_DIFF)\n");
                        sys.exit(1);

                # Initialize function values and gradients 
                nozzle.InitializeOutput(outputTags,value=-1,gradients=None);
	            
            else:

                # Intialize function values only
                nozzle.InitializeOutput(outputTags,value=-1,gradients=None);

        # Initialize output locations for QoI internal to nozzle flow
        for qoi in outputTags:
            if( qoi + '_LOCATIONS' in config ):
                
                if nozzle.dim == '3D':
                    sys.stderr.write('\nWARNING: %s is specified as an output will not be processed for high-fidelity.\n' % qoi);

                loc = config[qoi+'_LOCATIONS'].strip('()')
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
                    nozzle.outputLocations[qoi] = np.squeeze(np.array(locList))                        
                else: # filename given
                    nozzle.outputLocations[qoi] = np.loadtxt(loc[0])
                      
        # Initialize output locations for QoI on structural mesh
        for qoi in outputTags:
            # Nodes are specified
            if( qoi + '_NODES' in config ):
                loc = config[qoi+'_NODES'].strip('()')
                loc = loc.split(',')
                loc = [float(e) for e in loc]
                nozzle.outputLocations[qoi] = np.squeeze(np.array(loc))

        if len(nozzle.outputTags) == 0 :
            sys.stderr.write("\n  ## ERROR : No output function was given.\n\n");
            sys.exit(1);
        
        if output == 'verbose':
            sys.stdout.write('Setup Responses and Gradients complete\n');  


        
    def GetOutputFunctions (self,output='verbose'):
    
        nozzle = self;

        filename = nozzle.outputFile;
        
        # Output arrays
        tag_out = list();
        val_out = list();
        gra_out = list();
        gratag_out = list();
        
        # Stress and temperature outputs can have prefixes appended to them
        prefix = [];
        for i in range(len(nozzle.wall.layer)):
            prefix.append(nozzle.wall.layer[i].name);
        if nozzle.stringers.n > 0:
            prefix.append('STRINGERS');
        for i in range(nozzle.baffles.n):
            prefix.append('BAFFLE' + str(i+1)); 

        # Print function values first
        for i in range(0, len(nozzle.outputTags)):
            
            tag = nozzle.outputTags[i];  
			
            # 6 Get Hessian, gradient, and value Get Hessian and gradient
            # 5 Get Hessian and value
            # 4 Get Hessian
            # 3 Get gradient and value
            # 2 Get gradient
            # 1 Get value
            # 0 No data required, function is inactive
            code = nozzle.outputCode[i];
                            
            if code == 1 or code == 0:
            	pass;
            elif code == 2:
            	continue; # skip writing of function value to file
            elif code == 3:
            	pass;
            else:
                sys.stderr.write('\n ## ERROR : code %i in DV input file not available\n' % code);
                sys.exit(1);

            # Write response values to file
            if code == 1 or code == 3 or code == 5:
            
                if isinstance(nozzle.responses[tag],list):
                    for i in range(len(nozzle.responses[tag])):
                        if isinstance(nozzle.responses[tag][i],list): # i.e. nested list
                            for j in range(len(nozzle.responses[tag][i])):
                                
                                #fil.write('%0.16f %s_%i_%i\n' % (nozzle.responses[tag][i][j],tag,i,j));                                
                                tag_out.append("%s_%i_%i" % (tag,i,j));
                                val_out.append(nozzle.responses[tag][i][j]);
                                
                        else:
                            #fil.write('%0.16f %s_%i\n' % (nozzle.responses[tag][i],tag,i));
                            
                            tag_out.append("%s_%i" % (tag,i));
                            val_out.append(nozzle.responses[tag][i]);
                            
                elif isinstance(nozzle.responses[tag],np.ndarray):
                    arrayShape = nozzle.responses[tag].shape;
                    if len(arrayShape) == 1:
                        nr = arrayShape[0];
                        for i in range(nr):
                            #fil.write('%0.16f %s_%i\n' % (nozzle.responses[tag][i],tag,i));
                            tag_out.append("%s_%i" % (tag,i));
                            val_out.append(nozzle.responses[tag][i]);
                            
                    elif len(arrayShape) == 2:
                        nr, nc = arrayShape;      
                        for i in range(nr):
                            for j in range(nc):
                                #fil.write('%0.16f %s_%i_%i\n' % (nozzle.responses[tag][i,j],tag,i,j));
                                tag_out.append("%s_%i_%i" % (tag,i,j));
                                val_out.append(nozzle.responses[tag][i][j]);
                else:
                    #fil.write('%0.16f %s\n' % (nozzle.responses[tag],tag));
                    tag_out.append("%s" % (tag));
                    val_out.append(nozzle.responses[tag]);
                    

        # Print function gradients next
        for i in range(0, len(nozzle.outputTags)):
            
            tag = nozzle.outputTags[i];  
			
            # 6 Get Hessian, gradient, and value Get Hessian and gradient
            # 5 Get Hessian and value
            # 4 Get Hessian
            # 3 Get gradient and value
            # 2 Get gradient
            # 1 Get value
            # 0 No data required, function is inactive
            code = nozzle.outputCode[i];
                            
            if code == 1 or code == 0:
            	pass;
            elif code == 2:
            	continue; # skip writing of function value to file
            elif code == 3:
            	pass;
            else:
                sys.stderr.write('\n ## ERROR : code %i in DV input file not available\n' % code);
                sys.exit(1);
            
            # Write response gradients to file
            if code == 2 or code == 3 or code == 6:
                
                if output == 'verbose':
                    print nozzle.gradients[tag];
				
                if isinstance(nozzle.gradients[tag][0],list):
                
                    sys.stderr.write('Currently outputting gradients for a vector QOI '
                      'is not enabled\n');
                    sys.exit(1);
                        
                else:
                    
                    gra = list();
                    
                    for i in range(len(nozzle.gradients[tag])):
                        #fil.write('%0.16e ' % nozzle.gradients[tag][i]);
                        gra.append(nozzle.gradients[tag][i]);
		
                    gra_out.append(gra);
                    gratag_out.append(tag);
    
        return tag_out, val_out, gra_out, gratag_out;


    def WriteOutputFunctions_Plain (self, output='verbose'):
    
        nozzle = self;

        filename = nozzle.outputFile;
        
        if output == 'verbose':
            sys.stdout.write('\n');
            string = " Post-processing ";
            nch = (60-len(string))/2;
            sys.stdout.write('-' * nch);
            sys.stdout.write(string);
            sys.stdout.write('-' * nch);
            sys.stdout.write('\n\n');
        
        try:
            fil = open(filename, 'w');
        except:
            sys.stderr.write("  ## ERROR : Could not open output file %s\n" % filename);
            sys.exit(1);
  
        if output == 'verbose':
            sys.stdout.write('  -- Info : Output functions file : %s\n' % filename);
              
        # Stress and temperature outputs can have prefixes appended to them
        prefix = [];
        for i in range(len(nozzle.wall.layer)):
            prefix.append(nozzle.wall.layer[i].name);
        if nozzle.stringers.n > 0:
            prefix.append('STRINGERS');
        for i in range(nozzle.baffles.n):
            prefix.append('BAFFLE' + str(i+1)); 

        # For printing to screen
        prt_item = list();
        prt_comp = list();
        prt_val = list(); 

        # Print function values first
        for i in range(0, len(nozzle.outputTags)):
            
            tag = nozzle.outputTags[i];            
            
            # 6 Get Hessian, gradient, and value Get Hessian and gradient
            # 5 Get Hessian and value
            # 4 Get Hessian
            # 3 Get gradient and value
            # 2 Get gradient
            # 1 Get value
            # 0 No data required, function is inactive
            code = nozzle.outputCode[i];
            
            if code == 1:
            	pass;
            elif code == 2:
            	continue; # skip writing of function value to file
            elif code == 3:
            	pass;
            else:
                sys.stderr.write('\n ## ERROR : code %i in DV input file not available\n' % code);
                sys.exit(1);
            
            # Write response values to file
            if code == 1 or code == 3 or code == 5:

                if isinstance(nozzle.responses[tag],list):
                    for i in range(len(nozzle.responses[tag])):
                        if isinstance(nozzle.responses[tag][i],list): # i.e. nested list
                            for j in range(len(nozzle.responses[tag][i])):
                                fil.write('%0.16f %s_%i_%i\n' % (nozzle.responses[tag][i][j],tag,i,j));
                                # Do not print anything to the screen
                        else:
                            fil.write('%0.16f\n' % nozzle.responses[tag][i]);
                            prt_item.append('%s %i' % (nozzle.responses[tag][i],i));
                            prt_comp.append('%s' % nozzle.prefixLabels[i]);
                            prt_val.append('%0.16f' % nozzle.responses[tag][i]);
                elif isinstance(nozzle.responses[tag],np.ndarray):
                    arrayShape = nozzle.responses[tag].shape;
                    if len(arrayShape) == 1:
                        nr = arrayShape[0];
                        for i in range(nr):
                            fil.write('%0.16f\n' % nozzle.responses[tag][i]);
                            prt_item.append('%s loc %i' % (tag,i));
                            prt_comp.append('');
                            prt_val.append('%0.16f' % nozzle.responses[tag][i]);
                    elif len(arrayShape) == 2:
                        nr, nc = arrayShape;      
                        for i in range(nr):
                            for j in range(nc):
                                fil.write('%0.16f\n' % nozzle.responses[tag][i,j]);
                                prt_item.append('%s loc %i %i' % (tag,i,j));
                                prt_comp.append('');
                                prt_val.append('%0.16f' % nozzle.responses[tag][i,j]);                                  
                else:
                    fil.write('%0.16f\n' % nozzle.responses[tag]);
                    prt_item.append('%s' % tag);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.responses[tag]);                
            
            # Write response gradients to file
            if code == 2 or code == 3 or code == 6:
                
                if isinstance(nozzle.gradients[tag][0],list):
                
                    sys.stderr.write('Currently outputting gradients for a vector QOI '
                      'is not enabled\n');
                    sys.exit(1);
                    
                    for i in range(len(nozzle.responses[tag])):
                        fil.write('%0.16f %s_%i\n' % (nozzle.responses[tag][i],tag,i));
                        prt_item.append('%s %i' % (nozzle.responses[tag][i],i));
                        prt_comp.append('%s' % nozzle.prefixLabels[i]);
                        prt_val.append('%0.16f' % nozzle.responses[tag][i]);
                        
                else:
                
                    fil.write('[ ');
                    for i in range(len(nozzle.gradients[tag])):
                        fil.write('%0.16e ' % nozzle.gradients[tag][i]);
                        #prt_item.append('%s %i' % (nozzle.gradients[tag][i],i));
                        #prt_comp.append('%s' % nozzle.prefixLabels[i]);
                        #prt_val.append('%0.16f' % nozzle.gradients[tag][i]);
                    fil.write(']\n');		
                
        # Close output file
        if output == 'verbose':
            sys.stdout.write('\n');
        fil.close();
        
        # Print summary (function values only)
        if output == 'verbose':
            sys.stdout.write('-' * 85);
            sys.stdout.write('\n%s | %s | %s\n' % ("Item".ljust(40), "Component".ljust(20),"Value".ljust(20)));
            sys.stdout.write('-' * 85);
            sys.stdout.write('\n');
            for i in range(0,len(prt_item)):
                sys.stdout.write('%s | %s | %s\n' % (prt_item[i].ljust(40), prt_comp[i].ljust(20),prt_val[i].ljust(20)));
            sys.stdout.write('-' * 85);    
            sys.stdout.write('\n\n');            
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to "verbose" or "quiet" mode')          

        
    def WriteOutputFunctions_Dakota (self,output='verbose'):
    
        nozzle = self;

        filename = nozzle.outputFile;
        
        if output == 'verbose':
            sys.stdout.write('\n');
            string = " Post-processing ";
            nch = (60-len(string))/2;
            sys.stdout.write('-' * nch);
            sys.stdout.write(string);
            sys.stdout.write('-' * nch);
            sys.stdout.write('\n\n');
        
        try:
            fil = open(filename, 'w');
        except:
            sys.stderr.write("  ## ERROR : Could not open output file %s\n" % filename);
            sys.exit(1);
  
        if output == 'verbose':
            sys.stdout.write('  -- Info : Output functions file : %s\n' % filename);
              
        # Stress and temperature outputs can have prefixes appended to them
        prefix = [];
        for i in range(len(nozzle.wall.layer)):
            prefix.append(nozzle.wall.layer[i].name);
        if nozzle.stringers.n > 0:
            prefix.append('STRINGERS');
        for i in range(nozzle.baffles.n):
            prefix.append('BAFFLE' + str(i+1)); 

        # For printing to screen
        prt_item = list();
        prt_comp = list();
        prt_val  = list(); 

        # Print function values first
        for i in range(0, len(nozzle.outputTags)):
            
            tag = nozzle.outputTags[i];  
			
            # 6 Get Hessian, gradient, and value Get Hessian and gradient
            # 5 Get Hessian and value
            # 4 Get Hessian
            # 3 Get gradient and value
            # 2 Get gradient
            # 1 Get value
            # 0 No data required, function is inactive
            code = nozzle.outputCode[i];
            
            if output == 'verbose':
                print "TAG %s CODE %d" % (tag, code);
                
            if code == 1 or code == 0:
            	pass;
            elif code == 2:
            	continue; # skip writing of function value to file
            elif code == 3:
            	pass;
            else:
                sys.stderr.write('\n ## ERROR : code %i in DV input file not available\n' % code);
                sys.exit(1);

            # Write response values to file
            if code == 1 or code == 3 or code == 5:
            
                if isinstance(nozzle.responses[tag],list):
                    for i in range(len(nozzle.responses[tag])):
                        if isinstance(nozzle.responses[tag][i],list): # i.e. nested list
                            for j in range(len(nozzle.responses[tag][i])):
                                fil.write('%0.16f %s_%i_%i\n' % (nozzle.responses[tag][i][j],tag,i,j));
                                # Do not print anything to the screen
                        else:
                            fil.write('%0.16f %s_%i\n' % (nozzle.responses[tag][i],tag,i));
                            prt_item.append('%s %i' % (nozzle.responses[tag][i],i));
                            prt_comp.append('%s' % nozzle.prefixLabels[i]);
                            prt_val.append('%0.16f' % nozzle.responses[tag][i]);
                elif isinstance(nozzle.responses[tag],np.ndarray):
                    arrayShape = nozzle.responses[tag].shape;
                    if len(arrayShape) == 1:
                        nr = arrayShape[0];
                        for i in range(nr):
                            fil.write('%0.16f %s_%i\n' % (nozzle.responses[tag][i],tag,i));
                            prt_item.append('%s loc %i' % (tag,i));
                            prt_comp.append('');
                            prt_val.append('%0.16f' % nozzle.responses[tag][i]);
                    elif len(arrayShape) == 2:
                        nr, nc = arrayShape;      
                        for i in range(nr):
                            for j in range(nc):
                                fil.write('%0.16f %s_%i_%i\n' % (nozzle.responses[tag][i,j],tag,i,j));
                                prt_item.append('%s loc %i %i' % (tag,i,j));
                                prt_comp.append('');
                                prt_val.append('%0.16f' % nozzle.responses[tag][i,j]);                                  
                else:
                    fil.write('%0.16f %s\n' % (nozzle.responses[tag],tag));
                    prt_item.append('%s' % tag);
                    prt_comp.append('');
                    prt_val.append('%0.16f' % nozzle.responses[tag]);                

        # Print function gradients next
        for i in range(0, len(nozzle.outputTags)):
            
            tag = nozzle.outputTags[i];  
			
            # 6 Get Hessian, gradient, and value Get Hessian and gradient
            # 5 Get Hessian and value
            # 4 Get Hessian
            # 3 Get gradient and value
            # 2 Get gradient
            # 1 Get value
            # 0 No data required, function is inactive
            code = nozzle.outputCode[i];
            
            if output == 'verbose':
                print "TAG %s CODE %d" % (tag, code);
                
            if code == 1 or code == 0:
            	pass;
            elif code == 2:
            	continue; # skip writing of function value to file
            elif code == 3:
            	pass;
            else:
                sys.stderr.write('\n ## ERROR : code %i in DV input file not available\n' % code);
                sys.exit(1);
            
            # Write response gradients to file
            if code == 2 or code == 3 or code == 6:
                
                if output == 'verbose':
                    print nozzle.gradients[tag];
				
                if isinstance(nozzle.gradients[tag][0],list):
                
                    sys.stderr.write('Currently outputting gradients for a vector QOI '
                      'is not enabled\n');
                    sys.exit(1);
                    
                    for i in range(len(nozzle.responses[tag])):
                        fil.write('%0.16f %s_%i\n' % (nozzle.responses[tag][i],tag,i));
                        prt_item.append('%s %i' % (nozzle.responses[tag][i],i));
                        prt_comp.append('%s' % nozzle.prefixLabels[i]);
                        prt_val.append('%0.16f' % nozzle.responses[tag][i]);
                        
                else:
                
                    fil.write('[ ');
                    for i in range(len(nozzle.gradients[tag])):
                        fil.write('%0.16e ' % nozzle.gradients[tag][i]);
                        #prt_item.append('%s %i' % (nozzle.gradients[tag][i],i));
                        #prt_comp.append('%s' % nozzle.prefixLabels[i]);
                        #prt_val.append('%0.16f' % nozzle.gradients[tag][i]);
                    fil.write(']\n');		
                
        # Close output file
        if output == 'verbose':
            sys.stdout.write('\n');
        fil.close();
        
        # Print summary (function values only)
        if output == 'verbose':
            sys.stdout.write('-' * 85);
            sys.stdout.write('\n%s | %s | %s\n' % ("Item".ljust(40), "Component".ljust(20),"Value".ljust(20)));
            sys.stdout.write('-' * 85);
            sys.stdout.write('\n');
            for i in range(0,len(prt_item)):
                sys.stdout.write('%s | %s | %s\n' % (prt_item[i].ljust(40), prt_comp[i].ljust(20),prt_val[i].ljust(20)));
            sys.stdout.write('-' * 85);    
            sys.stdout.write('\n\n');       

        
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
        
        
        Box = [[nozzle.xinlet, nozzle.xoutlet],[ymin, ymax]];
                
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


    # Useful function for printing everything about nozzle
    def printNozzleData(self, output='verbose'):

        nozzle = self

        print '\n*****************NOZZLE INFORMATION*****************\n'

        print 'MISSION:'
        print nozzle.mission.__dict__
        print
        
        print 'CFD:'
        print nozzle.cfd.__dict__
        print
        
        print 'FLUID:'
        print nozzle.fluid.__dict__
        print
        
        print 'EXTERIOR ENVIRONMENT:'
        print nozzle.environment.__dict__
        print
        
        print 'INLET:'
        print nozzle.inlet.__dict__
        print
            
        print 'NOZZLE INTERIOR WALL:'
        if nozzle.dim == '3D':
            print 'The interior nozzle wall is parameterized with 3 B-splines.'
            print 'One B-spline parameterizes the centerline, and two more '
            print 'parameterize the major and minor axes of an ellipse centered '
            print 'on the centerline, drawn in the vertical plane. Two other '
            print 'parameters, the shovel exit height and shovel inlet start '
            print 'angle are used to define the flattening out of the geometry '
            print 'near the exit.'
            print 'Nozzle wall centerline geometry:'
            print nozzle.wall.centerline.geometry.__dict__
            print 'Nozzle wall major axis geometry:'
            print nozzle.wall.majoraxis.geometry.__dict__
            print 'Nozzle wall minor axis geometry:'
            print nozzle.wall.minoraxis.geometry.__dict__
            print 'Nozzle wall equivalent axisymmetric geometry:'
            print nozzle.wall.geometry.__dict__
            print 'Nozzle length: %f' % nozzle.length
            print 'Nozzle wall struct:'
            print nozzle.wall.__dict__
            print 'Nozzle wall shovel exit height: %f' % nozzle.wall.shovel_height
            print 'Nozzle wall shovel inlet start angle: %f' % nozzle.wall.shovel_start_angle
            print
        else:
            print 'Nozzle wall geometry:'
            print nozzle.wall.geometry.__dict__
            print 'Nozzle length: %f' % nozzle.length
            print 'Nozzle wall struct:'
            print nozzle.wall.__dict__
            print
        
        print 'WALL LAYERS:'
        if nozzle.dim == '3D':
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

        print 'MATERIALS:'
        print 'Not implemented yet for printing'
        
        print 'STRINGERS:'
        print 'Since the nozzle is nonaxisymmetric, stringers are defined by '
        print 'angular position around the centerline. Thicknesses are a '
        print 'function of global axial X-coordinate and stringer angle. '
        print 'A stringer is assumed to travel along the same angular coord.'
        print nozzle.stringers.__dict__
        if nozzle.dim == '3D':
            for i in range(len(nozzle.stringers.height)):
                print 'Nozzle stringer %i height:' % i
                print nozzle.stringers.height[i].__dict__
                print 'Nozzle stringer %i thickness:' % i
                print nozzle.stringers.thickness[i].__dict__
        
        print 'BAFFLES:'
        print nozzle.baffles.__dict__
        print
        
        print 'EXTERIOR:'
        print nozzle.exterior.__dict__
        print
        
        print 'Entire NOZZLE struct:'
        for k in nozzle.__dict__:
            print k
            print nozzle.__dict__[k]
            print

        return

    # END of Nozzle class.


def NozzleSetup( config, flevel, output='verbose'):
    
    nozzle = multif.nozzle.nozzle.Nozzle();

    # --- Begin setup of CFD structure for CFD related information
    nozzle.cfd = CFD();
    nozzle.cfd.mesh_name = 'nozzle.su2';
    nozzle.cfd.restart_name =  'nozzle.dat';
    nozzle.cfd.exit_mesh_name = 'nozzle_exit.mesh';
    nozzle.cfd.conv_filename = 'history';
    
    # --- General nozzle information
    
    if 'TEMP_RUN_DIR' in config and config['TEMP_RUN_DIR'] == 'YES':
        nozzle.runDir = tempfile.mkdtemp();    
    else:
        nozzle.runDir = '';

	# --- Mesh generation method
	if 'MESH_GENERATION_METHOD' in config : 
		if config['MESH_GENERATION_METHOD'] == 'DEFORM':
			nozzle.meshDeformationFlag = True; 
		elif config['MESH_GENERATION_METHOD'] == 'REGEN':
			nozzle.meshDeformationFlag = False;
		else :
			sys.stderr.write("  ## ERROR : Invalid option for MESH_GENERATION (DEFORM or REGEN)\n");
			sys.exit(1);
	else:
		nozzle.meshDeformationFlag = False;	
	
    # --- Path to SU2 exe
    
    if 'SU2_RUN' in config:
        nozzle.cfd.su2_run = config['SU2_RUN'];
    else:
        nozzle.cfd.su2_run = os.environ['SU2_RUN'];
        
    if 'SU2_OUTPUT_FORMAT' in config:
        nozzle.cfd.su2_output_format = config['SU2_OUTPUT_FORMAT'];
    else:
        nozzle.cfd.su2_output_format = 'TECPLOT';
        
    # --- Setup outputs
	
    nozzle.responses = {}
    nozzle.gradients = {}
    
    # --- Parse fidelity levels and parameterization
    
    nozzle.SetupFidelityLevels(config, flevel, output);
    
    # --- Set flight regime + fluid
    
    nozzle.SetupMission(config,output);
    
    # --- Setup inner wall & parameterization (B-spline)
    
    nozzle.SetupInnerWall(config,output);
    
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
    
    # --- Get function responses and gradients to be returned
    
    nozzle.SetupResponsesAndGradients(config,output);

    # --- Setup DV definition
    nozzle.dvList = [];
    nozzle.outputCode = [1] * len(nozzle.outputTags); # default: output values
    nozzle.SetupDV(config,output);
  
    # --- If input DV are provided, parse them and update nozzle
    
    if nozzle.NbrDVTot > 0 :

        # Parse DV from input DV file (plain or dakota format)
        nozzle.ParseDV(config,output);

        # Update DV using values provided in input DV file
        nozzle.UpdateDV(output);

    # --- Computer inner wall's B-spline and thermal and load layer thicknesses
    #     B-spline coefs, and thickness node arrays may have been updated by
    #     the design variables input file; update exterior geometry & baffles   
    nozzle.SetupWall(output);

    return nozzle;