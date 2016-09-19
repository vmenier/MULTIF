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
        
    def SetupDV (self, config):
        nozzle = self;
        
        coefs_size = nozzle.coefs_size;

        NbrDVTot = 0;        
                
        if 'DV_LIST' in config:

            dv_keys = ('WALL', 'THERMAL_LAYER_THICKNESS_LOCATIONS',          \
                      'THERMAL_LAYER_THICKNESS_VALUES',                      \
                      'LOAD_LAYER_THICKNESS_LOCATIONS',                      \
                      'LOAD_LAYER_THICKNESS_VALUES', 'THERMAL_LAYER_DENSITY', \
                      'THERMAL_LAYER_ELASTIC_MODULUS',                       \
                      'THERMAL_LAYER_POISSON_RATIO',                         \
                      'THERMAL_LAYER_THERMAL_CONDUCTIVITY',                  \
                      'THERMAL_LAYER_THERMAL_EXPANSION_COEF',                \
                      'LOAD_LAYER_DENSITY', 'LOAD_LAYER_ELASTIC_MODULUS',    \
                      'LOAD_LAYER_POISSON_RATIO',                            \
                      'LOAD_LAYER_THERMAL_CONDUCTIVITY',                     \
                      'LOAD_LAYER_THERMAL_EXPANSION_COEF',                   \
                      'INLET_PSTAG', 'INLET_TSTAG', 'ATM_PRES', 'ATM_TEMP',  \
                      'HEAT_XFER_COEF_TO_ENV');
            # Old keys include: UPPER_WALL_THICKNESS and LOWER_WALL_THICKNESS
            # DENSITY, ELASTIC_MODULUS, POISSON_RATIO, THERMAL_CONDUCTIVITY,
            # THERMAL_EXPANSION_COEF

            hdl = config['DV_LIST'].strip('()');
            hdl = hdl.split(",");

            dv_keys_size = len(hdl)/2;
            
            nozzle.DV_Tags = []; # tags: e.g. wall, tstag etc.
            nozzle.DV_Head = []; # correspondance btw DV_Tags and DV_List
            nozzle.wall_dv = []; # 
            
            for i in range(0,2*dv_keys_size,2):
                key = hdl[i].strip();
                NbrDV = int(hdl[i+1]);
                                
                if key == 'WALL' :

                    if 'GEOM_WALL_COEFS_DV' in config:

                        wall_dv = config['GEOM_WALL_COEFS_DV'].strip('()');
                        wall_dv = wall_dv.split(",");

                        wall_dv_size = len(wall_dv);
                        
                        if  wall_dv_size != coefs_size  :
                            sys.stderr.write('  ## ERROR : Inconsistent '     \
                              'number of design variables for the inner wall:'\
                              ' %d coefs provided, %d design variables.\n',   \
                              coefs_size,  wall_dv_size);
                            sys.stderr.exit(0);
                        
                        tag = np.zeros(NbrDV);
                        for iw in range(0,wall_dv_size):

                            val = int(wall_dv[iw]);

                            if ( val > NbrDV or val < 0 ):
                                sys.stderr.write('  ## ERROR : Inconsistent ' \
                                  'design variables were provided for the '   \
                                  'inner wall (idx=%d).\n' % val);
                                sys.exit(0);

                            nozzle.wall_dv.append(val-1);

                            if val-1 >= 0 :
                                tag[val-1] = 1;
                            
                        for iw in range(0,NbrDV):
                            if tag[iw] != 1 :
                                sys.stderr.write('  ## ERROR : Design '       \
                                  'variable %d of inner wall is not defined.' \
                                  '\n' % (iw+1));
                                sys.exit(0);
                    
                    else :
                        sys.stderr.write('\n ## ERROR : Expected '           \
                                         'GEOM_WALL_COEFS_DV.\n\n');
                        sys.exit(0);
                
                elif ( key == 'THERMAL_LAYER_THICKNESS_VALUES' 
                or key == 'THERMAL_LAYER_THICKNESS_LOCATIONS'):
                    
                    if NbrDV != len(nozzle.wall.thermal_layer.thicknessNodes) :
                        sys.stderr.write('\n ## ERROR : Inconsistent number ' \
                          'of DV for thermal layer thickness definition (%d ' \
                          'DV provided instead of %d).\n\n' % (NbrDV,         \
                          len(nozzle.wall.thermal_layer.thicknessNodes)));
                        sys.exit(0);

                elif ( key == 'LOAD_LAYER_THICKNESS_VALUES' 
                or key == 'LOAD_LAYER_THICKNESS_LOCATIONS' ):
                    
                    if NbrDV != len(nozzle.wall.load_layer.thicknessNodes) :
                        sys.stderr.write('\n ## ERROR : Inconsistent number ' \
                          'of DV for load layer thickness definition (%d '    \
                          'DV provided instead of %d).\n\n' % (NbrDV,         \
                          len(nozzle.wall.load_layer.thicknessNodes)));
                        sys.exit(0);
                    
                elif ( key == 'THERMAL_LAYER_DENSITY' 
                or key == 'LOAD_LAYER_DENSITY'
                or key == 'THERMAL_LAYER_POISSON_RATIO'
                or key == 'LOAD_LAYER_POISSON_RATIO' ):
                    
                    if NbrDV != 1 :
                        sys.stderr.write('\n ## ERROR : Inconsistent number ' \
                          'of DV for density or poisson ratio (%d '           \
                          'DV provided instead of 1).\n\n' % NbrDV);
                        sys.exit(0);
    
                elif ( key == 'THERMAL_LAYER_ELASTIC_MODULUS'
                or key == 'THERMAL_LAYER_THERMAL_CONDUCTIVITY'
                or key == 'THERMAL_LAYER_THERMAL_EXPANSION_COEF' ):
                    
                    if ( NbrDV == 0 
                    or NbrDV > nozzle.wall.thermal_layer.material.n ):
                        if nozzle.wall.thermal_layer.material.n == 1:
                            strMsg = '1'
                        else:
                            strMsg = '1 or 2'
                        sys.stderr.write('\n ## ERROR : Inconsistent number ' \
                          'of DV for thermal layer material definition (%d '  \
                          'DV provided, must be %s).\n\n' % (NbrDV,strMsg));
                        sys.exit(0);  
                    
                elif ( key == 'LOAD_LAYER_ELASTIC_MODULUS'
                or key == 'LOAD_LAYER_THERMAL_CONDUCTIVITY'
                or key == 'LOAD_LAYER_THERMAL_EXPANSION_COEF' ):
                    
                    if ( NbrDV == 0 
                    or NbrDV > nozzle.wall.load_layer.material.n ):
                        if nozzle.wall.load_layer.material.n == 1:
                            strMsg = '1'
                        else:
                            strMsg = '1 or 2'
                        sys.stderr.write('\n ## ERROR : Inconsistent number ' \
                          'of DV for load layer material definition (%d '     \
                          'DV provided, must be %s).\n\n' % (NbrDV,strMsg));
                        sys.exit(0);
                
                elif ( key == 'INLET_PSTAG' or key == 'INLET_TSTAG' 
                or key == 'ATM_PRES' or key == 'ATM_TEMP'
                or key == 'HEAT_XFER_COEF_TO_ENV' ):
                    
                    if NbrDV != 1 :
                        sys.stderr.write('\n ## ERROR : Only one design '     \
                          'variable expected for %s (%d given).\n' %          \
                          (key, NbrDV));
                        sys.exit(0);
                    
                else :
                    strMsg = '';
                    for k in dv_keys :
                        strMsg = "%s %s " % (strMsg,k);
                    sys.stderr.write('  ## ERROR : Unknown design variable '  \
                      'key : %s\n' % key);
                    sys.stderr.write('             Expected = %s\n\n' % strMsg);
                    sys.exit(0);

                nozzle.DV_Tags.append(key);
                nozzle.DV_Head.append(NbrDVTot);   
                 
                NbrDVTot = NbrDVTot + NbrDV;
                
            # for i in keys
        
        else :
            sys.stdout.write('\n  -- Info : No design variable set was '      \
              'defined. Running the baseline parameters.\n\n');
        
        nozzle.NbrDVTot = NbrDVTot;
        
        if NbrDVTot > 0 :
            nozzle.DV_Head.append(NbrDVTot);
            

    def SetupFidelityLevels (self, config, flevel, output='verbose'):
        
        nozzle = self;
        
        fidelity_tags = config['FIDELITY_LEVELS_TAGS'].strip('()');
        fidelity_tags = fidelity_tags.split(",");

        NbrFidLev = len(fidelity_tags);

        if NbrFidLev < 1 :
            print "  ## ERROR : No fidelity level was defined.\n"
            sys.exit(0)

        if output == 'verbose':
          sys.stdout.write('\n%d fidelity level(s) defined. Summary :\n'      \
            % (NbrFidLev));

          sys.stdout.write('-' * 90);
          sys.stdout.write('\n%s | %s | %s\n' % ("Level #".ljust(10),         \
            "Tag".ljust(10),"Description.".ljust(70)));
          sys.stdout.write('-' * 90);
          sys.stdout.write('\n'); 
        elif output == 'quiet':
          pass
        else:
          raise ValueError('keyword argument output can only be set to '      \
            '"verbose" or "quiet" mode')

        for i in range(NbrFidLev) :
            tag = fidelity_tags[i];
            kwd = "DEF_%s" % tag;

            if kwd not in config :
                sys.stderr.write("\n  ## ERROR : The fidelity level tagged '  \
                  '%s is not defined.\n\n" % kwd);
                sys.exit(0);

            cfgLvl = config[kwd].strip('()');
            cfgLvl = cfgLvl.split(",");

            method = cfgLvl[0];

            description = "";

            if i == flevel :
                nozzle.method = method;

            if method == 'NONIDEALNOZZLE' :

                tol = float(cfgLvl[1]);
                if tol < 1e-30 :
                    sys.stderr.write("\n ## ERROR : Wrong tolerance for "     \
                     "fidelity level %d (tagged %s)\n\n" % (i,tag));
                    sys.exit(0);

                if i == flevel :
                    nozzle.tolerance = tolerance.Tolerance();
                    nozzle.tolerance.setRelTol(tol);
                    nozzle.tolerance.setAbsTol(tol);
                    #nozzle.tolerance.exitTempPercentError = tol;
                description = "ODE solver relative and absolute tolerance "   \
                  "set to %le." % (tol);    

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
                description = "%s %s CFD simulation using the %s mesh level." \
                  % (dim, method, meshsize);

                if i == flevel :
                    nozzle.meshsize = meshsize;
                    
                    nozzle.bl_ds        = 0.000007;
                    nozzle.bl_ratio     = 1.3; 
                    nozzle.bl_thickness = 0.1;

                    scaleMesh = 1.0;
                    if meshsize == 'COARSE':
                        scaleMesh = 2.5;
                    elif meshsize == 'MEDIUM':
                        scaleMesh = 1.1;
                    elif meshsize == 'FINE':
                        scaleMesh = 0.5;

                    nozzle.meshhl = scaleMesh*np.asarray([0.1, 0.07, 0.06, 0.006, 0.0108]);

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
            sys.stdout.write('  -- Info : Fidelity level to be run : "        \
              "%d\n\n' % flevel);
        elif output == 'quiet':
            pass
        else:
            raise ValueError('keyword argument output can only be set to '    \
              '"verbose" or "quiet" mode')
        
    def SetupMission(self, config):
    
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
          
            mission_id = 0; # to denote custom mission
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
            hInf = 7.2; # W/m^2/K
        nozzle.environment = environment.Environment(altitude,hInf);

        # --- Setup fluid

        heatRatio = 1.4;
        gasCst    = 287.06;
        nozzle.fluid = fluid.Fluid(heatRatio, gasCst);
        
        # --- Setup convergence parameter
        
        if 'SU2_CONVERGENCE_ORDER' in config:
            nozzle.su2_convergence_order = config['SU2_CONVERGENCE_ORDER'];
        else:
            nozzle.su2_convergence_order = 3;
        
    def SetupBSplineCoefs(self, config):
        
        nozzle = self;
        
        # Set up nozzle wall
        nozzle.wall = component.AxisymmetricWall();
        
        wall_keys = ('GEOM_WALL_PARAM','GEOM_WALL_KNOTS','GEOM_WALL_COEFS');

        if all (key in config for key in wall_keys):

            # --- Get coefs

            hdl = config['GEOM_WALL_COEFS'].strip('()');
            hdl = hdl.split(",");

            coefs_size = len(hdl);

            coefs = [];

            for i in range(0,coefs_size):
                coefs.append(float(hdl[i]));

            # --- Get knots

            hdl = config['GEOM_WALL_KNOTS'].strip('()');
            hdl = hdl.split(",");

            knots_size = len(hdl);

            knots = [];

            for i in range(0,knots_size):
                knots.append(float(hdl[i]));

            #print "%d coefs, %d knots\n" % (coefs_size, knots_size);

        else:
            str = '';
            for key in wall_keys :
                str = "%s %s " % (str,key);
            sys.stderr.write('\n ## ERROR : NO INNER WALL DEFINITION IS '    \
                             'PROVIDED.\n\n');
            sys.stderr.write('             Expected = %s\n\n' % str);
            sys.exit(0);
        
        nozzle.coefs =  coefs;
        nozzle.knots =  knots;
        nozzle.coefs_size = coefs_size;
        
    
    def ParseThickness(self, config, key):
        
        loc_name = '%s_THICKNESS_LOCATIONS' % key;
        val_name = '%s_THICKNESS_VALUES' % key;
        
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
                         '_THERMAL_CONDUCTIVITY','_THERMAL_EXPANSION_COEF'];  
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
            if len(var) == 2:
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
            else:
                sys.stderr.write('\n ## ERROR: Key %s could not implemented' \
                ' for %s layer' % (key,name));
                sys.exit(0);        
        
    
    def SetupWallLayers(self, config):
  
        nozzle = self;
        
        # Inner thermal layer to take heat load
        nozzle.wall.thermal_layer = component.AxisymmetricWall('THERMAL_LAYER');
        
        # Outer structural layer to take structural load
        nozzle.wall.load_layer = component.AxisymmetricWall('LOAD_LAYER');
        
        # Setup thickness and material for thermal layer
        nozzle.SetupLayer(config,nozzle.wall.thermal_layer,'single')
        
        # Setup thickness and material for load layer
        nozzle.SetupLayer(config,nozzle.wall.load_layer,'single')
        

    def SetupWall (self, config):
        
        nozzle = self;
        
        # --- SHAPE OF INNER WALL (B-SPLINE)
        
        coefs = nozzle.coefs;
        knots = nozzle.knots;
        
        coefs_size = nozzle.coefs_size;
        
        nozzle.height = coefs[coefs_size-1];
        nozzle.length = coefs[coefs_size/2-1];
        
        if nozzle.method == 'RANS' or nozzle.method == 'EULER':
            x    = [];
            y    = [];

            nx = 100;
            _meshutils_module.py_BSplineGeo3 (knots, coefs, x, y, nx);

            nozzle.xwall = x;
            nozzle.ywall = y;
    
        coefsnp = np.empty([2, coefs_size/2]);
    
        for i in range(0, coefs_size/2):
            coefsnp[0][i] = coefs[i];
            coefsnp[1][i] = coefs[i+coefs_size/2];        
        
        nozzle.wall.geometry = geometry.Bspline(coefsnp);
        
        # --- THERMAL LAYER THICKNESS
        
        size_thickness = len(nozzle.wall.thermal_layer.thicknessNodes);
        thicknessNodeArray = np.zeros(shape=(2,size_thickness))
        
        for i in range(size_thickness):
            thicknessNodeArray[0][i] = nozzle.wall.thermal_layer.thicknessNodes[i][0]*nozzle.length;
            thicknessNodeArray[1][i] = nozzle.wall.thermal_layer.thicknessNodes[i][1];
        
        nozzle.wall.thermal_layer.thickness = geometry.PiecewiseLinear(thicknessNodeArray);

        # --- LOAD LAYER THICKNESS
        
        size_thickness = len(nozzle.wall.load_layer.thicknessNodes);
        thicknessNodeArray = np.zeros(shape=(2,size_thickness))
        
        for i in range(size_thickness):
            thicknessNodeArray[0][i] = nozzle.wall.load_layer.thicknessNodes[i][0]*nozzle.length;
            thicknessNodeArray[1][i] = nozzle.wall.load_layer.thicknessNodes[i][1];
        
        nozzle.wall.load_layer.thickness = geometry.PiecewiseLinear(thicknessNodeArray);
            

    def ParseDV (self, config):
        
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
                
                for iCoef in range(len(nozzle.wall_dv)):
                    id_dv = nozzle.DV_Head[iTag] + nozzle.wall_dv[iCoef];
                    
                    # --- Update coef iCoef if required
                    if id_dv >= 0 :
                        prt_name.append('Bspline coef #%d' % (iCoef+1));
                        prt_basval.append('%.4lf'% nozzle.coefs[iCoef]);
                        prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                        nozzle.coefs[iCoef] = nozzle.DV_List[id_dv];
                        NbrChanged = NbrChanged+1;
                        
            elif Tag == 'THERMAL_LAYER_THICKNESS_LOCATIONS':
                
                for i in range(NbrDV):
                    id_dv = nozzle.DV_Head[iTag]+i;
                    
                    prt_name.append('Thermal layer thickness location #%i' % i);
                    prt_basval.append('%.4lf'%                               \
                      nozzle.wall.thermal_layer.thicknessNodes[i][0]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    
                    nozzle.wall.thermal_layer.thicknessNodes[i][0] =              \
                      nozzle.DV_List[id_dv];
                        
            elif Tag == 'THERMAL_LAYER_THICKNESS_VALUES':
                
                for i in range(NbrDV):
                    id_dv = nozzle.DV_Head[iTag]+i;
                    
                    # if THERMAL_LAYER_THICKNESS_LOCATIONS is not updated yet,
                    # then the following print label will not be correct
                    prt_name.append('Thermal layer thickness t=%.3lf'        \
                      % nozzle.wall.thermal_layer.thicknessNodes[i][0]);
                    prt_basval.append('%.4lf'%                               \
                      nozzle.wall.thermal_layer.thicknessNodes[i][1]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    
                    nozzle.wall.thermal_layer.thicknessNodes[i][1] =         \
                      nozzle.DV_List[id_dv];
            
            elif Tag == 'LOAD_LAYER_THICKNESS_LOCATIONS':

                for i in range(NbrDV):
                    id_dv = nozzle.DV_Head[iTag]+i;
                    
                    prt_name.append('Load layer thickness location #%i' % i);
                    prt_basval.append('%.4lf'%                               \
                      nozzle.wall.load_layer.thicknessNodes[i][0]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    
                    nozzle.wall.load_layer.thicknessNodes[i][0] =            \
                      nozzle.DV_List[id_dv];
            
            elif Tag == 'LOAD_LAYER_THICKNESS_VALUES':
      
                for i in range(NbrDV):
                    id_dv = nozzle.DV_Head[iTag]+i;
                    
                    # if LOAD_LAYER_THICKNESS_LOCATIONS is not updated yet,
                    # then the following print label will not be correct
                    prt_name.append('Load layer thickness t=%.3lf'           \
                      % nozzle.wall.load_layer.thicknessNodes[i][0]);
                    prt_basval.append('%.4lf'%                               \
                      nozzle.wall.load_layer.thicknessNodes[i][1]);
                    prt_newval.append('%.4lf'% nozzle.DV_List[id_dv]);
                    
                    nozzle.wall.load_layer.thicknessNodes[i][1] =            \
                      nozzle.DV_List[id_dv];   
            
            elif Tag == 'THERMAL_LAYER_DENSITY':
                
                # Only single material enabled; panel material not implemented
                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Thermal layer density');
                prt_basval.append('%.2le'%                                   \
                  nozzle.wall.thermal_layer.material.getDensity());
                prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                
                var = nozzle.DV_List[id_dv];
                nozzle.wall.thermal_layer.material.setDensity(var);
            
            elif Tag == 'THERMAL_LAYER_ELASTIC_MODULUS':
                
                id_dv = nozzle.DV_Head[iTag];
                
                if NbrDV == 1:
                    prt_name.append('Thermal layer elastic modulus');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.E);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    var = nozzle.DV_List[id_dv]; 
                    nozzle.wall.thermal_layer.material.setElasticModulus(var);
                else: # NbrDV == 2
                    prt_name.append('Thermal layer elastic modulus (axial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.EAxial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    prt_name.append('Thermal layer elastic modulus (radial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.ERadial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                    
                    var = nozzle.DV_List[id_dv:id_dv+2]; 
                    nozzle.wall.thermal_layer.material.setElasticModulus(var);                
                
            elif Tag == 'THERMAL_LAYER_POISSON_RATIO':

                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Thermal layer Poisson ratio');
                prt_basval.append('%.4f'%                                   \
                  nozzle.wall.thermal_layer.material.getPoissonRatio());
                prt_newval.append('%.4f'% nozzle.DV_List[id_dv]);
                
                var = nozzle.DV_List[id_dv];
                nozzle.wall.thermal_layer.material.setPoissonRatio(var);
            
            elif Tag == 'THERMAL_LAYER_THERMAL_CONDUCTIVITY':
                
                id_dv = nozzle.DV_Head[iTag];
                
                if NbrDV == 1:
                    prt_name.append('Thermal layer thermal cond.');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.k);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    var = nozzle.DV_List[id_dv]; 
                    nozzle.wall.thermal_layer.material.setThermalConductivity(var);
                else: # NbrDV == 2
                    prt_name.append('Thermal layer thermal cond. (axial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.kAxial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    prt_name.append('Thermal layer thermal cond. (radial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.kRadial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                    
                    var = nozzle.DV_List[id_dv:id_dv+2]; 
                    nozzle.wall.thermal_layer.material.setThermalConductivity(var); 
                
            elif Tag == 'THERMAL_LAYER_THERMAL_EXPANSION_COEF':
                
                id_dv = nozzle.DV_Head[iTag];
                
                if NbrDV == 1:
                    prt_name.append('Thermal layer thermal expansion coef.');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.alpha);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    var = nozzle.DV_List[id_dv]; 
                    nozzle.wall.thermal_layer.material.setThermalExpansionCoef(var);
                else: # NbrDV == 2
                    prt_name.append('Thermal layer thermal expansion coef. (axial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.alphaAxial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    prt_name.append('Thermal layer thermal expansion coef. (radial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.thermal_layer.material.alphaRadial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                    
                    var = nozzle.DV_List[id_dv:id_dv+2]; 
                    nozzle.wall.thermal_layer.material.setThermalExpansionCoef(var);                 
            
            elif Tag == 'LOAD_LAYER_DENSITY':

                # Only single material enabled; panel material not implemented
                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Load layer density');
                prt_basval.append('%.2le'%                                   \
                  nozzle.wall.load_layer.material.getDensity());
                prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                
                var = nozzle.DV_List[id_dv];
                nozzle.wall.load_layer.material.setDensity(var);
            
            elif Tag == 'LOAD_LAYER_ELASTIC_MODULUS':

                id_dv = nozzle.DV_Head[iTag];
                
                if NbrDV == 1:
                    prt_name.append('Load layer elastic modulus');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.E);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    var = nozzle.DV_List[id_dv]; 
                    nozzle.wall.load_layer.material.setElasticModulus(var);
                else: # NbrDV == 2
                    prt_name.append('Load layer elastic modulus (axial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.EAxial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    prt_name.append('Load layer elastic modulus (radial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.ERadial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                    
                    var = nozzle.DV_List[id_dv:id_dv+2]; 
                    nozzle.wall.load_layer.material.setElasticModulus(var);
                
            elif Tag == 'LOAD_LAYER_POISSON_RATIO':

                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Load layer Poisson ratio');
                prt_basval.append('%.4f'%                                   \
                  nozzle.wall.load_layer.material.getPoissonRatio());
                prt_newval.append('%.4f'% nozzle.DV_List[id_dv]);
                
                var = nozzle.DV_List[id_dv];
                nozzle.wall.load_layer.material.setPoissonRatio(var);
            
            elif Tag == 'LOAD_LAYER_THERMAL_CONDUCTIVITY':

                id_dv = nozzle.DV_Head[iTag];
                
                if NbrDV == 1:
                    prt_name.append('Load layer thermal cond.');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.k);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    var = nozzle.DV_List[id_dv]; 
                    nozzle.wall.load_layer.material.setThermalConductivity(var);
                else: # NbrDV == 2
                    prt_name.append('Load layer thermal cond. (axial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.kAxial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    prt_name.append('Load layer thermal cond. (radial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.kRadial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                    
                    var = nozzle.DV_List[id_dv:id_dv+2]; 
                    nozzle.wall.load_layer.material.setThermalConductivity(var); 
                
            elif Tag == 'LOAD_LAYER_THERMAL_EXPANSION_COEF':

                id_dv = nozzle.DV_Head[iTag];
                
                if NbrDV == 1:
                    prt_name.append('Load layer thermal expansion coef.');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.alpha);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    var = nozzle.DV_List[id_dv]; 
                    nozzle.wall.load_layer.material.setThermalExpansionCoef(var);
                else: # NbrDV == 2
                    prt_name.append('Load layer thermal expansion coef. (axial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.alphaAxial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv]);
                    
                    prt_name.append('Load layer thermal expansion coef. (radial)');
                    prt_basval.append('%.2le'%                               \
                      nozzle.wall.load_layer.material.alphaRadial);
                    prt_newval.append('%.2le'% nozzle.DV_List[id_dv+1]);
                    
                    var = nozzle.DV_List[id_dv:id_dv+2]; 
                    nozzle.wall.load_layer.material.setThermalExpansionCoef(var);  

            elif Tag == 'INLET_PSTAG':
                
                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Inlet stagnation pressure');
                prt_basval.append('%.2lf'% nozzle.inlet.Pstag);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
                
                nozzle.inlet.setPstag(nozzle.DV_List[id_dv]);
                
            elif Tag == 'INLET_TSTAG':
                
                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Inlet stagnation temperature');
                prt_basval.append('%.2lf'% nozzle.inlet.Tstag);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
                
                nozzle.inlet.setTstag(nozzle.DV_List[id_dv]);

            elif Tag == 'ATM_PRES':
                
                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Atmospheric pressure');
                prt_basval.append('%.2lf'% nozzle.environment.P);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
                
                nozzle.environment.setPressure(nozzle.DV_List[id_dv]);
                
            elif Tag == 'ATM_TEMP':
                
                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Atmospheric temperature');
                prt_basval.append('%.2lf'% nozzle.environment.T);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
                
                nozzle.environment.setTemperature(nozzle.DV_List[id_dv]);

            elif Tag == 'HEAT_XFER_COEF_TO_ENV':

                id_dv = nozzle.DV_Head[iTag];
                
                prt_name.append('Heat xfer coef. to environment');
                prt_basval.append('%.2lf'% nozzle.environment.hInf);
                prt_newval.append('%.2lf'% nozzle.DV_List[id_dv]);
                
                nozzle.environment.setHeatTransferCoefficient(nozzle.DV_List[id_dv]);

            else :
                sys.stderr.write("\n ## Error : Unknown tag name %s\n\n" % Tag);
                sys.exit(0);
                
            if Tag != 'WALL':
                NbrChanged = NbrChanged + NbrDV;
        
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
        
    def SetupOutputFunctions (self, config):
        
        nozzle = self;
        
        nozzle.Output_Tags = [];
        
        nozzle.GetOutput = dict();
        
        if 'OUTPUT_NAME' in config:
            nozzle.Output_Name = config['OUTPUT_NAME'];
        else :
            sys.stderr.write("\n ## ERROR : Output function file name not "   \
              "specified in config file. (OUTPUT_NAME expected)\n\n");
            sys.exit(0);
        
        # --- Initialize outputs
        nozzle.thrust = -1;
        nozzle.volume = -1;
        nozzle.max_thermal_stress    = -1;
        nozzle.max_mechanical_stress = -1;
        
        if 'OUTPUT_FUNCTIONS' in config:

            dv_keys = ('VOLUME', 'THRUST');

            hdl = config['OUTPUT_FUNCTIONS'].strip('()');
            hdl = hdl.split(",");
            
            for i in range(len(hdl)):
                
                key = hdl[i].strip();
                        
                if( key == 'VOLUME' or key == 'THRUST'
                or 'MAX_MECHANICAL_STRESS' or 'MAX_THERMAL_STRESS' ):
                    nozzle.GetOutput[key] = 1;
                    nozzle.Output_Tags.append(key);
                
                else :
                    str = '';
                    for k in dv_keys :
                        str = "%s %s " % (str,k);
                    sys.stderr.write('\n ## ERROR : Unknown output '          \
                      'function name : %s\n' % key);
                    sys.stderr.write('             Expected = %s\n\n' % str);
                    sys.exit(0);
                
            # for i in keys
        
        if len(nozzle.Output_Tags) == 0 :
            sys.stderr.write("\n  ## Error : No output function was given.\n\n");
            sys.exit(0);


    def WriteOutputFunctions_Plain (self, output='verbose'):
      
        nozzle = self;

        filename = nozzle.Output_Name;
        
        if output == 'verbose':
            sys.stdout.write('\n');
            str = " Post-processing ";
            nch = (60-len(str))/2;
            sys.stdout.write('-' * nch);
            sys.stdout.write(str);
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

        for i in range(0, len(nozzle.Output_Tags)):
          
            tag = nozzle.Output_Tags[i];
          
            if tag == 'THRUST':
                fil.write('%0.16f\n' % nozzle.thrust);
                if output == 'verbose':
                    sys.stdout.write('      Thrust = %0.16f\n' % nozzle.thrust);
  
            if tag == 'VOLUME':
                fil.write('%0.16f\n' % nozzle.volume);
                if output == 'verbose':
                    sys.stdout.write('      Volume = %0.16f\n' % nozzle.volume);
                    
            if tag == 'MAX_MECHANICAL_STRESS':
                fil.write('%0.16f max_mechanical_stress\n' % nozzle.max_mechanical_stress);
                if output == 'verbose':
                    sys.stdout.write('      Max mechanical stress = %0.16f\n' % nozzle.max_mechanical_stress);
                    
            if tag == 'MAX_THERMAL_STRESS':
                fil.write('%0.16f max_thermal_stress\n' % nozzle.max_thermal_stress);
                if output == 'verbose':
                    sys.stdout.write('      Max thermal stress = %0.16f\n' % nozzle.max_thermal_stress);                    
      
        if output == 'verbose':
            sys.stdout.write('\n');
        fil.close();
        
    def WriteOutputFunctions_Dakota (self,output='verbose'):
    
        nozzle = self;
    
        filename = nozzle.Output_Name;
    
        try:
            fil = open(filename, 'w');
        except:
            sys.stderr.write("  ## ERROR : Could not open output file %s\n" % filename);
            sys.exit(0);
    
        sys.stdout.write('  -- Info : Output functions file : %s\n' % filename);
    
        for i in range(0, len(nozzle.Output_Tags)):
    
            tag = nozzle.Output_Tags[i];
    
            if tag == 'THRUST':
                fil.write('%0.16f thrust\n' % nozzle.thrust);
                if output == 'verbose':
                    sys.stdout.write('      Thrust = %0.16f\n' % nozzle.thrust);
                    
            if tag == 'VOLUME':
                fil.write('%0.16f volume\n' % nozzle.volume);
                if output == 'verbose':
                    sys.stdout.write('      Volume = %0.16f\n' % nozzle.volume);
                    
            if tag == 'MAX_MECHANICAL_STRESS':
                fil.write('%0.16f max_mechanical_stress\n' % nozzle.max_mechanical_stress);
                if output == 'verbose':
                    sys.stdout.write('      Max mechanical stress = %0.16f\n' % nozzle.max_mechanical_stress);
            
            if tag == 'MAX_THERMAL_STRESS':
                fil.write('%0.16f max_thermal_stress\n' % nozzle.max_thermal_stress);
                if output == 'verbose':
                    sys.stdout.write('      Max thermal stress = %0.16f\n' % nozzle.max_thermal_stress); 
                    
        sys.stdout.write('\n');
        fil.close();
        
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

def NozzleSetup( config_name, flevel, output='verbose' ):
    import tempfile
        
    if not os.path.isfile(config_name) :
        sys.stderr.write("  ## ERROR : could not find configuration file %s\n\ns" % config_name);
        sys.exit(0);
    
    nozzle = multif.nozzle.nozzle.Nozzle()
    
    config = SU2.io.Config(config_name)
    
    print config
    
    # --- File names
    
     #nozzle.mesh_name    =     'blabla.su2'
     #nozzle.restart_name =  'blabla.dat'

    #hdl, toto = tempfile.mkstemp(suffix='.su2');
     #
    #print "TMPNAME = %s" %  (toto);

    nozzle.mesh_name    =  'nozzle.su2'; #tempfile.mkstemp(suffix='.su2');
    nozzle.restart_name =  'nozzle.dat'; #tempfile.mkstemp(suffix='.dat');
    
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
    
    # --- Parse fidelity levels
    
    nozzle.SetupFidelityLevels(config, flevel, output);
    
    # --- Set flight regime + fluid
    
    nozzle.SetupMission(config);
    
    # --- Setup inner wall & parameterization (B-spline)
    
    nozzle.SetupBSplineCoefs(config);
    
    # --- Setup wall layer thickness(es) and material(s)
    
    nozzle.SetupWallLayers(config);
    
    # --- Setup DV definition
        
    nozzle.SetupDV(config);
    
    # --- If input DV are provided, parse them and update nozzle
    
    if nozzle.NbrDVTot > 0 :    
        
        # Parse DV from input DV file (plain or dakota format)
        nozzle.ParseDV(config);
        
        # Update DV using values provided in input DV file
        nozzle.UpdateDV(config,output);
    
    # --- Computer inner wall's B-spline and thermal and load layer thicknesses
    #     B-spline coefs, and thickness node arrays may have been updated by
    #     the design variables input file
    
    nozzle.SetupWall(config);
    
    # --- Get output functions to be returned
    
    nozzle.SetupOutputFunctions(config);
        
    #sys.exit(1);
    return nozzle;    
