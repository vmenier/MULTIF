import sys, os

import numpy as np

from .. import SU2
from hf_meshgeneration import *
from hf_runSU2 import *
import hf_postprocessing

try:
    from runAEROS import *
except ImportError:
    print 'Error importing all functions from runAEROS.\n'

def CheckOptions (nozzle):
    
    print "Check options"
    if nozzle.dim != '3D':
        sys.stderr.write("\n  ## ERROR : High-fidelity is 3D.\n\n");
        sys.exit(0);
    
def Run( nozzle, output = 'verbose', writeToFile=1 ):
	
	
	
	#HF_GenerateMesh_Deform(nozzle);
	#
	#return;
	
	# Introduction to this function and available data
	print
	print '*****'
	print
	print 'MISSION:'
	print nozzle.mission.__dict__
	print
	
	print 'CFD-related information:'
	print nozzle.cfd.__dict__
	print
	
	print 'Fluid-related information:'
	print nozzle.fluid.__dict__
	print
	
	print 'Exterior environment information:'
	print nozzle.environment.__dict__
	print
	
	print 'INLET information:'
	print nozzle.inlet.__dict__
	print
	
	print 'NOZZLE INTERIOR WALL:'
	print 'Nozzle wall centerline geometry:'
	print nozzle.wall.centerline.geometry.__dict__
	print 'Nozzle wall major axis geometry:'
	print nozzle.wall.majoraxis.geometry.__dict__
	print 'Nozzle wall minor axis geometry:'
	print nozzle.wall.minoraxis.geometry.__dict__
	print 'Nozzle length: %f' % nozzle.length
	print 'Nozzle wall struct:'
	print nozzle.wall.__dict__
	print 'Nozzle wall shovel exit height: %f' % nozzle.wall.shovel_height
	print 'Nozzle wall shovel inlet start angle: %f' % nozzle.wall.shovel_start_angle
	print
	
	print 'WALL LAYERS:'
	for i in range(len(nozzle.wall.layer)):
	    print 'Nozzle wall layer: %s' % nozzle.wall.layer[i].name
	    print nozzle.wall.layer[i].__dict__
	print
	
	print 'STRINGERS:'
	print nozzle.stringers.__dict__
	print 'Example access of stringers height definition for 2nd stringer:'
	print nozzle.stringers.height[1].__dict__
	print 'Example access of stringers thickness definition for 2nd stringer:'
	print nozzle.stringers.thickness[1].__dict__
	print
	
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
	    
	# Obtain mass and volume
	if 'MASS' in nozzle.responses or 'VOLUME' in nozzle.responses:
	    # XXX Calculate mass and volume here from CAD geometry & material prop.
	    # Perhaps modify existing function:
	    # volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle)
	    volume = [0,0];
	    mass = [0,0];
	    if 'MASS' in nozzle.responses:
	        nozzle.responses['MASS'] = np.sum(mass)
	    if 'VOLUME' in nozzle.responses:
	        nozzle.responses['VOLUME'] = np.sum(volume)
	
	# Calculate mass gradients if necessary
	if 'MASS' in nozzle.gradients and nozzle.gradients['MASS'] is not None:
	    if ( nozzle.gradientsMethod == 'ADJOINT' ):
	        # XXX Calculate derivative of mass here, using FD when adjoint is 
	        # used for fluid. Perhaps modify existing function:
	        # nozzle.gradients['MASS'] = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-8)
	        nozzle.gradients['MASS'] = [0.]*len(nozzle.derivativesDV);
	    elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
	        # XXX Calculate derivative of mass here using FD.
	        # Perhaps modify existing function:
	        # nozzle.gradients['MASS'] = \
	        #   nozzlemod.geometry.calcMassGradientsFD(nozzle,nozzle.fd_step_size)  
	        nozzle.gradients['MASS'] = [0.]*len(nozzle.derivativesDV);          
	    else:
		    sys.stderr.write('  ## ERROR : Unknown gradients computation '
		      'method.\n');
		    sys.exit(1);
	
	# Calculate volume gradients if necessary
	if 'VOLUME' in nozzle.gradients and nozzle.gradients['VOLUME'] is not None:
	    sys.stderr.write('\n ## ERROR : gradients for VOLUME are not supported\n\n');
	    sys.exit(1);
	    
	# Obtain mass of wall and gradients only if requested
	if 'MASS_WALL_ONLY' in nozzle.responses:
	    # XXX Calculate mass of nozzle wall only (no stringers, no baffles)
	    # from CAD geometry & material properties
	    # This QoI is not very important, so feel free to skip this.
	    nozzle.responses['MASS_WALL_ONLY'] = 0.
	    
	if 'MASS_WALL_ONLY' in nozzle.gradients and nozzle.gradients['MASS_WALL_ONLY'] is not None:
	    sys.stderr.write('\n ## ERROR : gradients for MASS_WALL_ONLY are not supported\n\n');
	    sys.exit(1);
	
	# Run aero-thermal-structural analysis if necessary
	runAeroThermalStructuralProblem = 0;
	for k in nozzle.responses:
	    if k not in ['MASS','VOLUME','MASS_WALL_ONLY']:
	        runAeroThermalStructuralProblem = 1;    
	
	# Run aero-thermal-structural gradient analysis if necessary
	if nozzle.gradientsFlag == 1:
	    runAeroThermalStructuralGradients = 0;        
	    for k in nozzle.gradients:
	        if k not in ['MASS','VOLUME','MASS_WALL_ONLY']:
	            if nozzle.gradients[k] is not None:
	                runAeroThermalStructuralGradients = 1;
	
	if runAeroThermalStructuralProblem:
	    
	    # Run aero analysis (and thrust adjoint if necessary)  
	    CheckSU2Version(nozzle);	
	    CheckOptions(nozzle);
	    
	    curDir = os.path.dirname(os.path.realpath(__file__));	
	    if nozzle.runDir != '':
	        os.chdir(nozzle.runDir);
	    # XXX Ensure high-fidelity SU2 analysis runs correctly
	    gradCalc = HF_runSU2(nozzle);
		
	   # Run thermal/structural analyses
	    if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
	        # XXX Ensure high-fidelity AEROS analysis runs correctly
	        try : 
	            runAEROS(nozzle, output);
	        except:
	           sys.stdout.write("  ## WARNING: CALL TO AERO-S IGNORED.\n");
	           
		
	    # Assign aero QoI if required
	    # XXX Fill in this PostProcess function to return correct outputs
	    hf_postprocessing.PostProcess(nozzle, output);
	    
	    # Assign thermal/structural QoI if required
	    if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
	        # XXX Ensure AEROS post-processing returns correct outputs
	        #AEROSpostprocessing.PostProcess(nozzle, output);      
	        try : 
	           AEROSpostprocessing.PostProcess(nozzle, output);    
	        except:
	           sys.stdout.write("  ## WARNING: CALL TO AERO-S IGNORED.\n");       
	    
	    # Calculate gradients if necessary
	    if nozzle.gradientsFlag == 1 and runAeroThermalStructuralGradients:
	
	        if ( nozzle.gradientsMethod == 'ADJOINT' ):
				
	            if gradCalc == 0: # i.e. failed adjoint calculation, use finite differences
	                # XXX Ensure FD gradients work for 3-D nozzle. If this
	                # function needs to be changed, ensure it is also changed
	                # in the other 2 spots below too.
	                multif.gradients.calcGradientsFD(nozzle,nozzle.fd_step_size,output);
	            else:
	                # Check for other required gradients
	                otherRequiredGradients = 0;
	                for k in nozzle.gradients:
	                    if k not in ['MASS','VOLUME','MASS_WALL_ONLY','THRUST']:
	                        otherRequiredGradients = 1;
	                        sys.stderr.write(' ## WARNING: QoI gradients desired using ADJOINT '
	                          'method which do not have an associated adjoint calculation.\n'
	                          ' Namely: %s. The current implementation requires finite '
	                          'differencing across the aero analysis, so this method is '
	                          'equivalent in computational cost to choosing the FINITE_DIFF'
	                          ' method\n' % k);
	                # Do finite difference for other QoI if there are any, but use
	                # adjoint gradients for thrust
	                if otherRequiredGradients:
	                    saveThrustGradients = nozzle.gradients['THRUST'];
	                    nozzle.gradients['THRUST'] = None;
	                    multif.gradients.calcGradientsFD(nozzle,nozzle.fd_step_size,output);
	                    nozzle.gradients['THRUST'] = saveThrustGradients;
	                    
	        elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
	        
	            multif.gradients.calcGradientsFD(nozzle,nozzle.fd_step_size,output);  
	                     
	        else:
			    sys.stderr.write('  ## ERROR : Unknown gradients computation '
			      'method.\n');
			    sys.exit(1);
			    
		 # Write separate gradients file
	        gradFile = open(nozzle.gradientsFile,'w');
	        for k in nozzle.outputTags:
			    np.savetxt(gradFile,nozzle.gradients[k]);
	        gradFile.close();   			
	
	# Write data
	if writeToFile:
	    if nozzle.outputFormat == 'PLAIN':
	        nozzle.WriteOutputFunctions_Plain();
	    else:
	        nozzle.WriteOutputFunctions_Dakota();
	
	return 0;
    
    
