import sys, os

import numpy as np

from .. import SU2
from hf_meshgeneration import *
from hf_runSU2 import *
import hf_postprocessing

from .. import nozzle as nozzlemod

try:
    from multif.MEDIUMF.runAEROS import *
except ImportError as e:
    print 'Error importing all functions from runAEROS in hf_run.py.'
    print e
    print


def CheckOptions (nozzle):
    
    print "Check options"
    if nozzle.dim != '3D':
        sys.stderr.write("\n  ## ERROR : High-fidelity is 3D.\n\n");
        sys.exit(0);
    
def Run( nozzle, output = 'verbose', writeToFile=1 ):

    # Obtain mass and volume
	if 'MASS' in nozzle.responses or 'VOLUME' in nozzle.responses or 'MASS_WALL_ONLY' in nozzle.responses:
		volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle)
        if 'MASS' in nozzle.responses:
            nozzle.responses['MASS'] = np.sum(mass)
        if 'VOLUME' in nozzle.responses:
            nozzle.responses['VOLUME'] = np.sum(volume)
        if 'MASS_WALL_ONLY' in nozzle.responses:
		    nozzle.responses['MASS_WALL_ONLY'] = np.sum(mass[:len(nozzle.wall.layer)])

    # Calculate mass gradients if necessary
	if 'MASS' in nozzle.gradients and nozzle.gradients['MASS'] is not None:
		if ( nozzle.gradientsMethod == 'ADJOINT' ):
			# Convergence study using B-spline coefs show finite difference mass gradients
			# converge. Most accurate gradients use absolute step size 1e-8. RWF 5/10/17
			# The above conclusion is based on 2D axisymmetric mass calculations.
			nozzle.gradients['MASS'] = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-8);
		elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
			nozzle.gradients['MASS'] = nozzlemod.geometry.calcMassGradientsFD( \
				nozzle,nozzle.fd_step_size);
		else:
			sys.stderr.write('  ## ERROR : Unknown gradients computation '
				'method.\n');
			sys.exit(1);
	
    # Calculate volume gradients if necessary
	if 'VOLUME' in nozzle.gradients and nozzle.gradients['VOLUME'] is not None:
		sys.stderr.write('\n ## ERROR : gradients for VOLUME are not supported\n\n');
		sys.exit(1);
        
    # Obtain gradients of mass of wall if requested        
	if 'MASS_WALL_ONLY' in nozzle.gradients and nozzle.gradients['MASS_WALL_ONLY'] is not None:
		if ( nozzle.gradientsMethod == 'ADJOINT' ):
			# Convergence study using B-spline coefs show finite difference mass gradients
			# converge. Most accurate gradients use absolute step size 1e-8. RWF 5/10/17
			# The above conclusion is based on 2D axisymmetric mass calculations.
			nozzle.gradients['MASS_WALL_ONLY'] = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-8,components='wall-only');
		elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
			nozzle.gradients['MASS_WALL_ONLY'] = nozzlemod.geometry.calcMassGradientsFD(\
				nozzle,nozzle.fd_step_size,components='wall-only');
		else:
			sys.stderr.write('  ## ERROR : Unknown gradients computation '
				'method.\n');
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
			runAEROS(nozzle, output);
			
	    # Assign aero QoI if required
	    # XXX Fill in this PostProcess function to return correct outputs
	    hf_postprocessing.PostProcess(nozzle, output);
	    
	    # Assign thermal/structural QoI if required
	    if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
	        # XXX Ensure AEROS post-processing returns correct outputs
	        AEROSpostprocessing.PostProcess(nozzle, output);           
	    
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
    
    
