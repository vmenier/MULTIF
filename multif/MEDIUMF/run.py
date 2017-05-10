from .. import SU2
from meshgeneration import *
from runSU2 import *
from postprocessing import *

try:
    from runAEROS import *
except ImportError:
    print 'Error importing all functions from runAEROS.\n'

def CheckOptions (nozzle):
    
		print "Check options"
    #if nozzle.Dim == 3 :
    #    sys.stderr.write("\n  ## ERROR : Only 2D axisymmetric simulations are available for now.\n\n");
    #    sys.exit(0);
    
    #if nozzle.method == 'RANS':
    #    sys.stderr.write("\n  ## ERROR : Only Euler simulations are available for now.\n\n");
    #    sys.exit(0);
    

def Run( nozzle, output = 'verbose' ):

    # Obtain mass and volume
    if nozzle.GetOutput['MASS'] == 1 or nozzle.GetOutput['VOLUME'] == 1:
        volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle) 
        nozzle.mass = np.sum(mass)
        nozzle.volume = np.sum(volume)
        
        # Calculate mass gradients if necessary
        if nozzle.mass_gradients == 'YES' or nozzle.output_gradients == 'YES':
            if ( nozzle.gradients_method == 'ADJOINT' ):
                nozzle.mass_grad = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-6);
            elif ( nozzle.gradients_method == 'FINITE_DIFF' ):
                sys.stderr.write('\n ## ERROR : No user-defined finite difference step has been defined\n');
                sys.exit(1);
                user_step = 1e-3;
                nozzle.mass_grad = nozzlemod.geometry.calcMassGradientsFD(nozzle,user_step);
            else:
			    sys.stderr.write("  ## ERROR : Unknown gradients computation method.\n");
			    sys.exit(1);
			    
            if nozzle.output_gradients == 'YES':
			    np.savetxt(nozzle.output_gradients_filename, nozzle.mass_grad, delimiter='\n')
        
    # Obtain mass of wall only if requested
    if nozzle.GetOutput['MASS_WALL_ONLY'] == 1:
        n_layers = len(nozzle.wall.layer);
        nozzle.mass_wall_only = np.sum(performanceTuple[1][:n_layers]);    

    # Run aero-thermal-structural analysis if other QoI are requested
    otherQoI = ['MAX_TOTAL_STRESS','KS_TEMPERATURE','KS_TOTAL_STRESS','MAX_TEMP_RATIO', \
        'KS_FAILURE_CRITERIA','WALL_PRESSURE','PN_TOTAL_STRESS','PN_TEMP_RATIO', \
        'PN_TEMPERATURE','VELOCITY','MAX_THERMAL_STRESS','MAX_TEMPERATURE', \
        'WALL_TEMPERATURE','MAX_FAILURE_CRITERIA','PRESSURE','PN_FAILURE_CRITERIA', \
        'THRUST','KS_TEMP_RATIO','MAX_MECHANICAL_STRESS']
    nRequested = 0
    for qoi in otherQoI:
        nRequested += np.sum(nozzle.GetOutput[qoi])
        
    # --- Check SU2 version	
    
    CheckSU2Version(nozzle);
        	
    if nRequested > 0:	
        	    
        	# --- Check fidelity level
        	
        	CheckOptions (nozzle);
        	
        	curDir = os.path.dirname(os.path.realpath(__file__));
        	
        	if nozzle.runDir != '':
        	    os.chdir(nozzle.runDir);
        	
        	# --- Run CFD
        	
        	runSU2 (nozzle);
        	
        	## --- Thrust gradient verification?
        	#
        	#thrust_grad = Compute_Thrust_Gradients_FD (nozzle);
        	#np.savetxt('thrust_grad_FD.dat', thrust_grad, delimiter='\n')  
        	
        	# --- Run AEROS  
        	
        	nozzle.runAEROS = 0;
        	if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
        	    nozzle.runAEROS = 1;
        	
        	    try:
        	        #from runAEROS import *
        	        if output == 'verbose':      
        	            print "SUCCESS IMPORTING AEROS"
        	    except ImportError:
        	            nozzle.runAEROS = 0;
        	            pass;
        	
        	if output == 'verbose':
        	    print "RUNAEROS = %d" % nozzle.runAEROS;
        	
        	print output;
        	if nozzle.runAEROS == 1:
        	    runAEROS(nozzle, output);
        	else :
        	    sys.stdout.write('  -- Info: Skip call to AEROS.\n');
        	
        	# --- Postprocessing
        	
        	PostProcessing(nozzle);
        	
        	if output == 'verbose':
        	    sys.stdout.write("\n  -- Info : Result directory :  %s\n\n" % nozzle.runDir);
        	
    # --- Output results
    	
    if nozzle.outputFormat == 'PLAIN':
        nozzle.WriteOutputFunctions_Plain();
    else:
        nozzle.WriteOutputFunctions_Dakota();
	
	
	
