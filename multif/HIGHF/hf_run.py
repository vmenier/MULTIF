from .. import SU2
from hf_meshgeneration import *
from hf_runSU2 import *
from hf_postprocessing import *


def CheckOptions (nozzle):
	
	print "Check options"
	if nozzle.Dim != '3D':
	    sys.stderr.write("\n  ## ERROR : High-fidelity is 3D.\n\n");
	    sys.exit(0);
	
def Run( nozzle ):
	
	# --- Check SU2 version
	
	#CheckSU2Version(nozzle);
	    
	# --- Check fidelity level
	
	CheckOptions (nozzle);
	
	curDir = os.path.dirname(os.path.realpath(__file__));
	
	if nozzle.runDir != '':
	    os.chdir(nozzle.runDir);
	
	# --- Generate mesh
	
	#HF_GenerateMesh(nozzle);
	
	# --- Run CFD
	
	#HF_runSU2 (nozzle);
	
	# --- Run AEROS  
	# TODO
	
	# --- Postprocessing
	
	nozzle.thrust = HF_Compute_Thrust (nozzle);
	
	#print "THRUST = %lf\n" % Thrust;
	
	sys.stdout.write("\n  -- Info : Result directory :  %s\n\n" % nozzle.runDir);
	
	
	
	