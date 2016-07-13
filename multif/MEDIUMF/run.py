#from .. import CFD
from .. import SU2
from meshgeneration import *
from runSU2 import *

def Run( nozzle ):
	
	# --- Check SU2 version
	
	CheckSU2Version(nozzle);
	
	print "-- Mesh generation"
	
	GenerateNozzleMesh(nozzle);
	
	config = SetupConfig(nozzle);
	
	SU2.run.CFD(config);
	
	