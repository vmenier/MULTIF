from .. import CFD

def Run( nozzle ):
	
	print "Mesh generation"
	
	CFD.GenerateNozzleMesh(nozzle);
	
	
	