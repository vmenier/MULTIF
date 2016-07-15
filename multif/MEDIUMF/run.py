#from .. import CFD
from .. import SU2
from meshgeneration import *
from runSU2 import *
from postprocessing import *

class Solver_Options:
	def __init__(self):
		pass


def Run( nozzle ):
	
	# --- Check SU2 version
	
	CheckSU2Version(nozzle);
		
	# --- 
	
	solver_options = Solver_Options();
	
	solver_options.Mach = nozzle.mission.mach;
	solver_options.Pres = nozzle.environment.P;
	solver_options.Temp = nozzle.environment.T;
	
	solver_options.InletPstag = nozzle.inlet.Pstag;
	solver_options.InletTstag = nozzle.inlet.Tstag;
	
	solver_options.LocalRelax = nozzle.LocalRelax;
	
	solver_options.NbrIte = 500;
	
	GenerateNozzleMesh(nozzle);
	
	config = SetupConfig(solver_options);
	
	info = SU2.run.CFD(config);
	
	CheckConvergence(config);