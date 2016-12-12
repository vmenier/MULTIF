from multif import _mshint_module
import numpy as np
import sys

def HF_Compute_Thrust (nozzle):
	
	info = [];
	Crd  = [];
	Tri  = [];
	Tet  = [];
	Sol  = [];
	Header = [];
	
	out = _mshint_module.py_Interpolation (nozzle.exit_mesh_name, nozzle.mesh_name, nozzle.restart_name,\
	info, Crd, Tri, Tet, Sol, Header);
	
	dim    = info[0];
	NbrVer = info[1]; 
	NbrTri = info[2];
	NbrTet = info[3];
	SolSiz = info[4];
	
	
	#for i in range(1,59):
	#	idx = (i-1)*dim;
	#	idxs = (i-1)*SolSiz;
	#	print "Ver %d : %lf %lf %lf / idxs %d  Sol = %le %le ..." % (i, Crd[idx], Crd[idx+1], Crd[idx+2], idxs, Sol[idxs+1], Sol[idxs+2])
	#
	#for i in range(0,40):
	#	print "%d %le" % (i,Sol[i])
	
	# --- Get solution field indices
  
	iMach  = -1;
	iTem   = -1;
	iCons1 = -1;
	iCons2 = -1;
	iCons3 = -1;
	iCons4 = -1;
	iPres  = -1;
  
	for iFld in range(0,len(Header)):
		if Header[iFld] == 'Mach':
			iMach = iFld;
		elif Header[iFld] == 'Temperature':
			iTem = iFld;
		elif Header[iFld] == 'Conservative_1':
			iCons1 = iFld;
		elif Header[iFld] == 'Conservative_2':
			iCons2 = iFld;
		elif Header[iFld] == 'Conservative_3':
			iCons3 = iFld;
		elif Header[iFld] == 'Conservative_4':
			iCons4 = iFld;
		elif Header[iFld] == 'Pressure':
			iPres = iFld;
			
	#print "iPres %lf iCons1 %lf iCons2 %lf" % (iPres, iCons1, iCons2)
  
	# --- Compute thrust	
  
	Thrust = 0.0;
  
	P0  = nozzle.environment.P;
	M0  = nozzle.mission.mach;
	Gam = nozzle.fluid.gam;
	Rs  = nozzle.fluid.R;
	T0  = nozzle.environment.T;
	U0  = M0*np.sqrt(Gam*Rs*T0);
	
	v = np.zeros([3,3]);
	a = np.zeros(3);
	b = np.zeros(3);
	
	areatot = 0;
	
	for iTri in range(1, NbrTri) :
		idt = 3*(iTri-1);
		
		rho  = 0.0;
		rhoU = 0.0;
		Pres = 0.0;
		Mach = 0.0;
		Temp = 0.0;
		
		#print "TRI %d : %d %d %d" % (iTri,  Tri[idt+0], Tri[idt+1], Tri[idt+2])
		
		
		
		for j in range(0,3):
			iVer = int(Tri[idt+j]);
			
			idv = 3*(iVer-1);
			ids = SolSiz*(iVer-1);
			
			for d in range(0,3):
				v[j][d] = Crd[idv+d];
			
			
			#if iTri == 1:
			#	print "rho %d (%d) = %lf (ids %d)" % (j,iVer,Sol[ids+iCons1], ids)
				
			rho  = rho  + Sol[ids+iCons1];
			rhoU = rhoU + Sol[ids+iCons2];
			Pres = Pres + Sol[ids+iPres];
			Mach = Mach + Sol[ids+iMach];
			Temp = Temp + Sol[ids+iTem];	
		
		us3 = 1./3.;
		rho  = us3*rho; 
		rhoU = us3*rhoU; 
		Pres = us3*Pres; 
		Mach = us3*Mach; 
		Temp = us3*Temp; 
		
		# --- Compute triangle area
		
		for d in range(0,3):
			a[d] = v[1][d] - v[0][d];
			b[d] = v[2][d] - v[0][d];
		
		area = 0.;
		
		area = area + (a[1]*b[2]-a[2]*b[1])*(a[1]*b[2]-a[2]*b[1]);
		area = area + (a[2]*b[0]-a[0]*b[2])*(a[2]*b[0]-a[0]*b[2]);
		area = area + (a[0]*b[1]-a[1]*b[0])*(a[0]*b[1]-a[1]*b[0]);
		
		area = 0.5*np.sqrt(area);
		
		areatot = areatot+area;
				
		# --- Compute thrust
					
		U = rhoU/rho;
		
		#if iTri < 10:
		#	print "Tri %d Pres %lf rho %lf U %lf rhoU %lf U0 %lf P0 %lf area %lf" % (iTri, Pres, rho, U, rhoU, U0, P0, area)
		
		try :
			Thrust = Thrust + area*(rhoU*(U-U0)+Pres-P0);
		except :
			Thrust = -1;
		#Thrust = Thrust + area*(rhoU*(U*U0-U0)+P0*Pres-P0);
	
	print "TOTAL AREA = %lf\n" % areatot;
	return Thrust;
	