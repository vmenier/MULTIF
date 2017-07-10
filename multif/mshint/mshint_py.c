#include "mshint.h"
#include "Python.h"
#include "compil.date"

Param    par_py;



/* read mesh */
int copyMesh(Mesh *mesh, VMesh *VMesh) {
  pPoint       ppt;
  pTetra       pt;
  pTria        pt1;
  double       dd;
  float        fp1,fp2,fp3;
  int          inm,i,k,ref;
  char        *ptr,data[256];

	Info        *info;

	mesh->dim = VMesh->Dim;

	mesh->np = VMesh->NbrVer;
	mesh->nt = VMesh->NbrTri;
	mesh->ne = VMesh->NbrTet;
	
	mesh->bin = 0;
	
  /* memory allocation */
  mesh->point = (pPoint)calloc(mesh->np+1,sizeof(Point));
  assert(mesh->point);
  if ( mesh->ne ) {
    mesh->tetra = (pTetra)calloc(mesh->ne+1,sizeof(Tetra));
    assert(mesh->tetra);
    mesh->adja = (int*)calloc(4*mesh->ne+5,sizeof(int));
    assert(mesh->adja);
    mesh->nt = 0;
  }
  else if ( mesh->nt ) {
    mesh->tria  = (pTria)calloc(mesh->nt+1,sizeof(Tria));
    assert(mesh->tria);
    mesh->adja = (int*)calloc(3*mesh->nt+5,sizeof(int));
    assert(mesh->adja);
  }


	
  /* read mesh vertices */
	
  if ( mesh->dim == 2 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
			ppt->c[0] = VMesh->Ver[k][0];
      ppt->c[1] = VMesh->Ver[k][1];	
    }
  }
  else {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
			
			ppt->c[0] = VMesh->Ver[k][0];
      ppt->c[1] = VMesh->Ver[k][1];
			ppt->c[2] = VMesh->Ver[k][2];
			
    }
  }
	
  /* copy mesh triangles */
  if ( mesh->nt > 0 ) {
		
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
			pt1->v[0] =  VMesh->Tri[k][0];
			pt1->v[1] =  VMesh->Tri[k][1];
			pt1->v[2] =  VMesh->Tri[k][2];
			ref  =  VMesh->Tri[k][3];
		}
		
		if ( mesh->dim == 2 ) {
			for (k=1; k<=mesh->nt; k++) {
			  for (i=0; i<3; i++) {    
    	    ppt = &mesh->point[pt1->v[i]];
			  	pt1->g[0] += ppt->c[0];
			  	pt1->g[1] += ppt->c[1];
    	    if ( !ppt->s )  ppt->s = k;
    	  }
			  pt1->g[0] /= 3.0;
			  pt1->g[1] /= 3.0;
			}	
		}

		

  }

	/* copy mesh tetra */
  if ( mesh->ne > 0 ) {
		
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];

			for (i=0; i<4; i++) {
				pt->v[i] = VMesh->Tet[k][i];
			}
			
			ref = VMesh->Tet[k][4];
    
			for (i=0; i<4; i++) {
	      ppt = &mesh->point[pt->v[i]];
				pt->g[0] += ppt->c[0];
				pt->g[1] += ppt->c[1];
				pt->g[2] += ppt->c[2];
	      if ( !ppt->s )  ppt->s = k;
	    }
			pt->g[0] *= 0.25;
			pt->g[1] *= 0.25; 
			pt->g[2] *= 0.25;
		}
		
  }

  /* find bounding box */
  info = &mesh->info;
  for (i=0; i<mesh->dim; i++) {
    info->min[i] =  1.e30;
    info->max[i] = -1.e30;
  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++) {
      if ( ppt->c[i] > info->max[i] )  info->max[i] = ppt->c[i];
      if ( ppt->c[i] < info->min[i] )  info->min[i] = ppt->c[i];
    }
  }
  info->delta = 0.0;
  for (i=0; i<mesh->dim; i++) {
    dd = fabs(info->max[i]-info->min[i]);
    if ( dd > info->delta )  info->delta = dd;
  }
  if ( info->delta < EPS1 ) {
    fprintf(stdout,"  ## Unable to scale mesh\n");
    return(0);
  }
	
  return(1);
}



/* load solution */
int copySol(Sol *sol, VMesh *Msh) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,i,inm,keyword;
  char        *ptr,data[128];

  if ( !sol->name )  return(-1);

	if ( !Msh->Sol ) {
		printf("  ## ERROR copySol : No solution in source file.\n");
		return(-1);
	}
	
	sol->dim = Msh->Dim;
	sol->ver = GmfDouble;
	
	sol->size[0] = Msh->SolSiz;
	sol->type[0] = Msh->NbrFld;
	for (i=0; i<sol->type[0]; i++) 
		sol->typtab[0][i] = Msh->FldTab[i];
	
	sol->ne = 0;

	sol->np = Msh->NbrVer;
  if ( !sol->np )  return(-1);
	
	sol->valp1 = (double*)calloc(sol->np+1,sol->size[0]*sizeof(double));
  assert(sol->valp1);
	
	int ia;
  for (k=1; k<sol->np; k++) {
		ia = (k-1)*sol->size[0] + 1;
		for (i=0; i<sol->size[0]; i++)
      sol->valp1[ia+i] = Msh->Sol[(k)*sol->size[0]+i];
  }
		
  return(1);
}


static int parsarPy( char *MshNam, char *BakMshNam, char *BakSolNam, pMesh mesh1,pSol sol1,pMesh mesh2,pSol sol2)
{
	
  char    *ptr;	


	mesh1->name = BakMshNam;
	sol1->name = BakSolNam;
	
	mesh2->name = MshNam;
  

	sol2->name = (char *)calloc(128,sizeof(char));
  assert(sol2->name);
  strcpy(sol2->name,mesh2->name);
	
  ptr = strstr(sol2->name,".mesh");
  if ( ptr )  *ptr = '\0';

  ptr = strstr(sol2->name,".su2");
  if ( ptr )  *ptr = '\0';
  //strcat(sol2->name,".solb");

}



int returnValuesToPython(Mesh *mesh, Sol *sol, PyObject *pyInfo, PyObject *pyCrd, PyObject *pyTri, PyObject *pyTet, PyObject *pySol )
{
	
	int i, j;
 	pPoint       ppt;
  pTetra       pt;
  pTria        pt1;
	
	// --- Infos : nbrver etc.
	
	PyList_Append(pyInfo, PyInt_FromLong((long) mesh->dim));
	PyList_Append(pyInfo, PyInt_FromLong((long) mesh->np));
	PyList_Append(pyInfo, PyInt_FromLong((long) mesh->nt));
	PyList_Append(pyInfo, PyInt_FromLong((long) mesh->ne));
	PyList_Append(pyInfo, PyInt_FromLong((long) sol->size[0]));
	
	// --- Ver coordinates
	
	for (i=1; i<=mesh->np; i++) {
		ppt = &mesh->point[i];
		
		PyList_Append(pyCrd, PyFloat_FromDouble(ppt->c[0]));
		PyList_Append(pyCrd, PyFloat_FromDouble(ppt->c[1]));
		
		if ( mesh->dim == 2 ) 
			PyList_Append(pyCrd, PyFloat_FromDouble(0));	
		else 
			PyList_Append(pyCrd, PyFloat_FromDouble(ppt->c[2]));	
	}
	
	// --- Triangles
	
	for (i=1; i<=mesh->nt; i++) {
		pt1 = &mesh->tria[i];
		
		PyList_Append(pyTri, PyInt_FromLong(pt1->v[0]));
		PyList_Append(pyTri, PyInt_FromLong(pt1->v[1]));
		PyList_Append(pyTri, PyInt_FromLong(pt1->v[2]));
	}
	
	// --- Tetras
	
	for (i=1; i<=mesh->ne; i++) {
		pt = &mesh->tetra[i];
		
		for (j=0; j<4; j++)
			PyList_Append(pyTet, PyInt_FromLong(pt->v[j]));
	}
	
	
	// --- Solution
	
	for (i=0; i<((sol->np+1)*sol->size[0]); i++) {
		PyList_Append(pySol, PyFloat_FromDouble(sol->valp1[i]));
	}
	
}




int py_Interpolation( char *MshNam, char *BakMshNam, char *BakSolNam, PyObject *pyInfo, PyObject *pyCrd, PyObject *pyTri, PyObject *pyTet, PyObject *pySol, PyObject *pyHeader)
{
	Mesh       mesh1,mesh2;
	Sol        sol1,sol2;
	int        i;
	char       stim[16];
	
	VMesh *Msh = NULL;
	VMesh *MshBak = NULL;
	int ier, FilTyp;
	
	fprintf(stdout,"  -- MSHINT Python module\n");
	
  
	tminit(par_py.ctim,TIMEMAX);
	chrono(ON,&par_py.ctim[0]);
  
	/* default values */
	memset(&mesh1,0,sizeof(Mesh));
	memset(&sol1,0,sizeof(Sol));
	memset(&mesh2,0,sizeof(Mesh));
	memset(&sol2,0,sizeof(Sol));
	par_py.imprim = -99;
	par_py.ddebug = 0;
	par_py.option = 1;
	par_py.ray    = 1.2;
  
	parsarPy( MshNam, BakMshNam, BakSolNam, &mesh1, &sol1, &mesh2, &sol2);
	
	// --- Source (back) mesh 
	
	FilTyp = GetInputFileType(BakMshNam);
	if ( FilTyp == FILE_GMF ) {
  	if ( !loadMesh(&mesh1,mesh1.name) )  return(1);
		if ( !loadSol(&sol1,sol1.name) )     return(1);
	}
	else if ( FilTyp == FILE_SU2 ) {

		MshBak = SetupMeshAndSolution (BakMshNam, BakSolNam);
		
		if ( !copyMesh(&mesh1, MshBak) ) return(1);		
		if ( !copySol(&sol1, MshBak) ) return(1);
		
		for (i=0; i<MshBak->SolSiz; i++){
			PyList_Append(pyHeader, PyString_FromString(MshBak->SolTag[i]));
		}
		
		if ( MshBak )
	 		FreeMesh(MshBak);
	}
	
	
	// --- Target mesh
	
	FilTyp = GetInputFileType(MshNam);
	if ( FilTyp == FILE_GMF ) {
		if ( !loadMesh(&mesh2,mesh2.name) ) return(1);
	}
	else if ( FilTyp == FILE_SU2 ) {
		Msh    = SetupMeshAndSolution (MshNam, "");
		if ( !copyMesh(&mesh2, Msh) )  return(1);
		if ( Msh )
	 		FreeMesh(Msh);
	}
	
	sol2.dim  = mesh2.dim;
  if ( sol1.np )
    sol2.np = mesh2.np;
  if ( sol1.ne ) {
		sol2.ne = mesh2.dim == 2 ? mesh2.nt : mesh2.ne;
	}
  sol2.ver  = sol1.ver;
  sol2.iter = sol1.iter;
  sol2.time = sol1.time;
  memcpy(sol2.type,sol1.type,2*sizeof(int));
  memcpy(sol2.size,sol1.size,2*sizeof(int));
  memcpy(&sol2.typtab[0],&sol1.typtab[0],sol1.type[0]*sizeof(int));
  memcpy(&sol2.typtab[1],&sol1.typtab[1],sol1.type[1]*sizeof(int));
  
	
	if ( sol2.np ) {
	  sol2.valp1 = (double*)calloc(sol2.np+1,sol2.size[0]*sizeof(double));
	  assert(sol2.valp1);
	}
	if ( sol2.ne ) {
	  sol2.valp0 = (double*)calloc(sol2.ne+1,sol2.size[1]*sizeof(double));
	  assert(sol2.valp0);
	}
	chrono(OFF,&par_py.ctim[1]);
	//if ( par_py.imprim )  stats(&mesh1,&sol1,&mesh2,&sol2);
	printim(par_py.ctim[1].gdif,stim);
	fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
    
	//if ( !saveSol(&sol1,"test") )   return(1);
	
	/* analysis */
	chrono(ON,&par_py.ctim[2]);
	if ( par_py.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
	for (i=0; i<mesh1.dim; i++) {
	  mesh1.info.min[i] = I_MIN(mesh1.info.min[i],mesh2.info.min[i]);
	  mesh1.info.max[i] = I_MAX(mesh1.info.max[i],mesh2.info.max[i]);
	}
	mesh1.info.delta = I_MAX(mesh1.info.delta,mesh2.info.delta);
	memcpy(&mesh2.info,&mesh1.info,sizeof(Info));
	if ( !scaleMesh(&mesh1,&sol1) )  return(1);
	if ( !scaleMesh(&mesh2,0) )      return(1);
	if ( mesh1.dim == 2 ) {
		if ( !hashelt_2d(&mesh1) )          return(1);
	}
	else {
		if ( !hashelt_3d(&mesh1) )          return(1);
	}
	chrono(OFF,&par_py.ctim[2]);
	printim(par_py.ctim[2].gdif,stim);
	if ( par_py.imprim )
	  fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  
	/* interpolation */
	if ( par_py.imprim )   fprintf(stdout,"\n  -- PHASE 2 : INTERPOLATION\n");
	chrono(ON,&par_py.ctim[3]);
	if ( !mshin1(&mesh1,&sol1,&mesh2,&sol2) )  return(1);
	chrono(OFF,&par_py.ctim[3]);
	printim(par_py.ctim[3].gdif,stim);
	if ( par_py.imprim )
	  fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s sec.\n",stim);
  
	if ( !unscaleMesh(&mesh2,&sol2) )  return(1);
	
	if ( !saveSol(&sol2,sol2.name) )   return(1);
	
	
	
  
	/*
		Return values to python
	*/
	
	returnValuesToPython(&mesh2, &sol2, pyInfo, pyCrd, pyTri, pyTet, pySol);
	

	/*
		Free memory
	*/

	free(mesh1.point);
	free(mesh1.adja);
	free(mesh2.point);
	free(mesh2.adja);
	if ( mesh1.dim == 3 ) {
	  free(mesh1.tetra);
	  free(mesh2.tetra);
	}
	else {
	  free(mesh1.tria);
	  free(mesh2.tria);
	}
	if ( sol2.np )
	  free(sol2.valp1);
	if ( sol2.ne )
		free(sol2.valp0);
  
	return(0);
	
}
