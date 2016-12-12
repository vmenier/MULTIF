#include "mshint.h"
#include "mi_calls.h"
#include "Python.h"


/* read mesh */
int copyMesh(Mesh *mesh, VMesh *VMesh) {
  pPoint       ppt;
  pTetra       pt;
  pTria        pt1;
  double       dd;
  float        fp1,fp2,fp3;
  int          inm,i,k,ref;
  char        *ptr,data[256];

	mesh->dim = VMesh->Dim;

	mesh->np = VMesh->NbrVer;
	mesh->nt = VMesh->NbrTri;
	mesh->ne = VMesh->NbrTet;
	
  /* memory allocation */
  mesh->point = (Point*)calloc(mesh->np+1,sizeof(Point));
  assert(mesh->point);

	
  /* read mesh vertices */
	
  if ( mesh->dim == 2 ) {
    /* 2d mesh */
    mesh->min[0] = mesh->min[1] =  FLT_MAX;
    mesh->max[0] = mesh->max[1] = -FLT_MAX;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
			
			ppt->c[0] = VMesh->Ver[k][0];
      ppt->c[1] = VMesh->Ver[k][1];
			
			for (i=0; i<2; i++) {
        mesh->min[i] = MI_MIN(ppt->c[i],mesh->min[i]);
        mesh->max[i] = MI_MAX(ppt->c[i],mesh->max[i]);
      }
			
    }
  }
  else {
    /* 3d mesh */
    mesh->min[0] = mesh->min[1] = mesh->min[2] =  FLT_MAX;
    mesh->max[0] = mesh->max[1] = mesh->max[2] = -FLT_MAX;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
			
			ppt->c[0] = VMesh->Ver[k][0];
      ppt->c[1] = VMesh->Ver[k][1];
			ppt->c[2] = VMesh->Ver[k][2];
			
      for (i=0; i<3; i++) {
        mesh->min[i] = MI_MIN(ppt->c[i],mesh->min[i]);
        mesh->max[i] = MI_MAX(ppt->c[i],mesh->max[i]);
      }
    }
  }
	
  /* copy mesh triangles */
  if ( mesh->nt > 0 ) {
    mesh->tria = (Tria*)calloc(mesh->nt+1,sizeof(Tria));
    assert(mesh->tria);
		
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
			pt1->v[0] =  VMesh->Tri[k][0];
			pt1->v[1] =  VMesh->Tri[k][1];
			pt1->v[2] =  VMesh->Tri[k][2];
			ref  =  VMesh->Tri[k][3];
    }
  }

	/* copy mesh tetra */
  if ( mesh->ne > 0 ) {
    mesh->tetra = (Tetra*)calloc(mesh->ne+1,sizeof(Tetra));
    assert(mesh->tetra);
		
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];

			for (i=0; i<4; i++) {
				pt->v[0] = VMesh->Tet[k][i];
			}
			
			ref = VMesh->Tet[k][4];
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
	
	sol->np = Msh->NbrVer;
  if ( !sol->np )  return(-1);
	
	
	
  sol->u = (double*)calloc(sol->np*sol->size[0],sizeof(double));
  assert(sol->u);
	
  for (k=0; k<sol->np; k++) {
		for (i=0; i<sol->size[0]; i++)
      sol->u[k*sol->size[0]+i] = Msh->Sol[(k+1)*sol->size[0]+i];
  }
		
  return(1);
}





static int parsarPy( char *MshNam, char *BakMshNam, char *BakSolNam, MIst *mist) {
	
	char    *ptr,*data;
	
	mist->info.verb  = 0;
	
	data = (char*)calloc(strlen(MshNam)+10,sizeof(char));
  strcpy(data,MshNam);
  ptr = strstr(data,".mesh");
  if ( !ptr )  strcat(data,".mesh");
  mist->mtgt.name = data;
	
  data = (char*)calloc(strlen(BakMshNam)+10,sizeof(char));
  strcpy(data,BakMshNam);
  ptr = strstr(data,".mesh");
  if ( !ptr )  strcat(data,".mesh");
  mist->msrc.name = data;
  
  data = (char*)calloc(strlen(BakSolNam)+10,sizeof(char));
  strcpy(data,BakMshNam);
  ptr = strstr(data,".sol");
  if ( !ptr )  strcat(data,".sol");
  mist->ssrc.name = data;
	
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
	
	printf("np %d, size %d\n", sol->np, sol->size[0]);
	
	for (i=0; i<(sol->np*sol->size[0]); i++) {
		printf("%lf\n", sol->u[i]);
		
		if (i==10 ) exit(1);
		PyList_Append(pySol, PyFloat_FromDouble(sol->u[i]));
	}
	
	
}


int py_Interpolation( char *MshNam, char *BakMshNam, char *BakSolNam, PyObject *pyInfo, PyObject *pyCrd, PyObject *pyTri, PyObject *pyTet, PyObject *pySol, PyObject *pyHeader) 
{
 	MIst       mist;
  int        ier, i;
	char       stim[32];
	int 			 FilTyp = -1;
	
	VMesh *Msh = NULL;
	VMesh *MshBak = NULL;
	
	char OutSol[1024];

	FILE *FilHdl=NULL;
	char *tok=NULL, *lin=NULL;

	int skip=0, SolSiz=0;
	size_t  len = 0;
	
  memset(&mist,0,sizeof(MIst));
  tminit(mist.info.ctim,TIMEMAX);
  chrono(ON,&mist.info.ctim[0]);

  /* init structure */
  memset(&mist.msrc,0,sizeof(Mesh));
  memset(&mist.mtgt,0,sizeof(Mesh));
  memset(&mist.ssrc,0,sizeof(Sol));
  memset(&mist.stgt,0,sizeof(Sol));
	
	parsarPy(MshNam, BakMshNam, BakSolNam, &mist);
	
	// --- Target mesh
	
	FilTyp = GetInputFileType(MshNam);
	if ( FilTyp == FILE_GMF ) {
  	ier = loadMesh(&mist.mtgt,mist.info.verb);
  	if ( ier <= 0 )  return(1);
	}
	else if ( FilTyp == FILE_SU2 ) {
		Msh    = SetupMeshAndSolution (MshNam, "");
		ier = copyMesh(&mist.mtgt, Msh);
		if ( ier <= 0 )  return(1);
		if ( Msh )
	 		FreeMesh(Msh);
	}
	
	// --- Source (back) mesh 
	
	FilTyp = GetInputFileType(BakMshNam);
	if ( FilTyp == FILE_GMF ) {
  	ier = loadMesh(&mist.msrc,mist.info.verb);
  	if ( ier <= 0 )  return(1);

		ier = loadSol(&mist.ssrc,mist.info.verb);
		if ( ier <= 0 )  return(1);
	}
	else if ( FilTyp == FILE_SU2 ) {

		MshBak = SetupMeshAndSolution (BakMshNam, BakSolNam);
		
		ier = copyMesh(&mist.msrc, MshBak);
		if ( ier <= 0 )  return(1);
		
		copySol(&mist.ssrc, MshBak);
		if ( ier <= 0 )  return(1);
		
		for (i=0; i<MshBak->SolSiz; i++){
			PyList_Append(pyHeader, PyString_FromString(MshBak->SolTag[i]));
		}
		
		if ( MshBak )
	 		FreeMesh(MshBak);
	}

	//sprintf(mist.msrc.name, "toto.meshb");
	//ier = saveMesh(&mist.msrc,0);
	
	sprintf(mist.ssrc.name, "toto.solb");
	ier = saveSol(&mist.ssrc,0);
	
  /* allocating memory */
  mist.stgt.dim  = mist.ssrc.dim;
  mist.stgt.ver  = mist.ssrc.ver;
  mist.stgt.np   = mist.mtgt.np;
  mist.stgt.it   = mist.ssrc.it;
  mist.stgt.time = mist.ssrc.time;
  memcpy(mist.stgt.type,mist.ssrc.type,2*sizeof(int));
  memcpy(mist.stgt.size,mist.ssrc.size,2*sizeof(int));
  memcpy(&mist.stgt.typtab[0],&mist.ssrc.typtab[0],mist.ssrc.type[0]*sizeof(int));
  memcpy(&mist.stgt.typtab[1],&mist.ssrc.typtab[1],mist.ssrc.type[1]*sizeof(int));
  if ( mist.stgt.np > 0 ) {
    mist.stgt.u = (double*)calloc(mist.stgt.np,mist.stgt.size[0]*sizeof(double));
    assert(mist.stgt.u);
  }
  
  /* set up adjacencies */
  if ( mist.msrc.dim == 2 )
    hashel_2d(&mist);
  else
    hashel_3d(&mist);
  
  chrono(OFF,&mist.info.ctim[1]);
	printim(mist.info.ctim[1].gdif,stim);
  if ( mist.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);
  
  if ( !mist.stgt.name ) {
    mist.stgt.name = (char *)calloc(128,sizeof(char));
    assert(mist.stgt.name);
    strcpy(mist.stgt.name,mist.mtgt.name);
  }
  
  /* L2 projection */
  chrono(ON,&mist.info.ctim[2]);
  if ( mist.info.verb != '0' )  fprintf(stdout,"\n ** MODULE MSHINT: %s\n",MI_VER);
  ier = MI_mshint(&mist);
  chrono(OFF,&mist.info.ctim[2]);
  if ( mist.info.verb != '0' ) {
		printim(mist.info.ctim[2].gdif,stim);
    if ( ier > 0 )  
      fprintf(stdout," ** COMPLETED: %s\n\n",stim);
    else
      fprintf(stdout," ** NOT COMPLETED!: %s\n\n",stim);
	}
  
  /* save file */
  if ( mist.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
  chrono(ON,&mist.info.ctim[3]);
  
  ier = saveSol(&mist.stgt,0);
	if ( !ier )   return(1);
  chrono(OFF,&mist.info.ctim[3]);
  if ( mist.info.verb != '0' ) {
    printim(mist.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }
  
  /* free mem */
  free(mist.msrc.point);
  if ( mist.msrc.nt )  free(mist.msrc.tria);
	if ( mist.msrc.ne )  free(mist.msrc.tetra);
  free(mist.mtgt.point);
  if ( mist.mtgt.nt )  free(mist.mtgt.tria);
	if ( mist.mtgt.ne )  free(mist.mtgt.tetra);
  if ( mist.ssrc.u )   free(mist.ssrc.u);
  
  chrono(OFF,&mist.info.ctim[0]);
  if ( mist.info.verb != '0' ) {
	  printim(mist.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

	/*
		Return values to python
	*/
	
	returnValuesToPython(&mist.mtgt, &mist.stgt, pyInfo, pyCrd, pyTri, pyTet, pySol);
  


  return(0);
	
}