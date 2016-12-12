#include "mshint.h"
#include "Python.h"
#include "compil.date"

Param    par;


//

//
///* read mesh */
//int copyMesh(Mesh *mesh, VMesh *VMesh) {
//  pPoint       ppt;
//  pTetra       pt;
//  pTria        pt1;
//  double       dd;
//  float        fp1,fp2,fp3;
//  int          inm,i,k,ref;
//  char        *ptr,data[256];
//
//	mesh->dim = VMesh->Dim;
//
//	mesh->np = VMesh->NbrVer;
//	mesh->nt = VMesh->NbrTri;
//	mesh->ne = VMesh->NbrTet;
//	
//  /* memory allocation */
//  mesh->point = (Point*)calloc(mesh->np+1,sizeof(Point));
//  assert(mesh->point);
//
//	
//  /* read mesh vertices */
//	
//  if ( mesh->dim == 2 ) {
//    /* 2d mesh */
//    //mesh->min[0] = mesh->min[1] =  FLT_MAX;
//    //mesh->max[0] = mesh->max[1] = -FLT_MAX;
//    for (k=1; k<=mesh->np; k++) {
//      ppt = &mesh->point[k];
//			
//			ppt->c[0] = VMesh->Ver[k][0];
//      ppt->c[1] = VMesh->Ver[k][1];
//			
//			//for (i=0; i<2; i++) {
//      //  mesh->min[i] = MI_MIN(ppt->c[i],mesh->min[i]);
//      //  mesh->max[i] = MI_MAX(ppt->c[i],mesh->max[i]);
//      //}
//			
//    }
//  }
//  else {
//    /* 3d mesh */
//    //mesh->min[0] = mesh->min[1] = mesh->min[2] =  FLT_MAX;
//    //mesh->max[0] = mesh->max[1] = mesh->max[2] = -FLT_MAX;
//    for (k=1; k<=mesh->np; k++) {
//      ppt = &mesh->point[k];
//			
//			ppt->c[0] = VMesh->Ver[k][0];
//      ppt->c[1] = VMesh->Ver[k][1];
//			ppt->c[2] = VMesh->Ver[k][2];
//			
//      //for (i=0; i<3; i++) {
//      //  mesh->min[i] = MI_MIN(ppt->c[i],mesh->min[i]);
//      //  mesh->max[i] = MI_MAX(ppt->c[i],mesh->max[i]);
//      //}
//    }
//  }
//	
//  /* copy mesh triangles */
//  if ( mesh->nt > 0 ) {
//    mesh->tria = (Tria*)calloc(mesh->nt+1,sizeof(Tria));
//    assert(mesh->tria);
//		
//    for (k=1; k<=mesh->nt; k++) {
//      pt1 = &mesh->tria[k];
//			pt1->v[0] =  VMesh->Tri[k][0];
//			pt1->v[1] =  VMesh->Tri[k][1];
//			pt1->v[2] =  VMesh->Tri[k][2];
//			ref  =  VMesh->Tri[k][3];
//    }
//  }
//
//	/* copy mesh tetra */
//  if ( mesh->ne > 0 ) {
//    mesh->tetra = (Tetra*)calloc(mesh->ne+1,sizeof(Tetra));
//    assert(mesh->tetra);
//		
//    for (k=1; k<=mesh->ne; k++) {
//      pt = &mesh->tetra[k];
//
//			for (i=0; i<4; i++) {
//				pt->v[0] = VMesh->Tet[k][i];
//			}
//			
//			ref = VMesh->Tet[k][4];
//    }
//  }
//	
//  return(1);
//}
//
//
//
///* load solution */
//int copySol(Sol *sol, VMesh *Msh) {
//  double       dbuf[GmfMaxTyp];
//  float        fbuf[GmfMaxTyp];
//  int          k,i,inm,keyword;
//  char        *ptr,data[128];
//
//  if ( !sol->name )  return(-1);
//
//	if ( !Msh->Sol ) {
//		printf("  ## ERROR copySol : No solution in source file.\n");
//		return(-1);
//	}
//	
//	sol->dim = Msh->Dim;
//	sol->ver = GmfDouble;
//	
//	sol->size[0] = Msh->SolSiz;
//	sol->type[0] = Msh->NbrFld;
//	for (i=0; i<sol->type[0]; i++) 
//		sol->typtab[0][i] = Msh->FldTab[i];
//	
//	sol->np = Msh->NbrVer;
//  if ( !sol->np )  return(-1);
//	
//	
//	
//  sol->u = (double*)calloc(sol->np*sol->size[0],sizeof(double));
//  assert(sol->u);
//	
//  for (k=0; k<sol->np; k++) {
//		for (i=0; i<sol->size[0]; i++)
//      sol->u[k*sol->size[0]+i] = Msh->Sol[(k+1)*sol->size[0]+i];
//  }
//		
//  return(1);
//}
//
//


//
//static int parsarPy( char *MshNam, char *BakMshNam, char *BakSolNam, MIst *mist) {
//	
//	char    *ptr,*data;
//	
//	mist->info.verb  = 0;
//	
//	data = (char*)calloc(strlen(MshNam)+10,sizeof(char));
//  strcpy(data,MshNam);
//  ptr = strstr(data,".mesh");
//  if ( !ptr )  strcat(data,".mesh");
//  mist->mtgt.name = data;
//	
//  data = (char*)calloc(strlen(BakMshNam)+10,sizeof(char));
//  strcpy(data,BakMshNam);
//  ptr = strstr(data,".mesh");
//  if ( !ptr )  strcat(data,".mesh");
//  mist->msrc.name = data;
//  
//  data = (char*)calloc(strlen(BakSolNam)+10,sizeof(char));
//  strcpy(data,BakMshNam);
//  ptr = strstr(data,".sol");
//  if ( !ptr )  strcat(data,".sol");
//  mist->ssrc.name = data;
//	
//}
//
//
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
		printf("%lf\n", sol->valp1[i]);
		
		if (i==10 ) exit(1);
		PyList_Append(pySol, PyFloat_FromDouble(sol->valp1[i]));
	}
	
	
}


//int py_Interpolation( char *MshNam, char *BakMshNam, char *BakSolNam, PyObject *pyInfo, PyObject *pyCrd, PyObject *pyTri, PyObject *pyTet, PyObject *pySol, PyObject *pyHeader) 
//{
// 	MIst       mist;
//  int        ier, i;
//	char       stim[32];
//	int 			 FilTyp = -1;
//	
//	VMesh *Msh = NULL;
//	VMesh *MshBak = NULL;
//	
//	char OutSol[1024];
//
//	FILE *FilHdl=NULL;
//	char *tok=NULL, *lin=NULL;
//
//	int skip=0, SolSiz=0;
//	size_t  len = 0;
//	
//  memset(&mist,0,sizeof(MIst));
//  tminit(mist.info.ctim,TIMEMAX);
//  chrono(ON,&mist.info.ctim[0]);
//
//  /* init structure */
//  memset(&mist.msrc,0,sizeof(Mesh));
//  memset(&mist.mtgt,0,sizeof(Mesh));
//  memset(&mist.ssrc,0,sizeof(Sol));
//  memset(&mist.stgt,0,sizeof(Sol));
//	
//	parsarPy(MshNam, BakMshNam, BakSolNam, &mist);
//	
//	// --- Target mesh
//	
//	FilTyp = GetInputFileType(MshNam);
//	if ( FilTyp == FILE_GMF ) {
//  	ier = loadMesh(&mist.mtgt,mist.info.verb);
//  	if ( ier <= 0 )  return(1);
//	}
//	else if ( FilTyp == FILE_SU2 ) {
//		Msh    = SetupMeshAndSolution (MshNam, "");
//		ier = copyMesh(&mist.mtgt, Msh);
//		if ( ier <= 0 )  return(1);
//		if ( Msh )
//	 		FreeMesh(Msh);
//	}
//	
//	// --- Source (back) mesh 
//	
//	FilTyp = GetInputFileType(BakMshNam);
//	if ( FilTyp == FILE_GMF ) {
//  	ier = loadMesh(&mist.msrc,mist.info.verb);
//  	if ( ier <= 0 )  return(1);
//
//		ier = loadSol(&mist.ssrc,mist.info.verb);
//		if ( ier <= 0 )  return(1);
//	}
//	else if ( FilTyp == FILE_SU2 ) {
//
//		MshBak = SetupMeshAndSolution (BakMshNam, BakSolNam);
//		
//		ier = copyMesh(&mist.msrc, MshBak);
//		if ( ier <= 0 )  return(1);
//		
//		copySol(&mist.ssrc, MshBak);
//		if ( ier <= 0 )  return(1);
//		
//		for (i=0; i<MshBak->SolSiz; i++){
//			PyList_Append(pyHeader, PyString_FromString(MshBak->SolTag[i]));
//		}
//		
//		if ( MshBak )
//	 		FreeMesh(MshBak);
//	}
//
//	//sprintf(mist.msrc.name, "toto.meshb");
//	//ier = saveMesh(&mist.msrc,0);
//	
//	sprintf(mist.ssrc.name, "toto.solb");
//	ier = saveSol(&mist.ssrc,0);
//	
//  /* allocating memory */
//  mist.stgt.dim  = mist.ssrc.dim;
//  mist.stgt.ver  = mist.ssrc.ver;
//  mist.stgt.np   = mist.mtgt.np;
//  mist.stgt.it   = mist.ssrc.it;
//  mist.stgt.time = mist.ssrc.time;
//  memcpy(mist.stgt.type,mist.ssrc.type,2*sizeof(int));
//  memcpy(mist.stgt.size,mist.ssrc.size,2*sizeof(int));
//  memcpy(&mist.stgt.typtab[0],&mist.ssrc.typtab[0],mist.ssrc.type[0]*sizeof(int));
//  memcpy(&mist.stgt.typtab[1],&mist.ssrc.typtab[1],mist.ssrc.type[1]*sizeof(int));
//  if ( mist.stgt.np > 0 ) {
//    mist.stgt.u = (double*)calloc(mist.stgt.np,mist.stgt.size[0]*sizeof(double));
//    assert(mist.stgt.u);
//  }
//  
//  /* set up adjacencies */
//  if ( mist.msrc.dim == 2 )
//    hashel_2d(&mist);
//  else
//    hashel_3d(&mist);
//  
//  chrono(OFF,&mist.info.ctim[1]);
//	printim(mist.info.ctim[1].gdif,stim);
//  if ( mist.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);
//  
//  if ( !mist.stgt.name ) {
//    mist.stgt.name = (char *)calloc(128,sizeof(char));
//    assert(mist.stgt.name);
//    strcpy(mist.stgt.name,mist.mtgt.name);
//  }
//  
//  /* L2 projection */
//  chrono(ON,&mist.info.ctim[2]);
//  if ( mist.info.verb != '0' )  fprintf(stdout,"\n ** MODULE MSHINT: %s\n",MI_VER);
//  ier = MI_mshint(&mist);
//  chrono(OFF,&mist.info.ctim[2]);
//  if ( mist.info.verb != '0' ) {
//		printim(mist.info.ctim[2].gdif,stim);
//    if ( ier > 0 )  
//      fprintf(stdout," ** COMPLETED: %s\n\n",stim);
//    else
//      fprintf(stdout," ** NOT COMPLETED!: %s\n\n",stim);
//	}
//  
//  /* save file */
//  if ( mist.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
//  chrono(ON,&mist.info.ctim[3]);
//  
//  ier = saveSol(&mist.stgt,0);
//	if ( !ier )   return(1);
//  chrono(OFF,&mist.info.ctim[3]);
//  if ( mist.info.verb != '0' ) {
//    printim(mist.info.ctim[3].gdif,stim);
//    fprintf(stdout," - COMPLETED: %s\n",stim);
//  }
//  
//  /* free mem */
//  free(mist.msrc.point);
//  if ( mist.msrc.nt )  free(mist.msrc.tria);
//	if ( mist.msrc.ne )  free(mist.msrc.tetra);
//  free(mist.mtgt.point);
//  if ( mist.mtgt.nt )  free(mist.mtgt.tria);
//	if ( mist.mtgt.ne )  free(mist.mtgt.tetra);
//  if ( mist.ssrc.u )   free(mist.ssrc.u);
//  
//  chrono(OFF,&mist.info.ctim[0]);
//  if ( mist.info.verb != '0' ) {
//	  printim(mist.info.ctim[0].gdif,stim);
//    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
//  }
//
//	/*
//		Return values to python
//	*/
//	
//	returnValuesToPython(&mist.mtgt, &mist.stgt, pyInfo, pyCrd, pyTri, pyTet, pySol);
//  
//
//
//  return(0);
//	
//}
//


int py_Interpolation( char *MshNam, char *BakMshNam, char *BakSolNam, PyObject *pyInfo, PyObject *pyCrd, PyObject *pyTri, PyObject *pyTet, PyObject *pySol, PyObject *pyHeader)
{
	Mesh       mesh1,mesh2;
	Sol        sol1,sol2;
	int        i;
	char       stim[16];
	
	fprintf(stdout,"  -- MSHINT, Release %s (%s) \n",I_VER,I_REL);
	fprintf(stdout,"     %s\n",I_CPY);
	fprintf(stdout,"    %s\n",COMPIL);
	
	///* trap exceptions */
	//signal(SIGABRT,excfun);
	//signal(SIGFPE,excfun);
	//signal(SIGILL,excfun);
	//signal(SIGSEGV,excfun);
	//signal(SIGTERM,excfun);
	//signal(SIGINT,excfun);
	//atexit(endcod);
	
	tminit(par.ctim,TIMEMAX);
	chrono(ON,&par.ctim[0]);
	
	/* default values */
	memset(&mesh1,0,sizeof(Mesh));
	memset(&sol1,0,sizeof(Sol));
	memset(&mesh2,0,sizeof(Mesh));
	memset(&sol2,0,sizeof(Sol));
	par.imprim = -99;
	par.ddebug = 0;
	par.option = 1;
	par.ray    = 1.2;
	
	mesh1.name = BakMshNam;
	sol1.name  = BakSolNam;
	mesh2.name = MshNam;
	
	
	
	//if ( !parsar(argc,argv,&mesh1,&sol1,&mesh2,&sol2) )  return(1);
	
	/* load data */
	if ( par.imprim )   fprintf(stdout,"\n  -- INPUT DATA\n");
	chrono(ON,&par.ctim[1]);
	if ( !loadMesh(&mesh1,mesh1.name) )  return(1);
	if ( !loadMesh(&mesh2,mesh2.name) )  return(1);
	if ( !setfunc(mesh1.dim) )  return(1);
	
	if ( par.option == 1 ) {
	  if ( !loadSol(&sol1,sol1.name) )  return(1);
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
	}
	else {
	  sol2.dim  = mesh2.dim;
	  sol2.np   = mesh2.np;
	  sol2.ver  = GmfFloat;
	  sol2.type[0]   = 1;
	  sol2.size[0]   = 1;
	  sol2.typtab[0][0] = GmfSca;
	}
	if ( sol2.np ) {
	  sol2.valp1 = (double*)calloc(sol2.np+1,sol2.size[0]*sizeof(double));
	  assert(sol2.valp1);
	}
	if ( sol2.ne ) {
	  sol2.valp0 = (double*)calloc(sol2.ne+1,sol2.size[1]*sizeof(double));
	  assert(sol2.valp0);
	}
	chrono(OFF,&par.ctim[1]);
	if ( par.imprim )  stats(&mesh1,&sol1,&mesh2,&sol2);
	printim(par.ctim[1].gdif,stim);
	fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
	
	fprintf(stdout,"\n  %s\n   MODULE MSHINT-LJLL : %s (%s)\n  %s\n",
	        I_STR,I_VER,I_REL,I_STR);
	
	/* analysis */
	chrono(ON,&par.ctim[2]);
	if ( par.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
	for (i=0; i<mesh1.dim; i++) {
	  mesh1.info.min[i] = I_MIN(mesh1.info.min[i],mesh2.info.min[i]);
	  mesh1.info.max[i] = I_MAX(mesh1.info.max[i],mesh2.info.max[i]);
	}
	mesh1.info.delta = I_MAX(mesh1.info.delta,mesh2.info.delta);
	memcpy(&mesh2.info,&mesh1.info,sizeof(Info));
	if ( !scaleMesh(&mesh1,&sol1) )  return(1);
	if ( !scaleMesh(&mesh2,0) )      return(1);
	if ( !hashelt(&mesh1) )          return(1);
	chrono(OFF,&par.ctim[2]);
	printim(par.ctim[2].gdif,stim);
	if ( par.imprim )
	  fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
	
	/* interpolation */
	if ( par.imprim )   fprintf(stdout,"\n  -- PHASE 2 : INTERPOLATION\n");
	chrono(ON,&par.ctim[3]);
	if ( !mshin1(&mesh1,&sol1,&mesh2,&sol2) )  return(1);
	chrono(OFF,&par.ctim[3]);
	printim(par.ctim[3].gdif,stim);
	if ( par.imprim )
	  fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s sec.\n",stim);
	
	/* save file */
	if ( par.imprim )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh2.name);
	chrono(ON,&par.ctim[1]);
	if ( !unscaleMesh(&mesh2,&sol2) )  return(1);
	if ( !saveSol(&sol2,sol2.name) )   return(1);
	chrono(OFF,&par.ctim[1]);
	if ( par.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");
	
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
