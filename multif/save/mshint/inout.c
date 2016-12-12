#include "mshint.h"
#include "mi_calls.h"


/* read mesh */
int loadMesh(Mesh *mesh,char verb) {
  pPoint       ppt;
  pTetra       pt;
  pTria        pt1;
  double       dd;
  float        fp1,fp2,fp3;
  int          inm,i,k,ref;
  char        *ptr,data[256];

  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        if ( verb != '0' )  fprintf(stderr," # %s: file not found.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    if ( verb != '0' )  fprintf(stderr," # %s: file not found.\n",data);
    return(0);
  }

  if ( verb != '0' )  fprintf(stdout,"    %s:",data);

  mesh->np = GmfStatKwd(inm,GmfVertices);
  mesh->nt = GmfStatKwd(inm,GmfTriangles);
  mesh->ne = GmfStatKwd(inm,GmfTetrahedra);
  if ( !mesh->np ) {
    if ( verb != '0' )  fprintf(stdout,"\n # missing data\n");
    return(0);
  }

  /* memory allocation */
  mesh->point = (Point*)calloc(mesh->np+1,sizeof(Point));
  assert(mesh->point);
  /* read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  if ( mesh->dim == 2 ) {
    /* 2d mesh */
    mesh->min[0] = mesh->min[1] =  FLT_MAX;
    mesh->max[0] = mesh->max[1] = -FLT_MAX;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ref);
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
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ref);
      for (i=0; i<3; i++) {
        mesh->min[i] = MI_MIN(ppt->c[i],mesh->min[i]);
        mesh->max[i] = MI_MAX(ppt->c[i],mesh->max[i]);
      }
    }
  }
  /* read mesh triangles */
  if ( mesh->nt > 0 ) {
    mesh->tria = (Tria*)calloc(mesh->nt+1,sizeof(Tria));
    assert(mesh->tria);
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
    }
  }
  if ( mesh->ne > 0 ) {
    mesh->tetra = (Tetra*)calloc(mesh->ne+1,sizeof(Tetra));
    assert(mesh->tetra);
    GmfGotoKwd(inm,GmfTetrahedra);
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
    }
  }
  GmfCloseMesh(inm);

  if ( verb != '0' ) {
    fprintf(stdout," %d vertices",mesh->np);
    if ( mesh->nt )  fprintf(stdout,", %d triangles",mesh->nt);
    if ( mesh->ne )  fprintf(stdout,", %d tetrahedra",mesh->ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* load solution */
int loadSol(Sol *sol,char verb) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,i,inm,keyword;
  char        *ptr,data[128];

  if ( !sol->name )  return(-1);
  strcpy(data,sol->name);

  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  
  /* look for data file */
  ptr = strstr(data,".sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".solb");
    inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim);
    if ( !inm ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim);
    }
  }
  if ( !inm )  return(-1);

  sol->np = GmfStatKwd(inm,GmfSolAtVertices,sol->type,sol->size,sol->typtab[0]);
  if ( !sol->np )  return(-1);

  if ( verb != '0' )  fprintf(stdout,"    %s:",data);

  sol->u = (double*)calloc(sol->np*sol->size[0],sizeof(double));
  assert(sol->u);

  GmfGotoKwd(inm,GmfSolAtVertices);
  for (k=0; k<sol->np; k++) {
    if ( sol->ver == GmfFloat ) {
      GmfGetLin(inm,GmfSolAtVertices,fbuf);
      for (i=0; i<sol->size[0]; i++)
        sol->u[k*sol->size[0]+i] = fbuf[i];
    }
    else {
      GmfGetLin(inm,GmfSolAtVertices,dbuf);
      for (i=0; i<sol->size[0]; i++)
        sol->u[k*sol->size[0]+i] = dbuf[i];
    }
  }

  if ( GmfStatKwd(inm,GmfIterations) ) {
    GmfGotoKwd(inm,GmfIterations);
    GmfGetLin(inm,GmfIterations,&sol->it);
	}
  if ( GmfStatKwd(inm,GmfTime) ) {
		GmfGotoKwd(inm,GmfTime);
    if ( sol->ver == GmfFloat ) {
      GmfGetLin(inm,GmfTime,fbuf);
      sol->time = (double)fbuf[0];
    }
    else {
      GmfGetLin(inm,GmfTime,&sol->time);
    }
	}
  GmfCloseMesh(inm);

  if ( verb != '0' ) {
    fprintf(stdout," %d data\n",sol->np);
  }

  return(1);
}


/* save solution */
int saveSol(Sol *sol,char verb) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,i,ia,outm,keyword;
  char        *ptr,data[128];

  strcpy(data,sol->name);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".sol");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".solb");
  }

  if ( !(outm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    if ( verb != '0' )  fprintf(stderr," # unable to open %s\n",data);
    return(0);
  }
  if ( verb != '0' )  fprintf(stdout,"    %s:",data);

  /* write sol */
  GmfSetKwd(outm,GmfSolAtVertices,sol->np,sol->type[0],sol->typtab[0]);
  for (k=0; k<sol->np; k++) {
    if ( sol->ver == GmfFloat ) {
      for (i=0; i<sol->size[0]; i++)
        fbuf[i] = sol->u[k*sol->size[0]+i];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
    else {
      for (i=0; i<sol->size[0]; i++)
        dbuf[i] = sol->u[k*sol->size[0]+i];
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }

  if ( sol->it > 0 ) {
    GmfSetKwd(outm,GmfIterations);
    GmfSetLin(outm,GmfIterations,sol->it);
  }
	if ( sol->time > 0.0 ) {
    GmfSetKwd(outm,GmfTime);
    if ( sol->ver == GmfFloat ) {
			fbuf[0] = sol->time;
      GmfSetLin(outm,GmfTime,fbuf[0]);
		}
		else
			GmfSetLin(outm,GmfTime,sol->time);
	}
  GmfCloseMesh(outm);

  if ( verb != '0' )  fprintf(stdout," %d data\n",sol->np);

  return(1);
}

