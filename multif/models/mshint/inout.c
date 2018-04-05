#include "mshint.h"

extern Param par;


/* read mesh data */
int loadMesh(pMesh mesh,char *filename) {
  pPoint       ppt;
  pTetra       pt;
  pTria        pt1;
  Info        *info;
  double       dd;
  float        fp1,fp2,fp3;
  int          inm,i,k,ref;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  mesh->bin = 0;
  if ( !ptr ) {
    strcat(data,".meshb");
    if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
      mesh->bin = 0;
    }
  }
  else if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(par.imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  mesh->np = GmfStatKwd(inm,GmfVertices);
  mesh->nt = GmfStatKwd(inm,GmfTriangles);
  mesh->ne = GmfStatKwd(inm,GmfTetrahedra);
  if ( !mesh->np ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }

  /* mem alloc */
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
  GmfGotoKwd(inm,GmfVertices);
  if ( mesh->dim == 2 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ref);
    }
  }
  else {
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
    }
  }

  /* read mesh triangles */
  GmfGotoKwd(inm,GmfTriangles);
  if ( mesh->dim == 2 ) {
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
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
  else {
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
    }
  }
  
  /* read mesh tetrahedra */
  GmfGotoKwd(inm,GmfTetrahedra);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
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

  GmfCloseMesh(inm);

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
int loadSol(pSol sol,char *filename) {
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  int          k,i,ia,inm,keyword;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".solb");
  if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
    ptr  = strstr(data,".solb");
    *ptr = '\0';
    strcat(data,".solb");
    if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
      return(0);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  sol->np = GmfStatKwd(inm,GmfSolAtVertices,&sol->type[0],&sol->size[0],sol->typtab[0]);
  keyword = sol->dim == 2 ? GmfSolAtTriangles : GmfSolAtTetrahedra;
  sol->ne = GmfStatKwd(inm,keyword,&sol->type[1],&sol->size[1],sol->typtab[1]);
  if ( sol->np + sol->ne == 0 ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    return(0);
  }

  /* read mesh solutions */
  if ( sol->np ) {
    sol->valp1 = (double*)calloc(sol->np+1,sol->size[0]*sizeof(double));
    assert(sol->valp1);

    GmfGotoKwd(inm,GmfSolAtVertices);
    for (k=1; k<=sol->np; k++) {
      if ( sol->ver == GmfFloat ) {
        GmfGetLin(inm,GmfSolAtVertices,fbuf);
        ia = (k-1)*sol->size[0] + 1;
        for (i=0; i<sol->size[0]; i++)
          sol->valp1[ia+i] = fbuf[i];
      }
      else {
        GmfGetLin(inm,GmfSolAtVertices,dbuf);
        ia = (k-1)*sol->size[0] + 1;
        for (i=0; i<sol->size[0]; i++)
          sol->valp1[ia+i] = dbuf[i];
      }
    }
  }

  if ( sol->ne ) {
    sol->valp0 = (double*)calloc(sol->ne+1,sol->size[1]*sizeof(double));
    assert(sol->valp0);

    GmfGotoKwd(inm,keyword);
    for (k=1; k<=sol->ne; k++) {
      if ( sol->ver == GmfFloat ) {
        GmfGetLin(inm,keyword,fbuf);
        ia = (k-1)*sol->size[1] + 1;
        for (i=0; i<sol->size[1]; i++)
          sol->valp0[ia+i] = fbuf[i];
      }
      else {
        GmfGetLin(inm,keyword,dbuf);
        ia = (k-1)*sol->size[1] + 1;
        for (i=0; i<sol->size[1]; i++)
          sol->valp0[ia+i] = dbuf[i];
      }
    }
  }

  if ( GmfStatKwd(inm,GmfIterations) ) {
    GmfGotoKwd(inm,GmfIterations);
    GmfGetLin(inm,GmfIterations,&sol->iter);
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
  return(1);
}


int saveSol(pSol sol,char *filename) {
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  int          k,i,ia,inm,keyword;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".meshb");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".mesh");
    if ( ptr ) {
      *ptr = '\0';
      strcat(data,".sol");
    }
    else
      strcat(data,".solb");
  }

  if (!(inm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* write sol */
  if ( sol->np ) {
    GmfSetKwd(inm,GmfSolAtVertices,sol->np,sol->type[0],sol->typtab[0]);
    for (k=1; k<=sol->np; k++) {
      ia = (k-1)*sol->size[0] + 1;
      if ( sol->ver == GmfFloat ) {
        for (i=0; i<sol->size[0]; i++)
          fbuf[i] = sol->valp1[ia+i];
        GmfSetLin(inm,GmfSolAtVertices,fbuf);
      }
      else {
        for (i=0; i<sol->size[0]; i++)
          dbuf[i] = sol->valp1[ia+i];
        GmfSetLin(inm,GmfSolAtVertices,dbuf);
      }
    }
  }

  if ( sol->ne ) {
    keyword = sol->dim == 2 ? GmfSolAtTriangles : GmfSolAtTetrahedra;
    GmfSetKwd(inm,keyword,sol->ne,sol->type[1],sol->typtab[1]);

    for (k=1; k<=sol->ne; k++) {
      ia = (k-1)*sol->size[1] + 1;
      if ( sol->ver == GmfFloat ) {
        for (i=0; i<sol->size[1]; i++)
          fbuf[i] = sol->valp0[ia+i];
        GmfSetLin(inm,keyword,fbuf);
      }
      else {
        for (i=0; i<sol->size[1]; i++)
          dbuf[i] = sol->valp0[ia+i];
        GmfSetLin(inm,keyword,dbuf);
      }
    }
  }

  if ( sol->iter ) {
    GmfSetKwd(inm,GmfIterations);
    GmfSetLin(inm,GmfIterations,sol->iter);
  }
	if ( sol->time > 0.0 ) {
    GmfSetKwd(inm,GmfTime);
    if ( sol->ver == GmfFloat ) {
			fbuf[0] = sol->time;
      GmfSetLin(inm,GmfTime,fbuf[0]);
		}
		else
			GmfSetLin(inm,GmfTime,sol->time);
	}
  GmfCloseMesh(inm);
  return(1);
}

