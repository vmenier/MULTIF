#include "mshint.h"

#define BUCKSIZ    32

extern Param par;


int mshin1(pMesh mesh1,pSol sol1,pMesh mesh2,pSol sol2) {
  pBucket  bucket;
  pPoint   ppt;
  pTria    pt;
  pTetra   ptt;
  double   cb[4],*sp,*sa,*se;
  int      k,iadr,iel,ret,fail;

  /* interpolation */
  bucket = newBucket(mesh1,BUCKSIZ);
  if ( !bucket )  return(0);

  fail = 0;
  if ( sol1->np ) {
    if ( par.imprim )  fprintf(stdout,"  %%%% Nodal interpolation (%d)\n",mesh2->np);
    for (k=1; k<=mesh2->np; k++) {
      ppt  = &mesh2->point[k];
      iel  = buckin(mesh1,bucket,ppt->c);
      iel  = locelt(mesh1,iel,ppt->c,cb);
      iadr = (k-1)*sol2->size[0] + 1;
      sp   = &sol2->valp1[iadr];
      if ( iel ) {
        if ( mesh1->dim == 2 )
          ret = intpp1(sol1,mesh1->tria[iel].v,sp,k,cb);
        else
          ret = intpp1(sol1,mesh1->tetra[iel].v,sp,k,cb);
        if ( !ret )  fail++;
      }
      else {
        iel = closept(mesh1,ppt->c);
        if ( iel ) {
          iadr = (iel-1)*sol1->size[0] + 1;
          sa   = &sol1->valp1[iadr];
          memcpy(sp,sa,sol2->size[0]*sizeof(double));
        }
        else  fail++;
      }
    }
  }

  if ( sol1->ne ) {
    if ( par.imprim )  fprintf(stdout,"  %%%% Element interpolation (%d)\n",sol2->ne);

    if ( mesh2->dim == 2 ) {
      for (k=1; k<=mesh2->nt; k++) {
        pt  = &mesh2->tria[k];
        iel = buckin(mesh1,bucket,pt->g);
        iel = locelt(mesh1,iel,pt->g,cb);
        iadr = (k-1)*sol2->size[1] + 1;
        se   = &sol2->valp0[iadr];
        if ( iel ) {
          ret = intpp0(mesh1,sol1,se,iel,pt->g,par.ray*pt->h);
          if ( !ret )  fail++;
        }
        else {
          if ( par.ddebug )  fprintf(stdout,"  ## No element found. Wrong solution\n");
          fail++;
        }
      }
    }
    else {
      for (k=1; k<=mesh2->ne; k++) {
        ptt = &mesh2->tetra[k];
        iel = buckin(mesh1,bucket,ptt->g);
        iel = locelt(mesh1,iel,ptt->g,cb);

        iadr = (k-1)*sol2->size[1] + 1;
        se   = &sol2->valp0[iadr];
        if ( iel ) {
          ret = intpp0(mesh1,sol1,se,iel,ptt->g,par.ray*ptt->h);
          if ( !ret )  fail++;
        }
        else {
          if ( par.ddebug )  fprintf(stdout,"  ## No element found. Wrong solution\n");
          fail++;
        }
      }
    }
  }
  if ( fail )  fprintf(stdout,"  ## WARNING: %d nodes failed. Wrong solution.\n",fail);

  free(bucket->head);
  free(bucket->link);
  free(bucket);
  return(1);
}
