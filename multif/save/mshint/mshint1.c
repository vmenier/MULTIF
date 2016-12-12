#include "mshint.h"


int mshint1_3d(MIst *mist) {
  Bucket  *bck;
  pPoint   ppt;
  pTetra   pt;
  double  *u,*ua,*ub,*uc,*ud,cb[4];
  int      i,k,iel,ip,nb,np;

  /* interpolation */
  bck = bucket_3d(&mist->msrc,BUCKSIZ);
  if ( !bck )  return(0);

  if ( mist->info.verb != '0' ) {
    fprintf(stdout,"    Nodal interpolation: ");
    fflush(stdout);
  }
 
  nb = np = 0;
  for (k=1; k<=mist->mtgt.np; k++) {
    ppt  = &mist->mtgt.point[k];
    /* get seed element */
    iel  = buckin_3d(&mist->msrc,bck,ppt->c);
    if ( iel ) {
      iel  = locelt_3d(&mist->msrc,iel,ppt->c,cb);
      if ( iel > 0 ) {
        pt = &mist->msrc.tetra[iel];
        u  = &mist->stgt.u[(k-1)*mist->stgt.size[0]];
        ua = &mist->ssrc.u[(pt->v[0]-1)*mist->ssrc.size[0]];
        ub = &mist->ssrc.u[(pt->v[1]-1)*mist->ssrc.size[0]];
        uc = &mist->ssrc.u[(pt->v[2]-1)*mist->ssrc.size[0]];
        ud = &mist->ssrc.u[(pt->v[3]-1)*mist->ssrc.size[0]];

        for (i=0; i<mist->stgt.size[0]; i++)
          u[i] = cb[0]*ua[i] + cb[1]*ub[i] + cb[2]*uc[i] + cb[3]*uc[i];
      }
    }
    if ( iel < 1 ) {
      /* exhaustive search */
      ip = closept_3d(&mist->msrc,ppt->c);
      if ( ip ) {
        ua = &mist->ssrc.u[(ip-1)*mist->ssrc.size[0]];
        u  = &mist->stgt.u[(k-1)*mist->stgt.size[0]];
        memcpy(u,ua,mist->stgt.size[0]*sizeof(double));
        np++;
      }
      else
        nb++;
    }
  }
  if ( mist->info.verb != '0' ) {
    fprintf(stdout,"%d, %d prox.\n",mist->mtgt.np,np);
    if ( nb > 0 )  fprintf(stdout," # Warning: %d nodes failed.\n",nb);
  }

  free(bck->cell);
  free(bck);

  return(nb == 0);
}


int mshint1_2d(MIst *mist) {
  Bucket  *bck;
  pPoint   ppt;
  pTria    pt;
  double  *u,*ua,*ub,*uc,cb[3];
  int      k,i,iel,ip,nb;

  /* interpolation */
  bck = bucket_2d(&mist->msrc,BUCKSIZ);
  if ( !bck )  return(0);

  if ( mist->info.verb != '0' ) {
    fprintf(stdout,"    Nodal interpolation: ");
    fflush(stdout);
  }
  puts("debut\n----------------------\n");
  nb = 0;
  for (k=1; k<=mist->mtgt.np; k++) {
    ppt  = &mist->mtgt.point[k];
    /* get seed element */
    iel  = buckin_2d(&mist->msrc,bck,ppt->c);
    if ( iel > 0 ) {
      iel  = locelt_2d(&mist->msrc,iel,ppt->c,cb);
      if ( iel > 0 ) {
        pt = &mist->msrc.tria[iel];
        u  = &mist->stgt.u[(k-1)*mist->stgt.size[0]];
        ua = &mist->ssrc.u[(pt->v[0]-1)*mist->ssrc.size[0]];
        ub = &mist->ssrc.u[(pt->v[1]-1)*mist->ssrc.size[0]];
        uc = &mist->ssrc.u[(pt->v[2]-1)*mist->ssrc.size[0]];
        for (i=0; i<mist->stgt.size[0]; i++)
          u[i] = cb[0]*ua[i] + cb[1]*ub[i] + cb[2]*uc[i];
      }
    }
    if ( iel < 1 ) {
      ip = closept_2d(&mist->msrc,ppt->c);
      if ( ip > 0 ) {
        ua = &mist->ssrc.u[(ip-1)*mist->ssrc.size[0]];
        u  = &mist->stgt.u[(k-1)*mist->stgt.size[0]];
        memcpy(u,ua,mist->stgt.size[0]*sizeof(double));
      }
      else
        nb++;
    }
  }
  if ( mist->info.verb != '0' ) {
    fprintf(stdout,"%d\n",mist->mtgt.np);
    if ( nb > 0 )  fprintf(stdout," # Warning: %d nodes failed.\n",nb);
  }

  free(bck->cell);
  free(bck);

  return(nb == 0);
}
