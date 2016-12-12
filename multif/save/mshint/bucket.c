#include "mshint.h"


/* create and manage bucket structure */
Bucket *bucket_3d(Mesh *mesh,int nmax) {
  pPoint   ppt;
  Bucket  *bck;
  double   dx,dy,dz;
  int      k,c,i,j,l;

  /* memory alloc */
  bck = (Bucket*)calloc(1,sizeof(Bucket));
  assert(bck);
  bck->size = nmax;
  bck->cell = (int*)calloc(nmax*nmax*nmax,sizeof(int));
  assert(bck->cell);

  /* insert vertices */
  dx = (nmax-1) / (mesh->max[0] - mesh->min[0]);
  dy = (nmax-1) / (mesh->max[1] - mesh->min[1]);
  dz = (nmax-1) / (mesh->max[2] - mesh->min[2]);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    i = (int)(dx * (ppt->c[0] - mesh->min[0]));
    j = (int)(dy * (ppt->c[1] - mesh->min[1]));
    l = (int)(dz * (ppt->c[2] - mesh->min[2]));
    c = (l*bck->size +j)*bck->size + i;

    if ( !bck->cell[c] )  bck->cell[c] = ppt->s;
  }

  return(bck);
}


/* create and manage bucket structure */
Bucket *bucket_2d(Mesh *mesh,int nmax) {
  pPoint   ppt;
  Bucket  *bck;
  double   dx,dy;
  int      k,c,i,j;

  /* memory alloc */
  bck = (Bucket*)calloc(1,sizeof(Bucket));
  assert(bck);
  bck->size = nmax;
  bck->cell = (int*)calloc(nmax*nmax,sizeof(int));
  assert(bck->cell);

  /* insert vertices */
  dx = (bck->size-1) / (mesh->max[0] - mesh->min[0]);
  dy = (bck->size-1) / (mesh->max[1] - mesh->min[1]);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    i = (int)(dx * (ppt->c[0] - mesh->min[0]));
    j = (int)(dy * (ppt->c[1] - mesh->min[1]));
    c = j*bck->size + i;

    if ( !bck->cell[c] )  bck->cell[c] = ppt->s;
  }

  return(bck);
}


/* return simplex in bucket cell close to c */
int buckin_3d(Mesh *mesh,Bucket *bck,double *c) {
  double    dx,dy,dz;
  int       k,i,j,ii,jj,ll,b,d;

  dx = (bck->size-1) / (mesh->max[0] - mesh->min[0]);
  dy = (bck->size-1) / (mesh->max[1] - mesh->min[1]);
  dz = (bck->size-1) / (mesh->max[2] - mesh->min[2]);
  ii = (int)(dx * (c[0] - mesh->min[0]));
  jj = (int)(dy * (c[1] - mesh->min[1]));
  ll = (int)(dz * (c[2] - mesh->min[2]));
  b = (ll*bck->size +jj)*bck->size + ii;

  /* check current cell */
  if ( bck->cell[b] )
    return(bck->cell[b]);
  /* explore neighboring cells */
  else {
    d = 1;
    do {
      for (k=MI_MAX(0,ll-d); k<MI_MIN(bck->size,ll+d); k++) {
        for (j=MI_MAX(0,jj-d); j<MI_MIN(bck->size,jj+d); j++) {
          for (i=MI_MAX(0,ii-d); i<MI_MIN(bck->size,ii+d); i++) {
            b = (k*bck->size +j)*bck->size + i;
            if ( bck->cell[b] )  return(bck->cell[b]);
          }
        }
      }
    }
    while ( ++d < bck->size/8 );
  }

  return(0);
}


/* return simplex in bucket cell close to c */
int buckin_2d(Mesh *mesh,Bucket *bck,double *c) {
  double    dx,dy;
  int       k,i,j,ii,jj,b,d;

  dx = (bck->size-1) / (mesh->max[0] - mesh->min[0]);
  dy = (bck->size-1) / (mesh->max[1] - mesh->min[1]);
  ii = (int)(dx * (c[0] - mesh->min[0]));
  jj = (int)(dy * (c[1] - mesh->min[1]));
  b  = jj*bck->size + ii;

  /* check current cell */
  if ( bck->cell[b] )
    return(bck->cell[b]);
  /* explore neighboring cells */
  else {
    d = 1;
    do {
      for (j=MI_MAX(0,jj-d); j<MI_MIN(bck->size,jj+d); j++) {
        for (i=MI_MAX(0,ii-d); i<MI_MIN(bck->size,ii+d); i++) {
          b = j*bck->size + i;
          if ( bck->cell[b] )  return(bck->cell[b]);
        }
      }
    }
    while ( ++d < bck->size/8 );
  }

  return(0);
}



