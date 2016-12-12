#include "mshint.h"


/* return closest point of p on mesh1 */
int closept_3d(Mesh *mesh,double *c) {
  pPoint   p1;
	double   ux,uy,uz,dd,dmin;
	int      i,ipmin;

  dmin  = MI_TGV;
  ipmin = 0;
  for (i=1; i<=mesh->np; i++) {
    p1 = &mesh->point[i];
    ux = p1->c[0] - c[0];
    uy = p1->c[1] - c[1];
    uz = p1->c[2] - c[2];
    dd = ux*ux + uy*uy + uz*uz;
    if ( dd < dmin) {
      ipmin = i;
      dmin  = dd;
    }
  }

	return(ipmin);
}


int closept_2d(pMesh mesh,double *c) {
	pPoint   p1;
	double   ux,uy,dd,dmin;
	int      i,ipmin;

  dmin  = MI_TGV;
  ipmin = 0;
  for (i=1; i<=mesh->np; i++) {
    p1 = &mesh->point[i];
    ux = p1->c[0] - c[0];
    uy = p1->c[1] - c[1];
    dd = ux*ux + uy*uy;
    if ( dd < dmin) {
      ipmin = i;
      dmin  = dd;
    }
  }

	return(ipmin);
}


/* find element containing c, starting from nsdep, return baryc. coord */
int locelt_3d(pMesh mesh,int nsd,double *c,double *cb) {
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
  double   eps,vto,vol1,vol2,vol3,vol4,dd; 
  int      i,nsf,nsp;

  nsf = nsd;
  nsp = nsd;
  ++mesh->mark;
  while ( nsf > 0 ) {
    pt = &mesh->tetra[nsf];
    if ( pt->mark == mesh->mark )  return(nsp);
    pt->mark = mesh->mark;
    
    /* measure of element */
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];

    /* barycentric and volume */
    bx  = p1->c[0] - p0->c[0];
    by  = p1->c[1] - p0->c[1];
    bz  = p1->c[2] - p0->c[2];
    cx  = p2->c[0] - p0->c[0];
    cy  = p2->c[1] - p0->c[1];
    cz  = p2->c[2] - p0->c[2];
    dx  = p3->c[0] - p0->c[0];
    dy  = p3->c[1] - p0->c[1];
    dz  = p3->c[2] - p0->c[2];

    /* test volume */
    vx  = cy*dz - cz*dy;
    vy  = cz*dx - cx*dz;
    vz  = cx*dy - cy*dx;
    vto = bx*vx + by*vy + bz*vz;
		eps = MI_EPS*vto;

    /* barycentric */
    apx = c[0] - p0->c[0];
    apy = c[1] - p0->c[1];
    apz = c[2] - p0->c[2];

    /* p in 2 */
    vol2  = apx*vx + apy*vy + apz*vz;
    if ( vol2 < eps ) {
      nsp = nsf;
      nsf = pt->adj[1] / 4;
      if ( !nsf ) {
        cb[1] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 3 */
    vx  = by*apz - bz*apy;
    vy  = bz*apx - bx*apz;
    vz  = bx*apy - by*apx;
    vol3 = dx*vx + dy*vy + dz*vz;
    if ( vol3 < eps ) {
			nsp = nsf;
      nsf = pt->adj[2] / 4;
      if ( !nsf ) {
        cb[2] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 4 */
    vol4 = -cx*vx - cy*vy - cz*vz;
    if ( vol4 < eps ) {
			nsp = nsf;
      nsf = pt->adj[3] / 4;
      if ( !nsf ) {
        cb[3] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 1 */
    vol1 = vto - vol2 - vol3 - vol4;
    if ( vol1 < eps ) {
      nsp = nsf;
      nsf = pt->adj[0] / 4;
      if ( !nsf ) {
        cb[0] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    dd = fabs(vol1+vol2+vol3+vol4);
    if ( dd > MI_EPSD ) {
	    dd = 1.0 / dd;
      cb[0] = vol1 * dd;
      cb[1] = vol2 * dd;
      cb[2] = vol3 * dd;
      cb[3] = vol4 * dd;
    }
    if ( (cb[0]>=0.0) && (cb[1]>=0.0) && (cb[2]>=0.0) && (cb[3]>=0.0) )
      return(nsf);
    else
      return(-nsf);
  }

  return(nsp);
}


/* find simplex containing c, starting from nsd, return baryc. coord  */
int locelt_2d(Mesh *mesh,int nsd,double *c,double *cb) {
  pTria     pt;
  pPoint    p0,p1,p2;
  double    ax,ay,bx,by,cx,cy,eps;
  double    aire1,aire2,aire3,dd; 
  int       i,sgn,nsf,nsp;

  nsf  = nsd;
  nsp  = 0;
  mesh->mark = ++mesh->mark;
  while ( nsf > 0 ) {
    pt = &mesh->tria[nsf];
    if ( pt->mark == mesh->mark )  return(nsp);
    pt->mark = mesh->mark;

    /* area of triangle */
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    bx = p2->c[0] - p0->c[0];
    by = p2->c[1] - p0->c[1];
    dd = ax*by - ay*bx;
    sgn = dd > 0.0 ? 1 : -1;
    eps = sgn == 1 ? -MI_EPS*dd : MI_EPS*dd;
    /* barycentric */
    bx = p1->c[0] - c[0];
    by = p1->c[1] - c[1];
    cx = p2->c[0] - c[0];
    cy = p2->c[1] - c[1];

    /* p in half-plane lambda_0 > 0 */
    aire1 = sgn*(bx*cy - by*cx);
    if ( aire1 < eps ) {
      nsp = nsf;
      nsf = pt->adj[0] / 3;
      if ( !nsf ) {
        cb[0] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    ax = p0->c[0] - c[0];
    ay = p0->c[1] - c[1];
    aire2 = sgn*(cx*ay - cy*ax);
    if ( aire2 < eps ) {
      nsp = nsf;
      nsf = pt->adj[1] / 3;
      if ( !nsf ) {
        cb[1] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    aire3 = sgn*dd - aire1 - aire2;
    if ( aire3 < eps ) {
      nsp = nsf;
      nsf = pt->adj[2] / 3;
      if ( !nsf ) {
        cb[2] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    dd = fabs(aire1+aire2+aire3);
    if ( dd > MI_EPSD ) {
      dd = 1.0 / dd;
      cb[0] = aire1 * dd;
      cb[1] = aire2 * dd;
      cb[2] = aire3 * dd;
    }
    if ( (cb[0]>=0.0) && (cb[1]>=0.0) && (cb[2]>=0.0) )
      return(nsf);
    else
      return(-nsf);
  }

  return(nsp);
}




