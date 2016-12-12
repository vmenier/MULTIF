#include "mshint.h"


/* P1 interpolation */
int intpp1_3d(pSol sol1,int *v,double *sp,int ip,double *cb) {
  double  *sa,*sb,*sc,*sd;
  int      i,iada,iadb,iadc,iadd;

  iada = (v[0]-1)*sol1->size[0] + 1;
  iadb = (v[1]-1)*sol1->size[0] + 1;
  iadc = (v[2]-1)*sol1->size[0] + 1;
  iadd = (v[3]-1)*sol1->size[0] + 1;

  sa = &sol1->valp1[iada];
  sb = &sol1->valp1[iadb];
  sc = &sol1->valp1[iadc];
  sd = &sol1->valp1[iadd];
  for (i=0; i<sol1->size[0]; i++)
    sp[i] = cb[0]*sa[i] + cb[1]*sb[i] + cb[2]*sc[i] + cb[3]*sd[i];

  return(1);
}


/* P1 interpolation */
int intpp1_2d(pSol sol1,int *v,double *sp,int ip,double *cb) {
  double  *sa,*sb,*sc;
  int      i,iada,iadb,iadc;

  iada = (v[0]-1)*sol1->size[0] + 1;
  iadb = (v[1]-1)*sol1->size[0] + 1;
  iadc = (v[2]-1)*sol1->size[0] + 1;

  sa = &sol1->valp1[iada];
  sb = &sol1->valp1[iadb];
  sc = &sol1->valp1[iadc];

  for (i=0; i<sol1->size[0]; i++)
    sp[i] = cb[0]*sa[i] + cb[1]*sb[i] + cb[2]*sc[i];

  return(1); 
}


/* return closest point of p on mesh1 */
int closept_3d(pMesh mesh,double *p) {
  pPoint   p1;
	double   ux,uy,uz,dd,dmin;
	int      i,ipmin;

  dmin  = 1.e30;
  ipmin = 0;
  for (i=1; i<=mesh->np; i++) {
    p1 = &mesh->point[i];
    ux = p1->c[0] - p[0];
    uy = p1->c[1] - p[1];
    uz = p1->c[2] - p[2];
    dd = ux*ux + uy*uy + uz*uz;
    if ( dd < dmin) {
      ipmin = i;
      dmin  = dd;
    }
  }
	return(ipmin);
}


int closept_2d(pMesh mesh,double *p) {
	pPoint   p1;
	double   ux,uy,dd,dmin;
	int      i,ipmin;

  dmin  = 1.e30;
  ipmin = 0;
  for (i=1; i<=mesh->np; i++) {
    p1 = &mesh->point[i];
    ux = p1->c[0] - p[0];
    uy = p1->c[1] - p[1];
    dd = ux*ux + uy*uy;
    if ( dd < dmin) {
      ipmin = i;
      dmin  = dd;
    }
  }
	return(ipmin);
}


int intpp0_3d(pMesh mesh1,pSol sol1,double *se,int iel,double *g,double r) {
	pTetra  pt;
	double *sk,dd,d1,ux,uy,uz;
	int     list[LONMAX+3],i,j,k,ne,nel,iadr;

	nel = 0;
  for (i=0; i<4; i++) {
		ne = boulep(mesh1,iel,i,list);
		nel += abs(ne);
		if ( ne < 0 )  break;
  }
	if ( nel < 0 )  return(0);

	/* compute barycenter */
	d1 = 0.0;
	for (i=0; i<=nel; i++) {
		k  = list[i] >> 2;
		pt = &mesh1->tetra[k];
		ux = pt->g[0] - g[0];
		uy = pt->g[1] - g[1];
		uz = pt->g[2] - g[2];
		dd = sqrt(ux*ux + uy*uy + uz*uz);
		if ( i > 0 && dd >= r )  continue;
		dd  = dd > 1.0e-15 ? 1.0 / dd : 1.0;
	  d1 += dd;
	  iadr = (k-1)*sol1->size[1] + 1;
	  sk   = &sol1->valp0[iadr];
	  for (j=0; j<sol1->size[1]; j++) {
	    se[j] += dd * sk[j];
    }
	}
  for (j=0; j<sol1->size[1]; j++)
		se[j] /= d1;

  return(1);
}


int intpp0_2d(pMesh mesh1,pSol sol1,double *se,int iel,double *g,double r) {
	pTria   pt;
	double *sk,dd,d1,ux,uy;
	int     list[LONMAX+3],i,j,k,ne,nel,iadr;

	nel = 0;
  for (i=0; i<3; i++) {
		ne = boulep(mesh1,iel,i,list);
		nel += abs(ne);
		if ( ne < 0 )  break;
  }
	if ( nel < 0 )  return(0);
	
	/* compute barycenter */
	d1 = 0.0;
	for (i=0; i<=nel; i++) {
		k  = list[i] / 3;
		pt = &mesh1->tria[k];
		ux = pt->g[0] - g[0];
		uy = pt->g[1] - g[1];
		dd = sqrt(ux*ux + uy*uy);
		if ( i > 0 && dd > r )  continue;
		dd  = dd > 1.0e-15 ? 1.0 / dd : 1.0;
	  d1 += dd;
	  iadr = (k-1)*sol1->size[1] + 1;
	  sk   = &sol1->valp0[iadr];
	  for (j=0; j<sol1->size[1]; j++)
			se[j] += dd * sk[j];
	}
  for (j=0; j<sol1->size[1]; j++)
		se[j] /= d1;

	return(1);
}



