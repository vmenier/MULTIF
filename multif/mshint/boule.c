#include "mshint.h"

unsigned char inxt[3] = {1,2,0};
unsigned char iprv[3] = {2,0,1};


/* find all vertices connected to P
   in:  start : tetrahedron containing p 
        ip    : index of p in start
        list  : dynamic list structure (allocated outside)
   out: list  : list of tets */
int boulep_3d(pMesh mesh,int start,int ip,int *list) {
  pTetra  pt,pt1;
  int    *adja,vois[4],adj,i,j,indp,iel,iadr,base,ilist,nump;

  pt   = &mesh->tetra[start];
  base = ++mesh->mark;
  pt->mark = base;
  list[0] = 4*start + ip;
  ilist   = 1;

  iadr = (start-1)*4 + 1;
  adja = &mesh->adja[iadr];
  vois[0]  = adja[0] >> 2;
  vois[1]  = adja[1] >> 2;
  vois[2]  = adja[2] >> 2;
  vois[3]  = adja[3] >> 2;

  /* store 3 neighbors sharing P */
  for (i=0; i<4; i++) {
    if ( i == ip )  continue;
    adj = vois[i];
		if ( !adj )  continue;
    pt1 = &mesh->tetra[adj];
		if ( pt1->mark != base ) {
      /* not stored yet */
      pt1->mark = base;
      for (j=0; j<4; j++) {
        if ( pt1->v[j] == pt->v[ip] ) {
          ilist++;
          list[ilist] = 4*adj + j;
	        break;
	      }
	    }
    }
  }
  if ( ilist < 1 )  return(ilist);

  /* explore list of neighbors */
  indp = 2;
	nump = pt->v[ip];
  do {
    iel  = list[indp] >> 2;
    pt   = &mesh->tetra[iel];
    iadr = (iel-1)*4 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0] >> 2;
    vois[1]  = adja[1] >> 2;
    vois[2]  = adja[2] >> 2;
    vois[3]  = adja[3] >> 2;

    for (i=0; i<4; i++) {
      if ( pt->v[i] == nump )  continue;
      adj = vois[i];
			if ( !adj )  continue;
      pt1 = &mesh->tetra[adj];
      if ( pt1->mark != base ) {
        pt1->mark = base;
        for (j=0; j<4; j++) {
          if ( pt1->v[j] == nump ) {
            list[ilist] = 4*adj + j;
            ilist++;
	          break;
          }
        }
      }
    }  
    /* overflow */
    if ( ilist > LONMAX-3 )  return(-ilist);
  }
  while ( ++indp <= ilist );

  return(ilist);
}


/* store neighboring vertices, return < 0 if overflow */
int boulep_2d(pMesh mesh,int start,int ip,int *list) {
  pTria    pt;
  int     *adja,iadr,base,ilist,iel;
	char     ii,i1,voy;

  pt = &mesh->tria[start];
  base = ++mesh->mark;
  pt->mark = base;
  list[0] = 3*start + ip;
  ilist   = 1;

  iadr = 3*(start-1) + 1;
  adja = &mesh->adja[iadr];

  /* store neighbors */
  i1  = inxt[ip];
  iel = adja[(int)i1] / 3;
  while ( iel && (iel != start) ) {
	  pt  = &mesh->tria[iel]; 
    voy = adja[(int)i1] % 3;
    i1  = iprv[(int)voy];
		ii  = inxt[(int)voy];
	  if ( ilist > LONMAX-2 )  return(-ilist);
    list[ilist] = 3*iel + ii;
	  ++ilist;
    iadr = 3*(iel-1) + 1;
    adja = &mesh->adja[iadr];
    iel  = adja[(int)i1] / 3;
  }

  /* reverse loop */
  if ( iel != start ) {
    i1   = iprv[ip];
    iadr = (start-1)*3 + 1;
    adja = &mesh->adja[(int)iadr];
    iel =  adja[(int)i1] / 3;
    while ( iel && (iel != start) ) {
      pt  = &mesh->tria[iel];
      voy = adja[(int)i1] % 3;
			i1  = inxt[(int)voy];
			ii  = iprv[(int)voy];
	    if ( ilist > LONMAX-2 )  return(-ilist);
      list[ilist] = 3*iel + ii;
   	  ++ilist;
      iadr = 3*(iel-1) + 1;
      adja = &mesh->adja[iadr];
      iel  = adja[(int)i1] / 3;
    }
  }

  return(ilist);
}



