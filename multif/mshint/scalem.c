#include "mshint.h"

extern Param  par;
extern unsigned char inxt[3];
unsigned char iedg[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
 

int scaleMesh(pMesh mesh,pSol sol) {
  pPoint    ppt;
	pTetra    ptt;
	pTria     pt;
  Info     *info;
  double    dd,d2,ux,uy,uz;
	int       i,j,k,kk,iadr,ia,ib;

  /* normalize coordinates */
  info = &mesh->info;
  dd = (double)PRECI / info->delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++)
      ppt->c[i] = dd * (ppt->c[i] - info->min[i]);
  }

  if ( mesh->dim == 2 ) {
		for (k=1; k<=mesh->nt; k++) {
			pt = &mesh->tria[k];
			pt->g[0] = dd * (pt->g[0] - info->min[0]);
			pt->g[1] = dd * (pt->g[1] - info->min[1]);
			for (i=0; i<3; i++) {
				ia = pt->v[i];
				ib = pt->v[inxt[i]];
				ux = mesh->point[ib].c[0] - mesh->point[ia].c[0];
				uy = mesh->point[ib].c[1] - mesh->point[ia].c[1];
				d2 = sqrt(ux*ux + uy*uy);
				pt->h = I_MAX(d2,pt->h);
			}
		}
	}
	else {
		for (k=1; k<=mesh->ne; k++) {
			ptt = &mesh->tetra[k];
			ptt->g[0] = dd * (ptt->g[0] - info->min[0]);
			ptt->g[1] = dd * (ptt->g[1] - info->min[1]);
			ptt->g[2] = dd * (ptt->g[2] - info->min[2]);
			for (i=0; i<6; i++) {
				ia = ptt->v[iedg[i][0]];
				ib = ptt->v[iedg[i][1]];
				ux = mesh->point[ib].c[0] - mesh->point[ia].c[0];
				uy = mesh->point[ib].c[1] - mesh->point[ia].c[1];
				uz = mesh->point[ib].c[2] - mesh->point[ia].c[2];
				d2 = sqrt(ux*ux + uy*uy + uz*uz);
				ptt->h = I_MAX(d2,ptt->h);
			}
		}
	}
	 
  /* normalize metric */
  if ( !sol )  return(1);

  d2 = 1.0 / (dd*dd);
  for (k=1; k<=sol->np; k++) {
    iadr = (k-1) * sol->size[0] + 1;
    for (i=0; i<sol->type[0]; i++) {
      switch(sol->typtab[0][i]){
      case GmfSca:
        sol->valp1[iadr] *= dd;
        iadr++;
        break;
      case GmfVec:
        for (j=0; j<mesh->dim; j++)
          sol->valp1[iadr+j] *= dd;
        iadr += j;
        break;
      case GmfSymMat:
        kk = mesh->dim * (mesh->dim+1) / 2;
        for (j=0; j<kk; j++)  
          sol->valp1[iadr+j] *= d2;
        iadr += kk;
        break;
      default:
        fprintf(stdout,"  ## SPECIFIC DATA IGNORED %d\n",sol->typtab[0][i]);
        return(0);
      }
    }
  }

  for (k=1; k<=sol->ne; k++) {
    iadr = (k-1) * sol->size[1] + 1;
    for (i=0; i<sol->type[1]; i++) {
      switch(sol->typtab[1][i]){
      case GmfSca:
        sol->valp0[iadr] *= dd;
        iadr++;
        break;
      case GmfVec:
        for (j=0; j<mesh->dim; j++)
          sol->valp0[iadr+j] *= dd;
        iadr += j;
        break;
      case GmfSymMat:
        kk = mesh->dim * (mesh->dim+1) / 2;
        for (j=0; j<kk; j++)  
          sol->valp0[iadr+j] *= d2;
        iadr += kk;
        break;
      default:
        fprintf(stdout,"  ## SPECIFIC DATA IGNORED %d\n",sol->typtab[1][i]);
        return(0);
      }
    }
  }
  return(1);
}


int unscaleMesh(pMesh mesh,pSol sol) {
  pPoint     ppt;
  Info      *info;
  double     dd,d2;
  int        i,j,k,kk,iadr;

  info = &mesh->info;

  /* de-normalize coordinates */
  dd = info->delta / (double)PRECI;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++)
      ppt->c[i] = ppt->c[i] * dd + info->min[i];
  }

  /* de-normalize metric */
  d2 = 1.0 / (dd*dd);
  for (k=1; k<=sol->np; k++) {
    iadr = (k-1) * sol->size[0] + 1;
    for (i=0; i<sol->type[0]; i++) {
      switch(sol->typtab[0][i]) {
      case GmfSca:
        sol->valp1[iadr] *= dd;
        iadr++;
        break;
      case GmfVec:
        for (j=0; j<mesh->dim; j++)
          sol->valp1[iadr+j] *= dd;
        iadr += mesh->dim;
        break;
      case GmfSymMat:
        kk = mesh->dim * (mesh->dim+1) / 2;
        for (j=0; j<kk; j++)  
          sol->valp1[iadr+j] *= d2;
        iadr += kk;
        break;
      }
    }
  }

  for (k=1; k<=sol->ne; k++) {
    iadr = (k-1) * sol->size[1] + 1;
    for (i=0; i<sol->type[1]; i++) {
      switch(sol->typtab[1][i]) {
      case GmfSca:
        sol->valp0[iadr] *= dd;
        iadr++;
        break;
      case GmfVec:
        for (j=0; j<mesh->dim; j++)
          sol->valp0[iadr+j] *= dd;
        iadr += mesh->dim;
        break;
      case GmfSymMat:
        kk = mesh->dim * (mesh->dim+1) / 2;
        for (j=0; j<kk; j++)  
          sol->valp0[iadr+j] *= d2;
        iadr += kk;
        break;
      }
    }
  }

  return(1);
}



