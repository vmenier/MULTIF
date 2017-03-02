#ifndef _INT_H
#define _INT_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "libmesh6.h"

#include "mesh.h"
#include "fproto.h"


#define I_VER   "2.1a"
#define I_REL   "Mar 16, 2007"
#define I_CPY   "Copyright (c) LJLL, 2007"
#define I_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

#define I_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define I_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define LONMAX   8192
#define EPS1     1.e-20
#define EPST    -1.e-2
#define EPSR     1.e+2
#define PRECI    1.0


typedef struct {
  double         c[3];
  int            s;
  unsigned char  b;
} Point;
typedef Point * pPoint;

typedef struct {
	double  g[2],h;
  int     v[3];
  int     mark;
} Tria;
typedef Tria * pTria;

typedef struct {
	double  g[3],h;
  int     v[4];
  int     mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  double   delta,min[3],max[3];
} Info;

typedef struct {
  int      np,nt,ne,ver,dim;
  int     *adja,mark;
  char    *name,bin;

  pPoint   point;
  pTria    tria;
  pTetra   tetra;
  Info     info;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
	int         np,ne,ver,dim,iter;
	int         type[2],size[2],typtab[2][GmfMaxTyp];
  double     *valp1,*valp0;
  float       time;
  char       *name;
} Sol;
typedef Sol * pSol;

typedef struct {
  int     size;
  int    *head;
  int    *link;
} Bucket;
typedef Bucket * pBucket;

typedef struct {
  double   dt,ray;
  char     imprim,ddebug,option; 
  mytime   ctim[TIMEMAX];
} Param;
	
/* prototypes */
int loadMesh(pMesh mesh,char *filename);
int loadSol(pSol sol,char *filename);
int saveSol(pSol sol,char *filename);
int mshin1(pMesh mesh1,pSol sol1,pMesh mesh2,pSol sol2);
int scaleMesh(pMesh mesh,pSol sol);
int unscaleMesh(pMesh mesh,pSol sol);


int locateTetra(pMesh mesh,int nsdep,int base,double *p,double *cb);
int locateTria(pMesh mesh,int nsdep,int base,double *p,double *cb);

pBucket newBucket_2d(pMesh ,int );
pBucket newBucket_3d(pMesh ,int );
int     buckin_2d(pMesh ,pBucket ,double *);
int     buckin_3d(pMesh ,pBucket ,double *);
int     intpp0_3d(pMesh ,pSol ,double *,int ,double *,double );
int     intpp0_2d(pMesh ,pSol ,double *,int ,double *,double );
int     intpp1_2d(pSol ,int *,double *,int ,double *);
int     intpp1_3d(pSol ,int *,double *,int ,double *);
int     locelt_2d(pMesh ,int ,double *,double *);
int     locelt_3d(pMesh ,int ,double *,double *);
int     closept_2d(pMesh ,double *);
int     closept_3d(pMesh ,double *);
int     hashelt_3d(pMesh );
int     hashelt_2d(pMesh );
int     boulep_2d(pMesh ,int ,int ,int *);
int     boulep_3d(pMesh ,int ,int ,int *);

/* function pointers */
//pBucket (*newBucket)(pMesh ,int );
//int     (*buckin)(pMesh ,pBucket ,double *);
//int     (*locelt)(pMesh ,int ,double *,double *);
//int     (*closept)(pMesh ,double *);
//int     (*intpp0)(pMesh ,pSol ,double *,int ,double *,double );
//int     (*intpp1)(pSol ,int *,double *,int ,double *);
//int     (*hashelt)(pMesh );
//int     (*boulep)(pMesh ,int ,int ,int *);


#endif
