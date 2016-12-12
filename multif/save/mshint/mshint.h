#ifndef _MSHINT_H
#define _MSHINT_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "libmesh6.h"
#include "mi_calls.h"


#include "mesh.h"
#include "fproto.h"

#define MI_VER   "4.1a"
#define MI_REL   "Mar. 9, 2016"
#define MI_CPY   "(C) Copyright 1997- , ICS-SU"

#define MI_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MI_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MI_MIN3(a,b,c) ( (a) < (b) ? ((a)<(c) ? (a) : (c)) : ((b)<(c) ? (b) : (c)) )
#define MI_MAX3(a,b,c) ( (a) > (b) ? ((a)>(c) ? (a) : (c)) : ((b)>(c) ? (b) : (c)) )

#define BUCKSIZ   32
#define MI_EPS    1.e-6
#define MI_EPSD   1.e-200
#define MI_TGV    1.e30


typedef struct {
  double         c[3];
  int            s;
} Point;
typedef Point * pPoint;

typedef struct {
  int     v[3],adj[3],mark;
} Tria;
typedef Tria * pTria;

typedef struct {
  int     v[4],adj[4],mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  char     verb;
  mytime   ctim[TIMEMAX];
} Info;

typedef struct {
  double   min[3],max[3];
  int      dim,ver,np,nt,ne,mark;
  char    *name;
  pPoint   point;
  pTria    tria;
  pTetra   tetra;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
	int      dim,ver,np,it,type[2],size[2],typtab[2][GmfMaxTyp];
  double  *u,time;
  char    *name;
} Sol;
typedef Sol * pSol;

struct _MIst {
  Mesh   msrc,mtgt;
  Sol    ssrc,stgt;
  Info   info;
};

typedef struct {
  int   *cell,size;
} Bucket;


/* prototypes */
int  loadMesh(Mesh *mesh,char verb);
int  loadSol(Sol *sol,char verb);
int  saveSol(Sol *sol,char verb);
int  hashel_2d(MIst *mist);
int  hashel_3d(MIst *mist);
int  mshint1_2d(MIst *mist);
int  mshint1_3d(MIst *mist);
Bucket *bucket_2d(Mesh *mesh,int nmax);
Bucket *bucket_3d(Mesh *mesh,int nmax);
int  buckin_2d(Mesh *mesh,Bucket *bck,double *c);
int  buckin_3d(Mesh *mesh,Bucket *bck,double *c);
int  closept_2d(pMesh mesh,double *c);
int  closept_3d(pMesh mesh,double *c);
int  locelt_2d(pMesh mesh,int nsd,double *c,double *cb);
int  locelt_3d(pMesh mesh,int nsd,double *c,double *cb);


#endif
