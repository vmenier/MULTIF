/*
 * main program file for mshint
 * (C) Copyright 1997 - , ICS-SU
 *
 * This file is part of mshint.
 *
 * mshint is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * mshint is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with AdaptTools.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "mshint.h"
#include "mi_calls.h"


static void excfun(int sigid) {
  fprintf(stdout,"\n # unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout," abnormal stop\n");  break;
    case SIGBUS:
      fprintf(stdout," code error...\n");  break;
    case SIGFPE:
      fprintf(stdout," floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout," illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout," segmentation fault.\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout," programm killed.\n");  break;
  }
  fprintf(stdout," # no data file saved.\n"); 
  exit(1);
}


static void usage(char *prog) {
  fprintf(stdout,"\n usage: %s [+/-v | -h | -e] source[.mesh] [-s data[.sol]] target[.mesh] [-o output[.sol]]\n",prog);

  fprintf(stdout,"\nOptions and flags:\n\
  --help       show the syntax and exit.\n\
  --version    show the version and date of release and exit.\n\n\
  -v           suppress any message (for use with function call).\n\
  +v           increase the verbosity level for output.\n\n\
  source.mesh    name of the source mesh\n\
  target.mesh    name of the target mesh\n\
  data.sol       name of file containing the solution\n\
  output.sol     name of the output file\n");
  exit(1);
}


/* parsing arguments on command line */
static int parsar(int argc,char *argv[],MIst *mist) {
  int      i;
  char    *ptr,*data;

  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i ]== '+') ) {
      switch(argv[i][1]) {
      case '-':  /* on-line help */
        if ( !strcmp(argv[i],"--help") )
          usage(argv[0]);
        else if ( !strcmp(argv[i],"--version") ) {
          fprintf(stdout,"%s: version: %s release: %s\n",argv[0],MI_VER,MI_REL);
          exit(1);
        }
        break;
      case 'h':  /* on-line help */
      case '?':
	      usage(argv[0]);
	      break;
      case 'v':
        if ( !strcmp(argv[i],"-v") )
          mist->info.verb = '0';
        else if ( !strcmp(argv[i],"+v") )
          mist->info.verb = '+';
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]);
        }
        break;
      default:
	      fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
	      usage(argv[0]);
      }
    }
    else {
      if ( mist->msrc.name == NULL ) {
        data = (char*)calloc(strlen(argv[i])+10,sizeof(char));
        strcpy(data,argv[i]);
        ptr = strstr(data,".mesh");
        if ( !ptr )  strcat(data,".mesh");
        mist->msrc.name = data;
      }
      else if ( mist->mtgt.name == NULL ) {
        data = (char*)calloc(strlen(argv[i])+10,sizeof(char));
        strcpy(data,argv[i]);
        ptr = strstr(data,".mesh");
        if ( !ptr )  strcat(data,".mesh");
        mist->mtgt.name = data;
      }
      else {
        fprintf(stdout,"%s: illegal option %s\n",argv[0],argv[i]);
	      usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */
  if ( (mist->msrc.name == NULL) || (mist->mtgt.name == NULL) ) {
    if ( mist->info.verb != '0' )  fprintf(stderr,"%s: missing argument\n",argv[0]);
    usage(argv[0]);
  }
  if ( mist->ssrc.name == NULL) {
    data = (char*)calloc(strlen(mist->msrc.name)+10,sizeof(char));
    strcpy(data,mist->msrc.name);
    ptr = strstr(data,".mesh");
    if ( ptr )  *ptr = '\0';
    strcat(data,".sol");
    mist->ssrc.name = data;
  }
  
  return(1);
}


int main(int argc,char **argv) {
  MIst       mist;
  int        ier;
	char       stim[32];

  memset(&mist,0,sizeof(MIst));
  tminit(mist.info.ctim,TIMEMAX);
  chrono(ON,&mist.info.ctim[0]);
  
  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  signal(SIGBUS,excfun);

  /* init structure */
  memset(&mist.msrc,0,sizeof(Mesh));
  memset(&mist.mtgt,0,sizeof(Mesh));
  memset(&mist.ssrc,0,sizeof(Sol));
  memset(&mist.stgt,0,sizeof(Sol));

  /* default values */
  mist.info.verb  = '1';

  /* parse command line */
  if ( !parsar(argc,argv,&mist) )  return(1);

  /* loading data */
  chrono(ON,&mist.info.ctim[1]);

  if ( mist.info.verb != '0' ) {
    fprintf(stdout," - MSHINT, Release %s, %s\n   %s\n\n",MI_VER,MI_REL,MI_CPY);
    fprintf(stdout," - LOADING DATA\n");
  }

  /* loading meshes and solution */
  ier = loadMesh(&mist.msrc,mist.info.verb);
  if ( ier <= 0 )  return(1);
  ier = loadMesh(&mist.mtgt,mist.info.verb);
  if ( ier <= 0 )  return(1);
  ier = loadSol(&mist.ssrc,mist.info.verb);
  if ( ier <= 0 )  return(1);
  
  /* allocating memory */
  mist.stgt.dim  = mist.ssrc.dim;
  mist.stgt.ver  = mist.ssrc.ver;
  mist.stgt.np   = mist.mtgt.np;
  mist.stgt.it   = mist.ssrc.it;
  mist.stgt.time = mist.ssrc.time;
  memcpy(mist.stgt.type,mist.ssrc.type,2*sizeof(int));
  memcpy(mist.stgt.size,mist.ssrc.size,2*sizeof(int));
  memcpy(&mist.stgt.typtab[0],&mist.ssrc.typtab[0],mist.ssrc.type[0]*sizeof(int));
  memcpy(&mist.stgt.typtab[1],&mist.ssrc.typtab[1],mist.ssrc.type[1]*sizeof(int));
  if ( mist.stgt.np > 0 ) {
    mist.stgt.u = (double*)calloc(mist.stgt.np,mist.stgt.size[0]*sizeof(double));
    assert(mist.stgt.u);
  }

  /* set up adjacencies */
  if ( mist.msrc.dim == 2 )
    hashel_2d(&mist);
  else
    hashel_3d(&mist);

  chrono(OFF,&mist.info.ctim[1]);
	printim(mist.info.ctim[1].gdif,stim);
  if ( mist.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);

  if ( !mist.stgt.name ) {
    mist.stgt.name = (char *)calloc(128,sizeof(char));
    assert(mist.stgt.name);
    strcpy(mist.stgt.name,mist.mtgt.name);
  }

  /* L2 projection */
  chrono(ON,&mist.info.ctim[2]);
  if ( mist.info.verb != '0' )  fprintf(stdout,"\n ** MODULE MSHINT: %s\n",MI_VER);
  ier = MI_mshint(&mist);
  chrono(OFF,&mist.info.ctim[2]);
  if ( mist.info.verb != '0' ) {
		printim(mist.info.ctim[2].gdif,stim);
    if ( ier > 0 )  
      fprintf(stdout," ** COMPLETED: %s\n\n",stim);
    else
      fprintf(stdout," ** NOT COMPLETED!: %s\n\n",stim);
	}

  /* save file */
  if ( mist.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
  chrono(ON,&mist.info.ctim[3]);

  ier = saveSol(&mist.stgt,0);
	if ( !ier )   return(1);
  chrono(OFF,&mist.info.ctim[3]);
  if ( mist.info.verb != '0' ) {
    printim(mist.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }

  /* free mem */
  free(mist.msrc.point);
  if ( mist.msrc.nt )  free(mist.msrc.tria);
	if ( mist.msrc.ne )  free(mist.msrc.tetra);
  free(mist.mtgt.point);
  if ( mist.mtgt.nt )  free(mist.mtgt.tria);
	if ( mist.mtgt.ne )  free(mist.mtgt.tetra);
  if ( mist.ssrc.u )   free(mist.ssrc.u);

  chrono(OFF,&mist.info.ctim[0]);
  if ( mist.info.verb != '0' ) {
	  printim(mist.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}


