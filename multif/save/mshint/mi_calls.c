#include "mshint.h"
#include "mi_calls.h"


MIst *MI_init(int dim,int ver) {
  MIst  *mist;
  
  /* default values */
  mist = (MIst*)calloc(1,sizeof(MIst));
  memset(&mist->msrc,0,sizeof(Mesh));
  memset(&mist->mtgt,0,sizeof(Mesh));
  memset(&mist->ssrc,0,sizeof(Sol));
  memset(&mist->stgt,0,sizeof(Sol));
  
  mist->info.verb  = '1';

  /* init timer */
  tminit(mist->info.ctim,TIMEMAX);
  chrono(ON,&mist->info.ctim[0]);

  return(mist);
}


/* free global data structure */
int MI_stop(MIst *mist) {
	char   stim[32];

	/* release memory */
  free(mist->msrc.point);
  if ( mist->msrc.nt )  free(mist->msrc.tria);
	if ( mist->msrc.ne ) free(mist->msrc.tetra);
  free(mist->mtgt.point);
  if ( mist->mtgt.nt )  free(mist->mtgt.tria);
	if ( mist->mtgt.ne )  free(mist->mtgt.tetra);
  if ( mist->ssrc.u )  free(mist->ssrc.u);

  chrono(OFF,&mist->info.ctim[0]);
  if ( mist->info.verb != '0' ) {
	  printim(mist->info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s sec.\n",stim);
  }

	return(1);
}


int MI_mshint(MIst *mist) {
  int   ier;

  if ( mist->msrc.dim == 2 )
    ier = mshint1_2d(mist);
  else
    ier = mshint1_3d(mist);

  return(ier);
}

