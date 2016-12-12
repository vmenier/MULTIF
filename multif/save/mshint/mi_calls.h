#ifndef __MI_CALLS_H
#define __MI_CALLS_H

/* data structure */
typedef struct _MIst MIst;

/* prototypes */
MIst *MI_init(int dim, int ver);
int   MI_stop(MIst *mist);

int   MI_mshint(MIst *mist);


#endif
