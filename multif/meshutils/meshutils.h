#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <setjmp.h>

#define max(a,b) (a>=b?a:b)
#define min(a,b) (a<=b?a:b)

#include "libmesh6.h"
#include "mesh.h"

#include "cad.h"
#include "option.h"
#include "fproto.h"

//#include "global.h"

#define NOZMAXNBRPTS 1000
#define NOZMAXNBRTRI 100000
#define NOZMAXNBRVER 300000
#define NOZMAXNBREDG 50000
