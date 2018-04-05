#include "meshutils.h"

/*
Victorien Menier Feb 2016
*/

int Geo2Egads ()
{
	
	int CadSiz[CadKwdSize];
	Cad *cad = NULL;
	
	if ( ! GetGeoSize (mshopt->InpNam, CadSiz) ) {
		return 0;
	}	
	
	cad = AllocGeo (CadSiz);
	
	printf(" -- Info %s : %d vertices, %d lines, %d loops, %d surfaces\n", mshopt->InpNam,\
		CadSiz[CadVertices], CadSiz[CadLines], CadSiz[CadLoops], CadSiz[CadSurfaces]);
	
	ReadGeo(mshopt->InpNam, cad);
	
	//int iVer;
	//CadVertex *ver;
	//for (iVer=1; iVer<10; iVer++){
	//	ver = &cad->Ver[iVer];
	//	printf("Ver %d (vidold %d) : %lf %lf %lf\n", iVer, cad->VidOld[iVer], ver->Crd[0],ver->Crd[1],ver->Crd[2]);
	//}
	
	//int iLin, i;
	//CadLine *lin;
	//for (iLin=1; iLin<10; iLin++) {
	//	lin = &cad->Lin[iLin];
	//	printf("Bspline %d (%d): ", iLin, cad->LidOld[iLin]);
	//	for (i=0; i<lin->NbrCtr; i++) 
	//		printf("%d ", cad->VidOld[cad->CtrVer[lin->HeadCtr+i]]);
	//	printf("\n");
	//}
	
	//int iLoo, i;
	//CadLoop *loo;
	//for (iLoo=1; iLoo<5; iLoo++) {
	//	loo = &cad->Loo[iLoo];
	//	printf("Loop %d : ", iLoo);
	//	for (i=0; i<loo->NbrLin; i++) 
	//		printf("%d ", loo->Lin[i]);
	//	printf("\n");
	//}
	
	//WriteEgadsSurf ("cad.igs", cad);
	WriteEgads ("cad.igs", cad);
	
	cad = FreeGeo (cad);
	
	return 1;
	
}


/*
spline1d

Copied from the ESP source code
Cf OpenCSM/udpFreeForm.c

*/
static int spline1d(ego    context, int    imax, double *x, double *y, double *z, ego    *ecurv)
{
    int    status = EGADS_SUCCESS;

    int    i, kk, iknot, icp, iter, niter, header[4];
    double *cp=NULL, du, dx, dy, dz, data[18], dxyzmax;
    double *knotu=NULL, *CP=NULL;
    double dxyztol = 1.0e-7, relax = 0.10;

    icp   = imax + 2;
    iknot = imax + 6;

    /* indices associated with various arrays:
             x,y,z      knot       cp
               0         0*         0
                         1*
                                    1       (set by bc)
                         2*
                         3*
               1         4          2
               2         5          3

                        ...

             imax-1    imax       imax-2
             imax-2    imax+1     imax-1
                       imax+2*
                       imax+3*
                                  imax      (set by bc)
                       imax+4*
             imax-1    imax+5*    imax+1

        note: there are 4 repeateed knots at beginning and
              4 repeated knots at end */

    cp = (double*) EG_alloc((iknot+3*icp)*sizeof(double));
    if (cp == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    knotu = &(cp[0    ]);
    CP    = &(cp[iknot]);

    /* create spline curve */
    header[0] = 0;
    header[1] = 3;
    header[2] = icp;
    header[3] = iknot;

    kk = 0;

    /* knots (equally spaced) */
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;

    for (i = 1; i < imax; i++) {
        cp[kk++] = i;
    }

    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;

    /* initial control point */
    cp[kk++] = x[0];
    cp[kk++] = y[0];
    cp[kk++] = z[0];

    /* initial interior control point (for slope) */
    cp[kk++] = (3 * x[0] + x[1]) / 4;
    cp[kk++] = (3 * y[0] + y[1]) / 4;
    cp[kk++] = (3 * z[0] + z[1]) / 4;

    /* interior control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = x[i];
        cp[kk++] = y[i];
        cp[kk++] = z[i];
    }

    /* penultimate interior control point (for slope) */
    cp[kk++] = (3 * x[imax-1] + x[imax-2]) / 4;
    cp[kk++] = (3 * y[imax-1] + y[imax-2]) / 4;
    cp[kk++] = (3 * z[imax-1] + z[imax-2]) / 4;

    /* final control point */
    cp[kk++] = x[imax-1];
    cp[kk++] = y[imax-1];
    cp[kk++] = z[imax-1];

    /* make the original BSPLINE (based upon the assumed control points) */
    status = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                             header, cp, ecurv);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* iterate to have knot evaluations match data points */
    niter = 1000;
    for (iter = 0; iter < niter; iter++) {
        dxyzmax = 0;

        /* beginning point is fixed */

        /* match finite-differenced slope d/du at beginning */
        i = 1;
        status = EG_evaluate(*ecurv, &(knotu[3]), data);
        if (status != EGADS_SUCCESS) goto cleanup;

        du = knotu[4] - knotu[3];
        dx = x[1] - x[0] - du * data[3];
        dy = y[1] - y[0] - du * data[4];
        dz = z[1] - z[0] - du * data[5];


        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

        CP[3*(i)  ] += relax * dx;
        CP[3*(i)+1] += relax * dy;
        CP[3*(i)+2] += relax * dz;

        /* match interior points */
        for (i = 2; i < imax; i++) {
            status = EG_evaluate(*ecurv, &(knotu[i+2]), data);
            if (status != EGADS_SUCCESS) goto cleanup;

            dx = x[i-1] - data[0];
            dy = y[i-1] - data[1];
            dz = z[i-1] - data[2];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*(i)  ] += dx;
            CP[3*(i)+1] += dy;
            CP[3*(i)+2] += dz;
        }
				
        /* match finite-differenced slope d/du at end */
        i = imax;
        status = EG_evaluate(*ecurv, &(knotu[imax+2]), data);
        if (status != EGADS_SUCCESS) goto cleanup;
				 
        du = knotu[imax+2] - knotu[imax+1];
        dx = x[imax-1] - x[imax-2] - du * data[3];
        dy = y[imax-1] - y[imax-2] - du * data[4];
        dz = z[imax-1] - z[imax-2] - du * data[5];

        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

        CP[3*(i)  ] -= relax * dx;
        CP[3*(i)+1] -= relax * dy;
        CP[3*(i)+2] -= relax * dz;

        /* convergence check */
        if (dxyzmax < dxyztol) break;

    }

cleanup:
    if (cp != NULL) EG_free(cp);

    if (status != EGADS_SUCCESS) {
        if (ecurv != NULL) {
            EG_free(ecurv);
            *ecurv = NULL;
        }
    }

    return status;
}



/*
 ************************************************************************
 *                                                                      *
 *   spline2d - create 2d cubic spline (uniform spacing, fixed ends)    *
 *                                                                      *
 ************************************************************************
 */

static int
spline2d(ego context, int imax, int jmax, double *x, double *y, double *z, ego *esurf)
{
    int    status = EGADS_SUCCESS;

    int    i, j, kk, iknot, jknot, icp, jcp, iter, niter, header[7];
    double *cp=NULL, du, dv, dx, dy, dz, data[18], dxyzmax, parms[2];
    double *knotu=NULL, *knotv=NULL, *CP=NULL;
    double dxyztol = 1.0e-7, relax = 0.10;

    icp   = imax + 2;
    iknot = imax + 6;
    jcp   = jmax + 2;
    jknot = jmax + 6;

    cp = (double*) EG_alloc((iknot+jknot+3*icp*jcp)*sizeof(double));
    if (cp == NULL) {
        status  = EGADS_MALLOC;
        goto cleanup;
    }

    knotu = &(cp[0          ]);
    knotv = &(cp[iknot      ]);
    CP    = &(cp[iknot+jknot]);

    /* create spline surface */
    header[0] = 0;
    header[1] = 3;
    header[2] = icp;
    header[3] = iknot;
    header[4] = 3;
    header[5] = jcp;
    header[6] = jknot;

    kk = 0;

    /* knots in i-direction (equally spaced) */
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;

    for (i = 1; i < imax; i++) {
        cp[kk++] = i;
    }

    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;

    /* knots in j-direction (equally spaced) */
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;

    for (j = 1; j < jmax; j++) {
        cp[kk++] = j;
    }

    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;

    /* map of point ID for imax=9 and jmax=5 (used in comments below)

             4   nw O  n  n  n  n  n  n  n  P ne
                 J  K  L  L  L  L  L  L  L  M  N
             3   w  H  *  *  *  *  *  *  *  I  e
             2   w  H  *  *  *  *  *  *  *  I  e
             1   w  H  *  *  *  *  *  *  *  I  e
                 C  D  E  E  E  E  E  E  E  F  G
             0   sw A  s  s  s  s  s  s  s  B se

                 0     1  2  3  4  5  6  7     8
    */

    /* southwest control point */
    cp[kk++] = x[(0)+(0)*imax];
    cp[kk++] = y[(0)+(0)*imax];
    cp[kk++] = z[(0)+(0)*imax];

    /* point A */
    cp[kk++] = (3 * x[(0)+(0)*imax] + x[(1)+(0)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(0)*imax] + y[(1)+(0)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(0)*imax] + z[(1)+(0)*imax]) / 4;

    /* south control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = x[(i)+(0)*imax];
        cp[kk++] = y[(i)+(0)*imax];
        cp[kk++] = z[(i)+(0)*imax];
    }

    /* point B */
    cp[kk++] = (3 * x[(imax-1)+(0)*imax] + x[(imax-2)+(0)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(0)*imax] + y[(imax-2)+(0)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(0)*imax] + z[(imax-2)+(0)*imax]) / 4;

    /* southeast control point */
    cp[kk++] = x[(imax-1)+(0)*imax];
    cp[kk++] = y[(imax-1)+(0)*imax];
    cp[kk++] = z[(imax-1)+(0)*imax];

    /* point C */
    cp[kk++] = (3 * x[(0)+(0)*imax] + x[(0)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(0)*imax] + y[(0)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(0)*imax] + z[(0)+(1)*imax]) / 4;

    /* point D */
    cp[kk++] = (3 * x[(0)+(0)*imax] + x[(1)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(0)*imax] + y[(1)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(0)*imax] + z[(1)+(1)*imax]) / 4;

    /* points E */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = (3 * x[(i)+(0)*imax] + x[(i)+(1)*imax]) / 4;
        cp[kk++] = (3 * y[(i)+(0)*imax] + y[(i)+(1)*imax]) / 4;
        cp[kk++] = (3 * z[(i)+(0)*imax] + z[(i)+(1)*imax]) / 4;
    }

    /* point F */
    cp[kk++] = (3 * x[(imax-1)+(0)*imax] + x[(imax-2)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(0)*imax] + y[(imax-2)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(0)*imax] + z[(imax-2)+(1)*imax]) / 4;

    /* point G */
    cp[kk++] = (3 * x[(imax-1)+(0)*imax] + x[(imax-1)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(0)*imax] + y[(imax-1)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(0)*imax] + z[(imax-1)+(1)*imax]) / 4;

    /* loop through interior j lines */
    for (j = 1; j < jmax-1; j++) {

        /* west control point */
        cp[kk++] = x[(0)+(j)*imax];
        cp[kk++] = y[(0)+(j)*imax];
        cp[kk++] = z[(0)+(j)*imax];

        /* point H */
        cp[kk++] = (3 * x[(0)+(j)*imax] + x[(1)+(j)*imax]) / 4;
        cp[kk++] = (3 * y[(0)+(j)*imax] + y[(1)+(j)*imax]) / 4;
        cp[kk++] = (3 * z[(0)+(j)*imax] + z[(1)+(j)*imax]) / 4;

        /* interior points */
        for (i = 1; i < imax-1; i++) {
            cp[kk++] = x[(i)+(j)*imax];
            cp[kk++] = y[(i)+(j)*imax];
            cp[kk++] = z[(i)+(j)*imax];
        }

        /* point I */
        cp[kk++] = (3 * x[(imax-1)+(j)*imax] + x[(imax-2)+(j)*imax]) / 4;
        cp[kk++] = (3 * y[(imax-1)+(j)*imax] + y[(imax-2)+(j)*imax]) / 4;
        cp[kk++] = (3 * z[(imax-1)+(j)*imax] + z[(imax-2)+(j)*imax]) / 4;

        /* east control point */
        cp[kk++] = x[(imax-1)+(j)*imax];
        cp[kk++] = y[(imax-1)+(j)*imax];
        cp[kk++] = z[(imax-1)+(j)*imax];
    }

    /* point J */
    cp[kk++] = (3 * x[(0)+(jmax-1)*imax] + x[(0)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(jmax-1)*imax] + y[(0)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(jmax-1)*imax] + z[(0)+(jmax-2)*imax]) / 4;

    /* point K */
    cp[kk++] = (3 * x[(0)+(jmax-1)*imax] + x[(1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(jmax-1)*imax] + y[(1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(jmax-1)*imax] + z[(1)+(jmax-2)*imax]) / 4;

    /* points L */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = (3 * x[(i)+(jmax-1)*imax] + x[(i)+(jmax-2)*imax]) / 4;
        cp[kk++] = (3 * y[(i)+(jmax-1)*imax] + y[(i)+(jmax-2)*imax]) / 4;
        cp[kk++] = (3 * z[(i)+(jmax-1)*imax] + z[(i)+(jmax-2)*imax]) / 4;
    }

    /* point M */
    cp[kk++] = (3 * x[(imax-1)+(jmax-1)*imax] + x[(imax-2)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(jmax-1)*imax] + y[(imax-2)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(jmax-1)*imax] + z[(imax-2)+(jmax-2)*imax]) / 4;

    /* point N */
    cp[kk++] = (3 * x[(imax-1)+(jmax-1)*imax] + x[(imax-1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(jmax-1)*imax] + y[(imax-1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(jmax-1)*imax] + z[(imax-1)+(jmax-2)*imax]) / 4;

    /* northwest control point */
    cp[kk++] = x[(0)+(jmax-1)*imax];
    cp[kk++] = y[(0)+(jmax-1)*imax];
    cp[kk++] = z[(0)+(jmax-1)*imax];

    /* point O */
    cp[kk++] = (3 * x[(0)+(jmax-1)*imax] + x[(1)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(jmax-1)*imax] + y[(1)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(jmax-1)*imax] + z[(1)+(jmax-1)*imax]) / 4;

    /* north control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = x[(i)+(jmax-1)*imax];
        cp[kk++] = y[(i)+(jmax-1)*imax];
        cp[kk++] = z[(i)+(jmax-1)*imax];
    }

    /* point P */
    cp[kk++] = (3 * x[(imax-1)+(jmax-1)*imax] + x[(imax-2)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(jmax-1)*imax] + y[(imax-2)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(jmax-1)*imax] + z[(imax-2)+(jmax-1)*imax]) / 4;

    /* northeast control point */
    cp[kk++] = x[(imax-1)+(jmax-1)*imax];
    cp[kk++] = y[(imax-1)+(jmax-1)*imax];
    cp[kk++] = z[(imax-1)+(jmax-1)*imax];

    /* make the original BSPLINE (based upon the assumed control points) */
    status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                             header, cp, esurf);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* iterate to have knot evaluations match data points */
    niter = 1000;
    for (iter = 0; iter < niter; iter++) {
        dxyzmax = 0;

        /* south boundary */
        {
            j = 0;

            /* sw corner is fixed */

            /* point A to match FD d/du */
            i = 1;
            parms[0] = knotu[3];
            parms[1] = knotv[3];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[4] - knotu[3];
            dx = x[(1)+(0)*imax] - x[(0)+(0)*imax] - du * data[3];
            dy = y[(1)+(0)*imax] - y[(0)+(0)*imax] - du * data[4];
            dz = z[(1)+(0)*imax] - z[(0)+(0)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* interior points along south boundary */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i+2];
                parms[1] = knotv[  3];
                status = EG_evaluate(*esurf, parms, data);
                if (status != EGADS_SUCCESS) goto cleanup;

                dx = x[(i-1)+(0)*imax] - data[0];
                dy = y[(i-1)+(0)*imax] - data[1];
                dz = z[(i-1)+(0)*imax] - data[2];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += dx;
                CP[3*((i)+(j)*(imax+2))+1] += dy;
                CP[3*((i)+(j)*(imax+2))+2] += dz;
            }

            /* point B to match FD d/du */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[     3];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[imax+2] - knotu[imax+1];
            dx = x[(imax-1)+(0)*imax] - x[(imax-2)+(0)*imax] - du * data[3];
            dy = y[(imax-1)+(0)*imax] - y[(imax-2)+(0)*imax] - du * data[4];
            dz = z[(imax-1)+(0)*imax] - z[(imax-2)+(0)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* se corner is fixed */

        }

        /* line just above south boundary */
        {
            j = 1;

            /* point C to match FD d/dv */
            i = 0;
            parms[0] = knotu[3];
            parms[1] = knotv[3];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            dv = knotv[4] - knotv[3];
            dx = x[(0)+(1)*imax] - x[(0)+(0)*imax] - dv * data[6];
            dy = y[(0)+(1)*imax] - y[(0)+(0)*imax] - dv * data[7];
            dz = z[(0)+(1)*imax] - z[(0)+(0)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* point D to anhihilate d2/du/dv */
            i = 1;
            parms[0] = knotu[3];
            parms[1] = knotv[3];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[4] - knotu[3];
            dv = knotv[4] - knotv[3];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* points E to match FD d/dv */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i+2];
                parms[1] = knotv[  3];
                status = EG_evaluate(*esurf, parms, data);
                if (status != EGADS_SUCCESS) goto cleanup;

                dv = knotv[4] - knotv[3];
                dx = x[(i-1)+(1)*imax] - x[(i-1)+(0)*imax] - dv * data[6];
                dy = y[(i-1)+(1)*imax] - y[(i-1)+(0)*imax] - dv * data[7];
                dz = z[(i-1)+(1)*imax] - z[(i-1)+(0)*imax] - dv * data[8];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
                CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
                CP[3*((i)+(j)*(imax+2))+2] += relax * dz;
            }

            /* point F to annihilate d2/du/dv */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[     3];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[imax+2] - knotu[imax+1];
            dv = knotv[     4] - knotv[     3];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* point G to match FD d/dv */
            i = imax + 1;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[     3];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            dv = knotv[4] - knotv[3];
            dx = x[(imax-1)+(1)*imax] - x[(imax-1)+(0)*imax] - dv * data[6];
            dy = y[(imax-1)+(1)*imax] - y[(imax-1)+(0)*imax] - dv * data[7];
            dz = z[(imax-1)+(1)*imax] - z[(imax-1)+(0)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;
        }

        /* interior j lines */
        for (j = 2; j < jmax; j++) {

            /* interior point along west boundary */
            i = 0;
            parms[0] = knotu[  3];
            parms[1] = knotv[j+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            dx = x[(0)+(j-1)*imax] - data[0];
            dy = y[(0)+(j-1)*imax] - data[1];
            dz = z[(0)+(j-1)*imax] - data[2];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += dx;
            CP[3*((i)+(j)*(imax+2))+1] += dy;
            CP[3*((i)+(j)*(imax+2))+2] += dz;

            /* point H to match FD d/du */
            i = 1;
            parms[0] = knotu[  3];
            parms[1] = knotv[j+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[4] - knotu[3];
            dx = x[(1)+(j-1)*imax] - x[(0)+(j-1)*imax] - du * data[3];
            dy = y[(1)+(j-1)*imax] - y[(0)+(j-1)*imax] - du * data[4];
            dz = z[(1)+(j-1)*imax] - z[(0)+(j-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* interior points */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i+2];
                parms[1] = knotv[j+2];
                status = EG_evaluate(*esurf, parms, data);
                if (status != EGADS_SUCCESS) goto cleanup;

                dx = x[(i-1)+(j-1)*imax] - data[0];
                dy = y[(i-1)+(j-1)*imax] - data[1];
                dz = z[(i-1)+(j-1)*imax] - data[2];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += dx;
                CP[3*((i)+(j)*(imax+2))+1] += dy;
                CP[3*((i)+(j)*(imax+2))+2] += dz;
            }

            /* point I to match FD d/du */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[j   +2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[imax+2] - knotu[imax+1];
            dx = x[(imax-1)+(j-1)*imax] - x[(imax-2)+(j-1)*imax] - du * data[3];
            dy = y[(imax-1)+(j-1)*imax] - y[(imax-2)+(j-1)*imax] - du * data[4];
            dz = z[(imax-1)+(j-1)*imax] - z[(imax-2)+(j-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* interior point along east boundary */
            i = imax + 1;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[j   +2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            dx = x[(imax-1)+(j-1)*imax] - data[0];
            dy = y[(imax-1)+(j-1)*imax] - data[1];
            dz = z[(imax-1)+(j-1)*imax] - data[2];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += dx;
            CP[3*((i)+(j)*(imax+2))+1] += dy;
            CP[3*((i)+(j)*(imax+2))+2] += dz;
        }

        /* line just below north boundary */
        {
            j = jmax;

            /* point J to match FD d/dv */
            i = 0;
            parms[0] = knotu[     3];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = x[(0)+(jmax-1)*imax] - x[(0)+(jmax-2)*imax] - dv * data[6];
            dy = y[(0)+(jmax-1)*imax] - y[(0)+(jmax-2)*imax] - dv * data[7];
            dz = z[(0)+(jmax-1)*imax] - z[(0)+(jmax-2)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* point K to annihilate d2/du/dv */
            i = 1;
            parms[0] = knotu[     3];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[     4] - knotu[     3];
            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* points L to match FD d/dv */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i   +2];
                parms[1] = knotv[jmax+2];
                status = EG_evaluate(*esurf, parms, data);
                if (status != EGADS_SUCCESS) goto cleanup;

                dv = knotv[jmax+2] - knotv[jmax+1];
                dx = x[(i-1)+(jmax-1)*imax] - x[(i-1)+(jmax-2)*imax] - dv * data[6];
                dy = y[(i-1)+(jmax-1)*imax] - y[(i-1)+(jmax-2)*imax] - dv * data[7];
                dz = z[(i-1)+(jmax-1)*imax] - z[(i-1)+(jmax-2)*imax] - dv * data[8];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
                CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
                CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;
            }

            /* point M to annihilate d2/du/dv */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[imax+2] - knotu[imax+1];
            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* point N to match FD d/dv */
            i = imax + 1;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = x[(imax-1)+(jmax-1)*imax] - x[(imax-1)+(jmax-2)*imax] - dv * data[6];
            dy = y[(imax-1)+(jmax-1)*imax] - y[(imax-1)+(jmax-2)*imax] - dv * data[7];
            dz = z[(imax-1)+(jmax-1)*imax] - z[(imax-1)+(jmax-2)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;
        }

        /* north boundary */
        {
            j = jmax + 1;

            /* nw corner is fixed */

            /* point O to match FD d/du */
            i = 1;
            parms[0] = knotu[     3];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[4] - knotu[3];
            dx = x[(1)+(jmax-1)*imax] - x[(0)+(jmax-1)*imax] - du * data[3];
            dy = y[(1)+(jmax-1)*imax] - y[(0)+(jmax-1)*imax] - du * data[4];
            dz = z[(1)+(jmax-1)*imax] - z[(0)+(jmax-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* interior points along north boundary */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i   +2];
                parms[1] = knotv[jmax+2];
                status = EG_evaluate(*esurf, parms, data);
                if (status != EGADS_SUCCESS) goto cleanup;

                dx = x[(i-1)+(jmax-1)*imax] - data[0];
                dy = y[(i-1)+(jmax-1)*imax] - data[1];
                dz = z[(i-1)+(jmax-1)*imax] - data[2];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += dx;
                CP[3*((i)+(j)*(imax+2))+1] += dy;
                CP[3*((i)+(j)*(imax+2))+2] += dz;
            }

            /* point P to match FD d/du */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            if (status != EGADS_SUCCESS) goto cleanup;

            du = knotu[imax+2] - knotu[imax+1];
            dx = x[(imax-1)+(jmax-1)*imax] - x[(imax-2)+(jmax-1)*imax] - du * data[3];
            dy = y[(imax-1)+(jmax-1)*imax] - y[(imax-2)+(jmax-1)*imax] - du * data[4];
            dz = z[(imax-1)+(jmax-1)*imax] - z[(imax-2)+(jmax-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* se corner is fixed */

        }

        /* convergence check */
        if (dxyzmax < dxyztol) break;

        /* make the new curve (after deleting old one) */
        status = EG_deleteObject(*esurf);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                                 header, cp, esurf);
        if (status != EGADS_SUCCESS) goto cleanup;
    }

cleanup:
    if (cp != NULL) EG_free(cp);

    if (status != EGADS_SUCCESS) {
        if (esurf != NULL) {
            EG_free(esurf);
            *esurf = NULL;
        }
    }
    return status;
}



int WriteEgads (char *Nam, Cad *cad)
{
	
  int    senses[300];
  double xyz[3], dum[3], range[2];
  ego    context, objs[10];
  ego    model = NULL;
	
	char cmd[512];
		
	int stat, i, iLin, vid, NbrLin=0, d, NbrNod=0, nid[2], Vid[2], iLoo, lid;
	CadLine   *lin = NULL;
	CadVertex *ver = NULL;
	CadLoop   *loo = NULL;
	double *x=NULL, *y=NULL, *z=NULL;
	
	ego *Bsp=NULL, *edges=NULL, *nodes=NULL, *loops=NULL, *bodies=NULL, *faces=NULL, edgbuf[100], surface;
	
	int *Tag=NULL;
	
	x = (double*) malloc(sizeof(double) * cad->NbrCtrVer);
	y = (double*) malloc(sizeof(double) * cad->NbrCtrVer);
	z = (double*) malloc(sizeof(double) * cad->NbrCtrVer);
	
	Tag = (int*) malloc(sizeof(int)*(cad->NbrVer+1));
	memset(Tag, 0, sizeof(int)*(cad->NbrVer+1));
	
  stat = EG_open(&context);

	if ( stat != 0 ) {
		printf("  ## ERROR : EG_open failed.\n ");
		exit(1);
	}
	
	//--- Count and alloc BSplines
	Bsp    = (ego*) malloc(sizeof(ego)*(cad->NbrLin+1));
	edges  = (ego*) malloc(sizeof(ego)*(cad->NbrLin+1));
	nodes  = (ego*) malloc(sizeof(ego)*2*(cad->NbrLin+1));
	loops  = (ego*) malloc(sizeof(ego)*(cad->NbrLin+1));
	bodies = (ego*) malloc(sizeof(ego)*(cad->NbrLoo+1));	
	faces  = (ego*) malloc(sizeof(ego)*(cad->NbrLoo+1));	
	
	//--- Create Lines
	NbrLin = 0;
	for (iLin=1; iLin<=cad->NbrLin; iLin++) {
		
		lin = &cad->Lin[iLin];
		
		NbrLin++;
		
		//--- Create nodes at the two ends if they don't exist already
		
		nid[0] = nid[1] = 0; //indices of both nodes in nodes[]
		
		Vid[0] = cad->CtrVer[lin->HeadCtr];
		Vid[1] = cad->CtrVer[lin->HeadCtr+lin->NbrCtr-1];
		
		for (i=0; i<2; i++) {
			
			if ( Tag[Vid[i]] > 0 ) {
				// node has already been created
				nid[i] = Tag[Vid[i]];
				continue;
			}
			
			// Create node
			NbrNod++;
			
			ver = &cad->Ver[Vid[i]];
			stat = EG_makeTopology(context, NULL, NODE, 0, ver->Crd, 0, NULL, NULL, &nodes[NbrNod]);
						
			if ( stat != 0 ) {
				printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
				exit(1);
			}
			
			Tag[Vid[i]] = NbrNod;
			nid[i] = NbrNod;
		}
		
		if ( lin->Typ == CADBSPLINE || lin->Typ == CADLINE ){
			
			//--- Create spline	
			
			for (i=0; i<lin->NbrCtr; i++) {			
				vid = cad->CtrVer[lin->HeadCtr+i];
				ver = &cad->Ver[vid];
				
				x[i] = ver->Crd[0];
				y[i] = ver->Crd[1];
				z[i] = ver->Crd[2];
			}
			
			
			//status = spline2d(context, KMAX, JMAX, x2d, y2d, z2d,               \
	    //                  &(esurfs[IFACE]));
			//    status = spline1d(context, IMAX, x2d, y2d, z2d, &(ecurvs[IEDGE]));  
			
			//stat = spline2d(context,  lin->NbrCtr, 1, x, y, z, &surface);
			//printf("STAT = %d\n", stat);
			//exit(1);
			
			stat = spline1d(context, lin->NbrCtr, x, y, z, &Bsp[NbrLin]);
			
			if ( stat != 0 ) {
				printf("  ## ERROR : spline1d failed.\n ");
				exit(1);
			}
			
		}
		//else if ( lin->Typ == CADLINE ) {	
		//	vid = cad->CtrVer[lin->HeadCtr];
		//	ver = &cad->Ver[vid];
		//	for (d=0; d<3; d++) {
		//		data[d] = ver->Crd[d];
		//	}
		//	
		//	vid = cad->CtrVer[lin->HeadCtr+lin->NbrCtr-1];
		//	ver = &cad->Ver[vid];
		//	for (d=0; d<3; d++) {
		//		data[d+3] = ver->Crd[d];
		//	}
		//	
		//  stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &Bsp[NbrLin]);
		//	
		//	if ( stat != 0 ) {
		//		printf("  ## ERROR : EG_makeGeometry failed (error %d).\n ", stat);
		//		exit(1);
		//	}
		//}
		else {
			printf("  ## ERROR WriteEgads: UNKNOWN LINE TYPE\n");
			exit(1);
		}
		
		//--- Make topo : edge
		
		for (i=0; i<2; i++) {
			
			ver = &cad->Ver[Vid[i]];
			for (d=0; d<3; d++)
				xyz[d] = ver->Crd[d];
			
			stat = EG_invEvaluate(Bsp[NbrLin], xyz, &range[i],dum);
			
			if ( stat != 0 ) {
				printf("  ## ERROR : EG_invEvaluate failed (error %d).\n ", stat);
				exit(1);
			}
		}
		
		objs[0] = nodes[nid[0]];
	  objs[1] = nodes[nid[1]];
		
		stat = EG_makeTopology(context, Bsp[NbrLin], EDGE,TWONODE, range, 2, objs, NULL, &edges[NbrLin]);
		
		if ( stat != 0 ) {
			printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
			exit(1);
		}
		
		////--- Create loop (one loop per edge for now)
		//senses[0] = 1;
	  //stat = EG_makeTopology(context, NULL, LOOP, OPEN, NULL, 1, &edges[NbrLin], senses, &loops[NbrLin-1]);
	  //
		//if ( stat != 0 ) {
		//	printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
		//	exit(1);
		//}
		//
		////--- Create body
		//
		//stat = EG_makeTopology(context, NULL, BODY,WIREBODY, NULL, 1, &loops[NbrLin-1], NULL, &bodies[NbrLin-1]);
	  // 
	  //if ( stat != 0 ) {
	  //	printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
	  //	exit(1);
	  //}
	   
	}
	
	//--- Create loops
	
	for (iLoo=1; iLoo<=cad->NbrLoo; iLoo++) {
		loo = &cad->Loo[iLoo];
		
		for (i=0; i<loo->NbrLin; i++) {
			lid = abs(loo->Lin[i]);
						
			senses[i] = loo->Lin[i] / lid;
			
			lid = cad->LidNew[lid];
			edgbuf[i] = edges[lid];
			
		}
		
		stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, loo->NbrLin, edgbuf, senses, &loops[iLoo-1]);
		if ( stat != 0 ) {
			printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
			exit(1);
		}
		
		//printf(" EG_isoCline        = %d\n", EG_isoCline(loops[iLoo-1], 0, 0.0, &surface));
		//exit(1);
		
		//stat =  EG_makeFace(loops[iLoo-1], SFORWARD, NULL, &faces[iLoo-1]);
		//if ( stat != 0 ) {
		//	printf("  ## ERROR : EG_makeFace failed (error %d).\n ", stat);
		//	exit(1);
		//}
		//if ( iLoo == 3 )
		//{
		//	stat = EG_makeTopology(context, NULL, BODY,WIREBODY, NULL, 1, &loops[iLoo-1], NULL, &bodies[iLoo-1]);
		//	if ( stat != 0 ) {
	  //		printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
	  //		exit(1);
	  //	}
		//}
		
		stat = EG_makeTopology(context, NULL, BODY,WIREBODY, NULL, 1, &loops[iLoo-1], NULL, &bodies[iLoo-1]);
		if ( stat != 0 ) {
	  	printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
	  	exit(1);
	  }
	
		//ego body;
		//objs[0] = bodies[iLoo-1];
		//printf(" EG_ruled           = %d\n", EG_ruled(1, objs, &body));
		//exit(1);
	
	}
	
	//---  model
	
	printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0, NULL, cad->NbrLoo, bodies, NULL, &model));
	
	//--- Save and close
	
	sprintf(cmd, "rm -f %s", Nam);
	system(cmd);																			
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, Nam));
  printf(" EG_close           = %d\n", EG_close(context));

	//--- Free memory
	
	if (x) 
		free(x);
	
	if (y) 
		free(y);
		
	if (z) 
		free(z);
		
	if (Bsp)
		free(Bsp);
		
	if (edges)
		free(edges);
		
	if (Tag)
		free(Tag);
	
	if (nodes)
		free(nodes);
		
	if (loops)
		free(loops);
		
	if (bodies)
		free(bodies);
		
	if (faces)
		free(faces);
	
  return 0;
	
}


//
//int WriteEgadsSurf (char *Nam, Cad *cad)
//{
//	
//  int    senses[1000], sens;
//  double xyz[3], dum[3], range[2];
//  ego    context, objs[10];
//  ego    model = NULL;
//	double data[10];
//	char cmd[512];
//		
//	int stat, i, j, iLin, vid, NbrLin=0, d, NbrNod=0, nid[2], Vid[2], iLoo, lid, nbe, e, NbrEdg;
//	CadLine   *lin = NULL;
//	CadVertex *ver = NULL;
//	CadLoop   *loo = NULL;
//	double *x=NULL, *y=NULL, *z=NULL;
//	
//	ego *Bsp=NULL, *edges=NULL, *nodes=NULL, *loops=NULL, *bodies=NULL, *faces=NULL, edgbuf[1000], surface;
//	
//	int *Tag=NULL;
//	
//	x = (double*) malloc(sizeof(double) * cad->NbrCtrVer);
//	y = (double*) malloc(sizeof(double) * cad->NbrCtrVer);
//	z = (double*) malloc(sizeof(double) * cad->NbrCtrVer);
//	
//	Tag = (int*) malloc(sizeof(int)*(cad->NbrVer+1));
//	memset(Tag, 0, sizeof(int)*(cad->NbrVer+1));
//	
//  stat = EG_open(&context);
//
//	if ( stat != 0 ) {
//		printf("  ## ERROR : EG_open failed.\n ");
//		exit(1);
//	}
//	
//	//--- Count and alloc BSplines
//	Bsp    = (ego*) malloc(sizeof(ego)*(cad->NbrLin+1));
//	edges  = (ego*) malloc(sizeof(ego)*(cad->NbrLin+1));
//	nodes  = (ego*) malloc(sizeof(ego)*2*(cad->NbrCtrVer+2));
//	loops  = (ego*) malloc(sizeof(ego)*(cad->NbrLin+1));
//	bodies = (ego*) malloc(sizeof(ego)*(cad->NbrLoo+1));	
//	faces  = (ego*) malloc(sizeof(ego)*(cad->NbrLoo+1));	
//	
//	int NbrIntEdg = 0;
//
//	int *EdgIdx  = (int*) malloc(sizeof(int)*(cad->NbrLin+1));
//	int *HeadEdg = (int*) malloc(sizeof(int)*(cad->NbrCtrVer+2));
//	ego *IntEdg  = (ego*) malloc(sizeof(ego)*(cad->NbrCtrVer+1));
//	ego *IntLin  = (ego*) malloc(sizeof(ego)*(cad->NbrCtrVer+1));
//	
//	//--- Create Lines
//	NbrLin = 0;
//	for (iLin=1; iLin<=cad->NbrLin; iLin++) {
//		
//		
//		HeadEdg[iLin] = NbrIntEdg;
//		
//		lin = &cad->Lin[iLin];
//		
//		NbrLin++;
//		
//		//--- Create nodes at the two ends if they don't exist already
//		
//		if ( lin->Typ == CADBSPLINE || lin->Typ == CADLINE ){
//			
//			//--- Create edges out of the bspline	
//			
//			for (i=1; i<lin->NbrCtr; i++) {
//				
//				//--- Create node
//				Vid[0] = cad->CtrVer[lin->HeadCtr+i-1];
//				Vid[1] = cad->CtrVer[lin->HeadCtr+i];
//				
//				for (j=0; j<2; j++) {
//
//					ver = &cad->Ver[Vid[j]];
//					
//					if ( Tag[Vid[j]] > 0 ) {
//						// node has already been created
//						nid[j] = Tag[Vid[j]];
//						continue;
//					}
//					
//					NbrNod++;
//					
//					stat = EG_makeTopology(context, NULL, NODE, 0, ver->Crd, 0, NULL, NULL, &nodes[NbrNod]);
//
//					if ( stat != 0 ) {
//						printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
//						exit(1);
//					}
//
//					Tag[Vid[j]] = NbrNod;
//					nid[j] = NbrNod;
//				}
//				
//				//--- create line segment
//				ver = &cad->Ver[Vid[0]];
//				for (d=0; d<3; d++)
//					data[d] = ver->Crd[d];
//				
//				ver = &cad->Ver[Vid[1]];
//				for (d=0; d<3; d++)
//					data[d+3] = ver->Crd[d]-data[d];
//				
//				
//				if ( data[3] == data[4] && data[4] == data[5] && data[5] == 0 ) {
//					// null edge -> skip
//					printf("  ## WARNING : Null edge length -> skip.\n");
//					continue;
//				}
//								
//				stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &IntLin[NbrIntEdg]);
//				if ( stat != 0 ) {
//					printf("  ## ERROR : EG_makeTopology (LINE) failed (error %d).\n ", stat);
//					exit(1);
//				}
//				
//				//--- Create topo edge
//				for (j=0; j<2; j++) {
//					ver = &cad->Ver[Vid[j]];
//					for (d=0; d<3; d++)
//						xyz[d] = ver->Crd[d];
//					
//					stat = EG_invEvaluate(IntLin[NbrIntEdg], xyz, &range[j],dum);
//
//					if ( stat != 0 ) {
//						printf("  ## ERROR : EG_invEvaluate failed (error %d).\n ", stat);
//						exit(1);
//					}
//				}
//				
//
//				objs[0] = nodes[nid[0]];
//			  objs[1] = nodes[nid[1]];
//				
//				stat = EG_makeTopology(context, IntLin[NbrIntEdg], EDGE,TWONODE, range, 2, objs, NULL, &IntEdg[NbrIntEdg]);
//				
//				printf("EDG %d created: [%d, %d]\n", NbrIntEdg, Vid[0], Vid[1]);
//				
//				if ( stat != 0 ) {
//					printf("  ## ERROR : EG_makeTopology failed for EDGE (error %d).\n ", stat);
//					exit(1);
//				}
//				
//				NbrIntEdg++;
//			}
//			
//		}
//		else {
//			printf("  ## ERROR WriteEgads: UNKNOWN LINE TYPE\n");
//			exit(1);
//		}
//		
//	}
//	
//	
//	HeadEdg[iLin] = NbrIntEdg;
//	
//	for (iLin=1; iLin<=cad->NbrLin; iLin++) {
//		printf("HeadEdg[%d] = %d\n", iLin, HeadEdg[iLin]);
//	}
//	
//	for (iLoo=1; iLoo<=cad->NbrLoo; iLoo++) {
//		
//		printf("LOOP %d\n", iLoo);
//		
//		loo = &cad->Loo[iLoo];
//		
//		nbe = 0;
//			
//		for (iLin=0; iLin<loo->NbrLin; iLin++) {
//			
//			lid = abs(loo->Lin[iLin]);
//			sens = loo->Lin[iLin] / lid;
//			
//			printf("lid = %d -> %d\n", lid, cad->LidNew[lid]);
//			
//			lid = cad->LidNew[lid];
//
//			//--- Number of edges (created above) contained by lin
//			NbrEdg = HeadEdg[lid+1]-HeadEdg[lid];
//			
//			printf("Lin %d : NbrEdg = %d - %d\n", lid, HeadEdg[lid+1], HeadEdg[lid]);
//			
//			for (e=0; e<NbrEdg; e++) {
//				
//				edgbuf[nbe] = IntEdg[HeadEdg[lid]+e];
//				
//				printf("edg %d = %p\n", nbe, edgbuf[nbe]);
//				
//				senses[nbe] = sens;
//				nbe++;
//			}
//			
//		}
//		
//		printf("nbe = %d\n", nbe);
//		
//			
//		
//		stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, nbe, edgbuf, senses, &loops[iLoo-1]);
//		if ( stat != 0 ) {
//			printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
//			exit(1);
//		}
//		
//		printf(" EG_isoCline        = %d\n", EG_isoCline(loops[iLoo-1], 0, 0.0, &surface));
//		
//	}
//	
//	////--- Create loops
//	//
//	//for (iLoo=1; iLoo<=cad->NbrLoo; iLoo++) {
//	//	loo = &cad->Loo[iLoo];
//	//	
//	//	for (i=0; i<loo->NbrLin; i++) {
//	//		lid = abs(loo->Lin[i]);
//	//					
//	//		senses[i] = loo->Lin[i] / lid;
//	//		
//	//		lid = cad->LidNew[lid];
//	//		edgbuf[i] = edges[lid];
//	//		
//	//	}
//	//	
//	//	stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, loo->NbrLin, edgbuf, senses, &loops[iLoo-1]);
//	//	if ( stat != 0 ) {
//	//		printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
//	//		exit(1);
//	//	}
//	//	
//  //
//	//	//printf(" EG_isoCline        = %d\n", EG_isoCline(loops[iLoo-1], 0, 0.0, &surface));
//	//	//exit(1);
//	//	
//	//	//stat =  EG_makeFace(loops[iLoo-1], SFORWARD, NULL, &faces[iLoo-1]);
//	//	//if ( stat != 0 ) {
//	//	//	printf("  ## ERROR : EG_makeFace failed (error %d).\n ", stat);
//	//	//	exit(1);
//	//	//}
//	//	if ( iLoo == 3 )
//	//	{
//	//		stat = EG_makeTopology(context, NULL, BODY,WIREBODY, NULL, 1, &loops[iLoo-1], NULL, &bodies[iLoo-1]);
//	//		if ( stat != 0 ) {
//	//  		printf("  ## ERROR : EG_makeTopology failed (error %d).\n ", stat);
//	//  		exit(1);
//	//  	}
//	//	}
//	//}
//	//
//	////---  model
//	//
//	//printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0, NULL, cad->NbrLoo, bodies, NULL, &model));
//	//
//	////--- Save and close
//	//
//	//sprintf(cmd, "rm -f %s", Nam);
//	//system(cmd);																			
//  //printf(" EG_saveModel       = %d\n", EG_saveModel(model, Nam));
//  //printf(" EG_close           = %d\n", EG_close(context));
//  //
//	//--- Free memory
//	
//	if (x) 
//		free(x);
//	
//	if (y) 
//		free(y);
//		
//	if (z) 
//		free(z);
//		
//	if (Bsp)
//		free(Bsp);
//		
//	if (edges)
//		free(edges);
//		
//	if (Tag)
//		free(Tag);
//	
//	if (nodes)
//		free(nodes);
//		
//	if (loops)
//		free(loops);
//		
//	if (bodies)
//		free(bodies);
//		
//	if (faces)
//		free(faces);
//		
//	if (IntEdg)
//		free(IntEdg);
//		
//	if (IntLin)
//		free(IntLin);
//		
//	if (EdgIdx)
//		free(EdgIdx);
//		
//	if (HeadEdg)
//		free(HeadEdg);
//	
//  return 0;
//	
//}
//