#include "meshutils.h"



int FreeCadBspline (CadBspline *Bsp) 
{
	
	if ( Bsp->Coefs )
		free(Bsp->Coefs);
		
 	if ( Bsp->Knots )
		free(Bsp->Knots); 
		
	return 1;
	
}

int FreeCadNozzle (CadNozzle *Noz) 
{
	
	if ( Noz )
		free(Noz);
	
	return 1;
	
}

int Evaluate_Nozzle ( CadNozzle *Noz, double *x, double *r1, double *r2, double *zcenter )
{
	
	double dydx=0.0;
	
	CadBspline * Bsp_center = Noz->Bsp_center;
	CadBspline * Bsp_r1     = Noz->Bsp_r1;
	CadBspline * Bsp_r2     = Noz->Bsp_r2;
	
	bSplineGeo3 (Bsp_center->Knots, Bsp_center->Coefs, x, zcenter, &dydx, 1, Bsp_center->NbrKnots, Bsp_center->NbrCoefs/2);
	bSplineGeo3 (Bsp_r1->Knots, Bsp_r1->Coefs, x, r1, &dydx, 1, Bsp_r1->NbrKnots, Bsp_r1->NbrCoefs/2);
	bSplineGeo3 (Bsp_r2->Knots, Bsp_r2->Coefs, x, r2, &dydx, 1, Bsp_r2->NbrKnots, Bsp_r2->NbrCoefs/2);
	
	return 1;
	
}



CadNozzle * AllocCadNozzle (int * SizCad)
{
	
	CadNozzle * Noz = malloc(sizeof(CadNozzle));
	
	Noz->ThetaCutIn  = 0.0;
	Noz->ThetaCutOut = 0.0;
	
	Noz->Bsp_center = NULL;
	Noz->Bsp_r1     = NULL;
	Noz->Bsp_r2     = NULL;
	
	//--- Centerline
	if ( SizCad[KwdnCoefs_center] > 0 && SizCad[KwdnKnots_center] > 0  ) {
		Noz->Bsp_center = AllocCadBspline (SizCad[KwdnCoefs_center], SizCad[KwdnKnots_center]);
	}
	else {
		printf("  ## ERROR  AllocCadNozzle : Incorrect coefs or knots number.\n");
		exit(1);
	}
	
	//--- R1
	if ( SizCad[KwdnCoefs_r1] > 0 && SizCad[KwdnKnots_r1] > 0  ) {
		Noz->Bsp_r1 = AllocCadBspline (SizCad[KwdnCoefs_r1], SizCad[KwdnKnots_r1]);
	}
	else {
		printf("  ## ERROR  AllocCadNozzle : Incorrect coefs or knots number.\n");
		exit(1);
	}
	
	//--- R2
	if ( SizCad[KwdnCoefs_r2] > 0 && SizCad[KwdnKnots_r2] > 0  ) {
		Noz->Bsp_r2 = AllocCadBspline (SizCad[KwdnCoefs_r2], SizCad[KwdnKnots_r2]);
	}
	else {
		printf("  ## ERROR  AllocCadNozzle : Incorrect coefs or knots number.\n");
		exit(1);
	}
	
	return Noz;
	
}

CadBspline * AllocCadBspline (int NbrCoefs, int NbrKnots)
{
	
	CadBspline * Bsp = (CadBspline*) malloc(sizeof(CadBspline));
	
	Bsp->NbrCoefs = NbrCoefs;
	Bsp->NbrKnots = NbrKnots;
	
	Bsp->Coefs = (double *) malloc(sizeof(double)*Bsp->NbrCoefs);
	Bsp->Knots = (double *) malloc(sizeof(double)*Bsp->NbrKnots);
	
	return Bsp;
	
}

int GetCadBsplineSize ( char *BasNam , int *NbrCoefs, int *NbrKnots )
{
	
	char KnotsNam[1024], CoefsNam[1024];
	
	sprintf(KnotsNam, "%s_knots.dat", BasNam);
	FILE *KnotsHdl = fopen(KnotsNam, "r");
  if ( !KnotsHdl ) {
    printf("  ## ERROR : Couldn't open %s\n", KnotsNam);
		exit(1);
  }

	fscanf(KnotsHdl, "%d" , NbrKnots);
	
	sprintf(CoefsNam, "%s_coefs.dat", BasNam);
	FILE *CoefsHdl = fopen(CoefsNam, "r");
  if ( !CoefsHdl ) {
    printf("  ## ERROR : Couldn't open %s\n", CoefsNam);
		exit(1);
  }
	
	fscanf(CoefsHdl, "%d" , NbrCoefs);

	return 1;
	
}

int LoadCadBspline (char *BasNam, CadBspline *Bsp, int verbose)
{
	
	int i;
	
	char KnotsNam[1024], CoefsNam[1024];
	
	sprintf(KnotsNam, "%s_knots.dat", BasNam);
	FILE *KnotsHdl = fopen(KnotsNam, "r");
  if ( !KnotsHdl ) {
    printf("  ## ERROR : Couldn't open %s\n", KnotsNam);
		exit(1);
  }

	printf ("%s opened.\n", KnotsNam);
	
	sprintf(CoefsNam, "%s_coefs.dat", BasNam);
	FILE *CoefsHdl = fopen(CoefsNam, "r");
  if ( !CoefsHdl ) {
    printf("  ## ERROR : Couldn't open %s\n", CoefsNam);
		exit(1);
  }

	printf ("%s opened.\n", CoefsNam);

	//--- Load coefs
	
	fscanf(CoefsHdl, "%d", &Bsp->NbrCoefs);
	Bsp->Coefs = (double *) malloc(sizeof(double)*Bsp->NbrCoefs);
	
	for (i=0; i<Bsp->NbrCoefs; i++) {
		fscanf(CoefsHdl, "%lf", &Bsp->Coefs[i]);
	}
	
	//--- Load knots
	
	fscanf(KnotsHdl, "%d", &Bsp->NbrKnots);
	Bsp->Knots = (double *) malloc(sizeof(double)*Bsp->NbrKnots);
	
	for (i=0; i<Bsp->NbrKnots; i++) {
		fscanf(KnotsHdl, "%lf", &Bsp->Knots[i]);
	}
	
	return 1;
	
}

int SetCadBspline (CadBspline *Bsp, double *Knots, int NbrKnots, double *Coefs, int NbrCoefs)
{
	
	int i;
	
	//--- Set coefs
		
	Bsp->NbrCoefs = NbrCoefs;	
	Bsp->Coefs = (double *) malloc(sizeof(double)*Bsp->NbrCoefs);
	
	for (i=0; i<Bsp->NbrCoefs; i++) {
		Bsp->Coefs[i] =  Coefs[i];
	}
	
	//--- Set knots
		
	Bsp->NbrKnots = NbrKnots;	
	Bsp->Knots = (double *) malloc(sizeof(double)*Bsp->NbrKnots);
	
	for (i=0; i<Bsp->NbrKnots; i++) {
		Bsp->Knots[i] =  Knots[i];
	}
	
	return 1;
	
}


double fzcut (double x,  double *BasParam)
{

	double alp = 0.0, z=0.0;

	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]}; 

	double zcut_in  = BasParam[10];
	double zcut_out = BasParam[BasOutletzcut];

	if ( x < cen_inlet[0]-1e-6 || x > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR fz_bas : x out of range!\n");
		exit(1);
	}

	alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);

	z = alp*zcut_out + (1.-alp)*zcut_in;
	
	return z;
}


double fycut (double x,  double *BasParam)
{

	double alp = 0.0, y=0.0;

	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]}; 

	double ycut_in  = BasParam[12];
	double ycut_out = BasParam[11];
	

	if ( x < cen_inlet[0]-1e-6 || x > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR fz_bas : x out of range!\n");
		exit(1);
	}

	alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);

 	y = alp*ycut_out + (1.-alp)*ycut_in;
	
	return y;
}


double fz_bas (double x,  double *BasParam)
{
	
	double alp = 0.0, z=0.0;
	
	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]}; 
	
	if ( x < cen_inlet[0]-1e-6 || x > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR fz_bas : x out of range!\n");
		exit(1);
	}
	
	alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	z = alp*cen_outlet[2] + (1.-alp)*cen_inlet[2];
	
	return z;
}

double fr1_bas (double x,  double *BasParam)
{
	
	double alp = 0.0, r1=0.0;
	
	double r1_inlet  = BasParam[3];
	double r1_outlet = BasParam[7];
	
	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]}; 
	
	if ( x < cen_inlet[0]-1e-6 || x > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR fz_bas : x out of range!\n");
		exit(1);
	}
	
	alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	r1 = alp*r1_outlet + (1.-alp)*r1_inlet;
	
	return r1;
}

double fr2_bas (double x,  double *BasParam)
{
	
	double alp = 0.0, r2=0.0;
	
	double r2_inlet  = BasParam[3];
	double r2_outlet = BasParam[8];
	
	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]}; 
	
	if ( x < cen_inlet[0]-1e-6 || x > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR fz_bas : x out of range!\n");
		exit(1);
	}
	
	alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	r2 = alp*r2_outlet + (1.-alp)*r2_inlet;
	
	return r2;
}


int WriteCadBspline(char *BasNam, CadBspline *Bsp, int verbose) {
	
	char CadNam [1024];
	double len = Bsp->Coefs[Bsp->NbrCoefs/2-1];
	int i;
	
	sprintf(CadNam, "%s_spl.dat", BasNam);
	
	FILE *outHdl = fopen(CadNam, "wb");
	
	if ( outHdl ) {
		if ( verbose > 0 )
			printf(" -- %s opened.\n", CadNam);
	}
	else {
		if ( verbose > 0 )
			printf(" ## WARNING: WriteCadBspline : Unable to open %s. Return.\n", CadNam);
		return 0;
	}
	
	double x[100], y[100], dydx[100];
	for (i=0; i<100; i++){
		x[i] = i*len/99;
	}	
	
	bSplineGeo3 (Bsp->Knots, Bsp->Coefs, x, y, dydx, 100, Bsp->NbrKnots, Bsp->NbrCoefs/2);
	
	for (i=0; i<100; i++){
		fprintf(outHdl, "%lf %lf \n", x[i], y[i]);
	}
	fclose(outHdl);
	
	return 0;
	
}

