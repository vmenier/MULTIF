#include "meshutils.h"

/*
Victorien Menier May 2017
*/


int ProjectNozzleWall_Up (double *CrdOld, double *CrdNew, double *BasParam)
{
	double r_inlet      = BasParam[3];              // r inlet
	
	double r_outlet[2] = {BasParam[7], BasParam[8]};            // r1 r2 outlet
	
	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]};
	
	double zbas = fz_bas (CrdOld[0],  BasParam);
	
	double cosTheta = 0.0, theta = 0.0, r1=0.0, r2=0.0, zcen=0.0;
	
	double r1_bas = fr1_bas (CrdOld[0],  BasParam);
	double r2_bas = fr2_bas (CrdOld[0],  BasParam);
	
	double rbas = 0.0;
	
	double alp = 0.0;
	
	double thetaMax = alp*BasParam[BasOutletTheta] + (1.0-alp)*BasParam[BasInletTheta];
	
	rbas += CrdOld[1]*CrdOld[1]/(r1_bas*r1_bas);
	rbas += (CrdOld[2]-zbas)*(CrdOld[2]-zbas)/(r2_bas*r2_bas);
	rbas = sqrt(rbas);
	
	cosTheta = (CrdOld[2]-zbas)/r2_bas/rbas;
	
	cosTheta = min(1.0, cosTheta);
	cosTheta = max(-1.0, cosTheta);
	
	theta = acos(cosTheta);
	
	r1 = 1.0;
	r2 = 1.0;
	zcen = 0.0;
	
	alp = (CrdOld[0]-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcen+r2*cos(theta);
	
	return 0;
	
}

double GetThetaProj(double y, double z)
{
	
	double  len, cosTheta, theta;
	
	len = sqrt(y*y+z*z);
	
	cosTheta = min(max(z/len, -1.0), 1.0);
	theta = acos(cosTheta);
	
	return theta;	
}

int ProjectNozzleWall_Down_Save (double *CrdOld, double *CrdNew, double *BasParam)
{
	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]}; 
	
	
	if ( CrdOld[0] < cen_inlet[0]-1e-6 || CrdOld[0] > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	double alp = (CrdOld[0]-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	double zbas = fz_bas (CrdOld[0],  BasParam);
	
	double cosTheta = 0.0, theta = 0.0, r1=0.0, r2=0.0, zcen=0.0, sinTheta=0.0;
	double y0, alp0, theta0, theta1;
	double r1_bas = fr1_bas (CrdOld[0],  BasParam);
	double r2_bas = fr2_bas (CrdOld[0],  BasParam);
	double zcut = fzcut (CrdOld[0],  BasParam);
	
	double z1, z2, z3;
	double nor[4];
	
	// ---- 
	
	theta1 = GetThetaProj(CrdOld[1], CrdOld[2]-zbas); 
	
	sinTheta = CrdOld[1]/r1_bas;
	sinTheta = min(1.0, sinTheta);
	sinTheta = max(-1.0, sinTheta);
	double theta2 = PI_NUMBER-asin(sinTheta);
	
	theta = alp*alp*theta2+(1-alp*alp)*theta1;
	
	r1 = 1.0;
	r2 = 1.0;
	zcen = 0.0;
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcen+r2*cos(theta);
	
	return 1;
	
}


double ComputeLength (double *crd0, double *crd1) 
{
	int j;
	double len=0.0;
	
	for (j=0; j<3; j++) {
		len += (crd1[j]-crd0[j])*(crd1[j]-crd0[j]);
	}
	
	return sqrt(len);
}

int GetVecCur (double *crd0, double *crd1, double *vec) 
{
	int j;
	double nrm=0.0;
	
	for (j=0; j<3; j++) {
		vec[j] = crd1[j]-crd0[j];
		nrm += vec[j]*vec[j];
	}
	
	nrm = sqrt(nrm);
	for (j=0; j<3; j++)
		vec[j] /= nrm;	
	
	return 1;
}


int ProjectNozzleWall_Down (double *CrdOld, double *CrdNew, double *BasParam)
{
	
	
	double cen_inlet[3]  = { BasParam[BasInletx], \
													 BasParam[BasInlety], \
													 BasParam[BasInletz]}; 
												
	double cen_outlet[3]  = { BasParam[BasOutletx], \
	                     		  BasParam[BasOutlety], \
	                     		  BasParam[BasOutletz]}; 
	
	
	if ( CrdOld[0] < cen_inlet[0]-1e-6 || CrdOld[0] > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	double alp = (CrdOld[0]-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	double zbas   = fz_bas  (CrdOld[0],  BasParam);
	double r1_bas = fr1_bas (CrdOld[0],  BasParam);
	double r2_bas = fr2_bas (CrdOld[0],  BasParam);
	double zcut   = fzcut   (CrdOld[0],  BasParam);
	
	//--- Perform dichotomy to find theta
	// ---- 
	
	double theta1 = GetThetaProj(CrdOld[1], CrdOld[2]-zbas); 
	
	double sinTheta = CrdOld[1]/r1_bas;
	sinTheta = min(1.0, sinTheta);
	sinTheta = max(-1.0, sinTheta);
	double theta2 = PI_NUMBER-asin(sinTheta);
	
	double theta = alp*alp*theta2+(1-alp*alp)*theta1;
	
	double thetaMin = alp*BasParam[BasOutletTheta] + (1.0-alp)*BasParam[BasInletTheta];
	
	//--- Compute thetaMin based on zcut, not thetacut
	
	double cosTheta = max(-1.0, min(1.0,(zcut-zbas)/r2_bas));
	thetaMin = acos(cosTheta);
	
	//printf("thetaMin %lf New %lf\n", thetaMin, theta);
	
	
	double thetaMax = PI_NUMBER-1e-6;
	
	double dTheta = 0.03*PI_NUMBER;
	
	int ite=0, maxIte=30, j;
	
	
	
	double t0 = thetaMin, t1 = thetaMax, t01 = 0.0;
	double l0=0.0, l1=0.0, n0=0.0, n1=0.0, nrm=0.0;
	
	t0 = max(thetaMin,theta-dTheta);
	t1 = min(thetaMax,theta+dTheta);
	
	printf("ini : t0 %lf t1 %lf\n", t0, t1);
	
	double crdCur[3] = {0.0,0.0,0.0};
	double crd[3] = {0.0,0.0,0.0};
	double vec[3] = {0.0,0.0,0.0};
	double nor[3]    = {0.0,0.0,0.0};
	
	crdCur[0] = CrdOld[0];
	crdCur[1] = CrdOld[1];
	crdCur[2] = CrdOld[2] - zbas;
	
	double t[2] = {0.0,0.0};
	double n[2] = {0.0,0.0};
	
	t[0] = max(thetaMin,theta-dTheta);
	t[1] = min(thetaMax,theta+dTheta);
	
	
	FILE *ellHdl = fopen("ellipse.dat", "wb");
	
	fprintf(ellHdl, "%le %le %le %le\n", r1_bas, r2_bas, crdCur[1], crdCur[2]);
	fclose(ellHdl);
	
	FILE *cvgHdl = fopen("cvgNor.dat", "wb");
	
	while( ite < maxIte ) {
		
		
		printf(" -- ITE %d : t0 %lf t1 %lf\n", ite, t[0], t[1]);
		
		for (j=0; j<2; j++) {
			
			//--- 
			
			crd[0] = CrdOld[0];
			crd[1] = r1_bas*sin(t[j]);
			crd[2] = r2_bas*cos(t[j]);
			
			//--- Normal to ellipse
			
			nor[0] = 0.0;
			nor[1] = -r1_bas*sin(t[j]);
			nor[2] = -r2_bas*cos(t[j]);
			nrm = sqrt(nor[1]*nor[1]+nor[2]*nor[2]);
			nor[1] /= nrm;
			nor[2] /= nrm;
			
			//--- Vector
			
			GetVecCur (crd, crdCur, vec);
			
			printf ("vec %d = %lf %lf\n", j, vec[1], vec[2]);
			printf ("nor %d = %lf %lf\n", j, nor[1], nor[2]);
			
			n[j] = (vec[1]-nor[1])*(vec[1]-nor[1])+(vec[2]-nor[2])*(vec[2]-nor[2]);
			fprintf(cvgHdl, "%le %le %le %le %le %le ", crd[1], crd[2], nor[1], nor[2], vec[1], vec[2]);
		}
		fprintf(cvgHdl, "\n");
		
		if ( fabs(n[0]-n[1]) < 1e-4 ) {
			theta = t[0];
			break;
		}
		
		if ( n[0] > n[1] ) {
			t[0] = 0.5*(t[0]+t[1]);
			theta = t[1];
		}
		else {
			t[1] = 0.5*(t[0]+t[1]);
			theta = t[0];
		}
		
		

		printf("n0 %lf n1 %lf\n", n[0], n[1]);
		
		ite++;
		
	}
	
	fclose(cvgHdl);
	
	
	//printf("thetamin %lf thetamax %lf theta %lf\n", thetaMin, thetaMax, theta);
	//
	//if ( theta < 0.21 )
	//	exit(1);
	

	if ( theta < thetaMin ) {
		printf("ERROR\n");
		exit(1);
	}
	
	//exit(1);
	
	// --- Project to cylinder
	
	double r1, r2, zcen;
	
	r1 = 1.0;
	r2 = 1.0;
	zcen = 0.0;
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcen+r2*cos(theta);
	
	return 1;
	
}


int ProjectToDV(double *Crd, CadNozzle *Noz, double *BasParam, int Patch) 
{
	
	int j=0;
	double  CrdNew[3] = {0.0,0.0,0.0}, alp=0.0, theta=0.0, zcut=0.0, thetaBas=0.0, ycut=0.0;
	
	CadBspline * Bsp_center = Noz->Bsp_center;
	CadBspline * Bsp_r1     = Noz->Bsp_r1;
	CadBspline * Bsp_r2     = Noz->Bsp_r2;
		
	double len = Bsp_center->Coefs[Bsp_center->NbrCoefs/2-1];
	double x = Bsp_center->Coefs[0]+Crd[0]*len, zcenter=0.0, r1=0.0, r2=0.0, dydx=0.0;
	
	//--- Centerline
	bSplineGeo3 (Bsp_center->Knots, Bsp_center->Coefs, &x, &zcenter, &dydx, 1, Bsp_center->NbrKnots, Bsp_center->NbrCoefs/2);
	
	//--- R1
	bSplineGeo3 (Bsp_r1->Knots, Bsp_r1->Coefs, &x, &r1, &dydx, 1, Bsp_r1->NbrKnots, Bsp_r1->NbrCoefs/2);
	
	//--- R2
	bSplineGeo3 (Bsp_r2->Knots, Bsp_r2->Coefs, &x, &r2, &dydx, 1, Bsp_r2->NbrKnots, Bsp_r2->NbrCoefs/2);
	
	//printf("x %lf ycenter %lf r1 %lf r2 %lf\n", x, zcenter, r1, r2);
	
	//--- Theta correction
	double thetaMaxBas, thetaMaxDV;
	
	alp = Crd[0];
	
	double xbas = Crd[0]*(BasParam[BasOutletx]-BasParam[BasInletx]);
	double zbas   = fz_bas  (xbas,  BasParam);
	double r2_bas = fr2_bas (xbas,  BasParam);
	zcut   = fzcut   (xbas,  BasParam);
	double cosTheta = max(-1.0, min(1.0,(zcut-zbas)/r2_bas));
	thetaMaxBas = acos(cosTheta);
	
	thetaMaxDV  = alp*Noz->ThetaCutOut + (1.0-alp)*Noz->ThetaCutIn;
	
	//printf("thetaMaxBas %lf thetaMaxDV %lf\n", thetaMaxBas, thetaMaxDV);
	
	CrdNew[0] = x;
	CrdNew[1] = r1*Crd[1];
	CrdNew[2] = zcenter+r2*Crd[2];
	
	if ( Patch == NOZZLEDOWN ) {
		
		theta = thetaMaxBas;//alp*Noz->ThetaCutOut + (1.0-alp)*Noz->ThetaCutIn;
		zcut  = zcenter+r2*cos(theta);		
		
	//	printf("alp %lf in %lf out %lf cut %lf zcut %lf theta(crd) = %lf\n", alp, Noz->ThetaCutIn/PI_NUMBER, Noz->ThetaCutOut/PI_NUMBER, theta/PI_NUMBER, zcut, acos(Crd[2])/PI_NUMBER);
		alp = Crd[0];
		//CrdNew[2] = zcut;
		CrdNew[2] = alp*zcut + (1.0-alp)*CrdNew[2];
	}
	
	//---
	
	Crd[0] = CrdNew[0];
	Crd[1] = CrdNew[1];
	Crd[2] = CrdNew[2];
	
	return 0;
	
}



int NozzleWallProjection (Options *mshopt, Mesh *Msh, CadNozzle * CadNoz, int refUp, int refDown, char *OutNam)
{
		
	//---
	
	int WrtMsh = 1, patch=-1;
		
	int iTri, ref=0, i, j, iVer, vid=-1, is[3];
	
	double CrdNew[3];
	
	int NbrVerPrj=0 , NbrTriPrj=0;
	
	int *Tag = (int*) malloc(sizeof(int)*(Msh->NbrVer+1));
	memset(Tag, 0, sizeof(int)*(Msh->NbrVer+1));
	
	//--- Parameters
	
	double BasParam[BasMaxKwd];
	memset(BasParam, 0.0, sizeof(double)*BasMaxKwd);
	
	BasParam[BasInletx    ]  = 0.0       ;
	BasParam[BasInlety    ]  = 0.0       ;
	BasParam[BasInletz    ]  = 0.0999079 ;
	BasParam[BasInletr    ]  = 0.438972  ;
	BasParam[BasOutletx   ]  = 2.33702   ;
	BasParam[BasOutlety   ]  = 0.0       ;
	BasParam[BasOutletz   ]  = 0.19      ;
	BasParam[BasOutletr1  ]  = 0.92      ;
	BasParam[BasOutletr2  ]  = 0.24      ;
	BasParam[BasOutletzcut]  = 0.122638  ;
	BasParam[BasInletzcut ]  = 0.099     ;
	BasParam[BasOutletycut]  = 0.879257  ;
	BasParam[BasInletycut ]  = 0.438972  ;
	BasParam[BasInletTheta]  = 1.572865  ;
	BasParam[BasOutletTheta] = 1.855294  ;
	
	//--- Tag one ref after the other for consistent line re-projection
	
	NbrTriPrj = 0;
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refDown )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	} 
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refUp  )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	}	
	
	//--- Create mesh for visu
	NbrVerPrj = 0;
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] > 0 )
			NbrVerPrj++;
	}
	
	int SizMsh[GmfMaxKwd+1];
	
	for (i=0; i<GmfMaxKwd; i++)
		SizMsh[i] = 0;
	
	SizMsh[GmfDimension]  = 3;
	SizMsh[GmfVertices]   = NbrVerPrj;
	SizMsh[GmfTriangles]  = NbrTriPrj;
	
	//printf("NbrVerPrj %d NbrRTriPrj %d\n", NbrVerPrj, NbrTriPrj);
	
	Mesh *MshViz = NULL;
	
	if ( WrtMsh == 1 ) {
		MshViz = AllocMesh(SizMsh);
		MshViz->Dim = 3;		
	}	
		
	FILE *motionHdl = fopen(OutNam, "wb");
	
	if ( motionHdl ) {
		printf("%%%% %s OPENED (WRITE).\n", OutNam);
	}
	else
	{
		printf("  ## ERROR : Could not open %s.\n", OutNam);
		exit(1);
	}
	
	//--- Project points
		
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
		if ( Tag[iVer] <= 0 )
			continue;
		
		CrdNew[0] = Msh->Ver[iVer][0];
		CrdNew[1] = Msh->Ver[iVer][1];
		CrdNew[2] = Msh->Ver[iVer][2];
		
		patch = -1;
		
		if ( Tag[iVer] == refUp ) {
			//--- 
			//fprintf(motionHdl, "%d %lf %lf %lf\n", iVer-1, Msh->Ver[iVer][0], Msh->Ver[iVer][1], Msh->Ver[iVer][2]);
			ProjectNozzleWall_Up (Msh->Ver[iVer], CrdNew, BasParam);			
			patch = NOZZLEUP;
			//printf("%d %lf %lf %lf\n", iVer-1, Msh->Ver[iVer][0], Msh->Ver[iVer][1], Msh->Ver[iVer][2]);
			
		}
		else if ( Tag[iVer] == refDown ) {
			//---
			ProjectNozzleWall_Down_Save (Msh->Ver[iVer], CrdNew, BasParam);
			patch = NOZZLEDOWN;
		}
		
		if ( patch == -1 ) {
			printf("  ## ERROR : Wrong nozzle CAD patch.\n");
			exit(1);
		}
		
		ProjectToDV(CrdNew, CadNoz, BasParam, patch);
		
		if ( WrtMsh == 1 ) {
			MshViz->NbrVer++;
			AddVertex(MshViz, MshViz->NbrVer, CrdNew);
			Tag[iVer] = MshViz->NbrVer;
		}
		
		fprintf(motionHdl, "%d %le %le %le\n", iVer-1, CrdNew[0], CrdNew[1], CrdNew[2]);
		
	}
	
	
	//--- Write visu mesh
	if ( WrtMsh == 1 ) {
		

		
		for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
			ref = Msh->Tri[iTri][3];
			if ( ref != refUp && ref != refDown )
				continue;
			
			for (j=0; j<3; j++) {
				vid = Msh->Tri[iTri][j];
				is[j] = Tag[vid];
			}
			
			MshViz->NbrTri++;
			AddTriangle(MshViz,MshViz->NbrTri,is,ref);
		}
		
		printf("%d ver, %d tri\n", MshViz->NbrVer, MshViz->NbrTri);
		WriteGMFMesh("visu", MshViz, 1);
	}
	
	if (motionHdl)
		fclose(motionHdl);
	
	if ( Tag )
		free(Tag);
		
	if ( MshViz )
		FreeMesh(MshViz);
	
	if ( CadNoz )
		FreeCadNozzle (CadNoz);
	
	return 0;
	
}
