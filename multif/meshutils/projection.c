#include "meshutils.h"

/*
Victorien Menier May 2017
*/

double GetThetaProj(double y, double z);
int GetEllipseNormal (double r1, double r2, double theta, double *nor);
int GetVecCur (double *crd0, double *crd1, double *vec);
int GetVecCur_theta (double *crd0, double r1, double r2, double theta, double *vec);
double NormalEllipseProjection (double r1, double r2, double *Crd, double t0);


int ProjectToDV_test(double *Crd, CadNozzle *Noz, double *BasParam, int Patch) 
{
	
	//--- Baseline parameters
	
	double cen_inlet[3]  = { BasParam[BasInletx], \
													 BasParam[BasInlety], \
													 BasParam[BasInletz]}; 
												
	double cen_outlet[3]  = { BasParam[BasOutletx], \
	                     		  BasParam[BasOutlety], \
	                     		  BasParam[BasOutletz]}; 
	
	double cut_inlet[3]  = {0,  0.439,  0.0999};
	double cut_outlet[3] = {2.33702,  0.88272, 0.122373};
	
	if ( Crd[0] < cen_inlet[0]-1e-6 || Crd[0] > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
		
	double zbas=0, r1_bas=0, r2_bas=0, alp=0;
	double ycut_bas=0, zcut_bas=0;
	double x=0, y=0, z=0;
	double r1=0, r2=0, dydx=0, zcenter=0;
	double ycut=0, zcut=0;
	double alpy=0, alpz=0;
	
	//--- Evaluate baseline at current x location
	
	x = Crd[0];
	y = Crd[1];
	z = Crd[2];
	
	alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	r1_bas = fr1_bas (x,  BasParam);
	r2_bas = fr2_bas (x,  BasParam);
	zbas   = fz_bas  (x,  BasParam);
	
	ycut_bas = alp*cut_outlet[1] + (1.0-alp)*cut_inlet[1];
	zcut_bas = alp*cut_outlet[2] + (1.0-alp)*cut_inlet[2];
	
	//--- Evaluate bsplines to project on
	
	CadBspline * Bsp_center = Noz->Bsp_center;
	CadBspline * Bsp_r1     = Noz->Bsp_r1;
	CadBspline * Bsp_r2     = Noz->Bsp_r2;
	
	//--- Centerline
	bSplineGeo3 (Bsp_center->Knots, Bsp_center->Coefs, &x, &zcenter, &dydx, 1, Bsp_center->NbrKnots, Bsp_center->NbrCoefs/2);
	
	//--- R1
	bSplineGeo3 (Bsp_r1->Knots, Bsp_r1->Coefs, &x, &r1, &dydx, 1, Bsp_r1->NbrKnots, Bsp_r1->NbrCoefs/2);
	
	//--- R2
	bSplineGeo3 (Bsp_r2->Knots, Bsp_r2->Coefs, &x, &r2, &dydx, 1, Bsp_r2->NbrKnots, Bsp_r2->NbrCoefs/2);
	
	Crd[1] = y/r1_bas*r1;
	Crd[2] = zcenter+(z-zbas)/r2_bas*r2;
	
	//if ( Patch == NOZZLEDOWN ) {
	//				
	//	ycut = ycut_bas/r1_bas*r1;
	//	zcut = zcenter+(zcut_bas-zbas)/r2_bas*r2;
	//	
	//	alpy = y/ycut_bas;
	//	alpz = z/zcut_bas;
	//	
	//	Crd[1] = alp*alpy*ycut + (1.0-alp)*Crd[1];
	//	Crd[2] = alp*alpz*zcut + (1.0-alp)*Crd[2];
	//	
	//}
	
	return 0;
}



int ProjectToDV_DV_test(double *Crd, CadNozzle *Noz, CadNozzle *NozBas, double *BasParam, int Patch) 
{
	
		
	double zbas=0, r1_bas=0, r2_bas=0, alp=0;
	double ycut_bas=0, zcut_bas=0;
	double x=0, y=0, z=0;
	double r1=0, r2=0, dydx=0, zcenter=0;
	double ycut=0, zcut=0;
	double alpy=0, alpz=0;
	double x_in, x_out;
	
	double zrou, r1_rou, r2_rou, ycut_rou, zcut_rou;
	
	//--- "rounded" baseline parameters
	
	double cen_inlet[3]  = { BasParam[BasInletx], \
													 BasParam[BasInlety], \
													 BasParam[BasInletz]}; 
												
	double cen_outlet[3]  = { BasParam[BasOutletx], \
	                     		  BasParam[BasOutlety], \
	                     		  BasParam[BasOutletz]}; 
	
	double cut_inlet[3]  = {0,  0.439,  0.0999};
	double cut_outlet[3] = {2.33702,  0.88272, 0.122373};
	
	//--- Evaluate baseline at current x location
	
	x = Crd[0];
	y = Crd[1];
	z = Crd[2];
	
	////alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	////
	////r1_bas = fr1_bas (x,  BasParam);
	////r2_bas = fr2_bas (x,  BasParam);
	////zbas   = fz_bas  (x,  BasParam);
		
	////ycut_bas = alp*cut_outlet[1] + (1.0-alp)*cut_inlet[1];
	////zcut_bas = alp*cut_outlet[2] + (1.0-alp)*cut_inlet[2];
	
	x_in   = NozBas->Bsp_center->Coefs[0];
	x_out  = NozBas->Bsp_center->Coefs[NozBas->Bsp_center->NbrCoefs/2-1];
	alp = (x-x_in)/(x_out-x_in);
		
	//--- Evaluate baseline nozzle
	Evaluate_Nozzle ( NozBas, &x, &r1_bas, &r2_bas, &zbas );
	
	//--- Evaluate current nozzle
	Evaluate_Nozzle ( Noz, &x, &r1, &r2, &zcenter );
	
	//--- Compute "blending" parameters (flat shovel)
	
	//alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	r1_rou = fr1_bas (x,  BasParam);
	r2_rou = fr2_bas (x,  BasParam);
	zrou   = fz_bas  (x,  BasParam);
	
	ycut_rou = alp*cut_outlet[1] + (1.0-alp)*cut_inlet[1];
	zcut_rou = alp*cut_outlet[2] + (1.0-alp)*cut_inlet[2];	
	
	ycut_bas = ycut_rou/r1_rou*r1;
	zcut_bas = zcenter+(zcut_rou-zrou)/r2_rou*r2;
		
	Crd[1] = y/r1_bas*r1;
	Crd[2] = zcenter+(z-zbas)/r2_bas*r2;
		
	if ( Patch == NOZZLEDOWN ) {
		
					
		ycut = ycut_bas/r1_bas*r1;
		zcut = zcenter+(zcut_bas-zbas)/r2_bas*r2;
		
		alpy = y/ycut_bas;
		alpz = z/zcut_bas;
		
		//Crd[1] = alp*alpy*ycut + (1.0-alp)*Crd[1];
		//Crd[2] = alp*alpz*zcut + (1.0-alp)*Crd[2];
					
	}
	
	return 0;
}


int NozzleWallProjection_test (Options *mshopt, Mesh *Msh, CadNozzle * CadNoz, int refUp, int refDown, char *OutNam)
{
	
	//---
	
	int WrtMsh = 1, patch=-1, WrtFullMsh=0;
	
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
	BasParam[BasInletz    ]  = 0.099908 ;
	BasParam[BasInletr    ]  = 0.439461   ;
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
			//ProjectNozzleWall_Up (Msh->Ver[iVer], CrdNew, BasParam);			
			patch = NOZZLEUP;
		}
		else if ( Tag[iVer] == refDown ) {
			
			//---
			//ProjectNozzleWall_Down_Save (Msh->Ver[iVer], CrdNew, BasParam);
			patch = NOZZLEDOWN;
		}
		
		if ( patch == -1 ) {
			printf("  ## ERROR : Wrong nozzle CAD patch.\n");
			exit(1);
		}
				
		ProjectToDV_test(CrdNew, CadNoz, BasParam, patch);
		
		if ( WrtMsh == 1 ) {
			MshViz->NbrVer++;
			AddVertex(MshViz, MshViz->NbrVer, CrdNew);
			Tag[iVer] = MshViz->NbrVer;
		}
		
		if (WrtFullMsh == 1) {
			Msh->Ver[iVer][0] = CrdNew[0];
			Msh->Ver[iVer][1] = CrdNew[1];
			Msh->Ver[iVer][2] = CrdNew[2];
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
		WriteGMFMesh("visu2", MshViz, 1);
	}
	
	if (WrtFullMsh == 1) {
		WriteGMFMesh("visu_full", Msh, 1);
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



int NozzleVolumeProjection (Options *mshopt, Mesh *Msh, CadNozzle * CadNoz, CadNozzle *CadNoz_bas,  int refUp, int refDown, char *OutNam, int verbose)
{
	//---
	
	int WrtMsh = 1, patch=-1;
		
	int iTri, ref=0, i, j, iVer, vid=-1, is[3];
	
	double CrdNew[3];
	
	int NbrVerPrj=0 , NbrTriPrj=0;
	
	int *Tag = (int*) malloc(sizeof(int)*(Msh->NbrVer+1));
	memset(Tag, 0, sizeof(int)*(Msh->NbrVer+1));
	
	double BasParam[BasMaxKwd];
	memset(BasParam, 0.0, sizeof(double)*BasMaxKwd);
	
	int refVol=2;
	int iTet;
	
	BasParam[BasInletx    ]  = 0.0       ;
	BasParam[BasInlety    ]  = 0.0       ;
	BasParam[BasInletz    ]  = 0.099908 ;
	BasParam[BasInletr    ]  = 0.439461   ;
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
	
	//--- Tag vertices that belong to volume 2
	
	refVol = 2;
	
	NbrTriPrj = 0;
	for (iTet=1; iTet<=Msh->NbrTet; iTet++) {
		ref = Msh->Tet[iTet][4];
		if ( ref != refVol )
			continue;
		for (j=0; j<4; j++) {
			vid = Msh->Tet[iTet][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	}
	
	printf("%d vertices belong to volume %d\n", NbrTriPrj, refVol);
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] <= 0 )
			continue;
		
		CrdNew[0] = Msh->Ver[iVer][0];
		CrdNew[1] = Msh->Ver[iVer][1];
		CrdNew[2] = Msh->Ver[iVer][2];
		
		if ( !CadNoz_bas )
			ProjectToDV_test(CrdNew, CadNoz, BasParam, patch);
		else
			ProjectToDV_DV_test(CrdNew, CadNoz, CadNoz_bas, BasParam, patch);
				
		Msh->Ver[iVer][0] = CrdNew[0];
		Msh->Ver[iVer][1] = CrdNew[1];
		Msh->Ver[iVer][2] = CrdNew[2];
	}
	
	//--- Write Mesh
	
	
	int FilTyp = GetInputFileType(OutNam);
  char *ptr = NULL;
	char BasNam[1024], OutSol[1024];
	
	// --- Get BasNam
	
  strcpy(BasNam,OutNam);
	
  ptr = strstr(BasNam,".su2");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
  ptr = strstr(BasNam,".meshb");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
	
	if ( FilTyp != FILE_SU2 ) {
		WriteGMFMesh(BasNam, Msh, 1);
	}
	else {
		WriteSU2Mesh(BasNam, Msh);
	}
	
	
	if ( Tag )
		free(Tag);
	
	if ( CadNoz )
		FreeCadNozzle (CadNoz);
	
	return 0;
	
}


