#include "meshutils.h"

/*
Victorien Menier Feb 2016
*/

int ConvertGMFtoSU2Sol (Options *mshopt)
{
	Mesh *Msh = NULL;
	char OutSol[1024];

	FILE *FilHdl=NULL;
	char *tok=NULL, *lin=NULL;

	int skip=0, SolSiz=0;
	size_t  len = 0;

	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);

	if ( Msh->FilTyp != FILE_GMF ) {
		printf("  ## ERROR : Input mesh file must be a .mesh (GMF) (FilTyp=%d)\n", Msh->FilTyp);
		return 0;
	}

	PrintMeshInfo (Msh);

	//if ( !Msh->Sol ) {
	//	printf(" ## ERROR : A solution must be provided\n");
	//	exit(1);
	//}



	//--- Get header information from reference file

	if ( strcmp(mshopt->HeaderNam, "") ) {

		FilHdl = fopen(mshopt->HeaderNam,"r");

		if ( !FilHdl ) {
			printf("  ## ERROR: ConvertGMFtoSU2Sol: UNABLE TO OPEN %s. \n", mshopt->HeaderNam);
	    return 0;
		}

		if ( getline(&lin, &len, FilHdl) != -1 ) {
			tok = strtok (lin, "	,");
			skip = 0;
			SolSiz = 0;
			while ( tok ) {
				if ( !strcmp(tok,"\"PointID\"") || !strcmp(tok,"\"x\"") || !strcmp(tok,"\"y\"") || !strcmp(tok,"\"z\"")   ) {
					tok = strtok (NULL, "	,");
					skip++;
					continue;
				}

				strcpy(Msh->SolTag[SolSiz], tok);
				//Str2Lower(Msh->SolTag[SolSiz]);
				StrRemoveChars(Msh->SolTag[SolSiz], '\"');
				StrRemoveChars(Msh->SolTag[SolSiz], '\n');
				SolSiz++;

				if ( SolSiz > Msh->SolSiz ) {
					printf("  ## ERROR: ConvertGMFtoSU2Sol: Provided header size does not match.\n");
			    return 0;
				}

				tok = strtok (NULL, "	,");
			}
	  }

	}

	if ( mshopt->clean == 1 )
		RemoveUnconnectedVertices(Msh);
	
	if ( mshopt->Dim == 2 )
		Msh->Dim = 2;
	
	WriteSU2Mesh(mshopt->OutNam, Msh);
	
	if ( Msh->Sol ) {
		sprintf(OutSol, "%s.dat", mshopt->OutNam);
		WriteSU2Solution (OutSol, Msh, Msh->Sol, Msh->NbrVer,  Msh->SolSiz, Msh->SolTag);
	}

	if ( Msh )
 		FreeMesh(Msh);

	return 1;
}

int ConvertSU2SolToGMF (Options *mshopt)
{
	Mesh *Msh = NULL;
	char OutSol[1024];

	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);

	if ( Msh->FilTyp != FILE_SU2 ) {
		printf("  ## ERROR : Input mesh file must be a .su2.\n");
		return 0;
	}

	PrintMeshInfo (Msh);

	WriteGMFMesh(mshopt->OutNam, Msh, 1);

	if ( Msh->Sol ) {
		sprintf(OutSol, "%s.solb", mshopt->OutNam);
		if ( ! WriteGMFSolutionItf(OutSol, Msh) ) {
			printf("  ## ERROR : outputmach FAILED.\n");
		}
	}

	if ( Msh )
 		FreeMesh(Msh);

	return 1;
}


int ConvertSU2ToGMSH (Options *mshopt)
{
	Mesh *Msh = NULL;
	char OutSol[1024];

	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);

	if ( Msh->FilTyp != FILE_SU2 ) {
		printf("  ## ERROR : Input mesh file must be a .su2\n");
		return 0;
	}

	PrintMeshInfo (Msh);

	WriteGMSHMesh(mshopt->OutNam, Msh);

	if ( Msh )
 		FreeMesh(Msh);

	return 1;
}



int Extraction (Options *mshopt)
{
	
	//--- Open mesh/solution file
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	PrintMeshInfo (Msh);
	
	
	if ( !Msh->Sol ) {
		printf("  ## ERROR SolutionExtraction : A solution must be provided.\n");
		return 0;
	}
		
	//SolutionExtraction(mshopt,Msh);
	
	double box[4] = {0.67,0.67,0,0.3048};
	
	int NbrRes=0, Siz=0;
	double *result = ExtractAlongLine(mshopt,Msh, box, &NbrRes, &Siz);
	
	//printf("NbrREs = %d\n", NbrRes);
	//
	//for (int i=0; i<NbrRes; i++) {
	//	printf ("res %d : ", i);
	//	for (int j=0; j<Siz; j++) {
	//		printf(" %lf", result[i*Siz+j]);
	//	}
	//	printf("\n");
	//}
	
	if ( result )
		free(result);
	
	if ( Msh )
 		FreeMesh(Msh);
	
	return 1;
	
}


int OutputMach (Options *mshopt)
{
	int NbrFld = 1, i, iVer;
	int FldTab[1] = {1};
	double *Mach = NULL;
	double *Pres = NULL;
	
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	PrintMeshInfo (Msh);
	WriteGMFMesh("mach", Msh, 1);	
	
	Mach = (double*)malloc(sizeof(double)*(Msh->NbrVer+1));
	Pres = (double*)malloc(sizeof(double)*(Msh->NbrVer+1));

	int iMach = -1;
	int iPres = -1;
	
	if ( Msh->FilTyp != FILE_SU2 ) {
		printf("  ## ERROR OutputMach : Input mesh file must be a .su2.\n");
		return 0;
	}
	
	for (i=0; i<Msh->NbrFld; i++) {
		if ( !strcmp(Msh->SolTag[i], "Mach") ) {
			iMach = i;
		}
		if ( !strcmp(Msh->SolTag[i], "Pressure") ) {
			iPres = i;
		}
	}
	
	if ( iMach < 0 ) {
		printf("  ## ERROR OutputMach : Mach index not found.\n");
		return 0;
	}
	
	if ( iPres < 0 ) {
		printf("  ## ERROR OutputMach : Pres index not found.\n");
		return 0;
	}
	
	printf("  -- Info : Mach number is field %d\n", iMach);
	printf("  -- Info : Pressure is field %d\n", iPres);
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {		
		Mach[iVer] = Msh->Sol[iVer*Msh->SolSiz+iMach];
		Pres[iVer] = Msh->Sol[iVer*Msh->SolSiz+iPres];
	}
		
	if ( ! WriteGMFSolution("mach.solb", Mach, 1, Msh->NbrVer, Msh->Dim, NbrFld, FldTab) ) {
		printf("  ## ERROR : outputmach FAILED.\n");
	}
	
	if ( ! WriteGMFSolution("pres.solb", Pres, 1, Msh->NbrVer, Msh->Dim, NbrFld, FldTab) ) {
		printf("  ## ERROR : output pressure FAILED.\n");
	}
		
	if ( Mach )
		free(Mach);
		
	if ( Pres )
		free(Pres);
	
	if ( Msh )
 		FreeMesh(Msh);
	
	return 1;
}



int CleanMesh (Options *mshopt)
{
	int NbrFld = 1, i, iVer;
	int FldTab[1] = {1};
	double *Mach = NULL;
	
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	PrintMeshInfo (Msh);

	RemoveUnconnectedVertices(Msh);
	
	Msh->NbrEfr = 0;
	
	WriteGMFMesh(mshopt->OutNam, Msh, 1);
	
	
	if ( Msh )
 		FreeMesh(Msh);
	
	return 1;
}


int ProjectNozzleWall (Options *mshopt)
{
	int NbrFld = 1, i, iVer;
	int FldTab[1] = {1};
	double *Mach = NULL;
	
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	PrintMeshInfo (Msh);

	//--- Load bsplines from data files (temporary)
	
	double BasParam[BasMaxKwd];
	memset(BasParam, 0.0, sizeof(double)*BasMaxKwd);
	// { 
	//	0.0, 0.0, 0.0999079,          // Crd of inlet circle center
	//	0.438972,                     // r inlet
	//	2.33702, 0.0, 0.19,           // Crd of outlet ellipse center
	//	0.92, 0.24,                   // r1 r2 outlet
	//	0.122638,                      // zcut out
	//	0.099,                         // zcut in
	//	0.879257,                     // ycut out
	//	0.438972                      // ycut in
	//};
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
	
	
	int SizCad[MaxKwdNozzleSize];
	memset(SizCad, 0, sizeof(int)*MaxKwdNozzleSize);
	
	CadNozzle *CadNoz = NULL;
	
	GetCadBsplineSize ("centerline", &SizCad[KwdnCoefs_center], &SizCad[KwdnKnots_center]);
	GetCadBsplineSize ("r1", &SizCad[KwdnCoefs_r1], &SizCad[KwdnKnots_r1]);
	GetCadBsplineSize ("r2", &SizCad[KwdnCoefs_r2], &SizCad[KwdnKnots_r2]);
	
	CadNoz = AllocCadNozzle (SizCad);
	
	LoadCadBspline ("centerline", CadNoz->Bsp_center);
	LoadCadBspline ("r1", CadNoz->Bsp_r1);
	LoadCadBspline ("r2", CadNoz->Bsp_r2);
	
	//double cosTheta;
	//cosTheta = max(-1.0, min(1.0,(BasParam[BasInletzcut]-BasParam[BasInletz])/BasParam[BasInletr]));
	//CadNoz->ThetaCutIn  = acos(cosTheta);
	//cosTheta = max(-1.0, min(1.0,(BasParam[BasOutletzcut]-BasParam[BasOutletz])/BasParam[BasOutletr2]));
	//CadNoz->ThetaCutOut  = acos(cosTheta);
	
	CadNoz->ThetaCutIn  = 1.572865;
	CadNoz->ThetaCutOut = 1.855294;
	
	WriteCadBspline ("centerline", CadNoz->Bsp_center);
	WriteCadBspline ("r1", CadNoz->Bsp_r1);
	WriteCadBspline ("r2", CadNoz->Bsp_r2);

	//--- End load bsplines from files

	int RefUp = 9;
	int RefDown = 10;
	
	NozzleWallProjection (mshopt ,Msh, CadNoz, RefUp, RefDown, "mesh_motion.dat");
	
	if ( Msh )
 		FreeMesh(Msh);
	
	return 1;
}

//
//int ConvertGMFtoSegMesh (Options *mshopt)
//{
//	Mesh *Msh = NULL;
//	char OutSol[1024];
//	
//	FILE *FilHdl=NULL;
//	char *tok=NULL, *lin=NULL;
//
//	int skip=0, SolSiz=0;
//	size_t  len = 0;
//
//	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
//
//	if ( Msh->FilTyp != FILE_GMF ) {
//		printf("  ## ERROR : Input mesh file must be a .mesh (GMF) (FilTyp=%d)\n", Msh->FilTyp);
//		return 0;
//	}
//
//	PrintMeshInfo (Msh);
//
//	if ( Msh->Dim != 2 ) {
//		printf("  ## ERROR : Mesh dimension must be 2.\n");
//		return 0;
//	}
//	
//	WriteSegMesh(mshopt->OutNam, Msh);
//	
//
//	if ( Msh )
// 		FreeMesh(Msh);
//
//	return 1;
//}
//


