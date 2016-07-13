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
		
	SolutionExtraction(mshopt,Msh);
	
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





