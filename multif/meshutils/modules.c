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

Mesh * ExtractSurfacePatches(Mesh *Msh, int *Ref, int NbrRef) 
{
	
	Mesh *MshOut = NULL;
	
	int iVer, iTri, j, i, flag=0, is[4];
	int SizMsh[GmfMaxKwd+1];
	int idx, idxOut;
	int *Tag = (int*)malloc(sizeof(int)*(Msh->NbrVer+1));
	int *TagTri = (int*)malloc(sizeof(int)*(Msh->NbrTri+1));
	
	memset(Tag,0,sizeof(int)*(Msh->NbrVer+1));
	memset(TagTri,0,sizeof(int)*(Msh->NbrTri+1));
	memset(SizMsh,0,sizeof(int)*(GmfMaxKwd+1));
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {	
		flag = 0;
		for (i=0; i<NbrRef; i++) {
			if ( Msh->Tri[iTri][3] == Ref[i] ) {
				flag = 1;
				break;
			}
		}
	
		if ( flag == 0 )
			continue;
		
		TagTri[iTri] = 1;
		
		SizMsh[GmfTriangles]++;
		
		for (j=0; j<3; j++) {
			iVer = Msh->Tri[iTri][j];
			Tag[iVer] = 1;
		}	
	}
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] < 1 )
			continue;
		SizMsh[GmfVertices]++;
		Tag[iVer] = SizMsh[GmfVertices];
	}
	
	//--- Alloc mesh
	
	MshOut = AllocMesh(SizMsh);
	MshOut->Dim = 3;
	
	//--- Fill mesh
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] < 1 )
			continue;
		MshOut->NbrVer++;
		AddVertex(MshOut, MshOut->NbrVer,  Msh->Ver[iVer]);
		Tag[iVer] = MshOut->NbrVer;
	}
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		if ( TagTri[iTri] < 1 )
			continue;
		
		for (j=0; j<3; j++) {
			is[j] = Tag[Msh->Tri[iTri][j]];
		}
		
		MshOut->NbrTri++;
		AddTriangle(MshOut, MshOut->NbrTri, is, Msh->Tri[iTri][3]);
	}
	
	//--- Solution
	
	MshOut->SolSiz = Msh->SolSiz;
	MshOut->NbrFld = Msh->NbrFld;
	MshOut->FldTab = (int*) malloc(sizeof(int)*Msh->SolSiz);
	for (j=0; j<Msh->NbrFld; j++){
		MshOut->FldTab[j] = Msh->FldTab[j];
		strcpy(MshOut->SolTag[j],Msh->SolTag[j]);
	}
	MshOut->Sol = (double*) malloc(sizeof(double)*(MshOut->NbrVer+1)*MshOut->SolSiz);
	memset(MshOut->Sol, 0, sizeof(double)*(MshOut->NbrVer+1)*MshOut->SolSiz);
	
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] < 1 )
			continue;
		
		idx = iVer*Msh->SolSiz;
		idxOut = Tag[iVer]*MshOut->SolSiz;
		
		for (j=0; j<MshOut->SolSiz; j++) 
			MshOut->Sol[idxOut+j] = Msh->Sol[idx+j];
	}
	
	return MshOut;
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


