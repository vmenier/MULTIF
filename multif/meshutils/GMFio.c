#include "meshutils.h"

/*
Victorien Menier Feb 2016
*/

int AddGMFMeshSize (char *MshNam, int *SizMsh)
{
	int dim, FilVer, InpMsh, i; 
	
	for (i=0; i<GmfMaxKwd; i++)
		SizMsh[i] = 0;

  if ( !(InpMsh = GmfOpenMesh(MshNam,GmfRead,&FilVer,&dim)) ) {
    fprintf(stderr,"  ## ERROR: Mesh data file %s.mesh[b] not found ! ",MshNam);
		return 0;
  }
	
  //--- get number of entities
	SizMsh[GmfDimension]  = dim;
  SizMsh[GmfVertices]   = GmfStatKwd(InpMsh, GmfVertices);
  SizMsh[GmfTriangles]  = GmfStatKwd(InpMsh, GmfTriangles);
	SizMsh[GmfTetrahedra] = GmfStatKwd(InpMsh, GmfTetrahedra);
  SizMsh[GmfEdges]      = GmfStatKwd(InpMsh, GmfEdges);
	SizMsh[GmfPrisms]     = GmfStatKwd(InpMsh, GmfPrisms);
  SizMsh[GmfPyramids]   = GmfStatKwd(InpMsh, GmfPyramids);

  if ( SizMsh[GmfVertices] <= 0 ) {
    fprintf(stderr,"\n  ## ERROR: NO VERTICES. IGNORED\n");
		return 0;
  }
	
	return 1;
}

int LoadGMFMesh (char *MshNam, Mesh *Msh)
{
	int i;
	int dim, FilVer, InpMsh, ref; 
	double bufDbl[3];
	int bufInt[8], is[8];
	
	int NbrVer, NbrTri, NbrEfr, NbrTet, NbrHex, NbrQua;
	
  if ( !(InpMsh = GmfOpenMesh(MshNam,GmfRead,&FilVer,&dim)) ) {
    printf("  ## ERROR: Mesh data file %s.mesh[b] not found ! \n",MshNam);
		return 0;
  }

	Msh->Dim = dim;
	strcpy(Msh->MshNam, MshNam);
	Msh->FilTyp = FILE_GMF;

  printf("  %%%% %s OPENED (READ)\n",MshNam);
		
	Msh->NbrVer = Msh->NbrTri = Msh->NbrEfr = 0;
	Msh->NbrTet = Msh->NbrHex = Msh->NbrQua = 0;
	
	NbrVer = NbrTri = NbrEfr = 0;
	NbrTet = NbrHex = NbrQua = 0;
	
	//--- Read vertices
	NbrVer = GmfStatKwd(InpMsh, GmfVertices);	
	
	GmfGotoKwd(InpMsh, GmfVertices);
	
	if ( Msh->Dim == 2 ) {
		for (i=1; i<=NbrVer; ++i) {
			GmfGetLin(InpMsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
			Msh->NbrVer++;
			AddVertex(Msh, Msh->NbrVer, bufDbl);
  	}		
	}
	else {
		for (i=1; i<=NbrVer; ++i) {
			GmfGetLin(InpMsh, GmfVertices, &bufDbl[0], &bufDbl[1], &bufDbl[2], &ref);
			Msh->NbrVer++;
			AddVertex(Msh, Msh->NbrVer, bufDbl);
  	}	
	}


	//--- Read Triangles
	NbrTri = GmfStatKwd(InpMsh, GmfTriangles);	
	GmfGotoKwd(InpMsh, GmfTriangles);
  for (i=1; i<=NbrTri; ++i) {
    GmfGetLin(InpMsh, GmfTriangles, &bufInt[0], &bufInt[1], &bufInt[2], &ref);
		Msh->NbrTri++;
		switchTriIdx(bufInt,is);
		AddTriangle(Msh,Msh->NbrTri,is,ref);
  }
	
	//--- Read boundary edges
	NbrEfr = GmfStatKwd(InpMsh, GmfEdges);	
	GmfGotoKwd(InpMsh, GmfEdges);
  for (i=1; i<=NbrEfr; ++i) {
		GmfGetLin(InpMsh, GmfEdges, &bufInt[0], &bufInt[1], &ref);
		Msh->NbrEfr++;
		AddEdge(Msh,Msh->NbrEfr,bufInt,ref);
  }

	//--- Read tetrahedra
	NbrTet = GmfStatKwd(InpMsh, GmfTetrahedra);	
	GmfGotoKwd(InpMsh, GmfTetrahedra);
  for (i=1; i<=NbrTet; ++i) {
		GmfGetLin(InpMsh, GmfTetrahedra, &bufInt[0], &bufInt[1], &bufInt[2], &bufInt[3], &ref);
		Msh->NbrTet++;
		AddTetrahedron(Msh,Msh->NbrTet,bufInt,ref);
  }
	
	if ( !GmfCloseMesh(InpMsh) ) {
    printf("  ## ERROR: Cannot close solution file %s ! \n",MshNam);
		return 0;
  }
	
	printf("NbrVer %d NbrTri %d NbrEfr %d\n", Msh->NbrVer, Msh->NbrTri, Msh->NbrEfr);
	
	return 1;
}

int LoadGMFSolution(char *SolNam, Mesh *Msh)
{
  int    SolMsh,FilVer=0,dim=0,SolTyp,iVer,i, idxVer;
  int    NbrLin,NbrTyp,SolSiz,TypTab[ GmfMaxTyp ];
  double *bufDbl = NULL;
		
	if ( Msh->Sol )
	{
		printf("  ## ERROR LoadGMFSolution : Msh->Sol already allocated.\n");
		return 0;
	}
	
  if ( (Msh == NULL) || (SolNam == NULL) ) {
    printf("  ## ERROR: LoadGMFSolution : MESH/FILE NAME NOT ALLOCATED \n");
    return 0; 
  }
	
	printf("SolNam %s GmfRead %d \n", SolNam, GmfRead);
  if ( !(SolMsh = GmfOpenMesh(SolNam,GmfRead,&FilVer,&dim)) ) {
    printf(" Solution data file %s.sol[b] not found ! \n",SolNam);
		return 0;
  }

  printf("  %%%% %s OPENED\n",SolNam);
	
	strcpy(Msh->SolNam, SolNam);

  if ( dim != 2 && dim != 3 ) {
    printf("  ## ERROR: WRONG DIMENSION NUMBER. IGNORED\n");
		return 0;
  }
	
  SolTyp = GmfSolAtVertices;
  NbrLin = GmfStatKwd(SolMsh, SolTyp, &NbrTyp, &SolSiz, TypTab);
	
	if ( NbrLin == 0 ) {
		printf("  ## ERROR LoadGMFSolution : No SolAtVertices in the solution file !\n");
		return 0;
	}
	
	if ( NbrLin != Msh->NbrVer ) {
		printf("  ## ERROR LoadGMFSolution : The number of vertices does not match (NbrLin %d, NbrVer %d).\n", NbrLin, Msh->NbrVer);
		return 0;
	} 
  
	//--- Allocate Msh->Sol
	
	Msh->Sol = (double*) malloc(sizeof(double)*(Msh->NbrVer+1)*SolSiz);
	memset(Msh->Sol, 0, sizeof(double)*(Msh->NbrVer+1)*SolSiz);
	
	Msh->SolSiz = SolSiz;
  
	Msh->NbrFld = NbrTyp;
	Msh->FldTab = (int*) malloc(sizeof(int)*SolSiz);
	for (i=0; i<NbrTyp; i++) {
		Msh->FldTab[i] = TypTab[i];
		switch (TypTab[i]) {
			case GmfSca:
			sprintf(Msh->SolTag[i], "scalar_%d", i);
			break;
						
			case GmfVec:
			sprintf(Msh->SolTag[i], "vector_%d", i);
			break;
			
			case GmfSymMat:
			sprintf(Msh->SolTag[i], "SymMatrix_%d", i);
			break;
			
			case GmfMat:
			sprintf(Msh->SolTag[i], "Matrix_%d", i);
			break;
			
			default:
			sprintf(Msh->SolTag[i], "field_%d", i);	
		}
	}
  //--- Read solution
  GmfGotoKwd(SolMsh, SolTyp);
  
  bufDbl = (double*)malloc(sizeof(double)*SolSiz);
  for (iVer=1; iVer<=Msh->NbrVer; ++iVer) {
    GmfGetLin(SolMsh, SolTyp, bufDbl);
		idxVer = iVer*SolSiz;
    for (i=0; i<SolSiz; i++)
			Msh->Sol[idxVer+i] = bufDbl[i];
  }
	free(bufDbl);
  
  if ( !GmfCloseMesh(SolMsh) ) {
    printf("  ## ERROR: Cannot close solution file %s ! \n",SolNam);
		return 0;
  }
	
	return 1;
}



int WriteGMFMesh(char *nam, Mesh *Msh, int OptBin)
{
  int       OutMsh,FilVer,i, j;
  int       iVer,iTri,iEfr,iTet; 
  long long idx[6];
  char      OutFil[512];
  
  int Dim = Msh->Dim;
	int NbrVer  = Msh->NbrVer;
	int NbrTri  = Msh->NbrTri;
	int NbrEfr  = Msh->NbrEfr;
	double3*Ver = Msh->Ver;
	int4*Tri    = Msh->Tri;
	int3*Efr    = Msh->Efr;
	int5*Tet    = Msh->Tet;
	
	
  //--- Define file name extension 
  strcpy(OutFil,nam);

  if ( OptBin == 1 )
    strcat(OutFil,".meshb");
  else
    strcat(OutFil,".mesh");
 
  FilVer = GmfDouble;
  
  //--- Open file
  if ( !(OutMsh = GmfOpenMesh(OutFil,GmfWrite,FilVer,Dim)) ) {
    printf("  ## ERROR: Cannot open mesh file %s ! \n",OutFil);
		return 0;
  }
  printf("  %%%% %s OPENED (WRITE)\n",OutFil);
  
  //--- Write vertices
  GmfSetKwd(OutMsh, GmfVertices, NbrVer);

	if ( Msh->Dim == 2 ) {
  	for (iVer=1; iVer<=NbrVer; ++iVer) {
    	GmfSetLin(OutMsh, GmfVertices,Ver[iVer][0],Ver[iVer][1],0);  
  	}
  }
	else {
		for (iVer=1; iVer<=NbrVer; ++iVer) {
    	GmfSetLin(OutMsh, GmfVertices,Ver[iVer][0],Ver[iVer][1],Ver[iVer][2],0);  
  	}
	}

  //--- Write triangles
  GmfSetKwd(OutMsh, GmfTriangles, NbrTri);
  for (iTri=1; iTri<=NbrTri; ++iTri) {
    for (i=0; i<3; ++i) {
      idx[i] = (long long)(Tri[iTri][i]);
    }
    GmfSetLin(OutMsh, GmfTriangles,idx[0],idx[1],idx[2],Tri[iTri][3]);  
  }

  //--- Write tetrahedra
  GmfSetKwd(OutMsh, GmfTetrahedra, Msh->NbrTet);
  for (iTet=1; iTet<=Msh->NbrTet; ++iTet) {
    for (i=0; i<4; ++i) {
      idx[i] = (long long)(Tet[iTet][i]);
    }
    GmfSetLin(OutMsh, GmfTetrahedra,idx[0],idx[1],idx[2],idx[3],Tet[iTet][4]);  
  }

  //--- Write prisms
  GmfSetKwd(OutMsh, GmfPrisms, Msh->NbrPri);
  for (i=1; i<=Msh->NbrPri; ++i) {
    for (j=0; j<6; ++j) {
      idx[j] = (long long)(Msh->Pri[i][j]);
    }
    GmfSetLin(OutMsh, GmfPrisms,idx[0],idx[1],idx[2],idx[3],idx[4],idx[5],Msh->Pri[i][6]);  
  }

  //--- Write pyr
  GmfSetKwd(OutMsh, GmfPyramids, Msh->NbrPyr);
  for (i=1; i<=Msh->NbrPyr; ++i) {
    for (j=0; j<5; ++j) {
      idx[j] = (long long)(Msh->Pyr[i][j]);
    }
    GmfSetLin(OutMsh, GmfPyramids,idx[0],idx[1],idx[2],idx[3],idx[4],Msh->Pyr[i][5]);  
  }

  //--- Write Edges
  GmfSetKwd(OutMsh, GmfEdges, NbrEfr);
  for (iEfr=1; iEfr<=NbrEfr; ++iEfr) {
    for (i=0; i<2; ++i) {
      idx[i] = (long long)(Efr[iEfr][i]);
    }
    GmfSetLin(OutMsh, GmfEdges,idx[0],idx[1],Efr[iEfr][2]);  
  }
  
  //--- close mesh file
  if ( !GmfCloseMesh(OutMsh) ) {
    printf("  ## ERROR: Cannot close mesh file %s ! \n",OutFil);
		return 0;
  }
  
  return 1;
}


int WriteGMFSolution(char *SolNam, double *Sol, int SolSiz, int NbrVer, int Dim, int NbrFld, int* FldTab)
{
  int       OutSol, iVer;
	double   *dbl=NULL;
	
	if ( !Sol ) {
		printf("  ## ERROR WriteGMFSolution : Sol not allocated.\n");
		return 0;	
	}
	
	if ( SolSiz < 1 ) {
		printf("  ## ERROR WriteGMFSolution : SolSiz < 1.\n");
		return 0;
	}
	
	//--- Open solution file
	if ( !(OutSol = GmfOpenMesh(SolNam, GmfWrite, GmfDouble, Dim)) ) {
    fprintf(stderr,"  ## ERROR: Cannot open solution file %s ! \n",SolNam);
    exit(1);
  }
  printf("  %%%% %s OPENED (WRITE)\n",SolNam);
	
	printf("OUTSOL=%d\n", OutSol);
	
  GmfSetKwd(OutSol, GmfSolAtVertices, NbrVer, NbrFld, FldTab);
	
  for (iVer=1; iVer<=NbrVer; ++iVer) {
		dbl = &Sol[iVer*SolSiz];
    GmfSetLin(OutSol, GmfSolAtVertices, dbl);
  }
		
	if ( !GmfCloseMesh(OutSol) ) {
	  printf("  ## ERROR: Cannot close solution file %s ! \n",SolNam);
		return 0;
	}
	
	return 1;
}


/*
	Interface to function WriteGMFSolution()
*/
int WriteGMFSolutionItf(char *SolNam, Mesh *Msh)
{
	double *Sol      = Msh->Sol;
	int     SolSiz   = Msh->SolSiz;
	int     NbrVer   = Msh->NbrVer;
	int     Dim      = Msh->Dim; 
	int     NbrFld   = Msh->NbrFld; 
	int    *FldTab   = Msh->FldTab; 
	
	return WriteGMFSolution(SolNam, Sol, SolSiz, NbrVer, Dim, NbrFld, FldTab);
}
