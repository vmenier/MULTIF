#include "meshutils.h"

/*
Victorien Menier Aug 2016
*/

void WriteGMSHMesh(char *nam, Mesh *Msh)
{
  int       i, s, iElt;
  int       iVer,iTri,iEfr, iTet, NbrElt=0;
  char      OutNam[512];
  
  int Dim = Msh->Dim;
	
	int2 *BdrTag=NULL;
	int NbrBdr, NbrTag, start, iTag, cpt;
	
	FILE *OutFil=NULL;
	
	sprintf(OutNam, "%s.msh", nam);
 
	OutFil = fopen(OutNam, "wb");
	
	if ( !OutFil ) {
		printf("  ## ERROR Write GMSH: Can't open %s\n", OutNam);
	}
 	
  printf("  %%%% %s OPENED (WRITE)\n",OutNam);
  
	fprintf(OutFil, "$MeshFormat    \n");
	fprintf(OutFil, "2.0 0 8        \n");
	fprintf(OutFil, "$EndMeshFormat \n");
	
	
	//--- Write vertices
  fprintf(OutFil, "$Nodes\n%d\n", Msh->NbrVer);

	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
  	fprintf(OutFil, "%d  %le %le %le \n", iVer, Msh->Ver[iVer][0], Msh->Ver[iVer][1], Msh->Ver[iVer][2]);
	}
	fprintf(OutFil, "$EndNodes\n");

	//--- Write elements

	if ( Msh->Dim == 2 ) {
		NbrElt = Msh->NbrTri+Msh->NbrEfr;
	}
	else {
		NbrElt = Msh->NbrTet;
	}

  fprintf(OutFil, "$Elements\n%d\n", NbrElt);
	
	iElt = 0;
	
	//--- Write triangles
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		iElt++;
    fprintf(OutFil, "%d %d 2 0 %d ", iElt, GMSH_TRIANGLE, Msh->Tri[iTri][3]); 
		for (i=0; i<3; ++i) {
      fprintf(OutFil, "%d ",Msh->Tri[iTri][i]);
    }
    fprintf(OutFil, "\n"); 
  }

	//--- Write edges
	for (iEfr=1; iEfr<=Msh->NbrEfr; iEfr++) {
		iElt++;
    fprintf(OutFil, "%d %d 2 0 %d ", iElt, GMSH_EDGE, Msh->Efr[iEfr][2]); 
		for (i=0; i<2; ++i) {
      fprintf(OutFil, "%d ",Msh->Efr[iEfr][i]);
    }
    fprintf(OutFil, "\n"); 
  }
	
	fprintf(OutFil, "$EndElements\n");
	
  //--- close mesh file
	if ( OutFil )
		fclose(OutFil);
	
  return;
}

