#include "meshutils.h"

/*
Victorien Menier Aug 2016
*/

void WriteFLUENTMesh(char *nam, Mesh *Msh)
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
		printf("  ## ERROR Write FLUENT: Can't open %s\n", OutNam);
	}
	
	//--- Write header
	
	fprintf(OutFil, "(1 \"File generated using MULTIF\")\n");
	
	//--- Dimension
	
	fprintf(OutFil, "(0 \"Dimension\")\n");
	fprintf(OutFil, "(2 %d)\n", Msh->Dim);
	
	//--- Cells
	
	fprintf(OutFil, "(0 \"Cells\")\n");
	
	//--- Vertices
	
	
	
	
	if ( OutFil )
		fclose(OutFil);
	
}