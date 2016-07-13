#include "meshutils.h"

/*
Victorien Menier March 2016
*/


int SolutionExtraction(Options *mshopt,Mesh *Msh)
{
	int iVer, i, idxVer;
	double x, y;
	
	double eps = 1e-15;
	
	FILE *FilHdl = fopen ("yextract.dat", "wb");
	
	if ( !FilHdl ) {
		printf("  ## ERROR SolutionExtraction: Could not open file.\n");
		return 0;
	}
	
	//--- Write header
	fprintf(FilHdl, "# x y ");
	for (i=0; i<Msh->SolSiz; i++) {
		fprintf(FilHdl, "%s ", Msh->SolTag[i]);
	}
	fprintf(FilHdl, "\n");
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
		x = Msh->Ver[iVer][0];
		y = Msh->Ver[iVer][1];
		
		if ( fabs(y) > eps ) continue;
		
		fprintf(FilHdl, "%lf %lf ", x, y);
		
		idxVer = iVer*Msh->SolSiz;
		
		for (i=0; i<Msh->SolSiz; i++) {
			fprintf(FilHdl, "%lf ", Msh->Sol[idxVer+i]);
		}
		
		fprintf(FilHdl, "\n");
		
		
	}
	
	if ( FilHdl )
		fclose(FilHdl);
	
	return 1;
	
}