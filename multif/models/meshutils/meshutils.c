#include "meshutils.h"

/*
Victorien Menier Feb 2016
*/

int main(int argc, char *argv[])
{
		
	Usage();
	
	Options *mshopt = AllocOptions();
	ParseCommandLine(mshopt, argc, argv);
	
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	} 
	
	PrintOptions(mshopt);
	
	
	
	switch ( mshopt->Mod ) {
		
		case 1:
			return OutputMach (mshopt);
		break;
		
		case 2:
			return ConvertSU2SolToGMF (mshopt);
		break;
		
		case 3:
			return ConvertGMFtoSU2Sol (mshopt);
		break;
		
		case 4:
			return Extraction(mshopt);
		break;
		
		//case 5:
		//	Geo2Egads ();
		//break;
		
		
		//case 6:
		//	CleanMesh ();
		//break;
		
		//case 7:
		//	return ProjectNozzleWall (mshopt);
		//break;
		
		//case 8:
		//	return ConvertGMFtoSegMesh (mshopt);
		//break;
		
		default:
		printf("  ## ERROR : Unknown module number.\n");
	}
	
	return 0;
}