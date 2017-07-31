#include "meshutils.h"


int P2SurfaceMesh()
{

	//--- Open mesh/solution file
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	PrintMeshInfo (Msh);
	
	
}