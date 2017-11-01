#include "amgio.h"
#include "Python.h"


int py_ConvertSU2toInria( char *MshNam, char *SolNam, char *OutNam ) 
{
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	mshopt->clean = 0; // remove unconnected vertices
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	return ConvertSU2SolToGMF (mshopt);
	
	
	return 1;
}


int py_ConvertInriatoSU2( char *MshNam, char *SolNam, char *OutNam ) 
{
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	mshopt->clean = 0; // remove unconnected vertices
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	return ConvertGMFtoSU2Sol (mshopt);
	
	
	return 1;
}


int py_SplitSolution(char *SolNam, int dim, char *prefix, char *adap_sensor)
{
	
	int SizMsh[GmfMaxSizMsh+1];
	memset(SizMsh,0,sizeof(int)*(GmfMaxSizMsh+1));
	
	Mesh *Msh = AllocMesh(SizMsh);
	
	Msh->NbrVer = GetSU2SolSize(SolNam);
	
	LoadSU2Solution(SolNam, Msh);
	
	Msh->Dim = dim;
	SplitSolution(Msh, prefix, adap_sensor);
	
}