#include "meshutils.h"
#include "Python.h"




int py_ProjectNozzleWall3D( char *MshNam,  
 PyObject *pyRefUp,  PyObject *pyRefDown,
 PyObject *pyKnots_center, PyObject *pyCoefs_center,
 PyObject *pyKnots_r1, PyObject *pyCoefs_r1,
 PyObject *pyKnots_r2, PyObject *pyCoefs_r2, 
 char *OutNam ) 
{
	
	int i;
	int *RefUp = NULL, *RefDown = NULL;
	int sizRefUp=0, sizRefDown=0;
	
	double  *Knots_center=NULL, *Coefs_center=NULL;
	double  *Knots_r1=NULL, *Coefs_r1=NULL;
	double  *Knots_r2=NULL, *Coefs_r2=NULL;
	
	int nKnots_center=0, nCoefs_center=0;
	int nKnots_r1=0, nCoefs_r1=0;
	int nKnots_r2=0, nCoefs_r2=0;
	
	CadNozzle *CadNoz = NULL;
					
	//--- Get nozzle patch references to re-project
	
	if ( PyList_Check(pyRefUp) )
  {
			sizRefUp = PyList_Size(pyRefUp);
      RefUp = (int*) malloc( sizRefUp * sizeof(int));
			for (i=0; i<sizRefUp; i++)
      {
       	PyObject *oo = PyList_GetItem(pyRefUp,i);
       	if ( PyInt_Check(oo) )
       	{
					RefUp[i] = (int) PyInt_AS_LONG(oo);
       	}
      }
  }

	if ( PyList_Check(pyRefDown) )
  {
			sizRefDown = PyList_Size(pyRefDown);
      RefDown = (int*) malloc( sizRefDown * sizeof(int));
			for (i=0; i<sizRefDown; i++)
      {
       	PyObject *oo = PyList_GetItem(pyRefDown,i);
       	if ( PyInt_Check(oo) )
       	{
					RefDown[i] = (int) PyInt_AS_LONG(oo);
       	}
      }
  }

	//---
	//--- Get coefs and knots of the 3 b-splines defining the nozzle
	//---
	
	//--- center knots and coefs
	
	if ( PyList_Check(pyKnots_center) )
  {
			nKnots_center = PyList_Size(pyKnots_center);
						
      Knots_center = (double*) malloc( nKnots_center * sizeof(double));
			for (i=0; i<nKnots_center; i++)
      {
       	PyObject *oo = PyList_GetItem(pyKnots_center,i);
       	if ( PyFloat_Check(oo) )
       	{
					Knots_center[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }

	if ( PyList_Check(pyCoefs_center) )
  {
			nCoefs_center = PyList_Size(pyCoefs_center);
      Coefs_center = (double*) malloc( nCoefs_center * sizeof(double));
			for (i=0; i<nCoefs_center; i++)
      {
       	PyObject *oo = PyList_GetItem(pyCoefs_center,i);
       	if ( PyFloat_Check(oo) )
       	{
					Coefs_center[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }

	//--- r1 coefs and knots

	if ( PyList_Check(pyKnots_r1) )
  {
			nKnots_r1 = PyList_Size(pyKnots_r1);
      Knots_r1 = (double*) malloc( nKnots_r1 * sizeof(double));
			for (i=0; i<nKnots_r1; i++)
      {
       	PyObject *oo = PyList_GetItem(pyKnots_r1,i);
       	if ( PyFloat_Check(oo) )
       	{
					Knots_r1[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }

	if ( PyList_Check(pyCoefs_r1) )
  {
			nCoefs_r1 = PyList_Size(pyCoefs_r1);
      Coefs_r1 = (double*) malloc( nCoefs_r1 * sizeof(double));
			for (i=0; i<nCoefs_r1; i++)
      {
       	PyObject *oo = PyList_GetItem(pyCoefs_r1,i);
       	if ( PyFloat_Check(oo) )
       	{
					Coefs_r1[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }

	//--- r2 coefs and knots

	if ( PyList_Check(pyKnots_r2) )
  {
			nKnots_r2 = PyList_Size(pyKnots_r2);
      Knots_r2 = (double*) malloc( nKnots_r2 * sizeof(double));
			for (i=0; i<nKnots_r2; i++)
      {
       	PyObject *oo = PyList_GetItem(pyKnots_r2,i);
       	if ( PyFloat_Check(oo) )
       	{
					Knots_r2[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }

	if ( PyList_Check(pyCoefs_r2) )
  {
			nCoefs_r2 = PyList_Size(pyCoefs_r2);
      Coefs_r2 = (double*) malloc( nCoefs_r2 * sizeof(double));
			for (i=0; i<nCoefs_r2; i++)
      {
       	PyObject *oo = PyList_GetItem(pyCoefs_r2,i);
       	if ( PyFloat_Check(oo) )
       	{
					Coefs_r2[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }

	//--- Alloc nozzle CAD data structure
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,"");
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	//--- Open mesh/solution file
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	//--- Setup nozzle CAD
	
	int SizCad[MaxKwdNozzleSize];
	memset(SizCad, 0, sizeof(int)*MaxKwdNozzleSize);
	
	SizCad[KwdnCoefs_center] = nCoefs_center;
	SizCad[KwdnKnots_center] = nKnots_center;
	SizCad[KwdnCoefs_r1]     = nCoefs_r1;
	SizCad[KwdnKnots_r1]     = nKnots_r1;
	SizCad[KwdnCoefs_r2]     = nCoefs_r2;
	SizCad[KwdnKnots_r2]     = nKnots_r2;

	CadNoz = AllocCadNozzle (SizCad);
	
	SetCadBspline (CadNoz->Bsp_center, Knots_center,  nKnots_center, Coefs_center, nCoefs_center);
	SetCadBspline (CadNoz->Bsp_r1, Knots_r1, nKnots_r1, Coefs_r1, nCoefs_r1);
	SetCadBspline (CadNoz->Bsp_r2, Knots_r2, nKnots_r2, Coefs_r2, nCoefs_r2);
	
	WriteCadBspline ("centerline", CadNoz->Bsp_center);
	WriteCadBspline ("r1", CadNoz->Bsp_r1);
	WriteCadBspline ("r2", CadNoz->Bsp_r2);
	
	//--- Project
	
	NozzleWallProjection (mshopt, Msh, CadNoz,  RefUp[0], RefDown[0], OutNam);
	
	//--- Free memory
	
	if ( Msh )
 		FreeMesh(Msh);
	
	if ( RefUp )
		free(RefUp);
	
	if ( RefDown )
		free(RefDown);
	
}

int py_ConvertGMFToSU2( char *MshNam, char *SolNam, char *OutNam ) 
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
}

int py_ConvertSU2toGMSH( char *MshNam, char *SolNam, char *OutNam ) 
{
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	mshopt->clean = 0; // remove unconnected vertices
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	return ConvertSU2ToGMSH (mshopt);
	
	
	return 1;
}


int py_MeshPrepro2D( char *InpNam, char *OutNam ) 
{
	
	Options *mshopt = AllocOptions();
		
	char BasNam[1024];
  char *ptr = NULL;
	
	// --- Get BasNam
	
  strcpy(BasNam,OutNam);
  
  ptr = strstr(BasNam,".su2");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
  
	printf("BASNAM = %s\n", BasNam);

	strcpy(mshopt->OutNam,BasNam);	
	
	//strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,InpNam);
	
	mshopt->clean = 1; // remove unconnected vertices
	mshopt->Dim   = 2; // force dim to be 2
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	return ConvertGMFtoSU2Sol (mshopt);
}


int py_BSplineGeo3 (PyObject *pyknots, PyObject *pycoefs, PyObject *pyx, PyObject *pyy, int nx)
{	
	int i, k, c;
	double len, hx;
	int size_knots = 0;
	int size_coefs = 0;
	int size_x     = 0;
	int size_y     = 0;
	int size_dydx  = 0;
	
	double *knots = NULL;
	double *coefs = NULL;
	double *x     = NULL;
	double *y     = NULL;
	double *dydx  = NULL;
		
	//--- Knots
	
  if ( PyList_Check(pyknots) )
  {
      size_knots = PyList_Size(pyknots);
      knots = malloc( size_knots * sizeof(double));
			
			for (i=0; i<size_knots; i++)
      {
       	PyObject *oo = PyList_GetItem(pyknots,i);
       	if ( PyFloat_Check(oo) )
       	{
					knots[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }
	
	//--- Coefs

  if ( PyList_Check(pycoefs) )
  {
      size_coefs = PyList_Size(pycoefs);
      coefs = malloc( size_coefs * sizeof(double));
			
			for (i=0; i<size_coefs; i++)
      {
       	PyObject *oo = PyList_GetItem(pycoefs,i);
       	if ( PyFloat_Check(oo) )
       	{
					coefs[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }
	
	//--- Call function
	
	x    = (double*)malloc(nx*sizeof(double));
	y    = (double*)malloc(nx*sizeof(double));
	dydx = (double*)malloc(nx*sizeof(double));
	
	k = size_knots;
	c = size_coefs/2;
	
	len = coefs[c-1];
	hx = len/(double)(nx-1);
	
	x[0] = 0.;
	for (i=1; i<nx; i++) {
		x[i] = x[i-1]+hx;
		//printf("x[%d] = %lf\n", i, x[i]);
	}
	
	bSplineGeo3(knots, coefs, x, y, dydx, nx, k, c);
	
	for (i=0; i<nx; i++){
		PyList_Append(pyx, PyFloat_FromDouble(x[i]));
	}
	
	for (i=0; i<nx; i++){
		PyList_Append(pyy, PyFloat_FromDouble(y[i]));
	}
	
	//--- Free memory
	
	if (x)
		free(x);
	if (y)
		free(y);
	if (dydx)
		free(dydx);
}



int py_BSplineGeo3LowF (PyObject *pyknots, PyObject *pycoefs, PyObject *pyx, PyObject *pyy, PyObject *pydydx)
{	
	int i, k, c, nx=0;
	double len, hx;
	int size_knots = 0;
	int size_coefs = 0;
	int size_x     = 0;
	int size_y     = 0;
	int size_dydx  = 0;
	
	double *knots = NULL;
	double *coefs = NULL;
	double *x     = NULL;
	double *y     = NULL;
	double *dydx  = NULL;
	
	//--- Knots
		
  if ( PyList_Check(pyknots) )
  {
      size_knots = PyList_Size(pyknots);
      knots = malloc( size_knots * sizeof(double));
			
			for (i=0; i<size_knots; i++)
      {
       	PyObject *oo = PyList_GetItem(pyknots,i);
       	if ( PyFloat_Check(oo) )
       	{
					knots[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }
	
	//--- Coefs

  if ( PyList_Check(pycoefs) )
  {
      size_coefs = PyList_Size(pycoefs);
      coefs = malloc( size_coefs * sizeof(double));
			
			for (i=0; i<size_coefs; i++)
      {
       	PyObject *oo = PyList_GetItem(pycoefs,i);
       	if ( PyFloat_Check(oo) )
       	{
					coefs[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }
	
	//--- x
	
	if ( PyList_Check(pyx) )
  {
      nx = PyList_Size(pyx);
			x  = (double*)malloc(nx*sizeof(double));
			
			for (i=0; i<nx; i++)
      {
       	PyObject *oo = PyList_GetItem(pyx,i);
       	if ( PyFloat_Check(oo) )
       	{
					x[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }
	
	//--- Call function
	
	y    = (double*)malloc(nx*sizeof(double));
	dydx = (double*)malloc(nx*sizeof(double));
	
	k = size_knots;
	c = size_coefs/2;
	
	len = coefs[c-1];
	
	bSplineGeo3(knots, coefs, x, y, dydx, nx, k, c);
	
	//for (i=0; i<nx; i++){
	//	PyList_Append(pyx, PyFloat_FromDouble(x[i]));
	//}
	
	for (i=0; i<nx; i++){
		PyList_Append(pyy, PyFloat_FromDouble(y[i]));
	}
	
	for (i=0; i<nx; i++){
		PyList_Append(pydydx, PyFloat_FromDouble(dydx[i]));
	}
	
	//--- Free memory
	
	if (x)
		free(x);
	if (y)
		free(y);
	if (dydx)
		free(dydx);
}




void py_ReadMesh (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pySol)
{
	int i, j, d;
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	//--- Open mesh/solution file
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	//PrintMeshInfo (Msh);
  //if ( !Msh->Sol ) {
	//	printf("  ## ERROR SolutionExtraction : A solution must be provided.\n");
	//	return;
	//}
	
	for (i=1; i<=Msh->NbrVer; i++){
		for (d=0; d<3; d++)
			PyList_Append(pyVer, PyFloat_FromDouble(Msh->Ver[i][d]));
	}
	
	for (i=1; i<=Msh->NbrTri; i++){
		for (j=0; j<4; j++)
			PyList_Append(pyTri, PyFloat_FromDouble(Msh->Tri[i][j]));
	}
	
	for (i=1; i<=Msh->NbrTet; i++){
		for (j=0; j<4; j++)
			PyList_Append(pyTet, PyFloat_FromDouble(Msh->Tet[i][j]));
	}
	
	for (i=1; i<=Msh->NbrEfr; i++){
		for (j=0; j<2; j++)
			PyList_Append(pyEdg, PyFloat_FromDouble(Msh->Efr[i][j]));
	}
	
	if ( Msh->Sol ) {
		
		//--- Output solution
		int iVer;
		for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
			for (i=0; i<Msh->SolSiz; i++) {
				PyList_Append(pySol, PyFloat_FromDouble(Msh->Sol[iVer*Msh->SolSiz+i]));
			}
		}
		
	}
	
	if ( Msh )
 		FreeMesh(Msh);
	
}




void py_ExtractAlongLine (char *MshNam, char *SolNam, PyObject *pyBox,  PyObject *pyResult, PyObject *PyInfo, PyObject *pyHeader)
{
	int i;
	double box[4] = {0.67,0.67,0,0.3048};
	
	//--- Get box
	
	if ( PyList_Check(pyBox) )
  {
		
			for (i=0; i<4; i++)
      {
       	PyObject *oo = PyList_GetItem(pyBox,i);
       	if ( PyFloat_Check(oo) )
       	{
					box[i] = (double) PyFloat_AS_DOUBLE(oo);
       	}
      }
  }
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	
	//--- Open mesh/solution file
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	//PrintMeshInfo (Msh);
	
	if ( !Msh->Sol ) {
		printf("  ## ERROR SolutionExtraction : A solution must be provided.\n");
		return;
	}
	
	int NbrRes=0, Siz=0;
	double *result = ExtractAlongLine(mshopt,Msh, box, &NbrRes, &Siz);
	
	for (i=0; i<NbrRes*Siz; i++){
		PyList_Append(pyResult, PyFloat_FromDouble(result[i]));
	}
	
	
	PyList_Append(PyInfo, PyInt_FromLong(NbrRes));
	PyList_Append(PyInfo, PyInt_FromLong(Siz));
	
	
	for (i=0; i<Msh->SolSiz; i++){
		PyList_Append(pyHeader, PyString_FromString(Msh->SolTag[i]));
	}
	
	
	if (result)
		free(result);
	
	if ( Msh )
 		FreeMesh(Msh);
	
}


//void py_NozzleWallProjection (char *MshNam, char *SolNam, PyObject *pyMeshMotion,  PyObject *pyDV)
//{
//	int i;
//	
//	int *Ref = NULL;
//	int size_Ref = 0;
//	
//	////--- Get box
//	//
//	//if ( PyList_Check(pyRefs) )
//  //{
//	//
//	//		size_Ref = PyList_Size(pyRefs);
//  //    Ref = (int*) malloc( size_Ref * sizeof(int));
//	//		
//	//		for (i=0; i<size_Ref; i++)
//  //    {
//  //     	PyObject *oo = PyList_GetItem(pyRefs,i);
//  //     	if ( PyInt_Check(oo) )
//  //     	{
//	//				Ref[i] = (int) PyInt_AS_LONG(oo);
//  //     	}
//  //    }
//  //}
//  //
//	Options *mshopt = AllocOptions();
//	//
//	strcpy(mshopt->InpNam,MshNam);
//	strcpy(mshopt->SolNam,SolNam);
//	//
//	//
//	//--- Open mesh/solution file
//	Mesh *Msh = NULL;
//	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
//	
//	PrintMeshInfo(Msh);
//	
//	int RefUp = 9;
//	int RefDown = 10;
//	
//	NozzleWallProjection (mshopt,Msh, &RefUp, &RefDown);
//	
//	////double * ExtractSolutionAtRef (Options *mshopt, Mesh *Msh, int *Ref, int NbrRef,  int *NbrRes, int *Siz)
//	//int NbrRes=0, Siz=0;
//	//double *result = ExtractSolutionAtRef(mshopt,Msh, Ref, size_Ref,  &NbrRes, &Siz);
//	//
//	//for (i=0; i<NbrRes*Siz; i++){
//	//	PyList_Append(pyResult, PyFloat_FromDouble(result[i]));
//	//}
//	//
//	//PyList_Append(PyInfo, PyInt_FromLong(NbrRes));
//	//PyList_Append(PyInfo, PyInt_FromLong(Siz));
//	//
//	//for (i=0; i<Msh->SolSiz; i++){
//	//	PyList_Append(pyHeader, PyString_FromString(Msh->SolTag[i]));
//	//}
//	//
//	
//	
//	if ( Msh )
// 		FreeMesh(Msh);
//	
//}



void py_ExtractAtRef (char *MshNam, char *SolNam, PyObject *pyRefs,  PyObject *pyResult, PyObject *PyInfo, PyObject *pyHeader)
{
	int i;
	
	int *Ref = NULL;
	int size_Ref = 0;
	
	//--- Get box
	
	if ( PyList_Check(pyRefs) )
  {
	
			size_Ref = PyList_Size(pyRefs);
      Ref = (int*) malloc( size_Ref * sizeof(int));
			
			for (i=0; i<size_Ref; i++)
      {
       	PyObject *oo = PyList_GetItem(pyRefs,i);
       	if ( PyInt_Check(oo) )
       	{
					Ref[i] = (int) PyInt_AS_LONG(oo);
       	}
      }
  }

	Options *mshopt = AllocOptions();
	//
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	//
	//
	//--- Open mesh/solution file
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	//PrintMeshInfo (Msh);
	
	if ( !Msh->Sol ) {
		printf("  ## ERROR SolutionExtraction : A solution must be provided.\n");
		return;
	}
	
	//double * ExtractSolutionAtRef (Options *mshopt, Mesh *Msh, int *Ref, int NbrRef,  int *NbrRes, int *Siz)
	int NbrRes=0, Siz=0;
	double *result = ExtractSolutionAtRef(mshopt,Msh, Ref, size_Ref,  &NbrRes, &Siz);
	
	for (i=0; i<NbrRes*Siz; i++){
		PyList_Append(pyResult, PyFloat_FromDouble(result[i]));
	}
	
	PyList_Append(PyInfo, PyInt_FromLong(NbrRes));
	PyList_Append(PyInfo, PyInt_FromLong(Siz));
	
	for (i=0; i<Msh->SolSiz; i++){
		PyList_Append(pyHeader, PyString_FromString(Msh->SolTag[i]));
	}
	
	
	if (result)
		free(result);
	
	if ( Msh )
 		FreeMesh(Msh);
	
}



void py_Extract_Vertices_Ref (char *MshNam, PyObject * pyRefs , PyObject * PyCrd_Out, PyObject * PyVid_Out, PyObject * PyRef_Tab)
{
	int i, iRef, j, idx;
	
	int *Ref = NULL;
	int size_Ref = 0;
	
	//--- Get box
	
	if ( PyList_Check(pyRefs) )
  {
	
			size_Ref = PyList_Size(pyRefs);
      Ref = (int*) malloc( size_Ref * sizeof(int));
			
			for (i=0; i<size_Ref; i++)
      {
       	PyObject *oo = PyList_GetItem(pyRefs,i);
       	if ( PyInt_Check(oo) )
       	{
					Ref[i] = (int) PyInt_AS_LONG(oo);
       	}
      }
  }

	Options *mshopt = AllocOptions();
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,"");

	//--- Open mesh/solution file
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, "");
	
	double * Crd_Out = (double *) malloc(sizeof(double)*3*Msh->NbrVer);
	int    * Vid_Out = (int * ) malloc(sizeof(int)*Msh->NbrVer);
	
	int Nbv=0;
	
	for (iRef=0; iRef<size_Ref; iRef++)
  {
		
		Extract_Vertices_Ref (mshopt, Msh, Ref[iRef], Crd_Out, Vid_Out, &Nbv);
		
		for (i=0; i<Nbv; i++) {
			idx = 3*i;
			for (j=0; j<3; j++)
				PyList_Append(PyCrd_Out, PyFloat_FromDouble(Crd_Out[idx+j]));
			PyList_Append(PyVid_Out, PyInt_FromLong(Vid_Out[i]));
		}
		
		
		PyList_Append(PyRef_Tab, PyInt_FromLong(Ref[iRef]));
		PyList_Append(PyRef_Tab, PyInt_FromLong(Nbv));
		
  }

	//
	////double * ExtractSolutionAtRef (Options *mshopt, Mesh *Msh, int *Ref, int NbrRef,  int *NbrRes, int *Siz)
	//int NbrRes=0, Siz=0;
	//double *result = ExtractSolutionAtRef(mshopt,Msh, Ref, size_Ref,  &NbrRes, &Siz);
	//
	//for (i=0; i<NbrRes*Siz; i++){
	//	PyList_Append(pyResult, PyFloat_FromDouble(result[i]));
	//}
	//
	//PyList_Append(PyInfo, PyInt_FromLong(NbrRes));
	//PyList_Append(PyInfo, PyInt_FromLong(Siz));
	//
	//for (i=0; i<Msh->SolSiz; i++){
	//	PyList_Append(pyHeader, PyString_FromString(Msh->SolTag[i]));
	//}
	
	
	if (Crd_Out)
		free(Crd_Out);
		
	if (Vid_Out)
		free(Vid_Out);
	
	if ( Msh )
 		FreeMesh(Msh);
	
}


