#include "meshutils.h"
#include "python.h"

int py_ConvertGMFToSU2( char *MshNam, char *SolNam, char *OutNam ) 
{
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	mshopt->clean = 1; // remove unconnected vertices
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	return ConvertGMFtoSU2Sol (mshopt);
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


