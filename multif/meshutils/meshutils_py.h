
//int hello_meshutils(int toto);

int py_ConvertGMFToSU2( char *MshNam, char *SolNam, char *OutNam ) ;


int py_ConvertSU2toGMSH( char *MshNam, char *SolNam, char *OutNam ) ;


/*
	Return nx points (pyx,pyy) of the bspline defined by knots and coefs
*/
int py_BSplineGeo3 (PyObject *pyknots, PyObject *pycoefs, PyObject *pyx, PyObject *pyy, int nx);

/*
	Mesh preprocessing: 
		remove unconnected vertices
		force 2D
		convert to su2 format
*/
int py_MeshPrepro2D (char *InpNam, char *OutNam);


/*
	Compute bspline for the low f level
	Similar to py_BSplineGeo3, only different i/o
*/
int py_BSplineGeo3LowF (PyObject *pyknots, PyObject *pycoefs, PyObject *pyx, PyObject *pyy, PyObject *pydydx);


void py_ExtractAlongLine (char *MshNam, char *SolNam, PyObject *pyBox,  PyObject *pyResult, PyObject *PyInfo, PyObject *pyHeader);

void py_ExtractAtRef (char *MshNam, char *SolNam, PyObject *pyRefs,  PyObject *pyResult, PyObject *PyInfo, PyObject *pyHeader);


void py_ReadMesh (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg);