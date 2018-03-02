
//int hello_meshutils(int toto);

int py_ProjectNozzleWall3D( char *MshNam,  
 PyObject *pyRefUp,  PyObject *pyRefDown,
 PyObject *pyKnots_center, PyObject *pyCoefs_center,
 PyObject *pyKnots_r1, PyObject *pyCoefs_r1,
 PyObject *pyKnots_r2, PyObject *pyCoefs_r2, 
 char *OutNam );
 
 int py_ProjectNozzleWall3D_DV( char *MshNam,  
  PyObject *pyRefUp,  PyObject *pyRefDown,
  PyObject *pyKnots_center_from, PyObject *pyCoefs_center_from,
  PyObject *pyKnots_r1_from, PyObject *pyCoefs_r1_from,
  PyObject *pyKnots_r2_from, PyObject *pyCoefs_r2_from, 
  PyObject *pyKnots_center_to, PyObject *pyCoefs_center_to,
  PyObject *pyKnots_r1_to, PyObject *pyCoefs_r1_to,
  PyObject *pyKnots_r2_to, PyObject *pyCoefs_r2_to, 
  char *OutNam ) ;

//void py_NozzleWallProjection (char *MshNam, char *SolNam, PyObject *pyMeshMotion,  PyObject *pyDV);

int py_ConvertGMFToSU2( char *MshNam, char *SolNam, char *OutNam );

int py_ConvertSU2toGMSH( char *MshNam, char *SolNam, char *OutNam );

int py_ConvertSU2toGMF( char *MshNam, char *SolNam, char *OutNam );

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

int py_PiecewiseLinear (PyObject *pyxnodes, PyObject *pyynodes, PyObject *pyx, PyObject *pyy, PyObject *pydydx);

void py_ExtractAlongLine (char *MshNam, char *SolNam, PyObject *pyBox,  PyObject *pyResult, PyObject *PyInfo, PyObject *pyHeader);

void py_ExtractAtRef (char *MshNam, char *SolNam, PyObject *pyRefs,  PyObject *pyResult, PyObject *PyInfo, PyObject *pyHeader);

void py_Extract_Vertices_Ref (char *MshNam, PyObject * pyRefs , PyObject * PyCrd_Out, PyObject * PyVid_Out, PyObject * PyRef_Tab);

void py_ReadMesh (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pySol);

void py_ReadMesh2 (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pySol);

void py_WriteMesh(char *OutNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pySol);

void py_ExtractSurfacePatches (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pySol, PyObject *pyRefs);
