//--- GMFio.c
int AddGMFMeshSize (char *MshNam, int *SizMsh);
int LoadGMFMesh (char *MshNam, VMesh *Msh);
int LoadGMFSolution(char *SolNam, VMesh *Msh);
int WriteGMFMesh(char *nam, VMesh *Msh, int OptBin);
int WriteGMFSolution(char *SolNam, double *Sol, int SolSiz, int NbrVer, int Dim, int NbrFld, int* FldTab);
int WriteGMFSolutionItf(char *SolNam, VMesh *Msh);

//---- mesh.c
VMesh* AllocMesh (int * SizMsh);
int   cmp_int2(const void *a, const void *b);
int   FreeMesh (VMesh *Msh);
int   GetMeshSize (char *MshNam, int *SizMsh);
VMesh *SetupMeshAndSolution (char *MshNam, char *SolNam);
int   RemoveUnconnectedVertices (VMesh *Msh);
void  AddEdge(VMesh *Msh, int idx, int *is, int ref);
void  AddHexahedron(VMesh *Msh, int idx, int *is, int ref);
void  AddQuadrilateral(VMesh *Msh, int idx, int *is, int ref);
void  AddTetrahedron(VMesh *Msh, int idx, int *is, int ref);
void  AddPyramid(VMesh *Msh, int idx, int *is, int ref);
void  AddPrism(VMesh *Msh, int idx, int *is, int ref);
void  AddTriangle(VMesh *Msh, int idxTri, int *is, int ref);
void  AddVertex(VMesh *Msh, int idxVer, double *Crd);
int   imin(int n, int *idx);
void  PrintMeshInfo (VMesh *Msh);
void  switchHexIdx(int *idx, int *swi);
void  switchQuaIdx(int *idx, int *swi);
void  switchTetIdx(int *idx, int *swi);
void  switchTriIdx(int *idx, int *swi);
int GetInputFileType (char *FilNam) ;
int  Str2Lower(char *buff);
void StrRemoveChars (char* str, char c);


//--- SU2io.c
int  AddSU2MeshSize (char *FilNam, int *SizMsh) ;
int  GetSU2KeywordValue (FILE *FilHdl, char *Kwd);
int  GetSU2KeywordValueStr (FILE *FilHdl, char *Kwd, char *StrVal);
int  LoadSU2Elements(FILE *FilHdl, VMesh *Msh);
int  LoadSU2Mesh (char *FilNam, VMesh *Msh);
int  LoadSU2Solution(char *SolNam, VMesh *Msh);
int  LoadSU2Vertices (FILE *FilHdl, VMesh *Msh);
void WriteSU2Mesh(char *nam, VMesh *Msh);
int  WriteSU2Solution (char *SolNam, VMesh *Msh, double *Sol, int NbrVer, int SolSiz, char SolTag[100][256]);

//--- GMFio.c
int AddGMFMeshSize (char *MshNam, int *SizMsh);
int LoadGMFMesh (char *MshNam, VMesh *Msh);
int LoadGMFSolution(char *SolNam, VMesh *Msh);
int WriteGMFMesh(char *nam, VMesh *Msh, int OptBin);
int WriteGMFSolution(char *SolNam, double *Sol, int SolSiz, int NbrVer, int Dim, int NbrFld, int* FldTab);
int WriteGMFSolutionItf(char *SolNam, VMesh *Msh);

//// --- mshint.c
//static void endcod() ;
//static void stats(pMesh mesh1,pSol sol1,pMesh mesh2,pSol sol2);
//static void excfun(int sigid);
//int setfunc(int dim) ;