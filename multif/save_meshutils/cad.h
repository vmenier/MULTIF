#define NBRCTRMAX 1000 //Max number of control vertices

#define CADBSPLINE 1 
#define CADLINE    2 
#define CADBEZIER  3 
#define CADLOOP  3 

enum CadKwdCod
{
	CadLines,    \
	CadVertices, \
	CadSurfaces, \
	CadLoops,    \
	CadVerMax,   \
	CadSurfMax,  \
	CadLineMax,  \
	CadLoopMax,  \
	CadCtrVer,   \
	CadKwdSize   
};


typedef struct CadLine
{
	int  Typ; // straight line, bspline, bez ?..
  int  Tag;
  int  Ver[2];

	int  HeadCtr; //  index of 1st ctr point in cad->CtrVer
	int  NbrCtr;  //  number of control points

} CadLine, *pCadLine;


typedef struct CadVertex
{
  int       Tag;
	double    Crd[3];
	short int active;
} CadVertex, *pCadVertex;


typedef struct CadLoop
{
	int NbrLin;
	int Lin[10];
} CadLoop, *pCadLoop;

typedef struct CadSurface
{
	int Typ;
	int Loop; // Loop id
} CadSurface, *pCadSurface;


typedef struct S_Cad
{
	
	int NbrLin;
	int LinMax; // Max index for a vertex
	CadLine *Lin;
	int *LidOld;
	int *LidNew;
	
	int  NbrCtrVer;
	int* CtrVer;

	int NbrVer;
	int VerMax;
	CadVertex *Ver;    
	int *VidNew;      // vertex renumbering
	int *VidOld;      // store indices such as put in the .geo file
	
	int NbrLoo;
	int LooMax; 
	CadLoop *Loo;
	
	int NbrSrf;
	int SrfMax;
	CadSurface *Srf; 
	
} Cad;