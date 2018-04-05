
// --- Flags for nozzle CAD patches
#define NOZZLEUP      1 
#define NOZZLEDOWN    2

enum KwdNozzleSize
{
	KwdnCoefs_center,   \
	KwdnKnots_center,   \
	KwdnCoefs_r1,       \
	KwdnKnots_r1,       \
	KwdnCoefs_r2,       \
	KwdnKnots_r2,       \
	MaxKwdNozzleSize
};


enum KwdBaselineParam
{
	BasInletx      , \
	BasInlety      , \
	BasInletz      , \
	BasInletr      , \
	BasOutletx     , \
	BasOutlety     , \
	BasOutletz     , \
	BasOutletr1    , \
	BasOutletr2    , \
	BasOutletzcut  , \
	BasInletzcut   , \
	BasOutletycut  , \
	BasInletycut   , \
	BasInletTheta  , \
	BasOutletTheta , \
	BasMaxKwd
};



typedef struct CadBspline
{
	int NbrCoefs;
	double *Coefs;

	int NbrKnots;
	double *Knots;

} CadBspline, *pCadBspline;


typedef struct CadNozzle
{
	CadBspline *Bsp_center;
	CadBspline *Bsp_r1;
	CadBspline *Bsp_r2;
	
	double ThetaCutIn;
	double ThetaCutOut;

} CadNozzle, *pCadNozzle;



