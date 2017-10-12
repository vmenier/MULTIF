double *allocateDoubleVector(int n);

double dynamicViscosity(double T);

double prair(double T);

double kair(double T);

double cpair(double T);

double areaMachFunc(double M, double g);

double dM2dx(double M2, double g, double A, double dAdx, double D, double Cf, 
    double Ts, double dTsdx);

double evaldM2dx(double xeval, double M2, double* xgeo, double* rgeo, int ngeo, 
        double g, double* x, double* Cf, double* Ts, double* dTs, int nb);

double findApparentThroat(double xstart, double h0, double M2, 
    double* xgeo, double* rgeo, int ngeo, double g, 
    double* x, double* Cf, double* Ts, double* dTs, int nb);

int analyzeNozzle(double* xgeo, double* rgeo, int ngeo, int nbreaks,
    double* xwalltemp, double* walltemp, int nwalltemp, 
    double* xlayer1, double* tlayer1, int nlayer1, double k1,
    double* xlayer2, double* tlayer2, int nlayer2, double k2,
    double* xlayer3, double* tlayer3, int nlayer3, double k3,
    double* xlayer4, double* tlayer4, int nlayer4, double k4,
    double* xlayer5, double* tlayer5, int nlayer5, double k5,
    double tsi, double dtsi, double psi, double cfi,
    double missionmach, double g, double gasconstant,
    double hinf, double tenv, double cenv, double penv,
    double eps1, int maxiter, int maxstep, 
    double eps2, double ns, double himag, double hminmag, double hmaxmag, 
    double singularitydy,
    double* x, double* temp, double* p, double* rho, double* u, double* mach,
    double* tempinside, double* tempoutside, double* netthrust);