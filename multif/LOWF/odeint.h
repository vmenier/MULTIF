#include <stdio.h>

void rkckstep(double* x, double* y, double dydx, double h, double *yout, 
    double *yerr, 
    double(*fcn)(double, double, double*, double*, int, double,
    double*, double*, double*, double*, int),
    double* xgeo, double* rgeo, int ngeo, double g, double* xn, double*cf,
    double* ts, double* dts, int nb);

void rkstepper(double* x, double* y, double dydx, 
    double htry, double hmin, double hmax, double* hnext,
    double eps, 
    double(*fcn)(double, double, double*, double*, int, double,
    double*, double*, double*, double*, int), 
    double* xgeo, double* rgeo, int ngeo, double g, double* xn, double* cf,
    double* ts, double* dts, int nb );

double odeint(double xi, double xe, double yi, 
    int maxstep, double hi, double hmin, double hmax, double eps,
    double* xsave, double* ysave, double dxsave, int ns, int* nsave,
    double(*fcn)(double, double, double*, double*, int, double,
    double*, double*, double*, double*, int), 
    double* xgeo, double* rgeo, int ngeo, double g, double* xn, double* cf,
    double* ts, double* dts, int nb );

