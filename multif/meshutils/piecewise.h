double findPiecwiseLinearMinimumLocation(double* xgeo, double* rgeo, int ngeo);

void interp1(double *xn, double *yn, int nn, double *x, double *y, int nx, int js);

void interp1grad(double *xn, double *yn, int nn, double *x, double *dydx, int nx, int js);

void cumtrapint(double *x, double *y, double *yint, int n);