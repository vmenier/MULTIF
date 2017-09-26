// Find x at which piecewise linear function obtains a minimum. If multiple
// minima are found, the last one is returned.
double findPiecwiseLinearMinimumLocation(double *xgeo, double *rgeo, int ngeo)
{
  double xloc, rmin;
  xloc = xgeo[0];
  rmin = rgeo[0];
  for (int i = 1; i < ngeo; i++)
  {
    if (rgeo[i] <= rmin)
    {
      xloc = xgeo[i];
      rmin = rgeo[i];
    }
  }
  return xloc;
}

/* Linear interpolation with linear extrapolation for value(s) y at x using 
arrays (xn, yn). xn, yn = [0 ... nn-1] and x, y = [0 ... nx-1]. Monotonically
increasing x is assumed. Use js to jumpstart interpolation at index js. */
void interp1(double *xn, double *yn, int nn, double *x, double *y, int nx, int js)
{

  int iL = 0;      // index for node on left
  int iR = 0;      // index for node on right
  int jStart = js; // begin interpolation search at index jjStart

  for (int i = 0; i < nx; i++)
  {

    // Determine which two values x[i] is between
    if (x[i] < xn[0])
    {
      iL = -1;
      iR = 0;
    }
    else
    {
      iL = nn - 1;
      iR = nn;
      for (int j = jStart; j < nn; j++)
      {
        if (x[i] < xn[j])
        {
          iL = j - 1;
          iR = j;
          jStart = j - 1;
          break;
        }
      }
    }

    // Calculate value at x[i] with linear extrapolation
    if (iL == -1)
    {
      y[i] = yn[0] + (yn[1] - yn[0]) * (x[i] - xn[0]) / (xn[1] - xn[0]);
    }
    else if (iR == nn)
    {
      y[i] = yn[nn - 2] + (yn[nn - 1] - yn[nn - 2]) * (x[i] - xn[nn - 2]) / (xn[nn - 1] - xn[nn - 2]);
    }
    else
    {
      y[i] = yn[iL] + (yn[iR] - yn[iL]) * (x[i] - xn[iL]) / (xn[iR] - xn[iL]);
    }
  }

  return;
}


/* Gradient estimation for linearly interpolated function dydx at x using 
arrays (xn, yn). xn, yn = [0 ... nn-1] and x, y = [0 ... nx-1]. Monotonically 
increasing x is assumed. Use js to jumpstart interpolation at index js. */
void interp1grad(double *xn, double *yn, int nn, double *x, double *dydx, int nx, int js)
{

  int iL = 0;      // index for node on left
  int iR = 0;      // index for node on right
  int jStart = js; // begin interpolation search at index jjStart

  for (int i = 0; i < nx; i++)
  {

    // Determine which two values x[i] is between
    if (x[i] < xn[0])
    {
      iL = -1;
      iR = 0;
    }
    else
    {
      iL = nn - 1;
      iR = nn;
      for (int j = jStart; j < nn; j++)
      {
        if (x[i] < xn[j])
        {
          iL = j - 1;
          iR = j;
          jStart = j - 1;
          break;
        }
      }
    }

    // Calculate gradient at x[ii] with linear extrapolation
    if (iL == -1)
    {
      dydx[i] = (yn[1] - yn[0]) / (xn[1] - xn[0]);
    }
    else if (iR == nn)
    {
      dydx[i] = (yn[nn - 1] - yn[nn - 2]) / (xn[nn - 1] - xn[nn - 2]);
    }
    else
    {
      dydx[i] = (yn[iR] - yn[iL]) / (xn[iR] - xn[iL]);
    }
  }

  return;
}

/* Cumulative trapezoidal integration of y over x using n steps. Return result 
in yint. An average slope is used to calculate y values at midpoints between
x values. Linear extrapolation is used at the ends of the interval. */
void cumtrapint(double *x, double *y, double *yint, int n) {

  double dxleft, dxright, dx;
  double dyleft, dyright;
  double x1, x2, y1, y2;
  double mleft, mright, m;
  double ycum = 0.;

  for(int i = 0; i < n; i++) {

    if(i == 0) {
      dxleft = x[1] - x[0];
      dxright = x[1] - x[0];
      dyleft = y[1] - y[0];
      dyright = y[1] - y[0];
    } else if(i == n-1) {
      dxleft = x[n-1] - x[n-2];
      dxright = x[n-1] - x[n-2];
      dyleft = y[n-1] - y[n-2];
      dyright = y[n-1] - y[n-2];
    } else {
      dxleft = x[i] - x[i-1];
      dxright = x[i+1] - x[i];
      dyleft = y[i] - y[i-1];
      dyright = y[i+1] - y[i];
    }

    mleft = dyleft/dxleft;
    mright = dyright/dxright;
    m = (mleft + mright)/2.;

    dx = dxright/2. + dxleft/2.;

    x1 = x[i] - dxleft/2.;
    y1 = y[i] - m*(x[i] - x1);
    x2 = x[i] + dxright/2.;
    y2 = y[i] + m*(x2 - x[i]);

    ycum += y1*y2*dx/2.;
    yint[i] = ycum;

  }

  return;

}