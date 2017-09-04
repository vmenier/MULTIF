#include "meshutils.h"

/* given a 1-D array of x-coordinates and y-coordinates of nodes for a 
piecewise linear function, as well as a 1-D array of x-values, calculate
value and gradient at each x-value, nx = # of x-values, nn = # of nodes.
Assume monotonic increasing x array. */
void piecewiseLinear(double *xnodes, double *ynodes, double *x, double *y,
		     double *dydx, int nx, int nn)
{

  int iL = 0; // index for node on left
  int iR = 0; // index for node on right
  int jjStart = 1; // index to start searching (speeds up search),
                   // assumes x is monotonic increasing

  for( int ii = 0; ii < nx; ii++ ) {
  
    // Determine which two values x[ii] is between
    if( x[ii] < xnodes[0] ) {
      iL = -1;
      iR = 0;
    } else {
      iL = nn-1;
      iR = nn;
      for( int jj = jjStart; jj < nn; jj++ ) {
        if( x[ii] < xnodes[jj] ) {
          iL = jj-1;
          iR = jj;
          jjStart = jj-1;
          break;
        }
      }
/*      for( int jj = nn-1; jj > -1; jj-- ) {
        if( x[ii] >= xnodes[jj] ) {
          iL = jj;
          iR = jj+1;
          break;
        }*/
    }
    
    // Calculate value at x[ii] with linear extrapolation 
    if( iL == -1 ) {
      y[ii] = ynodes[0] + (ynodes[1]-ynodes[0])*(x[ii]-xnodes[0])/(xnodes[1]-xnodes[0]);
    } else if( iR == nn ) {
      y[ii] = ynodes[nn-2] + (ynodes[nn-1]-ynodes[nn-2])*(x[ii]-xnodes[nn-2])/(xnodes[nn-1]-xnodes[nn-2]);
    } else {
      y[ii] = ynodes[iL] + (ynodes[iR]-ynodes[iL])*(x[ii]-xnodes[iL])/(xnodes[iR]-xnodes[iL]);
    }

    // Calculate gradient at x[ii] with linear extrapolation
    if( iL == -1 ) {
      dydx[ii] = (ynodes[1]-ynodes[0])/(xnodes[1]-xnodes[0]);
    } else if( iR == nn ) {
      dydx[ii] = (ynodes[nn-1]-ynodes[nn-2])/(xnodes[nn-1]-xnodes[nn-2]);
    } else {
      dydx[ii] = (ynodes[iR]-ynodes[iL])/(xnodes[iR]-xnodes[iL]);
    }    
  
  }

  return;

}
