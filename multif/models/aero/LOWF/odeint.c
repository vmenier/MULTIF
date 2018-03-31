#include <stdio.h>
#include <math.h>

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=a,dmaxarg2=b, dmaxarg1 > dmaxarg2 ? dmaxarg1 : dmaxarg2)
static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=a,dminarg2=b, dminarg1 < dminarg2 ? dminarg1 : dminarg2)

/* Takes a Cash-Karp Runge-Kutta step starting at x and given function value y,
derivative dydx, and stepsize h. The function fcn and associated arguments are 
used to provide values of dydx at the various steps required by the Cash-Karp
scheme. Finally, the function value y at x+h is provided in yout, and an
estimate of the error in y at x+h is provided in yerr. */
void rkckstep(double* x, double* y, double dydx, double h, double *yout, 
    double *yerr, 
    double(*fcn)(double, double, double*, double*, int, double,
    double*, double*, double*, double*, int),
    double* xgeo, double* rgeo, int ngeo, double g, double* xn, double*cf,
    double* ts, double* dts, int nb) {

    static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1., a6 = 0.875, 
    b21 = 0.2, b31 = 3./40., b41 = 0.3, b51 = -11./54., b61 = 1631./55296.,
    b32 = 9./40., b42 = -0.9, b52 = 2.5, b62 = 175./512.,
    b43 = 1.2, b53 = -70./27., b63 = 575./13824.,
    b54 = 35./27., b64 = 44275./110592., b65 = 253./4096.,
    c1 = 37./378., c3 = 250./621., c4 = 125./594., c6 = 512./1771.,
    dc5 = -277./14336.;
    double dc1 = c1 - 2825./27648., dc3 = c3 - 18575./48384.,
    dc4 = c4 - 13525./55296., dc6 = c6 - 0.25;
    double k1, k2, k3, k4, k5, k6, ytemp;
    
    k1 = h*dydx;

    ytemp = *y + b21*k1;
    k2 = h*(*fcn)(*x+a2*h, ytemp, xgeo, rgeo, ngeo, g, xn, cf, ts, dts, nb); 

    ytemp = *y + b31*k1 + b32*k2;
    k3 = h*(*fcn)(*x+a3*h, ytemp, xgeo, rgeo, ngeo, g, xn, cf, ts, dts, nb); 

    ytemp = *y + b41*k1 + b42*k2 + b43*k3;
    k4 = h*(*fcn)(*x+a4*h, ytemp, xgeo, rgeo, ngeo, g, xn, cf, ts, dts, nb); 

    ytemp = *y + b51*k1 + b52*k2 + b53*k3 + b54*k4;
    k5 = h*(*fcn)(*x+a5*h, ytemp, xgeo, rgeo, ngeo, g, xn, cf, ts, dts, nb);  

    ytemp = *y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5;
    k6 = h*(*fcn)(*x+a6*h, ytemp, xgeo, rgeo, ngeo, g, xn, cf, ts, dts, nb); 

    *yout = *y + c1*k1 + c3*k3 + c4*k4 + c6*k6;
    *yerr = dc1*k1 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6;

}


/* Stepper function for 5th order Runge-Kutta with error monitoring. Provided
a location x, function value y and function derivative dydx, attempts 
Cash-Karp Runge-Kutta step of size htry. If error is larger than tolerance eps,
the stepsize is decreased, possibly until step size becomes smaller than hmin.
Provides an estimate of next stepsize in hnext. The function fcn and associated
arguments are used to provide values of dydx at various values of x. x is 
updated with the new location post-step, and y is updated with the new function
value post-step. */
void rkstepper(double* x, double* y, double dydx, 
    double htry, double hmin, double hmax, double* hnext,
    double eps, 
    double(*fcn)(double, double, double*, double*, int, double,
    double*, double*, double*, double*, int), 
    double* xgeo, double* rgeo, int ngeo, double g, double* xn, double* cf,
    double* ts, double* dts, int nb ) {

    // Prepare stepping
    double yout, yerr;
    double h = htry;
    double yerratio, h2, hnext2;
    while(1) {

        // Take a step
        rkckstep(x, y, dydx, h, &yout, &yerr, fcn, xgeo, rgeo, ngeo, g, xn, cf, 
            ts, dts, nb);

        yerratio = fabs(yerr)/eps;
        //printf("yerratio: %f\n",yerratio);

        if( yerratio < 1.0 ) { break; } // tolerance met

        h2 = 0.9*h*pow(yerratio,-0.25); // try this smaller step
        h  = ( h >= 0. ? DMAX(h2,0.1*h) : DMIN(h2,0.1*h) ); // factor of 10 only

        if( fabs(h) < fabs(hmin) ) { 
            printf("Minimum step size reached.\n"); 
            break;
        }

        if( *x + h == *x ) {
            printf("Underflow error due to step size.\n");
            break;
        }
        
    }

    // Update step size for next iteration
    if( yerratio > 1.89e-4 ) {
        hnext2 = 0.9*h*pow(yerratio,-0.2);
        *hnext = ( fabs(hnext2) > fabs(hmax) ? hmax : hnext2 );
    } else {
        *hnext = ( fabs(5.*h) > fabs(hmax) ? hmax : 5.*h ); // maximum factor of 5 increase or hmax
    }

    // Update x location
    *x += h;

    // Update function value 
    *y = yout;

}


/* odeint integrates 1 ode starting at xi and ending at xe using a
4th order Runge-Kutta method with adaptive stepsize and initial condition yi.
Initial step size is hi, and min/max stepsize is hmin/hmax. eps specifies absolute
tolerance on error. Function fcn and associated arguments provides derivative
dydx. xsave is a 1-D array large enough (of size ns) to hold values of x 
approximately spaced by dxsave from xi to xe. ysave is a 1-D array of the same
size which saves values of y corresponding to the values in xsave. */
double odeint(double xi, double xe, double yi, 
    int maxstep, double hi, double hmin, double hmax, double eps,
    double* xsave, double* ysave, double dxsave, int ns, int* nsave,
    double(*fcn)(double, double, double*, double*, int, double,
    double*, double*, double*, double*, int), 
    double* xgeo, double* rgeo, int ngeo, double g, double* xn, double* cf,
    double* ts, double* dts, int nb ) {

    //printf("Beginning RK4 ODE solver.\n"); 
    double x, y, dydx, h, hnext;
    int count = -1;
    y = yi;
    x = xi;
    h = hi;
    double xprev = x + dxsave*2.;
    for(int i = 0; i < maxstep; i++) {

        // Derivative can be obtained here if scaling is desired
        dydx = (*fcn)(x, y, xgeo, rgeo, ngeo, g, xn, cf, ts, dts, nb); 
        //printf("\n%i %0.6f %0.6f %f\n",i,x,y,dydx);

        // Store intermediate results
        if( fabs(x - xprev) > fabs(dxsave) && count < ns-1 ) {
            xsave[++count] = x;
            ysave[count] = y;
            *nsave = count+1;
            xprev = x;
        }

        // Ensure step size does not overshoot
        if( (x+h-xe)*(x+h-xi) > 0.0 ) {
            h = xe - x;
        }

        // Call stepper routine
        rkstepper(&x, &y, dydx, h, hmin, hmax, &hnext, eps, fcn, xgeo, rgeo, 
            ngeo, g, xn, cf, ts, dts, nb);

        // Terminate if stepsize is too small
        if( fabs(hnext) < fabs(hmin) ) {
            printf("Step size too small. Terminating odeint at x=%0.8f.\n",x);
            return x;
        }

        // Check that xe has been reached
        if( (x-xe)*(xe-xi) >= 0. ) {
            // Save data at last point
            xsave[++count] = x;
            ysave[count] = y;
            *nsave = count+1;
            //printf("Successful completion of odeint.\n");
            return x;
        }
        
        // Update step for next time
        h = hnext;

    }

    printf("Maximum number of steps reached. Terminating odeint.\n");
    return x;

}


