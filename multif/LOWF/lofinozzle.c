#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../meshutils/piecewise.h"
#include "odeint.h"

#define PI 3.14159265358979323846


double *allocateDoubleVector(int n) {
    double *vector = (double*) malloc((n+1)*sizeof(double));
    if( vector == NULL ) {
        perror("malloc failed");
    }
    return vector;
}


// Sutherland's Law of dynamic viscosity for air
double dynamicViscosity(double T) {
    return 1.716e-5*pow((T/273.15),1.5)*(273.15 + 110.4)/(T + 110.4);
}


// Prandtl number for air
double prair(double T) {
    static double prval[31] = {0.744, 0.736, 0.728, 0.72, 0.713, 0.707, 0.701,
        0.697, 0.692, 0.688, 0.684, 0.68, 0.68, 0.68, 0.682, 0.684, 0.687,
        0.69, 0.693, 0.696, 0.699, 0.702, 0.704, 0.707, 0.709, 0.711,
        0.713, 0.715, 0.717, 0.719, 0.722};
    static double temp[31] = {175., 200., 225., 250., 275., 300., 325., 350., 
        375., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 
        900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 
        1400., 1500.};
    double pr;

    interp1(temp, prval, 31, &T, &pr, 1, 0);

    return pr;
    
}


// Thermal conductivity for air
double kair(double T) {
    static double kval[31] = {0.01593, 0.01809, 0.0202, 0.02227, 0.02428, 
        0.02624, 0.02816, 0.03003, 0.03186, 0.03365, 0.0371, 0.04041, 0.04357,
        0.04661, 0.04954, 0.05236, 0.05509, 0.05774, 0.0603, 0.06276, 0.0652, 
        0.06754, 0.06985, 0.07209, 0.07427, 0.0764, 0.07849, 0.08054, 0.08253, 
        0.0845, 0.08831};
    static double temp[31] = {175., 200., 225., 250., 275., 300., 325., 350., 
        375., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 
        900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 
        1400., 1500.};
    double k;

    interp1(temp, kval, 31, &T, &k, 1, 0);

    return k; // W/m-K
    
}


// Specific heat at constant pressure for air
double cpair(double T) {
    static double cpval[31] = {1002.3, 1002.5, 1002.7, 1003.1, 1003.8,
        1004.9, 1006.3, 1008.2, 1010.6, 1013.5, 1020.6, 1029.5, 1039.8,
        1051.1, 1062.9, 1075.0, 1087.0, 1098.7, 1110.1, 1120.9, 1131.3,
        1141.1, 1150.2, 1158.9, 1167.0, 1174.6, 1181.7, 1188.4, 1194.6,
        1200.5, 1211.2};
    static double temp[31] = {175., 200., 225., 250., 275., 300., 325., 350., 
        375., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 
        900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 
        1400., 1500.};
    double cp;

    interp1(temp, cpval, 31, &T, &cp, 1, 0);

    return cp; // J/kg-K
    
}


// Area-Mach function from 1-D mass conservation equations
double areaMachFunc(double M, double g) {
    return pow(((g+1.)/2.),((g+1.)/(2.*(g-1.))))*M/pow((1.+(g-1.)*pow(M,2)/2.),
        ((g+1.)/(2.*(g-1.))));
}


// Governing equation of motion for quasi-1D flow
double dM2dx(double M2, double g, double A, double dAdx, double D, double Cf, 
    double Ts, double dTsdx) {
    double delta = 1e-6; // to keep denominator from becoming zero
    double denom;
    if( fabs(1-M2) <= delta ) { // it is unclear whether this actually helps
        denom = delta;
    } else {
        denom = 1. - M2;
    }
    //return 2.*M2*(1.+(g-1.)*M2/2.)/(1.-M2+delta)*(-dAdx/A + 2.*g*M2*Cf/D + (1.+g*M2)*dTsdx/(2.*Ts));
    return 2.*M2*(1.+(g-1.)*M2/2.)/(denom)*(-dAdx/A + 2.*g*M2*Cf/D + (1.+g*M2)*dTsdx/(2.*Ts));
}


// Evaluate governing equation of motion at xeval & M2 for quasi-1D flow 
// assuming linear interpolation for xgeo, rgeo, x, Cf, Ts, dTs, nb
double evaldM2dx(double xeval, double M2, double* xgeo, double* rgeo, int ngeo, 
    double g, double* x, double* Cf, double* Ts, double* dTs, int nb) {
    
    double r, drdx;
    double a, da, d;
    double c, t, dt;
    
    interp1(xgeo, rgeo, ngeo, &xeval, &r, 1, 1);
    interp1grad(xgeo, rgeo, ngeo, &xeval, &drdx, 1, 1);
    a = PI*pow(r,2);
    da = 2*PI*r*drdx;
    d = 2*r;

    interp1(x, Cf, nb, &xeval, &c, 1, 1);
    interp1(x, Ts, nb, &xeval, &t, 1, 1);
    interp1(x, dTs, nb, &xeval, &dt, 1, 1);
    
    return dM2dx(M2, g, a, da, d, c, t, dt);    
}


// Find apparent throat of nozzle
double findApparentThroat(double xstart, double h0, double M2, 
    double* xgeo, double* rgeo, int ngeo, double g, 
    double* x, double* Cf, double* Ts, double* dTs, int nb) {

    // Initialize variables
    //double frac = 0.2; // fraction of nozzle length to search around min area point
    //double length = xgeo[ngeo-1]-xgeo[0]; // nozzle length
    //double h0 = 1e-3; // initial step size
    //double M2 = pow(1.0001,2);
    double xt; // throat location
    double h;

    // Determine approximate throat location by minimum area point
    double x1, f1, x2, f2, xstop;
    x2 = xstart;
    //x2 = findPiecwiseLinearMinimumLocation(xgeo, rgeo, ngeo);
    f2 = evaldM2dx(x2, M2, xgeo, rgeo, ngeo, g, x, Cf, Ts, dTs, nb);
    xt = x2;

    // Next, search around min area throat location for true apparent throat
    if( f2 < 0 ) {
        h = h0;
        //xstop = x2 + frac*length;
        xstop = xgeo[ngeo-1];
    } else {
        h = -h0;
        //xstop = x2 - frac*length;
        xstop = xgeo[0];
    }
    while(1) {
        f1 = f2;
        x1 = x2;
        x2 = x1 + h;
        f2 = evaldM2dx(x2, M2, xgeo, rgeo, ngeo, g, x, Cf, Ts, dTs, nb);
        //printf("%f, %f\n",x1,x2);
        //printf("%f, %f\n",f1,f2);
        if( f1*f2 <= 0 ) {
            break; // bracketing is successful
        }
        if( h > 0 && x2 > xstop ) {
            printf("Apparent throat set to nozzle exit.\n");
            //return xt;
            return xgeo[ngeo-1];
        } else if( h < 0 && x2 < xstop ) {
            printf("Bracketing for apparent throat failed (leftward search).\n");
            return xt;
        }
    }

    // Perform bisection until location of apparent throat is accurate to eps
    double x3, f3;
    double abserr = 1.;
    double eps = 1e-10;
    while(abserr > eps) {
        x3 = (x1 + x2)/2.;
        f3 = evaldM2dx(x3, M2, xgeo, rgeo, ngeo, g, x, Cf, Ts, dTs, nb);
        if(f1*f3 <= 0) { // throat b/w points 1 and 3
            x2 = x3;
        } else { // throat b/w points 3 and 2
            x1 = x3;
            f1 = f3;
        }
        abserr = fabs(x1-x2);
    }

    xt = x3;

    return xt;

}


/* Integrate M^2 across nozzle.
*/
int integrateM2(double* xgeo, double* rgeo, int ngeo, int nbreaks,
    double g, int maxstep, int maxattempts, int maxattempts2, double dxsave,
    double eps2, double ns, double himag, double hminmag, double hmaxmag, 
    double deltaMPrior, double deltaMPost, double xi, double xt, double xe,
    double* x, double* cf, double* ts, double* dts,
    double* xsave1, double* ysave1, double* xsave2, double* ysave2,
    int* nsave1, int* nsave2, double* xtguess) 

    {

    double yi, xs, xf, hi, hmin, hmax;
    double xterm;

    int converged = 0;

    // Integrate for M^2 backward from throat to inlet
    // Assume subsonic flow in this region for now
    yi = pow(1.+deltaMPrior,2);
    xs = xt;
    xf = xi;
    hi = -himag;
    hmin = -hminmag;
    hmax = -hmaxmag;
    xterm = odeint(xs, xf, yi, maxstep, hi, hmin, hmax, eps2, xsave1, ysave1, 
        dxsave, ns, nsave1, evaldM2dx, xgeo, rgeo, ngeo, g, x, cf, ts, dts, 
        nbreaks);
    if( xterm < xi + 1e-6 ) {
        converged = 1;
    } else {
        *xtguess = xterm; // update location for guess of throat
    }
    printf("completed LHS integration\n");

    // Integrate for M^2 forward from throat to exit
    // Assume supersonic flow for now
    if( xt >= xe ) {
        xsave2[0] = xe;
        ysave2[0] = 1.;
        *nsave2 = 1;
        printf("skipping RHS integration\n");
    } else {
        yi = pow(1.+deltaMPost,2);
        xs = xt;
        xf = xe;
        hi = himag;
        hmin = hminmag;
        hmax = hmaxmag;
        xterm = odeint(xs, xf, yi, maxstep, hi, hmin, hmax, eps2, xsave2, 
            ysave2, dxsave, ns, nsave2, evaldM2dx, xgeo, rgeo, ngeo, g, x, 
            cf, ts, dts, nbreaks);
        if( xterm > xe - 1e-6 ) {
            converged = converged == 0 ? 0 : 1;
        } else {
            converged = 0;
            *xtguess = xterm; // update location for guess of throat
        }
        printf("completed RHS integration\n");
    }

    return converged; // 1 = converged, 0 = not converged

    }


/* Analyze an axisymmetric nozzle with 5 wall layers defined by piecewise linear
functions using the quasi-1D Navier Stokes equation. Data is returned in
vectors x, temp, p, rho, u, mach, tempinside, and tempoutside, all of which are
length nbreaks. The net thrust is returned in netthrust. Nozzle geometry is
defined by xgeo and rgeo, of length ngeo. Wall thicknesses from the inside out 
are defined by xlayer1 through xlayer5 (tlayer1 to tlayer5, and nlayer1 to
nlayer5). k1 through k5 give the thermal conductivity in the radial direction
for each layer. Additional parameters include the inlet stagnation temperature
tsi, the (initially) constant gradient of tsi dtsi, inlet stag pressure psi, and
(initially) constant friction coef cfi. The flight mach number missionmach is 
used in the netthrust calculation. g is the ratio of specific heats and 
gasconstant is the ideal gas constant for air. Environmental temperature, 
speed of sound and pressure is given by tenv, cenv, and penv. The heat xfer 
coef to the environment from the nozzle exterior surface is given by hinf. The 
Gauss-Seidel conjugate heat xfer calculations terminate when either the relative
error in temperature at the outlet is less than eps1 or maxiter iterations is 
reached. The ODE integration terminates when maxstep steps is reached or the 
absolute error of integration is less than eps2. The initial step size magnitude
is himag, minimum allowable size is hminmag, and max allowable size is hmaxmag.
Integer ns gives the maximum length of the vectors where ODE integration data 
is saved. singularitydy is the value of Mach at which to start integrating
around the singularity M=1. */
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
    double* x, double* temp, double* p, double *rho, double* u, double* mach,
    double* tempinside, double* tempoutside, double* netthrust) 

    {

    // THINGS TO ADD:
    // 1) Consider implementing a shock case...
    
    printf("Beginning low-fi nozzle analysis.\n");

    // Initial parameters
    double xi = xgeo[0]; // inlet location
    double xe = xgeo[ngeo-1]; // outlet location
    double ri = rgeo[0]; // inlet radius
    double xt; // apparent throat location

    // Declare Gauss-Seidel fluid-thermal iteration properties
    double err;

    // Define ODE integration properties
    double dxsave = (xgeo[ngeo-1]-xgeo[0])/((double)ns);

    // Initialize loop variables
    double *r, *ts, *dts, *cf;
    double *ps, *re, *hf;
    double *xmach, *machtmp;
    r = allocateDoubleVector(nbreaks);
    ts = allocateDoubleVector(nbreaks);
    dts = allocateDoubleVector(nbreaks);
    cf = allocateDoubleVector(nbreaks);
    ps = allocateDoubleVector(nbreaks);
    re = allocateDoubleVector(nbreaks);
    hf = allocateDoubleVector(nbreaks);
    double *rwallprime = allocateDoubleVector(nbreaks);
    double *ro = allocateDoubleVector(nbreaks);
    double *rtotalprime = allocateDoubleVector(nbreaks);
    double *tsintegrand = allocateDoubleVector(nbreaks);
    double *tsintegral = allocateDoubleVector(nbreaks);
    double qwflux;
    double tprimeratio, reprimeratio, cfincomp;

    for (int i = 0; i < nbreaks; i++) {
        x[i] = xi + (xe-xi)*((double)i)/((double)(nbreaks-1));
        ts[i] = tsi + dtsi*(x[i]-xi);
        dts[i] = dtsi;
        cf[i] = cfi;
    }
    interp1(xgeo, rgeo, ngeo, x, r, nbreaks, 0);
    double tempold = ts[nbreaks-1];

    // Calculate wall thermal resistance & estimate outer wall radius
    double rotemp, ritemp, ttemp;
    for(int i = 0; i < nbreaks; i++) { // at each x station
        rwallprime[i] = 0.;
        // thermal layer
        interp1(xlayer1, tlayer1, nlayer1, &x[i], &ttemp, 1, 0);
        ritemp = r[i];
        rotemp = ritemp + ttemp;
        rwallprime[i] += log(rotemp/ritemp)/(2.*PI*k1);
        // air gap
        interp1(xlayer2, tlayer2, nlayer2, &x[i], &ttemp, 1, 0);
        ritemp = rotemp;
        rotemp = ritemp + ttemp;
        rwallprime[i] += log(rotemp/ritemp)/(2.*PI*k2);
        // inner load layer
        interp1(xlayer3, tlayer3, nlayer3, &x[i], &ttemp, 1, 0);
        ritemp = rotemp;
        rotemp = ritemp + ttemp;
        rwallprime[i] += log(rotemp/ritemp)/(2.*PI*k3); 
        // middle load layer
        interp1(xlayer4, tlayer4, nlayer4, &x[i], &ttemp, 1, 0);
        ritemp = rotemp;
        rotemp = ritemp + ttemp;
        rwallprime[i] += log(rotemp/ritemp)/(2.*PI*k4); 
        // outer load layer
        interp1(xlayer5, tlayer5, nlayer5, &x[i], &ttemp, 1, 0);
        ritemp = rotemp;
        rotemp = ritemp + ttemp;
        rwallprime[i] += log(rotemp/ritemp)/(2.*PI*k5);
        // assign outer wall radius
        ro[i] = rotemp;
    }

    // If wall temp is specified, perform these calculations outside the loop
    if( nwalltemp > 0 ) {

        // Inside wall temperature is fixed
        interp1(xwalltemp, walltemp, nwalltemp, x, tempinside, nbreaks, 1);
        
        // Estimate for stagnation temperature gradient (from old code)
        // Approx: Gradient of static wall temp = gradient of flow stagnation temp
        interp1grad(xwalltemp, walltemp, nwalltemp, x, dts, nbreaks, 1);
        
    }

    // Begin Gauss-Seidel iterations for aero-thermal analysis
    double xtguess;
    int nsave1 = 0;
    int nsave2 = 0; 
    int nm;
    int counter;
    int counter2;
    int maxattempts = 3;
    int maxattempts2 = 100;
    double deltaM, deltaMold;
    int converged = 0;

    // Allocate data arrays big enough to keep all integration results
    double *xsave1 = NULL;
    double *xsave2 = NULL;
    double *ysave1 = NULL;
    double *ysave2 = NULL;
    xsave1 = allocateDoubleVector(ns+2);
    ysave1 = allocateDoubleVector(ns+2);
    xsave2 = allocateDoubleVector(ns+2);
    ysave2 = allocateDoubleVector(ns+2);

    for(int i=0; i < maxiter; i++) {

        counter = 0;

        // Initial guess for throat location
        xtguess = findPiecwiseLinearMinimumLocation(xgeo, rgeo, ngeo);

        // Run integration until correct apparent throat is found and 
        // integration succeeds
        while( converged != 1 ) {

            xt = findApparentThroat(xtguess, fabs(himag), pow(1+singularitydy,2), 
            xgeo, rgeo, ngeo, g, x, cf, ts, dts, nbreaks);
            printf("\nLocation of minimum is: %f\n", xt);

            if( counter > 0 ) {
                printf("Attempting integration again from new apparent throat x = %0.6f\n",xt);            
            }

            converged = integrateM2(xgeo, rgeo, ngeo, nbreaks,
                g, maxstep, maxattempts, maxattempts2, dxsave,
                eps2, ns, himag, hminmag, hmaxmag, 
                -singularitydy, singularitydy, xi, xt, xe,
                x, cf, ts, dts,
                xsave1, ysave1, xsave2, ysave2,
                &nsave1, &nsave2, &xtguess);

            counter++;

            // If integration fails too many times (likely due to multiple 
            // throats), enter a failsafe mode where only subsonic flow is 
            // present in the nozzle with relaxed choking at selected throat.
            // The Mach number at the throat is determined through bisection.
            if( counter > maxattempts ) {

                printf("\n*****ODE integration failed after %d attempts.\n",maxattempts);
                printf("*****Entering fail safe mode with subsonic flow.\n");

                counter2 = 0;

                printf("\nRetaining same throat location at: %f\n", xt);

                deltaM = singularitydy;

                // First bracket the Mach number
                while( converged != 1 ) {

                    if( counter2 > 0 ) {
                        // Increase deltaM to relax the integration more
                        deltaM *= 2.;
                        printf("Attempting integration with new DeltaM = %0.6f\n",deltaM);          
                    }

                    converged = integrateM2(xgeo, rgeo, ngeo, nbreaks,
                        g, maxstep, maxattempts, maxattempts2, dxsave,
                        eps2, ns, himag, hminmag, hmaxmag, 
                        -deltaM, -deltaM, xi, xt, xe,
                        x, cf, ts, dts,
                        xsave1, ysave1, xsave2, ysave2,
                        &nsave1, &nsave2, &xtguess); 

                    if( counter2 > maxattempts2 ) {
                        printf("*****ODE integration failed after %d attempts.\n",maxattempts2);
                        break;
                    }

                    deltaMold = deltaM;
                    counter2++;

                }

                printf("\n*****Bracketing of Mach number succesful.\n");

                // Next bisect to obtain correct Mach number.
                double dMhi = deltaM;
                double dMlo = deltaM/2.;
                double dMmid;
                double relerrM = (dMhi-dMlo)/dMhi;
                printf("\n*****Beginning bisection.\n");
                printf("dMhi: %0.6f\n",dMhi);
                printf("dMlo: %0.6f\n",dMlo);
                while( relerrM > 1e-3 ) {

                    dMmid = (dMhi + dMlo)/2.;
                    printf("Bisecting Mach number. Using new DeltaM = %0.6f\n",dMmid);

                    converged = integrateM2(xgeo, rgeo, ngeo, nbreaks,
                        g, maxstep, maxattempts, maxattempts2, dxsave,
                        eps2, ns, himag, hminmag, hmaxmag, 
                        -dMmid, -dMmid, xi, xt, xe,
                        x, cf, ts, dts,
                        xsave1, ysave1, xsave2, ysave2,
                        &nsave1, &nsave2, &xtguess); 

                    if( converged == 0 ) {
                        dMlo = dMmid;
                    } else {
                        dMhi = dMmid;
                    }

                    relerrM = fabsf((dMhi-dMlo)/dMhi);

                }
                
                break;
            }

        }

        // Concatenate data for M^2 and combine into array for mach
        nm = nsave1+nsave2-1;
        xmach = allocateDoubleVector(nm);
        machtmp = allocateDoubleVector(nm);
        for(int j=0; j<nsave1-1; j++) {
            machtmp[j] = sqrt(ysave1[nsave1-j-1]);
            xmach[j] = xsave1[nsave1-j-1];
        }
        machtmp[nsave1-1] = sqrt( (ysave1[0] + ysave2[0])/2. );
        xmach[nsave1-1] = (xsave1[0] + xsave2[0])/2.;
        for(int j=1; j<nsave2; j++) {
            machtmp[nsave1-1+j] = sqrt(ysave2[j]);
            xmach[nsave1-1+j] = xsave2[j];
        }

        interp1(xmach, machtmp, nm, x, mach, nbreaks, 0);
        
        // Free memory for temporary M^2 storage
        free(xsave1);
        free(ysave1);
        free(xsave2);
        free(ysave2);
        free(xmach);
        free(machtmp);

/*        for(int j=0; j<nbreaks; j++) {
            printf("%0.8f %0.8f\n", x[j], mach[j]);
        } */

        // Update nozzle fluid variables
        for(int j=0; j<nbreaks; j++) {

            temp[j] = ts[j]/(1. + (g-1.)*pow(mach[j],2)/2.);
            ps[j] = psi*(pow(ri,2)/pow(r[j],2))*(areaMachFunc(mach[0],g)/
                areaMachFunc(mach[j],g))*sqrt(ts[j]/tsi);
            p[j] = ps[j]/pow( 1. + (g-1.)*pow(mach[j],2)/2., (g/(g-1.)) );
            rho[j] = p[j]/(gasconstant*temp[j]);
            u[j] = mach[j]*sqrt(g*gasconstant*temp[j]);
            re[j] = rho[j]*u[j]*2.*r[j]/dynamicViscosity(temp[j]);

            // Heat transfer coefficient from fluid to interior nozzle wall,
            // estimated using Chilton-Colburn analogy
            hf[j] = pow(prair(temp[j]),-2./3.)*rho[j]*cpair(temp[j])*
                u[j]*cf[j]/2.;

            // Total thermal resistance
            rtotalprime[j] = 1./(hf[j]*PI*2.*r[j]) + rwallprime[j] + 
                1./(hinf*PI*2.*ro[j]);

            // Integrand for integral used in calculation of stagnation temp
            tsintegrand[j] = 1./(rtotalprime[j]*rho[j]*u[j]*PI*pow(r[j],2)*
                cpair(temp[j]));
            
        }

        // Calculate other heat transfer quantities of interest
        if( nwalltemp > 0 ) { // if wall temperature is specified

            // Inside wall temp is already specified and set above.
            // Approximation for stagnation temperature gradient (see above)
            // makes this fixed as well. It is precalculated above.

            // Calculate heat flux, stagnation temperature, and external wall temp
            for(int j = 0; j<nbreaks; j++) {
                qwflux = (ts[j] - tenv)/rtotalprime[j]/(PI*2.*r[j]);
                ts[j] = tempinside[j] + qwflux/hf[j];
                tempoutside[j] = tenv + qwflux/hinf;
            }

        } else { // if wall temperature is calculated

            cumtrapint(x, tsintegrand, tsintegral, nbreaks);
            for(int j = 0; j<nbreaks; j++) {
    
                // Estimate stagnation temperature and its derivative
                ts[j] = tenv*(1. - exp(-tsintegral[j])) + tsi*exp(-tsintegral[j]);
                dts[j] = (tenv - ts[j])/(rtotalprime[j]*rho[j]*u[j]*PI*
                    pow(r[j],2)*cpair(temp[j]));
    
                // Calculate heat flux 
                qwflux = (ts[j] - tenv)/rtotalprime[j]/(PI*2.*r[j]);
                tempinside[j] = ts[j] - qwflux/hf[j];
                tempoutside[j] = tenv + qwflux/hinf;
            
            }

        }

        // Update friction coefficient
        for(int j = 0; j<nbreaks; j++) {
            tprimeratio = 1. + 0.035*pow(mach[j],2) + 0.45*(tempinside[j]/temp[j] - 1.);
            reprimeratio = 1./(tprimeratio*pow(tprimeratio,1.5)*
                (1. + 110.4/temp[j])/(tprimeratio + 110.4/temp[j]));
            cfincomp = 0.074/pow(re[j],0.2);
            cf[j] = cfincomp/tprimeratio/pow(reprimeratio,0.2);           
        }

        err = fabs(temp[nbreaks-1] - tempold)/temp[nbreaks-1];
        printf("Error in exit temperature: %0.16f\n", err);
        if( err < eps1 ) {
            printf("Converged.\n");
            break;
        }
        tempold = temp[nbreaks-1];

    }

    // Estimate nozzle thrust
    double drdxexit;
    interp1grad(xgeo, rgeo, ngeo, &x[nbreaks-1], &drdxexit, 1, 1);
    double mdot = rho[0]*u[0]*PI*pow(r[0],2);
    double exitangle = atan(drdxexit);
    double divfactor = (1. + cos(exitangle))/2.;
    *netthrust = divfactor*mdot*(u[nbreaks-1] - missionmach*
        cenv) + (p[nbreaks-1] - penv)*PI*pow(r[nbreaks-1],2);

    // Free memory
    free(r);
    free(ts);
    free(dts);
    free(cf);
    free(ps);
    free(re);
    free(hf);
    free(rwallprime);
    free(ro);
    free(rtotalprime);
    free(tsintegrand);
    free(tsintegral);

    return 0;

}

