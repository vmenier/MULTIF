#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "lofinozzle.h"
#include "Python.h"

double *allocateVectorFromPyList(PyObject *pylist, int *n) {

    double *vector = NULL;

    if( PyList_Check(pylist) ) {
        *n = PyList_Size(pylist);
        if( *n < 1 ) { // i.e. if Python list is empty
            return vector;
        } else {
            vector = malloc(*n*sizeof(double));
            if( vector == NULL ) {
                perror("malloc failed");
            }
            for(int i = 0; i<*n; i++) {
                PyObject *oo = PyList_GetItem(pylist,i);
                if( PyFloat_Check(oo) ) {
                    vector[i] = (double) PyFloat_AS_DOUBLE(oo);
                }
            }
        }
    }

    return vector;
}


int analyze(PyObject *pyxgeo, PyObject *pyrgeo, int nbreaks, 
    PyObject *pyxwalltemp, PyObject *pywalltemp,
    PyObject *pyxlayer1, PyObject *pytlayer1, double k1, 
    PyObject *pyxlayer2, PyObject *pytlayer2, double k2, 
    PyObject *pyxlayer3, PyObject *pytlayer3, double k3, 
    PyObject *pyxlayer4, PyObject *pytlayer4, double k4, 
    PyObject *pyxlayer5, PyObject *pytlayer5, double k5, 
    double tsi, double dtsi, double psi, double cfi,
    double missionmach, double g, double gasconstant, double hinf,
    double tenv, double cenv, double penv, double eps1, int maxiter, int maxstep,
    double eps2, int ns, double himag, double hminmag, double hmaxmag, 
    double singularitydy,
    PyObject *pyx, PyObject *pytemp, PyObject *pyp, PyObject *pyrho, PyObject *pyu, 
    PyObject *pymach, PyObject *pytempinside, PyObject *pytempoutside,
    PyObject *pynetthrust) {

    /*
    printf("%i\n",nbreaks);
    printf("%0.16f\n",k1);
    printf("%0.16f\n",k2);
    printf("%0.16f\n",k3);
    printf("%0.16f\n",k4);
    printf("%0.16f\n",k5);
    printf("%0.16f\n",tsi);
    printf("%0.16f\n",dtsi);
    printf("%0.16f\n",psi);
    printf("%0.16f\n",cfi);
    printf("%0.16f\n",missionmach);
    printf("%0.16f\n",g);
    printf("%0.16f\n",gasconstant);
    printf("%0.16f\n",hinf);
    printf("%0.16f\n",tenv);
    printf("%0.16f\n",cenv);
    printf("%0.16f\n",penv);
    printf("%0.16f\n",eps1);
    printf("%i\n",maxiter);
    printf("%0.16f\n",eps2);
    printf("%i\n",ns);
    printf("%0.16f\n",himag);
    printf("%0.16f\n",hminmag);
    printf("%0.16f\n",singularitydy);
    */

    double *xgeo = NULL;
    double *rgeo = NULL;
    double *xwalltemp = NULL;
    double *walltemp = NULL;
    double *xlayer1 = NULL;
    double *tlayer1 = NULL;
    double *xlayer2 = NULL;
    double *tlayer2 = NULL;
    double *xlayer3 = NULL;
    double *tlayer3 = NULL;
    double *xlayer4 = NULL;
    double *tlayer4 = NULL;
    double *xlayer5 = NULL;
    double *tlayer5 = NULL;

    int ngeo = 0;
    int nwalltemp = 0;
    int nlayer1 = 0;
    int nlayer2 = 0;
    int nlayer3 = 0;
    int nlayer4 = 0;
    int nlayer5 = 0;
    int nresult = 0;
    
    // Inner wall geometry
    xgeo = allocateVectorFromPyList(pyxgeo, &ngeo);
    rgeo = allocateVectorFromPyList(pyrgeo, &ngeo);

    // Wall temperature (if specified)
    xwalltemp = allocateVectorFromPyList(pyxwalltemp, &nwalltemp);
    walltemp = allocateVectorFromPyList(pywalltemp, &nwalltemp);

    // Thermal layer
    xlayer1 = allocateVectorFromPyList(pyxlayer1, &nlayer1);
    tlayer1 = allocateVectorFromPyList(pytlayer1, &nlayer1);

    // Air gap
    xlayer2 = allocateVectorFromPyList(pyxlayer2, &nlayer2);
    tlayer2 = allocateVectorFromPyList(pytlayer2, &nlayer2);

    // Inner load layer
    xlayer3 = allocateVectorFromPyList(pyxlayer3, &nlayer3);
    tlayer3 = allocateVectorFromPyList(pytlayer3, &nlayer3);

    // Middle load layer
    xlayer4 = allocateVectorFromPyList(pyxlayer4, &nlayer4);
    tlayer4 = allocateVectorFromPyList(pytlayer4, &nlayer4);

    // Outer load layer
    xlayer5 = allocateVectorFromPyList(pyxlayer5, &nlayer5);
    tlayer5 = allocateVectorFromPyList(pytlayer5, &nlayer5);

    // Allocate outputs
    double *x = allocateDoubleVector(nbreaks);
    double *temp = allocateDoubleVector(nbreaks);
    double *p = allocateDoubleVector(nbreaks);
    double *rho = allocateDoubleVector(nbreaks);
    double *u = allocateDoubleVector(nbreaks);
    double *mach = allocateDoubleVector(nbreaks);
    double *tempinside = allocateDoubleVector(nbreaks);
    double *tempoutside = allocateDoubleVector(nbreaks);
    double *netthrust = allocateDoubleVector(1);

    // Analyze nozzle
    analyzeNozzle(xgeo, rgeo, ngeo, nbreaks,
        xwalltemp, walltemp, nwalltemp,
        xlayer1, tlayer1, nlayer1, k1,
        xlayer2, tlayer2, nlayer2, k2,
        xlayer3, tlayer3, nlayer3, k3,
        xlayer4, tlayer4, nlayer4, k4,
        xlayer5, tlayer5, nlayer5, k5,
        tsi, dtsi, psi, cfi,
        missionmach, g, gasconstant,
        hinf, tenv, cenv, penv,
        eps1, maxiter, maxstep, 
        eps2, ns, himag, hminmag, hmaxmag, singularitydy,
        x, temp, p, rho, u, mach, tempinside, tempoutside, netthrust);

    // Return data to Python
    for(int i = 0; i < nbreaks; i++) {
        PyList_Append(pyx, PyFloat_FromDouble(x[i]));
        PyList_Append(pytemp, PyFloat_FromDouble(temp[i]));
        PyList_Append(pyp, PyFloat_FromDouble(p[i]));
        PyList_Append(pyrho, PyFloat_FromDouble(rho[i]));
        PyList_Append(pyu, PyFloat_FromDouble(u[i]));
        PyList_Append(pymach, PyFloat_FromDouble(mach[i]));
        PyList_Append(pytempinside, PyFloat_FromDouble(tempinside[i]));
        PyList_Append(pytempoutside, PyFloat_FromDouble(tempoutside[i]));
    }
    PyList_Append(pynetthrust,PyFloat_FromDouble(netthrust[0]));
    
    if(xgeo)
        free(xgeo);
    if(rgeo)
        free(rgeo);
    if(xlayer1)
        free(xlayer1);
    if(tlayer1)
        free(tlayer1);
    if(xlayer2)
        free(xlayer2);
    if(tlayer2)
        free(tlayer2);
    if(xlayer3)
        free(xlayer3);
    if(tlayer3)
        free(tlayer3);
    if(xlayer4)
        free(xlayer4);
    if(tlayer4)
        free(tlayer4);
    if(xlayer5)
        free(xlayer5);
    if(tlayer5)
        free(tlayer5);  

    free(x);
    free(temp);
    free(p);
    free(rho);
    free(u);
    free(mach);
    free(tempinside);
    free(tempoutside);
    free(netthrust);

    return 0;
}
