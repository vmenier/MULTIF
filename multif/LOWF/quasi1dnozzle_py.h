double *allocateVectorFromPyList(PyObject *pylist, int *n);

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
    PyObject *pynetthrust);