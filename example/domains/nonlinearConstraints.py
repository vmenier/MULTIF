"""
Functions used in the construction of nonlinear constraints for the 3-D nozzle
parameterization.

Rick Fenrich 6/19/18
"""

import numpy as np

def value(x):
    """
    Return values of nonlinear constraints for design variables x.
    
    Arguments:
    x: 1-D Numpy array, the first 18 values should correspond to the wall design
       variables defined in general-3d.cfg.
    
    Returns:
    f: 1-D Numpy array with values of nonlinear constraints
    """
    
    f = np.zeros((3,))
    # Quadratic constraint ensuring increase in slopes between 2nd and 3rd segments of 
    # major axis
    f[0] = x[10]*x[7] - 0.439461*x[7] + 0.439461*x[6] - x[6]*x[11] + \
           0.3*x[11] - 0.3*x[10]
    # Quadratic constraint ensuring decrease in slopes between 4th and 5th segments of
    # major axis
    f[1] = x[13]*x[8] - x[7]*x[13] + x[7]*x[12] - x[9]*x[12] + \
           x[9]*x[11] - x[8]*x[11]
    # Quadratic constraint ensuring decrease in slopes between 5th and 6th segments of
    # major axis
    f[2] = 0.92*x[9] - 0.92*x[8] + x[13]*x[8] - 2.33702*x[13] + \
           2.33702*x[12] - x[9]*x[12]
    
    return f
    

def gradient(x):
    """
    Return Jacobian of nonlinear constraints for design variables x.
    
    Arguments:
    x: 1-D Numpy array, the first 18 values should correspond to the wall design
       variables defined in general-3d.cfg.
       
    Returns:
    jac: 2-D Numpy array with Jacobian of nonlinear constraints
    """
    
    jac = np.zeros((3,len(x)))
    # First constraint
    jac[0,6] = 0.439461 - x[11]
    jac[0,7] = x[10] - 0.439461
    jac[0,10] = x[7] - 0.3
    jac[0,11] = -x[6] + 0.3
    # Second constraint
    jac[1,7] = -x[13] + x[12]
    jac[1,8] = x[13] - x[11]
    jac[1,9] = -x[12] + x[11]
    jac[1,11] = x[9] - x[8]
    jac[1,12] = x[7] - x[9]
    jac[1,13] = x[8] - x[7]
    # Third constraint
    jac[2,8] = -0.92 + x[13]
    jac[2,9] = 0.92 - x[12]
    jac[2,12] = 2.33702 - x[9]
    jac[2,13] = x[8] - 2.33702
    
    return jac
    
    
if __name__ == '__main__':
    """
    Test nonlinear constraint value and gradient functions.
    """
    
    x = np.zeros((25,))
    x[6] = 0.7
    x[7] = 1.1
    x[8] = 1.5
    x[9] = 1.8
    x[10] = 0.439460 #0.439462
    x[11] = 0.439457 #2.33696
    x[12] = 0.5 #2.33699
    x[13] = 0.6 #0.2.33701
    
    print("Constraint values:")
    fval = value(x)
    print(fval)
    
    jac = gradient(x)
    jac_fd = np.zeros((3,len(x)))
    delta = 1e-6
    for i in xrange(len(x)):
        xtmp = x + delta*np.squeeze(np.eye(N=1, M=len(x), k=i))
        ftmp = value(xtmp)
        jac_fd[:,i] = (ftmp - fval)/delta
    if np.max(np.max(np.abs(jac-jac_fd))) > 1e-6:
        print("Programmed Jacobian:")
        print(jac)
        print("Finite Difference Jacobian:")
        print(jac_fd)
    else:
        print("Programmed Jacobian and finite difference Jacobian within 1e-6")
    
    
    
    
    
    
