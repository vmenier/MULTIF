"""
Generate random samples from linear constraints of nozzle geometry and 
random parameters in nozzle problem. Also, generate the files necessary to run
random parameter sweeps in MULTI-F.

By changing the sampling controls below, any number of random samples, and 
sweeps with any number of points can be generated from the linear constraints.
These are then saved to the specified sample filename and sweep filename.

param = '3D'; # '2D' or '3D' nozzle parameterization
N = 100; # number of samples
P = 10; # number of sweeps
S = 50; # number of points in each sweep
samplefilename = '3d_samples.dat'; # filename to save samples in
sweepfilenameprefix = '3d_sweep'; # prefix of filename to save sweep points in
dN = 1000; # save every dN samples to the above samples file
dNplot = 10000; # plot every dNplot sampled geometries (avoid plot clutter)
output = 'verbose'; # print notifications to screen
plot = 'yes'; # plot a selection of sampled nozzle geometries

Rick Fenrich 9/28/17
Based on domain code supplied by Jeffrey Hokanson
"""

import sys
import numpy as np

from matplotlib import pyplot as plt

sys.path.append('../..');
from multif.nozzle import geometry

import multifDomains3D
import multifDomains2D

from domains import ComboDomain
from util import find_feasible_boundary

# ============================================================================
# Controls for sampling and sweeps
# ============================================================================
param = '3D'; # nozzle parameterization, either 2D or 3D
N = 10000; # number of samples
P = 0; # number of sweeps
S = 10; # number of points in each sweep
samplefilename = '2d_samples.dat'; # filename to save samples in
sweepfilenameprefix = '2d_sweep'; # prefix of filename to save sweep points in
dN = 100; # save every dN samples to the above samples file
dNplot = 100; # plot every dNplot sampled geometries (avoid plot clutter)
output = 'quiet'; # 'verbose' prints notifications to screen
plot = 'yes'; # plot sampled nozzle geometries

# ============================================================================
# Setup domains, linear constraints, and bounds
# ============================================================================

print('Parameterization is %s' % param);

if( param == '3D' ):
    designDomain = multifDomains3D.buildDesignDomain(output=output);
    randomDomain = multifDomains3D.buildRandomDomain(output=output);
else:
    designDomain = multifDomains2D.buildDesignDomain(output=output);
    randomDomain = multifDomains2D.buildRandomDomain(output=output);
fullDomain = ComboDomain([designDomain,randomDomain]);

A = fullDomain.A;
b = fullDomain.b;

# To save linear constraints
np.savetxt('linear_constraint_matrix.txt',A);
np.savetxt('linear_constraint_rhs.txt',b);
print('Linear constraint matrix saved to linear_constraint_matrix.txt.');
print('Linear constraint RHS saved to linear_constraint_rhs.txt.');
print;

# ============================================================================
# Randomly sample from constraints and bounds
# ============================================================================

# Sample directly from constraints in design domain
Z1 = designDomain.sample(draw=N);
Z2 = randomDomain.sample(draw=N);
Z = np.hstack((Z1,Z2));
print('Constraint matrix sampled %i times.' % N);
print('Random parameters sampled %i times.' % N);

for i in range(N):
    if(np.abs(Z[i,0]) < 1e-4):
        if(np.abs(Z[i,7]) < 1e-4):
            print i, Z[i,0], Z[i,7]

# To save samples
np.savetxt(samplefilename,Z[0:N:dN,:]);
print('Samples saved to %s.' % samplefilename);
print;

# ============================================================================
# Generate data for sweeps if desired
# ============================================================================
for i in range(P):

    # Randomly select 2 random samples and draw a line between them
    ind1 = np.random.randint(0,high=N-1);
    ind2 = np.random.randint(0,high=N-1);
    v = Z[ind2,:]-Z[ind1,:];
    
    # Find ends of feasible range on this line
    p2 = find_feasible_boundary(Z[ind2,:], v, 1e-6, 0., 1., A, b); # point 2
    p1 = find_feasible_boundary(Z[ind1,:],-v, 1e-6, 0., 1., A, b); # point 1
    d = np.linalg.norm(p2-p1,ord=2); # distance b/w point 1 and 2
    
    # Obtain samples along sweep line
    sweepsamples = np.zeros((S,len(p1)));
    for j in range(len(p1)):
        sweepsamples[:,j] = np.linspace(p1[j],p2[j],num=S);
    
    sweepfilename = sweepfilenameprefix + '_' + str(i) + '.dat';
    np.savetxt(sweepfilename,sweepsamples);
    print('Sweep %i samples saved to %s.' % (i,sweepfilename));
print;

# ============================================================================
# Plot geometries if desired
# ============================================================================
if( plot == 'yes' ):
    
    if( param == '3D' ):

        # ============================================================================
        # Wall design variables for free inlet
        # ============================================================================
        # # Centerline
        # WALL_COEFS1 = (0.0000, 0.0000, 0.3000, 0.5750, 1.1477, 1.1500, 1.1500, 1.1523, 1.7262, 2.0000, 2.3000, 2.3000, 
        #                0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)
        # WALL_COEFS1_DV= (1,    1,      2,      3,      4,      5,      5,      6,      7,      0,      0,      0,     
        #                  8,    8,      8,      9,      10,     10,     10,     10,     11,     0,      0,      0)
        
        # # Major Axis
        # WALL_COEFS2= (0.0000, 0.0000, 0.3000, 0.5000, 0.7000, 0.9000, 1.1477, 1.1500, 
        #               1.1500, 1.1523, 1.4000, 1.6500, 1.9000, 2.1000, 2.3000, 2.3000, 
        #               0.3255, 0.3255, 0.3255, 0.3195, 0.3046, 0.2971, 0.2956, 0.2956, 
        #               0.2956, 0.2956, 0.3065, 0.3283, 0.3611, 0.4211, 0.4265, 0.4265)
        # WALL_COEFS2_DV= (1,   1,      2,      12,     13,     14,     15,     5,
        #                  5,   16,     17,     18,     19,     20,     0,      0,
        #                  0,   0,      0,      21,     22,     23,     24,     24,
        #                  24,  24,     25,     26,     27,     28,     0,      0)
                         
        # # Minor axis
        # WALL_COEFS3= (0.0000, 0.0000, 0.3000, 0.5500, 0.9000, 
        #               1.1500, 1.8000, 2.1000, 2.3000, 2.3000, 
        #               0.3255, 0.3255, 0.3255, 0.3195, 0.2956, 
        #               0.2750, 0.2338, 0.2167, 0.2133, 0.2133)
        # WALL_COEFS3_DV= (1,   1,      2,      29,     30,
        #                  31,  32,     33,     0,      0, 
        #                  0,   0,      0,      34,     35,  
        #                  36,  37,     38,     0,      0)

        # ============================================================================
        # Wall design variables for fixed inlet
        # ============================================================================
        # Centerline
        WALL_COEFS1 = (0.0000,   0.0000,   0.3000,   0.5750, 1.1500, 1.7262, 2.0000, 2.33702, 2.33702, 
                    0.099908, 0.099908, 0.099908, 0.12,   0.14,   0.17,   0.19,   0.19,    0.19)
        WALL_COEFS1_DV= (0,      0,        0,        1,      2,      3,      0,      0,      0,     
                        0,      0,        0,        4,      5,      6,      0,      0,      0)
        
        # Major Axis
        WALL_COEFS2= (0.0000,   0.0000,   0.3000,   0.7000, 1.1500, 1.6000, 1.8,     2.33702, 2.33702, 
                    0.439461, 0.439461, 0.439461, 0.6,    0.7,    0.8,    0.85,     0.92,    0.92)
        WALL_COEFS2_DV= (0,     0,        0,        7,      8,      9,      10,      0,       0,
                        0,     0,        0,        11,     12,     13,     14,      0,       0)
                        
        # Minor axis
        WALL_COEFS3= (0.0000,   0.0000,   0.3000,   0.7000, 1.1500, 1.6000,   2.33702, 2.33702, 
                    0.439461, 0.439461, 0.439461, 0.3,    0.29,   0.26,     0.24,    0.24)
        WALL_COEFS3_DV= (0,     0,        0,        7,       8,     15,       0,       0,  
                        0,     0,        0,        16,      17,    18,       0,      0)     
   
        # Centerline
        plt.figure()
        for i in range(0,N,dNplot):
            # Build coefficient matrix
            coefs = list(WALL_COEFS1);
            for j in range(len(coefs)):
                if( WALL_COEFS1_DV[j] > 0 ):
                    coefs[j] = Z[i,WALL_COEFS1_DV[j]-1];
            spline = geometry.Bspline(coefs);
            xplot = np.linspace(min(coefs[0:len(coefs)/2]),max(coefs[0:len(coefs)/2]),200);
            plt.plot(xplot,spline.radius(xplot)); # B-spline
            #plt.plot(coefs[0:len(coefs)/2],coefs[len(coefs)/2:],'-o'); # bounding box
        plt.axis('equal');
        plt.grid();
        plt.title('Centerline');
        
        # Major axis
        plt.figure()
        for i in range(0,N,dNplot):
            # Build coefficient matrix
            coefs = list(WALL_COEFS2);
            for j in range(len(coefs)):
                if( WALL_COEFS2_DV[j] > 0 ):
                    coefs[j] = Z[i,WALL_COEFS2_DV[j]-1];
            spline = geometry.Bspline(coefs);
            xplot = np.linspace(min(coefs[0:len(coefs)/2]),max(coefs[0:len(coefs)/2]),200);
            plt.plot(xplot,spline.radius(xplot)); # B-spline
            #plt.plot(coefs[0:len(coefs)/2],coefs[len(coefs)/2:],'-o'); # bounding box
        plt.axis('equal');
        plt.grid();
        plt.title('Major axis');
            
        # Minor axis
        plt.figure()
        for i in range(0,N,dNplot):
            # Build coefficient matrix
            coefs = list(WALL_COEFS3);
            for j in range(len(coefs)):
                if( WALL_COEFS3_DV[j] > 0 ):
                    coefs[j] = Z[i,WALL_COEFS3_DV[j]-1];
            spline = geometry.Bspline(coefs);
            xplot = np.linspace(min(coefs[0:len(coefs)/2]),max(coefs[0:len(coefs)/2]),200);
            plt.plot(xplot,spline.radius(xplot)); # B-spline
            #plt.plot(coefs[0:len(coefs)/2],coefs[len(coefs)/2:],'-o'); # bounding box
        plt.axis('equal');
        plt.grid();
        plt.title('Minor axis');

        plt.show();

        print('Selected samples of major and minor axes and centerlines plotted.');
        
    else:

        # ============================================================================
        # Wall design variables for fixed inlet
        # ============================================================================        
        WALL_COEFS= (0.0000, 0.0000, 0.1, 0.3, 0.7, 1.0, 1.3, 
                     1.3500, 1.3500, 1.4000, 1.5000, 1.6000, 1.8000, 2.3371, 2.3371, 
                     0.4395, 0.4395, 0.4395, 0.4, 0.34, 0.31, 0.27, 
                     0.2700, 0.2700, 0.2700, 0.3, 0.33, 0.38, 0.3955, 0.3955)
        WALL_COEFS_DV= (0, 0, 0, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 0, 0, 
                        0, 0, 0, 10, 11, 12, 13, 13, 13, 13, 14, 15, 16, 0, 0)
                            
        # Nozzle inner wall
        plt.figure()
        for i in range(0,N,dNplot):
            # Build coefficient matrix
            coefs = list(WALL_COEFS);
            for j in range(len(coefs)):
                if( WALL_COEFS_DV[j] > 0 ):
                    coefs[j] = Z[i,WALL_COEFS_DV[j]-1];
            spline = geometry.Bspline(coefs);
            xplot = np.linspace(min(coefs[0:len(coefs)/2]),max(coefs[0:len(coefs)/2]),200);
            plt.plot(xplot,spline.radius(xplot)); # B-spline
            #plt.plot(coefs[0:len(coefs)/2],coefs[len(coefs)/2:],'-o'); # bounding box
        plt.axis('equal');
        plt.grid();
        plt.title('Inner wall');
        plt.show();
        
        print('Selected samples of nozzle inner wall shape plotted.');
    
