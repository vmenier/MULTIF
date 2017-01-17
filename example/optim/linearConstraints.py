import numpy as np

# 55 linear constraints for a B-spline nozzle geometry parameterized with 21 
# degrees of freedom. The function input x is a 21-D vector containing the 
# dofs of control points in the nozzle shape, ordered inlet to outlet, x-coord
# first, and then y-coord.
def wall(x):
    Rinlet = 0.3255
    Xinlet = 0.19
    delta = 1e-3 # for x-proximity for 2 control points controlling throat
    delta2 = 1e-3 # for proximity to throat
    delta3 = 1e-2 # for x-proximity of control points (except at throat)
    mMin = -1.8
    mMax = 0.5
    nLinearCon = 44 + 1 + 4 + 2 + 8
    A = np.zeros((nLinearCon, x.size)) # matrix for Ax <= b
    b = np.zeros(nLinearCon) # RHS vector for Ax <= b
    
    conNum = 0
    
    # Set constraints on x-position
    A[0,0] = -1; b[0] = -Xinlet - delta3; conNum += 1
    for i in range(conNum,conNum+10):
      A[i,i] = -1; A[i,i-1] = 1;
      if conNum == 4: # use different delta for at throat
        b[i] = -delta
      else:
        b[i] = -delta3
      conNum += 1
    
    # Set constraints on convergence of pre-throat area
    j = 11
    for i in range(conNum,conNum+4):
      A[i,j] = 1; b[i] = Rinlet; conNum += 1; j += 1
    
    # Set constraints making throat the lowest point
    j = 11
    for i in range(conNum,conNum+3):
      A[i,14] = 1; A[i,j] = -1; b[i] = -delta2; conNum += 1; j += 1
    j = 15
    for i in range(conNum,conNum+6):
      A[i,14] = 1; A[i,j] = -1; b[i] = -delta2; conNum +=1; j += 1
    
    # Set constraints for steepness of slopes
    A[conNum,0] = mMin; A[conNum,11] = -1; b[conNum] = mMin*Xinlet - Rinlet; conNum += 1
    A[conNum,11] = 1; A[conNum,0] = -mMax; b[conNum] = Rinlet - mMax*Xinlet; conNum += 1
    j = 1
    for i in range(conNum,conNum+3):
      A[i,j] = mMin; A[i,j-1] = -mMin; A[i,j+11] = -1; A[i,j+10] = 1; conNum += 1; j += 1
    j = 1
    for i in range(conNum,conNum+3):
      A[i,j+11] = 1; A[i,j+10] = -1; A[i,j] = -mMax; A[i,j-1] = mMax; conNum += 1; j += 1
    j = 5
    for i in range(conNum,conNum+6):
      A[i,j] = mMin; A[i,j-1] = -mMin; A[i,j+10] = -1; A[i,j+9] = 1; conNum += 1; j += 1
    j = 5
    for i in range(conNum,conNum+6):
      A[i,j+10] = 1; A[i,j+9] = -1; A[i,j] = -mMax; A[i,j-1] = mMax; conNum += 1; j += 1

    # Set upper bound constraints on slope in pre-throat region for all segments except first and last
    slopeCondition = [0.05, 0.05]
    j = 1
    c = 0
    for i in range(conNum,conNum+2):
        A[i,j] = -slopeCondition[c]; A[i,j-1] = slopeCondition[c]; A[i,j+11] = 1; A[i,j+10] = -1; conNum += 1; c += 1; j += 1;
    
    # Set lower bound on slope in segment just prior to throat
    slopeCondition = -0.2
    i = conNum
    j = 3
    A[i,j] = slopeCondition; A[i,j-1] = -slopeCondition; A[i,j+11] = -1; A[i,j+10] = 1; conNum += 1; j += 1
    
    # Set upper bound constraints on slope in post-throat region for all segments following throat
    slopeCondition = [0.16, 0.2, 0.3, 0.2, 0.1, 0.05]
    j = 5
    c = 0
    for i in range(conNum,conNum+6):
        A[i,j] = -slopeCondition[c]; A[i,j-1] = slopeCondition[c]; A[i,j+10] = 1; A[i,j+9] = -1; conNum += 1; c += 1; j += 1;

    # Set lower bound constraints on slope in post-throat region for all segments following post-throat segment
    slopeCondition = [0.04, 0.04, 0.07, 0., -0.05, -0.05]
    j = 5
    c = 0
    for i in range(conNum,conNum+6):
        A[i,j] = slopeCondition[c]; A[i,j-1] = -slopeCondition[c]; A[i,j+10] = -1; A[i,j+9] = 1; conNum += 1; c += 1; j += 1;
    return (A.dot(x) - b, A, b)
 
# linear constraints for a piecewise-linear shape with n nodes and 2*n - 2 dofs
# where the x-position of the first and last nodes is fixed at 0 and 1
def thermalLayer(x):
    delta = 1e-3 # for x-proximity of control points
    mMin = -0.2
    mMax = 0.2
    
    n = int((len(x) + 2)/2) # assume all node x-positions except first and last, and all node thicknesses are given
    nx = n - 2 # bumber of node x-positions
    nt = n # number of node thicknesses (y-values)
    nLinearCon = (n-3) + 2*(n-1)
    A = np.zeros((nLinearCon, x.size)) # matrix for Ax <= b
    b = np.zeros(nLinearCon) # RHS vector for Ax <= b
    
    conNum = 0
    
    # Set constraint(s) on x-position
    for i in range(n-3):
        A[conNum,i] = 1; A[conNum,i+1] = -1; b[conNum] = -delta; conNum += 1
        
    # Set constraints for steepness of slopes
    # First segment
    A[conNum,0] = mMin; A[conNum,nx] = 1; A[conNum,nx+1] = -1; b[conNum] = 0; conNum += 1
    A[conNum,0] = -mMax; A[conNum,nx] = -1; A[conNum,nx+1] = 1; b[conNum] = 0; conNum += 1
    
    # In-between segments
    for i in range(n-1-2):
        A[conNum,i] = -mMin; A[conNum,1+i] = mMin; A[conNum,nx+i+1] = 1; A[conNum,nx+i+2] = -1; b[conNum] = 0; conNum += 1
        A[conNum,i] = mMax; A[conNum,1+i] = -mMax; A[conNum,nx+i+1] = -1; A[conNum,nx+i+2] = 1; b[conNum] = 0; conNum += 1
    
    # Last segment
    A[conNum,nx-1] = -mMin; A[conNum,nx+nt-2] = 1; A[conNum,nx+nt-1] = -1; b[conNum] = -mMin; conNum += 1
    A[conNum,nx-1] = mMax; A[conNum,nx+nt-2] = -1; A[conNum,nx+nt-1] = 1; b[conNum] = mMax; conNum += 1
    
    return (A.dot(x) - b, A, b)        
    
def baffles(x):
    minDistance = 0.15 # non-dimensional minimum separation distance
    maxDistance = 0.25 # non-dimensional maximum separation distance
    nLinearCon = 10
    A = np.zeros((nLinearCon, x.size)) # matrix for Ax <= b
    b = np.zeros(nLinearCon) # RHS vector for Ax <= b
    
    conNum = 0     
    # Set constraints on x-position
    A[conNum,0] = 1; b[conNum] = maxDistance; conNum += 1
    A[conNum,0] = -1; b[conNum] = -minDistance; conNum += 1
    j = 1
    for i in range(conNum,conNum+3):
        A[i,j] = 1; A[i,j-1] = -1; b[conNum] = maxDistance; conNum += 1; j += 1
    j = 1
    for i in range(conNum,conNum+3):
        A[i,j] = -1; A[i,j-1] = 1; b[conNum] = -minDistance; conNum += 1; j += 1
    A[conNum,3] = -1; b[conNum] = maxDistance - 1; conNum += 1
    A[conNum,3] = 1; b[conNum] = 1 - minDistance; conNum += 1
    
    return (A.dot(x) - b, A, b)
    
      
