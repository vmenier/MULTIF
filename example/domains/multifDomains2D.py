import numpy as np

from linearConstraints import bspline
from linearConstraints import cleanupConstraintMatrix

from domains import LinIneqDomain
from domains import UniformDomain
from domains import LogNormalDomain
from domains import ComboDomain

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


def buildDesignDomain(output='verbose'):

    lb_perc = 0.8
    ub_perc = 1.2

    # ============================================================================
    # Wall design variables for fixed inlet
    # ============================================================================
    WALL_COEFS= (0.0000, 0.0000, 0.1, 0.3, 0.7, 1.0, 1.3, 
                 1.3500, 1.3500, 1.4000, 1.5000, 1.6000, 1.8000, 2.3371, 2.3371, 
                 0.4395, 0.4395, 0.4395, 0.4, 0.34, 0.31, 0.27, 
                 0.2700, 0.2700, 0.2700, 0.3, 0.33, 0.38, 0.3955, 0.3955)
    WALL_COEFS_DV= (0, 0, 0, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 0, 0, 
                    0, 0, 0, 10, 11, 12, 13, 13, 13, 13, 14, 15, 16, 0, 0)                  
    # DV_LIST = WALL, 21, 
    # Inner wall, 21-dimensions
#    x_wall = np.array([0.2124, 0.2269, 0.2734, 0.3218, 0.3230, 0.3343, 0.3474, 0.4392, 0.4828, 0.5673, 0.6700, 0.3238, 0.2981, 0.2817, 0.2787, 0.2797, 0.2807, 0.2936, 0.2978, 0.3049, 0.3048])
#    _, A, b = wall(x_wall)
#    inner_wall_domain = LinIneqDomain(A, b, lb = lb_perc * x_wall, ub = ub_perc * x_wall, center = x_wall)
    x_wall = np.array([0.542184, 0.861924, 1.072944, 1.211161, 1.311161, 
                       1.408983, 1.528983, 1.723828, 2.086573, 0.417017, 
                       0.365097, 0.301792, 0.267426, 0.277426, 0.332508, 
                       0.385631])
    A, b = bspline(WALL_COEFS, WALL_COEFS_DV, 7, (-0.3,0.005,-0.025,0.35), 
                     xLimits=[None,2.3], delta=0.1, throatIsLowest=1, 
                     minThroat=0.2, output=output) 
    A, b = cleanupConstraintMatrix(Alist=[A],blist=[b])
    inner_wall_domain = LinIneqDomain(A, np.squeeze(b), center = x_wall)
                     
    #             THERMAL_LAYER, 6
    # Thermal Layer 6-dimensions (first 2 correspond to the x-position between 2nd and 3rd break)
    x_thermal = np.array([0.3, 0.6, 0.03, 0.03, 0.03, 0.03])
    _, A, b  = thermalLayer(x_thermal)
    lb = np.hstack((0.8*x_thermal[0:2], 0.01*np.ones(4)))
    ub = np.hstack((1.2*x_thermal[0:2], 0.05*np.ones(4)))
    thermal_layer_domain = LinIneqDomain(A, b, lb = lb, ub = ub, center = x_thermal)

    #           AIR_GAP_THICKNESS, 1,
    # Air gap, 1-dimension
    air_gap_domain  = UniformDomain(0.003, 0.01, center = 0.005)

    # Load layer geometry, 3* 6-dimensions

    #           LOAD_LAYER_INSIDE, 6,
    x_load_layer_inner = np.array([0.3, 0.6, 0.002, 0.002, 0.002, 0.002])
    _, A, b = thermalLayer(x_load_layer_inner)
    lb = np.hstack((0.8*x_load_layer_inner[0:2], 0.0005*np.ones(4)))
    ub = np.hstack((1.2*x_load_layer_inner[0:2], 0.01*np.ones(4)))
    load_layer_inner_domain = LinIneqDomain(A, b, lb = lb, ub = ub, center = x_load_layer_inner)

    #             LOAD_LAYER_MIDDLE, 6,
    x_load_layer_middle = np.array([0.3, 0.6, 0.013, 0.013, 0.013, 0.013])
    _, A, b = thermalLayer(x_load_layer_middle)
    lb = np.hstack((0.8*x_load_layer_middle[0:2], 0.0064*np.ones(4)))
    ub = np.hstack((1.2*x_load_layer_middle[0:2], 0.0159*np.ones(4)))
    load_layer_middle_domain = LinIneqDomain(A, b, lb = lb, ub = ub, center = x_load_layer_middle)

    #            LOAD_LAYER_OUTSIDE, 6, 
    x_load_layer_outside = np.array([0.3, 0.6, 0.002, 0.002, 0.002, 0.002])
    _, A, b = thermalLayer(x_load_layer_outside)
    lb = np.hstack((0.8*x_load_layer_outside[0:2], 0.0005*np.ones(4)))
    ub = np.hstack((1.2*x_load_layer_outside[0:2], 0.01*np.ones(4)))
    load_layer_outside_domain = LinIneqDomain(A, b, lb = lb, ub = ub, center = x_load_layer_outside)
    
    #            BAFFLES, 10, 
    x_baffles = np.array([0.2, 0.4, 0.6, 0.8, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
    _, A, b = baffles(x_baffles)
    lb = np.hstack((lb_perc * x_baffles[0:4], 0.0074*np.ones(6)))
    ub = np.hstack((ub_perc * x_baffles[0:4], 0.0359*np.ones(6)))
    baffles_domain = LinIneqDomain(A, b, lb = lb, ub = ub, center = x_baffles)

    #            STRINGERS_THICKNESS_VALUES, 6,
    x_stringers = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
    lb = 0.002*np.ones((6,))
    ub = 0.01*np.ones((6,))
    stringer_domain = UniformDomain(lb, ub, center = x_stringers)


    design_domains = [inner_wall_domain, thermal_layer_domain, air_gap_domain, load_layer_inner_domain,
            load_layer_middle_domain, load_layer_outside_domain, baffles_domain, stringer_domain] 

    return ComboDomain(design_domains)


def buildRandomDomain(output='verbose', clip = None):

    random_domains = [
    #            CMC_DENSITY, 1,
        LogNormalDomain(7.7803, 0.0182**2, clip = clip),
    #            CMC_ELASTIC_MODULUS, 1,
        LogNormalDomain(4.2047, 0.0551**2, scaling = 1e9, clip = clip),
    #            CMC_POISSON_RATIO, 1,
        UniformDomain(0.23, 0.43),
    #             CMC_THERMAL_CONDUCTIVITY, 1,
        UniformDomain(1.37, 1.45),
    #            CMC_THERMAL_EXPANSION_COEF, 1, 
        UniformDomain(0.228e-6, 0.252e-6),     
    #            CMC_PRINCIPLE_FAILURE_STRAIN, 1, 
        LogNormalDomain(-2.6694, 0.1421**2,clip = clip),
    #            CMC_MAX_SERVICE_TEMPERATURE, 1, 
        UniformDomain(963, 983),
    #
    #
    #            GR-BMI_DENSITY, 1, 
        UniformDomain(1563, 1573), 
    #            GR-BMI_ELASTIC_MODULUS, 2,
        UniformDomain(57e9, 63e9),
        UniformDomain(57e9, 63e9),
    #            GR-BMI_SHEAR_MODULUS, 1, 
        UniformDomain(22.6e9, 24.0e9),
    #            GR-BMI_POISSON_RATIO, 1,
        UniformDomain(0.334, 0.354), 
    #            GR-BMI_MUTUAL_INFLUENCE_COEFS, 2, 
        UniformDomain(-0.1, 0.1),
        UniformDomain(-0.1, 0.1),
    #            GR-BMI_THERMAL_CONDUCTIVITY, 3,
        UniformDomain(3.208, 3.546),
        UniformDomain(3.208, 3.546),
        UniformDomain(3.243, 3.585),
    #            GR-BMI_THERMAL_EXPANSION_COEF, 3,
        UniformDomain(1.16e-6, 1.24e-6), 
        UniformDomain(1.16e-6, 1.24e-6), 
        UniformDomain(-0.04e-6, 0.04e-6),
    #            GR-BMI_LOCAL_FAILURE_STRAIN, 5,
        UniformDomain(0.675e-2, 0.825e-2, center = 0.75e-2),
        UniformDomain(-0.572e-2, -0.494e-2, center = -0.52e-2),
        UniformDomain(0.675e-2, 0.825e-2, center = 0.75e-2),
        UniformDomain(-0.572e-2, -0.494e-2, center = -0.52e-2),
        UniformDomain(0.153e-2, 0.187e-2, center = 0.17e-2),
    #            GR-BMI_MAX_SERVICE_TEMPERATURE, 1,
        UniformDomain(500, 510),
    #
    #
    #            TI-HC_DENSITY, 1, 
        UniformDomain(177.77, 181.37),
    #            TI-HC_ELASTIC_MODULUS, 1, 
        LogNormalDomain(0.6441, 0.0779**2, scaling = 1e9, clip = clip),
    #            TI-HC_POISSON_RATIO, 1, 
        UniformDomain(0.160, 0.196),
    #            TI-HC_THERMAL_CONDUCTIVITY, 1, 
        UniformDomain(0.680, 0.736),
    #            TI-HC_THERMAL_EXPANSION_COEF, 1, 
        UniformDomain(2.88e-6, 3.06e-6),
    #            TI-HC_YIELD_STRESS, 1,
        LogNormalDomain(2.5500, 0.1205**2, scaling = 1e6, clip = clip),
    #            TI-HC_MAX_SERVICE_TEMPERATURE, 1, 
        UniformDomain(745, 765),
    #
    #
    #            AIR_THERMAL_CONDUCTIVITY, 1, 
        UniformDomain(0.0320, 0.0530),
    #            PANEL_YIELD_STRESS, 1, 
        LogNormalDomain(5.7736, 0.1196**2, scaling = 1e6, clip = clip), 
    #            INLET_PSTAG, 1, 
        LogNormalDomain(11.5010, 0.0579**2, clip = clip),
    #            INLET_TSTAG, 1, 
        LogNormalDomain(6.8615, 0.0119**2, clip = clip),
    #            ATM_PRES, 1, 
        LogNormalDomain(9.8386, 0.0323**2, clip = clip),
    #            ATM_TEMP, 1, 
        LogNormalDomain(5.3781, 0.0282**2, clip = clip),
    #            HEAT_XFER_COEF_TO_ENV, 1
        LogNormalDomain(2.5090, 0.2285, clip = clip),
    ]

    return ComboDomain(random_domains)


