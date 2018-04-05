import sys, os
import numpy as np
from scipy.spatial import ConvexHull
from scipy.interpolate import griddata

# Global variables used for specifying KS and modified P-norm coefficients
KS_PARAM = 50.
PN_PARAM = 10.

# Kreselmeier-Steinhauser function
def ksFunction(x,p):
    return (1./p)*np.log(np.sum(np.exp(p*x)))

    
# Modified P-norm function
def pnFunction(x,p):
    return ((1./len(x))*np.sum(np.power(x,p)))**(1./p)
    
    
def assignTotalStress(k, filename, output='verbose'):
    """
    Read total stress data from Aero-S file and agglomerate.
    """
    
    # Load data
    try:
        data = np.loadtxt(filename,dtype=float,skiprows=3) # stresses in 4th column (0-indexed)
    except IOError as err:
        print('WARNING: could not open %s' % filename)
        print(err)
        return -2.
        
    stress = data[:,-1]
    
    # Agglomerate stresses and assign
    if 'MAX' in k:
        response = np.max(stress)
    elif 'KS' in k:
        stemp = np.mean(stress)
        response = ksFunction(stress/stemp,KS_PARAM)*stemp            
    elif 'PN' in k:
        stemp = np.mean(stress)
        response = pnFunction(stress/stemp,PN_PARAM)*stemp            
    else:
        raise ValueError("MAX, KS, or PN must be in QoI name: %s" % k)
        
    if output == 'verbose':
        print('%s assigned from Aero-S file %s' % (k,filename))

    return response
    
    
def assignFailureCriteria(nozzle, k, filesuffix, material, output='verbose'):
    """
    Read total stress data from Aero-S file, calculate failure criteria and 
    agglomerate.
    """

    # Assign failure criteria
    if not hasattr(material,'failureType'):    
        print('WARNING: no failure type for component %s' % k)
        return -2.
        
    elif material.failureType == 'VON_MISES':
    
        # Von Mises is read in through the STRESS files
        filename = ['STRESS.' + str(filesuffix)]
        try: 
            data = np.loadtxt(filename[0],dtype=float,skiprows=3) # stresses in 4th column (0-indexed)
        except IOError as err:
            print('WARNING: could not open %s' % filename[0])
            print(err)
            return -2.
            
        failureMeasure = data[:,-1]     
        failureLimit = material.yieldStress
                       
    elif material.failureType == 'PRINCIPLE_FAILURE_STRAIN':
                
        filename = ['STRAINP1.' + str(filesuffix), 'STRAINP3.' + str(filesuffix)]
        
        try:
            data = np.loadtxt(filename[0],dtype=float,skiprows=3) # strains in 4th column (0-indexed)
        except IOError as err:
            print('WARNING: could not open %s' % filename[0])
            print(err)
            return -2.    
        strainp1 = data[:,-1]

        try:
            data = np.loadtxt(filename[1],dtype=float,skiprows=3) # strains in 4th column (0-indexed)
        except IOError as err:
            print('WARNING: could not open %s' % filename[1])
            print(err)
            return -2.
        strainp3 = data[:,-1]
          
        # Assign failure criterion
        failureMeasure = np.max(np.vstack((strainp1,strainp3)),axis=0)
        failureLimit = material.getFailureLimit()
        
    elif material.failureType == 'LOCAL_FAILURE_STRAIN':

        filename = ['STRAINXX.' + str(filesuffix), 'STRAINYY.' + str(filesuffix)]
        
        try:
            data = np.loadtxt(filename[0],dtype=float,skiprows=3) # strains in 4th column (0-indexed)
        except IOError as err:
            print('WARNING: could not open %s' % filename[0])
            print(err)
            return -2.
        strainxx = data[:,-1]
        
        try:
            data = np.loadtxt(filename[1],dtype=float,skiprows=3) # strains in 4th column (0-indexed)
        except IOError as err:
            print('WARNING: could not open %s' % filename[1])
            print(err)
            return -2.  
        strainyy = data[:,-1]
        
        # Assign failure criterion
        failureStrain = material.getFailureLimit()
        failxx = np.empty((strainxx.size))
        failyy = np.empty((strainyy.size))
        for i in range(len(failxx)):
            if strainxx[i] >= 0.:
                failxx[i] = failureStrain[0]
            else:
                failxx[i] = failureStrain[1]
            if strainyy[i] >= 0.:
                failyy[i] = failureStrain[2]
            else:
                failyy[i] = failureStrain[3]
                
        failureMeasure = np.vstack((strainxx,strainyy))
        failureLimit = np.vstack((failxx,failyy))

    else:
        raise RuntimeError('Failure type %s not accepted.' % material.failureType)
    
    # Agglomerate stresses and assign
    if 'MAX' in k:
        response = np.max(failureMeasure/failureLimit)
    elif 'KS' in k:
        response = ksFunction(failureMeasure/failureLimit,KS_PARAM)                  
    elif 'PN' in k:
        response = pnFunction(failureMeasure/failureLimit,PN_PARAM)
    else: # pointwise stresses must be desired
        
        # Prepare failure ratio for each node
        failureRatio = failureMeasure/failureLimit # n x m
        frs = failureRatio.shape 
        if len(frs) > 1:
            n, m = frs # n = number of failure crit, m = number of nodes
        else:
            failureRatio = np.array([failureRatio])
            n, m = failureRatio.shape
            
        # Interpolate Radial Data on Convex Hull
        interpLoc = nozzle.qoi.getLocation(k)
        response = []
        for i in range(n):
            response.append(list(interpolateRadialDataOnConvexHull(nozzle, \
                       interpLoc, data[:,1:4],failureRatio[i,:].T)))
        response = np.squeeze(np.array(response))
    
    if output == 'verbose':
        for f in filename:
            print('%s assigned from Aero-S file %s' % (k,f))
    
    return response
    
    
def assignTemperature(k, filesuffix, output='verbose'):
    """
    Read temperature data from Aero-S file and agglomerate.
    """

    filename = 'TEMP.' + str(filesuffix)
    
    try:
        data = np.loadtxt(filename,dtype=float,skiprows=3) # temps in 4th column (0-indexed)
    except IOError as err:
        print('WARNING: could not open %s' % filename)
        print(err)
        return -2.
            
    # Agglomerate temperatures and assign
    if 'MAX' in k:
        response = np.max(data[:,-1])
    elif 'KS' in k:
        stemp = np.mean(data[:,-1])
        response = ksFunction(data[:,-1]/stemp,KS_PARAM)*stemp             
    elif 'PN' in k:
        stemp = np.mean(data[:,-1])
        response = pnFunction(data[:,-1]/stemp,PN_PARAM)*stemp         
    else:
        raise RuntimeError('MAX, KS, or PN must be in QoI name: %s' % k)
    
    if output == 'verbose':
        print('%s assigned from Aero-S file %s' % (k,filename))
    
    return response


def assignTempRatio(nozzle, k, filesuffix, material, output='verbose'):
    """
    Read temperature data from Aero-S file, calculate temperature ratio 
    (temperature divided by maximum allowable temperature) and agglomerate.
    """

    # Assign failure criteria
    if not hasattr(material,'Tmax'):    
        print('WARNING: no max temperature for component %s' % k)
        return -2.
        
    filename = 'TEMP.' + str(filesuffix)
    
    try:
        data = np.loadtxt(filename,dtype=float,skiprows=3) # temps in 4th column (0-indexed)
    except IOError as err:
        print('WARNING: could not open %s' % filename)
        print(err)
        return -2.
    
    # Agglomerate temperatures and assign
    if 'MAX' in k:
        response = np.max(data[:,-1])/material.Tmax
    elif 'KS' in k:
        response = ksFunction(data[:,-1]/material.Tmax,KS_PARAM)                   
    elif 'PN' in k:
        response = pnFunction(data[:,-1]/material.Tmax,PN_PARAM)           
    else: # pointwise temperature ratios must be desired
        
        # Prepare failure ratio for each node
        tempRatio = data[:,-1]/material.Tmax # n x m
        
        # Interpolate Radial Data on Convex Hull
        interpLoc = nozzle.qoi.getLocation(k)
        response = []
        response.append(list(interpolateRadialDataOnConvexHull(nozzle, \
                         interpLoc, data[:,1:4],tempRatio[:].T)))
        response = np.squeeze(np.array(response))
    
    if output == 'verbose':
        print('%s assigned from Aero-S file %s' % (k,filename))
    
    return response


def interpolateRadialDataOnConvexHull(nozzle, interpLoc, coord, data, output='verbose'):
    """
    Interpolate values at points specified by the intersection of radial rays
    defined in interpLoc array with convex hull defined by data. nozzle is used
    to provide information on centerlin defining data and provide a guess point.
    interpLoc is an n x 2 array depicting desired data at m points on the convex
    hull. An axial location (x) is specified, followed by an angle in degrees,
    depicting direction of ray in y-z plane as elevated from the x-y plane. Data
    is a Numpy array of m x 3 for m nodes containing x, y, and z-coordinates of 
    each node. data is an m x 1 array containing corresponding values for each 
    node. 
    """
    # Determine ray information for intersection with convex hull
    if( nozzle.dim == '3D' ): # reference from centerline
        print('WARNING: Pointwise stresses not implemented for 3D parameterization.')
        return [-1]*len(interpLoc)
    else:        
        x = interpLoc[:,0]
        r = nozzle.wall.geometry.radius(interpLoc[:,0])
        y = r*np.sin(interpLoc[:,1]*np.pi/180.)
        z = r*np.cos(interpLoc[:,1]*np.pi/180.)
        n = x.size
        rayDirections = np.vstack((np.zeros(n,),y,z))
        rayOrigins = np.vstack((x,np.zeros(n,),np.zeros(n,)))
    
    # Build convex hull of data
    hull = ConvexHull(coord) # convex hull of all nodes in 3D
    eq = hull.equations.T # transpose of hull equations
    V, b = eq[:-1], eq[-1] # normal vectors & offsets for hyperplanes

    # Calculate intersection of ray with convex hull
    denominator = np.dot(rayDirections.T,V)
    # XXX DO SOMETHING HERE TO GET RID OF DIVIDE BY ZERO ERROR
    alpha = (-b - np.dot(rayOrigins.T,V))/np.dot(rayDirections.T,V)
    #intersections = np.min(alpha[alpha>=0])*rayDirections + rayOrigins
    # Subtract a small number here to ensure point is inside convex hull within
    # rounding errors
    intersections = (np.min(alpha[alpha>=0])-1e-14)*rayDirections + rayOrigins
    
    # Now do interpolation for intersection points that have been found
    val = griddata(coord,data,intersections.T,method='linear')
    
    return val
   
    
def PostProcess(nozzle, runDir, output='verbose'):
    """
    Post-process Aero-S results after thermal/structural analysis. Assigns
    the following:
    'XXX_TOTAL_STRESS'
    'XXX_FAILURE_CRITERIA'
    'XXX_TEMPERATURE'
    'XXX_TEMP_RATIO'
    """

    cwd = os.getcwd()
    # Change to run directory if necessary
    if runDir:
        os.chdir(runDir)
    
    # ---- Indexing arrays (for use with specific material properties)
    mat = [nozzle.wall.layer[q].material for q in range(5)]
    mat.append(nozzle.stringers.material)
    for i in range(nozzle.baffles.n):
        mat.append(nozzle.baffles.material)
        
    # --- Determine labeling of files
    aerosSuffix = [0, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    multifPrefix = nozzle.prefixLabels
    # if output == 'verbose':
    #     print("Aero-S file labeling:")
    #     for i, j in zip(aerosSuffix, nozzle.prefixLabels):
    #         print(i, j)
        
    # --- Assign results as necessary
    for k in nozzle.qoi.names:
        
        # Assign total stress results
        if 'TOTAL_STRESS' in k:

            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to 
                    # calculate total stress in this component
                    filename = 'STRESS.' + \
                        str(aerosSuffix[multifPrefix.index(prefix)])
                    val = assignTotalStress(k, filename, output)
                    nozzle.qoi.setValue(k, val)
                    break # Go to next k in nozzle.responses
            else: # If k is not specific, provide stresses for all components
                for prefix in multifPrefix:
                    specificLabel = prefix + '_' + k
                    filename = 'STRESS.' + \
                        str(aerosSuffix[multifPrefix.index(prefix)])
                    val = assignTotalStress(k, filename, output)
                    nozzle.qoi.setValue(specificLabel, val)
        
        # Assign failure criteria results
        elif 'FAILURE_CRITERIA' in k:

            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to 
                    # calculate failure crit. in this component
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)]
                    val = assignFailureCriteria(nozzle, k, filesuffix, 
                        mat[multifPrefix.index(prefix)], output)
                    nozzle.qoi.setValue(k, val)
                    break # Go to next k in nozzle.responses 
            else: # If k is not specific, provide FC for all components
                for prefix in multifPrefix:
                    specificLabel = prefix + '_' + k
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)]
                    val = assignFailureCriteria(nozzle, k, filesuffix, 
                        mat[multifPrefix.index(prefix)], output)
                    nozzle.qoi.setValue(specificLabel, val)
        
        # Assign mechanical stress results
        elif 'MECHANICAL_STRESS' in k:
        
            raise RuntimeError('Output of MECHANICAL_STRESS is deprecated.')             
            
        # Assign thermal stress results
        elif 'THERMAL_STRESS' in k:
        
            raise RuntimeError('Output of THERMAL_STRESS is deprecated')

        # Assign temperature results
        elif 'TEMPERATURE' in k:

            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to 
                    # calculate temp ratio in this component
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)]
                    val = assignTemperature(k, filesuffix, output)
                    nozzle.qoi.setValue(k, val)
                    break # Go to next k in nozzle.responses
            else: # If k is not specific, provide temp ratio for all components
                # Only first four components have temperatures output
                for prefix in multifPrefix[0:4]: 
                    specificLabel = prefix + '_' + k
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)]
                    val = assignTemperature(k, filesuffix, output)
                    nozzle.qoi.setValue(specificLabel, val)
            
        # Assign temperature ratio results
        elif 'TEMP_RATIO' in k:

            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to 
                    # calculate temp ratio in this component
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)]
                    val = assignTempRatio(nozzle, k, filesuffix, 
                        mat[multifPrefix.index(prefix)], output)
                    nozzle.qoi.setValue(k, val)
                    break # Go to next k in nozzle.responses
            else: # If k is not specific, provide temp ratio for all components
                # Only first four components have temperatures output
                for prefix in multifPrefix[0:4]: 
                    specificLabel = prefix + '_' + k
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)]
                    val = assignTempRatio(k, filesuffix, 
                        mat[multifPrefix.index(prefix)], output)
                    nozzle.qoi.setValue(specificLabel, val)

    # Change to original directory
    if runDir:
        os.chdir(cwd)
    
    if output == 'verbose':
        print('Aero-S responses obtained')
    
    return 0                           

        
