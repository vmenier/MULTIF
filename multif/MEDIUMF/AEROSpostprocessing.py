import os, sys
import numpy as np

# Kreselmeier-Steinhauser function
def ksFunction(x,p):
    return (1./p)*np.log(np.sum(np.exp(p*x)));

    
# Modified P-norm function
def pnFunction(x,p):
    return ((1./len(x))*np.sum(np.power(x,p)))**(1./p)
    
    
def assignTotalStress(k, filename, output='verbose'):
 
    # KS and PN params for stresses
    ks_param = 50.;
    pn_param = 10.;
    
    # Load data
    try:
        data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
    except:
        sys.stdout.write('WARNING: could not open %s\n' % filename);
        return 0;
        
    stress = data[:,-1];
    
    # Agglomerate stresses and assign
    if 'MAX' in k:
        response = np.max(stress);
    elif 'KS' in k:
        stemp = np.mean(stress);
        response = ksFunction(stress/stemp,ks_param)*stemp;                   
    elif 'PN' in k:
        stemp = np.mean(stress);
        response = pnFunction(stress/stemp,pn_param)*stemp;               
    else:
        sys.stderr.write('  ## ERROR : MAX, KS, or PN must be in QoI name: %s' % k);
        sys.exit(1);
        
    if output == 'verbose':
        sys.stdout.write('%s assigned from Aero-S file %s\n' % (k,filename));

    return response;
    
    
def assignFailureCriteria(k, filesuffix, material, output='verbose'):
        
    # KS and PN params for failure criteria
    ks_param = 50.;
    pn_param = 10.;
    
    # Assign failure criteria
    if not hasattr(material,'failureType'):
    
        sys.stdout.write('WARNING: no failure type for component %i\n' % index);
        return 0;
        
    elif material.failureType == 'VON_MISES':
    
        # Von Mises is read in through the STRESS files
        filename = ['STRESS.' + str(filesuffix)];
        try: 
            data = np.loadtxt(filename[0],dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
        except IOError:
            sys.stdout.write('WARNING: could not open %s\n' % filename[0]);
            return 0;
            
        failureMeasure = data[:,-1];        
        failureLimit = material.yieldStress; 
                       
    elif material.failureType == 'PRINCIPLE_FAILURE_STRAIN':
                
        filename = ['STRAINP1.' + str(filesuffix), 'STRAINP3.' + str(filesuffix)];
        
        try:
            data = np.loadtxt(filename[0],dtype=float,skiprows=3); # strains in 4th column (0-indexed)
        except IOError:
            sys.stdout.write('WARNING: could not open %s\n' % filename[0]);
            return 0;        
        strainp1 = data[:,-1];   

        try:
            data = np.loadtxt(filename[1],dtype=float,skiprows=3); # strains in 4th column (0-indexed)
        except IOError:
            sys.stdout.write('WARNING: could not open %s\n' % filename[1]);
            return 0;        
        strainp3 = data[:,-1]; 
          
        # Assign failure criterion
        failureMeasure = np.max(np.vstack((strainp1,strainp3)),axis=0);
        failureLimit = material.getFailureLimit();
        
    elif material.failureType == 'LOCAL_FAILURE_STRAIN':

        filename = ['STRAINXX.' + str(filesuffix), 'STRAINYY.' + str(filesuffix)];
        
        try:
            data = np.loadtxt(filename[0],dtype=float,skiprows=3); # strains in 4th column (0-indexed)
        except IOError:
            sys.stdout.write('WARNING: could not open %s\n' % filename[0]);
            return 0;        
        strainxx = data[:,-1];
        
        try:
            data = np.loadtxt(filename[1],dtype=float,skiprows=3); # strains in 4th column (0-indexed)
        except IOError:
            sys.stdout.write('WARNING: could not open %s\n' % filename[1]);
            return 0;            
        strainyy = data[:,-1];
        
        # Assign failure criterion
        failureStrain = material.getFailureLimit();
        failxx = np.empty((strainxx.size));
        failyy = np.empty((strainyy.size));
        for i in range(len(failxx)):
            if strainxx[i] >= 0:
                failxx = failureStrain[0];
            else:
                failxx = failureStrain[1];
            if strainyy[i] >= 0:
                failyy = failureStrain[2];
            else:
                failyy = failureStrain[3];
                
        failureMeasure = np.vstack((strainxx,strainyy));
        failureLimit = np.vstack((failxx,failyy));

    else:
        sys.stderr.write('\n ## ERROR: failure type %s not accepted.\n\n' % material.failureType);
        sys.exit(0);    
    
    # Agglomerate stresses and assign
    if 'MAX' in k:
        response = np.max(failureMeasure/failureLimit);
    elif 'KS' in k:
        response = ksFunction(failureMeasure/failureLimit,ks_param);                   
    elif 'PN' in k:
        response = pnFunction(failureMeasure/failureLimit,pn_param);               
    else:
        sys.stderr.write('  ## ERROR : MAX, KS, or PN must be in QoI name: %s' % k);
        sys.exit(1);
    
    if output == 'verbose':
        for f in filename:
            sys.stdout.write('%s assigned from Aero-S file %s\n' % (k,f));
    
    return response;
    
    
def assignTemperature(k, filesuffix, output='verbose'):
        
    # KS and PN params for failure criteria
    ks_param = 50.;
    pn_param = 10.;
    
    filename = 'TEMP.' + str(filesuffix);
    
    try:
        data = np.loadtxt(filename,dtype=float,skiprows=3); # temps in 4th column (0-indexed)
    except IOError:
        sys.stdout.write('WARNING: could not open %s\n' % filename);
        return 0;
            
    # Agglomerate temperatures and assign
    if 'MAX' in k:
        response = np.max(data[:,-1]);
    elif 'KS' in k:
        stemp = np.mean(data[:,-1]);
        response = ksFunction(data[:,-1]/stemp,ks_param)*stemp;                   
    elif 'PN' in k:
        stemp = np.mean(data[:,-1]);
        response = pnFunction(data[:,-1]/stemp,pn_param)*stemp;            
    else:
        sys.stderr.write('  ## ERROR : MAX, KS, or PN must be in QoI name: %s' % k);
        sys.exit(1);
    
    if output == 'verbose':
        sys.stdout.write('%s assigned from Aero-S file %s\n' % (k,filename));
    
    return response;


def assignTempRatio(k, filesuffix, material, output='verbose'):
        
    # KS and PN params for failure criteria
    ks_param = 50.;
    pn_param = 10.;
    
    # Assign failure criteria
    if not hasattr(material,'Tmax'):    
        sys.stdout.write('WARNING: no max temperature for component %i\n' % index);
        return 0;
        
    filename = 'TEMP.' + str(filesuffix);
    
    try:
        data = np.loadtxt(filename,dtype=float,skiprows=3); # temps in 4th column (0-indexed)
    except IOError:
        sys.stdout.write('WARNING: could not open %s\n' % filename);
        return 0;
    
    # Agglomerate temperatures and assign
    if 'MAX' in k:
        response = np.max(data[:,-1])/material.Tmax;
    elif 'KS' in k:
        response = ksFunction(data[:,-1]/material.Tmax,ks_param);                   
    elif 'PN' in k:
        response = pnFunction(data[:,-1]/material.Tmax,pn_param);            
    else:
        sys.stderr.write('  ## ERROR : MAX, KS, or PN must be in QoI name: %s' % k);
        sys.exit(1);
    
    if output == 'verbose':
        sys.stdout.write('%s assigned from Aero-S file %s\n' % (k,filename));
    
    return response;
    
    
def PostProcess ( nozzle, output='verbose' ):
    
    # ---- KS and modified P-norm parameters
    ks_param = 50.;
    pn_param = 10.;
    
    # ---- Indexing arrays (for use with specific material properties)
    mat = [nozzle.wall.layer[q].material for q in range(5)];
    mat.append(nozzle.stringers.material);
    for i in range(nozzle.baffles.n):
        mat.append(nozzle.baffles.material);
        
    # --- Determine labeling of files
    aerosSuffix = [0, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
    multifPrefix = nozzle.prefixLabels;   
        
    # --- Assign results as necessary
    for k in nozzle.responses:
        
        # Assign total stress results
        if 'TOTAL_STRESS' in k:
        
            assigned = 0;
            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to calculate total stress 
                    # in this component
                    filename = 'STRESS.' + str(aerosSuffix[multifPrefix.index(prefix)]);
                    nozzle.responses[k] = assignTotalStress(k,filename,output);
                    assigned = 1;
                    break; # Go to next k in nozzle.responses
            
            # If k is not specific, provide stresses for all components
            if assigned == 0:
                for prefix in multifPrefix:
                    specificLabel = prefix + '_' + k;
                    filename = 'STRESS.' + str(aerosSuffix[multifPrefix.index(prefix)]);
                    nozzle.responses[k].append(assignTotalStress(k,filename,output));
        
        # Assign failure criteria results
        elif 'FAILURE_CRITERIA' in k:
        
            assigned = 0;
            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to calculate failure crit.
                    # in this component
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)];
                    nozzle.responses[k] = assignFailureCriteria(k,filesuffix,mat[multifPrefix.index(prefix)],output);
                    assigned = 1;
                    break; # Go to next k in nozzle.responses
            
            # If k is not specific, provide failure criteria for all components
            if assigned == 0:
                for prefix in multifPrefix:
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)];
                    nozzle.responses[k].append(assignFailureCriteria(k,filesuffix,mat[multifPrefix.index(prefix)],output));
        
        # Assign mechanical stress results
        elif 'MECHANICAL_STRESS' in k:
        
            sys.stderr.write('  ## ERROR: Output of MECHANICAL_STRESS is deprecated.\n\n')
            sys.exit(1);  
            
            # Example:
            #filename = 'MECHANICAL_STRESS.0';
            #data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
            #stemp = np.mean(data[:,-1]);
            #nozzle.max_mechanical_stress[0] = np.max(data[:,-1]);                  
            
        # Assign thermal stress results
        elif 'THERMAL_STRESS' in k:
        
            sys.stderr.write('  ## ERROR: Output of THERMAL_STRESS is deprecated\n\n')
            sys.exit(1);
            
            # Example:
            #filename = 'THERMAL_STRESS.1';
            #data = np.loadtxt(filename,dtype=float,skiprows=3); # stresses in 4th column (0-indexed)
            #stemp = np.mean(data[:,-1]);
            #nozzle.max_thermal_stress[2] = np.max(data[:,-1]);

        # Assign temperature results
        elif 'TEMPERATURE' in k:
        
            assigned = 0;
            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to calculate temp ratio
                    # in this component
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)];
                    nozzle.responses[k] = assignTemperature(k,filesuffix,output);
                    assigned = 1;
                    break; # Go to next k in nozzle.responses
            
            # If k is not specific, provide temp ratio for all components
            if assigned == 0:
                for prefix in multifPrefix[0:4]: # Only first four components have temperatures output
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)];
                    nozzle.responses[k].append(assignTemperature(k,filesuffix,output));
            
        # Assign temperature ratio results
        elif 'TEMP_RATIO' in k:

            assigned = 0;
            for prefix in multifPrefix:
                if prefix in k:
                    # Provide corresponding aeros filename to function to calculate temp ratio
                    # in this component
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)];
                    nozzle.responses[k] = assignTempRatio(k,filesuffix,mat[multifPrefix.index(prefix)],output);
                    assigned = 1;
                    break; # Go to next k in nozzle.responses
            
            # If k is not specific, provide temp ratio for all components
            if assigned == 0:
                filenames = []
                for prefix in multifPrefix[0:4]: # Only first four components have temperatures output
                    filesuffix = aerosSuffix[multifPrefix.index(prefix)];
                    nozzle.responses[k].append(assignTempRatio(k,filesuffix,mat[multifPrefix.index(prefix)],output));
            
    # END for k in nozzle.responses
    
    if output == 'verbose':
        sys.stdout.write('Aero-S responses obtained\n');
    
    return 0;                            

        
