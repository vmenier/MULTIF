import sys, os, copy
import numpy as np
import multiprocessing

import LOWF
import MEDIUMF

# Wrapping function for independent nozzle analysis in separate directory
def nozzleAnalysis(homedir, index, nozzle, output='verbose'):
    
    if output == 'verbose':
        sys.stdout.write('Entered separate nozzle analysis for index %i\n' % index);

    # Create and enter new directory
    dirname = os.path.join(homedir,'EVAL_' + str(index));
    if not os.path.exists(dirname):
        os.makedirs(dirname);
    os.chdir(dirname);
    
    if output == 'verbose':
        sys.stdout.write('Directory %s created and entered\n' % dirname);    
    
    # Run model analysis    
    if nozzle.dim == '1D':
        LOWF.Run(nozzle,output=output,writeToFile=1);
    elif nozzle.dim == '2D':
        MEDIUMF.Run(nozzle,output=output,writeToFile=1);
    else: # nozzle.dim == '3D'
        HIGHF.Run(nozzle, output=output, writeToFile=1);
                    
    if output == 'verbose':
        sys.stdout.write('Nozzle analysis completed in directory %s\n' % dirname);    
    
    # Exit directory
    os.chdir(homedir);
    
    # Return nozzle    
    return nozzle;


# Calculate and return forward finite difference gradients of nozzle QOI
# If rerun_center = 1 then the center point used for finite difference will
# be recalculated along with all the points corresponding to f.d. steps
def calcGradientsFD(nozzle,fd_step,rerun_center=0,output='verbose'):
    
    # Check that there are enough finite difference steps if a different step
    # size is provided for each variable
    if isinstance(fd_step,list):
        if len(nozzle.dvList) != len(fd_step):
            sys.stderr.write('  ## ERROR: %i finite difference step sizes have '
              'been provided for %i design variables. Please provide %i numbers'
              ' in the FD_STEP_SIZE list.\n\n' % (len(fd_step), \
              len(nozzle.dvList), len(nozzle.dvList)));
            sys.exit(1);
    
    # Gradients w.r.t. to all variables may not be required:
    # 0-index which design variables should have derivatives taken w.r.t. them
    derivativesDV = [i-1 for i in nozzle.derivativesDV];
    
    # Make local list of design variable 0-indexed indices
    dvList = []
    for i in derivativesDV:
        dvList.append(nozzle.dvList[i])

    # Setup each nozzle instance required for evaluation
    nozzleEval = [];
    for i in range(len(derivativesDV)):
        nozzleEval.append(copy.deepcopy(nozzle));
        if isinstance(fd_step,list):
            nozzleEval[i].dvList[derivativesDV[i]] += fd_step[derivativesDV[i]];
        else:
            nozzleEval[i].dvList[derivativesDV[i]] += fd_step;
        nozzleEval[i].UpdateDV(output='quiet');
        nozzleEval[i].SetupWall(output='quiet');
        nozzleEval[i].output_gradients = 0; # do not request gradients
        for k in nozzleEval[i].gradients:
            nozzleEval[i].gradients[k] = None; # do not request gradients
        for j in range(len(nozzleEval[i].outputCode)):
            nozzleEval[i].outputCode[j] = 1; # get value only
        nozzleEval[i].partitions = 1; # run evaluation on 1 core

    # Set up additional nozzle evaluation for centerpoint if desired
    if( rerun_center ):
        nozzleEval.append(copy.deepcopy(nozzle));
        # Design variables and wall do not change
        nozzleEval[-1].output_gradients = 0; # do not request gradients
        for k in nozzleEval[-1].gradients:
            nozzleEval[-1].gradients[k] = None; # do not request gradients
        for j in range(len(nozzleEval[-1].outputCode)):
            nozzleEval[-1].outputCode[j] = 1; # get value only
        nozzleEval[-1].partitions = 1; # run evaluation on 1 core
                    
    # Run gradient evaluations in serial
    if nozzle.partitions <= 1:   
    
        # For each design variable
        saveResponse = [];
        for i in range(len(nozzleEval)):

            # Create and enter new directory
            dirname = 'EVAL_' + str(i);
            if not os.path.exists(dirname):
                os.makedirs(dirname);
            os.chdir(dirname);
            
            if output == 'verbose':
                sys.stdout.write('Directory %s created and entered\n' % dirname);    
            
            # Run model analysis
            if nozzle.dim == '1D':
                LOWF.Run(nozzleEval[i], output=output, writeToFile=1);
            elif nozzle.dim == '2D':
                MEDIUMF.Run(nozzleEval[i], output=output, writeToFile=1);
            else: # nozzle.dim == '3D'
                HIGHF.Run(nozzleEval[i], output=output, writeToFile=1);

            if output == 'verbose':
                sys.stdout.write('Nozzle analysis completed in directory %s\n' % dirname);    
            
            # Exit directory
            os.chdir('..');  

            # Save reponses for assembly of gradients afterwards
            saveResponse.append(nozzleEval[i].responses)         

        # Calculate gradients here
        for k in nozzle.gradients:

            # Only calculate gradients that are requested, and avoid calculating 
            # mass gradient since it has already been calculated
            if nozzle.gradients[k] is not None and k not in ['MASS'] and \
                k not in ['MASS_WALL_ONLY']:

                if rerun_center == 0:
                    nozzleResponse = nozzle.responses[k];
                else:
                    nozzleResponse = nozzleEval[-1].responses[k];

                for i in range(len(derivativesDV)):
                    if isinstance(fd_step,list):
                        localGrad = (nozzleEval[i].responses[k] - nozzleResponse)/fd_step[derivativesDV[i]];
                    else:
                        localGrad = (nozzleEval[i].responses[k] - nozzleResponse)/fd_step;
                    nozzle.gradients[k].append(localGrad);
    
    # Run gradient evaluations in parallel            
    else:

        # Home directory where all subfolders are stored
        homedir = os.getcwd();
        
        # Start Python's multiprocessing pool
        if output == 'verbose':
            sys.stdout.write('Starting multiprocessing pool with %i '
              'processes\n' % nozzle.partitions);
        pool = multiprocessing.Pool(processes=nozzle.partitions);
        
        # Setup list to hold results for nozzle evaluation
        mEval = [];
        rEval = [];
        for i in range(len(nozzleEval)):
            mEval.append(-1);
            rEval.append(-1);
        
        for i in range(len(nozzleEval)):
            if output == 'verbose':
                sys.stdout.write('Adding analysis %i to the pool\n' % i);
            mEval[i] = pool.apply_async(nozzleAnalysis,(homedir,i,nozzleEval[i],output))
            
        pool.close();
        pool.join();
        
        # Obtain results of calculations
        for i in range(len(mEval)):
            # each entry in rEval contains a nozzle class instance 
            # with analysis responses
            rEval[i] = mEval[i].get();
            
        # Calculate gradients here
        for i in range(len(derivativesDV)):

            for k in nozzle.gradients:

                if rerun_center == 0:
                    nozzleResponse = nozzle.responses[k];
                else:
                    nozzleResponse = rEval[-1].responses[k];

                # Only calculate gradients that are requested, and avoid calculating 
                # mass gradient since it has already been calculated
                if nozzle.gradients[k] is not None and k not in ['MASS'] and  \
                  k not in ['MASS_WALL_ONLY']:
                    if isinstance(fd_step,list):
                        localGrad = (rEval[i].responses[k] - nozzleResponse)/fd_step[derivativesDV[i]];
                    else:
                        localGrad = (rEval[i].responses[k] - nozzleResponse)/fd_step;                    
                    nozzle.gradients[k].append(localGrad);
    	
    return nozzle.gradients
    


