import sys, os, copy
import numpy as np
import multiprocessing

import LOWF
import MEDIUMF
import HIGHF

def init(lock):

    global LOCK;
    LOCK = lock;

    return

# Wrapping function for independent nozzle analysis in separate directory
def nozzleAnalysis(homedir, index, nozzle, skipAero=0, output='verbose'):
    
    if output == 'verbose':
        sys.stdout.write('Entered separate nozzle analysis for index %i\n' % index);

    # Create and enter new directory
    dirname = os.path.join(homedir,'EVAL_' + str(index));
    if not os.path.exists(dirname):
        os.makedirs(dirname);
    os.chdir(dirname);
    
    if output == 'verbose':
        sys.stdout.write('Directory %s created and entered\n' % dirname);    

    # Write input file corresponding to this analysis for debugging
    if nozzle.inputDVformat == 'DAKOTA':
        sys.stdout.write('Note: PLAIN input format is used for writing of design variable input file\n');
        np.savetxt(nozzle.inputDVfilename,nozzle.dvList,fmt='%0.16f');
    else: # should be 'PLAIN'
        np.savetxt(nozzle.inputDVfilename,nozzle.dvList,fmt='%0.16f');

    # Link necessary files if aero analysis is to be skipped
    if skipAero == 1:
        LOCK.acquire() # Only one process should link these files at a time
        if os.path.exists(os.path.join(homedir,'nozzle.dat')):
            if os.path.exists('nozzle.dat'):
                os.remove('nozzle.dat');
            os.link(os.path.join(homedir,'nozzle.dat'),'nozzle.dat');
        if os.path.exists(os.path.join(homedir,'nozzle.su2')):
            if os.path.exists('nozzle.su2'):
                os.remove('nozzle.su2');
            os.link(os.path.join(homedir,'nozzle.su2'),'nozzle.su2');
        LOCK.release() # Release the lock so others can access these files
    
    # Run model analysis    
    if nozzle.dim == '1D':
        LOWF.Run(nozzle, output=output, writeToFile=1, skipAero=skipAero);
    elif nozzle.dim == '2D':
        MEDIUMF.Run(nozzle, output=output, writeToFile=1, skipAero=skipAero);
    else: # nozzle.dim == '3D'
        HIGHF.Run(nozzle, output=output, writeToFile=1, skipAero=skipAero);
                    
    if output == 'verbose':
        sys.stdout.write('Nozzle analysis completed in directory %s\n' % dirname);    
    
    # Exit directory
    os.chdir(homedir);
    
    # Return nozzle    
    return nozzle;  


# Calculate and return forward finite difference gradients of nozzle QOI
# If rerun_center = 1 then the center point used for finite difference will
# be recalculated along with all the points corresponding to f.d. steps
def calcGradientsFD(nozzle, fd_step, rerun_center=0, output='verbose'):
    
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
    skipAeroList = [];
    for i in range(len(derivativesDV)):

        # Determine whether aero analysis is required for this DV
        iDV_Head = [j for j in range(len(nozzle.DV_Head)) if i >= nozzle.DV_Head[j]][0]
        #tagtmp = nozzle.DV_Tags[iDV_Head]
        if nozzle.DV_Effect[iDV_Head] in [2, 3, 4]: # no need for aero analysis
            skipAeroList.append(1);
        else:
            skipAeroList.append(0);

        # Copy and setup nozzle
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
            if nozzleEval[i].outputCode <= 1:
                nozzleEval[i].outputCode[j] = 0; # no need to get any value
            else:
                nozzleEval[i].outputCode[j] = 1; # get value only
        nozzleEval[i].nTasks = 1; # run 1 task (this is it)
        # nozzleEval[i].cpusPerTask remains the same

    # Set up additional nozzle evaluation for centerpoint if desired
    if( rerun_center and sum(skipAeroList) > 0 ):
        # If some aero analyses are getting skipped and the center point is 
        # being rerun, we must rerun the center point first with the aero 
        # analysis so we can use this aero analysis for the runs which skip it
        sys.stderr.write("WARNING: Rerun center is requested, but center will"
                         " not be rerun since skipping the aero analysis " 
                         "is requested for some design variable finite "
                         "differences.\n");
        rerun_center = 0;
    elif( rerun_center ):
        # Rerun center is no longer applicable since the nTasks and cpusPerTask
        # are implemented to control parallelization of gradients. The 
        # center point for finite difference will use the same number of cores
        # as the points a difference away.
        sys.stderr.write("WARNING: Rerun center is not happening since it is no"
                         " longer applicable.\n")
        rerun_center = 0;
        # skipAeroList.append(0); #rerun aero analysis if QoI are effected by it
        # nozzleEval.append(copy.deepcopy(nozzle));
        # # Design variables and wall do not change
        # nozzleEval[-1].output_gradients = 0; # do not request gradients
        # for k in nozzleEval[-1].gradients:
        #     nozzleEval[-1].gradients[k] = None; # do not request gradients
        # for j in range(len(nozzleEval[-1].outputCode)):
        #     nozzleEval[-1].outputCode[j] = 1; # get value only

    # Home directory where all subfolders are stored
    homedir = os.getcwd();

    # Lock to prevent race condition when input files are linked at same time
    lock = multiprocessing.Lock()
                    
    # Run gradient evaluations in serial
    if nozzle.nTasks <= 1:

        # Initialize global LOCK
        init(lock);
    
        # For each design variable
        saveResponse = [];
        for i in range(len(nozzleEval)):

            # # Create and enter new directory
            # dirname = 'EVAL_' + str(i);
            # if not os.path.exists(dirname):
            #     os.makedirs(dirname);
            # os.chdir(dirname);
            
            # if output == 'verbose':
            #     sys.stdout.write('Directory %s created and entered\n' % dirname);    
            
            # # Write input file corresponding to this analysis for debugging
            # if nozzleEval[i].inputDVformat == 'DAKOTA':
            #     sys.stdout.write('Note: PLAIN input format is used for writing of design variable input file\n');
            #     np.savetxt(nozzleEval[i].inputDVfilename,nozzleEval[i].dvList,fmt='%0.16f');
            # else: # should be 'PLAIN'
            #     np.savetxt(nozzleEval[i].inputDVfilename,nozzleEval[i].dvList,fmt='%0.16f');

            # # Run model analysis
            # if nozzle.dim == '1D':
            #     LOWF.Run(nozzleEval[i], output=output, writeToFile=1);
            # elif nozzle.dim == '2D':
            #     MEDIUMF.Run(nozzleEval[i], output=output, writeToFile=1);
            # else: # nozzle.dim == '3D'
            #     HIGHF.Run(nozzleEval[i], output=output, writeToFile=1);

            # if output == 'verbose':
            #     sys.stdout.write('Nozzle analysis completed in directory %s\n' % dirname);    
            
            # # Exit directory
            # os.chdir('..');
            nozzleAnalysis(homedir, i, nozzleEval[i], skipAero=skipAeroList[i], output='verbose')

            # Save reponses for assembly of gradients afterwards
            saveResponse.append(nozzleEval[i].responses)         

        # Calculate gradients here
        for k in nozzle.gradients:

            # Only calculate gradients that are requested, and avoid calculating 
            # mass gradient since it has already been calculated
            # if nozzle.gradients[k] is not None and k not in ['MASS'] and \
            #     k not in ['MASS_WALL_ONLY']:
            if nozzle.gradients[k] is not None:

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
        
        # Start Python's multiprocessing pool
        if output == 'verbose':
            sys.stdout.write('Starting multiprocessing pool with %i '
              'processes\n' % nozzle.nTasks);
        pool = multiprocessing.Pool(initializer=init, initargs=(lock,), processes=nozzle.nTasks);
        
        # Setup list to hold results for nozzle evaluation
        mEval = [];
        rEval = [];
        for i in range(len(nozzleEval)):
            mEval.append(-1);
            rEval.append(-1);
        
        for i in range(len(nozzleEval)):
            if output == 'verbose':
                sys.stdout.write('Adding analysis %i to the pool\n' % i);
            mEval[i] = pool.apply_async(nozzleAnalysis, (homedir, i, nozzleEval[i], skipAeroList[i], output))
            
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
                # if nozzle.gradients[k] is not None and k not in ['MASS'] and  \
                #   k not in ['MASS_WALL_ONLY']:
                if nozzle.gradients[k] is not None:
                    if isinstance(fd_step,list):
                        localGrad = (rEval[i].responses[k] - nozzleResponse)/fd_step[derivativesDV[i]];
                    else:
                        localGrad = (rEval[i].responses[k] - nozzleResponse)/fd_step;                    
                    nozzle.gradients[k].append(localGrad);
    	
    return nozzle.gradients
    


