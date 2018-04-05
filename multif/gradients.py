import sys, os, copy
import numpy as np
import multiprocessing

import analysis

def init(lock):

    global LOCK
    LOCK = lock

    return

def getValue(qoi, q):
    """
    Wrapper function with warnings for obtaining function value.
    """

    value = qoi.getValue(q)

    if len(value) == 0:
        value = -1.
    elif len(value) == 1:
        value = value[0]
    else:
        print("WARNING: Multiple values for QoI %s will not be output." % q)
        value = value[0]

    return value


# Wrapping function for independent nozzle analysis in separate directory
def systemAnalysis(index, nozzle, post_process_only=False, skip_aero=False, 
    skip_post_processing=False, output='verbose'):

    cwd = os.getcwd()   
    
    if output == 'verbose':
        print('Entered separate nozzle analysis for index %i' % index)

    # Create and enter new directory
    dirname = os.path.join(cwd,'EVAL_' + str(index))
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    os.chdir(dirname)
    
    if output == 'verbose':
        print('Directory %s created and entered' % dirname)    

    # Write input file corresponding to this analysis for debugging
    if nozzle.inputDVformat == 'DAKOTA':
        print('Note: PLAIN input format is used for writing of design variable input file')
        np.savetxt(nozzle.inputDVfilename,nozzle.dvList,fmt='%0.16f')
    else: # should be 'PLAIN'
        np.savetxt(nozzle.inputDVfilename,nozzle.dvList,fmt='%0.16f')

    # Link necessary files if aero analysis is to be skipped
    
    if skip_aero:
        LOCK.acquire() # Only one process should link these files at a time
        if os.path.exists(os.path.join(cwd,'nozzle.dat')):
            if os.path.exists('nozzle.dat'):
                os.remove('nozzle.dat')
            os.link(os.path.join(cwd,'nozzle.dat'),'nozzle.dat')
        if os.path.exists(os.path.join(cwd,'nozzle.su2')):
            if os.path.exists('nozzle.su2'):
                os.remove('nozzle.su2')
            os.link(os.path.join(cwd,'nozzle.su2'),'nozzle.su2')
        LOCK.release() # Release the lock so others can access these files
    
    # Run model analysis
    analysis.run(nozzle, write_to_file=True, 
        post_process_only=post_process_only, skip_aero=skip_aero,
        skip_post_processing=skip_post_processing, output=output)
                    
    if output == 'verbose':
        print('Nozzle analysis completed in directory %s' % dirname)    
    
    # Exit directory
    os.chdir(cwd)
    
    # Return nozzle    
    return nozzle  


def calcGradientsFFD(system, fd_step, post_process_only=False, 
    skip_aero=False, skip_post_processing=False, output='verbose'):
    """
    Calculate and return forward finite difference gradients of system QOI.

    Arguments:
    system: a class instance that contains at least the following attributes:
        dvList
        DV_Head
        DV_Effect
        nTasks: integer specifying number of tasks to be concurrently run
        qoi: a Response() class instance
    fd_step: forward finite difference step size (float or list)
    """
    
    # Check that there are enough finite difference steps if a different step
    # size is provided for each variable
    if isinstance(fd_step,list):
        if len(system.dvList) != len(fd_step):
            raise RuntimeError('%i finite difference step sizes have '
              'been provided for %i design variables. Please provide %i numbers'
              ' in the FD_STEP_SIZE list.\n\n' % (len(fd_step), \
              len(system.dvList), len(system.dvList)))
    
    # Gradients w.r.t. to all variables may not be required:
    # 0-index which design variables should have derivatives taken w.r.t. them
    derivativesDV = [i-1 for i in system.derivativesDV]
    
    # Make local list of design variable 0-indexed indices
    dvList = []
    for i in derivativesDV:
        dvList.append(system.dvList[i])
    
    # Setup each system instance required for evaluation
    systemEval = []
    skipAeroList = []
    for i in range(len(derivativesDV)):

        # Determine whether aero analysis is required for this DV
        iDV_Head = [j for j in range(len(system.DV_Head)) if i >= system.DV_Head[j]][0]
        #tagtmp = system.DV_Tags[iDV_Head]
        if system.DV_Effect[iDV_Head] in [2, 3, 4] or skip_aero: # no need for aero analysis
            skipAeroList.append(1)
        else:
            skipAeroList.append(0)

        # Copy and setup system
        systemEval.append(copy.deepcopy(system))
        if isinstance(fd_step,list):
            systemEval[i].dvList[derivativesDV[i]] += fd_step[derivativesDV[i]]
        else:
            systemEval[i].dvList[derivativesDV[i]] += fd_step
        systemEval[i].UpdateDV(output='quiet')
        systemEval[i].SetupWall(output='quiet')
        systemEval[i].gradientsFlag = False # do not request gradients
        for k in systemEval[i].qoi.names:
            systemEval[i].qoi.initializeValue(k)
            systemEval[i].qoi.initializeGradient(k)
            systemEval[i].qoi.setGradient(k, None) # do not request gradients
        for j in range(len(systemEval[i].outputCode)):
            if systemEval[i].outputCode <= 1:
                systemEval[i].outputCode[j] = 0 # no need to get any value
            else:
                systemEval[i].outputCode[j] = 1 # get value only
        systemEval[i].nTasks = 1 # run 1 task (this is it)
        # systemEval[i].cpusPerTask remains the same

    # Lock to prevent race condition when input files are linked at same time
    lock = multiprocessing.Lock()
                    
    # Run gradient evaluations in serial
    if system.nTasks <= 1:

        # Initialize global LOCK
        init(lock)
    
        # Run analysis for each design variable
        for i in range(len(systemEval)):
            systemAnalysis(i, systemEval[i], post_process_only=post_process_only,
                skip_aero=skipAeroList[i], 
                skip_post_processing=skip_post_processing,
                output='verbose')

        # Calculate gradients here
        for q in system.qoi.names:

            if system.qoi.getGradient(q) is not None:
                
                q0 = getValue(system.qoi, q)


                if system.qoi.getKind(q) == 'field':
                    localGrad = np.zeros((len(q0),len(systemEval)))
                else:
                    localGrad = np.zeros((1,len(systemEval)))

                for i, r in enumerate(systemEval):

                    q1 = getValue(r.qoi, q)

                    if isinstance(fd_step,list):
                        localGrad[:,i] = (q1 - q0)/fd_step[derivativesDV[i]]
                    else:
                        localGrad[:,i] = (q1 - q0)/fd_step
                    
                system.qoi.setGradient(q, np.squeeze(localGrad))
    
    # Run gradient evaluations in parallel using Python multiprocessing.
    # This will only work on shared memory compute node.         
    else:
        
        # Start Python's multiprocessing pool
        if output == 'verbose':
            print('Starting multiprocessing pool with %i processes' % 
                system.nTasks)
        pool = multiprocessing.Pool(initializer=init, initargs=(lock,), 
            processes=system.nTasks)
        
        # Setup list to hold results for system evaluation
        mEval = []
        rEval = []
        for i in range(len(systemEval)):
            mEval.append(-1)
            rEval.append(-1)
        
        for i in range(len(systemEval)):
            if output == 'verbose':
                print('Adding analysis %i to the pool' % i)
            mEval[i] = pool.apply_async(systemAnalysis, (i, systemEval[i], 
                post_process_only, skipAeroList[i], skip_post_processing, output))
            
        pool.close()
        pool.join()
        
        # Obtain results of calculations
        for i in range(len(mEval)):
            rEval[i] = mEval[i].get() # each rEval[i] contains an analysis instance
            
        # Calculate gradients here
        for q in system.qoi.names:

            if system.qoi.getGradient(q) is not None:

                q0 = getValue(system.qoi, q)
                
                if system.qoi.getKind(q) == 'field':
                    localGrad = np.zeros((len(q0),len(rEval)))
                else:
                    localGrad = np.zeros((1,len(rEval)))

                for i, r in enumerate(rEval):

                    q1 = getValue(r.qoi, q)

                    if isinstance(fd_step,list):
                        localGrad[:,i] = (q1 - q0)/fd_step[derivativesDV[i]]
                    else:
                        localGrad[:,i] = (q1 - q0)/fd_step  

                system.qoi.setGradient(q, np.squeeze(localGrad))
    	
    return
    


