import os, sys
import copy
import numpy as np
import time

import gradients
import utils
import models


def runAeroAnalysis(nozzle, analysis, post_process_only=False, 
    skip_post_processing=False, run_analysis=True, output='verbose'):
    """
    Run aero analysis.
    """

    # Check if SST perturbations need to be performed 
    sst_perturb = None
    for s in nozzle.cfd.sst_tags:
        if s in analysis:
            sst_perturb = s
            break

    cwd = os.getcwd()

    # If separate directory is explicitly named, use it instead
    if nozzle.runDir:
        # XXX: Remove this option in the future
        runDir = os.path.join(cwd,nozzle.runDir)
        raise RuntimeError("nozzle.runDir is being deprecated")
    # If SST perturbation is being calculated, create separate directory
    elif sst_perturb is not None:
        runDir = os.path.join(cwd,'_'.join(("AERO",sst_perturb)))
    # Create separate directory for normal analyses too
    else:
        runDir = os.path.join(cwd,"AERO")

    # Run analysis
    if not post_process_only:

        # Run analysis
        utils.makeDirAndLink(runDir, [])
        os.chdir(runDir)
        models.aero.run.aeroAnalysis(nozzle, sst_perturbation=sst_perturb,
            run_analysis=run_analysis, output=output)
        os.chdir(cwd)

    # Postprocess
    if not skip_post_processing:
        
        models.aero.run.aeroPostprocessing(nozzle, runDir, output=output)

    return


def runThermalAnalysis(nozzle, analysis, post_process_only=False,
    skip_post_processing=False, run_analysis=True, output='verbose'):
    """
    Run thermal analysis.
    """

    # Check if SST perturbations need to be performed
    sst_perturb = None
    for s in nozzle.cfd.sst_tags:
        if s in analysis:
            sst_perturb = s
            break

    cwd = os.getcwd()

    # Determine directory for analysis
    if sst_perturb is not None:
        tsLocalDir = '_'.join(("THERMOSTRUCTURAL",sst_perturb))
        aeroLocalDir = '_'.join(("AERO",sst_perturb))
    else:
        tsLocalDir = "THERMOSTRUCTURAL"
        aeroLocalDir = "AERO"
    runDir = os.path.join(cwd,tsLocalDir)

    # Run analysis
    if not post_process_only:

        if nozzle.method != 'NONIDEALNOZZLE':
            linkFiles = [os.path.join(aeroLocalDir,"nozzle.su2"),
                        os.path.join(aeroLocalDir,"nozzle.dat")]
        else:
            linkFiles = [os.path.join(aeroLocalDir,"wall_results.txt")]
        utils.makeDirAndLink(runDir, linkFiles)
        os.chdir(runDir)
        models.thermstruct.run.thermalAnalysis(nozzle, 
            run_analysis=run_analysis, output=output)
        os.chdir(cwd)

    # Postprocess
    if not skip_post_processing:

        models.thermstruct.run.thermalPostprocessing(nozzle, runDir, output=output)

    return


def runStructuralAnalysis(nozzle, analysis, post_process_only=False,
    skip_post_processing=False, run_analysis=True, output='verbose'):
    """
    Run structural analysis.
    """

    # Check if SST perturbations need to be performed
    sst_perturb = None
    for s in nozzle.cfd.sst_tags:
        if s in analysis:
            sst_perturb = s
            break

    cwd = os.getcwd()

    # Determine directory for analysis
    if sst_perturb is not None:
        tsLocalDir = '_'.join(("THERMOSTRUCTURAL",sst_perturb))
        aeroLocalDir = '_'.join(("AERO",sst_perturb))
    else:
        tsLocalDir = "THERMOSTRUCTURAL"
        aeroLocalDir = "AERO"
    runDir = os.path.join(cwd,tsLocalDir)

    # Run analysis
    if not post_process_only:
        
        if nozzle.method != 'NONIDEALNOZZLE':
            linkFiles = [os.path.join(aeroLocalDir,"nozzle.su2"),
                        os.path.join(aeroLocalDir,"nozzle.dat")]
        else:
            linkFiles = [os.path.join(aeroLocalDir,"wall_results.txt")]
        utils.makeDirAndLink(runDir, linkFiles)
        os.chdir(runDir)
        models.thermstruct.run.structuralAnalysis(nozzle, 
            run_analysis=run_analysis, output=output)
        os.chdir(cwd)

    # Postprocess
    if not skip_post_processing:

        models.thermstruct.run.structuralPostprocessing(nozzle, runDir, output=output)

    return


def runMassAnalysis(nozzle, analysis, post_process_only=False, 
    skip_post_processing=False, run_analysis=True, output='verbose'):
    """
    Run mass analysis.
    """

    cwd = os.getcwd()

    runDir = os.path.join(cwd,"THERMOSTRUCTURAL")
    # XXX: Potential issue if only MASS is asked for? Do we need aero files?
    utils.makeDirAndLink(runDir, [])
    os.chdir(runDir)
    if post_process_only:
        models.thermstruct.run.massAnalysis(nozzle, run_analysis=False, 
            output=output)
    else:
        models.thermstruct.run.massAnalysis(nozzle, output=output)
    os.chdir(cwd)

    return


def runFullAnalysis(nozzle, output='verbose', post_process_only=False, 
    skip_aero=False, skip_post_processing=False, run_analysis=True):
    """
    Run all analyses required for a given nozzle instance and desired QoI.

    Keyword arguments:
    post_process_only: bool, only run postprocessing
    skip_aero: bool, skip aero setup, analysis, postprocessing entirely
    skip_post_processing: bool, skip postprocessing for each analysis
    run_analysis: bool, run analysis or not; if not, analysis setup and
        postprocessing will still be performed (currently does not affect
        aero analysis)
    """

    # List of analyses which are executed in that order
    analysis_base_list = ['AERO', 'THERMAL', 'STRUCTURAL', 'MASS']
    analysis_list = analysis_base_list
    for s in nozzle.cfd.sst_tags:
        analysis_list = analysis_list + ['_'.join((e,s)) for e in analysis_base_list]
    analysis_values = {} # Store outputs from each analysis
    analysis_gradients = {} # Store outputs from each analysis

    # Determine which analyses are requested
    requested = {key: 0 for key in analysis_list}
    for q in nozzle.qoi.names:
        dep = nozzle.qoi.getDependencies(q)
        for a in analysis_list:
            if a in dep:
                requested[a] += 1

    # Perform each requested analysis
    telapsed = {'AERO': 0., 'THERMAL': 0., 'STRUCTURAL': 0., 'MASS': 0.}
    for analysis in analysis_list:

        # At least 1 QoI uses this analysis
        if requested[analysis] > 0:
            
            # Deep copy the nozzle and remove all responses not generated by
            # this analysis
            noz = copy.deepcopy(nozzle)
            for q in nozzle.qoi.names:
                if analysis not in nozzle.qoi.getAnalysis(q):
                    noz.qoi.remove(q)

            # print analysis
            # print noz.qoi.names
            # print

            if output == 'verbose':
                print('Running %s analysis.' % analysis)

            # Run analysis
            if 'AERO' in analysis:
                if not skip_aero and noz.aeroFlag:
                    t0 = time.time()
                    runAeroAnalysis(noz, analysis, output=output, 
                        skip_post_processing=skip_post_processing,
                        post_process_only=post_process_only,
                        run_analysis=run_analysis)
                    t1 = time.time()
                    telapsed['AERO'] += t1-t0
                else:
                    if output == 'verbose':
                        print('Skipping %s analysis.' % analysis)
                    next
            elif 'THERMAL' in analysis:
                t0 = time.time()
                runThermalAnalysis(noz, analysis, output=output,
                    post_process_only=post_process_only,
                    run_analysis=run_analysis)
                t1 = time.time()
                telapsed['THERMAL'] += t1-t0
            elif 'STRUCTURAL' in analysis:
                t0 = time.time()
                runStructuralAnalysis(noz, analysis, output=output,
                    post_process_only=post_process_only,
                    run_analysis=run_analysis)
                t1 = time.time()
                telapsed['STRUCTURAL'] += t1-t0
            elif 'MASS' in analysis:
                t0 = time.time()
                runMassAnalysis(noz, analysis, 
                    post_process_only=post_process_only,
                    run_analysis=run_analysis)
                t1 = time.time()
                telapsed['MASS'] += t1-t0
            else:
                raise NotImplementedError("Analysis %s not implemented." 
                    % analysis)

            # Save analysis results
            analysis_values[analysis] = noz.qoi
            analysis_gradients[analysis] = noz.qoi

    # Assign responses
    for k, v in analysis_values.iteritems():
        for q in v.names:
            nozzle.qoi.setValue(q, v.getValue(q))

    if nozzle.gradientsFlag:
        for k, v in analysis_gradients.iteritems():
            for q in v.names:
                grad = v.getGradient(q)
                if grad is not None and len(grad) > 0:
                    nozzle.qoi.setGradient(q, v.getGradient(q))

    # Save timing information
    for k in nozzle.elapsedTime:
        nozzle.elapsedTime[k] += telapsed[k]
        
    return


def run(nozzle, output='verbose', write_to_file=True, post_process_only=False, 
    skip_aero=False, skip_post_processing=False):
    """
    Dispatch analyses to determine QoI values and/or gradients w.r.t. design
    variables.
    """

    # Timing information
    nozzle.elapsedTime = {'AERO': 0., 'THERMAL': 0., 'STRUCTURAL': 0., 'MASS': 0.}

    # Run analyses to obtain function values
    # If there are special methods for obtaining gradients (i.e. adjoints),
    # these should be detected and run during the following runFullAnalysis() call.
    if output == 'verbose':
        print("Beginning full nozzle analysis.")
    runFullAnalysis(nozzle, output=output, 
        post_process_only=post_process_only, 
        skip_aero=skip_aero, 
        skip_post_processing=skip_post_processing)

    # Run analyses to obtain function gradients
    # Calculate gradients if necessary
    if nozzle.gradientsFlag:

        # Deep copy the nozzle and remove all responses that already have had
        # gradients calculated
        noz = copy.deepcopy(nozzle)
        for q in nozzle.qoi.names:
            grad = nozzle.qoi.getGradient(q)
            if grad is None or len(grad) > 0:
                noz.qoi.remove(q)

        # All gradients not yet calculated will be approximated by forward
        # finite difference gradients
        if output == 'verbose':
            print("Beginning finite difference gradient calculations for nozzle analysis.")
        gradients.calcGradientsFFD(noz, noz.fd_step_size, output=output,
            post_process_only=post_process_only, 
            skip_aero=skip_aero, 
            skip_post_processing=skip_post_processing)

        # Assign finite difference gradients
        for q in noz.qoi.names:
            nozzle.qoi.setGradient(q, noz.qoi.getGradient(q))

    # Write data
    if write_to_file:
        if output == 'verbose':
            print("Writing output function values to file")
        nozzle.writeOutputFunctions(format=nozzle.outputFormat)  

    # Write timing information
    with open('timing.dat',mode='w') as f:
        for k in nozzle.elapsedTime: 
            f.write('%s, %0.16f\n' % (k,nozzle.elapsedTime[k]))
            if output == 'verbose':
                print("%s analysis took %0.2f seconds" % (k,nozzle.elapsedTime[k]))

    return
