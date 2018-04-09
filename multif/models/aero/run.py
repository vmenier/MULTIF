"""
Run aero analyses (any fidelity) and corresponding post-processing.
"""

import LOWF
import MEDIUMF
import HIGHF


def aeroAnalysis(nozzle, sst_perturbation=None, run_analysis=True, 
    output='verbose'):
    """
    Run aero analysis for provided nozzle.
    """

    if run_analysis:
        if nozzle.method == 'NONIDEALNOZZLE':
            LOWF.runQ1D.run(nozzle, output=output)
        elif nozzle.dim == '2D':
            #MEDIUMF.runSU2.CheckSU2Version(nozzle)
            MEDIUMF.runSU2.runSU2(nozzle, sst_perturbation=sst_perturbation, 
                output=output)
        else: # nozzle.dim == '3D':
            #HIGHF.runSU2.CheckSU2Version(nozzle)
            HIGHF.runSU2.runSU2(nozzle, sst_perturbation=sst_perturbation, 
                output=output)

    return
    

def aeroPostprocessing(nozzle, runDir, output='verbose'):
    """
    Run post-processing for provided nozzle, assuming aero analysis has already
    been performed.
    """

    if nozzle.method == 'NONIDEALNOZZLE':
        pass
    elif nozzle.dim == '2D':
        MEDIUMF.postprocessing.PostProcess(nozzle, runDir, output=output)
    else: # nozzle.dim == '3D':
        HIGHF.postprocessing.PostProcess(nozzle, runDir, output=output)

    return
