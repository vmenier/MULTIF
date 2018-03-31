import os
import numpy as np

from .. import _nozzle_module
import MEDIUMF

def thermalAnalysis(nozzle, run_analysis=True, output='verbose'):
    """
    Run thermal analysis for provided nozzle.
    """

    # Generate mesh and input files if they are not already generated
    if os.path.exists("wall_results.txt"):
        nozzle.wallResults = np.loadtxt("wall_results.txt") # for lo-fi
    if not os.path.exists('nozzle.aeroh'):
        MEDIUMF.runAEROS.prepareAEROS(nozzle, run_analysis=run_analysis,
            output=output)
    else:
        print("WARNING: Aero-S thermal input files and mesh already exist" + \
              "and will not be regenerated.")
    
    # Run thermal analysis
    if run_analysis:
        MEDIUMF.runAEROS.callThermalAnalysis()

    return


def thermalPostprocessing(nozzle, runDir, output='verbose'):
    """
    Run post-processing for provided nozzle, assuming thermal analysis has 
    already been performed.
    """

    MEDIUMF.postprocessing.PostProcess(nozzle, runDir, output=output)

    return


def structuralAnalysis(nozzle, run_analysis=True, output='verbose'):
    """
    Run structural analysis for provided nozzle.
    """

    # Generate mesh and input files if they are not already generated
    if os.path.exists("wall_results.txt"):
        nozzle.wallResults = np.loadtxt("wall_results.txt") # for lo-fi
    if not os.path.exists('nozzle.aeros'):
        MEDIUMF.runAEROS.prepareAEROS(nozzle, run_analysis=run_analysis,
            output=output)
    else:
        print("WARNING: Aero-S structural input files and mesh already exist" + \
              "and will not be regenerated.")
    
    # Run structural analysis
    if run_analysis:
        MEDIUMF.runAEROS.callStructuralAnalysis()

    return


def structuralPostprocessing(nozzle, runDir, output='verbose'):
    """
    Run post-processing for provided nozzle, assuming structural analysis has 
    already been performed.
    """

    MEDIUMF.postprocessing.PostProcess(nozzle, runDir, output=output)

    return


def massAnalysis(nozzle, run_analysis=True, output='verbose'):
    """
    Run mass analysis for provided nozzle.
    """

    total_mass, wall_mass = MEDIUMF.runAEROS.getMass(nozzle, 
        run_analysis=run_analysis, output=output)
    if 'MASS' in nozzle.qoi.names:
        nozzle.qoi.setValue('MASS', total_mass)
    if 'MASS_WALL_ONLY' in nozzle.qoi.names:
        nozzle.qoi.setValue('MASS_WALL_ONLY', wall_mass)

    return