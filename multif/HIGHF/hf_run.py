import sys, os

import numpy as np

from .. import SU2
from hf_meshgeneration import *
from hf_runSU2 import *
import hf_postprocessing

from hf_aeros import *

from .. import nozzle as nozzlemod

try:
    from multif.MEDIUMF.runAEROS import *
except ImportError as e:
    print 'Error importing all functions from runAEROS in hf_run.py.'
    print e
    print


def CheckOptions (nozzle):
    
    print "Check options"
    if nozzle.dim != '3D':
        sys.stderr.write("\n  ## ERROR : High-fidelity is 3D.\n\n");
        sys.exit(0);
    
def Run( nozzle, **kwargs ):

    output = 'verbose';
    writeToFile=1;
    postpro = 0;    
    
    if 'output' in kwargs:
        output = kwargs['output'];
        
    if 'writeToFile' in kwargs:
        writeToFile = int(kwargs['writeToFile']);
        
    if 'postpro' in kwargs:
        postpro = int(kwargs['postpro']);    

    # # Obtain mass and volume
    # if 'MASS' in nozzle.responses or 'VOLUME' in nozzle.responses or 'MASS_WALL_ONLY' in nozzle.responses:
    #     volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle)
    #     if 'MASS' in nozzle.responses:
    #         nozzle.responses['MASS'] = np.sum(mass)
    #     if 'VOLUME' in nozzle.responses:
    #         nozzle.responses['VOLUME'] = np.sum(volume)
    #     if 'MASS_WALL_ONLY' in nozzle.responses:
    #         nozzle.responses['MASS_WALL_ONLY'] = np.sum(mass[:len(nozzle.wall.layer)])

    # # Calculate mass gradients if necessary
    # if 'MASS' in nozzle.gradients and nozzle.gradients['MASS'] is not None:
    #     if ( nozzle.gradientsMethod == 'ADJOINT' ):
    #         # Convergence study using B-spline coefs show finite difference mass gradients
    #         # converge. Most accurate gradients use absolute step size 1e-8. RWF 5/10/17
    #         # The above conclusion is based on 2D axisymmetric mass calculations.
    #         nozzle.gradients['MASS'] = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-8);
    #     elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
    #         nozzle.gradients['MASS'] = nozzlemod.geometry.calcMassGradientsFD( \
    #             nozzle,nozzle.fd_step_size);
    #     else:
    #         sys.stderr.write('  ## ERROR : Unknown gradients computation '
    #             'method.\n');
    #         sys.exit(1);
    
    # # Calculate volume gradients if necessary
    # if 'VOLUME' in nozzle.gradients and nozzle.gradients['VOLUME'] is not None:
    #     sys.stderr.write('\n ## ERROR : gradients for VOLUME are not supported\n\n');
    #     sys.exit(1);
        
    # # Obtain gradients of mass of wall if requested        
    # if 'MASS_WALL_ONLY' in nozzle.gradients and nozzle.gradients['MASS_WALL_ONLY'] is not None:
    #     if ( nozzle.gradientsMethod == 'ADJOINT' ):
    #         # Convergence study using B-spline coefs show finite difference mass gradients
    #         # converge. Most accurate gradients use absolute step size 1e-8. RWF 5/10/17
    #         # The above conclusion is based on 2D axisymmetric mass calculations.
    #         nozzle.gradients['MASS_WALL_ONLY'] = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-8,components='wall-only');
    #     elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
    #         nozzle.gradients['MASS_WALL_ONLY'] = nozzlemod.geometry.calcMassGradientsFD(\
    #             nozzle,nozzle.fd_step_size,components='wall-only');
    #     else:
    #         sys.stderr.write('  ## ERROR : Unknown gradients computation '
    #             'method.\n');
    #         sys.exit(1);
    
    # Run aero-thermal-structural analysis if necessary
    runAeroThermalStructuralProblem = 0;
    for k in nozzle.responses:
        if k not in ['MASS','VOLUME','MASS_WALL_ONLY']:
            runAeroThermalStructuralProblem = 1;    
    
    # Run aero-thermal-structural gradient analysis if necessary
    if nozzle.gradientsFlag == 1:
        runAeroThermalStructuralGradients = 0;        
        for k in nozzle.gradients:
            if k not in ['MASS','VOLUME','MASS_WALL_ONLY']:
                if nozzle.gradients[k] is not None:
                    runAeroThermalStructuralGradients = 1;
    
    if runAeroThermalStructuralProblem:
        
        # Run aero analysis (and thrust adjoint if necessary)  
        CheckSU2Version(nozzle);    
        CheckOptions(nozzle);

        if postpro != 1:
            curDir = os.path.dirname(os.path.realpath(__file__));    
            if nozzle.runDir != '':
                os.chdir(nozzle.runDir);
            # XXX Ensure high-fidelity SU2 analysis runs correctly
            gradCalc = HF_runSU2(nozzle);
            
            # Run thermal/structural analyses
            if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
                # XXX Ensure high-fidelity AEROS analysis runs correctly
                
                try:
                    runAEROS(nozzle, output);
                except:
                    sys.stderr.write("\n  ## WARNING : runAEROS disabled.\n\n");
                
        # ------- BEGIN HACK VERIF EQUIV AREA
        # -> Creates N meshes of nozzle cross sections and compares areas with the analytic ones
        
        #import multif
        #
        #coefs_center = nozzle.wall.centerline.coefs;
        #
        #x_in   = coefs_center[0];
        #x_out  = coefs_center[len(coefs_center)/2-1];
        #
        #xtab = np.linspace(x_in,x_out,15);
        #
        #area_tab = [];
        #
        #for i in range(len(xtab)):
        #    
        #    multif.HIGHF.HF_GenerateExitMesh(nozzle, xtab[i]);
        #    area_msh = multif.HIGHF.Get2DMeshArea("nozzle_exit_hin.mesh");
        #    
        #    rad = MF_GetRadius (xtab[i],nozzle);
        #    area_rad = 0.5*rad*rad*math.pi;
        #    
        #    print "x %lf area_msh %le (rad %lf) area_rad %le (rad %lf)" % (xtab[i],area_msh,math.sqrt(2*area_msh/math.pi),area_rad, rad)
        #    
        #    area_tab.append([xtab[i],area_msh ,math.sqrt(2*area_msh/math.pi),area_rad, rad]);
        #    
        #
        #for i in range(len(area_tab)):
        #    print "x=%lf area_msh %le area_rad %le: " % (area_tab[i][0],area_tab[i][1],area_tab[i][2])
        #
        #hdl=open("equivarea.dat", "w");
        #for i in range(len(area_tab)):
        #    hdl.write("%le %le %le %le %le\n" % (area_tab[i][0],area_tab[i][1],area_tab[i][2],area_tab[i][3],area_tab[i][4]))
        #hdl.close();
        #sys.exit(1);
        
        # ------- END HACK VERIF EQUIV AREA
        
        
        #print "Interface AEROS"
        #
        #MshNam_str = "nozzle_struct.mesh"
        #MshNam_cfd = "nozzle.su2"
        #SolNam_cfd = "nozzle.dat"
        #
        #Crd, Tri, Pres, Temp = hf_FluidStructureInterpolation(MshNam_str, MshNam_cfd, SolNam_cfd);
        #        
        #for i in range(0,10):
        #    print "ver %d (%lf %lf %lf) : pres %lf temp %lf" % (i+1,Crd[i][0],Crd[i][1],Crd[i][2],Pres[i], Temp[i])
        #
        #sys.exit(1);
        #
        
        hf_postprocessing.PostProcess(nozzle, output);
        
        # Assign thermal/structural QoI if required
        if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
            # XXX Ensure AEROS post-processing returns correct outputs
            AEROSpostprocessing.PostProcess(nozzle, output);           
        
        # Calculate gradients if necessary
        if nozzle.gradientsFlag == 1 and runAeroThermalStructuralGradients:
    
            if ( nozzle.gradientsMethod == 'ADJOINT' ):
                
                if gradCalc == 0: # i.e. failed adjoint calculation, use finite differences
                    # XXX Ensure FD gradients work for 3-D nozzle. If this
                    # function needs to be changed, ensure it is also changed
                    # in the other 2 spots below too.
                    multif.gradients.calcGradientsFD(nozzle,nozzle.fd_step_size,output);
                else:
                    # Check for other required gradients
                    otherRequiredGradients = 0;
                    for k in nozzle.gradients:
                        if k not in ['MASS','VOLUME','MASS_WALL_ONLY','THRUST']:
                            otherRequiredGradients = 1;
                            sys.stderr.write(' ## WARNING: QoI gradients desired using ADJOINT '
                              'method which do not have an associated adjoint calculation.\n'
                              ' Namely: %s. The current implementation requires finite '
                              'differencing across the aero analysis, so this method is '
                              'equivalent in computational cost to choosing the FINITE_DIFF'
                              ' method\n' % k);
                    # Do finite difference for other QoI if there are any, but use
                    # adjoint gradients for thrust
                    if otherRequiredGradients:
                        saveThrustGradients = nozzle.gradients['THRUST'];
                        nozzle.gradients['THRUST'] = None;
                        multif.gradients.calcGradientsFD(nozzle,nozzle.fd_step_size,output);
                        nozzle.gradients['THRUST'] = saveThrustGradients;
                        
            elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
            
                multif.gradients.calcGradientsFD(nozzle,nozzle.fd_step_size,output);  
                         
            else:
                sys.stderr.write('  ## ERROR : Unknown gradients computation '
                  'method.\n');
                sys.exit(1);
                
         # Write separate gradients file
            gradFile = open(nozzle.gradientsFile,'w');
            for k in nozzle.outputTags:
                np.savetxt(gradFile,nozzle.gradients[k]);
            gradFile.close();               
    
    # Write data
    if writeToFile:
        if nozzle.outputFormat == 'PLAIN':
            nozzle.WriteOutputFunctions_Plain();
        else:
            nozzle.WriteOutputFunctions_Dakota();
    
    return 0;
    
    
