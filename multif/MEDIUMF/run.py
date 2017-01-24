from .. import SU2
from meshgeneration import *
from runSU2 import *
from postprocessing import *

try:
    from runAEROS import *
except ImportError:
    print 'Error importing all functions from runAEROS.\n'

#try:
#    from runAEROS import *
#else ImportError:
#    pass;

def CheckOptions (nozzle):
    
		print "Check options"
    #if nozzle.Dim == 3 :
    #    sys.stderr.write("\n  ## ERROR : Only 2D axisymmetric simulations are available for now.\n\n");
    #    sys.exit(0);
    
    #if nozzle.method == 'RANS':
    #    sys.stderr.write("\n  ## ERROR : Only Euler simulations are available for now.\n\n");
    #    sys.exit(0);
    

def Run( nozzle, output = 'verbose' ):
    
    # --- Check SU2 version
    
    CheckSU2Version(nozzle);
        
    # --- Check fidelity level
    
    CheckOptions (nozzle);
    
    curDir = os.path.dirname(os.path.realpath(__file__));
    
    if nozzle.runDir != '':
        os.chdir(nozzle.runDir);
    
    # --- Run CFD
    
    runSU2 (nozzle);
    
    # --- Run AEROS  
    
    nozzle.runAEROS = 0;
    if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
        nozzle.runAEROS = 1;
    
#        try:
#            #from runAEROS import *
#            if output == 'verbose':      
#                print "SUCCESS IMPORTING AEROS"
#        except ImportError:
#                nozzle.runAEROS = 0;
#                pass;
    
    if output == 'verbose':
        print "RUNAEROS = %d" % nozzle.runAEROS;
    
    if nozzle.runAEROS == 1:
        runAEROS (nozzle);
    else :
        sys.stdout.write('  -- Info: Skip call to AEROS.\n');

    # --- Postprocessing
    
    PostProcessing(nozzle);
    
    if output == 'verbose':
        sys.stdout.write("\n  -- Info : Result directory :  %s\n\n" % nozzle.runDir);
    
    # --- Output results

    if nozzle.outputFormat == 'PLAIN':
        nozzle.WriteOutputFunctions_Plain();
    else:
        nozzle.WriteOutputFunctions_Dakota();
    

