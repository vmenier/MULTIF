import aero

try:
    import thermstruct
except:
    print "## Warning: thermstruct module not found."
    
import _meshutils_module
import _mshint_module

try:
    import _nozzle_module
except:
    print "## Warning: _nozzle_module not found."