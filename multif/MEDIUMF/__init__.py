from .. import nozzle
from .. import SU2

from run import *
from meshgeneration import *
from runSU2 import *
import SU2postprocessing
import AEROSpostprocessing

try:
	from runAEROS import *
except ImportError:
	pass;
