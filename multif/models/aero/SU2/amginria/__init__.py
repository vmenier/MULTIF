# SU2/inria_amg/__init__.py

from tools     import *
from interface import *
#from setup import setup_amgio

try:
    import _amgio
except:
    raise;

#try:
#    import _amgio
#except:
#    try:
#        setup_amgio();
#    except:
#        raise
#    try:
#        import _amgio
#    except:
#        raise;
#