# -*- coding: utf-8 -*-
"""
Define nozzle component classes

Rick Fenrich 6/28/16
"""

class AxisymmetricWall:
    def __init__(self,*args):
        if len(args) == 1:
            self.name = args[0]
        else:
            self.name = 'NONDESCRIPT_WALL'
    
class Baffles:
    def __init__(self):
        pass
    
class Stringers:
    def __init__(self):
        pass