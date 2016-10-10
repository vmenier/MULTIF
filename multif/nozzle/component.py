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
        self.param = 'NONE'
    
class Baffles:
    def __init__(self,n):
        self.n = int(n)
        self.location = None
        self.height = None
        self.thickness = None
    
class Stringers:
    def __init__(self,n):
        self.n = int(n)
        self.height = None
        self.thickness = None