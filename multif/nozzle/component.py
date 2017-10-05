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
            self.name = 'NONDESCRIPT_AXISYMMETRIC_WALL'
        self.param = 'NONE'
        
class NonaxisymmetricWall:
    def __init__(self,*args):
        if len(args) == 1:
            self.name = args[0]
        else:
            self.name = 'NONDESCRIPT_NONAXISYMMETRIC_WALL'
        self.param = 'NONE'
        
class Wall:
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
        
class Distribution:
    def __init__(self,*args):
        if len(args) == 1:
            self.name = args[0]
        else:
            self.name = 'NONDESCRIPT_DISTRIBUTION'
        self.param = 'NONE'
        
class Spline:
    def __init__(self,*args):
        if len(args) == 1:
            self.name = args[0]
        else:
            self.name = 'NONDESCRIPT_SPLINE'
        self.coefs = []
        self.knots = []
        self.coefs_size = 0