# -*- coding: utf-8 -*-
"""
Material module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Inlet:
    def __init__(self, PstagInlet, TstagInlet):
        self.setPstag(PstagInlet)
        self.setTstag(TstagInlet)
        
    def setMach(self, mach):
        self.mach = mach
        
    def setPstag(self,Pstag):
        self.Pstag = Pstag # Pa, stagnation pressure of inlet
    
    def setTstag(self,Tstag):
        self.Tstag = Tstag # K, stagnation temp. of inlet