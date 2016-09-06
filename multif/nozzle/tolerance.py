# -*- coding: utf-8 -*-
"""
Mission module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Tolerance:
    def __init__(self):			
				self.exitTempPercentError         = 1e-10    # tolerance used only for quasi-1D model
				self.solverRelTol                 = 1e-10   # general solver relative tolerance
				self.solverAbsTol                 = 1e-10   # general solver absolute tolerance
				self.solverApparentThroatLocation = 1e-6    # tolerance usedy only for quasi-1D model
    
    # Set solver's absolute tolerance
    def setAbsTol(self, tol):
        self.solverAbsTol = tol;
    
    # Set solver's relative tolerance
    def setRelTol(self, tol):
        self.solverRelTol = tol;
