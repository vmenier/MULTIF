# -*- coding: utf-8 -*-
"""
Mission module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Tolerance:
    def __init__(self):			
				self.exitTempPercentError         = 1e-8
				self.solverRelTol                 = 1e-10
				self.solverAbsTol                 = 1e-10
				self.solverApparentThroatLocation = 1e-6
        
    def setTol(self, tol):
        self.solverRelTol = tol;
        self.solverAbsTol = tol;
