# -*- coding: utf-8 -*-
"""
Geometry module for axisymmetric nozzle. Currently only a B-spline geometry is 
implemented.

Rick Fenrich 9/4/17
"""

import copy
import numpy as np 
import scipy.optimize
import scipy.integrate   

#import matplotlib.pyplot as plt

from .. import _meshutils_module

class Bspline():
    def __init__(self, coefs): # assumes 3rd degree B-spline
        self.type = "B-spline"
        
        # coefs given as a list
        if isinstance(coefs,list): # convert to Numpy array
            if isinstance(coefs[0],list): # nested lists
                self.coefs = np.array(coefs)
            else:
                self.coefs = np.array([coefs[0:len(coefs)/2],coefs[len(coefs)/2:]])
                
                
        # coefs given as a Numpy array
        elif isinstance(coefs,np.ndarray):
            if len(coefs.shape) == 1: # 1D array
                self.coefs = np.array([coefs[0:coefs.size/2],coefs[coefs.size/2:]])
            else:
                self.coefs = coefs
        
        self.knots = np.hstack(([np.zeros(4), np.arange(1.,self.coefs.size/2-3),  \
          np.ones(4)*(self.coefs.size/2-3)])) # calculate here
        self.degree = self.knots.size - self.coefs.size/2 - 1
        self.xstart = self.coefs[0,0]
        self.xend = self.coefs[0,-1]
        self.length = self.coefs[0,-1] - self.coefs[0,0]
        self.inletRadius = self.coefs[1,0]
        self.n = self.coefs.size/2
        
    def findMinimumRadius(self):
        xSeg = np.zeros(self.knots.size)
        ySeg = np.zeros(self.knots.size)
        for ii in range(0,self.knots.size):
            (xTemp,yTemp,temp1,temp2) = uMap3(self.knots[ii],self)
            xSeg[ii] = xTemp
            ySeg[ii] = yTemp
        yMinKnotIndex = np.argmin(ySeg)
        minFunc = lambda x: self.diameter(x).item()
        DMinOld = 1e12
        for ii in range(max(yMinKnotIndex-2,0),min(yMinKnotIndex+3,          \
          self.knots.size-1)):
            lowerBound = xSeg[ii].item()
            upperBound = xSeg[ii+1].item()
            xMin = scipy.optimize.fminbound(minFunc,lowerBound,upperBound)  
            DMin = minFunc(xMin)
            if( DMin < DMinOld ):
                DMinOld = DMin
                xMinOld = xMin
        self.xThroat = xMinOld
        self.yThroat = DMinOld/2
        self.Ainlet2Athroat = (self.inletRadius)**2/self.yThroat**2
        self.Aexit2Athroat = (self.coefs[1,-1])**2/self.yThroat**2
        return (self.xThroat, self.yThroat)
        
    def radius(self, x): # r
        #y = bSplineGeometry(x,self)[0] # Python version (slower)
        y = bSplineGeometryC(x,self)[0] # uses dynamic C library
        return y     
        
    def diameter(self, x): # D
        #y = bSplineGeometry(x,self)[0] # Python version (slower)
        y = bSplineGeometryC(x,self)[0] # uses dynamic C library
        return y*2
        
    def area(self, x): # A
        #y = bSplineGeometry(x,self)[0] # Python version (slower)
        y = bSplineGeometryC(x,self)[0] # uses dynamic C library
        return np.pi*y**2
        
    def radiusGradient(self, x): # drdx
        (y, dydx) = bSplineGeometryC(x,self) # uses dynamic C library
        return dydx
        
    def areaGradient(self, x): # dAdx
        #(y, dydx) = bSplineGeometry(x,self) # Python version (slower)
        (y, dydx) = bSplineGeometryC(x,self) # uses dynamic C library
        return 2*np.pi*y*dydx
        
class PiecewiseLinear:
    def __init__(self,nodes):
        # nodes should be a Numpy array of n x 2, each row contains an
        # x-coordinate and thickness value for a node
        self.type = "piecewise-linear"
        self.nodes = nodes
        self.xstart = nodes[0,0]
        self.xend = nodes[-1,0]
        self.length = np.max(nodes[:,0])
        self.inletRadius = nodes[0,1]
        self.nx = nodes.size/2
        
    def findMinimumRadius(self):
        ii = np.argmin(self.nodes[:,1])
        self.xThroat = self.nodes[ii,0]
        self.rThroat = self.nodes[ii,1]
        self.Ainlet2Athroat = (self.inletRadius)**2/self.rThroat**2
        self.Aexit2Athroat = (self.nodes[-1,1])**2/self.rThroat**2
        return (self.xThroat, self.rThroat)
        
    def radius(self, x): # r
        #r = np.interp(x,self.nodes[:,0],self.nodes[:,1])
        r = piecewiseLinearGeometryC(x,self)[0] # uses dynamic C library
        return r
        
    def diameter(self, x): # D
        #r = self.radius(x)
        r = piecewiseLinearGeometryC(x,self)[0] # uses dynamic C library
        return 2*r
        
    def area(self, x): # A
        #r = self.radius(x)
        r = piecewiseLinearGeometryC(x,self)[0] # uses dynamic C library
        return np.pi*r**2
        
    def radiusGradient(self, x): # drdx
#        if( isinstance(x,float) ):
#            upperIndex = find(x,self.nodes[:,0])
#            if( upperIndex == self.nodes.size/2 ):
#                upperIndex = upperIndex - 1
#            lowerIndex = upperIndex - 1
#            drdx = (self.nodes[upperIndex,1] - self.nodes[lowerIndex,1])/    \
#              (self.nodes[upperIndex,0] - self.nodes[lowerIndex,0])
#        else: # x is an array
#            drdx = np.zeros(x.size)
#            for ii in range(0,x.size):
#                upperIndex = find(x[ii],self.nodes[:,0])
#                if( upperIndex == self.nodes.size/2 ):
#                    upperIndex = upperIndex - 1
#                lowerIndex = upperIndex - 1
#                drdx[ii] = (self.nodes[upperIndex,1] - 
#                  self.nodes[lowerIndex,1])/(self.nodes[upperIndex,0]        \
#                  - self.nodes[lowerIndex,0])
        (r, drdx) = piecewiseLinearGeometryC(x,self) # uses dynamic C library
        return drdx
        
    def areaGradient(self, x): # dAdx
        #r = self.radius(x)
        #drdx = self.radiusGradient(x)  #        
        (r, drdx) = piecewiseLinearGeometryC(x,self) # uses dynamic C library          
        return 2*np.pi*r*drdx

        
class PiecewiseBilinear:
    def __init__(self,nx,ny,nodes):
	    # nodes should be a Numpy array of nx*ny x 3, each row contains an
	    # x-coordinate, y-coordinate, and thickness value for a node, nodes
	    # should be arranged in a rectilinear grid
	    # nodes = np.array([[x1, y1, t1],
	    #                   [x2, y1, t2],
	    #                       ...
	    #                   [x1, y2, t_],
	    #                       etc.
	    self.type = "piecewise-bilinear"
	    self.nodes = nodes # nx*ny x 3 Numpy array of nodes & thicknesses
	    self.size = nx*ny # number of nodes
	    self.nx = nx # dimension of grid in x-direction
	    self.ny = ny # dimension of grid in y-direction
	    
	    # Build interpolating function
	    #self.finterp = interpolate.RectBivariateSpline(nodes[0:nx,0], \
	    #               np.array(list(nodes[:,1])[::nx]), \
	    #               np.reshape(nodes[:,2],(nx,ny),'F'),kx=1,ky=1,s=0)
	    
    def findNearestPoints(self,x,y):
	    # x and y should be scalar
	
	    # Check that point is valid. Extrapolation will not be performed.
	    if( x > np.max(self.nodes[:,0])):
	        raise ValueError('Requested interpolant (%f) > data range.' % x)
	    if( y > np.max(self.nodes[:,1])):
	        raise ValueError('Requested interpolant (%f) > data range.' % y)
	    if( x < np.min(self.nodes[:,0])):
	        raise ValueError('Requested interpolant (%f) < data range.' % x)
	    if( y < np.min(self.nodes[:,1])):
	        raise ValueError('Requested interpolant (%f) < data range.' % y)                    
	    
	    # Cycle through all nodes and find the four that bound point p:
	    # y ^        
	    #   |        
	    #  z12 -- z22
	    #   |  .p  |
	    #  z11 -- z21 ---> x
	    for i in range(self.size):
	        # Find grid point which has x and y coordinates > than that desired
	        if( self.nodes[i,0] >= x and self.nodes[i,1] >= y ):
	            z22 = self.nodes[i,2]
	            x2 = self.nodes[i,0]
	            y2 = self.nodes[i,1]
	            
	            # If we are on left edge of domain
	            if( self.nodes[i,0] == np.min(self.nodes[i,0]) ):
	                # Collapse dimension in x-direction
	                x1 = self.nodes[0,0]
	                z12 = self.nodes[i,2]
	            else:
	                x1 = self.nodes[i-1,0]
	                z12 = self.nodes[i-1,2]                        
	            
	            # If we are on bottom edge of domain
	            if( self.nodes[i,1] == np.min(self.nodes[i,1]) ):
	                # Collapse dimension in y-direction
	                y1 = self.nodes[0,1]
	                z21 = self.nodes[i,2]
	            else:
	                y1 = self.nodes[i-self.nx,1]
	                z21 = self.nodes[i-self.nx,2]
	            
	            # If we are in bottom left corner of domain
	            if( self.nodes[i,0] == np.min(self.nodes[i,0]) and \
	                self.nodes[i,1] == np.min(self.nodes[i,1])):
	                z11 = self.nodes[0,2]
	            else:
	                z11 = self.nodes[i-self.nx-1,2]
	                
	    return x1, x2, y1, y2, z11, z12, z21, z22        
	    
    def height(self,x,y):
	    #z = self.finterp(x,y,grid=False)
	
	    # Perform bilinear interpolation with custom script here
        if( isinstance(x,float) and isinstance(y,float) ):
	        
            x1, x2, y1, y2, z11, z12, z21, z22 = self.findNearestPoints(x,y)
            z = 1./((x2-x1)*(y2-y1))*(z11*(x2-x)*(y2-y) + z21*(x-x1)*(y2-y) + \
	            z12*(x2-x)*(y-y1) + z22*(x-x1)*(y-y1))
            return z
	            
        else: # assume array
            if isinstance(y,list):
                z = np.zeros(len(y))
            elif isinstance(y,np.ndarray):
                z = np.zeros(y.size)

            if isinstance(x,float):
                for i in range(0,z.size):
                    x1, x2, y1, y2, z11, z12, z21, z22 = self.findNearestPoints(x,y[i])
                    z[i] = 1./((x2-x1)*(y2-y1))*(z11*(x2-x)*(y2-y[i]) + \
                        z21*(x-x1)*(y2-y[i]) + \
                        z12*(x2-x)*(y[i]-y1) + z22*(x-x1)*(y[i]-y1)) 
            else:
                for i in range(0,z.size):
                    x1, x2, y1, y2, z11, z12, z21, z22 = self.findNearestPoints(x[i],y[i])
                    z[i] = 1./((x2-x1)*(y2-y1))*(z11*(x2-x[i])*(y2-y[i]) + \
                        z21*(x[i]-x1)*(y2-y[i]) + \
                        z12*(x2-x[i])*(y[i]-y1) + z22*(x[i]-x1)*(y[i]-y1))                

            return z
	    
	def gradient(self,x,y):
	    #temp = self.finterp(x,y,dx=1,dy=1,grid=False)
	
	    # Perform bilinear interpolation with custom script here
	    if( isinstance(x,float) and isinstance(y,float) ):
	        
	        x1, x2, y1, y2, z11, z12, z21, z22 = self.findNearestPoints(x,y)
	        dzdx = 1./((x2-x1)*(y2-y1))*(-z11*(y2-y) + z21*(y2-y) - z12*(y-y1) + \
	               z22*(y-y1))
	        dzdy = 1./((x2-x1)*(y2-y1))*(-z11*(x2-x) - z21*(x-x1) + z12*(x2-x) + \
	               z22*(x-x1))
	            
	    else: # assume array
	        if isinstance(x,list):
	            dzdx = np.zeros(len(x))
	            dzdy = np.zeros(len(x))
	        else:
	            dzdx = np.zeros(x.size)
	            dzdy = np.zeros(x.size)
	        for i in range(0,dzdx.size):
	            x1, x2, y1, y2, z11, z12, z21, z22 = self.findNearestPoints(x[i],y[i])
	            dzdx[i] = 1./((x2-x1)*(y2-y1))*(-z11*(y2-y[i]) + z21*(y2-y[i]) - \
	                      z12*(y[i]-y1) + z22*(y[i]-y1))
	            dzdy[i] = 1./((x2-x1)*(y2-y1))*(-z11*(x2-x[i]) - z21*(x[i]-x1) + \
	                      z12*(x2-x[i]) + z22*(x[i]-x1))
	            if dzdx[i] < pow(10,-16):
	                dzdx[i] = 0.
	            if dzdy[i] < pow(10,-16):
	                dzdy[i] = 0.
	                   
	    return dzdx, dzdy
        

class EllipticalExterior:

    # Elliptical exterior shape is defined by default for a nozzle with an 
    # elliptical exit shape: major axis = 0.92, minor axis = 0.24, centered at
    # (2.33702, 0.19) in the x-z plane. The bottom of the ellipse is truncated
    # at the exit by the line z = 0.122638. This is the bottom of the shovel.
    # For other nozzle exit shapes (such as when a 2d axisymmetric 
    # parameterization is used), the space between the exterior and top and 
    # bottom of the interior nozzle wall exit shape is kept constant at 0.1. 
    # Thus, the exterior is adjustable in such cases.

    def __init__(self,surface,xexit,zoutlettop=0.43,zoutletbottom=0.122638):

        self.xexit = xexit;
        self.surface = surface;

        if( surface == 'top' ):
            # Original parameterization
            # self.angle = 5.; # degrees
            # self.offset = 0.03; # m
            # self.a = 1.5; # m, major axis
            # self.b = 0.2; # m, minor axis

            # Parameterization for 44cm inlet, fixed inlet
            self.angle = 6.; # degrees
            self.a = 3.0; # m, major axis
            self.b = 0.4; # m, minor axis
            self.spacer = 0.2; # m, minimum space between interior surface of
                               # outlet and exterior in y-z plane. Maximum wall 
                               # thickness is 0.0879 m.
            self.offset = zoutlettop - self.b + self.spacer; # m

        elif( surface == 'bottom' ):
            # Original parameterization
            # self.angle = -7.; # degrees
            # self.offset = -0.06; # m
            # self.a = 1.; # m, major axis
            # self.b = 0.05; # m, minor axis

            # Parameterization for 44cm inlet, fixed inlet
            self.angle = -12.; # degrees
            self.a = 4.0; # m, major axis
            self.b = 0.12; # m, minor axis
            self.spacer = 0.2; # m, minimum space between interior surface of
                               # outlet and exterior in y-z plane. Maximum wall 
                               # thickness is 0.0879 m.
            self.offset = zoutletbottom + self.b - self.spacer; # m

        else:
            raise NotImplementedError('Only top or bottom can be used for ' + \
                'surfaces of elliptical exterior.');

    # Given x-coordinate and theta angle in radians measured from Y axis of 
    # ellipse, return global Y and Z coordinates of point on elliptical exterior.
    def coord(self,x,theta):
        
        r = self.a*self.b/np.sqrt(self.b**2*np.cos(theta)**2 + \
            self.a**2*np.sin(theta)**2);
        c = self.offset + (self.xexit - x)*np.tan(np.pi*self.angle/180.);
        y = r*np.cos(theta);
        z = r*np.sin(theta) + c;

        return y, z;

    # Given global x-coordinate and y-coordinate, return global z-coordinate of 
    # point on elliptical exterior.
    def z(self,x,y):

        if np.abs(y) > self.coord(x,0.)[0]:
            raise RuntimeError("Prescribed Y-coordinate lies outside exterior surface definition")

        yabs = np.abs(y)
        f = lambda theta: self.coord(x,theta)[0] - yabs;

        if self.surface == 'top':
            a = 0.;
            b = np.pi/2;
        elif self.surface == 'bottom':
            a = 0.;
            b = -np.pi/2;

        zLocal = scipy.optimize.brentq(f,a,b);

        z = zLocal + self.offset + (self.xexit - x)*np.tan(np.pi*self.angle/180.);

        return z;


#==============================================================================
# Find first 1-based index where scalar xFind < xVec[ii]
#==============================================================================
def find(xFind,xVec):
    counter = 1
    for ii in range(1,xVec.size):
        if(xFind < xVec[ii]):
            break # assume xVec is in ascending order
        counter += 1
    return counter

#==============================================================================
# Calculate x, y, dxdu, and dydu given u for 3rd degree B-spline
# Currently implemented only for scalar inputs of u
#==============================================================================
def uMap3(u,bSpline):
    x = 0.
    y = 0.
    dxdu = 0.
    dydu = 0.
    
    # Calculate the highest knot index that gives a value of u below the given
    # value (1-based index)
    hh = find(u,bSpline.knots)
    if( hh == bSpline.knots.size ): # at end of knots vector
        hh = bSpline.knots.size - 4
        
    nn = 0
    if( hh == 1 ):
        nn = 1
    elif( hh == 2 ):
        nn = 2
    elif( hh == 3 ):
        nn = 3
    else:
        nn = 4
        
    ii = 0
    while( ii > -nn ): # i.e. for each contributing basis
        jj = hh + ii - 1 # subtracted 2 for C++ compatibility
        
        # Redefine k1 through k5 here
        k1 = bSpline.knots[jj]
        k2 = bSpline.knots[jj+1]
        k3 = bSpline.knots[jj+2]
        k4 = bSpline.knots[jj+3]
        k5 = bSpline.knots[jj+4]
        
        if( ii == 0 ): # calculate basis N1
            if( abs(k1-k2) <= 1e-12 ):
                Ncurrent = 0.
                dNducurrent = 0.
            else:
                Ncurrent = (u-k1)/(k4-k1)*(u-k1)/(k3-k1)*(u-k1)/(k2-k1)
                dNducurrent = -(3*pow(k1 - u,2))/((k1 - k2)*(k1 - k3)*       \
                  (k1 - k4))
        elif( ii == -1 ): # calculate basis N2
            if( abs(k2-k3) <= 1e-12):
                Ncurrent = 0.
                dNducurrent = 0.
            else:
                Ncurrent = (u-k1)/(k4-k1)*((u-k1)/(k3-k1)*(k3-u)/(k3-k2) +   \
                  (k4-u)/(k4-k2)*(u-k2)/(k3-k2)) + (k5-u)/(k5-k2)*(u-k2)     \
                  /(k4-k2)*(u-k2)/(k3-k2)
                dNducurrent = (((k1 - u)*(k3 - u))/((k1 - k3)*(k2 - k3)) +   \
                  ((k2 - u)*(k4 - u))/((k2 - k3)*(k2 - k4)))/(k1 - k4) +     \
                  pow(k2 - u,2)/((k2 - k3)*(k2 - k4)*(k2 - k5)) + ((k5 - u)  \
                  *(2*k2 - 2*u))/((k2 - k3)*(k2 - k4)*(k2 - k5)) +           \
                  (2*(k1 - u)*(k1*k2 - k3*k4 - k1*u - k2*u + k3*u + k4*u))/  \
                  ((k1 - k3)*(k1 - k4)*(k2 - k3)*(k2 - k4))
        elif( ii == -2 ): # calculate basis N3
            if( abs(k3-k4) <= 1e-12 ):
                Ncurrent = 0.
                dNducurrent = 0.
            else:
                Ncurrent = (u-k1)/(k4-k1)*(k4-u)/(k4-k2)*(k4-u)/(k4-k3) +    \
                  (k5-u)/(k5-k2)*((u-k2)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/    \
                  (k5-k3)*(u-k3)/(k4-k3))
                dNducurrent = - (((k2 - u)*(k4 - u))/((k2 - k4)*(k3 - k4)) + \
                  ((k3 - u)*(k5 - u))/((k3 - k4)*(k3 - k5)))/(k2 - k5) -     \
                  pow(k4 - u,2)/((k1 - k4)*(k2 - k4)*(k3 - k4)) -            \
                  ((k1 - u)*(2*k4 - 2*u))/((k1 - k4)*(k2 - k4)*(k3 - k4)) -  \
                  (2*(k5 - u)*(k2*k3 - k4*k5 - k2*u - k3*u + k4*u + k5*u))/  \
                  ((k2 - k4)*(k2 - k5)*(k3 - k4)*(k3 - k5))
        else: # calculate basis N4
            if( abs(k4-k5) <= 1e-12 ):
                Ncurrent = 0.
                dNducurrent = 0.
            else:            
                Ncurrent = (k5-u)/(k5-k2)*(k5-u)/(k5-k3)*(k5-u)/(k5-k4)
                dNducurrent = (3*pow(k5 - u,2))/((k2 - k5)*(k3 - k5)*        \
                  (k4 - k5))
        
        x += bSpline.coefs[0,jj]*Ncurrent
        y += bSpline.coefs[1,jj]*Ncurrent
        dxdu += bSpline.coefs[0,jj]*dNducurrent
        dydu += bSpline.coefs[1,jj]*dNducurrent
        
        ii -= 1
        
        if( ii < -4):
            break
        
    # END OF while( ii > -nn )
        
    return (x,y,dxdu,dydu)
    
# END OF uMap3

#==============================================================================
# Return y given x for a 3rd degree B-spline
#==============================================================================
def bSplineGeometry(x,bSpline):
    
    if( isinstance(x,float) ):
        if( x > (bSpline.coefs[0][-1]) ):
            x = bSpline.coefs[0][-1]
        elif( x < (bSpline.coefs[0][0]) ):
            x = bSpline.coefs[0][0]
            #raise ValueError("x is outside bounds of B-spline")
    
    if( isinstance(x,np.ndarray) ):
        nx = x.size # number of x
    else: # convert to numpy array for calculations
        x = np.array(([x]))
        nx = 1
        
    k = bSpline.knots.size # number of knots
    
    y = np.zeros(nx)
    dydx = np.zeros(nx)
    
    # Determine x-value at breaks
    xKnot = np.empty(k)
    for ii in range(0,k):
        xKnot[ii] = uMap3(bSpline.knots[ii],bSpline)[0]
    xKnot[k-4] += 1e-6 # so finding upper bound on u works (assumes same 4
                       # knots at end of knot vector)
    
    tolerance = 1e-6 # tolerance for Newton solver
    
    for ii in range(0,nx):
        
        # Determine lower and upper bounds on u
        seg = find(x[ii],xKnot)
        uLower = bSpline.knots[seg-1]
        uUpper = bSpline.knots[seg]
        
        # Pick a guess for u (a linear interpolation)
        u = (uLower + uUpper)/2
        
        # Calculate x and dxdu corresponding to u
        (xEst, temp1, dxduEst, temp2) = uMap3(u,bSpline)
        
        # Perform 1 Newton iteration
        if( dxduEst < tolerance ):
            uNew = 0.
        else:
            uNew = u - (xEst - x[ii])/dxduEst
            
        # Perform remaining Newton iterations
        counter = 0
        errorMeasure = abs((uNew - u)/uNew)
        while( errorMeasure > tolerance ):
            u = uNew
            
            (xEst, temp1, dxduEst, temp2) = uMap3(u,bSpline)
            
            if( dxduEst < 1e-12 ):
                uNew = u
            else:
                uNew = u - (xEst - x[ii])/dxduEst
                
            counter += 1
            
            if( counter > 20 ):
                break
            
            errorMeasure = abs((uNew-u)/uNew)
            
        # END OF while( errorMeasure > tolerance )
            
        u = uNew
        
        (xTemp, yTemp, dxduTemp, dyduTemp) = uMap3(u,bSpline)
        
        y[ii] = yTemp
        
        if( dxduTemp < tolerance ):
            dydx[ii] = 0.
        else:
            dydx[ii] = dyduTemp/dxduTemp
            
    # END OF for ii in range(0,nx)
            
    return (y, dydx)

# END OF bSplineGeometry
        
#==============================================================================
# Return y given x for a 3rd degree B-spline
#==============================================================================
def bSplineGeometryC(x,bSpline):
	
	if( isinstance(x,float) ):
		if( x > (bSpline.coefs[0][-1]) ):
		    x = np.array([bSpline.coefs[0][-1]])
		elif( x < (bSpline.coefs[0][0]) ):
		    x = np.array([bSpline.coefs[0][0]])
		else:
		    x = np.array([x])
 
 	coefs = list(bSpline.coefs.flatten());
		
	knots = list(bSpline.knots.flatten());
	
	x = list(x);
	
	y    = [];
	dydx = [];

	_meshutils_module.py_BSplineGeo3LowF (knots, coefs, x, y, dydx);
	
	return (np.asarray(y), np.asarray(dydx))
 
#==============================================================================
# Return y given x for a piecewise linear function
#==============================================================================
def piecewiseLinearGeometryC(x,piecewiseLinear):
	
	if( isinstance(x,float) ):
		if( x > piecewiseLinear.nodes[-1,0] ):
		    x = np.array([piecewiseLinear.nodes[-1,0]])
		elif( x < (piecewiseLinear.nodes[0,0]) ):
		    x = np.array([piecewiseLinear.nodes[0,0]])
		else:
		    x = np.array([x])
	
	x_nodes = list(piecewiseLinear.nodes[:,0].flatten());
	y_nodes = list(piecewiseLinear.nodes[:,1].flatten());
	
	x = list(x);
	
	y    = [];
	dydx = [];
	
	_meshutils_module.py_PiecewiseLinear (x_nodes, y_nodes, x, y, dydx); 
	
	return (np.asarray(y), np.asarray(dydx))
    
#==============================================================================
# Calculate volume of axisymmetric nozzle wall using trapezoidal integration
#==============================================================================
def wallVolume(innerWall,thickness):
    
    xVolume = np.linspace(0,innerWall.length,2000)
    volumeIntegrand = np.pi*innerWall.diameter(xVolume)*                     \
      thickness.radius(xVolume) + np.pi*thickness.radius(xVolume)**2
    volume = scipy.integrate.trapz(volumeIntegrand,xVolume)
    
    return volume
    
def wallVolume2Layer(innerWall,thickness1,thickness2):

    xVolume = np.linspace(0,innerWall.length,2000)
    volumeIntegrand = np.pi*innerWall.diameter(xVolume)*                     \
      (thickness1.radius(xVolume)+thickness2.radius(xVolume)) +              \
      np.pi*(thickness1.radius(xVolume)+thickness2.radius(xVolume))**2
    volume = scipy.integrate.trapz(volumeIntegrand,xVolume)
    
    return volume

#==============================================================================
# Calculate volume and mass of axisymmetric nozzle wall 10/12/16 RWF
#==============================================================================
    
# Take two functions of x, a lower one in a global (x,r) coordinate frame and 
# an upper one in a local (x,n) coordinate frame, where n is defined as the 
# normal to the lower function. The function below takes both functions and 
# returns the (x-r) coordinates of the upper function.
def localToGlobalCoordConversion(x,rLower,nUpper,drdxLower):

    theta = np.arctan2(drdxLower,1)
    
    # Check theta for discontinuities and smooth out if necessary
    thetaSlope = (theta[1:] - theta[0:-1])/(x[1:] - x[0:-1])
    
    slopeLimit = 300.
    neighborhood = 2 # look for large slopes within this neighborhood of steps
    indPos = -100 # count index of first large positive slope
    indPosCount = 0 # count number of steps of large positive slope after first
    indNeg = -100 # count index of first large negative slope
    indNegCount = 0 # count number of steps of large negative slope after first
    for i in range(len(thetaSlope)):
        if thetaSlope[i] > slopeLimit:
            if i < indPos + indPosCount + 1 + neighborhood: # point follows previous large positive slopes w/i a neighborhood
                indPosCount += i - indPos
            else: # point is the first of its large positive slope type
                indPos = i
                indPosCount = 0
        if thetaSlope[i] < -slopeLimit:
            if i < indNeg + indNegCount + 1 + neighborhood: # point follows previous large negative slopes w/i a neighborhood
                indNegCount += i - indNeg
            else: # point is the first of its large negative slope type
                indNeg = i
                indNegCount = 0
        # At a distance from steep phenomena, check for spike, and then replace with linear interpolation
        if i == max(indPos + indPosCount + 1 + neighborhood, indNeg + indNegCount + 1 + neighborhood):
            # positive spike
            if indNeg > indPos and indNeg - indPos - indPosCount <= neighborhood:
                #print 'found positive spike'
                s = indPos
                e = indNeg + indNegCount + 1
                theta[s:e] = np.interp(x[s:e],[x[s],x[e]],[theta[s],theta[e]])
            elif indPos > indNeg and indPos - indNeg - indNegCount <= neighborhood:
                #print 'found negative spike'
                s = indNeg
                e = indPos + indPosCount + 1
                theta[s:e] = np.interp(x[s:e],[x[s],x[e]],[theta[s],theta[e]])
    
    xTransform = x - nUpper*np.sin(theta)
    rTransform = rLower + nUpper*np.cos(theta)

    # Clean up transformed coordinates if necessary (remove overlap caused by 
    # 'kinks' in the shape)
    searchComplete = 0
    xa = np.max(xTransform)+1.
    while not searchComplete:
        decFlag = 0
        sitFlag = 0
        xa = np.max(xTransform)+1.
        for i in range(1,len(xTransform)):
            xDiff = xTransform[i] - xTransform[i-1]
            if xDiff < 0 and decFlag == 0: # i.e. x starts to decrease
                xa = xTransform[i-1] # point where x starts to decrease
                ra = rTransform[i-1] # radius where x starts to decrease
                ia = i-1
                decFlag = 1 # flag we are in a decreasing x portion of the shape
            elif decFlag == 1 and xDiff >= 0: # i.e. x starts to increase after having decreased
                xb = xTransform[i-1] # point where x starts to increase after decreasing
                rb = rTransform[i-1]
                ib = i-1
                decFlag = 2 # flag we are in an increasing portion of the shape after having been in a decreasing portion
            if xTransform[i] > xa: # just cut out offending kink entirely
                xTransform = np.hstack((xTransform[0:ia+1],xTransform[i:]))
                rTransform = np.hstack((rTransform[0:ia+1],rTransform[i:]))
                #print 'remedied basic kink'
                break
                
            if decFlag == 2:
                xcp = xTransform[i]
                rcp = rTransform[i]
                icp = i
                rI = np.interp(xcp,xTransform[0:ia+1],rTransform[0:ia+1])
                rII = np.interp(xcp,xTransform[ia:ib+1][::-1],rTransform[ia:ib+1][::-1]) 
                #print xTransform[ia:ib+1]
                #print rTransform[ia:ib+1]
                if rI > rcp and rcp > rII: # situation 3 before crossover
                    sitFlag = 3
                elif rII > rcp and rcp > rI: # situation 4 before crossover
                    sitFlag = 4
                elif rcp > rI and rI > rII and sitFlag == 3: # situation 3 after crossover
                    x4 = xcp; r4 = rcp
                    x3 = xTransform[i-1]; r3 = xTransform[i-1]
                    x1 = x3; x2 = x4;
                    r1 = np.interp(x1,xTransform[0:ia+1],rTransform[0:ia+1])
                    r2 = np.interp(x2,xTransform[0:ia+1],rTransform[0:ia+1])
                    m1 = (r2-r1)/(x2-x1); m2 = (r4-r3)/(x4-x3)
                    xi = (m1*x1 - m2*x3 + r3 - r1)/(m1 - m2)
                    ri = m1*(xi-x1) + r1                    
                    ii = len([q for q in xTransform[0:ia+1] if q < xi]) - 1
                    xTransform = np.hstack((xTransform[0:ii+1],xi,xTransform[icp:]))
                    rTransform = np.hstack((rTransform[0:ii+1],ri,rTransform[icp:]))
                    #print 'found intersection xi for situation 3'
                    break
                elif rII > rI and rI > rcp and sitFlag == 4: # situation 4 after crossover
                    x4 = xcp; r4 = rcp
                    x3 = xTransform[i-1]; r3 = xTransform[i-1]
                    x1 = x3; x2 = x4;
                    r1 = np.interp(x1,xTransform[0:ia+1],rTransform[0:ia+1])
                    r2 = np.interp(x2,xTransform[0:ia+1],rTransform[0:ia+1])
                    m1 = (r2-r1)/(x2-x1); m2 = (r4-r3)/(x4-x3)
                    xi = (m1*x1 - m2*x3 + r3 - r1)/(m1 - m2)
                    ri = m1*(xi-x1) + r1                    
                    ii = len([q for q in xTransform[0:ia+1] if q < xi]) - 1
                    xTransform = np.hstack((xTransform[0:ii+1],xi,xTransform[icp:]))
                    rTransform = np.hstack((rTransform[0:ii+1],ri,rTransform[icp:]))
                    #print 'found intersection xi for situation 4'
                    break      

            if i == len(xTransform)-1:
                searchComplete = 1

    # Ensure bounds are right for linear extrapolation if necessary
    if xTransform[0] > x[0]:
        xStart = x[0] - xTransform[0]
    else:
        xStart = xTransform[0] - (x[0]-xTransform[0])
        
    if xTransform[-1] < x[-1]:
        xEnd = x[-1] + (x[-1] - xTransform[-1])
    else:
        xEnd = xTransform[-1] + (xTransform[-1] - x[-1])
    mStart = (rTransform[1] - rTransform[0])/(xTransform[1] - xTransform[0])
    mEnd = (rTransform[-1] - rTransform[-2])/(xTransform[-1] - xTransform[-2])
    rStart = rTransform[0] - mStart*(0 - xStart)
    rEnd = rTransform[-1] + mEnd*(xEnd - xTransform[-1])
    
    xTransform = np.hstack((xStart,xTransform,xEnd))
    rTransform = np.hstack((rStart,rTransform,rEnd))
    
    rUpper = np.interp(x,xTransform,rTransform)
    
    return rUpper


# Used for 2D axisymmetric volume and mass calculations.
# Calculate the (x,r) coordinates for the inside wall shape, as well each layer
# thickness and the stringers (if any) given vector x for interpolation. As x 
# decreases in length and the number of layers increases, this function will 
# degrade in accuracy due to approximations of gradients and function values. 
def layerCoordinatesInGlobalFrame(nozzle,x):

    rList = list()
    if nozzle.stringers.n > 0 and nozzle.stringers.heightDefinition != 'EXTERIOR' \
        and nozzle.stringers.heightDefinition != 'BAFFLES_HEIGHT':
        N = len(nozzle.wall.layer) + 1
    else:
        N = len(nozzle.wall.layer)

    for i in range(N):
        if i == 0:
            lower = nozzle.wall.geometry # original normal coordinates (x-r)
            rLower = lower.radius(x)
            drdxLower = lower.radiusGradient(x)
            rList.append(rLower)
        else:
            lower = nozzle.wall.layer[i-1].thickness
            rLower = rUpper # update from last time
            # Approximate gradient from vector
            drdxTemp = (rLower[1:] - rLower[:-1])/(x[1:] - x[:-1])
            xTemp = (x[1:] - x[:-1])/2 + x[:-1]
            xTemp = np.hstack((nozzle.xinlet,xTemp,nozzle.xoutlet))
            drdxTemp = np.hstack((drdxTemp[0],drdxTemp,drdxTemp[-1]))
            drdxLower = np.interp(x,xTemp,drdxTemp)
        
        if i < len(nozzle.wall.layer): # wall layer
            upper = nozzle.wall.layer[i].thickness
        else: # stringers
            upper = nozzle.stringers.height
        nUpper = upper.radius(x)
    
        rUpper = localToGlobalCoordConversion(x,rLower,nUpper,drdxLower)
        rList.append(rUpper)
        
    return rList


# Approximate derivative at N points using data from those N points
def approxDerivative(x,y):

    dydxTmp = (y[1:] - y[:-1])/(x[1:] - x[:-1])
    #print dydxTmp
    #dydx = np.hstack((dydxTmp,dydxTmp[-1])) EVEN WORSE
    xTmp = (x[1:] - x[:-1])/2 + x[:-1]
    xTmp = np.hstack((x[0],xTmp,x[-1]))
    dydxTmp = np.hstack((dydxTmp[0],dydxTmp,dydxTmp[-1]))
    dydx = (dydxTmp[1:] + dydxTmp[:-1])/2
    # ANOTHER POTENTIAL SOLUTION
    # *if |dydx| < 1e-2 set to zero? although this is a hack
    #dydx = np.interp(x,xTmp,dydxTmp) Cannot use interpolation if x is not monotonic ascending

    return dydx


# Used for 3D non-axisymmetric volume and mass calculations.
# Calculate the (theta,r) coordinates for the inside wall shape, as well as
# each layer exterior at a given x-coordinate. 
def radialCoordinatesInGlobalFrame(nozzle,x,theta):

    r1 = nozzle.wall.majoraxis.geometry.radius(x)[0]
    r2 = nozzle.wall.minoraxis.geometry.radius(x)[0]
    z0 = nozzle.wall.centerline.geometry.radius(x)[0]
    thetain = nozzle.wall.shovel_start_angle
    thetaout = nozzle.wall.shovel_end_angle
    xin = nozzle.wall.centerline.geometry.xstart
    xout = nozzle.wall.centerline.geometry.xend    
    alpha = (x - xin)/(xout - xin)
    thetacut = alpha*thetaout + (1-alpha)*thetain

    ic = np.searchsorted(theta,thetacut) # all theta indices prior to this one
                                         # correspond to a shape with no
                                         # shovel effect
    rList = list()
    N = len(nozzle.wall.layer)
    for i in range(N):

        if i == 0: # inner wall

            rLower = np.zeros((len(theta),)) # store r^2
            yTmp = np.zeros((len(theta),))
            zTmp = np.zeros((len(theta),))
            thetaTmp = np.zeros((len(theta),)) # store temporary theta

            # Upper half & lower half (y,z) coordinates
            yTmp[0:ic+1] = r1*np.sin(theta[0:ic+1])
            zTmp[0:ic+1] = r2*np.cos(theta[0:ic+1])
            yTmp[ic+1:]  = r1*np.sin(theta[ic+1:])
            zTmp[ic+1:]  = alpha*r2*np.cos(theta[ic+1:]) + \
                           (1-alpha)*r2*np.cos(theta[ic+1:])

            # Form radius of wall
            rLower[:] = np.sqrt(yTmp**2 + zTmp**2)
            thetaTmp[:] = np.arctan2(yTmp,zTmp)
            rLower = np.interp(theta,thetaTmp,rLower)

            # plt.plot(yTmp,zTmp)
            # plt.plot(rLower*np.sin(theta),rLower*np.cos(theta))
            # plt.axis('equal')
            # plt.show()

            rList.append(rLower)
        else: # exterior of other walls
            rLower = rUpper # update from last time

        yLower = rLower*np.sin(theta)
        zLower = rLower*np.cos(theta)

        # Approximate derivative for now
        dzdyLower = approxDerivative(yLower,zLower)
        normalAngle = np.arctan(-1./dzdyLower) # angle of normal vector
        # Artifically fix normal angles at top and bottom of nozzle to ensure 
        # angles point outward. This may occur near top and bottom of nozzle.
        # So long as there are no steep gradients in layer thickness 
        # parameterization, the below method does not pose a problem.
        normalAngle[0:int(len(theta)/6)] = np.abs(normalAngle[0:int(len(theta)/6)])
        normalAngle[len(theta)-int(len(theta)/6):] = -np.abs(normalAngle[len(theta)-int(len(theta)/6):])
        if nozzle.wall.layer[i].thickness.type == 'piecewise-bilinear':
            thickness = nozzle.wall.layer[i].thickness.height(x,theta*180/np.pi)
        else:
            thickness = nozzle.wall.layer[i].thickness.radius(x)
        yUpper = yLower + thickness*np.cos(normalAngle)
        zUpper = zLower + thickness*np.sin(normalAngle)
        rTmp = np.sqrt(yUpper**2 + zUpper**2)
        thetaTmp = np.arctan2(yUpper,zUpper)
        rUpper = np.interp(theta,thetaTmp,rTmp)

        rList.append(rUpper)

    return rList

    
# Calculate and return volume and mass of nozzle and structure (stringers &
# baffles) given fully parameterized
# nozzle. Assumes piecewise-linear layer thickness definitions. On the baseline
# geometry, using n = 1e4 results in ~0.02% error in gradients with respect
# to some inner wall B-spline coefficients (compared to the derivatives
# estimated using 1e7, where convergence of the mass calculation was observed).
# The rule of thumb above is for the 2D parameterization only.
def calcVolumeAndMass(nozzle):
    
    if nozzle.dim == '3D':
        n = 1000 # previously 1e4
    else:
        n = 1e4
    x = np.linspace(nozzle.xinlet,nozzle.xoutlet,n)
    # Pick x smartly
    xHit = set()
    if( nozzle.wall.geometry.type == 'B-spline' ):
        for i in range(len(nozzle.wall.geometry.coefs[0,:])):
            xHit.add(nozzle.wall.geometry.coefs[0,i])
    for i in range(len(nozzle.wall.layer)):
        for j in range(len(nozzle.wall.layer[i].thicknessNodes[:,0])):
            xHit.add(nozzle.wall.layer[i].thicknessNodes[j,0])
    for i in range(nozzle.baffles.n):
        xHit.add(nozzle.baffles.location[i])
    xHit = list(xHit)
    xHit.sort()
    x = np.array([])
    for i in range(len(xHit)-1):
        m = int(n*(xHit[i+1] - xHit[i])/xHit[-1])
        xTemp = np.linspace(xHit[i],xHit[i+1],m)
        if len(xTemp) <= 1:
            x = np.hstack((x[:-1],[xHit[i],xHit[i+1]]))
        else:
            x = np.hstack((x[:-1],xTemp))

    if nozzle.dim == '3D':

        # Pick theta for integration in Y-Z plane, theta = 0 corresponds
        # to Z axis, not Y axis
        n2 = 100
        theta = np.linspace(0.,np.pi,n2)
        deltatheta = theta[1]-theta[0] # equally-spaced

        # Now calculate volume of shell formed by inner wall and each layer's
        # exterior. A shell symmetric across X-Z plane is assumed.
        Vshell = np.array([0., 0., 0., 0., 0., 0.])
        AbaffleInsideEdge = list()
        zTop = np.zeros((len(x),)) # approximate Z-coord of top surface of nozzle
        zBot = np.zeros((len(x),)) # approximate Z-coord of bottom surface of nozzle

        for i in range(len(x)): # for each x-station

            if i == 0:
                deltax = (x[1]-x[0])/2.
            elif i < len(x)-1:
                deltax = (x[i+1]-x[i-1])/2.
            else: # i == len(x)-1
                deltax = (x[-1]-x[-2])/2.

            xc = x[i]

            rList = radialCoordinatesInGlobalFrame(nozzle,xc,theta)

            # Save approximate z-coordinate of top and bottom of nozzle for 
            # stringer volume and mass calculations
            zTop[i] = rList[-1][0]
            zBot[i] = -rList[-1][-1]

            for j in range(len(Vshell)):
                # Assumes symmetric shell across X-Z plane
                # Assumes vertical slices in Y-Z plane...a correction term for
                # deltax in terms of deltatheta should probably be used or a 
                # more accurate integration scheme should be used.
                Vshell[j] = Vshell[j] + np.sum(np.power(rList[j],2)*deltatheta)*deltax

            # Record area for baffle calculations
            if xc in nozzle.baffles.location:
                AbaffleInsideEdge.append(np.sum(np.power(rList[j],2)*deltatheta))

            # for item in rList:
            #     ytmp = item*np.sin(theta)
            #     ztmp = item*np.cos(theta)
            #     plt.plot(ytmp,ztmp)
            # plt.axis('equal')
            # plt.show()

        # Calculate layer volume
        V = list(Vshell[1:] - Vshell[0:-1])

        # Obtain layer mass
        m = list()
        for i in range(len(nozzle.wall.layer)):
            m.append(nozzle.wall.layer[i].material.getDensity()*V[i])

        # Calculate volume and mass for stringers
        for i in range(nozzle.stringers.n):
            angle = nozzle.stringers.thickness[i].angle # measured from Y axis CCW
            if np.abs(angle-90.) < 1e-6: # 90deg, top of nozzle
                dz = nozzle.exterior.geometry['top'].coord(x,np.pi/2)[1] - zTop
                for j in range(len(dz)):
                    if dz[j] < 0:
                        print "Nozzle wall thickness approximation intersects exterior at top"
                        dz[j] = 0.
            elif np.abs(angle-270.) < 1e-6: # 270deg, bottom of nozzle
                dz = zBot - nozzle.exterior.geometry['bottom'].coord(x,-np.pi/2)[1]
                for j in range(len(dz)):
                    if dz[j] < 0:
                        print "Nozzle wall thickness approximation intersects exterior at bottom"
                        dz[j] = 0.                
            else:
                raise NotImplementedError("Volume and mass calculation for " + \
                    "stringers in 3D param with angles other than 90 and 270" + \
                    "degrees is not implemented.")

            dt = nozzle.stringers.thickness[i].radius(x)
            dx = np.hstack(((x[1]-x[0])/2,(x[2:]-x[0:-2])/2.,(x[-1]-x[-2])/2.))

            V.append(np.sum(dz*dt*dx))
            m.append(V[-1]*nozzle.stringers.material.getDensity())

        # Calculate volume and mass for baffles
        for i in range(nozzle.baffles.n):
            xTmp = nozzle.baffles.location[i]
            yTmp = nozzle.baffles.halfWidth
            zTmp1 = nozzle.exterior.geometry['top'].z(xTmp,yTmp)
            zTmp2 = nozzle.exterior.geometry['bottom'].z(xTmp,yTmp)
            Atmp = (zTmp1-zTmp2)*yTmp*2 # approximate baffle outer edge as rectangle
                # This approximation is justified since area between rectangle
                # and true baffle outer edge is the same for each baffle by definition.
                # As a result this will just underapproximate baffle mass a slight amount,
                # by the same amount each time.
            Abaffle = Atmp - AbaffleInsideEdge[i]
            Vbaffle = Abaffle*nozzle.baffles.thickness[i]
            V.append(Vbaffle)
            m.append(Vbaffle*nozzle.baffles.material.getDensity())

    else: # 2D or 1D parameterization implies axisymmetric nozzle

        # Obtain radii of layers
        radiusList = layerCoordinatesInGlobalFrame(nozzle,x)
        
        # Calculate volume and mass for nozzle layers
        s = list()
        V = list()
        for i in range(len(nozzle.wall.layer)):
            midpoint = (radiusList[i+1] + radiusList[i])/2
            ds = np.sqrt( (midpoint[1:] - midpoint[:-1])**2 + (x[1:] - x[:-1])**2 )
            xMid = (x[1:] + x[:-1])/2
            mMid = np.interp(xMid,x,midpoint)
            dV = 2*np.pi*mMid*nozzle.wall.layer[i].thickness.radius(xMid)*ds
            s.append(np.sum(ds))
            V.append(np.sum(dV))
        
        m = list()
        for i in range(len(nozzle.wall.layer)):
            m.append(nozzle.wall.layer[i].material.getDensity()*V[i])

        # Calculate volume and mass for stringers
        # if nozzle.stringers.n > 0:
        #     deltaR = radiusList[-1] - radiusList[-2]
        #     xMid = (x[1:] + x[:-1])/2
        #     dr = np.interp(xMid,x,deltaR)
        #     dw = nozzle.stringers.thickness.radius(xMid)
        #     dx = x[1:] - x[:-1]
        #     dV = dr*dw*dx
        #     V.append(nozzle.stringers.n*np.sum(dV))
        #     m.append(V[-1]*nozzle.stringers.material.getDensity())
        if( nozzle.stringers.heightDefinition == 'EXTERIOR' \
            or nozzle.stringers.heightDefinition == 'BAFFLES_HEIGHT' ): 
            
            # Assume one stringer on top and one on bottom
            for i in range(2):
                if i == 0: # 90deg, top of nozzle
                    dz = nozzle.exterior.geometry['top'].coord(x,np.pi/2)[1] - radiusList[-1]
                    for j in range(len(dz)):
                        if dz[j] < 0:
                            print "Nozzle wall thickness approximation intersects exterior at top"
                            dz[j] = 0.
                else: # i == 1, 270deg, bottom of nozzle
                    dz = radiusList[-1] - nozzle.exterior.geometry['bottom'].coord(x,-np.pi/2)[1]
                    for j in range(len(dz)):
                        if dz[j] < 0:
                            print "Nozzle wall thickness approximation intersects exterior at bottom"
                            dz[j] = 0.           

                # THIS IS A HACK
                try:
                    dt = nozzle.stringers.thickness[i].radius(x)
                except:
                    dt = nozzle.stringers.thickness.radius(x)
                dx = np.hstack(((x[1]-x[0])/2,(x[2:]-x[0:-2])/2.,(x[-1]-x[-2])/2.))

                V.append(np.sum(dz*dt*dx))
                m.append(V[-1]*nozzle.stringers.material.getDensity())

        # Each stringer has individual height
        elif nozzle.stringers > 0:

            deltaR = radiusList[-1] - radiusList[-2]
            xMid = (x[1:] + x[:-1])/2
            dr = np.interp(xMid,x,deltaR)
            dw = nozzle.stringers.thickness.radius(xMid)
            dx = x[1:] - x[:-1]
            dV = dr*dw*dx
            V.append(nozzle.stringers.n*np.sum(dV))
            m.append(V[-1]*nozzle.stringers.material.getDensity())   

        # Calculate volume and mass for baffles
        # for i in range(nozzle.baffles.n):
        #     # baffles connect to outside of nozzle wall
        #     if nozzle.stringers.n > 0: # stringer radius is returned
        #         rInner = np.interp(nozzle.baffles.location[i],x,radiusList[-2])
        #     else: # no stringer radius is returned; last entry is last layer
        #         rInner = np.interp(nozzle.baffles.location[i],x,radiusList[-1])
        #     rOuter = rInner + nozzle.baffles.height[i]
        #     V.append(np.pi*(rOuter**2 - rInner**2)*nozzle.baffles.thickness[i])
        #     # All baffles have the same fixed ratio panel material
        #     m.append(V[-1]*nozzle.baffles.material.getDensity())
        AbaffleInsideEdge = list()
        for i in range(nozzle.baffles.n):
            if nozzle.stringers.n > 0: # stringer radius is returned
                rInner = np.interp(nozzle.baffles.location[i],x,radiusList[-2])
            else: # no stringer radius is returned; last entry is last layer
                rInner = np.interp(nozzle.baffles.location[i],x,radiusList[-1]) 
            AbaffleInsideEdge.append(np.pi*rInner**2)
        for i in range(nozzle.baffles.n):
            xTmp = nozzle.baffles.location[i]
            yTmp = nozzle.baffles.halfWidth
            zTmp1 = nozzle.exterior.geometry['top'].z(xTmp,yTmp)
            zTmp2 = nozzle.exterior.geometry['bottom'].z(xTmp,yTmp)
            Atmp = (zTmp1-zTmp2)*yTmp*2 # approximate baffle outer edge as rectangle
                # This approximation is justified since area between rectangle
                # and baffle outer edge is the same for each baffle by definition.
                # As a result this will just underapproximate baffle mass a slight amount.
            Abaffle = Atmp - AbaffleInsideEdge[i]
            Vbaffle = Abaffle*nozzle.baffles.thickness[i]
            V.append(Vbaffle)
            m.append(Vbaffle*nozzle.baffles.material.getDensity())

    
    return V, m
    

# Calculate and return forward finite difference gradients of nozzle mass
def calcMassGradientsFD(nozzle,fd_step,components='all'):
    
    if components == 'all': # calculate mass of nozzle wall & structure
        mass = nozzle.responses['MASS'];
    elif components == 'wall-only': # calculate mass of wall layers only
        mass = nozzle.responses['MASS_WALL_ONLY'];
    
    # For parallel computation
    #print nozzle.partitions
    
    # Perform serial forward finite difference on mass
    dmdx = [];
    
    # For each design variable
    
    #print len(nozzle.dvList);
    #sys.exit(1);
    

    for i in range(len(nozzle.derivativesDV)):
        
        nozzle2 = copy.deepcopy(nozzle);
        
        if isinstance(fd_step,list):
            dx = fd_step[nozzle.derivativesDV[i]-1];
        else:
            dx = fd_step;
            
        nozzle2.dvList[nozzle.derivativesDV[i]-1] += dx;
        nozzle2.UpdateDV(output='quiet');
        nozzle2.SetupWall(output='quiet');
        
        volume2, mass2 = calcVolumeAndMass(nozzle2);
        if components == 'all': # calcuate mass of nozzle wall & structure
            mass2 = np.sum(mass2);
            volume2 = np.sum(volume2);
        elif components == 'wall-only': # calculate mass of wall layers only
            n_layers = len(nozzle.wall.layer);
            mass2 = np.sum(mass2[:n_layers]);
            volume2 = np.sum(volume2[:n_layers]);
        
        dmdxLocal = (mass2-mass)/dx;
        
        dmdx.append((mass2-mass)/dx);

    	
		#print "DV %d : copy %lf sec , updateDV %lf sec , setupWall %lf sec , calcvol %lf sec , rest %lf " % (i, t1-t0, t2-t1, t3-t2, t4-t3, t5-t4 )

    return dmdx


