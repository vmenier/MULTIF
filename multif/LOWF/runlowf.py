"""
Perform quasi-1D area-averaged Navier-Stokes analysis on axisymmetric nozzle

Rick Fenrich 9/24/17
"""

import numpy as np
import scipy.optimize
import scipy.integrate
import sys, os

import quasi1dnozzle

from .. import nozzle as nozzlemod
from multif.MEDIUMF.AEROSpostprocessing import PostProcess as AEROSPostProcess

try:
    from multif.MEDIUMF.runAEROS import *
except ImportError as e: 
    print 'Error importing all functions from runAEROS in runlowf.py.'
    print e
    print

try:
	from matplotlib import pyplot as plt
except ImportError as e:
    print 'Error importing matplotlib in runlowf.py.'
    print e
    print

#==============================================================================
# Sutherland's Law of dynamic viscosity of air
#==============================================================================
def dynamicViscosity(T):
    mu = 1.716e-5*(T/273.15)**1.5*(273.15 + 110.4)/(T + 110.4)
    return mu
   
#==============================================================================
# Area-Mach function from 1-D mass conservation equations
#==============================================================================
def areaMachFunc(g,M):
    a = ((g+1)/2)**((g+1)/(2*(g-1)))*M/(1+(g-1)*M**2/2)**((g+1)/(2*(g-1)))
    return a

#==============================================================================
# Mass flow-rate
#==============================================================================
def massFlowRate(fluid,Pstag,Area,Tstag,M):
    gam = fluid.gam
    R = fluid.R;
    mdot = (gam/((gam+1)/2)**((gam+1)/(2*(gam-1))))*Pstag*Area*              \
      areaMachFunc(gam,M)/np.sqrt(gam*R*Tstag)
    return mdot

#==============================================================================
# Determine state of nozzle stagnation pressure and temperatures at throat and
# exit, nozzle geometry, and pressure ratio between reservoir and atmosphere
#==============================================================================
def nozzleState(nozzle,pressureRatio,PsT,TsT,PsE,TsE):
    gam = nozzle.fluid.gam
    Athroat = nozzle.wall.geometry.area(nozzle.wall.geometry.xThroat)
    Aexit = nozzle.wall.geometry.area(nozzle.wall.geometry.length)
    
    rootFunc = lambda x: ((gam+1)/2)**((gam+1)/(2*(gam-1)))*x/(1+(gam-1)*    \
      x**2/2)**((gam+1)/(2*(gam-1)))*np.sqrt(TsT/TsE)*(PsE/PsT) - Athroat/Aexit
      
    MsubsonicCritical = scipy.optimize.fsolve(func=rootFunc,x0=0.5,xtol=1e-12)
    MsupersonicCritical = scipy.optimize.fsolve(func=rootFunc,x0=2.,xtol=1e-12)
    
    PtRatioSubsonic = (1 + (gam-1)*MsubsonicCritical**2/2)**(gam/(gam-1))
    PtRatioSupersonic = (1 + (gam-1)*MsupersonicCritical**2/2)**(gam/(gam-1))
      
    MbehindShock = np.sqrt((1 + (gam-1)*MsupersonicCritical**2/2)/           \
      (gam*MsupersonicCritical**2 - (gam-1)/2))
    normalShockPtRatio = (Aexit/Athroat)*((gam+1)/2)**((gam+1)/(2*(gam-1)))  \
      *MbehindShock*np.sqrt(1 + (gam-1)*MbehindShock**2/2)
      
    deltaPtRatio = 0.05 # a tolerance to determine if fully expanded flow occur
    Mshock = 0.
    if( pressureRatio <= 1 ):
        status = "no flow"
    elif( pressureRatio < PtRatioSubsonic ):
        status = "subsonic"
        Mexit = np.sqrt(2/(gam-1))*np.sqrt(pressureRatio**((gam-1)/gam) - 1)
        Mshock = Mexit # return Mexit in place of Mshock
    elif( pressureRatio < normalShockPtRatio):
        status = "shock"
        
        rootFunc2 = lambda x: ((gam+1)/2)**((gam+1)/(2*(gam-1)))*x*          \
          np.sqrt(1 + (gam-1)*x**2/2) - pressureRatio*(Athroat/Aexit)*       \
          np.sqrt(TsE/TsT)
        Mexit = scipy.optimize.fsolve(func=rootFunc2,x0=0.5,xtol=1e-12)
        PtRatio = (Athroat/Aexit)*np.sqrt(TsE/TsT)/(((gam+1)/2)**((gam+1)/   \
          (2*(gam-1)))*Mexit/(1+(gam-1)*Mexit**2/2)**((gam+1)/(2*(gam-1))))
        
        rootFunc3 = lambda x: (((gam+1)*x**2/2)/(1 + (gam-1)*x**2/2))**      \
          (gam/(gam-1))*(((gam+1)/2)/(gam*x**2 - (gam-1)/2))**(1/(gam-1))    \
          - PtRatio
        Mshock = scipy.optimize.fsolve(func=rootFunc3,x0=2,xtol=1e-12)
        #MpostShock = np.sqrt((1 + ((gam-1)/2)*Mshock**2)/(gam*Mshock**2 -    \
        #  (gam-1)/2))
    elif( pressureRatio < PtRatioSupersonic - deltaPtRatio ):
        status = "overexpanded"
    elif( pressureRatio < PtRatioSupersonic + deltaPtRatio ):
        status = "fully expanded"
    else:
        status = "underexpanded"
        
    return (status,Mshock)
    
#==============================================================================
# Solve for location where M = 1 (location of apparent throat)
#==============================================================================
def findApparentThroat(nozzle,tol,params):

    n = 10000
    x = np.linspace(nozzle.xinlet,nozzle.xoutlet,n)
    
    # First, determine approximate throat location by minimum area point
    minInd = np.argmin(nozzle.wall.geometry.area(x))
    xMinArea = x[minInd]    
    
    # Next, search in area of throat for exact location.
    # Find axial location in nozzle where dM^2/dx first becomes positive given
    # M = 1+delta. This is the location of the apparent throat. If dM^2/dx < 0
    # when M = 1+delta, then M will have a tendency to decrease and will choke
    # the flow again. We require M to increase after this point.
    dM2dx = np.zeros((n,))
    #params = (xInterp,Cf,Tstag,dTstagdx);    
    xRight = nozzle.xoutlet
    
    # Search to the right:
    frac = 0.1 # fraction of nozzle length to search around min area point
    for i in range(minInd,int(minInd+frac*n)):
        dM2dx[i] = dM2dxForward(x[i], 1.0001**2, nozzle.fluid.gam, \
           nozzle.wall.geometry, params)
        if dM2dx[i] > 0: # reached point of increase for Mach at M = 1
            if i == 0:
                xRight = x[i]
                break
            else:
                xRight = (x[i]+x[i-1])/2
                break
            
    xThroat = xRight
    
#    plt.plot(x,dM2dx,'b')
#    plt.plot(x,np.zeros((n,)),'r')
#    plt.show()
#    sys.exit()

#    relTol = tol["solverApparentThroatLocation"]
#
#    dMdxCoeffFunc = lambda x: -nozzle.wall.geometry.areaGradient(x)/         \
#      nozzle.wall.geometry.area(x) + 2*gam*np.interp(x,xInterp,Cf)/          \
#      nozzle.wall.geometry.diameter(x) + (1+gam)*np.interp(x,xInterp,        \
#      dTstagdx)/(2*np.interp(x,xInterp,Tstag))
#      
#    # Find sign changes in dMdxCoeffFunc
#    xFind = np.linspace(0.,nozzle.wall.geometry.length-1e-4,1000.)
#    coeffFind = dMdxCoeffFunc(xFind)
#    coeffFindSign = np.sign(coeffFind)
#    coeffFindSignChange = ((np.roll(coeffFindSign,1) - coeffFindSign)        \
#      != 0).astype(int)
#    coeffFindSignChange[0] = 0
#    signChangeLocations = np.nonzero(coeffFindSignChange)[0]
#    
#    indKeep = 0
#    if( signChangeLocations.size == 0 ):
#        throatGuess = nozzle.wall.geometry.findMinimumRadius()[0]
#    else:
#        minInd = np.argmin(nozzle.wall.geometry.area(xFind[signChangeLocations]))
#        throatGuess = xFind[signChangeLocations[minInd]]
#        
#        # Check to make sure each following possible throat is far enough away
#        ind = signChangeLocations[minInd]
#        indKeep = ind
#        for ii in range(minInd+1,signChangeLocations.size):
#            currentInd = signChangeLocations[ii]
#            dx = xFind[currentInd] - xFind[ind]
#            
#            #dAdxbar = nozzle.wall.geometry.area(xFind[ind])
#            dAdxbarEst = (nozzle.wall.geometry.area(xFind[currentInd]) -     \
#              nozzle.wall.geometry.area(xFind[ind]))/dx
#            
#            dTstagdxbar = np.interp(xFind[currentInd],xInterp,dTstagdx)
#            Abar = nozzle.wall.geometry.area(xFind[ind])
#            Tstagbar = np.interp(xFind[ind],xInterp,Tstag)
#            Cfbar = np.interp(xFind[ind],xInterp,Cf)
#            Dbar = nozzle.wall.geometry.diameter(xFind[ind])
#            
#            RHS = dTstagdxbar*Abar*(1+gam)/(2*Tstagbar) + 2*gam*Abar*Cfbar/Dbar
#            
#            if( dAdxbarEst <= RHS ):
#                throatGuess = xFind[currentInd]
#                indKeep = currentInd
#            
#        # END OF for ii in range(minInd,signChangeLocations.size)
#            
#    xApparentThroat = scipy.optimize.fsolve(func=dMdxCoeffFunc,              \
#      x0=throatGuess,xtol=relTol)[0]
#    
#    # Perform some error checking for the apparent throat location
#    if( np.isnan(xApparentThroat) or                                         \
#      xApparentThroat > nozzle.wall.geometry.length):
#        xApparentThroat = nozzle.wall.geometry.length
#        print "Mach = 1 at exit\n"
#    
#    if( nozzle.wall.geometry.diameter(xApparentThroat) >                      \
#      nozzle.wall.geometry.diameter(nozzle.wall.geometry.length)):
#        xApparentThroat = nozzle.wall.geometry.length
#        
#    if( indKeep != len(xFind) and indKeep != 0 and 
#      (xApparentThroat > xFind[indKeep+1] or xApparentThroat < xFind[indKeep-1]) ):
#      print 'scipy fsolve found wrong apparent throat at %f' % xApparentThroat
#      # Reset apparent throat to a linear interpolation between the 2 points
#      # Sign change index of 1 implies sign change between indices 0 and 1
#      slope = (coeffFind[indKeep] - coeffFind[indKeep-1])/(xFind[indKeep] -   \
#        xFind[indKeep-1])
#      xApparentThroat = xFind[indKeep] - coeffFind[indKeep]/slope
#      print 'apparent throat interpolated and set to %f' % xApparentThroat
#    
#    print
#    print 'Throat originally at: %f' % xApparentThroat
#    minInd = np.argmin(nozzle.wall.geometry.area(xFind))
#    throatGuess = xFind[minInd]
#    xApparentThroat = throatGuess
#    print 'Throat at: %f' % xApparentThroat
#    print
    
    return xThroat

#==============================================================================
# Define non-ideal quasi-1D equations of motion for forward integration    
#==============================================================================
def dM2dxForward(x,M2,gam,geo,(xInterp,CfInterp,TstagInterp,dTstagdxInterp)):
    dAdx = geo.areaGradient(x)
    A = geo.area(x)
    D = geo.diameter(x)
    Cf = np.interp(x,xInterp,CfInterp)
    Tstag = np.interp(x,xInterp,TstagInterp)
    dTstagdx = np.interp(x,xInterp,dTstagdxInterp)
    dM2dx = (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx/A + 2*gam*M2*Cf/D +        \
      (1+gam*M2)*dTstagdx/(2*Tstag))    
    return dM2dx

#==============================================================================
# Define non-ideal quasi-1D equations of motion for backward integration    
#==============================================================================    
def dM2dxBackward(x,M2,xRef,gam,geo,(xInterp,CfInterp,TstagInterp,dTstagdxInterp)):
    dAdx = geo.areaGradient(xRef-x)
    A = geo.area(xRef-x)
    D = geo.diameter(xRef-x)
    Cf = np.interp(xRef-x,xInterp,CfInterp)
    Tstag = np.interp(xRef-x,xInterp,TstagInterp)
    dTstagdx = np.interp(xRef-x,xInterp,dTstagdxInterp)
    dM2dx = -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx/A + 2*gam*M2*Cf/D +       \
      (1+gam*M2)*dTstagdx/(2*Tstag))
    return dM2dx
    
#==============================================================================
# Check ODE integration for problems such as hitting a singularity
#==============================================================================
#def checkODEintegration(x,y):
#    r = 2. # maximum allowable multiplicative factor for slopes between 2
#            # adjacent steps
#    tol = 0.01 # tolerance to determine whether y is close enough to zero
#    
#    # calculate slopes for each step
#    m = (y - np.roll(y,1))/(x - np.roll(x,1))
#    m[0] = m[1] # assume problem doesn't occur at first step
#
#    v1 = [ii for ii in range(1,m.size) if abs(m[ii]) > r*abs(m[ii-1]) and (m[ii-1] > tol or m[ii] > tol)]
#    
#    if( len(v1) != 0 ):
#        minIndex = v1[0] - 1 # last good value before solver went crazy
#    else:
#        minIndex = -1
#
#    return minIndex
    
#==============================================================================
# Integrate ODE built with scipy.integrate.ode and perform simple event check.
# Provide solution at equally spaced points. Integration will end when event is
# detected.
#==============================================================================
def integrateODEwithEvents(ode,dt,tfinal,ycrit,direction):
    
    if( isinstance(dt,float) ): # uniform spacing
        if( direction == "b" ):
            x = np.linspace(tfinal,ode.t,np.round((ode.t-tfinal)/dt))
        else:
            x = np.linspace(ode.t,tfinal,np.round((tfinal-ode.t)/dt))
    else:
        raise TypeError(("integration for y(t) at uniform intervals not ",
                        "enabled, dt must be a float"))
        # if non-uniform x-spacing is desired, x must be flipped for backwards
        # integration
        # for backwards integration: x = abs(x - max(x)) + min(x)
    
    ii = 0
    y = np.zeros(x.size)
    delta = 1e-12 # gives distance away from event when solver terminates
    
    if( direction == "b" ): # backward integration
    
        # translate into a forward integration problem
        tTemp = ode.t
        ode.set_initial_value(ode.y,tfinal)
        tfinal = tTemp
    
        if( ode.y > ycrit ):
            
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                #print yTemp
                if( yTemp < ycrit + delta ):
                    y[x.size-ii-1] = yTemp
                    break # critical y reached
                    
                y[x.size-ii-1] = yTemp
                ii +=1
                    
        elif( ode.y <= ycrit ):
        
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                #print yTemp
                if( yTemp > ycrit - delta ):
                    y[x.size-ii-1] = yTemp
                    break # critical y reached
                
                y[x.size-ii-1] = yTemp
                ii +=1       
    
    else: # forward integration
    
        if( ode.y > ycrit ):
            
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                #print yTemp
                if( yTemp < ycrit + delta ):
                    y[ii] = yTemp
                    break # critical y reached
                    
                y[ii] = yTemp
                ii +=1
                    
        elif( ode.y <= ycrit ):
        
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                #print yTemp
                if( yTemp > ycrit - delta ):
                    y[ii] = yTemp
                    break # critical y reached
                
                y[ii] = yTemp
                ii +=1  
    
    return (x,y,ii)
    
#==============================================================================
# Integrate subsonic flow through an axial nozzle geometry
#==============================================================================
def integrateSubsonic(nozzle,tol,params,xThroat,nPartitions):
    
    # Use inlet Mach number provided by user, if available
    if( hasattr(nozzle.inlet, "mach") ):
        print "Using prescribed inlet Mach number\n"
        M0 = nozzle.inlet.mach
        
        f = scipy.integrate.ode(dM2dxForward,jac=None).set_integrator(       \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        f.set_initial_value(M0**2,0.)
        f.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
        dt = nozzle.wall.geometry.length/nPartitions
        tfinal = nozzle.wall.geometry.length
        (xIntegrate,M2,eventIndex) = integrateODEwithEvents(f,dt,tfinal,1.,"f")
        
        if( eventIndex != xIntegrate.size ):
            raise ValueError(("Integration terminated early: prescribed ",
                              "inlet Mach number is too large"))
    
    # Else if nozzle converges only, assume choked flow at the exit    
    elif( nozzle.wall.geometry.length - xThroat < 1e-12 ):
        M0 = 0.9999 # start integration from this Mach number
        
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(M0**2,nozzle.wall.geometry.length)
        b.set_f_params(nozzle.wall.geometry.length,nozzle.fluid.gam,         \
          nozzle.wall.geometry,params)
        dt = nozzle.wall.geometry.length/nPartitions
        tfinal = 0.
        (xIntegrate,M2,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
                
        if( eventIndex != xIntegrate.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating backwards from the exit"))
    
    # Else, assume nozzle is choked at throat, integrate forwards & backwards
    else:
        UpperM = 0.9999 # start integration at this Mach number for aft portion
        LowerM = 0.9999 # start integ. at this Mach number for fore portion
        dx = 1e-5 # 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001
        
        # Integrate forward from throat
        f = scipy.integrate.ode(dM2dxForward,jac=None).set_integrator(       \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        f.set_initial_value(UpperM**2,xThroat+dx/2)
        f.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round((1-xThroat/nozzle.wall.geometry.length)*nPartitions)
        dt = (nozzle.wall.geometry.length-xThroat-dx/2)/nP
        tfinal = nozzle.wall.geometry.length
        (xF,M2F,eventIndex) = integrateODEwithEvents(f,dt,tfinal,1.,"f")
        
        if( eventIndex != xF.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating forwards from the throat"))
        
        # Integrate backward from throat
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(LowerM**2,xThroat-dx/2)
        b.set_f_params(xThroat,nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round(xThroat/nozzle.wall.geometry.length*nPartitions)
        dt = (xThroat-dx/2)/nP
        tfinal = 0.
        (xB,M2B,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
                
        if( eventIndex != xB.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating backwards from throat"))
    
        xIntegrate = np.concatenate((xB,xF))
        M2 = np.concatenate((M2B,M2F))
    
    return (xIntegrate, M2)

#==============================================================================
# Integrate subsonic, shock post-throat, supersonic flow through an axial 
# nozzle geometry
#==============================================================================
def integrateShock(nozzle,tol,params,xThroat,nPartitions):
    raise RuntimeError("Shock in nozzle case not implemented yet")
    return -1

#==============================================================================
# Integrate supersonic flow through an axial nozzle geometry
#==============================================================================
def integrateSupersonic(nozzle,tol,params,xThroat,nPartitions):
        
    # If nozzle converges only, assume choked flow at the exit    
    if( nozzle.wall.geometry.length - xThroat < 1e-12 ):
        M0 = 0.999 # start integration from this Mach number
        
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(M0**2,nozzle.wall.geometry.length)
        b.set_f_params(xThroat,nozzle.wall.geometry.length,nozzle.fluid.gam,         \
          nozzle.wall.geometry,params)
        dt = nozzle.wall.geometry.length/nPartitions
        tfinal = 0.
        (xIntegrate,M2,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
        
#        plt.plot(xIntegrate,M2)
#        plt.show()
#        sys.exit()
                
        if( eventIndex != xIntegrate.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating backwards from the exit"))

    # Else, assume nozzle is choked at throat, integrate forwards & backwards
    else:
        
        # Integrate forward from throat
        UpperM = 1.001 # start integration at this Mach number for aft portion
        dx = 1e-4 # 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001        
        f = scipy.integrate.ode(dM2dxForward,jac=None).set_integrator(       \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        f.set_initial_value(UpperM**2,xThroat+dx/2)
        f.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round((1-xThroat/nozzle.wall.geometry.length)*nPartitions)
        dt = (nozzle.wall.geometry.length-xThroat-dx/2)/nP
        tfinal = nozzle.wall.geometry.length
        (xF,M2F,eventIndex) = integrateODEwithEvents(f,dt,tfinal,1.,"f")
        
        if( eventIndex != xF.size ): # Relax UpperM and try again
            print 'Integration terminated early while integrating forwards from throat, trying again...'
            UpperM = 1.01
            dx = 1e-3
            f = scipy.integrate.ode(dM2dxForward,jac=None).set_integrator(       \
              'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
            f.set_initial_value(UpperM**2,xThroat+dx/2)
            f.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
            nP = np.round((1-xThroat/nozzle.wall.geometry.length)*nPartitions)
            dt = (nozzle.wall.geometry.length-xThroat-dx/2)/nP
            tfinal = nozzle.wall.geometry.length
            (xF,M2F,eventIndex) = integrateODEwithEvents(f,dt,tfinal,1.,"f")    
            if( eventIndex != xF.size ):        
                raise RuntimeError(("Integration terminated early while ",
                                   "integrating forwards from the throat"))
        
        # Integrate backward from throat
        LowerM = 0.999 # start integ. at this Mach number for fore portion
        dx = 1e-4 # 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(LowerM**2,xThroat-dx/2)
        b.set_f_params(xThroat,nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round(xThroat/nozzle.wall.geometry.length*nPartitions)
        dt = (xThroat-dx/2)/nP
        tfinal = 0.
        (xB,M2B,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
                
        if( eventIndex != xB.size ): # Relax LowerM and try again
            print 'Integration terminated early while integrating backwards from throat, trying again...'
            LowerM = 0.99
            dx = 1e-3
            b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
              'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
            b.set_initial_value(LowerM**2,xThroat-dx/2)
            b.set_f_params(xThroat,nozzle.fluid.gam,nozzle.wall.geometry,params)
            nP = np.round(xThroat/nozzle.wall.geometry.length*nPartitions)
            dt = (xThroat-dx/2)/nP
            tfinal = 0.
            (xB,M2B,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")            
            if( eventIndex != xB.size ): 
                raise RuntimeError(("Integration terminated early while ",
                                    "integrating backwards from throat"))
    
        xIntegrate = np.concatenate((xB,xF))
        M2 = np.concatenate((M2B,M2F))
    
    return (xIntegrate, M2)

#==============================================================================
# Cumulative trapezoidal integration using linear interpolation to return n 
# results for a vector of size n    
#==============================================================================
def integrateTrapezoidal(y,x):
        integrand = np.empty(x.size+1)
        dx = np.empty(x.size+1)
        
        dx0 = x[1] - x[0]
        dxN = x[-1] - x[-2]
        
        integrand[0] = y[0] - (y[1] - y[0])/dx0*(x[0] - (x[0] - dx0/2))
        integrand[1:-1] = y[1:] - (y[1:] - y[:-1])/(x[1:] - x[:-1])*         \
          (x[1:] - (x[1:] + x[:-1])/2)
        integrand[-1] = y[-1] + (y[-1] - y[-1])/dxN*((x[-1] + dxN/2) - x[-1])
        
        dx[0] = x[0] - dx0/2
        dx[1:-1] = (x[1:] + x[:-1])/2
        dx[-1] = x[-1] + dxN/2
        
        integral = scipy.integrate.cumtrapz(integrand,dx)
        
        return (integral)

#==============================================================================
# Perform quasi-1D area-averaged Navier-Stokes analysis of axisymmetric nozzle.
#% Solve for flow along length of non-ideal nozzle given geometry, inlet
#% stagnation temperature and pressure, and freestream temperature and
#% pressure. Iterate for Cf and stagnation temperature. An ODE for M^2 is 
#% solved given A, Cf, and Tstag. Pstag is found from mass conservation. T 
#% and P are found from def'n of stag. temp. Density rho is found from ideal 
#% gas law. 
#%
#% Returns M, density, pressure P, temperature T, stagnation 
#% temp. Tstag, stagnation pressure Pstag, velocity U, Re, internal heat
#% transfer coefficient hf, friction coefficient Cf, interior wall temp. Tw,
#% exterior wall temp. Text, and approximate stress along length of nozzle.
#==============================================================================
def Quasi1D_old(nozzle,output='verbose'):
    
    # Initialize
    gam = nozzle.fluid.gam
    xApparentThroat = nozzle.wall.geometry.findMinimumRadius()[0]
    # Using 1000 equally-spaced ODE integration steps will give thrust an
    # accuracy of 7 digits, to 2 decimal places. Integration will take about 
    # 6 seconds. Increasing the number of ODE integration steps to 10,000 will
    # yield an additional 1 or 2 decimal places of accuracy at nearly 10 times
    # the computational expense. Unfortunately, adaptive timesteps are not
    # implemented.
    nOdeIntegrationSteps = 1000;
    
    tol = {};
    tol["exitTempPercentError"]         = nozzle.tolerance.exitTempPercentError;
    tol["solverRelTol"]                 = nozzle.tolerance.solverRelTol;
    tol["solverAbsTol"]                 = nozzle.tolerance.solverAbsTol;
    tol["solverApparentThroatLocation"] = nozzle.tolerance.solverApparentThroatLocation;
    
    # Determine state of nozzle assuming ideal conditions
    pressureRatio = nozzle.inlet.Pstag/nozzle.environment.P
    TstagThroat = nozzle.inlet.Tstag
    PstagThroat = nozzle.inlet.Pstag
    TstagExit = nozzle.inlet.Tstag
    PstagExit = nozzle.inlet.Pstag    
    (status,shock) = nozzleState(nozzle,pressureRatio,PstagThroat,           \
      TstagThroat,PstagExit,TstagExit)
    
    # Initialize loop variables
    if hasattr(nozzle.wall,'temperature'):
        xPositionOld = np.linspace(0.,nozzle.wall.geometry.length,4000)
        Tstag = nozzle.wall.temperature.geometry.radius(xPositionOld)
        dTstagdx = nozzle.wall.temperature.geometry.radiusGradient(xPositionOld)
        Cf = np.array([0.004]*4000)
    else:
        Cf = np.array(([0.004, 0.004]))
        Tstag = np.array(([nozzle.inlet.Tstag, nozzle.inlet.Tstag-6.*nozzle.wall.geometry.length]))
        dTstagdx = np.array(([-6., -6.]))
        xPositionOld = np.array(([0., nozzle.wall.geometry.length]))
    
    if nozzle.thermalFlag == 1:
        maxIterations = 12 # max number of iterations to solve for Cf and Tstag
    else: # do not perform thermal analysis
        maxIterations = 1
    counter = 0 # used to count b/w number of iterations
    tolerance = tol["exitTempPercentError"] # tolerance for % error in
                # exit static temperature between iterations

    Texit_old = 0 # save previous static temperature
    
    #------------------------------ Begin Solver -----------------------------
    
    if output == 'verbose':
        introString = " Begin Solver ";
        nch = (60-len(introString))/2;
        sys.stdout.write('-' * nch);
        sys.stdout.write(introString);
        sys.stdout.write('-' * nch);
        sys.stdout.write('\n\n');
    
        sys.stdout.write(" Running non ideal nozzle computation (target error in exit temp %.3le): \n\n" % tolerance);
        sys.stdout.write('\t %s %s\n' % ("Iter".ljust(10), "Error %".ljust(10)));
        
    while( 1 ):
        
        counter += 1
                
        # Parameters passed to functions called by ODE 
        params = (xPositionOld,Cf,Tstag,dTstagdx)
        
        # Find where M = 1
        xApparentThroat = findApparentThroat(nozzle,tol,params)
        
        if( status == "no flow" ):
            raise UserWarning(("Prescribed inputs result in flow reversal ",
                               "in nozzle"))
        elif( status == "subsonic" ):
            (xPosition,M2) = integrateSubsonic(nozzle,tol,params,            \
              xApparentThroat,nOdeIntegrationSteps)
        elif( status == "shock" ):
            (xPosition,M2) = integrateShock(nozzle,tol,params,               \
              xApparentThroat,nOdeIntegrationSteps)
        else: # supersonic flow
            (xPosition,M2) = integrateSupersonic(nozzle,tol,params,          \
              xApparentThroat,nOdeIntegrationSteps)
              
        # Check output
        if any(M2<0.):
            print 'Impossible flow scenario: Square of Mach number < 0'
            print 'No flow in nozzle.'
            sys.exit(0);
            
        if( np.isnan(M2.any()) or np.isinf(M2.any()) ):
            raise RuntimeError("Unrealistic Mach number calculated")    

        # Calculate geometric properties
        D = nozzle.wall.geometry.diameter(xPosition)
        A = nozzle.wall.geometry.area(xPosition)
        
        # Calculate other 1D flow properties
        M = np.sqrt(M2)
        Tstag = np.interp(xPosition,xPositionOld,Tstag)
        T = Tstag/(1. + (gam-1.)*M2/2.) # static temp. from stag. temp. def.
        Pstag = nozzle.inlet.Pstag*(A[0]/A)*(areaMachFunc(gam,M[0])/         \
          areaMachFunc(gam,M))*np.sqrt(Tstag/Tstag[0]) # from mass conserv.
        P = Pstag/(1. + (gam-1.)*M2/2.)**(gam/(gam-1)) # from stag. press. def.
        density = P/(nozzle.fluid.R*T)
        U = M*np.sqrt(gam*nozzle.fluid.R*T) # velocity
        Re = density*U*D/dynamicViscosity(T) # Reynolds number from definition
                
        Cf = np.interp(xPosition,xPositionOld,Cf) # friction coefficient
            
        # Recalculate friction and heat
        # heat transfer coefficient to interior nozzle wall, estimated using
        # Chilton-Colburn analogy
        hf = nozzle.fluid.Pr(T)**(-2./3.)*density*nozzle.fluid.Cp(T)*U*Cf/2
        
        # Total thermal resistance (*dx) from fluid stag. temp. to ambient temp.
        # For a 2-layer wall
        #RwallPrime = np.log((D/2.+ti+to)/(D/2.+ti))/(2.*np.pi*ki) +          \
        #  np.log((D/2.+ti)/(D/2.))/(2.*np.pi*ko)
        #RtotPrime = 1./(hf*np.pi*D) + RwallPrime +                           \
        #  1./(nozzle.environment.hInf*np.pi*(D+2.*ti+2.*to))
        RwallPrime = np.zeros(len(D))
        tTempLower = np.zeros(len(D))
        for i in range(len(nozzle.wall.layer)):
            kTemp = nozzle.wall.layer[i].material.getThermalConductivity(3)
            tTempUpper = tTempLower + nozzle.wall.layer[i].thickness.radius(xPosition)
            RwallPrime = RwallPrime + np.log((D/2.+tTempUpper)/(D/2.+tTempLower))/(2.*np.pi*kTemp)
            tTempLower = tTempUpper
        RtotPrime = 1./(hf*np.pi*D) + RwallPrime +                           \
          1./(nozzle.environment.hInf*np.pi*(D+2.*tTempUpper))
        
        # Redefine stagnation temperature distribution (for axisymmetric nozzle)
#        if hasattr(nozzle.wall,'temperature'):                    
#            dTstagdx = np.interp(xPosition,xPositionOld,dTstagdx)
#        else:
        TstagXIntegrand = 1./(RtotPrime*density*U*A*nozzle.fluid.Cp(T))
        TstagXIntegral = integrateTrapezoidal(TstagXIntegrand,xPosition)
        Tstag = nozzle.environment.T*(1. - np.exp(-TstagXIntegral)) +        \
          nozzle.inlet.Tstag*np.exp(-TstagXIntegral)
        dTstagdx = (nozzle.environment.T - Tstag)/(RtotPrime*density*U*A*    \
          nozzle.fluid.Cp(T))
        
        # Redefine stagnation temperature distribution (for flat plate)
#        t = ti+to
#        TstagXIntegrand = 4./(nozzle.fluid.Cp(T)*density*U*D*(1./hf +        \
#          t/ki + 1./nozzle.environment.hInf))   
#        TstagXIntegral = integrateTrapezoidal(TstagXIntegrand,xPosition)
#        Tstag = nozzle.environment.T*(1. - np.exp(-TstagXIntegral)) +        \
#          nozzle.inlet.Tstag*np.exp(-TstagXIntegral)
#        dTstagdx = (nozzle.environment.T - Tstag)*4./(nozzle.fluid.Cp(T)*    \
#          density*U*D*(1./hf + t/ki +                    \
#          1./nozzle.environment.hInf))
          
        # Estimate interior wall temperature
        #Qw = nozzle.fluid.Cp(T)*density*U*D*dTstagdx/4.
        QwFlux = (Tstag - nozzle.environment.T)/RtotPrime/(np.pi*D) # W/m
                
        #Tinside = Tstag + Qw/hf # interior wall temperature
        if hasattr(nozzle.wall,'temperature'):
            Tinside = nozzle.wall.temperature.geometry.radius(xPosition)
            dTstagdx = (Tinside[1:]-Tinside[0:-1])/(xPosition[1:]-xPosition[0:-1])
            dTstagdx = np.hstack((dTstagdx,dTstagdx[-1]))
            Tstag = Tinside + QwFlux/hf
        else:
            Tinside = Tstag - QwFlux/hf
        #recoveryFactor = (Tinside/T - 1)/((gam-1)*M2/2)
        
        # Estimate exterior wall temperature
        #Toutside = nozzle.environment.T - Qw/nozzle.environment.hInf
        Toutside = nozzle.environment.T + QwFlux/nozzle.environment.hInf
    
        # Redefine friction coefficient distribution (Sommer & Short's method)
        TPrimeRatio = 1. + 0.035*M2 + 0.45*(Tinside/T - 1.)
        RePrimeRatio = 1./(TPrimeRatio*(TPrimeRatio)**1.5*(1. + 110.4/T)/    \
          (TPrimeRatio + 110.4/T))
        CfIncomp = 0.074/Re**0.2
        Cf = CfIncomp/TPrimeRatio/RePrimeRatio**0.2
        
        # Save old solution x position for next iteration
        xPositionOld = xPosition    
        
        # Recalculate nozzle status
        TstagThroat = np.interp(xApparentThroat,xPosition,Tstag)
        PstagThroat = np.interp(xApparentThroat,xPosition,Pstag)
        TstagExit = Tstag[-1]
        if( status == "shock" ):
            # approximate calculation
            shockIndex = 10000 
            pressureRatio = Pstag[shockIndex-1]/nozzle.environment.P
            PstagExit = Pstag[shockIndex-1] # temp. fix
            (status,shock) = nozzleState(nozzle,pressureRatio,PstagThroat,   \
              TstagThroat,PstagExit,TstagExit)
            PstagExit = Pstag[-1]
        else: # no shock in nozzle; assume no Pstag loss
            pressureRatio = Pstag[-1]/nozzle.environment.P
            PstagExit = Pstag[-1]
            (status,shock) = nozzleState(nozzle,pressureRatio,PstagThroat,   \
              TstagThroat,PstagExit,TstagExit)
        
        #print "%i\n" % counter
        
        # Check tolerance on static temperature at nozzle exit
        percentError = abs(T[-1] - Texit_old)/T[-1]
        #print "Percent error: %e\n" % (percentError*100)
        Texit_old = T[-1]
        
        if( counter >= maxIterations ):
            sys.stdout.write("\n WARNING: Done (max number of iterations (%i) reached)\n" % maxIterations);
            sys.stdout.write("Terminated with error in exit temp: %le\n\n" % percentError);
            #print "Iteration limit for quasi-1D heat xfer & friction reached\n"
            break
                
        if output == 'verbose':
            sys.stdout.write("\t %s %s\n" % (("%d" % counter).ljust(10), ("%.3le" % percentError).ljust(10)));
                
        if( percentError < tolerance ):
            if output == 'verbose':
                sys.stdout.write("\n Done (converged)\n\n");
            #print "%i iterations to converge quasi-1D heat xfer & friction \
#calcs\n" % counter    
            break
    
    # END OF while( ~converged )
    
    dAdx = nozzle.wall.geometry.areaGradient(xPosition)
    maxSlope = max(dAdx/np.pi/D)
    minSlope = min(dAdx/np.pi/D)
    
    # Calculate mass flow rate
    mdot = massFlowRate(nozzle.fluid,Pstag,A,Tstag,M)
    
    # Calculate thrust
    exitAngle = np.arctan2(dAdx[-1],np.pi*D[-1])
    divergenceFactor = (1. + np.cos(exitAngle))/2.
    netThrust = divergenceFactor*mdot[0]*(U[-1] - nozzle.mission.mach*       \
      nozzle.environment.c) + (P[-1] - nozzle.environment.P)*A[-1]
    
    if nozzle.structuralFlag == 1:
#        # Simplified stress calculation (calculate stresses IN LOAD LAYER ONLY)
#        # Assumptions: nozzle is a cylinder; nozzle length is not constrained;
#        #              thermal mismatch at interface of both material layers is 
#        #              neglected; steady state
#    
#        # Determine hoop stress for outermost layer only
#        stressHoop = P*(D/2.+tTempUpper)/nozzle.wall.layer[-1].thickness.radius(xPosition)
#        
#        # Determine thermal stress for outermost layer only
#        ri = D/2. + tTempUpper - nozzle.wall.layer[-1].thickness.radius(xPosition) # inner radius
#        ro = D/2. + tTempUpper # outer radius
#        
#        E1 = nozzle.wall.layer[-1].material.getElasticModulus(1)
#        E2 = nozzle.wall.layer[-1].material.getElasticModulus(2)
#        alpha1 = nozzle.wall.layer[-1].material.getThermalExpansionCoef(1)
#        alpha2 = nozzle.wall.layer[-1].material.getThermalExpansionCoef(2)
#        v = nozzle.wall.layer[-1].material.getPoissonRatio()
#        
#        stressThermalRadial = E1*alpha1*(Tinside-Toutside)/            \
#          (2.*(1.-v))*(1./np.log(ro/ri))*(1. - 2.*ri**2./                        \
#          (ro**2. - ri**2.)*np.log(ro/ri))
#        stressThermalAxial = E2*alpha2*(Tinside-Toutside)/               \
#          (2.*(1.-v))*(1./np.log(ro/ri))*(1. - 2.*ri**2./                        \
#          (ro**2. - ri**2.)*np.log(ro/ri))      
#        
#        # THE EQUATIONS BELOW NEED TO BE CHECKED (ITEMS HAVE BEEN RENAMED)
#        # Estimate vonMises, even though not really valid for composites
#    #    stressVonMises = np.sqrt( (stressHoop+stressThermalAxial)**2 -           \
#    #       stressHoop*stressThermalRadial + stressThermalRadial**2 )    
#        stressMaxPrincipal = stressHoop + stressThermalAxial
#        stressPrincipal = (stressMaxPrincipal, stressThermalRadial,              \
#          np.zeros(xPosition.size))
#          
#        # Calculate cycles to failure, Nf
#        Nf = nozzlemod.lifetime.estimate(Tinside,stressMaxPrincipal,1)
    
        # --- Run AEROS
        nozzle.wallResults = np.transpose(np.array([xPosition,Tinside,P]))
        nozzle.runAEROS = 0;
        if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1:
            
            #try:
            #    runAEROS(nozzle, output);
            #    AEROSPostProcess(nozzle, output);
            #except:
            #    sys.stdout.write("  ## WARNING : CALL TO AEROS IGNORED.\n");
            runAEROS(nozzle, output);
            AEROSPostProcess(nozzle, output);

    else: # do not perform structural analysis
        pass;
    
    # Calculate volume of nozzle material (approximately using trap. integ.)
    volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle)    
    
    # Assign all data for output
    flowTuple = (M, U, density, P, Pstag, T, Tstag, Re)
    heatTuple = (Tinside, Toutside, Cf, hf, QwFlux)
    geoTuple = (D, A, dAdx, minSlope, maxSlope)
    performanceTuple= (volume, mass, netThrust, mdot, Pstag[-1]/Pstag[0], \
      Tstag[-1]/Tstag[0], status)
    
    return (xPosition, flowTuple, heatTuple, geoTuple, performanceTuple)
    
# END OF analysis(nozzle,tol)


#==============================================================================
# Perform quasi-1D area-averaged Navier-Stokes analysis of axisymmetric nozzle.
#% Solve for flow along length of non-ideal nozzle given geometry, inlet
#% stagnation temperature and pressure, and freestream temperature and
#% pressure. Iterate for Cf and stagnation temperature. An ODE for M^2 is 
#% solved given A, Cf, and Tstag. Pstag is found from mass conservation. T 
#% and P are found from def'n of stag. temp. Density rho is found from ideal 
#% gas law. 
#%
#% Returns M, density, pressure P, temperature T, stagnation 
#% temp. Tstag, stagnation pressure Pstag, velocity U, Re, internal heat
#% transfer coefficient hf, friction coefficient Cf, interior wall temp. Tw,
#% exterior wall temp. Text.
#==============================================================================
def Quasi1D(nozzle,output='verbose'):
    
    # Initialize inputs for nozzle analysis

    # Discretization
    nbreaks = 1000 # save data for nbreaks along length of nozzle
    nApproxPoints = 8000 # number of points to use for piecewise-linear
                         # approximation of non-piecewise-linear geometries

    # Nozzle interior wall geometry
    if nozzle.wall.geometry.type == 'piecewise-linear':
        xgeo = list(nozzle.wall.geometry.nodes[:,0].transpose())
        rgeo = list(nozzle.wall.geometry.nodes[:,1].transpose())
    else:
        xgeo = list(np.linspace(nozzle.wall.geometry.xstart,
            nozzle.wall.geometry.xend, nApproxPoints))
        rgeo = list(nozzle.wall.geometry.radius(xgeo))

    # Interior wall temperature (if necessary)
    if hasattr(nozzle.wall, 'temperature'):
        if nozzle.wall.temperature.geometry.type == 'piecewise-linear':
            xwalltemp = list(nozzle.wall.temperature.geometry.nodes[:,0].transpose())
            walltemp = list(nozzle.wall.temperature.geometry.nodes[:,1].transpose())
        else:
            xwalltemp = list(np.linspace(nozzle.wall.temperature.geometry.xstart,
                nozzle.wall.temperature.geometry.xend, nApproxPoints))
            walltemp = list(nozzle.wall.temperature.geometry.radius(xwalltemp))
    else:
        xwalltemp = []
        walltemp = []

    # Thermal layer
    if nozzle.wall.layer[0].thickness.type == 'piecewise-linear':
        xlayer1 = list(nozzle.wall.layer[0].thickness.nodes[:,0].transpose())
        tlayer1 = list(nozzle.wall.layer[0].thickness.nodes[:,1].transpose())
    else:
        xlayer1 = list(np.linspace(nozzle.wall.layer[0].thickness.xstart,
            nozzle.wall.layer[0].thickness.xend, nApproxPoints))
        tlayer1 = list(nozzle.wall.layer[0].thickness.radius(xlayer1))
    k1 = nozzle.wall.layer[0].material.getThermalConductivity(direction=3)

    # Air gap
    if nozzle.wall.layer[1].thickness.type == 'piecewise-linear':
        xlayer2 = list(nozzle.wall.layer[1].thickness.nodes[:,0].transpose())
        tlayer2 = list(nozzle.wall.layer[1].thickness.nodes[:,1].transpose())
    else:
        xlayer2 = list(np.linspace(nozzle.wall.layer[1].thickness.xstart,
            nozzle.wall.layer[1].thickness.xend, nApproxPoints))
        tlayer2 = list(nozzle.wall.layer[1].thickness.radius(xlayer2))
    k2 = nozzle.wall.layer[1].material.getThermalConductivity(direction=3)

    # Inner load layer
    if nozzle.wall.layer[2].thickness.type == 'piecewise-linear':
        xlayer3 = list(nozzle.wall.layer[2].thickness.nodes[:,0].transpose())
        tlayer3 = list(nozzle.wall.layer[2].thickness.nodes[:,1].transpose())
    else:
        xlayer3 = list(np.linspace(nozzle.wall.layer[2].thickness.xstart,
            nozzle.wall.layer[2].thickness.xend, nApproxPoints))
        tlayer3 = list(nozzle.wall.layer[2].thickness.radius(xlayer3))
    k3 = nozzle.wall.layer[2].material.getThermalConductivity(direction=3)

    # Middle load layer
    if nozzle.wall.layer[3].thickness.type == 'piecewise-linear':
        xlayer4 = list(nozzle.wall.layer[3].thickness.nodes[:,0].transpose())
        tlayer4 = list(nozzle.wall.layer[3].thickness.nodes[:,1].transpose())
    else:
        xlayer4 = list(np.linspace(nozzle.wall.layer[3].thickness.xstart,
            nozzle.wall.layer[3].thickness.xend, nApproxPoints))
        tlayer4 = list(nozzle.wall.layer[3].thickness.radius(xlayer4))
    k4 = nozzle.wall.layer[3].material.getThermalConductivity(direction=3)
    
    # Outer load layer
    if nozzle.wall.layer[4].thickness.type == 'piecewise-linear':
        xlayer5 = list(nozzle.wall.layer[4].thickness.nodes[:,0].transpose())
        tlayer5 = list(nozzle.wall.layer[4].thickness.nodes[:,1].transpose())
    else:
        xlayer5 = list(np.linspace(nozzle.wall.layer[4].thickness.xstart,
            nozzle.wall.layer[4].thickness.xend, nApproxPoints))
        tlayer5 = list(nozzle.wall.layer[4].thickness.radius(xlayer5))
    k5 = nozzle.wall.layer[4].material.getThermalConductivity(direction=3)
   
    # Inlet parameters
    tsi = nozzle.inlet.Tstag
    psi = nozzle.inlet.Pstag

    # Initial parameters
    dtsi = -6. # stagnation temp gradient guess
    cfi = 0.004 # friction coefficient guess

    # Fluid properties
    g = nozzle.fluid.gam
    gasconstant = nozzle.fluid.R

    # Mission properties
    missionmach = nozzle.mission.mach

    # Environment properties
    hinf = nozzle.environment.hInf
    tenv = nozzle.environment.T
    penv = nozzle.environment.P
    cenv = nozzle.environment.c

    # Gauss-Seidel fluid-thermal iteration properties
    eps1 = nozzle.tolerance.exitTempPercentError # relative tolerance for exit temp b/w iterations
    if nozzle.thermalFlag == 1:
        maxiter = 12 # max number of iterations for conjugate heat xfer
    else: # do not perform thermal analysis
        maxiter = 1

    # ODE integration properties
    maxstep = 100000 # max number of steps
    eps2 = nozzle.tolerance.solverAbsTol # absolute allowable error of RK4 integration
    ns = 1001 # approximate number of steps to save (governs closeness of integrated mach values)
    himag = 1e-4 # guess for starting step
    hminmag = 1e-12 # smallest allowable step
    hmaxmag = 5e-3 # largest allowable step
    singularitydy = 1e-3 # increment above/below M=1 to start integration

    # Outputs
    x = [] # x-coordinate along nozzle axis
    temp = [] # temperature
    p = [] # static pressure
    rho = [] # density
    u = [] # velocity
    mach = [] # mach
    tempinside = [] # inside wall temperature
    tempoutside = [] # outside wall temperature
    netthrust = [] # net thrust

    # Run thermo-fluid analysis
    quasi1dnozzle.analyze(xgeo, rgeo, nbreaks, xwalltemp, walltemp, 
        xlayer1, tlayer1, k1, xlayer2, 
        tlayer2, k2, xlayer3, tlayer3, k3, xlayer4, tlayer4, k4, xlayer5,
        tlayer5, k5, tsi, dtsi, psi, cfi, missionmach, g, gasconstant, hinf,
        tenv, cenv, penv, eps1, maxiter, maxstep, eps2, ns, himag, hminmag, hmaxmag, 
        singularitydy,
        x, temp, p, rho, u, mach, tempinside, tempoutside, netthrust)

    # Convert data to Numpy array
    x = np.array(x)
    temp = np.array(temp)
    p = np.array(p)
    rho = np.array(rho)
    u = np.array(u)
    mach = np.array(mach)
    tempinside = np.array(tempinside)
    tempoutside = np.array(tempoutside)
    netthrust = netthrust[0]

    # plt.plot(x,mach,'.')
    # plt.show()

    ps = p*(1+(g-1)*mach**2/2.)**(g/(g-1.)) # stagnation pressure
    
    # Run secondary thermal analysis and structural analysis if necessary
    nozzle.wallResults = np.transpose(np.array([x, tempinside, ps]))
    nozzle.runAEROS = 0
    if nozzle.thermalFlag == 1 or nozzle.structuralFlag == 1 or \
        'MASS' in nozzle.responses or 'MASS_WALL_ONLY' in nozzle.responses:
            runAEROS(nozzle, output)
            AEROSPostProcess(nozzle, output)

    return netthrust, x, tempinside, ps, p, u


    
def Run(nozzle, **kwargs):
        
    output = 'verbose';
    writeToFile = 1;
    skipAero = 0;
    
    if 'output' in kwargs:
        output = kwargs['output'];
        
    if 'writeToFile' in kwargs:
        writeToFile = int(kwargs['writeToFile']);

    if 'skipAero' in kwargs and kwargs['skipAero'] == 1:
        skipAero = 1;
    
    # Obtain mass and volume
    # if 'MASS' in nozzle.responses or 'VOLUME' in nozzle.responses:
    #     volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle)
    #     if 'MASS' in nozzle.responses:
    #         nozzle.responses['MASS'] = np.sum(mass)
    #     if 'VOLUME' in nozzle.responses:
    #         nozzle.responses['VOLUME'] = np.sum(volume)
        
    # Calculate mass gradients if necessary
    # if 'MASS' in nozzle.gradients and nozzle.gradients['MASS'] is not None:
    #     if ( nozzle.gradientsMethod == 'ADJOINT' ):
    #         # Convergence study using B-spline coefs show finite difference mass gradients
    #         # converge. Most accurate gradients use absolute step size 1e-8. RWF 5/10/17
    #         nozzle.gradients['MASS'] = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-8);
    #     elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
    #         nozzle.gradients['MASS'] = nozzlemod.geometry.calcMassGradientsFD(\
    #           nozzle,nozzle.fd_step_size);
    #     else:
	# 	    sys.stderr.write('  ## ERROR : Unknown gradients computation '
	# 	      'method.\n');
	# 	    sys.exit(1);
	
	# Calculate volume gradients if necessary
    # if 'VOLUME' in nozzle.gradients and nozzle.gradients['VOLUME'] is not None:
    #     sys.stderr.write('\n ## ERROR : gradients for VOLUME are not supported\n\n');
    #     sys.exit(1);
 
    # # Obtain mass of wall and gradients only if requested
    # if 'MASS_WALL_ONLY' in nozzle.responses:
    #     volume, mass = nozzlemod.geometry.calcVolumeAndMass(nozzle)
    #     n_layers = len(nozzle.wall.layer);
    #     nozzle.responses['MASS_WALL_ONLY'] = np.sum(mass[:n_layers]);

    # if 'MASS_WALL_ONLY' in nozzle.gradients and nozzle.gradients['MASS_WALL_ONLY'] is not None:
    #     if ( nozzle.gradientsMethod == 'ADJOINT' ):
    #         # Convergence study using B-spline coefs show finite difference mass gradients
    #         # converge. Most accurate gradients use absolute step size 1e-8. RWF 5/10/17
    #         nozzle.gradients['MASS_WALL_ONLY'] = nozzlemod.geometry.calcMassGradientsFD(nozzle,1e-8,components='wall-only');
    #     elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
    #         nozzle.gradients['MASS_WALL_ONLY'] = nozzlemod.geometry.calcMassGradientsFD(\
    #           nozzle,nozzle.fd_step_size,components='wall-only');
    #     else:
	# 	    sys.stderr.write('  ## ERROR : Unknown gradients computation '
	# 	      'method.\n');
	# 	    sys.exit(1);

    # Run aero-thermal-structural analysis if necessary
    runAeroThermalStructuralProblem = 0;
    for k in nozzle.responses:
        if k not in ['MASS','VOLUME','MASS_WALL_ONLY']:
            runAeroThermalStructuralProblem = 1;  
            break;  
    
    # Run aero-thermal-structural gradient analysis if necessary
    if nozzle.gradientsFlag == 1:
        runAeroThermalStructuralGradients = 0;        
        for k in nozzle.gradients:
            #if k not in ['MASS','VOLUME','MASS_WALL_ONLY']:
            if nozzle.gradients[k] is not None:
                runAeroThermalStructuralGradients = 1;  
                break;
   
    if runAeroThermalStructuralProblem:
        
        # # Run aero analysis (+ 1D thermal analysis if requested)
        # xPosition, flowTuple, heatTuple, geoTuple, \
        #     performanceTuple = Quasi1D_old(nozzle,output);

        # # Assign function values        
        # if 'THRUST' in nozzle.responses:
        #     nozzle.responses['THRUST'] = performanceTuple[2]  

        # if 'WALL_TEMPERATURE' in nozzle.responses:
        #     nozzle.responses['WALL_TEMPERATURE'] = np.interp(\
        #         nozzle.outputLocations['WALL_TEMPERATURE'], xPosition, heatTuple[0])
        
        # if 'WALL_PRESSURE' in nozzle.responses:
        #     nozzle.responses['WALL_PRESSURE'] = np.interp(\
        #         nozzle.outputLocations['WALL_PRESSURE'], xPosition, flowTuple[3])

        # if 'PRESSURE' in nozzle.responses:
        #     nozzle.responses['PRESSURE'] = np.interp(\
        #         nozzle.outputLocations['PRESSURE'][:,0], xPosition, flowTuple[3])

        # if 'VELOCITY' in nozzle.responses:
        #     nr, nc = nozzle.outputLocations['VELOCITY'].shape
        #     nozzle.responses['VELOCITY'] = np.zeros((nr,3))
        #     nozzle.responses['VELOCITY'][:,0] = np.interp(\
        #         nozzle.outputLocations['VELOCITY'][:,0], xPosition, flowTuple[1])

        # Run aero analysis, and thermal and structural analyses if necessary
        if skipAero:
            sys.stdout.write("WARNING: Skipping low-fidelity aero analysis.\n")
        else:
            thrust, x, walltemp, wallpress, press, velocity = Quasi1D(nozzle, output)
        
        # Assign function values        
        if 'THRUST' in nozzle.responses:
            nozzle.responses['THRUST'] = thrust
            
        if 'WALL_TEMPERATURE' in nozzle.responses:
            nozzle.responses['WALL_TEMPERATURE'] = np.interp(\
                nozzle.outputLocations['WALL_TEMPERATURE'], x, walltemp)

        if 'WALL_PRESSURE' in nozzle.responses:
            nozzle.responses['WALL_PRESSURE'] = np.interp(\
                nozzle.outputLocations['WALL_PRESSURE'], x, wallpress)

        if 'PRESSURE' in nozzle.responses:
            nozzle.responses['PRESSURE'] = np.interp(\
                nozzle.outputLocations['PRESSURE'][:,0], x, press)

        if 'VELOCITY' in nozzle.responses:
            nr, nc = nozzle.outputLocations['VELOCITY'].shape
            nozzle.responses['VELOCITY'] = np.zeros((nr,3))
            nozzle.responses['VELOCITY'][:,0] = np.interp(\
                nozzle.outputLocations['VELOCITY'][:,0], x, velocity)

    # Obtain mass (volume is currently not accepted as a nozzle response)
    if 'MASS' in nozzle.responses or 'MASS_WALL_ONLY' in nozzle.responses:
        total_mass, wall_mass = getMass(nozzle, output)
        if 'MASS' in nozzle.responses:
            nozzle.responses['MASS'] = total_mass;
        if 'MASS_WALL_ONLY' in nozzle.responses:
            nozzle.responses['MASS_WALL_ONLY'] = wall_mass;

    # Calculate gradients if necessary
    if nozzle.gradientsFlag == 1 and runAeroThermalStructuralGradients:
    
        if ( output == 'verbose' ):
            sys.stdout.write('Running gradient analysis\n');
    
        if ( nozzle.gradientsMethod == 'ADJOINT' ):
            sys.stderr.write('\n ## ERROR : Adjoint gradients for low-fidelity '
                'thrust calculation are not available.\n');
            sys.exit(1);
        elif ( nozzle.gradientsMethod == 'FINITE_DIFF' ):
            multif.gradients.calcGradientsFD(nozzle,nozzle.fd_step_size,output);           
        else:
            sys.stderr.write('  ## ERROR : Unknown gradients computation '
                'method.\n');
            sys.exit(1);
            
        # Write separate gradients file
        # gradFile = open(nozzle.gradientsFile,'w');
        # for k in nozzle.outputTags:
        #     np.savetxt(gradFile,nozzle.gradients[k]);
        # gradFile.close();
    
    # Write data
    if writeToFile:
        if nozzle.outputFormat == 'PLAIN':
            nozzle.WriteOutputFunctions_Plain();
        else:
            nozzle.WriteOutputFunctions_Dakota();
        
    return 0;
