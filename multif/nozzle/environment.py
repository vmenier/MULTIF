# -*- coding: utf-8 -*-
"""
Environment module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

import numpy as np

class Environment:
    def __init__(self, altitude, hInf):
        
        # Altitude measured in FEET, all other parameters in SI units
        self.reset(altitude)
        self.setHeatTransferCoefficient(hInf)
    
    # Set heat transfer coefficient from external surface of nozzle to
    # surrounding structure and eventually atmosphere
    def setHeatTransferCoefficient(self,hInf):
        self.hInf = hInf # W/m^2/K, hea transfer coef. to atmosphere
    
    # Manually set pressure
    def setPressure(self,P):
        self.P = P # Pa
        self.consistentWithStandardAtmosphere = 0
    
    # Manually set temperature
    def setTemperature(self,T):
        self.T = T # K
        self.consistentWithStandardAtmosphere = 0
        
    def reset(self,altitude):
        self.altitude = altitude
        (T, rho, P, mu, c) = StndAtm(altitude)
        self.P = P/0.020885 # Pa, atmospheric pressure
        self.T = T/1.8 # K, atmospheric temperature
        self.rho = rho/0.068521/0.028317 # kg/m^3, atmospheric density
        self.mu = mu/0.22481/0.092903 # N-s/m^2, dynamic viscosity
        self.c = c/3.2808 # m/s, speed of sound
        
        # record whether environment properties are consistent with the US
        # 1976 Standard Atmosphere
        self.consistentWithStandardAtmosphere = 1
				
				
				
				
#==============================================================================
# % StndAtm calculates properties of the 1976 U.S. Standard Atmosphere for
# % up to altitudes of 230,000 ft. This function is adapted from Ilan Kroo's
# % 1995 JavaScript code found on Stanford's AA 241 webpage. 
# %
# % StndAtm(H) returns a structure with atmospheric properties corresponding
# % to altitude H in English units. H can be a vector or scalar.
# %
# % Property                  Units (US)         Units (SI)
# % ------------------------------------------------------------------------
# % Altitude                  ft                 m
# % Velocity                  ft/sec             m/s
# % Characteristic Length     ft                 m
# % Temperature               degR               degK
# % Density                   sl/ft^3            kg/m^3
# % Pressure                  lb/ft^2            N/m^2
# % Kinematic Viscosity       ft/sec             m/s
# % Dynamic Viscosity         lb-sec/ft^2        N-s/m^2
# % Dynamic Pressure          lb/ft^2            N/m^2
# % Speed of Sound            ft/sec             m/s
# % Mach Number               -                  -
# % Reynolds Number           -                  -
# % Turbulent Friction Coeff. -                  -
# % Laminar Friction Coeff.   -                  -
# % Critical Pressure Coeff.  -                  -
# % Minimum Pressure Coeff.   -                  -
# %
# % Rick Fenrich 9/12/14
#============================================================================== 
def StndAtm(alt):
    
    T = 0. # temperature
    rho =  0. # density
    P  = 0. # pressure
    c  = 0. # speed of sound
    mu =  0. # dynamic viscosity
    
    # ================== CALCULATE ATMOSPHERE PROPERTIES =====================
       
    TEMPSL = 518.67 # degR
    RHOSL = 0.00237689 # slug/ft^3
    PRESSSL = 2116.22 # lbf/ft^2
    saTheta = 1.0
    saSigma = 1.0
    saDelta = 1.0
    
    if( alt > 232940 ):
        print "Atmospheric model only valid up to 232,294 ft"
        saTheta = 0.
        saSigma = 0.
        saDelta = 0.
    if( alt< 232940 ):
        saTheta = 1.434843 - alt /337634
        saSigma = (0.79899-alt /606330)**11.20114
        saDelta = (0.838263-alt /577922)**12.20114
    if( alt < 167323 ):
        saTheta = 0.939268
        saSigma = 0.00116533*np.exp( (alt - 154200)/-25992 )
        saDelta = 0.00109456*np.exp( (alt - 154200)/-25992 )
    if( alt < 154199 ):
        saTheta = 0.482561 + alt /337634
        saSigma = (0.857003+alt /190115)**-13.20114
        saDelta = (0.898309+alt /181373)**-12.20114
    if( alt < 104987 ):
        saTheta = 0.682457 + alt /945374
        saSigma = (0.978261+alt /659515)**-35.16319
        saDelta = (0.988626+alt /652600)**-34.16319
    if( alt < 65617 ):
        saTheta = 0.751865
        saSigma = 0.297076*np.exp((36089- alt)/20806 )
        saDelta = 0.223361*np.exp((36089- alt)/20806 )
    if( alt < 36089 ):
        saTheta = 1.0 - alt /145442
        saSigma = (1.0-alt /145442)**4.255876
        saDelta = (1.0-alt /145442)**5.255876
    
    T = TEMPSL*saTheta # degR
    rho = RHOSL*saSigma # sl/ft^3
    P = PRESSSL*saDelta # lb/ft^2
    mu = 0.0226968*T**1.5 / (T+198.72) / 1000000.0 # lb-sec/ft^2
    c = np.sqrt( 1.4*1716.56*T ) # ft/sec
		
    
    return (T, rho, P, mu, c)
