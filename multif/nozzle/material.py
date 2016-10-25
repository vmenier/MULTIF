# -*- coding: utf-8 -*-
"""
Material module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

import component

class Material:
    def __init__(self, name, matType, structure):
        if matType == 'ISOTROPIC' or matType == 'ANISOTROPIC_SHELL':
            self.name = name
            self.type = matType
        elif matType == 'FIXED_RATIO_PANEL':
            self.name = name
            self.type = matType
        else:
            raise NotImplementedError('Only ISOTROPIC or ANISOTROPIC_SHELL ' \
                  'materials are implemented')
                
        if structure == 'single':
            self.panel = 0 # material properties DO NOT depend on thickness
        elif structure == 'panel':
            self.panel = 1 # material properties DO depend on thickness
        else:
            raise NotImplementedError('Only "single" materials or "panel" '  \
                  'material structures are implemented')
                  
    def setPanelLayers(self,layerNames,layerRatios,layerMaterials):
        if self.panel == 0:
            raise RuntimeError('Panel layers cannot be set for a material '  \
              'with a non-panel structure.')
              
        if abs(sum(layerRatios) - 1) > 1e-12:
            raise RuntimeError('Panel layer ratios must add up to 1. Each '  \
              'layer ratio determines the fraction of the total thickness '  \
              'of the panel used by that material')
        
        self.layer = list()
        for i in range(len(layerNames)):
            self.layer.append(component.Wall(layerNames[i]));
            self.layer[i].ratio = layerRatios[i]
            self.layer[i].material = layerMaterials[i]
            
        if ( len(self.layer) != 3 or self.layer[0].ratio != self.layer[2].ratio
          or self.layer[0].material.name != self.layer[2].material.name ):
              raise RuntimeError('Currently FIXED_RATIO_PANEL only works for' \
                ' symmetric panels with 3 layers.')
            
    
    def setDensity(self,rho):
        if self.panel == 0: # single material
            try:
                self.rho = float(rho)
            except:
                raise RuntimeError('Single material must only have 1 value ' \
                      'provided for density')
        else: # panel material
            try:
                self.rhoSurface = float(rho[0])
                self.rhoMiddle = float(rho[1])
            except:
                raise RuntimeError('Panel structure material must have 2 '   \
                      'values provided for density')
        return 0

    def setElasticModulus(self,E):
        if self.type == 'ISOTROPIC':
            if hasattr(self,'G') and hasattr(self,'v'):
                raise RuntimeError('Specifying elastic modulus E when shear' \
                      ' modulus G and poisson ratio v have already been '    \
                      'specified violates the isotropic relation: G=E/'      \
                      '(2*(1+v)). Only specify 2 of the 3.')            
            try:
                self.E = float(E)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for elastic modulus')
        elif self.type == 'ANISOTROPIC_SHELL' and self.panel == 0:
            if isinstance(E,float): # 1 value provided
                self.E1 = E
                self.E2 = E        
            else:
                try:
                    self.E1 = float(E[0])
                    self.E2 = float(E[1])
                except:
                    raise RuntimeError('Anisotropic material must have 1-2 ' \
                          'values provided for elastic modulus')
        else: # anisotropic shell, panel material
            try:
                self.E1Surface = float(E[0])
                self.E2Surface = float(E[1])
                self.E1Middle = float(E[2])
                self.E2Middle = float(E[3])
            except:
                raise RuntimeError('Anisotropic panel material must have 4'  \
                      'values provided for elastic modulus')
        return 0
        
    def setShearModulus(self,G):
        if self.type == 'ISOTROPIC':
            if hasattr(self,'E') and hasattr(self,'v'):
                raise RuntimeError('Specifying shear modulus G when elastic' \
                      ' modulus E and poisson ratio v have already been '    \
                      'specified violates the isotropic relation: G=E/'      \
                      '(2*(1+v)). Only specify 2 of the 3.')
            try:
                self.G12 = float(G)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for shear modulus')
        elif self.type == 'ANISOTROPIC_SHELL' and self.panel == 0:
            if isinstance(G,float): # 1 value provided
                self.G12 = G      
            else:
                raise RuntimeError('Anisotropic shell material must have 1 ' \
                      'value provided for shear modulus')
        else: # anisotropic shell, panel material
            try:
                self.G12Surface = float(G[0])
                self.G12Middle = float(G[1])
            except:
                raise RuntimeError('Anisotropic shell panel material must '  \
                      'have 2 values provided for shear modulus')
        return 0        

    def setPoissonRatio(self,v):
        if self.type == 'ISOTROPIC':
            if hasattr(self,'E') and hasattr(self,'G'):
                raise RuntimeError('Specifying poisson ratio v when elastic' \
                      ' modulus E and shear modulus G have already been '    \
                      'specified violates the isotropic relation: G=E/'      \
                      '(2*(1+v)). Only specify 2 of the 3.')     
            try:
                self.v = float(v)
            except:
                raise NotImplementedError('Isotropic material must have 1 '  \
                      'value specified for poisson ratio')
        elif self.type == 'ANISOTROPIC_SHELL' and self.panel == 0:
            try:
                self.v12 = float(v)
            except:
                raise NotImplementedError('Anisotropic shell material must ' \
                      'have 1 value specified for poisson ratio')
        else: # anisotropic shell, panel material
            try:
                self.v12Surface = float(v[0])
                self.v12Middle = float(v[1])
            except:
                raise RuntimeError('Anisotropic shell panel material must '  \
                      'have 2 values provided for poisson ratio')                            
        return 0
        
    def setMutualInfluenceCoefs(self,nu):
        if self.type == 'ISOTROPIC':
            raise NotImplementedError('Isotropic material does not have '    \
                  'coefficients of mutual influence')
        elif self.type == 'ANISOTROPIC_SHELL' and self.panel == 0:
            if isinstance(nu,float): # 1 value provided
                self.nu1 = nu
                self.nu2 = nu       
            else:
                try:
                    self.nu1 = float(nu[0])
                    self.nu2 = float(nu[1])
                except:
                    raise RuntimeError('Anisotropic shell material must '    \
                          'have 1-2 values provided for coefficients of '    \
                          'mutual influence')
        else: # anisotropic shell, panel material
            try:
                self.nu1Surface = float(nu[0])
                self.nu2Surface = float(nu[1])
                self.nu1Middle = float(nu[2])
                self.nu2Middle = float(nu[3])
            except:
                raise RuntimeError('Anisotropic panel material must have 4'  \
                      'values provided for coefficients of mutual influence')                           
        return 0        
       
    def setThermalConductivity(self,k):
        if self.type == 'ISOTROPIC':
            try:
                self.k = float(k)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for thermal conductivity')
        elif self.type == 'ANISOTROPIC_SHELL' and self.panel == 0:
            if isinstance(k,float): # 1 value provided
                self.k1 = k
                self.k2 = k
                self.k3 = k
            elif len(k) == 2:
                self.k1 = float(k[0])
                self.k2 = float(k[1])
                self.k3 = (float(k[0]) + float(k[1]))/2
            elif len(k) == 3:
                self.k1 = float(k[0])
                self.k2 = float(k[1])
                self.k3 = float(k[2])
            else:
                raise RuntimeError('1 to 3 values must be input for the '    \
                      'anisotropic-shell material thermal conductivity')
        else: # anisotropic_shell, panel material
            try:
                self.k1Surface = float(k[0])
                self.k2Surface = float(k[1])
                self.k3Surface = float(k[2])
                self.k1Middle = float(k[3])
                self.k2Middle = float(k[4])
                self.k3Middle = float(k[5])
            except:
                raise RuntimeError('Anisotropic panel material must have 6 ' \
                      'values provided for thermal conductivity')
        return 0

    def setThermalExpansionCoef(self,alpha):
        if self.type == 'ISOTROPIC':
            try:
                self.alpha = float(alpha)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for thermal expansion coef')
        elif self.type == 'ANISOTROPIC_SHELL' and self.panel == 0:
            if isinstance(alpha,float): # 1 value provided
                self.alpha1 = alpha
                self.alpha2 = alpha
                self.alpha12 = alpha
            else:
                try:
                    self.alpha1 = float(alpha[0])
                    self.alpha2 = float(alpha[1])
                    self.alpha12 = float(alpha[2])
                except:
                    raise RuntimeError('Anisotropic material must have 3 '   \
                          'values provided for thermal expansion coef')
        else: # anisotropic, panel material
            try:
                self.alpha1Surface = float(alpha[0])
                self.alpha2Surface = float(alpha[1])
                self.alpha12Surface = float(alpha[2])
                self.alpha1Middle = float(alpha[3])
                self.alpha2Middle = float(alpha[4])
                self.alpha12Middle = float(alpha[5])
            except:
                raise RuntimeError('Anisotropic panel material must have 6'  \
                      'values provided for thermal expansion coef')
        return 0
     
    # Calculate (if necessary) and return density
    def getDensity(self,*args):
        if self.type == 'ISOTROPIC' or self.type == 'ANISOTROPIC_SHELL': # single material
            return self.rho
        elif self.type == 'FIXED_RATIO_PANEL':
            rho = 0.
            for i in range(len(self.layer)):
                rho = rho + self.layer[i].material.getDensity()*self.layer[i].ratio
            self.rho = rho
            return rho
        else: # panel material that does not have fixed ratio
            if len(args) == 3:
                x = args[0]
                tSurface = args[1]
                tMiddle = args[2]
                rho = (self.RhoSurface*2.*tSurface.radius(x) +               \
                      self.RhoMiddle*tMiddle.radius(x))/                     \
                      (2.*tSurface.radius(x) + tMiddle.radius(x))
                return rho
            else:
                raise RuntimeError('Panel material density calculation '     \
                      'takes 3 inputs: numpy 1-D array x, thickness dist. '  \
                      'tSurface, and thickness dist. tMiddle')

    def getElasticModulus(self,direction=0,*args):
        if self.type == 'ISOTROPIC':
            if hasattr(self,'E'):
                return self.E
            elif hasattr(self,'G') and hasattr(self,'v'):
                return self.G*(2*(1+self.v))
            else:                    
                raise RuntimeError('Not enough information to return E ' \
                      'for isotropic material')
        elif self.type == 'ANISOTROPIC_SHELL':
            if direction == 0:
                return [self.E1,self.E2]
            elif direction == 1:
                return self.E1
            elif direction == 2:
                return self.E2
            else:
                raise NotImplementedError('Only 0, 1, or 2 can be '   \
                      'specified as directions for elastic modulus' \
                      ' where 0 returns values for each direction and '  \
                      '1 and 2 return a value for the specified '    \
                      'direction')
        elif self.type == 'FIXED_RATIO_PANEL':
            E1Surface = self.layer[0].material.getElasticModulus(1)
            E2Surface = self.layer[0].material.getElasticModulus(2)
            E1Middle = self.layer[1].material.getElasticModulus(1)
            E2Middle = self.layer[1].material.getElasticModulus(2)
            ratioSurface = self.layer[0].ratio
            ratioMiddle = self.layer[1].ratio
            E1 = (E1Surface*2*ratioSurface + E1Middle*ratioMiddle)
            E2 = (E2Surface*E2Middle)/(2*ratioSurface*E2Middle + ratioMiddle*E2Surface)               
            if direction == 0:
                return [E1,E2]
            elif direction == 1:
                return E1
            elif direction == 2:
                return E2
            else:
                raise NotImplementedError('Only 0, 1, or 2 can be '   \
                      'specified as directions for elastic modulus' \
                      ' where 0 returns values for each direction and '  \
                      '1 and 2 return a value for the specified '    \
                      'direction')
        else: # panel material
            if len(args) == 3:
                x = args[0]
                tSurface = args[1]
                tMiddle = args[2]
                if direction == 1:
                    E = (self.E1Surface*2*tSurface.radius(x) +           \
                        self.E1Middle*tMiddle.radius(x))/                \
                        (2*tSurface.radius(x) + tMiddle.radius(x))
                    return E
                elif direction == 2:
                    E = (self.E2Surface*self.E2Middle*             \
                        (2*tSurface.radius(x) + tMiddle.radius(x)))/         \
                        (2*tSurface.radius(x)*self.E2Middle +           \
                        tMiddle.radius(x)*self.E2Surface)
                    return E
                else:
                    raise NotImplementedError('Only 1 or 2 can be specified ' \
                          'as directions for panel material elastic modulus')
            else:
                raise RuntimeError('Panel material elastic modulus '         \
                      'takes 3 inputs: numpy 1-D array x, thickness dist. '  \
                      'tSurface, and thickness dist. tMiddle')  

    def getShearModulus(self,direction=0):
        if self.type == 'ISOTROPIC':
            if hasattr(self,'G'):
                return self.G
            elif hasattr(self,'E') and hasattr(self,'v'):
                return self.E/(2*(1+self.v))
            else:
                raise RuntimeError('Not enough information to return G ' \
                      'for isotropic material')
        elif self.type == 'ANISOTROPIC_SHELL':
            return self.G12
        elif self.type == 'FIXED_RATIO_PANEL':
            G12 = self.layer[1].material.getShearModulus()
            print 'WARNING: FIXED_RATIO_PANEL assumes shear modulus of '     \
              'entire panel is equal to middle layer (%s)' % self.name
            return G12
        else: # panel material
            raise RuntimeError('Panel material shear modulus not implemented.')

    def getPoissonRatio(self):
        if self.type == 'ISOTROPIC':
            return self.v
        elif self.type == 'ANISOTROPIC_SHELL':
            return self.v12
        elif self.type == 'FIXED_RATIO_PANEL':
            vSurface = self.layer[0].material.getPoissonRatio()
            vMiddle = self.layer[1].material.getPoissonRatio()
            ratioSurface = self.layer[0].ratio
            ratioMiddle = self.layer[1].ratio
            v = vSurface*ratioSurface*2 + vMiddle*ratioMiddle
            print 'WARNING: FIXED_RATIO_PANEL approximates a global Poisson' \
              ' ratio (%s)' % self.name
            return v
        else: # anisotropic shell, panel material
            raise RuntimeError('Panel material poisson ratio not implemented.')

    def getMutualInfluenceCoefs(self,direction=0):
        if self.type == 'ISOTROPIC':
            #print 'WARNING: Isotropic material does not have coefficients '   \
            #  'of mutual influence (%s)' % self.name
            #return 0.      
            raise NotImplementedError('Isotropic material does not have '    \
                  'coefficients of mutual influence')
        elif self.type == 'ANISOTROPIC_SHELL':
            if direction == 0:
                return [self.nu1, self.nu2]
            elif direction == 1:
                return self.nu1
            elif direction == 2:
                return self.nu2     
            else:
                raise NotImplementedError('Only 1 or 2 can be specified ' \
                      'as directions for mutual influence coefficients')
        elif self.type == 'FIXED_RATIO_PANEL':
            print 'WARNING: FIXED_RATIO_PANEL assumes 0 for coefficients of ' \
              'mutual influence (%s)' % self.name
            if direction == 0:
                return [0.,0.]
            else:
                return 0.
            #raise NotImplementedError('coefficients of mutural influence '   \
            #  'not implemented for FIXED_RATIO_PANEL')            
        else: # anisotropic shell, panel material
            raise RuntimeError('Panel material mutual influence coefs not implemented.')                         
        return 0 
                      
    # Calculate (if necessary) and return thermal conductivity
    def getThermalConductivity(self,direction=0,*args):
        if self.type == 'ISOTROPIC':
            return self.k
        elif self.type == 'ANISOTROPIC_SHELL':
            if direction == 0:
                return [self.k1, self.k2, self.k3]
            if direction == 1:
                return self.k1
            elif direction == 2:
                return self.k2
            elif direction == 3:
                return self.k3
            else:
                raise NotImplementedError('Only 0, 1, 2, or 3 can be '   \
                      'specified as directions for thermal conductivity' \
                      ' where 0 returns values for each direction and '  \
                      '1, 2, and 3 return a value for the specified '    \
                      'direction')
        elif self.type == 'FIXED_RATIO_PANEL':
            k1Surface = self.layer[0].material.getThermalConductivity(1)
            k2Surface = self.layer[0].material.getThermalConductivity(2)
            k3Surface = self.layer[0].material.getThermalConductivity(3)
            k1Middle = self.layer[1].material.getThermalConductivity(1)
            k2Middle = self.layer[1].material.getThermalConductivity(2)
            k3Middle = self.layer[1].material.getThermalConductivity(3)
            ratioSurface = self.layer[0].ratio
            ratioMiddle = self.layer[1].ratio
            k1 = (2*k1Surface*ratioSurface + k1Middle*ratioMiddle)
            k2 = (2*k2Surface*ratioSurface + k2Middle*ratioMiddle)
            k3 = (k3Surface*k3Middle)/(2*ratioSurface*k3Middle + ratioMiddle*k3Surface)                      
            if direction == 0:
                return [k1, k2, k3]
            if direction == 1:
                return k1
            elif direction == 2:
                return k2
            elif direction == 3:
                return k3
            else:
                raise NotImplementedError('Only 0, 1, 2, or 3 can be '   \
                      'specified as directions for thermal conductivity' \
                      ' where 0 returns values for each direction and '  \
                      '1, 2, and 3 return a value for the specified '    \
                      'direction')            
        else: # panel material
            if len(args) == 3:
                x = args[0]
                tSurface = args[1]
                tMiddle = args[2]
                if direction == 1 or direction == 2:
                    kAverageSurf = (self.k1Surface+self.k2Surface)/2
                    kAverageMidd = (self.k1Middle+self.k2Middle)/2
                    k = (kAverageSurf*2*tSurface.radius(x) +           \
                        kAverageMidd*tMiddle.radius(x))/                \
                        (2*tSurface.radius(x) + tMiddle.radius(x))
                    return k
                elif direction == 2:
                    k = (self.k3Surface*self.k3Middle*             \
                        (2*tSurface.radius(x) + tMiddle.radius(x)))/         \
                        (2*tSurface.radius(x)*self.k3Middle +           \
                        tMiddle.radius(x)*self.k3Surface)
                    return k
                else:
                    raise NotImplementedError('Only 1, 2 or 3 '   \
                          'can be specified as directions for properties')
            else:
                raise RuntimeError('Panel material thermal cond. calc. '     \
                      'takes 3 inputs: numpy 1-D array x, thickness dist. '  \
                      'tSurface, and thickness dist. tMiddle')    
    
    def getThermalExpansionCoef(self,direction=0,*args):
        if self.type == 'ISOTROPIC':
            if direction == 12:
                print 'WARNING! Istropic material assumes alpha12 = 0'
                return 0.
            else:
                return self.alpha
        elif self.type == 'ANISOTROPIC_SHELL':
            if direction == 0:
                return [self.alpha1,self.alpha2,self.alpha12]
            elif direction == 1:
                return self.alpha1
            elif direction == 2:
                return self.alpha2
            elif direction == 12:
                return self.alpha12
            else:
                raise NotImplementedError('Only 0, 1, 2, or 12 can be '  \
                      'specified as directions for thermal conductivity' \
                      ' where 0 returns values for each direction and '  \
                      '1, 2, and 12 return a value for the specified '   \
                      'direction')
        elif self.type == 'FIXED_RATIO_PANEL':
            a1Surface = self.layer[0].material.getThermalExpansionCoef(1)
            a2Surface = self.layer[0].material.getThermalExpansionCoef(2)
            a12Surface = self.layer[0].material.getThermalExpansionCoef(12)
            a1Middle = self.layer[1].material.getThermalExpansionCoef(1)
            a2Middle = self.layer[1].material.getThermalExpansionCoef(2)
            a12Middle = self.layer[1].material.getThermalExpansionCoef(12)
            E1Surface = self.layer[0].material.getElasticModulus(1)
            E2Surface = self.layer[0].material.getElasticModulus(2)
            E1Middle = self.layer[1].material.getElasticModulus(1)
            E2Middle = self.layer[1].material.getElasticModulus(2)  
            G12Surface = self.layer[0].material.getShearModulus(12)
            G12Middle = self.layer[1].material.getShearModulus(12)
            ratioMiddle = self.layer[1].ratio
            ratio = 1/ratioMiddle - 1
            alpha1 = (ratio*E1Surface*a1Surface - E1Middle*a1Middle)/(ratio*E1Surface - E1Middle)
            alpha2 = (ratio*E2Surface*a2Surface - E2Middle*a2Middle)/(ratio*E2Surface - E2Middle)
            alpha12 = (ratio*G12Surface*a12Surface - G12Middle*a12Middle)/(ratio*G12Surface - G12Middle)
            #alpha3 = (2*a1Surface*ratioSurface + a1Middle*ratioMiddle)    
            if direction == 0:
                return [alpha1,alpha2,alpha12]
            elif direction == 1:
                return alpha1
            elif direction == 2:
                return alpha2
            elif direction == 12:
                return alpha12
            else:
                raise NotImplementedError('Only 0, 1, 2, or 12 can be '  \
                      'specified as directions for thermal conductivity' \
                      ' where 0 returns values for each direction and '  \
                      '1, 2, and 12 return a value for the specified '   \
                      'direction')                      
        else: # panel material
            if len(args) == 3:
                x = args[0]
                tSurface = args[1]
                tMiddle = args[2]
                if direction == 1:
                    alpha = (self.alpha1Surface*self.E1Surface*      \
                            2*tSurface.radius(x) - self.alpha1Middle*    \
                            self.E1Middle*tMiddle.radius(x))/            \
                            (self.E1Surface*2*tSurface.radius(x) -       \
                            self.E1Middle*tMiddle.radius(x))
                    return alpha
                elif direction == 2:
                    alpha = (self.alpha2Surface*2*tSurface.radius(x) +  \
                            self.alpha2Middle*tMiddle.radius(x))/       \
                            (2*tSurface.radius(x) + tMiddle.radius(x))
                    return alpha
                else:
                    raise NotImplementedError('Only 1 or 2 can be specified ' \
                          'as directions for panel material thermal '        \
                          'expansion coefficient')
            else:
                raise RuntimeError('Panel material thermal exp. coef. '      \
                      'takes 3 inputs: numpy 1-D array x, thickness dist. '  \
                      'tSurface, and thickness dist. tMiddle')  