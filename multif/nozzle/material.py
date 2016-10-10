# -*- coding: utf-8 -*-
"""
Material module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Material:
    def __init__(self, name, matType, structure):
        if matType == 'ISOTROPIC' or matType == 'ANISOTROPIC_SHELL':
            self.name = name
            self.type = matType
        else:
            raise NotImplementedError('Only ISOTROPIC or ANISOTROPIC_SHELL ' \
                  'materials are implemented')
                
        if structure == 'single':
            self.depends_on_t = 0 # material properties DO NOT depend on thickness
        elif structure == 'panel':
            self.depends_on_t = 1 # material properties DO depend on thickness
            if self.type == 'ISOTROPIC':
                raise Warning('Isotropic material for panel structure does'  \
                      'not make any sense')
        else:
            raise NotImplementedError('Only "single" materials or "panel" '  \
                  'material structures are implemented')
    
    def setDensity(self,rho):
        if self.depends_on_t == 0: # single material
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
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
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
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
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
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
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
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
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
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
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
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
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
        if ~self.depends_on_t: # single material
            return self.rho
        else: # panel material
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
        if ~self.depends_on_t: # single material
            if self.type == 'ISOTROPIC':
                if hasattr(self,'E'):
                    return self.E
                elif hasattr(self,'G') and hasattr(self,'v'):
                    return self.G*(2*(1+self.v))
                else:                    
                    raise RuntimeError('Not enough information to return E ' \
                          'for isotropic material')
            else: # ANISOTROPIC_SHELL
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
        if ~self.depends_on_t: # single material
            if self.type == 'ISOTROPIC':
                if hasattr(self,'G'):
                    return self.G
                elif hasattr(self,'E') and hasattr(self,'v'):
                    return self.E/(2*(1+self.v))
                else:
                    raise RuntimeError('Not enough information to return G ' \
                          'for isotropic material')
            else: # ANISOTROPIC_SHELL
                return self.G12
        else: # panel material
            raise RuntimeError('Panel material shear modulus not implemented.')

    def getPoissonRatio(self):
        if self.type == 'ISOTROPIC':
            return self.v
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
            return self.v12
        else: # anisotropic shell, panel material
            raise RuntimeError('Panel material poisson ratio not implemented.')

    def getMutualInfluenceCoefs(self,direction=0):
        if self.type == 'ISOTROPIC':
            #print 'WARNING: Isotropic material does not have coefficients '   \
            #  'of mutual influence (%s)' % self.name
            #return 0.;            
            raise NotImplementedError('Isotropic material does not have '    \
                  'coefficients of mutual influence')
        elif self.type == 'ANISOTROPIC_SHELL' and self.depends_on_t == 0:
            if direction == 0:
                return [self.nu1, self.nu2]
            elif direction == 1:
                return self.nu1
            elif direction == 2:
                return self.nu2     
            else:
                raise NotImplementedError('Only 1 or 2 can be specified ' \
                      'as directions for mutual influence coefficients')
        else: # anisotropic shell, panel material
            raise RuntimeError('Panel material mutual influence coefs not implemented.')                         
        return 0 
                      
    # Calculate (if necessary) and return thermal conductivity
    def getThermalConductivity(self,direction=0,*args):
        if ~self.depends_on_t: # single material
            if self.type == 'ISOTROPIC': # isotropic
                return self.k
            else: # ANISOTROPIC_SHELL
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
        if ~self.depends_on_t: # single material
            if self.type == 'ISOTROPIC':
                return self.alpha
            else: # anisotropic
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