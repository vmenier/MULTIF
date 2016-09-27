# -*- coding: utf-8 -*-
"""
Material module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Material:
    #def __init__(self, rho, k, alpha, E, v):
    def __init__(self, matType, structure):
        if matType == 'isotropic':
            self.n = 1 # number of property entries for each property
        elif matType == 'anisotropic':
            self.n = 2 # number of property entries for each property
        else:
            raise NotImplementedError('Only isotropic or anisotropic in 2 '  \
                  'direction materials are implemented')
                
        if structure == 'single':
            self.depends_on_t = 0 # material properties DO NOT depend on thickness
        elif structure == 'panel':
            self.depends_on_t = 1 # material properties DO depend on thickness
            if self.n == 1:
                self.n = 2
                raise Warning('Panel structure material reset from '         \
                      'isotropic to anisotropic')
        else:
            raise NotImplementedError('Only "single" materials or "panel" '  \
                  'material structures are implemented')
    
    def setDensity(self,rho):
        if self.depends_on_t == 0: # single material
            try:
                self.rho = float(rho)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for density')
        else: # panel material
            try:
                self.rhoSurface = float(rho[0])
                self.rhoMiddle = float(rho[1])
            except:
                raise RuntimeError('Panel structure material must have 2 '   \
                      'values provided for density')
        return 0
       
    def setThermalConductivity(self,k):
        if self.n == 1: # isotropic material
            try:
                self.k = float(k)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for thermal conductivity')
        elif self.n == 2 and self.depends_on_t == 0: # anisotropic, single material
            if isinstance(k,float): # 1 value provided
                self.kAxial = k
                self.kRadial = k
            else:
                try:
                    self.kAxial = float(k[0])
                    self.kRadial = float(k[1])
                except:
                    raise RuntimeError('Anisotropic material must have '     \
                          '2 values provided for thermal conductivity')
        else: # anisotropic, panel material
            try:
                self.kAxialSurface = float(k[0])
                self.kRadialSurface = float(k[1])
                self.kAxialMiddle = float(k[2])
                self.kRadialMiddle = float(k[3])
            except:
                raise RuntimeError('Anisotropic panel material must have 4 ' \
                      'values provided for thermal conductivity')
        return 0

    def setThermalExpansionCoef(self,alpha):
        if self.n == 1: # isotropic material
            try:
                self.alpha = float(alpha)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for thermal expansion coef')
        elif self.n == 2 and self.depends_on_t == 0: # anisotropic, single material
            if isinstance(alpha,float): # 1 value provided
                self.alphaAxial = alpha
                self.alphaRadial = alpha
            else:
                try:
                    self.alphaAxial = float(alpha[0])
                    self.alphaRadial = float(alpha[1])
                except:
                    raise RuntimeError('Anisotropic material must have 2 '   \
                          'values provided for thermal expansion coef')
        else: # anisotropic, panel material
            try:
                self.alphaAxialSurface = float(alpha[0])
                self.alphaRadialSurface = float(alpha[1])
                self.alphaAxialMiddle = float(alpha[2])
                self.alphaRadialMiddle = float(alpha[3])
            except:
                raise RuntimeError('Anisotropic panel material must have 4'  \
                      'values provided for thermal expansion coef')
        return 0
   
    def setElasticModulus(self,E):
        if self.n == 1: # isotropic material
            try:
                self.E = float(E)
            except:
                raise RuntimeError('Isotropic material must have 1 value '   \
                      'provided for elastic modulus')
        elif self.n == 2 and self.depends_on_t == 0: # anisotropic, single material
            if isinstance(E,float): # 1 value provided
                self.EAxial = E
                self.ERadial = E
            else:
                try:
                    self.EAxial = float(E[0])
                    self.ERadial = float(E[1])
                except:
                    raise RuntimeError('Anisotropic material must have 2 '   \
                          'values provided for elastic modulus')
        else: # anisotropic, panel material
            try:
                self.EAxialSurface = float(E[0])
                self.ERadialSurface = float(E[1])
                self.EAxialMiddle = float(E[2])
                self.ERadialMiddle = float(E[3])
            except:
                raise RuntimeError('Anisotropic panel material must have 4'  \
                      'values provided for elastic modulus')
        return 0
        
    def setPoissonRatio(self,v):
        try:
            self.v = float(v)
        except:
            raise NotImplementedError('Only isotropic single material '      \
                  'Poisson ratio implemented. Use 1 value for panel '        \
                  'materials and anisotropic materials')
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
    
    # Calculate (if necessary) and return thermal conductivity
    def getThermalConductivity(self,direction,*args):
        if ~self.depends_on_t: # single material
            if self.n == 1: # isotropic
                return self.k
            else: # anisotropic
                if direction == 'axial':
                    return self.kAxial
                elif direction == 'radial':
                    return self.kRadial
                else:
                    raise NotImplementedError('Only "axial" and "radial" '   \
                          'can be specified as directions for properties')
        else: # panel material
            if len(args) == 3:
                x = args[0]
                tSurface = args[1]
                tMiddle = args[2]
                if direction == 'axial':
                    k = (self.kAxialSurface*2*tSurface.radius(x) +           \
                        self.kAxialMiddle*tMiddle.radius(x))/                \
                        (2*tSurface.radius(x) + tMiddle.radius(x))
                    return k
                elif direction == 'radial':
                    k = (self.kRadialSurface*self.kRadialMiddle*             \
                        (2*tSurface.radius(x) + tMiddle.radius(x)))/         \
                        (2*tSurface.radius(x)*self.kRadialMiddle +           \
                        tMiddle.radius(x)*self.kRadialSurface)
                    return k
                else:
                    raise NotImplementedError('Only "axial" and "radial" '   \
                          'can be specified as directions for properties')
            else:
                raise RuntimeError('Panel material thermal cond. calc. '     \
                      'takes 3 inputs: numpy 1-D array x, thickness dist. '  \
                      'tSurface, and thickness dist. tMiddle')    
    
    def getThermalExpansionCoef(self,direction,*args):
        if ~self.depends_on_t: # single material
            if self.n == 1: # isotropic
                return self.alpha
            else: # anisotropic
                if direction == 'axial':
                    return self.alphaAxial
                elif direction == 'radial':
                    return self.alphaRadial
                else:
                    raise NotImplementedError('Only "axial" and "radial" '   \
                          'can be specified as directions for properties')
        else: # panel material
            if len(args) == 3:
                x = args[0]
                tSurface = args[1]
                tMiddle = args[2]
                if direction == 'axial':
                    alpha = (self.alphaAxialSurface*self.EAxialSurface*      \
                            2*tSurface.radius(x) - self.alphaAxialMiddle*    \
                            self.EAxialMiddle*tMiddle.radius(x))/            \
                            (self.EAxialSurface*2*tSurface.radius(x) -       \
                            self.EAxialMiddle*tMiddle.radius(x))
                    return alpha
                elif direction == 'radial':
                    alpha = (self.alphaRadialSurface*2*tSurface.radius(x) +  \
                            self.alphaRadialMiddle*tMiddle.radius(x))/       \
                            (2*tSurface.radius(x) + tMiddle.radius(x))
                    return alpha
                else:
                    raise NotImplementedError('Only "axial" and "radial" '   \
                          'can be specified as directions for properties')
            else:
                raise RuntimeError('Panel material thermal exp. coef. '      \
                      'takes 3 inputs: numpy 1-D array x, thickness dist. '  \
                      'tSurface, and thickness dist. tMiddle')  
                      
    def getElasticModulus(self,direction,*args):
        if ~self.depends_on_t: # single material
            if self.n == 1: # isotropic
                return self.E
            else: # anisotropic
                if direction == 'axial':
                    return self.EAxial
                elif direction == 'radial':
                    return self.ERadial
                else:
                    raise NotImplementedError('Only "axial" and "radial" '   \
                          'can be specified as directions for properties')
        else: # panel material
            if len(args) == 3:
                x = args[0]
                tSurface = args[1]
                tMiddle = args[2]
                if direction == 'axial':
                    E = (self.EAxialSurface*2*tSurface.radius(x) +           \
                        self.EAxialMiddle*tMiddle.radius(x))/                \
                        (2*tSurface.radius(x) + tMiddle.radius(x))
                    return E
                elif direction == 'radial':
                    E = (self.ERadialSurface*self.ERadialMiddle*             \
                        (2*tSurface.radius(x) + tMiddle.radius(x)))/         \
                        (2*tSurface.radius(x)*self.ERadialMiddle +           \
                        tMiddle.radius(x)*self.ERadialSurface)
                    return E
                else:
                    raise NotImplementedError('Only "axial" and "radial" '   \
                          'can be specified as directions for properties')
            else:
                raise RuntimeError('Panel material elastic modulus '         \
                      'takes 3 inputs: numpy 1-D array x, thickness dist. '  \
                      'tSurface, and thickness dist. tMiddle')  
                
    def getPoissonRatio(self):
        return self.v
        
#        self.matType = matType
#        self.rho = rho # kg/m^3, density
#        self.k = k # W/m*K, thermal conductivity of wall
#        self.alpha = alpha # 1/K, coeff. of thermal expansion
#        self.E = E # Pa, elastic modulus
#        self.v = v # Poisson's ratio
    
