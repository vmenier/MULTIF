"""
Response class for holding response values and gradients.

Rick Fenrich 3/18/18
"""

import numpy as np

class Response(object):
    """
    Generic response object.
    """

    def __init__(self): 
        self.names = []
        self.analysis = {} # analysis name which generates this response
        self.value = {}
        self.gradient = {}
        self.kind = {} # 'scalar' or 'field' for each response
        self.location = {} # record location of responses for field responses
        self.analyses = [] # list of analysis names available for generating
                           # responses
        self.analysisDependence = {} # list of analyses a response is dependent
                           # on for each response
        return

    def __len__(self):
        return len(self.names)

    # def __repr__(self):
    #     """
    #     Print information about responses.
    #     """
    #     string = 'Not implemented yet'
    #     # for n in self.names:
    #     #     string += '%s (%s) requires %s analysis.' % (n,self.kind[n],self.analysis[n])
    #     #     string += '    Dependent on:'
    #     #     string += self.analysisDependence[n]
    #     #     string += '    Value:'
    #     #     string += self.value[n]
    #     #     string += '    Gradients:'
    #     #     string += self.gradient[n]
    #         string += '\n'
    #     return string

    def addAnalyses(self, analyses):
        """
        Add analysis names in analyses list. These analyses can be used to 
        generate responses.
        """
        for a in analyses:
            if a not in self.analyses:
                self.analyses.append(a)
        return
    
    def addResponse(self, name, analysis=None, kind='scalar', dependencies=[], 
        location=[]):
        """
        Add a response with specified name.
        """
        if kind == 'field':
            assert len(location) > 0
        assert analysis in self.analyses
        self.names.append(name)
        self.analysis[name] = analysis
        self.value[name] = None
        self.gradient[name] = None
        self.kind[name] = kind
        self.location[name] = np.array(location)
        self.analysisDependence[name] = []
        for d in dependencies:
            if d in self.analyses:
                self.analysisDependence[name].append(d)
            else:
                raise KeyError("%s is not a valid analysis type." % d)
        return

    def exists(self, name):
        """
        Check if provided response exists. Return name if it does, otherwise 
        empty string.
        """
        if name in self.names:
            return name
        else:
            return ''

    def getAnalyses(self):
        """
        Get list of analyses that can generate responses.
        """
        return self.analyses

    def getAnalysis(self, name):
        """
        Get analysis name that generates named QoI.
        """
        return self.analysis[name]

    def getDependencies(self, name):
        """
        Get list of analysis dependencies for named response.
        """
        return self.analysisDependence[name]

    def getKind(self, name):
        """
        Get kind of named response ('scalar' or 'field').
        """
        return self.kind[name]

    def getLocation(self, name):
        """
        Get locations where named field variable should be evaluated.
        """
        assert self.kind[name] == 'field'
        return self.location[name]

    def getValue(self, name):
        """
        Get response value corresponding to name.
        """
        return self.value[name]

    def getGradient(self, name):
        """
        Get response gradient corresponding to name.
        """
        return self.gradient[name]

    def setValue(self, name, value):
        """
        Set value for given response.
        """
        if self.kind[name] == 'field':
            assert isinstance(value, np.ndarray) or value == None
        else:
            assert isinstance(value, float) or value == None
        if name in self.names:
            self.value[name] = value
        else:
            raise KeyError("%s is not a valid response." % name)
        return
    
    def setGradient(self, name, gradient):
        """
        Set gradient for given response. 
        """
        assert isinstance(gradient, np.ndarray) or gradient == None
        if name in self.names:
            self.gradient[name] = gradient
        else:
            raise KeyError("%s is not a valid response." % name)
        return

    def remove(self, name):
        """
        Remove response from list of responses.
        """
        self.value.pop(name, None)
        self.gradient.pop(name, None)
        self.kind.pop(name, None)
        self.location.pop(name, None)
        self.analysis.pop(name, None)
        self.analysisDependence.pop(name, None)
        self.names.remove(name)
        return

    @property
    def dependencyMatrix(self):
        raise NotImplementedError()


class NozzleResponse(Response):
    """
    Predefined set of responses for nozzle problem.
    """

    def __init__(self): 
        self.names = []
        self.value = {}
        self.gradient = {}
        self.kind = {} # 'scalar' or 'field' for each response
        self.analyses = [] # list of analysis names available for generating
                           # responses
        self.analysisDependence = {} # list of analyses a response is dependent
                           # on for each response

        # Categorized lists for allowable QoI base names & prefixes/suffixes
        self.sst_tags = ['SSTC1', 'SSTC2', 'SSTC3', 'SSTP1C1', 'SSTP1C2']
        self.agglom_tags = ['PN', 'KS', 'MAX', 'MIN']
        self.geo_scalar = ['MASS', 'MASS_WALL_ONLY', 'VOLUME']
        self.aero_scalar = ['THRUST', 'WALL_PRES_AVG', 'WALL_TEMP_AVG', 
                        'SU2_RESIDUAL']
        self.aero_field = ['WALL_PRESSURE', 'PRESSURE', 'VELOCITY']
        self.struct_scalar = []
        self.struct_field = ['TOTAL_STRESS', 'MECHANICAL_STRESS',
                        'THERMAL_STRESS', 'FAILURE_CRITERIA']
        self.therm_scalar = []
        self.therm_field = ['TEMPERATURE', 'TEMP_RATIO', 'WALL_TEMPERATURE']

        # Build more general lists for checking QoI names
        self.scalarNames = self.geo_scalar + self.aero_scalar + \
            self.struct_scalar + self.therm_scalar
        self.fieldNames = self.aero_field + self.struct_field + self.therm_field
        self.generalNames = self.geo_scalar + self.aero_scalar + \
            self.aero_field + self.struct_scalar + \
            self.struct_field + self.therm_scalar + self.therm_field
        self.thermstructNames = self.struct_scalar + self.therm_scalar + \
            self.struct_field + self.therm_field

        return

    def exists(self, name, comp_tags):
        """
        Check if provided response exists. Return standard response name if it 
        does (and name of analysis that generates it), otherwise return empty string.

        Arguments:
        comp_tags: list of component names which can be used in thermal/
            structural QoI
        """

        # Check for SST specific tags first
        sst_suffix = [e for e in self.sst_tags if e in name]
        if len(sst_suffix) > 1:
            #sst_suffix = max(sst_suffix, key=len) # pick longest match
            raise ValueError("Only one SST tag can be specified for response.")
        elif len(sst_suffix) == 1:
            sst_suffix = sst_suffix[0]
        else:
            sst_suffix = ''         

        # Check for agglomeration tags next
        agglom_prefix = [e for e in self.agglom_tags if e in name]
        if len(agglom_prefix) > 1:
            raise ValueError("Only one agglomeration tag can be specified for response.")
        elif len(agglom_prefix) == 1:
            agglom_prefix = agglom_prefix[0]
        else:
            agglom_prefix = ''   

        # Check for component specific tags next
        comp_prefix = [e for e in comp_tags if e in name]
        if len(comp_prefix) > 1:
            raise ValueError("Only one component name can be specified for response.")
        elif len(comp_prefix) == 1:
            comp_prefix = comp_prefix[0]
        else:
            comp_prefix = ''   

        # Remove prefixes and suffixes from name
        remainder = name.replace(sst_suffix,'')
        remainder = remainder.replace(agglom_prefix,'')
        remainder = remainder.replace(comp_prefix,'')

        # Check for base QoI name
        base_name = [e for e in self.generalNames if e in remainder]
        if len(base_name) > 1:
            base_name = sorted(base_name, key=len)[-1] # longest matching name
        elif len(base_name) == 1:
            base_name = base_name[0]
        else:
            raise ValueError("%s not accepted as output QoI." % name)

        # Ensure there is nothing extra in base_name
        remainder = remainder.replace(base_name,'')
        remainder = remainder.replace('_','')
        if len(remainder) > 0:
            raise ValueError("%s not accepted as output QoI." % name)

        # Check for improper mixing and matching of prefixes, suffixes, and base
        # QoI names
        if comp_prefix != '' and base_name not in self.thermstructNames:
            # Specifying component name only acceptable for thermal/structural QoI
            raise ValueError("%s not accepted as output QoI." % name)
        if sst_suffix != '' and base_name in self.geo_scalar:
            # Specifying MASS etc. with RANS perturbation tag not acceptable
            raise ValueError("%s not accepted as output QoI." % name)
        if agglom_prefix != '' and base_name not in self.fieldNames:
            # Specifying agglomeration tag for scalar variable not acceptable
            raise ValueError("%s not accepted as output QoI." % name)

        # Determine type of QoI
        if base_name in self.geo_scalar:
            analysis = 'MASS'
        elif base_name in self.aero_scalar + self.aero_field:
            analysis = '_'.join(('AERO',sst_suffix)) if sst_suffix != '' else 'AERO'
        elif base_name in self.struct_scalar + self.struct_field:
            analysis = '_'.join(('STRUCTURAL',sst_suffix)) if sst_suffix != '' else 'STRUCTURAL'
        elif base_name in self.therm_scalar + self.therm_field:
            analysis = '_'.join(('THERMAL',sst_suffix)) if sst_suffix != '' else 'THERMAL'

        # Build standard name and return standard name
        joinList = [comp_prefix, agglom_prefix, base_name, sst_suffix]
        standard_name = '_'.join([e for e in joinList if e != ''])

        return standard_name, analysis

    def getDependencies(self, name, method):
        """
        Assuming that the provided name has already been checked that it is
        acceptable for the nozzle, provide a list of names of analyses that
        the named QoI is dependent upon. Possible analyses include:

        MASS, AERO, THERMAL, STRUCTURAL
        
        The last 3 may also have SSTC1, SSTC2, SSTC3, SSTP1C1, or SSTP1C2 appended to them
        (spaced by an underscore) according to whether a RANS perturbation 
        analysis is desired.

        Arguments:
        method: 'EULER', 'RANS' or 'NONIDEALNOZZLE'
        """

        # Check for MASS analysis requirement
        mass_analysis = [e for e in self.geo_scalar if e in name]

        # Check for AERO analysis requirement
        aero_analysis = [e for e in self.aero_scalar+self.aero_field if e in name]

        # Check for THERMAL analysis requirement
        therm_analysis = [e for e in self.therm_scalar+self.therm_field if e in name]

        # Check for STRUCTURAL analysis requirement
        struct_analysis = [e for e in self.struct_scalar+self.struct_field if e in name]

        # Check for RANS perturbation requirements
        sst_analysis = [e for e in self.sst_tags if e in name]

        # Build list of dependencies
        dependencyList = set()

        assert len(sst_analysis) <= 1 # If longer, than multiple matches are found
        if len(sst_analysis) == 1:
            #sst_analysis = max(sst_analysis, key=len) # pick longest match
            sst_suffix = ''.join(('_',sst_analysis[0]))
        else:
            sst_suffix = ''

        if len(mass_analysis) > 0:
            dependencyList.add('MASS')
        
        if len(aero_analysis) > 0:
            dependencyList.add(''.join(('AERO',sst_suffix)))
            if method == 'NONIDEALNOZZLE':
                #dependencyList.add(''.join(('THERMAL',sst_suffix)))
                pass
            else:
                pass # no thermal coupling for other aero analyses
        
        if len(therm_analysis) > 0:
            dependencyList.add(''.join(('AERO',sst_suffix)))
            dependencyList.add(''.join(('THERMAL',sst_suffix)))
        
        if len(struct_analysis) > 0:
            dependencyList.add(''.join(('AERO',sst_suffix)))
            dependencyList.add(''.join(('THERMAL',sst_suffix)))  
            dependencyList.add(''.join(('STRUCTURAL',sst_suffix)))  

        return list(dependencyList)

    def getKind(self, name):
        """
        Assuming that the provided name has already been checked that it is
        acceptable for the nozzle, return whether the QoI is a 'scalar' or 
        'field' variable.
        """
        kind = [e for e in self.fieldNames if e in name]
        agglom = [e for e in self.agglom_tags if e in name]
        if len(kind) > 0 and len(agglom) == 0:
            return 'field'
        else:
            return 'scalar'

    def getAnalyses(self):
        """
        Return list of analyses.
        """
        analyses = ['MASS']
        ats = ['AERO', 'THERMAL', 'STRUCTURAL']
        analyses += ats
        for k in ats:
            for s in self.sst_tags:
                analyses.append('_'.join((k,s)))
        return analyses

