# genDakotaInput.py generates a Dakota input file for deterministic
# optimization of a reduced minimum mass nozzle problem given a 
# single thrust constraint at the top-of-climb operating condition.
# Derivatives are obtained via adjoint sensitivities.
#
# >> dakota -i climb-optim.in -o climb-optim.out
#
# Rick Fenrich 10/4/17

import numpy as np

# The imports below require access to the following files from
# MULTIF/example/domains:
#    domains.py
#    linearConstraints.py
#    multifDomains2D.py
#    util.py

from linearConstraints import bspline
from linearConstraints import cleanupConstraintMatrix

from domains import LinIneqDomain

filename = "det-opt"
maxLineWidth = 2000
ncores = 8 # number of cores available

## ---- Set up design variables (wall variables only)
designVar = dict()

# Set up initial values of design variables
designVar['initial'] = np.array([0.542184, 0.861924, 1.072944, 1.211161, 1.311161, 1.408983, 1.528983, 1.723828, 2.086573, 0.417017, 0.365097, 0.301792, 0.267426, 0.277426, 0.332508, 0.385631])

# Set up lower and upper bounds of design variables
lbPerc = -np.inf
ubPerc = np.inf
designVar['lower_bound'] = lbPerc*designVar['initial']
designVar['upper_bound'] = ubPerc*designVar['initial']

# Set descriptions for design variables
designVar['description'] = ['x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7']

# Useful indices for generating constraints, etc.
n = len(designVar['initial']) # number of design variables

## ---- Define nonlinear inequality constraints, g(x) <= d
conNLI = dict()
conNLI['lower_bound'] = np.array([21500])
conNLI['upper_bound'] = np.array([50000])
conNLI['description'] = ['thrust']

## ---- Define linear inequality constraints, Ax <= b

WALL_COEFS= (0.0000, 0.0000, 0.1, 0.3, 0.7, 1.0, 1.3, 
             1.3500, 1.3500, 1.4000, 1.5000, 1.6000, 1.8000, 2.3371, 2.3371, 
             0.4395, 0.4395, 0.4395, 0.4, 0.34, 0.31, 0.27, 
             0.2700, 0.2700, 0.2700, 0.3, 0.33, 0.38, 0.3955, 0.3955)
WALL_COEFS_DV= (0, 0, 0, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 0, 0, 
                0, 0, 0, 10, 11, 12, 13, 13, 13, 13, 14, 15, 16, 0, 0)
x_wall = np.array([0.542184, 0.861924, 1.072944, 1.211161, 1.311161, 
                   1.408983, 1.528983, 1.723828, 2.086573, 0.417017, 
                   0.365097, 0.301792, 0.267426, 0.277426, 0.332508, 
                   0.385631])
A, b = bspline(WALL_COEFS, WALL_COEFS_DV, 7, (-0.3,0.005,-0.025,0.35), 
                 xLimits=[None,2.3], delta=0.1, throatIsLowest=1, 
                 minThroat=0.2, output='verbose') 
A, b = cleanupConstraintMatrix(Alist=[A],blist=[b])
inner_wall_domain = LinIneqDomain(A, np.squeeze(b), center = x_wall)
b = np.squeeze(b)
shapea = A.shape
ncon = shapea[0]

# Check feasibility of linear constraints at starting point
b2 = np.dot(A,np.transpose(designVar['initial']))
bDiff = np.transpose(b)-b2
nConInfeasible = 0
for i in range(0,len(bDiff[:])):
  if bDiff[i] < 0:
    nConInfeasible += 1
if nConInfeasible > 0:
  print("{} linear inequality constraints infeasible at starting point\n".format(nConInfeasible))
  print(bDiff)
else:
  print("no linear inequality infeasibilities!\n")

## ---- Write Dakota input file

filein = "%s.in" %filename
fileout = "%s.out" %filename
f = open(filein,"w")

# write header information
f.write("# Deterministic optimization of reduced aero only design problem\n")
f.write("#     Minimize Mass of Wall only s.t. Thrust < 21500 N\n")
f.write("# 16 design variables related to shape of inner wall\n")
f.write("# Usage:\n")
f.write("# dakota -i {} -o {}\n".format(filein,fileout))
f.write("\n")
f.flush()

# write environment
f.write("environment\n")
#f.write("    graphics\n")
f.write("    output_precision = 16\n")
f.write("    method_pointer = 'OPTIM'\n")
f.write("\n")
f.flush()

# write method
f.write("method\n")
f.write("    id_method = 'OPTIM'\n")
f.write("    model_pointer = 'OPTIM_M'\n")
f.write("        npsol_sqp\n")
f.write("        max_iterations = 100\n")
f.write("        max_function_evaluations = 100\n")
f.write("        function_precision = 1.e-6\n")
f.write("        convergence_tolerance = 1.e-3\n")
f.write("        constraint_tolerance = 1.e-3\n")
#f.write("        linesearch_tolerance = 0.9\n")
f.write("\n")
f.flush()

# write model
f.write("model\n")
f.write("    id_model = 'OPTIM_M'\n")
f.write("    single\n") # for deterministic optimization only
f.write("        variables_pointer = 'OPTIM_V'\n")
f.write("        responses_pointer = 'OPTIM_R'\n")
f.write("        interface_pointer = 'OPTIM_I'\n")
f.write("\n")
f.flush()

# write variables
f.write("variables\n")
f.write("    id_variables = 'OPTIM_V'\n")
f.write("    continuous_design = {}\n".format(designVar['initial'].size))
f.write("        initial_point    {}\n".format(np.array_str(designVar['initial'],max_line_width=maxLineWidth)[1:-1]))
f.write("        lower_bounds     {}\n".format(np.array_str(designVar['lower_bound'],max_line_width=maxLineWidth,precision=4)[1:-1]))
f.write("        upper_bounds     {}\n".format(np.array_str(designVar['upper_bound'],max_line_width=maxLineWidth,precision=4)[1:-1]))
f.write("        descriptors      {}\n".format("'" + "' '".join(map(str,designVar['description'])) + "'"))
f.write("        scale_types = 'auto'\n")
f.write("    linear_inequality_constraint_matrix =\n")
for i in range(0,ncon):
  f.write("        {}\n".format(np.array_str(A[i,:],max_line_width=maxLineWidth)[1:-1]))
f.write("    linear_inequality_upper_bounds = %0.8f\n" % b[0])
for i in range(1,len(b)):
  f.write("                                     %0.8f\n" % b[i])
f.write("\n")
f.flush()

# write responses
f.write("responses\n")
f.write("    id_responses = 'OPTIM_R'\n")
f.write("    objective_functions = 1\n")
f.write("    primary_scale_types = 'value'\n") # scale objective function (mass) by 100
f.write("    primary_scales = 0.01\n")
f.write("    nonlinear_inequality_constraints = {}\n".format(len(conNLI['lower_bound'])))
f.write("        lower_bounds = {}\n".format(np.array_str(conNLI['lower_bound'],max_line_width=maxLineWidth,precision=4)[1:-1]))
f.write("        upper_bounds = {}\n".format(np.array_str(conNLI['upper_bound'],max_line_width=maxLineWidth,precision=4)[1:-1]))
f.write("        scale_types = 'value'\n")
f.write("        scales = 2.15e-5\n")
f.write("    analytic_gradients\n") # provided by user
f.write("    no_hessians\n")
f.write("\n")
f.flush()

# write interface
f.write("interface\n")
f.write("    id_interface = 'OPTIM_I'\n")
f.write("    analysis_driver = 'pythona runModel.py -l 0 -n %i -f det_optim_aero.cfg'\n" % ncores)
#f.write("        fork asynchronous evaluation_concurrency = 4\n")
f.write("        fork\n")
f.write("        work_directory named 'tmp'\n")
#f.write("        link_files '*'\n")
f.write("        link_files 'runModel.py' 'det_optim_aero.cfg'\n")
f.write("        directory_tag\n")
f.write("        directory_save\n")
#f.write("        file_tag\n")
f.write("        file_save\n")
f.write("        parameters_file = 'params.in'\n")
f.write("        results_file = 'results.out'\n")
#f.write("    failure_capture\n")
#f.write("        recover 1.e3 0.e0\n")

f.close()

