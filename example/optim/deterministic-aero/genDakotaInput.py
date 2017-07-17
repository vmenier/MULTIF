# genDakotaInput.py generates a Dakota input file for deterministic
# optimization of a reduced minimum mass nozzle problem given a 
# single thrust constraint at the top-of-climb operating condition.
# Derivatives are obtained via adjoint sensitivities.
#
# >> dakota -i climb-optim.in -o climb-optim.out
#
# Rick Fenrich 7/17/17

import numpy as np

import linearConstraints

filename = "det-opt"
maxLineWidth = 2000
ncores = 8 # number of cores available

## ---- Set up design variables (wall variables only)
designVar = dict()

# Set up initial values of design variables
designVar['initial'] = np.array([0.2124, 0.2269, 0.2734, 0.3218, 0.3230, 0.3343, 0.3474, 0.4392, 0.4828, 0.5673, 0.6700, 0.3238, 0.2981, 0.2817, 0.2787, 0.2797, 0.2807, 0.2936, 0.2978, 0.3049, 0.3048])

# Set up lower and upper bounds of design variables
lbPerc = 0.8
ubPerc = 1.2
designVar['lower_bound'] = lbPerc*designVar['initial']
designVar['upper_bound'] = ubPerc*designVar['initial']

# Set descriptions for design variables
designVar['description'] = ['c6x', 'c7x', 'c8x', 'c9x', 'c11x', 'c12x', 'c13x', 'c14x', 'c15x', 'c16x', 'c17x', 'c6y', 'c7y', 'c8y', 'c9y', 'c12y', 'c13y', 'c14y', 'c15y', 'c16y', 'c17y']

# Useful indices for generating constraints, etc.
n = len(designVar['initial']) # number of design variables

## ---- Define nonlinear inequality constraints, g(x) <= d
conNLI = dict()
conNLI['lower_bound'] = np.array([21500])
conNLI['upper_bound'] = np.array([50000])
conNLI['description'] = ['thrust']

## ---- Define linear inequality constraints, Ax <= b

(tmp1, A_w, b_w) = linearConstraints.wall(designVar['initial'])
A = A_w # matrix for Ax <= b
b = b_w # RHS vector for Ax <= b
ncon, tmp2 = A.shape

# Scale the linear constraints to unity on RHS
for row in range(ncon):
  if np.abs(b[row]) < 1e-12:
    scale = np.min(np.abs([a for a in A[row,:] if a != 0]))
    # do not scale b[row]
  else:
    scale = np.abs(b[row])
    b[row] = b[row]/scale
  A[row,:] = A[row,:]/scale

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
f.write("# 21 design variables related to shape of inner wall\n")
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
f.write("        scales = 4.6512e-5\n")
f.write("    analytic_gradients\n") # provided by user
f.write("    no_hessians\n")
f.write("\n")
f.flush()

# write interface
f.write("interface\n")
f.write("    id_interface = 'OPTIM_I'\n")
f.write("    analysis_driver = 'python runModel.py -l 2 -n %i -f det_optim_aero.cfg'\n" % ncores)
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

