# genDakotaInput.py generates a Dakota input file for deterministic
# optimization of the small climb problem.
#
# >> dakota -i climb-optim.in -o climb-optim.out
#
# Rick Fenrich 12/1/16

import numpy as np

from multifDomains3D import buildDesignDomain, buildRandomDomain

def concatenateDesignVariables(key,components):
	previousArray = components[0][key]
	for i in range(1,len(components)):
		previousArray = np.hstack((previousArray,components[i][key]))
	return previousArray		

filename = "det-opt"
maxLineWidth = 2000

## ---- Set up design variables
wall = dict()
thermalLayer = dict()
airGap = dict()
loadLayerIn = dict()
loadLayerMid = dict()
loadLayerOut = dict()
stringers = dict()
baffles = dict()

# Set up initial values of design variables
init_data = np.loadtxt("det_mdo_3d.init")
wall['initial'] = init_data[0:18]
thermalLayer['initial'] = init_data[18:31]
airGap['initial'] = init_data[31]
loadLayerIn['initial'] = init_data[32:45]
loadLayerMid['initial'] = init_data[45:58]
loadLayerOut['initial'] = init_data[58:71]
stringers['initial'] = init_data[71:87]
baffles['initial'] = init_data[87:96]

# Set up lower and upper bounds of design variables
# TO DO: set up a check to ensure certain bounds are not < 0 or > 1, otherwise set to delta and 1 - delta
lbPerc = 0.8
ubPerc = 1.2
wall['lower_bound'] = lbPerc*wall['initial']
wall['upper_bound'] = ubPerc*wall['initial']
thermalLayer['lower_bound'] = np.hstack((lbPerc*thermalLayer['initial'][0],12*[0.01]))
thermalLayer['upper_bound'] = np.hstack((ubPerc*thermalLayer['initial'][0],12*[0.05]))
airGap['lower_bound'] = np.array([0.003])
airGap['upper_bound'] = np.array([0.05])
loadLayerIn['lower_bound'] = np.hstack((lbPerc*loadLayerIn['initial'][0],12*[0.0005]))
loadLayerIn['upper_bound'] = np.hstack((ubPerc*loadLayerIn['initial'][0],12*[0.01]))
loadLayerMid['lower_bound'] = np.hstack((lbPerc*loadLayerMid['initial'][0],12*[0.0064]))
loadLayerMid['upper_bound'] = np.hstack((ubPerc*loadLayerMid['initial'][0],12*[0.0159]))
loadLayerOut['lower_bound'] = np.hstack((lbPerc*loadLayerOut['initial'][0],12*[0.0005]))
loadLayerOut['upper_bound'] = np.hstack((ubPerc*loadLayerOut['initial'][0],12*[0.01]))
stringers['lower_bound'] = np.hstack((lbPerc*stringers['initial'][0:4],12*[0.002]))
stringers['upper_bound'] = np.hstack((ubPerc*stringers['initial'][0:4],12*[0.01]))
baffles['lower_bound'] = np.hstack((lbPerc*baffles['initial'][0:4],5*[0.0074]))
baffles['upper_bound'] = np.hstack((ubPerc*baffles['initial'][0:4],5*[0.0359]))

# Set descriptions for design variables
wall['description'] = ['wall_' + str(i) for i in range(18)]
thermalLayer['description'] = ['thermal_' + str(i) for i in range(13)]
airGap['description'] = ['air_gap']
loadLayerIn['description'] = ['load_in_' + str(i) for i in range(13)]
loadLayerMid['description'] = ['load_mid_' + str(i) for i in range(13)]
loadLayerOut['description'] = ['load_out_' + str(i) for i in range(13)]
stringers['description'] = ['stringer_' + str(i) for i in range(16)]
baffles['description'] = ['baffle_' + str(i) for i in range(9)]

# Concatenate all design variables
designVar = dict()
componentList = [wall, thermalLayer, airGap, loadLayerIn, loadLayerMid, loadLayerOut, stringers, baffles]
designVar['initial'] = concatenateDesignVariables('initial', componentList)
designVar['lower_bound'] = concatenateDesignVariables('lower_bound', componentList)
designVar['upper_bound'] = concatenateDesignVariables('upper_bound', componentList)
designVar['description'] = concatenateDesignVariables('description', componentList)

## ---- Define uncertain parameters

## ---- Define nonlinear inequality constraints, g(x) <= d
conNLI = dict()
conNLI['lower_bound'] = np.array([21500] + 14*[0])
conNLI['upper_bound'] = np.array([100000] + 14*[1])
conNLI['description'] = ['thrust', 'thermal_layer_temp', 'load_layer_in_temp', 
  'load_layer_mid_temp', 'load_layer_out_temp', 'thermal_layer_fc', 
  'load_layer_in_fc', 'load_layer_mid_fc', 'load_layer_out_fc', 'stringers_fc', 
  'baffle1_fc', 'baffle2_fc', 'baffle3_fc', 'baffle4_fc', 'baffle5_fc']

## ---- Define linear inequality constraints, Ax <= b
domain = buildDesignDomain(output='verbose')
A = domain.A # linear inequality constraints
b = domain.b # linear inequality constraint RHS
nLinearCon = A.shape[0]

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
f.write("# Deterministic optimization of nozzle problem\n")
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
f.write("        max_function_evaluations = 10000\n")
f.write("        function_precision = 1.e-5\n")
f.write("        convergence_tolerance = 3.2e-3\n")
f.write("        constraint_tolerance = 3.2e-3\n")
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
for i in range(0,nLinearCon):
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
f.write("        scales = 4.6512e-5 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n")
f.write("    numerical_gradients\n")
f.write("        method_source dakota\n")
f.write("        interval_type forward\n")
f.write("        fd_gradient_step_size = 3.2e-3\n")
f.write("    no_hessians\n")
f.write("\n")
f.flush()

# write interface
f.write("interface\n")
f.write("    id_interface = 'OPTIM_I'\n")
f.write("    analysis_driver = 'python runModel.py -l 0 -f det_mdo_3d.cfg'\n")
f.write("        fork asynchronous evaluation_concurrency = 4\n")
f.write("        work_directory named 'tmp'\n")
#f.write("        link_files '*'\n")
f.write("        link_files 'runModel.py' 'det_mdo_3d.cfg'\n")
f.write("        directory_tag\n")
f.write("        directory_save\n")
#f.write("        file_tag\n")
f.write("        file_save\n")
f.write("        parameters_file = 'params.in'\n")
f.write("        results_file = 'results.out'\n")
#f.write("    failure_capture\n")
#f.write("        recover 1.e3 0.e0\n")

f.close()

