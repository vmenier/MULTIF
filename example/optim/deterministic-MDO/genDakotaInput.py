# genDakotaInput.py generates a Dakota input file for deterministic
# optimization of the small climb problem.
#
# >> dakota -i climb-optim.in -o climb-optim.out
#
# Rick Fenrich 12/1/16

import numpy as np

import linearConstraints

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
wall['initial'] = np.array([0.2124, 0.2269, 0.2734, 0.3218, 0.3230, 0.3343, 0.3474, 0.4392, 0.4828, 0.5673, 0.6700, 0.3238, 0.2981, 0.2817, 0.2787, 0.2797, 0.2807, 0.2936, 0.2978, 0.3049, 0.3048])
thermalLayer['initial'] = np.array([0.3, 0.6, 0.03, 0.03, 0.03, 0.03])
airGap['initial'] = np.array([0.005])
loadLayerIn['initial'] = np.array([0.3, 0.6, 0.002, 0.002, 0.002, 0.002])
loadLayerMid['initial'] = np.array([0.3, 0.6, 0.013, 0.013, 0.013, 0.013])
loadLayerOut['initial'] = np.array([0.3, 0.6, 0.002, 0.002, 0.002, 0.002])
stringers['initial'] = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
baffles['initial'] = np.array([0.2, 0.4, 0.6, 0.8, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

# Set up lower and upper bounds of design variables
# TO DO: set up a check to ensure certain bounds are not < 0 or > 1, otherwise set to delta and 1 - delta
lbPerc = 0.8
ubPerc = 1.2
wall['lower_bound'] = lbPerc*wall['initial']
wall['upper_bound'] = ubPerc*wall['initial']
thermalLayer['lower_bound'] = np.hstack((lbPerc*thermalLayer['initial'][0:2],4*[0.01]))
thermalLayer['upper_bound'] = np.hstack((ubPerc*thermalLayer['initial'][0:2],4*[0.05]))
airGap['lower_bound'] = np.array([0.003])
airGap['upper_bound'] = np.array([0.01])
loadLayerIn['lower_bound'] = np.hstack((lbPerc*loadLayerIn['initial'][0:2],4*[0.0005]))
loadLayerIn['upper_bound'] = np.hstack((ubPerc*loadLayerIn['initial'][0:2],4*[0.01]))
loadLayerMid['lower_bound'] = np.hstack((lbPerc*loadLayerMid['initial'][0:2],4*[0.0064]))
loadLayerMid['upper_bound'] = np.hstack((ubPerc*loadLayerMid['initial'][0:2],4*[0.0159]))
loadLayerOut['lower_bound'] = np.hstack((lbPerc*loadLayerOut['initial'][0:2],4*[0.0005]))
loadLayerOut['upper_bound'] = np.hstack((ubPerc*loadLayerOut['initial'][0:2],4*[0.01]))
stringers['lower_bound'] = np.array(6*[0.002])
stringers['upper_bound'] = np.array(6*[0.01])
baffles['lower_bound'] = np.hstack((lbPerc*baffles['initial'][0:4],6*[0.0074]))
baffles['upper_bound'] = np.hstack((ubPerc*baffles['initial'][0:4],6*[0.0359]))

# Set descriptions for design variables
wall['description'] = ['c6x', 'c7x', 'c8x', 'c9x', 'c11x', 'c12x', 'c13x', 'c14x', 'c15x', 'c16x', 'c17x', 'c6y', 'c7y', 'c8y', 'c9y', 'c12y', 'c13y', 'c14y', 'c15y', 'c16y', 'c17y']
thermalLayer['description'] = ['thermal_layer_x2', 'thermal_layer_x3', 'thermal_layer_t1', 'thermal_layer_t2', 'thermal_layer_t3', 'thermal_layer_t4']
airGap['description'] = ['air_gap_thickness']
loadLayerIn['description'] = ['load_layer_in_x2', 'load_layer_in_x3', 'load_layer_in_t1', 'load_layer_in_t2', 'load_layer_in_t3', 'load_layer_in_t4']
loadLayerMid['description'] = ['load_layer_mid_x2', 'load_layer_mid_x3', 'load_layer_mid_t1', 'load_layer_mid_t2', 'load_layer_mid_t3', 'load_layer_mid_t4']
loadLayerOut['description'] = ['load_layer_out_x2', 'load_layer_out_x3', 'load_layer_out_t1', 'load_layer_out_t2', 'load_layer_out_t3', 'load_layer_out_t4']
stringers['description'] = ['stringer_t1', 'stringer_t2', 'stringer_t3', 'stringer_t4', 'stringer_t5', 'stringer_t6']
baffles['description'] = ['baffle2_x', 'baffle3_x', 'baffle4_x', 'baffle5_x', 'baffle1_t', 'baffle2_t', 'baffle3_t', 'baffle4_t', 'baffle5_t', 'baffle6_t']

# Concatenate all design variables
designVar = dict()
componentList = [wall, thermalLayer, airGap, loadLayerIn, loadLayerMid, loadLayerOut, stringers, baffles]
designVar['initial'] = concatenateDesignVariables('initial', componentList)
designVar['lower_bound'] = concatenateDesignVariables('lower_bound', componentList)
designVar['upper_bound'] = concatenateDesignVariables('upper_bound', componentList)
designVar['description'] = concatenateDesignVariables('description', componentList)

# Useful indices for generating constraints, etc.
n = [len(wall['initial'])] # number of variables for each component
n.append(len(thermalLayer['initial']))
n.append(len(airGap['initial']))
n.append(len(loadLayerIn['initial']))
n.append(len(loadLayerMid['initial']))
n.append(len(loadLayerOut['initial']))
n.append(len(stringers['initial']))
n.append(len(baffles['initial']))
nDesignVar = len(designVar['initial'])
index = [0] # starting index for each component's constraints
index.append(index[-1]+n[0])
index.append(index[-1]+n[1])
index.append(index[-1]+n[2])
index.append(index[-1]+n[3])
index.append(index[-1]+n[4])
index.append(index[-1]+n[5])
index.append(index[-1]+n[6])

## ---- Define uncertain parameters

## ---- Define nonlinear inequality constraints, g(x) <= d
conNLI = dict()
conNLI['lower_bound'] = np.array([21500] + 12*[0])
conNLI['upper_bound'] = np.array([100000] + 12*[1])
conNLI['description'] = ['thrust', 'load_layer_out_temp', 'thermal_layer_stress', 'load_layer_in_stress', 'load_layer_mid_stress', 'load_layer_out_stress', 'stringers_stress', 'baffle1_stress', 'baffle2_stress', 'baffle3_stress', 'baffle4_stress', 'baffle5_stress', 'baffle6_stress']

## ---- Define linear inequality constraints, Ax <= b

(tmp1, A_w, b_w) = linearConstraints.wall(wall['initial'])
(tmp2, A_tl, b_tl) = linearConstraints.thermalLayer(thermalLayer['initial'])
# no linear constraints for air gap
(tmp3, A_lli, b_lli) = linearConstraints.thermalLayer(loadLayerIn['initial']) # setup is the same
(tmp4, A_llm, b_llm) = linearConstraints.thermalLayer(loadLayerMid['initial']) # setup is the same 
(tmp5, A_llo, b_llo) = linearConstraints.thermalLayer(loadLayerOut['initial']) # setup is the same
# no linear constraints for stringers since only thickness is used
(tmp6, A_b, b_b) = linearConstraints.baffles(baffles['initial']) # setup is the same

# Assemble the linear inequality constraints
nLinearCon = 55 + 7 + 7*3 + 2 + 12
A = np.zeros((nLinearCon,nDesignVar)) # matrix for Ax <= b
b = np.zeros(nLinearCon) # RHS vector for Ax <= b
A[0:59,index[0]:index[0]+n[0]] = A_w; b[0:59] = b_w
A[59:59+7,index[1]:index[1]+n[1]] = A_tl; b[59:59+7] = b_tl
A[59+7:59+7*2,index[3]:index[3]+n[3]] = A_lli; b[59+7:59+7*2] = b_lli
A[59+7*2:59+7*3,index[4]:index[4]+n[4]] = A_llm; b[59+7*2:59+7*3] = b_llm
A[59+7*3:59+7*4,index[5]:index[5]+n[5]] = A_llo; b[59+7*3:59+7*4] = b_llo
A[59+7*4:59+7*4+10,index[7]:index[7]+n[7]] = A_b; 
b[59+7*4:59+7*4+10] = b_b

# Scale the linear constraints to unity on RHS
for row in range(nLinearCon):
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
f.write("# Deterministic optimization of small climb problem\n")
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
f.write("        scales = 4.6512e-5 1 1 1 1 1 1 1 1 1 1 1 1\n")
f.write("    numerical_gradients\n")
f.write("        method_source dakota\n")
f.write("        interval_type forward\n")
f.write("        fd_gradient_step_size = 1.e-3\n")
f.write("    no_hessians\n")
f.write("\n")
f.flush()

# write interface
f.write("interface\n")
f.write("    id_interface = 'OPTIM_I'\n")
f.write("    analysis_driver = 'python runModel.py -l 0 -f det_optim_standard.cfg'\n")
f.write("        fork asynchronous evaluation_concurrency = 4\n")
f.write("        work_directory named 'tmp'\n")
#f.write("        link_files '*'\n")
f.write("        link_files 'runModel.py' 'det_optim_standard.cfg'\n")
f.write("        directory_tag\n")
f.write("        directory_save\n")
#f.write("        file_tag\n")
f.write("        file_save\n")
f.write("        parameters_file = 'params.in'\n")
f.write("        results_file = 'results.out'\n")
#f.write("    failure_capture\n")
#f.write("        recover 1.e3 0.e0\n")

f.close()

