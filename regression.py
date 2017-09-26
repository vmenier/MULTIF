"""
Run a series of regression tests for MULTI-F. Tests are broken down into
the following sets:

*Test everything
*Test low-fidelity analysis only (2D and 3D parameterizations)
*Test medium-fidelity analysis only (2D and 3D parameterizations)
*Test medium-fidelity adjoint cability (2D parameterization only)
*Test finite difference gradient capability
*Test linear and nonlinear FEA for 2D parameterization

Rick Fenrich 9/2/17
"""

import os
from optparse import OptionParser

from TestCase import TestCase

# =========================================================================== #
# Setup options etc.
# =========================================================================== #

# Set up options for regression testing
parser = OptionParser();
parser.add_option("-n", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS");
(options, args)=parser.parse_args();

nCores = int(options.partitions);

testNum = 1;

# Control which tests are run
test = {'all': 0, 'lofi': 0, 'medfi': 1, 'adjoint': 0, 'fd': 0, 'fea': 0};

# Set up necessary filepaths
#rootdir = os.getcwd()
rootdir = os.path.dirname(os.path.abspath(__file__));

# Start fresh log file by overwriting previous one
f = open('regression.out','w');
f.close();

# Add MULTI-F to PYTHONPATH environment variable
#os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + rootdir;
#print os.environ['PYTHONPATH'];
# DO THE ABOVE ONLY IF THIS IS NOT ALREADY DONE

test_list = [];

# =========================================================================== #
# Begin tests with low-fidelity analyses
# =========================================================================== #

# --------------------------- 2D param, lo-fi --------------------------------
if( test['all'] or test['lofi'] or test['fea'] ):

    lofi_2dparam_1 = TestCase('lofi_2dparam_1');
    lofi_2dparam_1.description = 'General, 2D, serial, low-fi analysis';
    lofi_2dparam_1.cfg_dir = 'example';
    lofi_2dparam_1.cfg_file = 'general.cfg';
    lofi_2dparam_1.input_file = 'general.in';
    lofi_2dparam_1.compare_file = 'example/regression/results_general_lofi_1.out';
    lofi_2dparam_1.fidelity = 0;
    lofi_2dparam_1.nproc = 1;
    lofi_2dparam_1.diff_tol = 1e-6;
    test_list.append(lofi_2dparam_1);

# --------------------------- 3D param, lo-fi --------------------------------
#if( test['all'] or test['lofi'] ):
#
#    lofi_3dparam = TestCase('lofi_3dparam');
#    lofi_3dparam.description = 'General, 3D, serial, low-fi analysis';
#    lofi_3dparam.cfg_dir = 'example';
#    lofi_3dparam.cfg_file = 'general-3d.cfg';
#    lofi_3dparam.input_file = 'general-3d.in';
#    lofi_3dparam.compare_file = 'example/regression/results_general-3d_lofi.out';
#    lofi_3dparam.fidelity = 0;
#    lofi_3dparam.nproc = 1;
#    lofi_3dparam.diff_tol = 1e-6;
#    #test_list.append(lofi_3dparam);

# ------------------- 2D param, lo-fi, f.d. gradients ------------------------
if( test['all'] or test['lofi'] or test['fd'] ):

    lofi_2dparam_fd = TestCase('lofi_2dparam_fd');
    lofi_2dparam_fd.description = 'Mass and thrust, 2D, parallel, low-fi analysis with finite difference gradients';
    lofi_2dparam_fd.cfg_dir = os.path.join('example','gradients');
    lofi_2dparam_fd.cfg_file = 'general_fd_gradients.cfg';
    lofi_2dparam_fd.input_file = 'params.in';
    lofi_2dparam_fd.compare_file = 'example/regression/results_general_fd_gradients_lofi.out';
    lofi_2dparam_fd.fidelity = 0;
    lofi_2dparam_fd.nproc = nCores;
    lofi_2dparam_fd.diff_tol = 1e-6;
    test_list.append(lofi_2dparam_fd);

# ----------------- 2D param, lo-fi, wall temp as input ----------------------
if( test['all'] or test['lofi'] ):

    lofi_2dparam_temp = TestCase('lofi_2dparam_temp');
    lofi_2dparam_temp.description = 'Temperature input, 2D, serial, low-fi analysis';
    lofi_2dparam_temp.cfg_dir = os.path.join('example','inference','infer-temp');
    lofi_2dparam_temp.cfg_file = 'inference_standard.cfg';
    lofi_2dparam_temp.input_file = 'inference_standard.in';
    lofi_2dparam_temp.dependencies = ['pressure_locations.in',
                                      'velocity_locations.in',
                                      'wall_pressure_locations.in',
                                      'wall_temperature.in',
                                      'wall_temperature_locations-100.in'];
    lofi_2dparam_temp.compare_file = 'example/regression/results_temp_input_lofi.out';
    lofi_2dparam_temp.fidelity = 0;
    lofi_2dparam_temp.nproc = 1;
    lofi_2dparam_temp.diff_tol = 1e-6;
    test_list.append(lofi_2dparam_temp);

# ------------------- 2D param, lo-fi, nonlinear FEA -------------------------
if( test['all'] or test['lofi'] or test['fea'] ):

    lofi_2dparam_2 = TestCase('lofi_2dparam_2');
    lofi_2dparam_2.description = 'General, 2D, serial, low-fi analysis with nonlinear FEA';
    lofi_2dparam_2.cfg_dir = 'example';
    lofi_2dparam_2.cfg_file = 'general.cfg';
    lofi_2dparam_2.input_file = 'general.in';
    lofi_2dparam_2.compare_file = 'example/regression/results_general_lofi_2.out';
    lofi_2dparam_2.fidelity = 1;
    lofi_2dparam_2.nproc = 1;
    lofi_2dparam_2.diff_tol = 1e-6;
    test_list.append(lofi_2dparam_2);

# =========================================================================== #
# Perform medium-fidelity analyses
# =========================================================================== #

# --------------------------- 2D param, med-fi -------------------------------
if( test['all'] or test['medfi'] ):

    medfi_2dparam = TestCase('medfi_2dparam');
    medfi_2dparam.description = 'General, 2D, parallel, med-fi analysis';
    medfi_2dparam.cfg_dir = 'example';
    medfi_2dparam.cfg_file = 'general.cfg';
    medfi_2dparam.input_file = 'general.in';
    medfi_2dparam.compare_file = 'example/regression/results_general_medfi_1.out';
    medfi_2dparam.fidelity = 2;
    medfi_2dparam.nproc = nCores;
    medfi_2dparam.diff_tol = 1e-6;
    test_list.append(medfi_2dparam);

# ----------------------- 3D param, dim 2D, med-fi ---------------------------
if( test['all'] or test['medfi'] ):

    medfi_3dparam = TestCase('medfi_3dparam');
    medfi_3dparam.description = 'General, 3D, parallel, med-fi analysis';
    medfi_3dparam.cfg_dir = 'example';
    medfi_3dparam.cfg_file = 'general-3d.cfg';
    medfi_3dparam.input_file = 'general-3d.in';
    medfi_3dparam.compare_file = 'example/regression/results_general-3d_medfi_1.out';
    medfi_3dparam.fidelity = 2;
    medfi_3dparam.nproc = nCores;
    medfi_3dparam.diff_tol = 1e-6;
    test_list.append(medfi_3dparam);

# ------------------ 2D param, med-fi, adjoint gradients ---------------------
if( test['all'] or test['medfi'] or test['adjoint'] ):

    medfi_2dparam_adjoint = TestCase('medfi_2dparam_adjoint');
    medfi_2dparam_adjoint.description = 'General, 2D, parallel, med-fi analysis with adjoint gradients';
    medfi_2dparam_adjoint.cfg_dir = os.path.join('example','gradients');
    medfi_2dparam_adjoint.cfg_file = 'general_gradients.cfg';
    medfi_2dparam_adjoint.input_file = 'params.in';
    medfi_2dparam_adjoint.compare_file = 'example/regression/results_general_gradients_medfi.out';
    medfi_2dparam_adjoint.fidelity = 0;
    medfi_2dparam_adjoint.nproc = nCores;
    medfi_2dparam_adjoint.diff_tol = 1e-6;
    test_list.append(medfi_2dparam_adjoint);

# --------- 2D param, med-fi, f.d. gradients w/ mesh deformation -------------
if( test['all'] or test['medfi'] or test['fd'] ):
   medfi_2dparam_fd = TestCase('medfi_2dparam_fd');
   medfi_2dparam_fd.description = 'Mass and thrust, 2D, parallel, med-fi analysis with finite #difference gradients using mesh deformation';
   medfi_2dparam_fd.cfg_dir = os.path.join('example','gradients');
   medfi_2dparam_fd.cfg_file = 'general_fd_gradients.cfg';
   medfi_2dparam_fd.input_file = 'params.in';
   medfi_2dparam_fd.compare_file = 'example/regression/results_general_fd_gradients_medfi.out';
   medfi_2dparam_fd.fidelity = 1;
   medfi_2dparam_fd.nproc = nCores;
   medfi_2dparam_fd.diff_tol = 1e-6;
   test_list.append(medfi_2dparam_fd);

# Remaining tests get placed below as they come online

pass_list = [ t.run_test() for t in test_list ];
print '\n==================================================\n';
print 'Summary of the tests'
for i, test in enumerate(test_list):
    if( pass_list[i] == 1 ):
        print '  passed - %s' % test.name;
    elif( pass_list[i] == 2 ):
        print '  DIFFED - %s' % test.name;
    else:
        print '* FAILED - %s' % test.name;
