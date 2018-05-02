"""
Run a series of regression tests for MULTI-F's 3D nozzle parameterization.
Tests are broken down into the following sets:

*Test everything
*Test low-fidelity analysis only
*Test medium-fidelity analysis only
*Test high-fidelity analysis only
*Test finite-difference gradient code
*Test nonlinear FEA

Rick Fenrich 5/1/17
"""

import os
from optparse import OptionParser

from TestCase import TestCase

# =========================================================================== #
# Setup options etc.
# =========================================================================== #

# Set up options for regression testing
parser = OptionParser()
parser.add_option("-n", "--ntasks", dest="nTasks", default=1,
                    help="number of tasks", metavar="NTASKS")    
parser.add_option("-c", "--cpus-per-task", dest="cpusPerTask",
                    default=1, help="cpus requested per task",
                    metavar="CPUS_PER_TASK")
(options, args)=parser.parse_args()

nTasks = int( options.nTasks )
cpusPerTask = int( options.cpusPerTask )

testNum = 1

# Control which tests are run
test = {'all': 1, 'lofi': 0, 'medfi': 0, 'hifi': 0, 'euler': 0, 'rans': 0, 
        'gradient': 0, 'inference': 0, 'fea': 0}

# Set up necessary filepaths
#rootdir = os.getcwd()
rootdir = os.path.dirname(os.path.abspath(__file__))

# Start fresh log file by overwriting previous one
f = open('regression.out','w')
f.close()

test_list = []

# =========================================================================== #
# Low-fidelity analysis tests
# =========================================================================== #

# --------------------------- 3D param, 1D dim lo-fi --------------------------
if( test['all'] or test['lofi'] ):

    lofi_3dparam = TestCase('lofi_3dparam')
    lofi_3dparam.description = 'General, 3D->2D, serial, low-fi analysis'
    lofi_3dparam.cfg_dir = 'example'
    lofi_3dparam.cfg_file = 'general-3d.cfg'
    lofi_3dparam.input_file = 'general-3d.in'
    lofi_3dparam.compare_file = 'example/regression/lofi_3dparam.out'
    lofi_3dparam.fidelity = 0
    lofi_3dparam.ntasks = 1
    lofi_3dparam.cpus_per_task = 1
    lofi_3dparam.diff_tol = 1e-6
    test_list.append(lofi_3dparam)

# =========================================================================== #
# Medium-fidelity analysis tests (2D axisymmetric downscaled from 3D)
# =========================================================================== #

# ----------------------- 3D param, 2D dim, med-fi ----------------------------
if( test['all'] or test['medfi'] or test['euler'] ):

    medfi_3dparam = TestCase('medfi_3dparam')
    medfi_3dparam.description = 'General, 3D->2D, parallel, med-fi (Euler) analysis, coarse'
    medfi_3dparam.cfg_dir = 'example'
    medfi_3dparam.cfg_file = 'general-3d.cfg'
    medfi_3dparam.input_file = 'general-3d.in'
    medfi_3dparam.compare_file = 'example/regression/medfi_3dparam.out'
    medfi_3dparam.fidelity = 2
    medfi_3dparam.ntasks = nTasks
    medfi_3dparam.cpus_per_task = cpusPerTask
    medfi_3dparam.diff_tol = 1e-6
    test_list.append(medfi_3dparam)

# ----------------------- 3D param, 2D dim, med-fi ----------------------------
if( test['all'] or test['medfi'] or test['euler'] ):

    medfi_3dparam_med = TestCase('medfi_3dparam_med')
    medfi_3dparam_med.description = 'General, 3D->2D, parallel, med-fi (Euler) analysis, medium'
    medfi_3dparam_med.cfg_dir = 'example'
    medfi_3dparam_med.cfg_file = 'general-3d.cfg'
    medfi_3dparam_med.input_file = 'general-3d.in'
    medfi_3dparam_med.compare_file = 'example/regression/medfi_3dparam_med.out'
    medfi_3dparam_med.fidelity = 3
    medfi_3dparam_med.ntasks = nTasks
    medfi_3dparam_med.cpus_per_task = cpusPerTask
    medfi_3dparam_med.diff_tol = 1e-6
    test_list.append(medfi_3dparam_med)

# ----------------------- 3D param, 2D dim, med-fi ----------------------------
if( test['all'] or test['medfi'] or test['euler'] ):

    medfi_3dparam_fine = TestCase('medfi_3dparam_fine')
    medfi_3dparam_fine.description = 'General, 3D->2D, parallel, med-fi (Euler) analysis, fine'
    medfi_3dparam_fine.cfg_dir = 'example'
    medfi_3dparam_fine.cfg_file = 'general-3d.cfg'
    medfi_3dparam_fine.input_file = 'general-3d.in'
    medfi_3dparam_fine.compare_file = 'example/regression/medfi_3dparam_fine.out'
    medfi_3dparam_fine.fidelity = 4
    medfi_3dparam_fine.ntasks = nTasks
    medfi_3dparam_fine.cpus_per_task = cpusPerTask
    medfi_3dparam_fine.diff_tol = 1e-6
    test_list.append(medfi_3dparam_fine)

# ----------------------- 3D param, 2D dim, med-fi ----------------------------
if( test['all'] or test['medfi'] or test['rans'] ):

    medfi_3dparam_rans = TestCase('medfi_3dparam_rans')
    medfi_3dparam_rans.description = 'General, 3D->2D, parallel, med-fi (RANS) analysis, coarse'
    medfi_3dparam_rans.cfg_dir = 'example'
    medfi_3dparam_rans.cfg_file = 'general-3d.cfg'
    medfi_3dparam_rans.input_file = 'general-3d.in'
    medfi_3dparam_rans.compare_file = 'example/regression/medfi_3dparam_rans.out'
    medfi_3dparam_rans.fidelity = 8
    medfi_3dparam_rans.ntasks = nTasks
    medfi_3dparam_rans.cpus_per_task = cpusPerTask
    medfi_3dparam_rans.diff_tol = 1e-6
    test_list.append(medfi_3dparam_rans)

# ----------------------- 3D param, 2D dim, med-fi ----------------------------
if( test['all'] or test['medfi'] or test['rans'] ):

    medfi_3dparam_rans_med = TestCase('medfi_3dparam_rans_med')
    medfi_3dparam_rans_med.description = 'General, 3D->2D, parallel, med-fi (RANS) analysis, medium'
    medfi_3dparam_rans_med.cfg_dir = 'example'
    medfi_3dparam_rans_med.cfg_file = 'general-3d.cfg'
    medfi_3dparam_rans_med.input_file = 'general-3d.in'
    medfi_3dparam_rans_med.compare_file = 'example/regression/medfi_3dparam_rans_med.out'
    medfi_3dparam_rans_med.fidelity = 9
    medfi_3dparam_rans_med.ntasks = nTasks
    medfi_3dparam_rans_med.cpus_per_task = cpusPerTask
    medfi_3dparam_rans_med.diff_tol = 1e-6
    test_list.append(medfi_3dparam_rans_med)    

# ----------------------- 3D param, 2D dim, med-fi ----------------------------
if( test['all'] or test['medfi'] or test['rans'] ):

    medfi_3dparam_rans_fine = TestCase('medfi_3dparam_rans_fine')
    medfi_3dparam_rans_fine.description = 'General, 3D->2D, parallel, med-fi (RANS) analysis, fine'
    medfi_3dparam_rans_fine.cfg_dir = 'example'
    medfi_3dparam_rans_fine.cfg_file = 'general-3d.cfg'
    medfi_3dparam_rans_fine.input_file = 'general-3d.in'
    medfi_3dparam_rans_fine.compare_file = 'example/regression/medfi_3dparam_rans_fine.out'
    medfi_3dparam_rans_fine.fidelity = 10
    medfi_3dparam_rans_fine.ntasks = nTasks
    medfi_3dparam_rans_fine.cpus_per_task = cpusPerTask
    medfi_3dparam_rans_fine.diff_tol = 1e-6
    test_list.append(medfi_3dparam_rans_fine) 

# =========================================================================== #
# High-fidelity analysis tests
# =========================================================================== #

# ----------------------- 3D param, 3D dim, Euler ---------------------------
if( test['all'] or test['hifi'] or test['euler'] ):

    hifi_euler = TestCase('hifi_euler')
    hifi_euler.description = 'General, 3D, parallel, hi-fi (Euler) analysis, coarse'
    hifi_euler.cfg_dir = 'example'
    hifi_euler.cfg_file = 'general-3d.cfg'
    hifi_euler.input_file = 'general-3d.in'
    hifi_euler.compare_file = 'example/regression/hifi_euler.out'
    hifi_euler.fidelity = 5
    hifi_euler.ntasks = nTasks
    hifi_euler.cpus_per_task = cpusPerTask
    hifi_euler.diff_tol = 1e-6
    test_list.append(hifi_euler)

# ----------------------- 3D param, 3D dim, Euler ---------------------------
if( test['all'] or test['hifi'] or test['euler'] ):

    hifi_euler_med = TestCase('hifi_euler_med')
    hifi_euler_med.description = 'General, 3D, parallel, hi-fi (Euler) analysis, medium'
    hifi_euler_med.cfg_dir = 'example'
    hifi_euler_med.cfg_file = 'general-3d.cfg'
    hifi_euler_med.input_file = 'general-3d.in'
    hifi_euler_med.compare_file = 'example/regression/hifi_euler_med.out'
    hifi_euler_med.fidelity = 6
    hifi_euler_med.ntasks = nTasks
    hifi_euler_med.cpus_per_task = cpusPerTask
    hifi_euler_med.diff_tol = 1e-6
    test_list.append(hifi_euler_med)

# ----------------------- 3D param, 3D dim, Euler ---------------------------
if( test['all'] or test['hifi'] or test['euler'] ):

    hifi_euler_fine = TestCase('hifi_euler_fine')
    hifi_euler_fine.description = 'General, 3D, parallel, hi-fi (Euler) analysis, fine'
    hifi_euler_fine.cfg_dir = 'example'
    hifi_euler_fine.cfg_file = 'general-3d.cfg'
    hifi_euler_fine.input_file = 'general-3d.in'
    hifi_euler_fine.compare_file = 'example/regression/hifi_euler_fine.out'
    hifi_euler_fine.fidelity = 7
    hifi_euler_fine.ntasks = nTasks
    hifi_euler_fine.cpus_per_task = cpusPerTask
    hifi_euler_fine.diff_tol = 1e-6
    test_list.append(hifi_euler_fine)

# ----------------------- 3D param, 3D dim, RANS ----------------------------
if( test['all'] or test['hifi'] or test['rans'] ):

    hifi_rans = TestCase('hifi_rans')
    hifi_rans.description = 'General, 3D, parallel, hi-fi (RANS) analysis, coarse'
    hifi_rans.cfg_dir = 'example'
    hifi_rans.cfg_file = 'general-3d.cfg'
    hifi_rans.input_file = 'general-3d.in'
    hifi_rans.compare_file = 'example/regression/hifi_rans.out'
    hifi_rans.fidelity = 11
    hifi_rans.ntasks = nTasks
    hifi_rans.cpus_per_task = cpusPerTask
    hifi_rans.diff_tol = 1e-6
    test_list.append(hifi_rans)

# ----------------------- 3D param, 3D dim, RANS ----------------------------
if( test['all'] or test['hifi'] or test['rans'] ):

    hifi_rans_med = TestCase('hifi_rans_med')
    hifi_rans_med.description = 'General, 3D, parallel, hi-fi (RANS) analysis, medium'
    hifi_rans_med.cfg_dir = 'example'
    hifi_rans_med.cfg_file = 'general-3d.cfg'
    hifi_rans_med.input_file = 'general-3d.in'
    hifi_rans_med.compare_file = 'example/regression/hifi_rans_med.out'
    hifi_rans_med.fidelity = 12
    hifi_rans_med.ntasks = nTasks
    hifi_rans_med.cpus_per_task = cpusPerTask
    hifi_rans_med.diff_tol = 1e-6
    test_list.append(hifi_rans_med)

# ----------------------- 3D param, 3D dim, RANS ----------------------------
if( test['all'] or test['hifi'] or test['rans'] ):

    hifi_rans_fine = TestCase('hifi_rans_fine')
    hifi_rans_fine.description = 'General, 3D, parallel, hi-fi (RANS) analysis, fine'
    hifi_rans_fine.cfg_dir = 'example'
    hifi_rans_fine.cfg_file = 'general-3d.cfg'
    hifi_rans_fine.input_file = 'general-3d.in'
    hifi_rans_fine.compare_file = 'example/regression/hifi_rans_fine.out'
    hifi_rans_fine.fidelity = 13
    hifi_rans_fine.ntasks = nTasks
    hifi_rans_fine.cpus_per_task = cpusPerTask
    hifi_rans_fine.diff_tol = 1e-6
    test_list.append(hifi_rans_fine)

# =========================================================================== #
# SST eigen perturbation analysis tests
# =========================================================================== #

# ------------ 3D param, 2D dim, med-fi with SST perturbation -----------------
if( test['all'] or test['medfi'] or test['rans'] ):

    medfi_3dparam_sst = TestCase('medfi_3dparam_sst')
    medfi_3dparam_sst.description = 'General, 3D->2D, parallel, med-fi (RANS) analysis with SST perturbation, coarse'
    medfi_3dparam_sst.cfg_dir = 'example'
    medfi_3dparam_sst.cfg_file = 'general-3d_sst.cfg'
    medfi_3dparam_sst.input_file = 'general.in'
    medfi_3dparam_sst.compare_file = 'example/regression/medfi_3dparam_sst.out'
    medfi_3dparam_sst.fidelity = 0
    medfi_3dparam_sst.ntasks = nTasks
    medfi_3dparam_sst.cpus_per_task = cpusPerTask
    medfi_3dparam_sst.diff_tol = 1e-6
    test_list.append(medfi_3dparam_sst)

# ------------ 3D param, 2D dim, med-fi with SST perturbation -----------------
if( test['all'] or test['medfi'] or test['rans'] ):

    medfi_3dparam_sst_med = TestCase('medfi_3dparam_sst_med')
    medfi_3dparam_sst_med.description = 'General, 3D->2D, parallel, med-fi (RANS) analysis with SST perturbation, medium'
    medfi_3dparam_sst_med.cfg_dir = 'example'
    medfi_3dparam_sst_med.cfg_file = 'general-3d_sst.cfg'
    medfi_3dparam_sst_med.input_file = 'general.in'
    medfi_3dparam_sst_med.compare_file = 'example/regression/medfi_3dparam_sst_med.out'
    medfi_3dparam_sst_med.fidelity = 1
    medfi_3dparam_sst_med.ntasks = nTasks
    medfi_3dparam_sst_med.cpus_per_task = cpusPerTask
    medfi_3dparam_sst_med.diff_tol = 1e-6
    test_list.append(medfi_3dparam_sst_med)

# ------------ 3D param, 2D dim, med-fi with SST perturbation -----------------
if( test['all'] or test['medfi'] or test['rans'] ):

    medfi_3dparam_sst_fine = TestCase('medfi_3dparam_sst_fine')
    medfi_3dparam_sst_fine.description = 'General, 3D->2D, parallel, med-fi (RANS) analysis with SST perturbation, fine'
    medfi_3dparam_sst_fine.cfg_dir = 'example'
    medfi_3dparam_sst_fine.cfg_file = 'general-3d_sst.cfg'
    medfi_3dparam_sst_fine.input_file = 'general.in'
    medfi_3dparam_sst_fine.compare_file = 'example/regression/medfi_3dparam_sst_fine.out'
    medfi_3dparam_sst_fine.fidelity = 2
    medfi_3dparam_sst_fine.ntasks = nTasks
    medfi_3dparam_sst_fine.cpus_per_task = cpusPerTask
    medfi_3dparam_sst_fine.diff_tol = 1e-6
    test_list.append(medfi_3dparam_sst_fine)

# ------------ 3D param, 3D dim, hi-fi with SST perturbation -----------------
if( test['all'] or test['hifi'] or test['rans'] ):

    hifi_3dparam_sst = TestCase('hifi_3dparam_sst')
    hifi_3dparam_sst.description = 'General, 3D, parallel, hi-fi (RANS) analysis with SST perturbation, coarse'
    hifi_3dparam_sst.cfg_dir = 'example'
    hifi_3dparam_sst.cfg_file = 'general-3d_sst.cfg'
    hifi_3dparam_sst.input_file = 'general.in'
    hifi_3dparam_sst.compare_file = 'example/regression/hifi_3dparam_sst.out'
    hifi_3dparam_sst.fidelity = 3
    hifi_3dparam_sst.ntasks = nTasks
    hifi_3dparam_sst.cpus_per_task = cpusPerTask
    hifi_3dparam_sst.diff_tol = 1e-6
    test_list.append(hifi_3dparam_sst)

# ------------ 3D param, 3D dim, hi-fi with SST perturbation -----------------
if( test['all'] or test['hifi'] or test['rans'] ):

    hifi_3dparam_sst_med = TestCase('hifi_3dparam_sst_med')
    hifi_3dparam_sst_med.description = 'General, 3D, parallel, hi-fi (RANS) analysis with SST perturbation, med'
    hifi_3dparam_sst_med.cfg_dir = 'example'
    hifi_3dparam_sst_med.cfg_file = 'general-3d_sst.cfg'
    hifi_3dparam_sst_med.input_file = 'general.in'
    hifi_3dparam_sst_med.compare_file = 'example/regression/hifi_3dparam_sst_med.out'
    hifi_3dparam_sst_med.fidelity = 4
    hifi_3dparam_sst_med.ntasks = nTasks
    hifi_3dparam_sst_med.cpus_per_task = cpusPerTask
    hifi_3dparam_sst_med.diff_tol = 1e-6
    test_list.append(hifi_3dparam_sst_med)

# ------------ 3D param, 3D dim, hi-fi with SST perturbation -----------------
if( test['all'] or test['hifi'] or test['rans'] ):

    hifi_3dparam_sst_fine = TestCase('hifi_3dparam_sst_fine')
    hifi_3dparam_sst_fine.description = 'General, 3D, parallel, hi-fi (RANS) analysis with SST perturbation, fine'
    hifi_3dparam_sst_fine.cfg_dir = 'example'
    hifi_3dparam_sst_fine.cfg_file = 'general-3d_sst.cfg'
    hifi_3dparam_sst_fine.input_file = 'general.in'
    hifi_3dparam_sst_fine.compare_file = 'example/regression/hifi_3dparam_sst_fine.out'
    hifi_3dparam_sst_fine.fidelity = 5
    hifi_3dparam_sst_fine.ntasks = nTasks
    hifi_3dparam_sst_fine.cpus_per_task = cpusPerTask
    hifi_3dparam_sst_fine.diff_tol = 1e-6
    test_list.append(hifi_3dparam_sst_fine)

# =========================================================================== #
# Gradient calculation tests
# =========================================================================== #

# # ------------------- 2D param, 1D dim lo-fi, f.d. gradients ------------------
# if( test['all'] or test['gradient'] ):

#     lofi_3dparam_fd = TestCase('lofi_3dparam_fd')
#     lofi_3dparam_fd.description = 'Mass and thrust, 2D, parallel, low-fi analysis with finite difference gradients'
#     lofi_3dparam_fd.cfg_dir = os.path.join('example','gradients')
#     lofi_3dparam_fd.cfg_file = 'general_fd_gradients.cfg'
#     lofi_3dparam_fd.input_file = 'params.in'
#     lofi_3dparam_fd.compare_file = 'example/regression/lofi_3dparam_fd.out'
#     lofi_3dparam_fd.fidelity = 0
#     lofi_3dparam_fd.ntasks = nTasks
#     lofi_3dparam_fd.cpus_per_task = cpusPerTask
#     lofi_3dparam_fd.diff_tol = 1e-6
#     test_list.append(lofi_3dparam_fd)

# # ------------------ 2D param, 2D dim med-fi, adjoint gradients ---------------
# if( test['all'] or test['medfi'] or test['adjoint'] ):

#     medfi_2dparam_adjoint = TestCase('medfi_2dparam_adjoint')
#     medfi_2dparam_adjoint.description = 'General, 2D, parallel, med-fi (Euler) analysis with adjoint gradients'
#     medfi_2dparam_adjoint.cfg_dir = os.path.join('example','gradients')
#     medfi_2dparam_adjoint.cfg_file = 'general_gradients.cfg'
#     medfi_2dparam_adjoint.input_file = 'params.in'
#     medfi_2dparam_adjoint.compare_file = 'example/regression/medfi_2dparam_adjoint.out'
#     medfi_2dparam_adjoint.fidelity = 0
#     medfi_2dparam_adjoint.ntasks = nTasks
#     medfi_2dparam_adjoint.cpus_per_task = cpusPerTask
#     medfi_2dparam_adjoint.diff_tol = 1e-6
#     test_list.append(medfi_2dparam_adjoint)

# # --------- 2D param, 2D dim, med-fi, f.d. gradients w/ mesh deformation ------
# if( test['all'] or test['medfi'] or test['fd'] ):
#    medfi_2dparam_fd = TestCase('medfi_2dparam_fd')
#    medfi_2dparam_fd.description = 'Mass and thrust, 2D, parallel, med-fi (Euler) analysis with finite difference gradients using mesh deformation'
#    medfi_2dparam_fd.cfg_dir = os.path.join('example','gradients')
#    medfi_2dparam_fd.cfg_file = 'general_fd_gradients.cfg'
#    medfi_2dparam_fd.input_file = 'params.in'
#    medfi_2dparam_fd.compare_file = 'example/regression/medfi_2dparam_fd.out'
#    medfi_2dparam_fd.fidelity = 1
#    medfi_2dparam_fd.ntasks = nTasks
#    medfi_2dparam_fd.cpus_per_task = cpusPerTask
#    medfi_2dparam_fd.diff_tol = 1e-6
#    test_list.append(medfi_2dparam_fd)

# =========================================================================== #
# Inference tests
# =========================================================================== #

# # ----------------- 2D param, 1D dim lo-fi, wall temp as input ----------------
# if( test['all'] or test['inference'] ):

#     lofi_2dparam_temp = TestCase('lofi_2dparam_temp')
#     lofi_2dparam_temp.description = 'Temperature input, 2D, serial, low-fi analysis'
#     lofi_2dparam_temp.cfg_dir = os.path.join('example','inference','infer-temp')
#     lofi_2dparam_temp.cfg_file = 'inference_standard.cfg'
#     lofi_2dparam_temp.input_file = 'inference_standard.in'
#     lofi_2dparam_temp.dependencies = ['pressure_locations.in',
#                                       'velocity_locations.in',
#                                       'wall_pressure_locations.in',
#                                       'wall_temperature.in',
#                                       'wall_temperature_locations-100.in']
#     lofi_2dparam_temp.compare_file = 'example/regression/lofi_2dparam_temp.out'
#     lofi_2dparam_temp.fidelity = 0
#     lofi_2dparam_temp.ntasks = 1
#     lofi_2dparam_temp.cpus_per_task = 1
#     lofi_2dparam_temp.diff_tol = 1e-6
#     test_list.append(lofi_2dparam_temp)

# =========================================================================== #
# Begin nonlinear FEA tests
# =========================================================================== #

# # ------------------- 2D param, 1D dim lo-fi, nonlinear FEA -------------------
# if( test['all'] or test['fea'] ):

#     lofi_2dparam_2 = TestCase('lofi_2dparam_2')
#     lofi_2dparam_2.description = 'General, 2D, serial, low-fi analysis with nonlinear FEA'
#     lofi_2dparam_2.cfg_dir = 'example'
#     lofi_2dparam_2.cfg_file = 'general.cfg'
#     lofi_2dparam_2.input_file = 'general.in'
#     lofi_2dparam_2.compare_file = 'example/regression/lofi_2dparam_2.out'
#     lofi_2dparam_2.fidelity = 1
#     lofi_2dparam_2.ntasks = 1
#     lofi_2dparam_2.cpus_per_task = 1
#     lofi_2dparam_2.diff_tol = 1e-6
#     test_list.append(lofi_2dparam_2)

# Remaining tests get placed below as they come online

pass_list = [ t.run_test() for t in test_list ]
print '\n==================================================\n'
print 'Summary of the tests'
for i, test in enumerate(test_list):
    if( pass_list[i] == 1 ):
        print '  passed - %s' % test.name
    elif( pass_list[i] == 2 ):
        print '  DIFFED - %s' % test.name
    else:
        print '* FAILED - %s' % test.name
