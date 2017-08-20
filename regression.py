"""
Run a series of regression tests for MULTI-F

Rick Fenrich 8/7/17
"""

import os, sys
import shutil
from subprocess import call
from optparse import OptionParser
import datetime

def compareResponses(absfilename,resultsFile1,resultsFile2):

    f = open(absfilename,'a');
    f.write('Response file 1: %s\n' % resultsFile1);
    f.write('Response file 2: %s\n' % resultsFile2);
    f.write('Responses will be compared here.\n\n');
    f.close();

    return 0;

# Set up options for regression testing
parser = OptionParser();
parser.add_option("-n", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS");
(options, args)=parser.parse_args();

nCores = int(options.partitions);

testNum = 1;

# Set up necessary filepaths

rootdir = os.getcwd()
print 'MULTI-F root directory: %s' % rootdir;

testdir = os.path.join('local','regression_tests');

filename = 'regression.out';
absfilename = os.path.join(rootdir,testdir,filename);

# Add MULTI-F to PYTHONPATH environment variable
#os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + rootdir;
#print os.environ['PYTHONPATH'];
# DO THE ABOVE ONLY IF THIS IS NOT ALREADY DONE

# Create and clean directory for tests
if not os.path.isdir(testdir):
    print 'Making directory: %s' % testdir;
    os.mkdir(testdir);
else:
    print 'Removing existing directory: %s' % testdir;
    shutil.rmtree(testdir);
    print 'Making directory: %s' % testdir;
    os.mkdir(testdir);
    
os.chdir(testdir);

# Setup file for recording test information
f = open(absfilename,'w');
f.write('%s\n' % str(datetime.datetime.now()));
f.write('Regression directory: %s\n\n' % os.path.join(rootdir,testdir));
f.close();

# =========================================================================== #
# Build file structure for tests
# =========================================================================== #
categories = ['general','general_3d','gradients','temp_input'];

for subdir in categories:
    os.mkdir(subdir);
    os.chdir(subdir);
    os.mkdir('low_fi');
    if(subdir != 'temp_input'):
        os.mkdir('med_fi');
    if(subdir == 'general'):
        os.mkdir('low_fi2');
    os.chdir('..');

# =========================================================================== #
# Begin tests with low-fidelity analyses
# =========================================================================== #

f = open(absfilename,'a');
f.write('\nBeginning low-fidelity analyses\n\n');
f.close();

# --------------------------- 2D param, lo-fi --------------------------------

os.chdir(os.path.join(rootdir,testdir,'general'));
os.chdir('low_fi');

f = open(absfilename,'a');
f.write('\nTEST %i\n' % testNum);
f.write('General, 2D, serial low-fi analysis\n');
f.write('Directory: %s\n\n' % os.getcwd());
f.close();
testNum += 1;

shutil.copyfile(os.path.join(rootdir,'example','general.cfg'),'general.cfg');
shutil.copyfile(os.path.join(rootdir,'example','general.in'),'general.in');
shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');

f = open(absfilename,'a');
f.write('Running...\n\n');
f.close();

call(['pythona','runModel.py','-f','general.cfg','-l','0']);

compareResponses(absfilename,'results_general_lofi_1.out','results.out');

# --------------------------- 3D param, lo-fi --------------------------------

#os.chdir(os.path.join(rootdir,testdir,'general_3d'));
#os.chdir('low_fi');
#
#f = open(absfilename,'a');
#f.write('\nTEST %i\n' % testNum);
#f.write('General, 3D, serial, low-fi analysis\n');
#f.write('Directory: %s\n\n' % os.getcwd());
#f.close();
#testNum += 1;
#
#shutil.copyfile(os.path.join(rootdir,'example','general-3d.cfg'),'general-3d.cfg');
#shutil.copyfile(os.path.join(rootdir,'example','general-3d.in'),'general-3d.in');
#shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');
#
#f = open(absfilename,'a');
#f.write('Running...\n\n');
#f.close();
#
#call(['pythona','runModel.py','-f','general-3d.cfg','-l','0']);
#
#compareResponses(absfilename,'results_general-3d_lofi_1.out','results.out');

# ------------------- 2D param, lo-fi, f.d. gradients ------------------------

os.chdir(os.path.join(rootdir,testdir,'gradients'));
os.chdir('low_fi');

f = open(absfilename,'a');
f.write('\nTEST %i\n' % testNum);
f.write('Mass and thrust, 2D, parallel, low-fi F.D. gradient analysis\n');
f.write('Directory: %s\n\n' % os.getcwd());
f.close();
testNum += 1;

shutil.copyfile(os.path.join(rootdir,'example','gradients','general_fd_gradients.cfg'),'general_fd_gradients.cfg');
shutil.copyfile(os.path.join(rootdir,'example','gradients','params.in'),'params.in');
shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');

f = open(absfilename,'a');
f.write('Running...\n\n');
f.close();

call(['pythona','runModel.py','-f','general_fd_gradients.cfg','-l','0','-n',str(nCores)]);

compareResponses(absfilename,'results_general_fd_gradients_lofi.out','results.out');

# ----------------- 2D param, lo-fi, wall temp as input ----------------------

os.chdir(os.path.join(rootdir,testdir,'temp_input'));
os.chdir('low_fi');

f = open(absfilename,'a');
f.write('\nTEST %i\n' % testNum);
f.write('Temperature input, 2D, serial, low-fi analysis\n');
f.write('Directory: %s\n\n' % os.getcwd());
f.close();
testNum += 1;

shutil.copyfile(os.path.join(rootdir,'example','inference','infer-temp','inference_standard.cfg'),'inference_standard.cfg');
shutil.copyfile(os.path.join(rootdir,'example','inference','infer-temp','inference_standard.in'),'inference_standard.in');
shutil.copyfile(os.path.join(rootdir,'example','inference','infer-temp','pressure_locations.in'),'pressure_locations.in');
shutil.copyfile(os.path.join(rootdir,'example','inference','infer-temp','velocity_locations.in'),'velocity_locations.in');
shutil.copyfile(os.path.join(rootdir,'example','inference','infer-temp','wall_pressure_locations.in'),'wall_pressure_locations.in');
shutil.copyfile(os.path.join(rootdir,'example','inference','infer-temp','wall_temperature.in'),'wall_temperature.in');
shutil.copyfile(os.path.join(rootdir,'example','inference','infer-temp','wall_temperature_locations-100.in'),'wall_temperature_locations-100.in');
shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');

f = open(absfilename,'a');
f.write('Running...\n\n');
f.close();

call(['pythona','runModel.py','-f','inference_standard.cfg','-l','0']);

compareResponses(absfilename,'results_temp_input_lofi.out','results.out');

# ------------------- 2D param, lo-fi, nonlinear FEA -------------------------

os.chdir(os.path.join(rootdir,testdir,'general'));
os.chdir('low_fi2');

f = open(absfilename,'a');
f.write('\nTEST %i\n' % testNum);
f.write('General, 2D, serial low-fi analysis with nonlinear FEA\n');
f.write('Directory: %s\n\n' % os.getcwd());
f.close();
testNum += 1;

shutil.copyfile(os.path.join(rootdir,'example','general.cfg'),'general.cfg');
shutil.copyfile(os.path.join(rootdir,'example','general.in'),'general.in');
shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');

f = open(absfilename,'a');
f.write('Running...\n\n');
f.close();

call(['pythona','runModel.py','-f','general.cfg','-l','1']);

compareResponses(absfilename,'results_general_lofi_2.out','results.out');

# =========================================================================== #
# Perform medium-fidelity analyses
# =========================================================================== #

f = open(absfilename,'a');
f.write('\nBeginning medium-fidelity analyses\n\n');
f.close();

# --------------------------- 2D param, med-fi -------------------------------

os.chdir(os.path.join(rootdir,testdir,'general'));
os.chdir('med_fi');

f = open(absfilename,'a');
f.write('\nTEST %i\n' % testNum);
f.write('General, 2D, parallel med-fi analysis\n');
f.write('Directory: %s\n\n' % os.getcwd());
f.close();
testNum += 1;

shutil.copyfile(os.path.join(rootdir,'example','general.cfg'),'general.cfg');
shutil.copyfile(os.path.join(rootdir,'example','general.in'),'general.in');
shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');

f = open(absfilename,'a');
f.write('Running...\n\n');
f.close();

call(['pythona','runModel.py','-f','general.cfg','-l','2','-n',str(nCores)]);

compareResponses(absfilename,'results_general_medfi_1.out','results.out');

# --------------------------- 3D param, med-fi -------------------------------

os.chdir(os.path.join(rootdir,testdir,'general_3d'));
os.chdir('med_fi');

f = open(absfilename,'a');
f.write('\nTEST %i\n' % testNum);
f.write('General, 3D, parallel, med-fi analysis\n');
f.write('Directory: %s\n\n' % os.getcwd());
f.close();
testNum += 1;

shutil.copyfile(os.path.join(rootdir,'example','general-3d.cfg'),'general-3d.cfg');
shutil.copyfile(os.path.join(rootdir,'example','general-3d.in'),'general-3d.in');
shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');

f = open(absfilename,'a');
f.write('Running...\n\n');
f.close();

call(['pythona','runModel.py','-f','general-3d.cfg','-l','2','-n',str(nCores)]);

compareResponses(absfilename,'results_general-3d_medfi_1.out','results.out');

# --------- 2D param, med-fi, f.d. gradients w/ mesh deformation -------------

os.chdir(os.path.join(rootdir,testdir,'gradients'));
os.chdir('med_fi');

f = open(absfilename,'a');
f.write('\nTEST %i\n' % testNum);
f.write('Mass and thrust 2D, parallel, med-fi F.D. gradients with mesh deformation\n');
f.write('Directory: %s\n\n' % os.getcwd());
f.close();
testNum += 1;

shutil.copyfile(os.path.join(rootdir,'example','gradients','general_fd_gradients.cfg'),'general_fd_gradients.cfg');
shutil.copyfile(os.path.join(rootdir,'example','gradients','params.in'),'params.in');
shutil.copyfile(os.path.join(rootdir,'runModel.py'),'runModel.py');

f = open(absfilename,'a');
f.write('Running...\n\n');
f.close();

call(['pythona','runModel.py','-f','general_fd_gradients.cfg','-l','1','-n',str(nCores)]);

compareResponses(absfilename,'results_general_fd_gradients_medfi.out','results.out');
