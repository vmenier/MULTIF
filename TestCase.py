"""
Python class for running regression test cases for MULTI-F

Rick Fenrich 8/29/17
"""

import time, os, shutil, subprocess, datetime, sys
import numpy as np

class TestCase:
    
    def __init__(self,name):
        
        self.name = name
        self.description = ''
        self.log_file = os.path.join(os.getcwd(),'regression.out')
        
        # Configuration and input file information
        self.multif_dir = os.path.dirname(os.path.abspath(__file__))
        self.cfg_dir = os.getcwd()
        
        self.working_dir = os.getcwd()
        
        self.cfg_file = 'general.cfg'
        self.input_file = 'general.in'
        self.dependencies = [] # list of strings giving other required MULTI-F files
        self.fidelity = 0
        self.ntasks = 1 # number of concurrent tasks (i.e. model runs)
                         # only used when finite difference gradients are used
        self.cpus_per_task = 1 # number of cores to be used for each model run
        
        # Termination information
        self.timeout = 15*60 # s
        self.fail_tol = 1e-2 # relative percent difference for determining failure
        self.diff_tol = 1e-6 # relative percent difference for determing diff
        
        # Comparison information
        self.compare_file = 'results_general.out'
     
    # return -1 for failure tolerance exceeded
    # return -2 for difference tolerance exceeded
    # return -3 for wrong number of results included
    # return -4 for no results file
    def compare_responses(self):

        return_flag = 0    
        
        f = open(self.log_file,'a')
        f.write('Response file 1: %s\n' % os.path.join(self.multif_dir,self.compare_file))
        f.write('Response file 2: %s\n' % os.path.join(os.getcwd(),'results.out'))
        f.close()
                
        # Obtain data from example results file
        r1 = open(os.path.join(self.multif_dir,self.compare_file),'r')
        val1 = []
        grad1 = []
        val_name = []
        for line in r1:
            linelist = line.split()
            if( len(linelist) > 2 ): # gradients are recorded
                grad1.append([float(i) for i in linelist[1:-1]])
            else: # values are recorded
                val1.append(float(linelist[0]))
                if( len(linelist) > 1 ):
                    val_name.append(linelist[1])
                else:
                    val_name.append('')
        r1.close()
        
        # Obtain data from results file just obtained
        if( not os.path.isfile('results.out') ):
            return -4
        r2 = open('results.out','r')
        val2 = []
        grad2 = []
        for line in r2:
            linelist = line.split()
            if( len(linelist) > 2 ): # gradients are recorded
                grad2.append([float(i) for i in linelist[1:-1]])
            else: # values are recorded
                val2.append(float(linelist[0]))
        r2.close()
        
        # Find differences in response values
        valDiff = []
        if( len(val1) != len(val2) ):
            print 'Obtained results.out file does not contain same number of entries as example results file.'
            return -3
        for i in range(len(val1)):
            if( val1[i] != 0. ):
                valDiff.append((val2[i]-val1[i])/val1[i])
            else:
                valDiff.append((val2[i]-val1[i]))
        
        # Find differences in response gradients
        gradDiff = []
        if( len(grad1) != len(grad2) ):
            print 'Obtained results.out file does not contain same number of gradient entries as example results file.'
            return -3
        for i in range(len(grad1)):
            if len(grad1[i]) != len(grad2[i]):
                print 'Obtained results.out file does not contain same number gradient variables as example results file.'
                return -3
        for i in range(len(grad1)):
            for j in range(len(grad1[i])):
                if( grad1[i][j] != 0. ):
                    gradDiff.append((grad2[i][j]-grad1[i][j])/grad1[i][j])
                else:
                    gradDiff.append((grad2[i][j]-grad1[i][j]))
        
        # Check for exceedance
        maxValDiff = np.max(np.abs(np.array(valDiff)))
        if( maxValDiff > self.diff_tol or np.isnan(maxValDiff) or np.isinf(maxValDiff) ):
            if( maxValDiff > self.fail_tol ):
                return_flag = -1
            else:
                return_flag = -2
            maxValDiffLoc = np.argmax(np.abs(np.array(valDiff)))
            f = open(self.log_file,'a')
            f.write('Max relative difference in response value is %f at %i\n' % (maxValDiff, maxValDiffLoc))
            f.write('{0:6} {1:4} {2:25} {3:25} {4}\n'.format('Line #', '', 'Expected', 'Obtained', 'Name'))           
            for i in range(len(valDiff)):
                if( np.max(np.abs(valDiff[i])) > self.fail_tol or np.isnan(valDiff[i]) or np.isinf(valDiff[i]) ):
                    f.write('{0:6} {1:4} {2:25.16f} {3:25.16f} {4}\n'.format(i, '*', val1[i], val2[i], val_name[i]))
                elif( np.max(np.abs(valDiff[i])) > self.diff_tol ):
                    f.write('{0:6} {1:4} {2:25.16f} {3:25.16f} {4}\n'.format(i, 'DIFF', val1[i], val2[i], val_name[i]))
                else:
                    f.write('{0:6} {1:4} {2:25.16f} {3:25.16f} {4}\n'.format(i, '', val1[i], val2[i], val_name[i]))
            f.close()
        
        if( len(gradDiff) > 0 ):
            maxGradDiff = np.max(np.abs(np.array(gradDiff)))
            if( maxGradDiff > self.diff_tol or np.isnan(maxGradDiff) or np.isinf(maxValDiff) ):
                if( maxGradDiff > self.fail_tol ):
                    return_flag = -1
                else:
                    return_flag = -2
                f = open(self.log_file,'a')
                f.write('Max relative difference in gradient value is %f\n' % maxGradDiff)
                f.write('{0:3} {1:3} {2:4} {3:25} {4:25}\n'.format('Row', 'Col', '', 'Expected', 'Obtained'))           
                k = 0
                for i in range(len(grad1)):
                    for j in range(len(grad1[i])):
                        if( np.max(np.abs(gradDiff[k])) > self.fail_tol or np.isnan(gradDiff[k]) or np.isinf(gradDiff[k]) ):
                            f.write('{0:3} {1:3} {2:4} {3:25.16f} {4:25.16f}\n'.format(i+len(valDiff), j, '*', grad1[i][j], grad2[i][j]))
                        elif( np.max(np.abs(gradDiff[k])) > self.diff_tol ):
                            f.write('{0:3} {1:3} {2:4} {3:25.16f} {4:25.16f}\n'.format(i+len(valDiff), j, 'DIFF', grad1[i][j], grad2[i][j]))                          
                        else:
                            f.write('{0:3} {1:3} {2:4} {3:25.16f} {4:25.16f}\n'.format(i+len(valDiff), j, '', grad1[i][j], grad2[i][j]))
                        k = k + 1
                f.close()            
    
        return return_flag
        
    def run_test(self):
        
        passed = 1
        timed_out = 0
        
        # Go to MULTI-F directory and run all tests from here
        #os.chdir(self.multif_dir)
        
        # Create directory to run test in 
        if not os.path.isdir('local'):
            os.mkdir('local')
            os.chdir('local')
        else:
            os.chdir('local')
            
        if os.path.isdir(self.name):
            shutil.rmtree(self.name)
        os.mkdir(self.name)
        os.chdir(self.name)
        
        workdir = os.getcwd()
                
        # Write preliminary information to log file
        f = open(self.log_file,'a')
        f.write('======================= Test: %s =======================\n' % self.name)
        f.write('%s\n' % self.description)
        f.write('Directory: %s\n\n' % os.getcwd())
        f.close()
            
        # Copy all necessary files to directory
        shutil.copyfile(os.path.join(self.multif_dir,self.cfg_dir,self.cfg_file),self.cfg_file)
        shutil.copyfile(os.path.join(self.multif_dir,self.cfg_dir,self.input_file),self.input_file)
        shutil.copyfile(os.path.join(self.multif_dir,'runModel.py'),'runModel.py')
        for f in self.dependencies:
            shutil.copyfile(os.path.join(self.multif_dir,self.cfg_dir,f),f)
            
        # Build shell command to run MULTI-F
        command = ['python','%s/runModel.py'%self.multif_dir,'-f',self.cfg_file,'-l',str(self.fidelity),'-n',str(self.ntasks),'-c',str(self.cpus_per_task)]
        
        # Write command line in file
        f = open(self.log_file,'a')
        f.write('Command line: %s\n' % ' '.join(command))
        f.close()
        
        # Run MULTI-F
        start = datetime.datetime.now()
        
#        process = subprocess.Popen(command, shell = True) # start MULTI-F
        try:
            so = open('stdout.txt','w')
            se = open('stderr.txt','w')
            sys.stdout.write('\nRunning test case %s\n' % self.name)
            sys.stdout.write('Working directory: %s\n' % os.getcwd())
            sys.stdout.write('Comparison data file: %s\n' % os.path.join(self.multif_dir,self.compare_file))
            subprocess.call(command,stdout=so,stderr=se)
            so.close()
            se.close()
        except:
            print 'Error encountered during execution'
            f = open(self.log_file,'a')
            f.write('Test failed -- MULTI-F error encountered\n')
            f.write('===================== End Test: %s =====================\n\n' % self.name)
            f.close()
            # Return to MULTI-F directory
            #os.chdir(self.multif_dir)
            os.chdir(self.working_dir)
            return 0
        
        # IMPLEMENT PROCESS KILL HERE FOR TIMEOUT
        
        end = datetime.datetime.now()
        running_time = (end - start).seconds
        
#        # Check for timeout
#        while process.poll() is None:
#            time.sleep(0.1)
#            now = datetime.datetime.now()
#            running_time = (now - start).seconds
#            if running_time > self.timeout:
#                try:
#                    process.kill()
#                    os.system('killall %s' % 'python')
#                except AttributeError:
#                    pass
#                timed_out = 1
#                passed = 0
        
        # Examine and compare output
        
        return_flag = self.compare_responses()        
        if( return_flag == -1 ):
            print 'Test failed -- tolerance exceeded (%e)' % self.fail_tol
            passed = 0
        elif( return_flag == -2 ):
            print 'Test diffed -- tolerance exceeded (%e)' % self.diff_tol
            f = open(self.log_file,'a')
            f.write('Test diffed\n')
            f.close()
            passed = 2            
        elif( return_flag == -3 ):
            print 'Test failed -- wrong number of results in results file'
            f = open(self.log_file,'a')
            f.write('Test failed -- wrong number of results in results file\n')
            f.close()
            passed = 0
        elif( return_flag == -4 ):
            print 'Test failed -- no results file'
            f = open(self.log_file,'a')
            f.write('Test failed -- no results file\n')
            f.close()
            passed = 0
        elif( timed_out ):
            print 'Test failed -- timed out'
            f = open(self.log_file,'a')
            f.write('Test failed -- timed out\n')
            f.close()
            passed = 0
        else:
            print 'Test passed'
            f = open(self.log_file,'a')
            f.write('Test passed\n')
            f.close()              
        
        # Closing remarks
        f = open(self.log_file,'a')
        f.write('Test duration: %.2f min\n' % (running_time/60.) )
        f.write('===================== End Test: %s =====================\n\n' % self.name)
        f.close()

        # Return to MULTI-F directory
        #os.chdir(self.multif_dir)    
        os.chdir(self.working_dir)
        
        return passed
