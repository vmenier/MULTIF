"""
Take a file containing design and random variable samples and produce MULTI-F 
input files. It is assumed the sample file already has the data in the right order.

Rick Fenrich 9/20/17
"""

import numpy as np;

sampleFile = '3d_samples.dat';
outputFilePrefix = '3d_samples';

# Read sample file
samples = np.loadtxt(sampleFile);
n,m = samples.shape;
print('%i samples read of %i points each.' % (n,m));
        
# Print input files
for i in range(n):
    ftmp = outputFilePrefix + '_' + '%i' % i + '.in';
    np.savetxt(ftmp,samples[i,:]);
print('%i input files written.' % n);