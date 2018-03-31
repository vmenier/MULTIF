"""
A collection of useful utility functions.
"""

import os, sys

def makeDirAndLink(dirname, link_files):
    '''
    Make a new directory called dirname and link the files from the current
    directory listed in the link_files list to the new directory.
    '''
          
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    else:
        print('Directory %s already exists.' % dirname)

    # Link necessary files
    for f in link_files:
        ffrom = f # file linking from
        fto = f.split('/')[-1] # file linking to (may be joined path)
        if( os.path.exists(os.path.join(dirname,fto)) ):
            pass
        else:
            try:
                os.link(ffrom, os.path.join(dirname,fto))
            except OSError as err:
                print(err)
                print(ffrom)
                print(fto)
                print("ERROR: %s failed to link from %s." % (ffrom,os.getcwd()))
                #raise err

    return