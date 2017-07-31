from distutils.core import setup, Extension
import numpy
#meshutils.o clean.o extraction.o GMFio.o libmesh6.o mesh.o modules.o option.o parser.o SU2io.o utils.o
# define the extension module
meshutils_py = Extension('meshutils_py', sources=['meshutils_py.c', 'meshutils.c', 'meshutils.h', 'fproto.h', 'option.c', 'option.h'],
                          include_dirs=[numpy.get_include()])

# run the setup
setup(ext_modules=[meshutils_py])