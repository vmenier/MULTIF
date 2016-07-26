from distutils.core import setup, Extension

setup(ext_modules=[Extension("_meshutils_module",
      sources=["./meshutils/meshutils_py.c", "./meshutils/meshutils.c", "./meshutils/GMFio.c", "./meshutils/SU2io.c", "./meshutils/clean.c", "./meshutils/extraction.c", "./meshutils/libmesh6.c", "./meshutils/modules.c", "./meshutils/utils.c", "./meshutils/parser.c", "./meshutils/option.c", "./meshutils/mesh.c", "./meshutils/bspline3.c", "./meshutils/meshutils_py.i"],
       extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result"]),
       Extension('_nozzle_module',
       sources = ['./meshutils/nozzle.cpp'],
#       include_dirs=['/home/avery/Projects/Gmsh/gmsh-2.13.1-install/include/gmsh'],
#       library_dirs=['/home/avery/Projects/Gmsh/gmsh-2.13.1-install/lib'],
       libraries=['Gmsh'])])

