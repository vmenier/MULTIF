from distutils.core import setup, Extension

setup(ext_modules=[Extension("_meshutils_module",
      sources=["./meshutils/meshutils_py.c", "./meshutils/meshutils.c", "./meshutils/GMFio.c", "./meshutils/SU2io.c", "./meshutils/clean.c", "./meshutils/extraction.c", "./meshutils/libmesh6.c", "./meshutils/modules.c", "./meshutils/utils.c", "./meshutils/parser.c", "./meshutils/option.c", "./meshutils/mesh.c", "./meshutils/bspline3.c", "./meshutils/meshutils_py.i", "./meshutils/GMSHio.c"],
       extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result"])#,
#       Extension('_nozzle_module',
#       sources = ['./meshutils/nozzle.cpp'],
##       include_dirs=['/home/avery/Projects/Gmsh/gmsh-2.13.1-install/include/gmsh'],
##       library_dirs=['/home/avery/Projects/Gmsh/gmsh-2.13.1-install/lib'],
#       libraries=['Gmsh'])
])

#setup(ext_modules=[Extension("_mshint_module",
#      sources=["./mshint/bucket.c"   , \
#			"./mshint/hash.c"     , \
#			"./mshint/inout.c"    , \
#			"./mshint/locelt.c"   , \
#			"./mshint/mi_calls.c" , \
#			"./mshint/mshint.c"   , \
#			"./mshint/chrono.c"   , \
#			"./mshint/libmesh5.c"   , \
#			"./mshint/mshint_py.i"   , \
#			"./mshint/mshint_py.c"   , \
#			"./mshint/mesh.c"   , \
#			"./mshint/SU2io.c"   , \
#			"./mshint/GMFio.c"   , \
#			"./mshint/mshint1.c" ] , \
#			extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result"])
#])
#

#
#setup(ext_modules=[Extension("_mshint_2_module",
#      sources=[ \
##"./mshint_2/boule.c"          , \
##			"./mshint_2/bucket.c"         , \
##			"./mshint_2/chrono.c"         , \
##			"./mshint_2/hash.c"           , \
##			"./mshint_2/inout.c"          , \
##			"./mshint_2/intelt.c"         , \
##			"./mshint_2/libmesh6.c"       , \
##			"./mshint_2/locelt.c"         , \
##			"./mshint_2/mesh.c"           , \
##			"./mshint_2/mshin1.c"         , \
#			"./mshint_2/mshint_py.c"      , \
##			"./mshint_2/scalem.c"         , \
##			"./mshint_2/GMFio.c"          , \
###			"./mshint_2/SU2io.c"          , \
##			"./mshint_2/mshint_py.i"   ,\
##x			"./mshint_2/mshint.c"         ] , \
#			extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result", \
#			"-O3","-c","-Wuninitialized","-Wshadow","-fexpensive-optimizations","-funroll-loops"])
#])
#
#


#setup(ext_modules=[Extension("_mshint_module",
#      sources=["./mshint/boule.c"      , \
#			"./mshint/bucket.c"     , \
#			"./mshint/chrono.c"     , \
#			"./mshint/hash.c"       , \
#			"./mshint/inout.c"      , \
#			"./mshint/intelt.c"     , \
#			"./mshint/libmesh5.c"   , \
#			"./mshint/locelt.c"     , \
#			"./mshint/mshin1.c"     , \
#			"./mshint/mshint.c"     , \
#			"./mshint/mshint_py.c"  , \
#			"./mshint/mshint_py.i"  , \
#			"./mshint/GMFio.c" , \
#			"./mshint/SU2io.c" , \
#			"./mshint/mesh.c" ,\
#			"./mshint/scalem.c"] , \
#			extra_compile_args=["-O3","-c","-Wuninitialized","-Wunused","-Winline","-Wshadow","-fexpensive-optimizations","-funroll-loops"])
#		#	extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result"])
#])
#
#
#