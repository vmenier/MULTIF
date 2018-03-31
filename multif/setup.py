from distutils.core import setup, Extension

setup(ext_modules=[ \
      Extension("./models/_meshutils_module",
      sources=["./models/meshutils/meshutils_py.c", \
			   "./models/meshutils/meshutils.c", \
			   "./models/meshutils/GMFio.c", \
			   "./models/meshutils/SU2io.c", \
			   "./models/meshutils/clean.c", \
			   "./models/meshutils/extraction.c", \
			   "./models/meshutils/libmesh6.c", \
			   "./models/meshutils/modules.c", \
			   "./models/meshutils/utils.c", \
			   "./models/meshutils/parser.c", \
			   "./models/meshutils/option.c", \
			   "./models/meshutils/mesh.c", \
			   "./models/meshutils/nozzle.c", \
			   "./models/meshutils/bspline3.c", \
			   "./models/meshutils/piecewise.c", \
			   "./models/meshutils/meshutils_py.i", \
			   "./models/meshutils/projection.c", \
			   "./models/meshutils/GMSHio.c"],
       extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result"]),
       
       Extension('./models/_nozzle_module',
       sources = ['./models/meshutils/nozzle.cpp'],
       extra_compile_args=["-Wno-maybe-uninitialized","-std=c++11","-Wno-sign-compare","-Wno-unused-but-set-variable","-DHAVE_NO_OCC_CONFIG_H"],
       libraries=['Gmsh','TKOffset', 'TKTopAlgo', 'TKGeomAlgo', 'TKBRep', 'TKGeomBase', 'TKG3d', 'TKG2d',
                  'TKMath', 'TKernel', 'TKBool', 'TKFeat', 'TKShHealing']),
       
       Extension("./models/_mshint_module",
       sources=["./models/mshint/boule.c"      , \
			    "./models/mshint/bucket.c"     , \
			    "./models/mshint/chrono.c"     , \
			    "./models/mshint/hash.c"       , \
			    "./models/mshint/inout.c"      , \
			    "./models/mshint/intelt.c"     , \
			    "./models/mshint/libmesh5.c"   , \
			    "./models/mshint/locelt.c"     , \
			    "./models/mshint/mshin1.c"     , \
			    "./models/mshint/mshint.c"     , \
			    "./models/mshint/mshint_py.c"  , \
			    "./models/mshint/mshint_py.i"  , \
			    "./models/mshint/GMFio.c" , \
			    "./models/mshint/SU2io.c" , \
			    "./models/mshint/mesh.c" ,\
			    "./models/mshint/scalem.c"] , \
	   extra_compile_args=["-O3","-c","-Wuninitialized","-Wunused",
	   "-Winline","-Wshadow","-fexpensive-optimizations","-funroll-loops"]),   
       
       Extension("./models/aero/SU2/amginria/_amgio",
       sources=[ "./models/aero/SU2/amginria/amgio/amgio_py.c", \
       			 "./models/aero/SU2/amginria/amgio/mesh.c", \
       			 "./models/aero/SU2/amginria/amgio/GMFio.c", \
       			 "./models/aero/SU2/amginria/amgio/SU2io.c", \
       			 "./models/aero/SU2/amginria/amgio/option.c", \
       			 "./models/aero/SU2/amginria/amgio/libmesh6.c", \
                 "./models/aero/SU2/amginria/amgio/amgio_py.i", \
       			 "./models/aero/SU2/amginria/amgio/convert.c"],
        extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result"]), 
       	   
	  Extension("./models/aero/LOWF/_quasi1dnozzle",
      sources=["./models/aero/LOWF/quasi1dnozzle_py.c", \
               "./models/aero/LOWF/lofinozzle.c", \
               "./models/aero/LOWF/odeint.c", \
               "./models/meshutils/piecewise.c", \
               "./models/aero/LOWF/quasi1dnozzle_py.i"],
      extra_compile_args=["-std=c99","-Wno-unused-variable","-Wno-unused-result"])
    
       ])
       
       
       
