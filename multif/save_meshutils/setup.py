from distutils.core import setup, Extension

setup(ext_modules=[Extension("_meshutils_module",
      sources=["meshutils_py.c", "meshutils.c", "GMFio.c", "SU2io.c", "clean.c", "extraction.c", "libmesh6.c", "modules.c", "utils.c", "parser.c", "option.c", "mesh.c", "meshutils_py.i", "bspline3.c"])])

