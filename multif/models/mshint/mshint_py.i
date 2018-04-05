/*  Example of wrapping cos function from math.h using SWIG. */

%module mshint_module
%{
    /* the resulting C file should be built as a python extension */
    #define SWIG_FILE_WITH_INIT
    /*  Includes the header in the wrapper code */
		#include "mshint_py.h" 
%}

/*  Parse the header file to generate wrappers */
%include "mshint_py.h"

