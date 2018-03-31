/*  Example of wrapping cos function from math.h using SWIG. */

%module meshutils_module
%{
    /* the resulting C file should be built as a python extension */
    #define SWIG_FILE_WITH_INIT
    /*  Includes the header in the wrapper code */
	#include "./models/meshutils/meshutils_py.h" 
%}

%include "numpy.i"


/*  Parse the header file to generate wrappers */
%include "./models/meshutils/meshutils_py.h"

/*
%include "typemaps.i"
%apply int *OUTPUT { int *pyNbrRes, int *pyResSiz };
%inlne %{
extern void py_ExtractAlongLine (char *MshNam, char *SolNam, PyObject *pyResult, int *pyNbrRes, int *pyResSiz)
%}
*/
