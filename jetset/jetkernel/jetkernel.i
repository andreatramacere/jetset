/* jetkernel.i */
 %module jetkernel
 %{
    #define SWIG_FILE_WITH_INIT
    /* Put header files here or function declarations like below */
    #include "../../jetkernel_src/include/Blazar_SED.h"
%}

/*
%include "numpy.i"
%init %{
import_array();
%}
*/

 /* Parse the header file to generate wrappers */
 %include "../../jetkernel_src/include/Blazar_SED.h"