/* jetkernel.i */
 %module jetkernel
 %{
 /* Put header files here or function declarations like below */
 #include "../../jetkernel_src/include/Blazar_SED.h"
 %}

 /* Parse the header file to generate wrappers */
 %include "../../jetkernel_src/include/Blazar_SED.h"