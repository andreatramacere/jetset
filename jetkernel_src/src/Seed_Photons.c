//-------------------------------------------------------------------------------------------------------
//
//  Seed Photons for IC process
//
//-------------------------------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"


double  I_nu_to_n(double I_nu, double nu){
    //----  phtons/(Hz ster cm^3) ------------------
    return   I_nu/(HPLANCK*nu*vluce_cm);
}

