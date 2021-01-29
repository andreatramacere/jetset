//================================================================
//                                                                                                    
//   Flux and Freq Transformations                                                          
//                                                                                                    
//================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"
//#include "nrutil.h"

/**
 * \file funz_math.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief funzioni matematiche
 *
 */

//=======================================
// Frequencies Transformations
//=======================================
//the disk rest frame is the source rest frame

double nu_disk_to_nu_obs_disk(double nu, double z) {
    return nu / (1 + z);
}

double nu_obs_disk_to_nu_disk(double nu, double z){
    return nu*(1+z);
}

//from blob to observer rest frame

double nu_blob_to_nu_obs(double nu, double delta, double z) {
    return nu * delta / (1 + z);
}
//

double nu_obs_to_nu_blob(double nu, double delta, double z) {
    return nu * (1 + z) / delta;
}

double nu_obs_to_nu_src(double nu, double z) {
    return nu * (1 + z);
}

//from blob to source rest frame

double nu_blob_to_nu_src(double nu, double delta, double z) {
    return nu*delta;
}
//=========================================================================================


//======================================
// j_nu to L_nu 4pi in n_seed() is multiplied here
// L_nu_blob is the luminosity in the blob
// rest frame so there is no beming
//======================================

double j_nu_to_L_nu_blob(double j_nu, double V) {
    return j_nu * four_pi*V;
}
//=========================================================================================


//==================================
// I_nu to L_nu
// here the angular dependence is pi not 4pi
// L_nu_blob is the luminosity in the blob
// rest frame so there is no beming
//==================================

double I_nu_to_L_nu_blob(double I_nu, double S) {
    //Kataoka Thesis pag. 195
    return I_nu * (pi * S); /*erg s^-1  Hz^-1 */
}
//=========================================================================================


//======================================
// j_nu to L_nu 4pi in n_seed() is multiplied here
// L_nu_src is the luminosity in the source
// rest frame so there is  the beaming 
//======================================

double j_nu_to_L_nu_src(double j_nu, double V, double beam_obj) {
    return beam_obj * beam_obj * beam_obj * j_nu * four_pi*V;
}

double j_nu_src_to_L_nu_src(double j_nu, double V, double beam_obj)
{
    return beam_obj  * j_nu * four_pi * V;
}
//=========================================================================================


//==================================
// I_nu to L_nu
// here the angular dependence is pi not 4pi
// L_nu_src is the luminosity in the source
// rest frame so there is  the beaming 
//==================================

double I_nu_to_L_nu_src(double I_nu, double S, double beam_obj) {
    //Kataoka Thesis pag. 195
  
    return beam_obj * beam_obj * beam_obj * I_nu * (pi * S); /*erg s^-1  Hz^-1 */
}
//=========================================================================================


//===============================
// L_nu to F_nu
//F(nu_obs)=L(nu_blob)*(1+z)/(4*pi*dl^2)*(delta^3)
//===============================

double L_nu_blob_to_F_nu(double L_nu, double beam_obj, double z, double dist) {
    return L_nu * (beam_obj * beam_obj * beam_obj)*(1 + z) / (dist * dist * four_pi);
}

double L_nu_src_to_F_nu(double L_nu, double beam_obj, double z, double dist) {
    return L_nu * (1 + z) / (dist * dist * four_pi);
}

double L_nu_Disk_to_F_nu(double L_nu, double z, double dist) {
    return L_nu * (1 + z) / (dist * dist * four_pi);
}

//=========================================================================================

//===============================
// nuFnu_nu to nuL_nu
//nu_obsF(nu_obs)=nu_blob*L(nu_blob)/(4*pi*dl^2)*dellta^4
//nu_blob*L(nu_blob)=nu_obsF(nu_obs)*(4*pi*ld^2)/(delta^4)
//===============================

double nuFnu_obs_to_nuLnu_src(double nuFnu_obs, double beam_obj, double z, double dist) {
    return nuFnu_obs * (dist * dist * four_pi);
}

double nuFnu_obs_to_nuLnu_blob(double nuFnu_obs, double beam_obj, double z, double dist) {
    return nuFnu_obs / (beam_obj * beam_obj * beam_obj * beam_obj) * (dist * dist * four_pi);
}
//=========================================================================================
