//=========================================================================================
//                   CALCOLO DELLO SPETTRO Gamma PP
//=========================================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file spettro_Compton.c
 * \author Andrea Tramacere
 * \date 04-05-2010
 * \brief CALCOLO DELLO SPETTRO COMPTON
 *
 */


void spettro_pp_gamma(int Num_file, struct blob *pt) {
    double  k, nu_1;
    double L_nu_pp, F_nu_pp_obs;
    double log_nu_start, gmax;
    unsigned int NU_INT, i, I_MAX, stop;
    
    //============================================================
    //         inizio  loop sulle freq per spettro  pp
    //============================================================
    stop = 0;


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!RICODATI DI CAMBIARE la distr e- con p
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // massima e minima freq pp
    gmax=Find_gmax(pt,pt->Np,pt->griglia_gamma_Np_log);
    pt->nu_stop_pp_gamma_pred = gmax * MPC2 / HPLANCK * 100;
    pt->nu_start_pp_gamma = E_th_pp * 1E12 * ev_to_erg / HPLANCK / 100;
    pt->nu_start_pp_gamma_obs = nu_blob_to_nu_obs(pt->nu_start_pp_gamma, pt->beam_obj, pt->z_cosm);
    pt->nu_stop_pp_gamma_obs = nu_blob_to_nu_obs(pt->nu_stop_pp_gamma_pred, pt->beam_obj, pt->z_cosm);
   
    NU_INT = 0;
    k = (log10(pt->nu_stop_pp_gamma_pred) - log10(pt->nu_start_pp_gamma));
    I_MAX = pt->nu_IC_size;
    if (pt->verbose)
    {
        printf("**********************  CALCOLO DELLO SPETTRO pp   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->nu_start_pp_gamma,
               pt->nu_stop_pp_gamma_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }
    

    log_nu_start = log10(pt->nu_start_pp_gamma);
    for (i = 0; i < I_MAX; i++) {

        // calcola nu con incremento logaritmico
        nu_1 = pow(10, log_nu_start + k * (double) i / (double) I_MAX);
        if (i == 0) {
            // !!! ATTENZIONE POICHE' nu per errori
            // di arrot. numerico puo' avere valore
            // iniziale <nu_start_comp allora impongo
            nu_1 = pt->nu_start_pp_gamma;
        }
        pt->nu_1 = nu_1;
        pt->nu_pp_gamma[NU_INT]=nu_1;
        pt->nu_pp_gamma_obs[NU_INT] = nu_blob_to_nu_obs(nu_1, pt->beam_obj, pt->z_cosm);
       
        if ((nu_1 >= pt->nu_start_pp_gamma) && (nu_1 <= pt->nu_stop_pp_gamma_pred)) {
            //printf("hi\n");
            if (!stop) {
                //rate_gamma_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
                //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
                //then you have to multipli by HPLANCK_TeV*nu->TeV*(TeV^-1 cm^-3 s^-1)
                //the you multiply by HPLANCK in  to get erg/( cm^3 Hz s) 
                //and then you divide by 4pi to get erg/( cm^3 Hz s setard) that are j_nu units 

                pt->j_pp_gamma[NU_INT] = vluce_cm * pt->NH_pp * bn_to_cm2 *
                        (HPLANCK)* (HPLANCK_TeV * pt->nu_1) *one_by_four_pi*
                        rate_gamma_pp(pt);
                if (pt->verbose) {
                    printf("#-> NU_INT=%d j[NU_INT]=%e nu_1=%e i=%d \n",
                            NU_INT,
                            pt->j_pp_gamma[NU_INT],
                            nu_1, i);
                }
                //nu_src = nu_blob_to_nu_src(nu_1, pt->beam_obj, pt->z_cosm);
                L_nu_pp = j_nu_to_L_nu_src(pt->j_pp_gamma[NU_INT], pt->Vol_sphere, pt->beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->beam_obj, pt->z_cosm, pt->dist);
                pt->nuFnu_pp_gamma_obs[NU_INT] = F_nu_pp_obs * pt->nu_pp_gamma_obs[NU_INT];

                pt->nu_stop_pp_gamma = nu_1;
                pt->NU_INT_STOP_PP_GAMMA = NU_INT;
                if (pt->verbose) {
                    printf("nu_stop_pp_pred=%e nu_stop_pp=%e NU_INT=%d\n ",
                            pt->nu_stop_pp_gamma_pred,
                            pt->nu_stop_pp_gamma,
                            NU_INT);
                }
            }
            if (pt->j_pp_gamma[NU_INT] < pt->emiss_lim) {
                stop = 1;
                pt->j_pp_gamma[NU_INT] = pt->emiss_lim;
                pt->nuFnu_pp_gamma_obs[NU_INT] = pt->emiss_lim;
                F_nu_pp_obs = pt->emiss_lim;
                if (pt->verbose) {
                    printf("%e %d\n ", nu_1, NU_INT);
                }
            }

            if (pt->verbose) {
                printf("#-> ********************************\n\n");
            }
            //==========================  END of Loop ove frequencies ====================================
        }
        NU_INT++;
    }

    //Se ancora non ha trovato nu_stop
    if (!stop) {
        pt->NU_INT_STOP_PP_GAMMA = NU_INT - 1;
    }
    //printf("nu_stop_pp=%e NU_INT_STOP_PP=%d\n", pt->nu_stop_pp, pt->NU_INT_STOP_PP);
    pt->nu_stop_pp_gamma_obs = nu_blob_to_nu_obs(pt->nu_stop_pp_gamma, pt->beam_obj, pt->z_cosm);
    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

        FindEpSp(pt->nu_pp_gamma, pt->nuFnu_pp_gamma_obs,   pt->NU_INT_STOP_PP_GAMMA, pt,
                &(pt->nu_peak_PP_gamma_obs),
                &(pt->nu_peak_PP_gamma_src),
                &(pt->nu_peak_PP_gamma_blob),
                &(pt->nuFnu_peak_PP_gamma_obs),
                &(pt->nuLnu_peak_PP_gamma_src),
                &(pt->nuLnu_peak_PP_gamma_blob));

        if (pt->verbose)
        {
            printf("nu_PP_blob peak=%e\n", pt->nu_peak_PP_gamma_blob);
            printf("nu_PP_src   peak=%e\n", pt->nu_peak_PP_gamma_src);
            printf("nu_PP_obs  peak=%e\n", pt->nu_peak_PP_gamma_obs);

            printf("nuFnu PP  blob    peak=%e\n", pt->nuFnu_peak_PP_gamma_obs);
            printf("nuLnu PP  src      peak=%e\n", pt->nuLnu_peak_PP_gamma_src);
            printf("nuLnu PP  obs     peak=%e\n", pt->nuLnu_peak_PP_gamma_blob);
        }
        return;
}
//=========================================================================================

