//=========================================================================================
//                   CALCOLO DELLO SPETTRO bremsstrahlung ep
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
 * \date 24-12-2010
 * \brief CALCOLO DELLO SPETTRO COMPTON
 *
 */


void spettro_bremss_ep(int Num_file, struct blob *pt) {
    double  k, nu_1, nu_check;
    double L_nu_bremss_ep, F_nu_bremss_ep_obs;
    double log_nu_start,gmax;
    unsigned int NU_INT, i, I_MAX, stop;
    
    //============================================================
    //         inizio  loop sulle freq per spettro  pp
    //============================================================
    stop = 0;



    // massima e minima freq bremss
    gmax=Find_gmax(pt,pt->Ne,pt->griglia_gamma_Ne_log);
    pt->nu_stop_bremss_ep_pred = gmax*MEC2/HPLANCK*10;
    pt->nu_start_bremss_ep = pt->gmin_griglia*MEC2/HPLANCK/100;
    pt->nu_start_bremss_ep_obs = nu_blob_to_nu_obs(pt->nu_start_bremss_ep, pt->beam_obj, pt->z_cosm);
    pt->nu_stop_bremss_ep_obs = nu_blob_to_nu_obs(pt->nu_stop_bremss_ep_pred, pt->beam_obj, pt->z_cosm);
    nu_check=(pt->nu_start_bremss_ep)*0.5;
    NU_INT = 0;
    k = (log10(pt->nu_stop_bremss_ep_pred) - log10(pt->nu_start_bremss_ep));
    I_MAX = pt->nu_IC_size;
    if (pt->verbose)
    {
        printf("**********************  CALCOLO DELLO SPETTRO bremss ep   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->nu_start_bremss_ep,
               pt->nu_stop_bremss_ep_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }
    

    log_nu_start = log10(pt->nu_start_bremss_ep);
    for (i = 0; i < I_MAX; i++) {

        // calcola nu con incremento logaritmico
        nu_1 = pow(10, log_nu_start + k * (double) i / (double) I_MAX);
        if (i == 0) {
            // !!! ATTENZIONE POICHE' nu per errori
            // di arrot. numerico puo' avere valore
            // iniziale <nu_start_comp allora impongo
            nu_1 = pt->nu_start_bremss_ep;
        }
        pt->nu_1 = nu_1;
        pt->nu_bremss_ep[NU_INT]=nu_1;
        pt->nu_bremss_ep_obs[NU_INT] = nu_blob_to_nu_obs(nu_1, pt->beam_obj, pt->z_cosm);
       
        if ((nu_1 >= pt->nu_start_bremss_ep) && (nu_1 <= pt->nu_stop_bremss_ep_pred)) {
            //printf("hi\n");
            if (!stop) {
                //rate_gamma_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
                //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
                //then you have to multipli by HPLANCK_TeV*nu->TeV*(TeV^-1 cm^-3 s^-1)
                //the you multiply by HPLANCK in  to get erg/( cm^3 Hz s) 
                //and then you divide by 4pi to get erg/( cm^3 Hz s setard) that are j_nu units 

                pt->j_bremss_ep[NU_INT] = pt->NH_pp*j_nu_bremss_ep(pt,pt->nu_bremss_ep[NU_INT]);
                if (pt->verbose) {
                    printf("#-> NU_INT=%d j[NU_INT]=%e nu_1=%e i=%d \n",
                            NU_INT,
                            pt->j_bremss_ep[NU_INT],
                            nu_1, i);
                }
                //nu_src = nu_blob_to_nu_src(nu_1, pt->beam_obj, pt->z_cosm);
                L_nu_bremss_ep = j_nu_to_L_nu_src(pt->j_bremss_ep[NU_INT], pt->Vol_sphere, pt->beam_obj);
                //nuL_nu_ep_brem = L_nu_bremss_ep*nu_src;
                F_nu_bremss_ep_obs = L_nu_src_to_F_nu(L_nu_bremss_ep, pt->beam_obj, pt->z_cosm, pt->dist);
                pt->nuFnu_bremss_ep_obs[NU_INT] = F_nu_bremss_ep_obs * pt->nu_bremss_ep_obs[NU_INT];

                pt->nu_stop_bremss_ep = nu_1;
                pt->NU_INT_STOP_BREMSS_EP = NU_INT;
                if (pt->verbose) {
                    printf("nu_stop_brems_ep_pred=%e nu_stop_bremss_ep=%e NU_INT=%d\n ",
                            pt->nu_stop_bremss_ep_pred,
                            pt->nu_stop_bremss_ep,
                            NU_INT);
                }
            }
            if (pt->j_bremss_ep[NU_INT] < pt->emiss_lim) {
                pt->j_bremss_ep[NU_INT] = pt->emiss_lim;
                pt->nuFnu_bremss_ep_obs[NU_INT] = pt->emiss_lim;
                if (nu_1>nu_check){
                    stop = 1;

                    F_nu_bremss_ep_obs = pt->emiss_lim;
                }
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
        pt->NU_INT_STOP_BREMSS_EP = NU_INT - 1;
    }
    //printf("nu_stop_pp=%e NU_INT_STOP_PP=%d\n", pt->nu_stop_pp, pt->NU_INT_STOP_PP);
    pt->nu_start_bremss_ep_obs = nu_blob_to_nu_obs(pt->nu_start_bremss_ep, pt->beam_obj, pt->z_cosm);
    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

        FindEpSp(pt->nu_bremss_ep_obs, pt->nuFnu_bremss_ep_obs,   pt->NU_INT_STOP_BREMSS_EP, pt,
                &(pt->nu_peak_bremss_ep_obs),
                &(pt->nu_peak_bremss_ep_src),
                &(pt->nu_peak_bremss_ep_blob),
                &(pt->nuFnu_peak_bremss_ep_obs),
                &(pt->nuLnu_peak_bremss_ep_src),
                &(pt->nuLnu_peak_bremss_ep_blob));

        if (pt->verbose)
        {
            printf("nu_bremss_ep_blob peak=%e\n", pt->nu_peak_bremss_ep_blob);
            printf("nu_bremss_ep_src   peak=%e\n", pt->nu_peak_bremss_ep_src);
            printf("nu_bremss_ep_obs  peak=%e\n", pt->nu_peak_bremss_ep_src);

            printf("nuFnu bremss_ep  blob    peak=%e\n", pt->nuFnu_peak_bremss_ep_obs);
            printf("nuLnu bremss_ep  src      peak=%e\n", pt->nuLnu_peak_bremss_ep_src);
            printf("nuLnu bremss_ep  obs     peak=%e\n", pt->nuLnu_peak_bremss_ep_blob);
        }
        return;
}
//=========================================================================================

