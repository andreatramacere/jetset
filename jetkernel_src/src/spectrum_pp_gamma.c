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
    double L_nu_pp, F_nu_pp_obs;
    double  gmax;
    unsigned int NU_INT, i, I_MAX, stop;
    void *(*eval_j_ptr)(void * args);
    //============================================================
    //         inizio  loop sulle freq per spettro  pp
    //============================================================
    stop = 0;


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!RICODATI DI CAMBIARE la distr e- con p
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // massima e minima freq pp
    gmax=Find_gmax(pt,pt->emitters.Np,pt->emitters.griglia_gamma_Np_log);
    pt->PP_gamma.nu_stop_pp_gamma_pred = gmax * MPC2 / HPLANCK * 100;
    pt->PP_gamma.spec.nu_min = E_th_pp * 1E12 * ev_to_erg / HPLANCK / 100/10;
    pt->PP_gamma.spec.nu_min_obs = nu_blob_to_nu_obs(pt->PP_gamma.spec.nu_min, pt->core.beam_obj, pt->core.z_cosm);
    pt->PP_gamma.spec.nu_max_obs = nu_blob_to_nu_obs(pt->PP_gamma.nu_stop_pp_gamma_pred, pt->core.beam_obj, pt->core.z_cosm);
   
    NU_INT = 0;
    //k = (log10(pt->PP_gamma.nu_stop_pp_gamma_pred) - log10(pt->PP_gamma.spec.nu_min));

    build_log_grid(pt->PP_gamma.spec.nu_min,  pt->PP_gamma.nu_stop_pp_gamma_pred, pt->core.nu_IC_size, pt->PP_gamma.spec.nu);
    build_log_grid(pt->PP_gamma.spec.nu_min_obs,  pt->PP_gamma.spec.nu_max_obs, pt->core.nu_IC_size, pt->PP_gamma.spec.nu_obs);


    I_MAX = pt->core.nu_IC_size -1;
    eval_j_ptr = &eval_j_pp_gamma;
    pt->PP_gamma.pp_racc_gamma=rate_gamma_pp(pt ,pt->PP_gamma.spec.nu_min,1);
    threaded_j_evaluation(pt, eval_j_ptr, pt->PP_gamma.spec.j_nu,pt->PP_gamma.spec.nu,pt->PP_gamma.spec.nu_min, pt->PP_gamma.nu_stop_pp_gamma_pred,I_MAX,pt->core.N_THREADS);
    if (pt->core.verbose){
        printf("**********************  CALCOLO DELLO SPETTRO pp   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->PP_gamma.spec.nu_min,
               pt->PP_gamma.nu_stop_pp_gamma_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }

    for (i = 0; i <= I_MAX; i++) {
        if ((pt->PP_gamma.spec.nu[i] >= pt->PP_gamma.spec.nu_min) && (pt->PP_gamma.spec.nu[i] <= pt->PP_gamma.nu_stop_pp_gamma_pred)) {
            //printf("hi\n");
            if (!stop) {

                L_nu_pp = j_nu_to_L_nu_src(pt->PP_gamma.spec.j_nu[NU_INT], pt->core.Vol_region, pt->core.beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);
                pt->PP_gamma.spec.nuFnu_obs[NU_INT] = F_nu_pp_obs * pt->PP_gamma.spec.nu_obs[NU_INT];

                pt->PP_gamma.spec.nu_max = pt->PP_gamma.spec.nu[i];
                pt->PP_gamma.NU_INT_STOP_PP_GAMMA = NU_INT;
                if (pt->core.verbose) {
                    printf("nu_stop_pp_pred=%e nu_stop_pp=%e NU_INT=%d\n ",
                            pt->PP_gamma.nu_stop_pp_gamma_pred,
                            pt->PP_gamma.spec.nu_max,
                            NU_INT);
                }
            }
            if (pt->PP_gamma.spec.j_nu[NU_INT] < pt->core.emiss_lim) {
                //stop = 1;
                pt->PP_gamma.spec.j_nu[NU_INT] = pt->core.emiss_lim;
                pt->PP_gamma.spec.nuFnu_obs[NU_INT] = pt->core.emiss_lim;
                F_nu_pp_obs = pt->core.emiss_lim;
                if (pt->core.verbose) {
                    printf("%e %d\n ", pt->PP_gamma.spec.nu[i], NU_INT);
                }
            }

            if (pt->core.verbose) {
                printf("nuFnu_pp_gamma_obs= %e j=%e nu_stop_pp_pred=%e nu_stop_pp=%e NU_INT=%d\n ",
                pt->PP_gamma.spec.nuFnu_obs[NU_INT], 
                pt->PP_gamma.spec.j_nu[NU_INT],
                pt->PP_gamma.nu_stop_pp_gamma_pred,
                pt->PP_gamma.spec.nu_max,
                NU_INT);
                printf("#-> ********************************\n\n");
            }
            //==========================  END of Loop ove frequencies ====================================
        }
        NU_INT++;
    }

    //Se ancora non ha trovato nu_stop
    if (!stop) {
        pt->PP_gamma.NU_INT_STOP_PP_GAMMA = NU_INT - 1;
    }
    //printf("nu_stop_pp=%e NU_INT_STOP_PP=%d\n", pt->nu_stop_pp, pt->NU_INT_STOP_PP);
    pt->PP_gamma.spec.nu_max_obs = nu_blob_to_nu_obs(pt->PP_gamma.spec.nu_max, pt->core.beam_obj, pt->core.z_cosm);
    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

        FindEpSp(pt->PP_gamma.spec.nu, pt->PP_gamma.spec.nuFnu_obs,   pt->PP_gamma.NU_INT_STOP_PP_GAMMA, pt,
                &(pt->PP_gamma.spec.nu_peak_obs),
                &(pt->PP_gamma.spec.nu_peak_src),
                &(pt->PP_gamma.spec.nu_peak_blob),
                &(pt->PP_gamma.spec.nuFnu_peak_obs),
                &(pt->PP_gamma.spec.nuLnu_peak_src),
                &(pt->PP_gamma.spec.nuLnu_peak_blob));

        if (pt->core.verbose)
        {
            printf("nu_PP_blob peak=%e\n", pt->PP_gamma.spec.nu_peak_blob);
            printf("nu_PP_src   peak=%e\n", pt->PP_gamma.spec.nu_peak_src);
            printf("nu_PP_obs  peak=%e\n", pt->PP_gamma.spec.nu_peak_obs);

            printf("nuFnu PP  blob    peak=%e\n", pt->PP_gamma.spec.nuFnu_peak_obs);
            printf("nuLnu PP  src      peak=%e\n", pt->PP_gamma.spec.nuLnu_peak_src);
            printf("nuLnu PP  obs     peak=%e\n", pt->PP_gamma.spec.nuLnu_peak_blob);
        }
        return;
}
//=========================================================================================

void * eval_j_pp_gamma(void *data){
    struct j_args *thread_args = data;
    unsigned int NU_INT;
    double nu_out;
    for (NU_INT = thread_args->NU_INT_START; NU_INT <= thread_args->NU_INT_STOP; NU_INT++) {
        nu_out=thread_args->nu_array[NU_INT];
        thread_args->blob_pt->PP_gamma.spec.j_nu[NU_INT] = 0.;
       
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->1 in eval_j_pp_gamma   NU_INT=%d   nu_out=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }
        //rate_gamma_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
        //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
        //then you have to multipli by HPLANCK_TeV*nu->TeV*(TeV^-1 cm^-3 s^-1)
        //the you multiply by HPLANCK in  to get erg/( cm^3 Hz s) 
        //and then you divide by 4pi to get erg/( cm^3 Hz s setard) that are j_nu units 

        thread_args->blob_pt->PP_gamma.spec.j_nu[NU_INT] = vluce_cm * thread_args->blob_pt->PP_gamma.NH_pp * bn_to_cm2 *
                        (HPLANCK)* (HPLANCK_TeV * nu_out) *one_by_four_pi* rate_gamma_pp(thread_args->blob_pt,nu_out,-1);
        if (thread_args->blob_pt->core.verbose > 1) {
                 printf("#-> NU_INT=%d j[NU_INT]=%e nu_out=%e  \n",
                            NU_INT,
                            thread_args->blob_pt->PP_gamma.spec.j_nu[NU_INT],
                            nu_out);
        }
    }
    return NULL; 
}
   
