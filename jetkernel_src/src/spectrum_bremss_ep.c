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
 * \brief CALCOLO DELLO SPETTRO Bremsstrahlung
 *
 */


void spettro_bremss_ep(int Num_file, struct blob *pt) {
    double  nu_check;
    double L_nu_bremss_ep, F_nu_bremss_ep_obs;
    double gmax;
    unsigned int NU_INT, i, I_MAX, stop;
    void *(*eval_j_bremss)(void * args);

    //============================================================
    //         inizio  loop sulle freq per spettro  pp
    //============================================================
    stop = 0;



    // massima e minima freq bremss
    gmax=Find_gmax(pt,pt->Ne,pt->griglia_gamma_Ne_log);
    pt->nu_stop_bremss_ep_pred = gmax*MEC2/HPLANCK*10;
    pt->nu_start_bremss_ep = pt->gmin_griglia*MEC2/HPLANCK/100/10;
    pt->nu_start_bremss_ep_obs = nu_blob_to_nu_obs(pt->nu_start_bremss_ep, pt->beam_obj, pt->z_cosm);
    pt->nu_stop_bremss_ep_obs = nu_blob_to_nu_obs(pt->nu_stop_bremss_ep_pred, pt->beam_obj, pt->z_cosm);
    nu_check=(pt->nu_start_bremss_ep)*0.5;
    NU_INT = 0;
    I_MAX = pt->nu_IC_size -1;
    
    build_log_grid(pt->nu_start_bremss_ep,  pt->nu_stop_bremss_ep_pred, pt->nu_IC_size, pt->nu_bremss_ep);
    build_log_grid(pt->nu_start_bremss_ep_obs,  pt->nu_stop_bremss_ep_obs, pt->nu_IC_size, pt->nu_bremss_ep_obs);

    eval_j_bremss = &eval_j_pp_bremss_ep;
    threaded_j_evaluation(pt, eval_j_bremss, pt->j_bremss_ep,pt->nu_bremss_ep,pt->nu_start_bremss_ep, pt->nu_stop_bremss_ep_pred,I_MAX,pt->N_THREADS);
    if (pt->verbose)
    {
        printf("**********************  CALCOLO DELLO SPETTRO bremss ep   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->nu_start_bremss_ep,
               pt->nu_stop_bremss_ep_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }
    
    for (i = 0; i <= I_MAX; i++) {
       
        if ((pt->nu_bremss_ep[i] >= pt->nu_start_bremss_ep) && (pt->nu_bremss_ep[i] <= pt->nu_stop_bremss_ep_pred)) {
            //printf("hi\n");
            if (!stop) {
                
                //nu_src = nu_blob_to_nu_src(nu_1, pt->beam_obj, pt->z_cosm);
                L_nu_bremss_ep = j_nu_to_L_nu_src(pt->j_bremss_ep[NU_INT], pt->Vol_region, pt->beam_obj);
                //nuL_nu_ep_brem = L_nu_bremss_ep*nu_src;
                F_nu_bremss_ep_obs = L_nu_src_to_F_nu(L_nu_bremss_ep, pt->beam_obj, pt->z_cosm, pt->dist);
                pt->nuFnu_bremss_ep_obs[NU_INT] = F_nu_bremss_ep_obs * pt->nu_bremss_ep_obs[NU_INT];

                pt->nu_stop_bremss_ep = pt->nu_bremss_ep[i];
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
                if (pt->nu_bremss_ep[i]>nu_check){
                    stop = 1;

                    F_nu_bremss_ep_obs = pt->emiss_lim;
                }
                if (pt->verbose) {
                    printf("%e %d\n ", pt->nu_bremss_ep[i], NU_INT);
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

void * eval_j_pp_bremss_ep(void *data){
    unsigned int NU_INT;
    double nu_out;
    struct j_args *thread_args = data;
    for (NU_INT = thread_args->NU_INT_START; NU_INT <= thread_args->NU_INT_STOP; NU_INT++) {
        nu_out=thread_args->nu_array[NU_INT];
        thread_args->blob_pt->j_bremss_ep[NU_INT] = 0.;
       
        if (thread_args->blob_pt->verbose > 1) {
                printf("#->1 in eval_j_pp_bremss_ep   NU_INT=%d   nu_out=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }
        //rate_gamma_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
        //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
        //then you have to multipli by HPLANCK_TeV*nu->TeV*(TeV^-1 cm^-3 s^-1)
        //the you multiply by HPLANCK in  to get erg/( cm^3 Hz s) 
        //and then you divide by 4pi to get erg/( cm^3 Hz s setard) that are j_nu units 

                
        thread_args->blob_pt->j_bremss_ep[NU_INT] = thread_args->blob_pt->NH_pp*j_nu_bremss_ep(thread_args->blob_pt,thread_args->blob_pt->nu_bremss_ep[NU_INT]);

        if (thread_args->blob_pt->verbose > 1) {
                 printf("#-> NU_INT=%d j[NU_INT]=%e nu_out=%e  \n",
                            NU_INT,
                            thread_args->blob_pt->j_bremss_ep[NU_INT],
                            nu_out);
        }
    }
    return NULL; 
}