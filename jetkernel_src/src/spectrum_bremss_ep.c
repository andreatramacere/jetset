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
    gmax=Find_gmax(pt,pt->emitters.Ne,pt->emitters.griglia_gamma_Ne_log);
    pt->Bremss_ep.nu_stop_bremss_ep_pred = gmax*MEC2/HPLANCK*10;
    pt->Bremss_ep.spec.nu_min = pt->emitters.gmin_griglia*MEC2/HPLANCK/100/10;
    pt->Bremss_ep.spec.nu_min_obs = nu_blob_to_nu_obs(pt->Bremss_ep.spec.nu_min, pt->core.beam_obj, pt->core.z_cosm);
    pt->Bremss_ep.spec.nu_max_obs = nu_blob_to_nu_obs(pt->Bremss_ep.nu_stop_bremss_ep_pred, pt->core.beam_obj, pt->core.z_cosm);
    nu_check=(pt->Bremss_ep.spec.nu_min)*0.5;
    NU_INT = 0;
    I_MAX = pt->core.nu_IC_size -1;
    
    build_log_grid(pt->Bremss_ep.spec.nu_min,  pt->Bremss_ep.nu_stop_bremss_ep_pred, pt->core.nu_IC_size, pt->Bremss_ep.spec.nu);
    build_log_grid(pt->Bremss_ep.spec.nu_min_obs,  pt->Bremss_ep.spec.nu_max_obs, pt->core.nu_IC_size, pt->Bremss_ep.spec.nu_obs);

    eval_j_bremss = &eval_j_pp_bremss_ep;
    threaded_j_evaluation(pt, eval_j_bremss, pt->Bremss_ep.spec.j_nu,pt->Bremss_ep.spec.nu,pt->Bremss_ep.spec.nu_min, pt->Bremss_ep.nu_stop_bremss_ep_pred,I_MAX,pt->core.N_THREADS);
    if (pt->core.verbose)
    {
        printf("**********************  CALCOLO DELLO SPETTRO bremss ep   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->Bremss_ep.spec.nu_min,
               pt->Bremss_ep.nu_stop_bremss_ep_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }
    
    for (i = 0; i <= I_MAX; i++) {
       
        if ((pt->Bremss_ep.spec.nu[i] >= pt->Bremss_ep.spec.nu_min) && (pt->Bremss_ep.spec.nu[i] <= pt->Bremss_ep.nu_stop_bremss_ep_pred)) {
            //printf("hi\n");
            if (!stop) {
                
                //nu_src = nu_blob_to_nu_src(nu_1, pt->core.beam_obj, pt->core.z_cosm);
                L_nu_bremss_ep = j_nu_to_L_nu_src(pt->Bremss_ep.spec.j_nu[NU_INT], pt->core.Vol_region, pt->core.beam_obj);
                //nuL_nu_ep_brem = L_nu_bremss_ep*nu_src;
                F_nu_bremss_ep_obs = L_nu_src_to_F_nu(L_nu_bremss_ep, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);
                pt->Bremss_ep.spec.nuFnu_obs[NU_INT] = F_nu_bremss_ep_obs * pt->Bremss_ep.spec.nu_obs[NU_INT];

                pt->Bremss_ep.spec.nu_max = pt->Bremss_ep.spec.nu[i];
                pt->Bremss_ep.NU_INT_STOP_BREMSS_EP = NU_INT;
                if (pt->core.verbose) {
                    printf("nu_stop_brems_ep_pred=%e nu_stop_bremss_ep=%e NU_INT=%d\n ",
                            pt->Bremss_ep.nu_stop_bremss_ep_pred,
                            pt->Bremss_ep.spec.nu_max,
                            NU_INT);
                }
            }
            if (pt->Bremss_ep.spec.j_nu[NU_INT] < pt->core.emiss_lim) {
                pt->Bremss_ep.spec.j_nu[NU_INT] = pt->core.emiss_lim;
                pt->Bremss_ep.spec.nuFnu_obs[NU_INT] = pt->core.emiss_lim;
                if (pt->Bremss_ep.spec.nu[i]>nu_check){
                    stop = 1;

                    F_nu_bremss_ep_obs = pt->core.emiss_lim;
                }
                if (pt->core.verbose) {
                    printf("%e %d\n ", pt->Bremss_ep.spec.nu[i], NU_INT);
                }
            }

            if (pt->core.verbose) {
                printf("#-> ********************************\n\n");
            }
            //==========================  END of Loop ove frequencies ====================================
        }
        NU_INT++;
    }

    //Se ancora non ha trovato nu_stop
    if (!stop) {
        pt->Bremss_ep.NU_INT_STOP_BREMSS_EP = NU_INT - 1;
    }
    //printf("nu_stop_pp=%e NU_INT_STOP_PP=%d\n", pt->nu_stop_pp, pt->NU_INT_STOP_PP);
    pt->Bremss_ep.spec.nu_min_obs = nu_blob_to_nu_obs(pt->Bremss_ep.spec.nu_min, pt->core.beam_obj, pt->core.z_cosm);
    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

        FindEpSp(pt->Bremss_ep.spec.nu_obs, pt->Bremss_ep.spec.nuFnu_obs,   pt->Bremss_ep.NU_INT_STOP_BREMSS_EP, pt,
                &(pt->Bremss_ep.spec.nu_peak_obs),
                &(pt->Bremss_ep.spec.nu_peak_src),
                &(pt->Bremss_ep.spec.nu_peak_blob),
                &(pt->Bremss_ep.spec.nuFnu_peak_obs),
                &(pt->Bremss_ep.spec.nuLnu_peak_src),
                &(pt->Bremss_ep.spec.nuLnu_peak_blob));

        if (pt->core.verbose)
        {
            printf("nu_bremss_ep_blob peak=%e\n", pt->Bremss_ep.spec.nu_peak_blob);
            printf("nu_bremss_ep_src   peak=%e\n", pt->Bremss_ep.spec.nu_peak_src);
            printf("nu_bremss_ep_obs  peak=%e\n", pt->Bremss_ep.spec.nu_peak_src);

            printf("nuFnu bremss_ep  blob    peak=%e\n", pt->Bremss_ep.spec.nuFnu_peak_obs);
            printf("nuLnu bremss_ep  src      peak=%e\n", pt->Bremss_ep.spec.nuLnu_peak_src);
            printf("nuLnu bremss_ep  obs     peak=%e\n", pt->Bremss_ep.spec.nuLnu_peak_blob);
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
        thread_args->blob_pt->Bremss_ep.spec.j_nu[NU_INT] = 0.;
       
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->1 in eval_j_pp_bremss_ep   NU_INT=%d   nu_out=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }
        //rate_gamma_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
        //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
        //then you have to multipli by HPLANCK_TeV*nu->TeV*(TeV^-1 cm^-3 s^-1)
        //the you multiply by HPLANCK in  to get erg/( cm^3 Hz s) 
        //and then you divide by 4pi to get erg/( cm^3 Hz s setard) that are j_nu units 

                
        thread_args->blob_pt->Bremss_ep.spec.j_nu[NU_INT] = thread_args->blob_pt->PP_gamma.NH_pp*j_nu_bremss_ep(thread_args->blob_pt,thread_args->blob_pt->Bremss_ep.spec.nu[NU_INT]);

        if (thread_args->blob_pt->core.verbose > 1) {
                 printf("#-> NU_INT=%d j[NU_INT]=%e nu_out=%e  \n",
                            NU_INT,
                            thread_args->blob_pt->Bremss_ep.spec.j_nu[NU_INT],
                            nu_out);
        }
    }
    return NULL; 
}