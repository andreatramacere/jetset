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


void spettro_pp_neutrino(int Num_file, struct blob *pt) {
    double gmax;
    double L_nu_pp, F_nu_pp_obs;
    unsigned int NU_INT, i, I_MAX, stop;
    double gamma_e,j_neutrino_mu_1,j_neutrino_mu_2,j_neutrino_tot;
    void *(*eval_j_neutrio_ptr)(void * args);

    stop = 0;


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!RICODATI DI CAMBIARE la distr e- con p
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // massima e minima freq pp
    gmax=Find_gmax(pt,pt->Np,pt->griglia_gamma_Np_log);
    pt->nu_stop_pp_neutrino_pred = gmax * MPC2 / HPLANCK * 100;
    pt->nu_start_pp_neutrino = E_th_pp * 1E12 * ev_to_erg / HPLANCK / 100;
    pt->nu_start_pp_neutrino_obs = nu_blob_to_nu_obs(pt->nu_start_pp_neutrino, pt->beam_obj, pt->z_cosm);
    pt->nu_stop_pp_neutrino_obs = nu_blob_to_nu_obs(pt->nu_stop_pp_neutrino_pred, pt->beam_obj, pt->z_cosm);
   
    NU_INT = 0;
   
    build_log_grid(pt->nu_start_pp_neutrino,  pt->nu_stop_pp_neutrino_pred, pt->nu_IC_size, pt->nu_pp_neutrino_tot);
    build_log_grid(pt->nu_start_pp_neutrino_obs,  pt->nu_stop_pp_neutrino_obs, pt->nu_IC_size, pt->nu_pp_neutrino_tot_obs);

    build_log_grid(pt->nu_start_pp_neutrino,  pt->nu_stop_pp_neutrino_pred, pt->nu_IC_size, pt->nu_pp_neutrino_mu);
    build_log_grid(pt->nu_start_pp_neutrino_obs,  pt->nu_stop_pp_neutrino_obs, pt->nu_IC_size, pt->nu_pp_neutrino_mu_obs);

    build_log_grid(pt->nu_start_pp_neutrino,  pt->nu_stop_pp_neutrino_pred, pt->nu_IC_size, pt->nu_pp_neutrino_e);
    build_log_grid(pt->nu_start_pp_neutrino_obs,  pt->nu_stop_pp_neutrino_obs, pt->nu_IC_size, pt->nu_pp_neutrino_e_obs);
    I_MAX = pt->nu_IC_size -1;

    eval_j_neutrio_ptr = &eval_j_pp_neutrino;
    pt->pp_racc_nu_mu=rate_neutrino_mu_1_pp(pt ,pt->nu_start_pp_neutrino,1);
    threaded_j_evaluation(pt, eval_j_neutrio_ptr, pt->j_pp_neutrino_tot,pt->nu_pp_neutrino_tot,pt->nu_start_pp_neutrino, pt->nu_stop_pp_neutrino_pred,I_MAX,pt->N_THREADS);

    
    if (pt->verbose)
    {
        printf("**********************  CALCOLO DELLO SPETTRO pp   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->nu_start_pp_neutrino,
               pt->nu_stop_pp_neutrino_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }
    
    for (i = 0; i <= I_MAX; i++) {        
        if ((pt->nu_pp_neutrino_tot[NU_INT] >= pt->nu_start_pp_neutrino) && (pt->nu_pp_neutrino_tot[NU_INT] <= pt->nu_stop_pp_neutrino_pred)) {
            //printf("hi\n");
            if (!stop) {
                //tot
                //nu_src = nu_blob_to_nu_src(nu_1, pt->beam_obj, pt->z_cosm);
                L_nu_pp = j_nu_to_L_nu_src(pt->j_pp_neutrino_tot[NU_INT], pt->Vol_region, pt->beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->beam_obj, pt->z_cosm, pt->dist);
                
                pt->nuFnu_pp_neutrino_tot_obs[NU_INT] = F_nu_pp_obs * pt->nu_pp_neutrino_tot_obs[NU_INT];
                
                //mu
                L_nu_pp = j_nu_to_L_nu_src(pt->j_pp_neutrino_mu[NU_INT], pt->Vol_region, pt->beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->beam_obj, pt->z_cosm, pt->dist);

                pt->nuFnu_pp_neutrino_mu_obs[NU_INT] = F_nu_pp_obs * pt->nu_pp_neutrino_mu_obs[NU_INT];

                //e-
                L_nu_pp = j_nu_to_L_nu_src(pt->j_pp_neutrino_e[NU_INT], pt->Vol_region, pt->beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->beam_obj, pt->z_cosm, pt->dist);

                pt->nuFnu_pp_neutrino_e_obs[NU_INT] = F_nu_pp_obs * pt->nu_pp_neutrino_e_obs[NU_INT];


                pt->nu_stop_pp_neutrino = pt->nu_pp_neutrino_tot[NU_INT];
                pt->NU_INT_STOP_PP_NUETRINO = NU_INT;
                if (pt->verbose) {
                    printf("nu_stop_pp_pred=%e nu_stop_pp=%e NU_INT=%d\n ",
                            pt->nu_stop_pp_neutrino_pred,
                            pt->nu_stop_pp_neutrino,
                            NU_INT);
                }
            }
            if (pt->j_pp_neutrino_tot[NU_INT] <pt->emiss_lim) {
                //stop = 1;
                pt->j_pp_neutrino_tot[NU_INT] = pt->emiss_lim;
                pt->nuFnu_pp_neutrino_tot_obs[NU_INT] = pt->emiss_lim;

                pt->j_pp_neutrino_mu[NU_INT] = pt->emiss_lim;
                pt->nuFnu_pp_neutrino_mu_obs[NU_INT] = pt->emiss_lim;

                pt->j_pp_neutrino_e[NU_INT] = pt->emiss_lim;
                pt->nuFnu_pp_neutrino_e_obs[NU_INT] = pt->emiss_lim;

                F_nu_pp_obs = pt->emiss_lim;
                
                if (pt->verbose) {
                    printf("%e %d\n ", pt->nu_pp_neutrino_tot[NU_INT], NU_INT);
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
        pt->NU_INT_STOP_PP_NUETRINO = NU_INT - 1;
    }
    pt->nu_stop_pp_neutrino_obs = nu_blob_to_nu_obs(pt->nu_stop_pp_neutrino, pt->beam_obj, pt->z_cosm);
  
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

        FindEpSp(pt->nu_pp_neutrino_tot, pt->nuFnu_pp_neutrino_tot_obs,   pt->NU_INT_STOP_PP_NUETRINO, pt,
                &(pt->nu_peak_PP_neutrino_obs),
                &(pt->nu_peak_PP_neutrino_src),
                &(pt->nu_peak_PP_neutrino_blob),
                &(pt->nuFnu_peak_PP_neutrino_obs),
                &(pt->nuLnu_peak_PP_neutrino_src),
                &(pt->nuLnu_peak_PP_neutrino_blob));

        if (pt->verbose)
        {
            printf("nu_PP_blob peak=%e\n", pt->nu_peak_PP_neutrino_blob);
            printf("nu_PP_src   peak=%e\n", pt->nu_peak_PP_neutrino_src);
            printf("nu_PP_obs  peak=%e\n", pt->nu_peak_PP_neutrino_obs);

            printf("nuFnu PP  blob    peak=%e\n", pt->nuFnu_peak_PP_neutrino_obs);
            printf("nuLnu PP  src      peak=%e\n", pt->nuLnu_peak_PP_neutrino_src);
            printf("nuLnu PP  obs     peak=%e\n", pt->nuLnu_peak_PP_neutrino_blob);
        }
        return;
}
//=========================================================================================

void  * eval_j_pp_neutrino(void *data){
    unsigned int NU_INT;
    struct j_args *thread_args = data;
    double nu_out,j_neutrino_mu_1,j_neutrino_mu_2,j_neutrino_tot,gamma_e;

    for (NU_INT = thread_args->NU_INT_START; NU_INT <= thread_args->NU_INT_STOP; NU_INT++) {
        nu_out=thread_args->nu_array[NU_INT];
        thread_args->blob_pt->j_pp_neutrino_e[NU_INT] = 0.;
        thread_args->blob_pt->j_pp_neutrino_tot[NU_INT] = 0.;
        thread_args->blob_pt->j_pp_neutrino_mu[NU_INT] = 0.;

       
        if (thread_args->blob_pt->verbose > 1) {
                printf("#->1 in eval_j_pp_neutrino   NU_INT=%d   nu_out=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }
        //rate_neutrino_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
        //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
        //then you have to multipli by HPLANCK_TeV*nu->TeV*(TeV^-1 cm^-3 s^-1)
        //the you multiply by HPLANCK in  to get erg/( cm^3 Hz s) that is our units
        //and then you divide by 4pi to get erg/( cm^3 Hz s setard) that are j_nu units 

        //contribution from Eq. 66 neutirno_mu_1 Kenler 2006
        j_neutrino_mu_1 = rate_neutrino_mu_1_pp(thread_args->blob_pt ,nu_out,-1)*vluce_cm * thread_args->blob_pt->NH_pp * bn_to_cm2 *
                (HPLANCK)* (HPLANCK_TeV * nu_out)*one_by_four_pi;
                
        //contribution from Eq. 62  neutirno_mu_2 and neutirno_2 Kenler 2006
        //assuminf F_neutrino_e = F_e and F_neutirno_mu_2~F_neutrino_e

        //F_neutrino_e  is obteined from  injetcted e-
        //vluce_cm * pt->NH_pp * bn_to_cm2 alredy done in N_distr
        gamma_e=nu_out*HPLANCK/MEC2;
        thread_args->blob_pt->j_pp_neutrino_e[NU_INT]= HPLANCK*gamma_e*N_distr_interp(thread_args->blob_pt->gamma_grid_size, gamma_e, thread_args->blob_pt->griglia_gamma_Ne_log, thread_args->blob_pt->Q_inj_e_second)*one_by_four_pi;
        
        // F_neutirno_mu_2~F_neutrino_e
        j_neutrino_mu_2=thread_args->blob_pt->j_pp_neutrino_e[NU_INT];
        
        j_neutrino_tot=j_neutrino_mu_1+j_neutrino_mu_2+thread_args->blob_pt->j_pp_neutrino_e[NU_INT];
        
        thread_args->blob_pt->j_pp_neutrino_tot[NU_INT]=j_neutrino_tot;
        thread_args->blob_pt->j_pp_neutrino_mu[NU_INT]=j_neutrino_mu_1+j_neutrino_mu_2;
    
        if (thread_args->blob_pt->verbose > 1) {
                 printf("#-> NU_INT=%d j[NU_INT]=%e nu_out=%e  \n",
                            NU_INT,
                            thread_args->blob_pt->j_pp_neutrino_tot[NU_INT],
                            nu_out);
        }
    }
    return NULL; 
}