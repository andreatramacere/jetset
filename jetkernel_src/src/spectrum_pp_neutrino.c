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
    gmax=Find_gmax(pt,pt->emitters.Np,pt->emitters.griglia_gamma_Np_log);
    pt->PP_neutrino.nu_stop_pp_neutrino_pred = gmax * MPC2 / HPLANCK * 100;
    pt->PP_neutrino.spec_tot.nu_min = E_th_pp * 1E12 * ev_to_erg / HPLANCK / 100;
    pt->PP_neutrino.spec_tot.nu_min_obs = nu_blob_to_nu_obs(pt->PP_neutrino.spec_tot.nu_min, pt->core.beam_obj, pt->core.z_cosm);
    pt->PP_neutrino.spec_tot.nu_max_obs = nu_blob_to_nu_obs(pt->PP_neutrino.nu_stop_pp_neutrino_pred, pt->core.beam_obj, pt->core.z_cosm);
   
    NU_INT = 0;
   
    build_log_grid(pt->PP_neutrino.spec_tot.nu_min,  pt->PP_neutrino.nu_stop_pp_neutrino_pred, pt->core.nu_IC_size, pt->PP_neutrino.spec_tot.nu);
    build_log_grid(pt->PP_neutrino.spec_tot.nu_min_obs,  pt->PP_neutrino.spec_tot.nu_max_obs, pt->core.nu_IC_size, pt->PP_neutrino.spec_tot.nu_obs);

    build_log_grid(pt->PP_neutrino.spec_tot.nu_min,  pt->PP_neutrino.nu_stop_pp_neutrino_pred, pt->core.nu_IC_size, pt->PP_neutrino.spec_mu.nu);
    build_log_grid(pt->PP_neutrino.spec_tot.nu_min_obs,  pt->PP_neutrino.spec_tot.nu_max_obs, pt->core.nu_IC_size, pt->PP_neutrino.spec_mu.nu_obs);

    build_log_grid(pt->PP_neutrino.spec_tot.nu_min,  pt->PP_neutrino.nu_stop_pp_neutrino_pred, pt->core.nu_IC_size, pt->PP_neutrino.spec_e.nu);
    build_log_grid(pt->PP_neutrino.spec_tot.nu_min_obs,  pt->PP_neutrino.spec_tot.nu_max_obs, pt->core.nu_IC_size, pt->PP_neutrino.spec_e.nu_obs);
    I_MAX = pt->core.nu_IC_size -1;

    eval_j_neutrio_ptr = &eval_j_pp_neutrino;
    pt->PP_gamma.pp_racc_nu_mu=rate_neutrino_mu_1_pp(pt ,pt->PP_neutrino.spec_tot.nu_min,1);
    threaded_j_evaluation(pt, eval_j_neutrio_ptr, pt->PP_neutrino.spec_tot.j_nu,pt->PP_neutrino.spec_tot.nu,pt->PP_neutrino.spec_tot.nu_min, pt->PP_neutrino.nu_stop_pp_neutrino_pred,I_MAX,pt->core.N_THREADS);

    
    if (pt->core.verbose)
    {
        printf("**********************  CALCOLO DELLO SPETTRO pp   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->PP_neutrino.spec_tot.nu_min,
               pt->PP_neutrino.nu_stop_pp_neutrino_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }
    
    for (i = 0; i <= I_MAX; i++) {        
        if ((pt->PP_neutrino.spec_tot.nu[NU_INT] >= pt->PP_neutrino.spec_tot.nu_min) && (pt->PP_neutrino.spec_tot.nu[NU_INT] <= pt->PP_neutrino.nu_stop_pp_neutrino_pred)) {
            //printf("hi\n");
            if (!stop) {
                //tot
                //nu_src = nu_blob_to_nu_src(nu_1, pt->core.beam_obj, pt->core.z_cosm);
                L_nu_pp = j_nu_to_L_nu_src(pt->PP_neutrino.spec_tot.j_nu[NU_INT], pt->core.Vol_region, pt->core.beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);
                
                pt->PP_neutrino.spec_tot.nuFnu_obs[NU_INT] = F_nu_pp_obs * pt->PP_neutrino.spec_tot.nu_obs[NU_INT];
                
                //mu
                L_nu_pp = j_nu_to_L_nu_src(pt->PP_neutrino.spec_mu.j_nu[NU_INT], pt->core.Vol_region, pt->core.beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);

                pt->PP_neutrino.spec_mu.nuFnu_obs[NU_INT] = F_nu_pp_obs * pt->PP_neutrino.spec_mu.nu_obs[NU_INT];

                //e-
                L_nu_pp = j_nu_to_L_nu_src(pt->PP_neutrino.spec_e.j_nu[NU_INT], pt->core.Vol_region, pt->core.beam_obj);
                //nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);

                pt->PP_neutrino.spec_e.nuFnu_obs[NU_INT] = F_nu_pp_obs * pt->PP_neutrino.spec_e.nu_obs[NU_INT];


                pt->PP_neutrino.spec_tot.nu_max = pt->PP_neutrino.spec_tot.nu[NU_INT];
                pt->PP_neutrino.NU_INT_STOP_PP_NUETRINO = NU_INT;
                if (pt->core.verbose) {
                    printf("nu_stop_pp_pred=%e nu_stop_pp=%e NU_INT=%d\n ",
                            pt->PP_neutrino.nu_stop_pp_neutrino_pred,
                            pt->PP_neutrino.spec_tot.nu_max,
                            NU_INT);
                }
            }
            if (pt->PP_neutrino.spec_tot.j_nu[NU_INT] <pt->core.emiss_lim) {
                //stop = 1;
                pt->PP_neutrino.spec_tot.j_nu[NU_INT] = pt->core.emiss_lim;
                pt->PP_neutrino.spec_tot.nuFnu_obs[NU_INT] = pt->core.emiss_lim;

                pt->PP_neutrino.spec_mu.j_nu[NU_INT] = pt->core.emiss_lim;
                pt->PP_neutrino.spec_mu.nuFnu_obs[NU_INT] = pt->core.emiss_lim;

                pt->PP_neutrino.spec_e.j_nu[NU_INT] = pt->core.emiss_lim;
                pt->PP_neutrino.spec_e.nuFnu_obs[NU_INT] = pt->core.emiss_lim;

                F_nu_pp_obs = pt->core.emiss_lim;
                
                if (pt->core.verbose) {
                    printf("%e %d\n ", pt->PP_neutrino.spec_tot.nu[NU_INT], NU_INT);
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
        pt->PP_neutrino.NU_INT_STOP_PP_NUETRINO = NU_INT - 1;
    }
    pt->PP_neutrino.spec_tot.nu_max_obs = nu_blob_to_nu_obs(pt->PP_neutrino.spec_tot.nu_max, pt->core.beam_obj, pt->core.z_cosm);
  
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

        FindEpSp(pt->PP_neutrino.spec_tot.nu, pt->PP_neutrino.spec_tot.nuFnu_obs,   pt->PP_neutrino.NU_INT_STOP_PP_NUETRINO, pt,
                &(pt->PP_neutrino.spec_tot.nu_peak_obs),
                &(pt->PP_neutrino.spec_tot.nu_peak_src),
                &(pt->PP_neutrino.spec_tot.nu_peak_blob),
                &(pt->PP_neutrino.spec_tot.nuFnu_peak_obs),
                &(pt->PP_neutrino.spec_tot.nuLnu_peak_src),
                &(pt->PP_neutrino.spec_tot.nuLnu_peak_blob));

        if (pt->core.verbose)
        {
            printf("nu_PP_blob peak=%e\n", pt->PP_neutrino.spec_tot.nu_peak_blob);
            printf("nu_PP_src   peak=%e\n", pt->PP_neutrino.spec_tot.nu_peak_src);
            printf("nu_PP_obs  peak=%e\n", pt->PP_neutrino.spec_tot.nu_peak_obs);

            printf("nuFnu PP  blob    peak=%e\n", pt->PP_neutrino.spec_tot.nuFnu_peak_obs);
            printf("nuLnu PP  src      peak=%e\n", pt->PP_neutrino.spec_tot.nuLnu_peak_src);
            printf("nuLnu PP  obs     peak=%e\n", pt->PP_neutrino.spec_tot.nuLnu_peak_blob);
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
        thread_args->blob_pt->PP_neutrino.spec_e.j_nu[NU_INT] = 0.;
        thread_args->blob_pt->PP_neutrino.spec_tot.j_nu[NU_INT] = 0.;
        thread_args->blob_pt->PP_neutrino.spec_mu.j_nu[NU_INT] = 0.;

       
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->1 in eval_j_pp_neutrino   NU_INT=%d   nu_out=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }
        //rate_neutrino_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
        //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
        //then you have to multipli by HPLANCK_TeV*nu->TeV*(TeV^-1 cm^-3 s^-1)
        //the you multiply by HPLANCK in  to get erg/( cm^3 Hz s) that is our units
        //and then you divide by 4pi to get erg/( cm^3 Hz s setard) that are j_nu units 

        //contribution from Eq. 66 neutirno_mu_1 Kenler 2006
        j_neutrino_mu_1 = rate_neutrino_mu_1_pp(thread_args->blob_pt ,nu_out,-1)*vluce_cm * thread_args->blob_pt->PP_gamma.NH_pp * bn_to_cm2 *
                (HPLANCK)* (HPLANCK_TeV * nu_out)*one_by_four_pi;
                
        //contribution from Eq. 62  neutirno_mu_2 and neutirno_2 Kenler 2006
        //assuminf F_neutrino_e = F_e and F_neutirno_mu_2~F_neutrino_e

        //F_neutrino_e  is obteined from  injetcted e-
        //vluce_cm * pt->PP_gamma.NH_pp * bn_to_cm2 alredy done in N_distr
        gamma_e=nu_out*HPLANCK/MEC2;
        thread_args->blob_pt->PP_neutrino.spec_e.j_nu[NU_INT]= HPLANCK*gamma_e*N_distr_interp(thread_args->blob_pt->emitters.gamma_grid_size, gamma_e, thread_args->blob_pt->emitters.griglia_gamma_Ne_log, thread_args->blob_pt->emitters.Q_inj_e_second)*one_by_four_pi;
        
        // F_neutirno_mu_2~F_neutrino_e
        j_neutrino_mu_2=thread_args->blob_pt->PP_neutrino.spec_e.j_nu[NU_INT];
        
        j_neutrino_tot=j_neutrino_mu_1+j_neutrino_mu_2+thread_args->blob_pt->PP_neutrino.spec_e.j_nu[NU_INT];
        
        thread_args->blob_pt->PP_neutrino.spec_tot.j_nu[NU_INT]=j_neutrino_tot;
        thread_args->blob_pt->PP_neutrino.spec_mu.j_nu[NU_INT]=j_neutrino_mu_1+j_neutrino_mu_2;
    
        if (thread_args->blob_pt->core.verbose > 1) {
                 printf("#-> NU_INT=%d j[NU_INT]=%e nu_out=%e  \n",
                            NU_INT,
                            thread_args->blob_pt->PP_neutrino.spec_tot.j_nu[NU_INT],
                            nu_out);
        }
    }
    return NULL; 
}