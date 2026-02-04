//==========================================================================
//
//                CALCOLO DELLO SPETTRO DI SINCROTRONE
//
//===========================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file spettro_sincrotrone.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief CALCOLO DELLO SPETTRO DI SINCROTRONE
 *
 */

void spettro_sincrotrone(int Num_file, struct blob * pt) {
    int stop;
    unsigned int NU_INT, I_MAX;
    double nu_src,gmax;
    double nu_p_ext;

    double suggested_nu_stop_Sync;


    //double N_tot_e_Sferic;

    //double tau_nu;

    double F_nu_Sync_obs,S_nu;
    double L_nu_Sync, nuL_nu_Sync;

    //double (*pf_norm) (struct spettro *, double x);


    //char f_Synch[static_file_name_max_legth];

    void *(*eval_j_ptr)(void * args);
    //FILE *fp_Synch;
    if (pt->emitters.Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }

    //*fpe_sinc,*fpf_sinc;
    stop = 0;


    

    //=============================================================================
    //     Starts loop over Synch frequencies and computes Synch stuff
    //=============================================================================
    //initialize the index of pt->Sync.spec.nu[]
    NU_INT = 0;
    I_MAX = pt->core.nu_seed_size -1;
    //nu_B and UB
    pt->Sync.nu_B = (q_esu * pt->core.B) / (2 * pi * me_g * vluce_cm);
    pt->Sync.UB = pow(pt->core.B, 2.0) / (8.0 * pi); /*dens. ener. B */
    
    //Check that  nu_Sync min is consistent with gmax
    //FindNe_NpGp(struct spettro *pt)
    gmax=Find_gmax(pt,pt->emitters.Ne,pt->emitters.griglia_gamma_Ne_log);
    FindNe_NpGp(pt);
    nu_p_ext=pt->emitters.Gamma_p3*pt->emitters.Gamma_p3*pt->Sync.nu_B;
    suggested_nu_stop_Sync = pt->Sync.nu_B * gmax*gmax * pt->Sync.sin_psi * 100.0;
    pt->Sync.spec.nu_max = suggested_nu_stop_Sync;
    if (pt->Sync.spec.nu_max < suggested_nu_stop_Sync) {
        if (pt->core.verbose) {
            printf("!!!!!Warning nu_stop_Sync %e should be at least %e\n",pt->Sync.spec.nu_max,suggested_nu_stop_Sync);
            printf("!!!!!Warning nu_stop_Sync changed to %e\n", suggested_nu_stop_Sync);        
        }
        pt->Sync.spec.nu_max = suggested_nu_stop_Sync;
    }

    


    
    pt->Sync.spec.nu_min_obs =nu_blob_to_nu_obs(pt->Sync.spec.nu_min, pt->core.beam_obj, pt->core.z_cosm);
    pt->Sync.spec.nu_max_obs = nu_blob_to_nu_obs(pt->Sync.spec.nu_max, pt->core.beam_obj, pt->core.z_cosm);

    build_log_grid(pt->Sync.spec.nu_min,  pt->Sync.spec.nu_max, pt->core.nu_seed_size, pt->Sync.spec.nu);
    build_log_grid(pt->Sync.spec.nu_min_obs,  pt->Sync.spec.nu_max_obs, pt->core.nu_seed_size, pt->Sync.spec.nu_obs);




    //========================================================
    // INFORMAZIONI GENERALI
    //========================================================
    if (pt->core.verbose>0) {
        printf("**********************  CALCOLO DELLO SPETTRO DI SINCROTRONE   ****************************\n");
        printf("informazioni generali sul Sync\n");
        printf("nu_B=%e\n", pt->Sync.nu_B);
        printf("gmin*sin_psi=%e\n", pt->emitters.gmin * pt->Sync.sin_psi);
        printf("gmax*sin_psi=%e\n", pt->emitters.gmax * pt->Sync.sin_psi);
        printf("2*gmin*nu_B=%e\n", 2 * pt->emitters.gmin * pt->Sync.nu_B);
        printf("2*gmax*nu_B=%e\n", 2 * pt->emitters.gmax * pt->Sync.nu_B);
        printf("nu_B/(gmin*sin_psi^2)=%e\n", (pt->Sync.nu_B / (pt->emitters.gmin * pt->Sync.sin_psi * pt->Sync.sin_psi)));
        printf("nu_B/(gmax*sin_psi^2)=%e\n", (pt->Sync.nu_B / (pt->emitters.gmax * pt->Sync.sin_psi * pt->Sync.sin_psi)));
        printf("gmin cooling time (s)=%e\n",Sync_tcool(pt,pt->emitters.gmin));
        printf("gmax cooling time (s)=%e\n",Sync_tcool(pt,pt->emitters.gmax));
        printf("gmax from Ne>0 = %e\n", gmax);
        printf("Power_Sync Total From e-=%e\n", Power_Sync_Electron(pt));
        printf("nu_start_Sync=%+-2.20e\n", pt->Sync.spec.nu_min);
        printf("nu_stop_Sync=%+-2.20e\n", pt->Sync.spec.nu_max);
        printf("nu_peak extim =%e\n ",nu_p_ext);
        printf("Number of freq to eval=%d\n", I_MAX);
        printf("out_file=%d\n", pt->core.OUT_FILE);
    }
    //========================================================
    eval_j_ptr = &eval_j_Sync;
    threaded_j_evaluation(pt, eval_j_ptr, pt->Sync.spec.j_nu,pt->Sync.spec.nu,pt->Sync.spec.nu_min, pt->Sync.spec.nu_max,I_MAX,pt->core.N_THREADS);
    for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++) {

    
        if (pt->core.verbose>1) {
            printf("nu=%+-2.20e NU_INT=%d\n", pt->Sync.spec.nu[NU_INT], NU_INT);
            printf("stop=%d\n", stop);
        }

        // Eval j_un and alpha_nu   and set to 0 all the arrays
   

        // Synch is evaluated as int as  nu_start_Sync <nu<nu_stop_Sync
        if ( pt->Sync.spec.nu[NU_INT] <= pt->Sync.spec.nu_max &&  pt->Sync.spec.nu[NU_INT] >= pt->Sync.spec.nu_min && stop != 1) {

            /* erg*s^-1*cm^-3*Hz^-1*sterad^-1 */
            pt->Sync.nu_stop_Sync_ssc =  pt->Sync.spec.nu[NU_INT];
            pt->Sync.NU_INT_STOP_Sync_SSC = NU_INT;
            
        }

        //If j_nu<1e-60 stops the Synch eval
        if ((pt->Sync.spec.j_nu[NU_INT] < pt->core.emiss_lim) && (pt->Sync.spec.nu[NU_INT]>nu_p_ext)){
            stop = 1;
            pt->Sync.spec.n_nu[NU_INT] =0.0;
            pt->Sync.spec.nuFnu_obs[NU_INT]=0.0;
            pt->Sync.spec.nu_max = pt->Sync.spec.nu[NU_INT];
            pt->Sync.spec.nu_max_obs = pt->Sync.spec.nu_obs[NU_INT];
            
        } else if (!stop) {
        	S_nu=solve_S_nu_Sync(pt,NU_INT);

        	//=============================
            //Fluxes transformations and n_Synch
            L_nu_Sync = I_nu_to_L_nu_src(S_nu, pt->core.Surf_region, pt->core.beam_obj); /*erg s^-1  Hz^-1 */
            nu_src = nu_blob_to_nu_src(pt->Sync.spec.nu[NU_INT], pt->core.beam_obj, pt->core.z_cosm);
            nuL_nu_Sync = L_nu_Sync*nu_src; /* erg*s^-1 */
            F_nu_Sync_obs = L_nu_src_to_F_nu(L_nu_Sync, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);
            pt->Sync.spec.nuFnu_obs[NU_INT] = F_nu_Sync_obs*pt->Sync.spec.nu_obs[NU_INT];
            //Gould correction factor
            pt->Sync.spec.n_nu[NU_INT] =pt->core.n_sync_corr_factor*I_nu_to_n(pt->Sync.spec.I_nu[NU_INT], pt->Sync.spec.nu[NU_INT]);
     

            //=============================
            if (pt->core.verbose>1) {
                printf("nuL_nu_Sync=%e\n", nuL_nu_Sync);
            }
        }

       
    }

    //Se ancora non ha trovato nu_stop
    if (!stop) {
        pt->Sync.NU_INT_STOP_Sync_SSC = NU_INT - 1;
    }

    //==========================  END of Loop ove frequencies ====================================
    
   
    //==============================================================
    // se a cauasa di qualche arrotondamento
    // l'ultima nu calcolata nu_Sync[NU_INT-1]
    // e' minore di
    // nu_stop_Sync riaggiorna nu_stop_Sync
    //==============================================================
    if (pt->Sync.spec.nu_max > pt->Sync.spec.nu[NU_INT - 1]) {
        if (pt->core.verbose>1) {
            printf("#-> per arrot. sulla nu_seed_size ho aggioranato\n");
            printf("#-> nu_stop_Sync da=%e a=%e\n", pt->Sync.spec.nu_max, pt->Sync.spec.nu[NU_INT - 1]);
            printf("#-> NU_INT_STOP_Sync_SSC da=%d a=%d\n", pt->Sync.NU_INT_STOP_Sync_SSC, NU_INT - 1);
        }
        pt->Sync.spec.nu_max = pt->Sync.spec.nu[NU_INT - 1];
        pt->Sync.NU_INT_STOP_Sync_SSC = NU_INT - 1;
    }

    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

    FindEpSp(pt->Sync.spec.nu, pt->Sync.spec.nuFnu_obs, pt->Sync.NU_INT_STOP_Sync_SSC, pt,
            &(pt->Sync.spec.nu_peak_obs),
            &(pt->Sync.spec.nu_peak_src),
            &(pt->Sync.spec.nu_peak_blob),
            &(pt->Sync.spec.nuFnu_peak_obs),
            &(pt->Sync.spec.nuLnu_peak_src),
            &(pt->Sync.spec.nuLnu_peak_blob));


    if (pt->core.verbose>0) {
        printf("nu_stop_Sync_ssc =%e NU_INT_STOP_Sync_SSC=%d\n",
            pt->Sync.nu_stop_Sync_ssc, pt->Sync.NU_INT_STOP_Sync_SSC);

        printf("nu_Synch_blob peak=%e\n", pt->Sync.spec.nu_peak_blob);
        printf("nu_Synch_src   peak=%e\n", pt->Sync.spec.nu_peak_src);
        printf("nu_Synch_obs  peak=%e\n", pt->Sync.spec.nu_peak_obs);
    
        printf("nuFnu Synch  blob    peak=%e\n", pt->Sync.spec.nuFnu_peak_obs);
        printf("nuLnu Synch  src      peak=%e\n", pt->Sync.spec.nuLnu_peak_src);
        printf("nuLnu Synch  obs     peak=%e\n", pt->Sync.spec.nuLnu_peak_blob);
    }
    return;
}
//=========================================================================================

void  * eval_j_Sync(void *data){
    unsigned int NU_INT;
	struct j_args *thread_args = data;
    double nu_sync;
    //printf("sono qui, eval_j_SSC \n");
    for (NU_INT = thread_args->NU_INT_START; NU_INT <= thread_args->NU_INT_STOP; NU_INT++) {    
        nu_sync=thread_args->nu_array[NU_INT];
        thread_args->blob_pt->Sync.spec.j_nu[NU_INT] = 0.0;
        thread_args->blob_pt->Sync.alfa_Sync[NU_INT] = 0.0;
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->1 in eval_j_sync   NU_INT=%d   nu_1=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }
        thread_args->blob_pt->Sync.spec.j_nu[NU_INT] = j_nu_Sync(thread_args->blob_pt, nu_sync);
        if (thread_args->blob_pt->core.do_Sync == 2) {
            /* cm^-1 */
            thread_args->blob_pt->Sync.alfa_Sync[NU_INT] = alfa_nu_Sync(thread_args->blob_pt, nu_sync);

        }
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->2 in  eval_j_sync NU_INT=%d j[%d]=%e nu_1=%e \n", NU_INT,NU_INT,
                        thread_args->blob_pt->Sync.spec.j_nu[NU_INT],
                        thread_args->nu_array[NU_INT]);
        }
    }
    //}
    return NULL;   
}