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


    //FILE *fp_Synch;
    if (pt->Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }

    //*fpe_sinc,*fpf_sinc;
    stop = 0;


    

    //=============================================================================
    //     Starts loop over Synch frequencies and computes Synch stuff
    //=============================================================================
    //initialize the index of pt->nu_Sync[]
    NU_INT = 0;
    I_MAX = pt->nu_seed_size-1;
    //nu_B and UB
    pt->nu_B = (q_esu * pt->B) / (2 * pi * me_g * vluce_cm);
    pt->UB = pow(pt->B, 2.0) / (8.0 * pi); /*dens. ener. B */
    
    //Check that  nu_Sync min is consistent with gmax
    //FindNe_NpGp(struct spettro *pt)
    gmax=Find_gmax(pt,pt->Ne,pt->griglia_gamma_Ne_log);
    FindNe_NpGp(pt);
    nu_p_ext=pt->Gamma_p3*pt->Gamma_p3*pt->nu_B;
    suggested_nu_stop_Sync = pt->nu_B * gmax*gmax * pt->sin_psi * 100.0;
    pt->nu_stop_Sync = suggested_nu_stop_Sync;
    if (pt->nu_stop_Sync < suggested_nu_stop_Sync) {
        if (pt->verbose) {
            printf("!!!!!Warning nu_stop_Sync %e should be at least %e\n",pt->nu_stop_Sync,suggested_nu_stop_Sync);
            printf("!!!!!Warning nu_stop_Sync changed to %e\n", suggested_nu_stop_Sync);        
        }
        pt->nu_stop_Sync = suggested_nu_stop_Sync;
    }

    


    
    pt->nu_start_Sync_obs =nu_blob_to_nu_obs(pt->nu_start_Sync, pt->beam_obj, pt->z_cosm);
    pt->nu_stop_Sync_obs = nu_blob_to_nu_obs(pt->nu_stop_Sync, pt->beam_obj, pt->z_cosm);

    build_log_grid(pt->nu_start_Sync,  pt->nu_stop_Sync, pt->nu_seed_size, pt->nu_Sync);
    build_log_grid(pt->nu_start_Sync_obs,  pt->nu_stop_Sync_obs, pt->nu_seed_size, pt->nu_Sync_obs);




    //========================================================
    // INFORMAZIONI GENERALI
    //========================================================
    if (pt->verbose>0) {
        printf("**********************  CALCOLO DELLO SPETTRO DI SINCROTRONE   ****************************\n");
        printf("informazioni generali sul Sync\n");
        printf("nu_B=%e\n", pt->nu_B);
        printf("gmin*sin_psi=%e\n", pt->gmin * pt->sin_psi);
        printf("gmax*sin_psi=%e\n", pt->gmax * pt->sin_psi);
        printf("2*gmin*nu_B=%e\n", 2 * pt->gmin * pt->nu_B);
        printf("2*gmax*nu_B=%e\n", 2 * pt->gmax * pt->nu_B);
        printf("nu_B/(gmin*sin_psi^2)=%e\n", (pt->nu_B / (pt->gmin * pt->sin_psi * pt->sin_psi)));
        printf("nu_B/(gmax*sin_psi^2)=%e\n", (pt->nu_B / (pt->gmax * pt->sin_psi * pt->sin_psi)));
        printf("gmin cooling time (s)=%e\n",Sync_tcool(pt,pt->gmin));
        printf("gmax cooling time (s)=%e\n",Sync_tcool(pt,pt->gmax));
        printf("gmax from Ne>0 = %e\n", gmax);
        printf("Power_Sync Total From e-=%e\n", Power_Sync_Electron(pt));
        printf("nu_start_Sync=%+-2.20e\n", pt->nu_start_Sync);
        printf("nu_stop_Sync=%+-2.20e\n", pt->nu_stop_Sync);
        printf("nu_peak extim =%e\n ",nu_p_ext);
        printf("Number of freq to eval=%d\n", I_MAX);
        printf("out_file=%d\n", pt->OUT_FILE);
    }
    //========================================================


    for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++) {

        // !!!  ATT! NON TOGLIERE QUESTA ISTRUZIONE
        pt->nu = pt->nu_Sync[NU_INT];

        if (pt->verbose>1) {
            printf("nu=%+-2.20e NU_INT=%d\n", pt->nu_Sync[NU_INT], NU_INT);
            printf("stop=%d\n", stop);
        }

        // Eval j_un and alpha_nu   and set to 0 all the arrays
        pt->j_Sync[NU_INT] = 0.0;
        pt->alfa_Sync[NU_INT] = 0.0;

        // Synch is evaluated as int as  nu_start_Sync <nu<nu_stop_Sync
        if ( pt->nu_Sync[NU_INT] <= pt->nu_stop_Sync &&  pt->nu_Sync[NU_INT] >= pt->nu_start_Sync && stop != 1) {

            /* erg*s^-1*cm^-3*Hz^-1*sterad^-1 */
            pt->j_Sync[NU_INT] = j_nu_Sync(pt);

            if (pt->do_Sync == 2) {
                /* cm^-1 */
                pt->alfa_Sync[NU_INT] = alfa_nu_Sync(pt);

            }

            // Update the maximum frequency of the Synch
            // to use in the IC computation
            pt->nu_stop_Sync_ssc =  pt->nu_Sync[NU_INT];
            pt->NU_INT_STOP_Sync_SSC = NU_INT;
            
        }

        //If j_nu<1e-60 stops the Synch eval
        if ((pt->j_Sync[NU_INT] < pt->emiss_lim) && (pt->nu_Sync[NU_INT]>nu_p_ext)){
            stop = 1;
            pt->n_Sync[NU_INT] =0.0;
            pt->nuF_nu_Sync_obs[NU_INT]=0.0;
            pt->nu_stop_Sync = pt->nu_Sync[NU_INT];
            pt->nu_stop_Sync_obs = pt->nu_Sync_obs[NU_INT];
            
        } else if (!stop) {
        	S_nu=solve_S_nu_Sync(pt,NU_INT);

        	//=============================
            //Fluxes transformations and n_Synch
            L_nu_Sync = I_nu_to_L_nu_src(S_nu, pt->Surf_sphere, pt->beam_obj); /*erg s^-1  Hz^-1 */
            nu_src = nu_blob_to_nu_src(pt->nu_Sync[NU_INT], pt->beam_obj, pt->z_cosm);
            nuL_nu_Sync = L_nu_Sync*nu_src; /* erg*s^-1 */
            F_nu_Sync_obs = L_nu_src_to_F_nu(L_nu_Sync, pt->beam_obj, pt->z_cosm, pt->dist);
            pt->nuF_nu_Sync_obs[NU_INT] = F_nu_Sync_obs*pt->nu_Sync_obs[NU_INT];
            //Gould correction factor
            pt->n_Sync[NU_INT] =0.75*I_nu_to_n(pt->I_nu_Sync[NU_INT], pt->nu_Sync[NU_INT]);
     

            //=============================
            if (pt->verbose>1) {
                printf("nuL_nu_Sync=%e\n", nuL_nu_Sync);
            }
        }

       
    }

    //Se ancora non ha trovato nu_stop
    if (!stop) {
        pt->NU_INT_STOP_Sync_SSC = NU_INT - 1;
    }

    //==========================  END of Loop ove frequencies ====================================
    
   
    //==============================================================
    // se a cauasa di qualche arrotondamento
    // l'ultima nu calcolata nu_Sync[NU_INT-1]
    // e' minore di
    // nu_stop_Sync riaggiorna nu_stop_Sync
    //==============================================================
    if (pt->nu_stop_Sync > pt->nu_Sync[NU_INT - 1]) {
        if (pt->verbose>1) {
            printf("#-> per arrot. sulla nu_seed_size ho aggioranato\n");
            printf("#-> nu_stop_Sync da=%e a=%e\n", pt->nu_stop_Sync, pt->nu_Sync[NU_INT - 1]);
            printf("#-> NU_INT_STOP_Sync_SSC da=%d a=%d\n", pt->NU_INT_STOP_Sync_SSC, NU_INT - 1);
        }
        pt->nu_stop_Sync = pt->nu_Sync[NU_INT - 1];
        pt->NU_INT_STOP_Sync_SSC = NU_INT - 1;
    }

    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

    FindEpSp(pt->nu_Sync, pt->nuF_nu_Sync_obs, pt->NU_INT_STOP_Sync_SSC, pt,
            &(pt->nu_peak_Sync_obs),
            &(pt->nu_peak_Sync_src),
            &(pt->nu_peak_Sync_blob),
            &(pt->nuFnu_peak_Sync_obs),
            &(pt->nuLnu_peak_Sync_src),
            &(pt->nuLnu_peak_Sync_blob));


    if (pt->verbose>0) {
        printf("nu_stop_Sync_ssc =%e NU_INT_STOP_Sync_SSC=%d\n",
            pt->nu_stop_Sync_ssc, pt->NU_INT_STOP_Sync_SSC);

        printf("nu_Synch_blob peak=%e\n", pt->nu_peak_Sync_blob);
        printf("nu_Synch_src   peak=%e\n", pt->nu_peak_Sync_src);
        printf("nu_Synch_obs  peak=%e\n", pt->nu_peak_Sync_obs);
    
        printf("nuFnu Synch  blob    peak=%e\n", pt->nuFnu_peak_Sync_obs);
        printf("nuLnu Synch  src      peak=%e\n", pt->nuLnu_peak_Sync_src);
        printf("nuLnu Synch  obs     peak=%e\n", pt->nuLnu_peak_Sync_blob);
    }
    return;
}
//=========================================================================================
