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


void spettro_pp(int Num_file, struct spettro *pt) {
    double  k, nu_1, nu_src;
    double L_nu_pp, nuL_nu_pp, F_nu_pp_obs;
    double log_nu, log_nu_start;
    unsigned long l, NU_INT, i, I_MAX, stop;
    char f_pp[static_file_name_max_legth], f_pp_energy[static_file_name_max_legth];
    FILE *fp_pp, *fp_pp_energy;

    //=================================
    // apre i files dove scrive i dati di SSC
    // HEADER FILES
    //=================================
    if (pt->WRITE_TO_FILE==1){
        sprintf(f_pp, "%s%s-pp-gamma.dat",
                pt->path, pt->STEM);

        sprintf(f_pp_energy, "%s%s-pp-gamma-energy.dat",
                pt->path, pt->STEM);

        fp_pp = fopen(f_pp, "w");
        if (fp_pp == NULL) {
            printf("non posso aprire %s\n ", f_pp);
            exit(1);
        }
        fp_pp_energy = fopen(f_pp_energy, "w");
        if (fp_pp == NULL) {
            printf("non posso aprire %s\n ", f_pp);
            exit(1);
        }
        flux_header(fp_pp);
    }
    //==================================================================


    //============================================================
    //         inizio  loop sulle freq per spettro  pp
    //============================================================
    stop = 0;


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!RICODATI DI CAMBIARE la distr e- con p
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // massima e minima freq pp
    pt->nu_stop_pp_pred = pt->gmax * MPC2 / HPLANCK * 100;
    pt->nu_start_pp = E_th_pp * 1E12 * ev_to_erg / HPLANCK / 100;
    pt->nu_start_pp_obs = nu_blob_to_nu_obs(pt->nu_start_pp, pt->beam_obj, pt->z_cosm);
    pt->nu_stop_pp_obs = nu_blob_to_nu_obs(pt->nu_stop_pp_pred, pt->beam_obj, pt->z_cosm);
   
    NU_INT = 0;
    k = (log10(pt->nu_stop_pp_pred) - log10(pt->nu_start_pp));
    I_MAX = pt->nu_IC_size;
    if (pt->verbose)
    {
        printf("**********************  CALCOLO DELLO SPETTRO pp   ****************************\n");

        printf("nu_start_pp=%e nu_stop_pp=%e\n",
               pt->nu_start_pp,
               pt->nu_stop_pp_pred);
        printf("Number of freq to eval=%d\n", I_MAX);
    }
    

    log_nu_start = log10(pt->nu_start_pp);
    for (i = 0; i < I_MAX; i++) {

        // calcola nu con incremento logaritmico
        nu_1 = pow(10, log_nu_start + k * (double) i / (double) I_MAX);
        if (i == 0) {
            // !!! ATTENZIONE POICHE' nu per errori
            // di arrot. numerico puo' avere valore
            // iniziale <nu_start_comp allora impongo
            nu_1 = pt->nu_start_pp;
        }
        pt->nu_1 = nu_1;
        pt->nu_pp[NU_INT]=nu_1;
        pt->nu_pp_obs[NU_INT] = nu_blob_to_nu_obs(nu_1, pt->beam_obj, pt->z_cosm);
       
        if ((nu_1 >= pt->nu_start_pp) && (nu_1 <= pt->nu_stop_pp_pred)) {
            //printf("hi\n");
            if (!stop) {
                //rate_gamma_pp is (dN/dEg)/(c*NH_pp)   TeV^-1 cm^-3 s^-1/(c*NH_pp)
                //you have to multiply by (c*NH_pp) to get dN/dEg (TeV^-1 cm^-3 s^-1)
                //then you have to multipli by HPLANCK_TeV*nu*HPLANCK_TeV*nu*Tev_to_erg/HPLANCK
                //to get erg/(erg cm^3 Hz) that is our units
                pt->j_pp[NU_INT] = vluce_cm * pt->NH_pp * bn_to_cm2 *
                        (HPLANCK_TeV * Tev_to_erg)* (HPLANCK_TeV * pt->nu_1) *
                        rate_gamma_pp(pt);
                if (pt->verbose) {
                    printf("#-> NU_INT=%d j[NU_INT]=%e nu_1=%e i=%d \n",
                            NU_INT,
                            pt->j_pp[NU_INT],
                            nu_1, i);
                }
                nu_src = nu_blob_to_nu_src(nu_1, pt->beam_obj, pt->z_cosm);
                L_nu_pp = j_nu_to_L_nu_src(pt->j_pp[NU_INT], pt->Vol_sphere, pt->beam_obj);
                nuL_nu_pp = L_nu_pp*nu_src;
                F_nu_pp_obs = L_nu_src_to_F_nu(L_nu_pp, pt->beam_obj, pt->z_cosm, pt->dist);
                pt->nuFnu_pp_obs[NU_INT] = F_nu_pp_obs * pt->nu_pp_obs[NU_INT];

                pt->nu_stop_pp = nu_1;
                pt->NU_INT_STOP_PP = NU_INT;
                if (pt->verbose) {
                    printf("nu_stop_pp_pred=%e nu_stop_pp=%e NU_INT=%d\n ",
                            pt->nu_stop_pp_pred,
                            pt->nu_stop_pp,
                            NU_INT);
                }
            }
            if (pt->j_pp[NU_INT] < 1.0e-60) {
                stop = 1;
                pt->j_pp[NU_INT] = 1.0e-60;
                pt->nuFnu_pp_obs[NU_INT] = 1.0e-60;
                F_nu_pp_obs = 1.0e-60;
                if (pt->verbose) {
                    printf("%e %d\n ", nu_1, NU_INT);
                }
            }

            //printf("#-> nu_obs=%e  i=%e\n", pt->nu_pp_obs[NU_INT], pt->nuFnu_pp_obs[NU_INT]);

            //===========================================
            // FILES output nu dnu nuFnu dnuFnu
            //==========================================
            if (pt->WRITE_TO_FILE==1){
                if (!stop) {
                    fprintf(fp_pp, "%4.4e\t%4.4e\t%4.4e\t %4.4e\t%4.4e\t%4.4e\n",
                            log10(pt->nu_pp_obs[NU_INT]),
                            log10(pt->nuFnu_pp_obs[NU_INT]),
                            pt->nu_pp_obs[NU_INT],
                            pt->nuFnu_pp_obs[NU_INT],
                            nu_src,
                            nuL_nu_pp);

                    fprintf(fp_pp_energy, "%4.4e\t%4.4e\n",
                            nu_1*HPLANCK_TeV,
                            pt->j_pp[NU_INT] * nu_1 * erg_to_TeV);
                    // nuL_nu_pp);
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
        pt->NU_INT_STOP_PP = NU_INT - 1;
    }
    //printf("nu_stop_pp=%e NU_INT_STOP_PP=%d\n", pt->nu_stop_pp, pt->NU_INT_STOP_PP);
    pt->nu_stop_pp_obs = nu_blob_to_nu_obs(pt->nu_stop_pp, pt->beam_obj, pt->z_cosm);
    if (pt->WRITE_TO_FILE==1){
        fclose(fp_pp);
        fclose(fp_pp_energy);
    }

    //===========================================
    //    trova nu peak e Flux peak
    //===========================================

        FindEpSp(pt->nu_pp, pt->nuFnu_pp_obs,   pt->NU_INT_STOP_PP, pt,
                &(pt->nu_peak_PP_obs),
                &(pt->nu_peak_PP_src),
                &(pt->nu_peak_PP_blob),
                &(pt->nuFnu_peak_PP_obs),
                &(pt->nuLnu_peak_PP_src),
                &(pt->nuLnu_peak_PP_blob));

        if (pt->verbose)
        {
            printf("nu_PP_blob peak=%e\n", pt->nu_peak_PP_blob);
            printf("nu_PP_src   peak=%e\n", pt->nu_peak_PP_src);
            printf("nu_PP_obs  peak=%e\n", pt->nu_peak_PP_obs);

            printf("nuFnu PP  blob    peak=%e\n", pt->nuFnu_peak_PP_obs);
            printf("nuLnu PP  src      peak=%e\n", pt->nuLnu_peak_PP_src);
            printf("nuLnu PP  obs     peak=%e\n", pt->nuLnu_peak_PP_blob);
        }
        return;
}
//=========================================================================================

