//=================================================================
//  FP solver
//=================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file temp_evol.c
 * \author Andrea Tramacere
 * \date 04-05-2010
 * \brief perform temporal evolution of the plasma
 *
 */



void time_evolve_emitters(struct blob *pt_spec, 
                        struct temp_ev *pt_ev, 
                        int do_inj,
                        double t,
                        unsigned int T,
                        unsigned int E_SIZE,
                        unsigned int E_N_SIZE,
                        double E_acc,
                        double *T_esc,
                        double *N_escped,
                        double *N,
                        double *N_swap,
                        double *A,
                        double *B,
                        double *C,
                        double *R,
                        double *x, 
                        double *xm_p,
                        double *xm_m,
                        double *dxm_p,
                        double *dxm_m,
                        double *dxm){

    // if luminosity_distance is negative is evaluated internally
    // otherwise the passed value is used

    unsigned int TMP,Gamma;
    //double E_acc_pre,E_acc_post,delta_E_acc;
    double K, Cp, Cm, wm_p, wm_m;
    double WP_pm, WM_mm, WP_mm, WM_pm;
        
    
    if (pt_ev->do_Compton_cooling > 0 && T>0) {
        for (TMP = 0; TMP < pt_spec->gamma_grid_size; TMP++) {
            //iterpolate Ne from N
            pt_spec->Ne[TMP]=N_distr_interp(E_SIZE,pt_spec->griglia_gamma_Ne_log[TMP],pt_ev->gamma,N);
        }            
        spettro_sincrotrone(1, pt_spec);
        spectra_External_Fields(1,pt_spec);
    }
    
    

    //--ACC+COOLIG--------------------------
    K = pt_ev->deltat / (dxm_p[0]);
    
    if (pt_ev->T_acc_profile[T]>0 && E_acc<pt_ev->E_acc_max && do_inj==1){
        Cp = Cfp(xm_p[0], pt_ev) / dxm_p[0];
        wm_p = Bfp(xm_p[0], pt_ev, pt_spec) / (Cfp(xm_p[0], pt_ev)) * dxm_p[0];

        //W^{+}_{m+1/2} W^{-}_{m+1/2}
        Wm(wm_p, &WP_pm, &WM_pm);

        //-------A[0]---------
        A[0] = 0;

        //-------C[0]---------
        C[0] = -1 * K * Cp*WP_pm;

        //-------B[0]----------
        B[0] = 1 + K * Cp * WM_pm + pt_ev->deltat / (T_esc[0]);

        //-------R[0]----------
        R[0] =  N[0];
        if (do_inj==1){
            R[0] += pt_ev->deltat * pt_ev->Q_inj[0] * pt_ev->T_inj_profile[T];
        }
        if (do_inj==2){
            R[0] += N_escped[0];
        } 
        //---------------------------------------
    }           
    //----ONLY COOLING---------------------
    else {
        //printf("Cool\n");
        A[0] = 0;
        B[0] = 1 + K * Cooling(xm_m[0], pt_ev, pt_spec) + pt_ev->deltat / (T_esc[0]);
        C[0] = -1 * K * Cooling(xm_p[0], pt_ev, pt_spec);
        R[0] =  N[0];
        if (do_inj==1){
            R[0] += pt_ev->deltat * pt_ev->Q_inj[0] * pt_ev->T_inj_profile[T];
        }
        if (do_inj==2){
            R[0] += N_escped[0];
        } 
    }
    
    for (Gamma = 1; Gamma < E_N_SIZE; Gamma++) {
        //--------COOLING + ACC--------------
        K = pt_ev->deltat / dxm[Gamma];
        if (pt_ev->T_acc_profile[T]>0 && E_acc<pt_ev->E_acc_max && do_inj==1){
        //if (t >= pt_ev->TStart_Acc && t <= pt_ev->TStop_Acc) {
            //C_{m-1/2}
            Cm = Cfp(xm_m[Gamma], pt_ev) / dxm_m[Gamma];

            //C_{m+1/2]
            Cp = Cfp(xm_p[Gamma], pt_ev) / dxm_p[Gamma];

            //w_{m+1/2}
            wm_p = Bfp(xm_p[Gamma], pt_ev, pt_spec) / (Cfp(xm_p[Gamma], pt_ev)) * dxm_p[Gamma];

            //w_{m-1/2}
            wm_m = Bfp(xm_m[Gamma], pt_ev, pt_spec) / (Cfp(xm_m[Gamma], pt_ev)) * dxm_m[Gamma];

            //W^{+}_{m+1/2} W^{-}^{m+1/2}
            Wm(wm_p, &WP_pm, &WM_pm);

            // W^{+}_{m-1/2} W^{-}^{m-1/2}
            Wm(wm_m, &WP_mm, &WM_mm);

            A[Gamma] = -1 * (K * Cm * WM_mm);
            C[Gamma] = -1 * (K * Cp * WP_pm);
            B[Gamma] = 1 + K * (Cm * WP_mm + Cp * WM_pm) + pt_ev->deltat / (T_esc[Gamma]);
            R[Gamma] =  N[Gamma];
            if (do_inj==1){
                R[Gamma] += pt_ev->deltat * pt_ev->Q_inj[Gamma] * pt_ev->T_inj_profile[T];
            }
            if (do_inj==2){
                R[Gamma] += N_escped[Gamma];
            }    
        }
        //-----ONLY COOLING-----------------
        else {
            A[Gamma] = 0;
            B[Gamma] = 1 + pt_ev->deltat / (T_esc[Gamma]) + K * Cooling(xm_m[Gamma], pt_ev, pt_spec);
            C[Gamma] = -1 * K * Cooling(xm_p[Gamma], pt_ev, pt_spec);
            R[Gamma] = N[Gamma];
            if (do_inj==1){
                R[Gamma] += pt_ev->deltat * pt_ev->Q_inj[Gamma] * pt_ev->T_inj_profile[T];
            }
            if (do_inj==2){
                R[Gamma] += N_escped[Gamma];
            }   
        }

    }

    //--------COOLING + ACC--------------
    K = pt_ev->deltat / (dxm_m[E_N_SIZE]);
    if (pt_ev->T_acc_profile[T]>0  && E_acc<pt_ev->E_acc_max && do_inj==1){
    
        Cm = Cfp(xm_m[E_N_SIZE], pt_ev) / dxm_m[E_N_SIZE];
        wm_m = Bfp(xm_m[E_N_SIZE], pt_ev, pt_spec) / (Cfp(xm_m[E_N_SIZE], pt_ev)) * dxm_m[E_N_SIZE];

        //   W^{+}_{m-1/2} W^{-}^{m-1/2}
        Wm(wm_m, &WP_mm, &WM_mm);

        //WP_mm=Wm_m*exp(wm_p*0.5);

        //--------A[N_SIZE]------------
        A[E_N_SIZE] = -1 * K * Cm*WM_mm;
        //--------C[N_SIZE]------------
        C[E_N_SIZE] = 0;
        //--------B[N_SIZE]------------
        B[E_N_SIZE] = 1 + K * Cm * WP_mm + pt_ev->deltat / (T_esc[E_N_SIZE]);
        //--------R[N_SIZE]------------
        R[E_N_SIZE] =  N[E_N_SIZE];
        if (do_inj==1){
            R[E_N_SIZE] += pt_ev->deltat * pt_ev->Q_inj[E_N_SIZE] * pt_ev->T_inj_profile[T];
        } 
        if (do_inj==2){
            R[E_N_SIZE] += N_escped[E_N_SIZE];
        }  
    } else {
        A[E_N_SIZE] = 0;
        B[E_N_SIZE] = 1 + pt_ev->deltat / (T_esc[E_N_SIZE]) + K * Cooling(xm_m[E_N_SIZE], pt_ev, pt_spec);
        C[E_N_SIZE] = 0;
        R[E_N_SIZE] = N[E_N_SIZE];
        if (do_inj==1){
            R[E_N_SIZE] += pt_ev->deltat * pt_ev->Q_inj[E_N_SIZE] * pt_ev->T_inj_profile[T];
        }
        if (do_inj==2){
            R[E_N_SIZE] += N_escped[E_N_SIZE];
        }
    }

    //--------------TRIDIAG SYSTEM SOLVING--------------------
    for (TMP = 0; TMP < E_SIZE; TMP++) {
        N_swap[TMP] = N[TMP];
    }
    if (solve_sys1(A, B, C, R, N_swap, E_SIZE) > 0) {
        printf("errore nella soluzione del sistema condizione di positivita' non soddisfatta\n ");
        exit(1);
    }
    for (TMP = 0; TMP < E_SIZE; TMP++) {
        //iterpolate Ne from N1
        N[TMP] = N_swap[TMP];
    }

}   
                