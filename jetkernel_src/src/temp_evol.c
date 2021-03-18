//=================================================================
//  Time Evol
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

void Init_Q_inj(struct temp_ev *pt_ev ){
    unsigned int TMP;
    //printf("alloc Init_Q_inj start \n");
    alloc_temp_ev_array(&(pt_ev->Q_inj), pt_ev->gamma_grid_size);
    alloc_temp_ev_array(&(pt_ev->gamma), pt_ev->gamma_grid_size);
    alloc_temp_ev_array(&(pt_ev->Q_inj_jetset), pt_ev->Q_inj_jetset_gamma_grid_size);
    alloc_temp_ev_array(&(pt_ev->gamma_inj_jetset), pt_ev->Q_inj_jetset_gamma_grid_size);
    for (TMP = 0; TMP < pt_ev->gamma_grid_size; TMP++) {
        pt_ev->gamma[TMP] = 0;
        pt_ev->Q_inj[TMP] = 0;
    }
    for (TMP = 0; TMP < pt_ev->Q_inj_jetset_gamma_grid_size; TMP++) {
        pt_ev->gamma_inj_jetset[TMP] = 0;
        pt_ev->Q_inj_jetset[TMP] = 0;
    }
    //printf("alloc Init_Q_inj stop \n");
}

void Init_temp_evolution(struct blob *pt_spec, struct temp_ev *pt_ev, double luminosity_distance)
{
    //int grid_bounded_to_gamma;

    //unsigned int  gamma_grid_size;

    unsigned int i, TMP;

    double delta_t_eq_t_D, delta_t_eq_t_A, delta_t_eq_t_D_old, delta_t_eq_t_A_old;
    double delta_t_eq_t_DA, delta_t_eq_t_DA_old;


    double delta_log;
    double log_a, log_b;


    pt_ev->t_unit = (pt_spec->R / vluce_cm);
    pt_ev->t_unit_acc = (pt_ev->R_acc / vluce_cm);

    pt_ev->deltat = pt_ev->duration / (double)pt_ev->T_SIZE;

    //----- Diffusion coeff ----
    pt_ev->t_DA0 = pt_ev->t_D0 * 0.5;
    pt_ev->Diff_Coeff = 1.0 / (pt_ev->t_D0);
   
    //----- Acc coeff --------
    pt_ev->Acc_Coeff = 1.0 / (pt_ev->t_A0);
    

    pt_ev->Gamma_Max_Turb_L_max = Larmor_radius_to_gamma(pt_ev->Lambda_max_Turb, pt_spec->B, pt_spec->sin_psi);
    pt_ev->Gamma_Max_Turb_L_coher = Larmor_radius_to_gamma(pt_ev->Lambda_max_Turb * pt_ev->Lambda_choer_Turb_factor, pt_spec->B, pt_spec->sin_psi);
    pt_ev->Diff_coeff_CD = pt_ev->Diff_Coeff * pow((pt_ev->Gamma_Max_Turb_L_coher), pt_ev->Diff_Index);
    pt_ev->Diff_coeff_CA = pt_ev->Diff_Coeff * 2 * pow((pt_ev->Gamma_Max_Turb_L_coher), pt_ev->Diff_Index);

    //----- Esc coeff --------
    pt_ev->T_esc_Coeff_acc = pt_ev->T_esc_Coeff_R_by_c_acc * pt_ev->t_unit_acc;
    pt_ev->T_esc_Coeff_rad = pt_ev->T_esc_Coeff_R_by_c_rad * pt_ev->t_unit;

    pt_ev->gamma_eq_t_D = -1.0;
    pt_ev->gamma_eq_t_DA = -1.0;
    pt_ev->gamma_eq_t_A = -1.0;
    
    delta_t_eq_t_A_old = 1.0;
    delta_t_eq_t_D_old = 1.0;
    delta_t_eq_t_DA_old = 1.0;

    log_a = log10(pt_ev->gmin_griglia);
    log_b = log10(pt_ev->gmax_griglia);
    delta_log = (log_b - log_a) / ((double)static_ev_arr_grid_size - 1);

    for (i = 0; i <= static_ev_arr_grid_size - 1; i++)
    {
        pt_ev->g[i] = pow(10, (log_a + delta_log * (double)(i)));
        pt_ev->t_Sync_cool[i] = Sync_tcool(pt_spec, pt_ev->g[i]);
        pt_ev->t_D[i] = (pt_ev->g[i] * pt_ev->g[i]) / f_Dp(pt_ev->g[i], pt_ev);
        pt_ev->t_DA[i] = 0.5 * pt_ev->t_D[i];
        pt_ev->t_A[i] = pt_ev->g[i] / f_A(pt_ev->g[i], pt_ev);
        pt_ev->t_Esc_acc[i] = f_Tesc(pt_ev->g[i], pt_ev->T_esc_Coeff_acc, pt_ev->Esc_Index);
        pt_ev->t_Esc_rad[i] = f_Tesc(pt_ev->g[i], pt_ev->T_esc_Coeff_rad, pt_ev->Esc_Index);
        delta_t_eq_t_D = pt_ev->t_Sync_cool[i] - pt_ev->t_D[i];
        delta_t_eq_t_A = pt_ev->t_Sync_cool[i] - pt_ev->t_A[i];
        delta_t_eq_t_DA = pt_ev->t_Sync_cool[i] - pt_ev->t_DA[i];
        
        if (i > 0 && pt_ev->g[i] < pt_ev->Gamma_Max_Turb_L_max)
        {
            if (delta_t_eq_t_D <= 0 && delta_t_eq_t_D_old >= 0)
            {
                pt_ev->gamma_eq_t_D = pt_ev->g[i];
            }

            if (delta_t_eq_t_DA <= 0 && delta_t_eq_t_DA_old >= 0)
            {
                pt_ev->gamma_eq_t_DA = pt_ev->g[i];
            }

            if (delta_t_eq_t_A <= 0 && delta_t_eq_t_A_old >= 0)
            {
                pt_ev->gamma_eq_t_A = pt_ev->g[i];
            }


            delta_t_eq_t_A_old = delta_t_eq_t_A;
            delta_t_eq_t_D_old = delta_t_eq_t_D;
            delta_t_eq_t_DA_old = delta_t_eq_t_DA;

        }
    }

    // if gamma_eq still negative, then t_acc never crossed t_cool

    if (delta_t_eq_t_A_old >= 0 && pt_ev->gamma_eq_t_A < 0)
    {
        pt_ev->gamma_eq_t_A = pt_ev->Gamma_Max_Turb_L_max;
    }
    else if (delta_t_eq_t_A_old <= 0 && pt_ev->gamma_eq_t_A < 0)
    {
        pt_ev->gamma_eq_t_A = pt_spec->griglia_gamma_Ne_log[0];
    }
    
    if (delta_t_eq_t_D_old >= 0 && pt_ev->gamma_eq_t_D < 0)
    {
        pt_ev->gamma_eq_t_D = pt_ev->Gamma_Max_Turb_L_max;
    }
    else if (delta_t_eq_t_D_old <= 0 && pt_ev->gamma_eq_t_D < 0)
    {
        pt_ev->gamma_eq_t_D = pt_spec->griglia_gamma_Ne_log[0];
    }

    if (pt_ev->gamma_eq_t_A > pt_ev->Gamma_Max_Turb_L_max)
    {
        pt_ev->gamma_eq_t_A = pt_ev->Gamma_Max_Turb_L_max;
    }
    if (pt_ev->gamma_eq_t_D > pt_ev->Gamma_Max_Turb_L_max)
    {
        pt_ev->gamma_eq_t_D = pt_ev->Gamma_Max_Turb_L_max;
    }

     if (delta_t_eq_t_DA_old >= 0 && pt_ev->gamma_eq_t_DA < 0)
    {
        pt_ev->gamma_eq_t_DA = pt_ev->Gamma_Max_Turb_L_max;
    }
    else if (delta_t_eq_t_DA_old <= 0 && pt_ev->gamma_eq_t_DA < 0)
    {
        pt_ev->gamma_eq_t_DA = pt_spec->griglia_gamma_Ne_log[0];
    }
    if (pt_ev->gamma_eq_t_DA > pt_ev->Gamma_Max_Turb_L_max)
    {
        pt_ev->gamma_eq_t_DA = pt_ev->Gamma_Max_Turb_L_max;
    }

    
    log_a = log10(pt_ev->gmin_griglia);
    log_b = log10(pt_ev->gmax_griglia);
    delta_log = (log_b - log_a) / ((double)pt_ev->gamma_grid_size - 1);
    delta_log = (log_b - log_a) / ((double)pt_ev->gamma_grid_size - 1);
    //-----gamma grid-------
    for (TMP = 0; TMP < pt_ev->gamma_grid_size; TMP++) {
        pt_ev->gamma[TMP] = pow(10, (log_a + delta_log * (double) (TMP)));
    }
    alloc_temp_ev_array(&(pt_ev->T_inj_profile), pt_ev->T_SIZE);
    alloc_temp_ev_array(&(pt_ev->T_acc_profile), pt_ev->T_SIZE);
    for (TMP = 0; TMP < pt_ev->T_SIZE; TMP++) {
        pt_ev->T_inj_profile[TMP]=0.;
        pt_ev->T_acc_profile[TMP]=0.;
    }
    pt_ev->T_COUNTER=0;
 
}

void Run_temp_evolution(struct blob *pt_spec, struct temp_ev *pt_ev, int only_injection) {

        // if luminosity_distance is negative is evaluated internally
        // otherwise the passed value is used

        unsigned int i, E_SIZE, E_N_SIZE, Gamma, T, TMP,NUM_OUT;
        //double Q_scalig_factor;

        double STEP_FILE, COUNT_FILE, OUT_FILE;
        double *x, *N1, *N, *N_escaped;
        double  t;
        //double g, t_D, t_DA, t_A, t_Sync_cool;
        double *xm_p, *xm_m;
        double *dxm_p, *dxm_m, *dxm;
        double K, Cp, Cm, wm_p, wm_m;
        double WP_pm, WM_mm, WP_mm, WM_pm;
        double *A, *B, *C, *R;
        double E_acc_pre,E_acc_post,delta_E_acc,E_acc;
        double Vol_acc,V_rad;
        // *T_inj;
        //int grid_bounded_to_gamma;
        //unsigned int gamma_grid_size;

        
    

        //----SCELTA DELLA VARIABILE ADIMENSIONALE-------------------------
        //    E=mec2*gamma, mec2 rest electron energy
        //    beta=v/c
        //    E_kin=E-mec2=mec2(gamma-1)
        //    x=E_kin/mec2=(g-1)=beta*gamma
        //------------------------------------------------------------------



   

    
    //--------------ALLOCATION OF DYNAMICAL ARRAYS-----------------------
    //--------------SIZE DEFINITION------------------
    
    //---Useful size definition for 0..N-1----------
    E_SIZE = pt_ev->gamma_grid_size;
    E_N_SIZE = pt_ev->gamma_grid_size - 1;
    //------ALLOCATION-------------------------------
    
    //printf("alloc run1 start \n");
    alloc_temp_ev_array(&(pt_ev->N_gamma), pt_ev->NUM_SET * E_SIZE);
    alloc_temp_ev_array(&(pt_ev->N_escaped_gamma), pt_ev->NUM_SET * E_SIZE);
    alloc_temp_ev_array(&(pt_ev->N_time), pt_ev->T_SIZE);
    alloc_temp_ev_array(&(pt_ev->T_esc_acc),E_SIZE);
    alloc_temp_ev_array(&(pt_ev->T_esc_rad),E_SIZE);
    //printf("alloc run1 stop \n");

    
    x = (double *)calloc(E_SIZE, sizeof(double));
    //--------Gamma energia normalizzata-------------
    N1 = (double *) calloc(E_SIZE, sizeof (double));
    N = (double *) calloc(E_SIZE, sizeof (double));
    N_escaped = (double *) calloc(E_SIZE, sizeof (double));
    A = (double *) calloc(E_SIZE, sizeof (double));
    B = (double *) calloc(E_SIZE, sizeof (double));
    C = (double *) calloc(E_SIZE, sizeof (double));
    R = (double *) calloc(E_SIZE, sizeof (double));
    
    //T_inj = (double *) calloc(pt_ev->T_SIZE, sizeof (double));

    xm_p = (double *) calloc(E_SIZE, sizeof (double));
    xm_m = (double *) calloc(E_SIZE, sizeof (double));
    dxm_p = (double *) calloc(E_SIZE, sizeof (double));
    dxm_m = (double *) calloc(E_SIZE, sizeof (double));
    dxm = (double *) calloc(E_SIZE, sizeof (double));
    //--------------------------------------------------------------------



    //----------------Griglia X--------------------------------
    //x=gamma-1
    for (i = 0; i <= E_SIZE - 1; i++) {
        x[i] = pt_ev->gamma[i] - 1;
        //printf("x[%d]=%e\n",i,x[i]);
    }

    for (i = 1; i < E_N_SIZE; i++) {
        xm_p[i] = (x[i + 1] + x[i]) / 2;
        xm_m[i] = (x[i] + x[i - 1]) / 2;
        dxm_p[i] = (x[i + 1] - x[i]);
        dxm_m[i] = (x[i] - x[i - 1]);
        dxm[i] = (x[i + 1] - x[i - 1]) / 2;
        //printf("N_SIZE=%d i-1=%d,i+1=%d\n",N_SIZE,i-1,i+1);
    }
    xm_p[0] = (x[1] + x[0]) / 2;
    xm_m[E_N_SIZE] = (x[E_N_SIZE] + x[E_N_SIZE - 1]) / 2;
    dxm_p[0] = (x[1] - x[0]);
    dxm_m[E_N_SIZE] = (x[E_N_SIZE] - x[E_N_SIZE - 1]);

    //---- x grids -------------
    
    //---------------------------------------------------------
        if (pt_ev->LOG_SET >=1){
        STEP_FILE=log10((double) pt_ev->T_SIZE)/(double) pt_ev->NUM_SET;
        STEP_FILE =  pow(10,STEP_FILE);
    } else{
        STEP_FILE = (double)pt_ev->T_SIZE / (double)pt_ev->NUM_SET;
    }
    COUNT_FILE = STEP_FILE;
    OUT_FILE=-1.0;
    //------------------------------------------


    //-----Injection Term-------
    for (TMP = 0; TMP < pt_ev->gamma_grid_size; TMP++)
    {   
        //iterpolate Q_inj from Q_inj_jetset
        pt_ev->Q_inj[TMP] =  N_distr_interp(pt_ev->Q_inj_jetset_gamma_grid_size,pt_ev->gamma[TMP],pt_ev->gamma_inj_jetset,pt_ev->Q_inj_jetset);
    }

    t = 0;
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        t = t + pt_ev->deltat;
        //T_inj[T] = Inj_temp_prof(t, pt_ev);
        //printf("t=%e,T_inj=%e\n",t,T_inj[T]);
    }
    
    for (TMP = 0; TMP < E_SIZE; TMP++) {
        N_escaped[TMP]=0;
        pt_ev->T_esc_acc[TMP] = f_Tesc(x[TMP], pt_ev->T_esc_Coeff_acc, pt_ev->Esc_Index);
        pt_ev->T_esc_rad[TMP] = f_Tesc(x[TMP], pt_ev->T_esc_Coeff_rad, pt_ev->Esc_Index);
    }

    if (only_injection>0){
        for (TMP = 0; TMP < E_SIZE; TMP++) {
                //iterpolate N from Ne
                N[TMP] =0;                
        }
    }else{
        for (TMP = 0; TMP < E_SIZE; TMP++) {
               //iterpolate N from Ne
               N[TMP] =N_distr_interp(pt_spec->gamma_grid_size,pt_ev->gamma[TMP],pt_spec->griglia_gamma_Ne_log,pt_spec->Ne);
       }
    }
    

    
    //--------------Loop over Gamma------------------------------------------
    t = 0;
    NUM_OUT=0;
    delta_E_acc=0;
    E_acc_pre=0;
    E_acc_post=0;
    E_acc=0;
    Vol_acc=four_by_three_pi*(pt_ev->R_acc*pt_ev->R_acc*pt_ev->R_acc);
    V_rad=four_by_three_pi*(pt_spec->R*pt_spec->R*pt_spec->R);
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        //printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>  temporal evolution step=%d,%e,%e <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n", T, OUT_FILE,COUNT_FILE);
        
        
        
        if (pt_ev->do_Compton_cooling > 0 && T>0) {
            //printf(">>>>>>>>>>>>>>>> Eval Sync \n");            
            spettro_sincrotrone(1, pt_spec);
            spectra_External_Fields(1,pt_spec);
        }
        t = t + pt_ev->deltat;
        pt_ev->T_COUNTER=T;
        //-------G[0]---------------------------

        //--ACC+COOLIG--------------------------
        K = pt_ev->deltat / (dxm_p[0]);
        E_acc_pre=0;
        if (pt_ev->T_acc_profile[T]>0 && E_acc<pt_ev->E_acc_max){
            E_acc_pre=eval_E_acc(pt_ev->gamma, N, E_SIZE, Vol_acc);
           //if (t >= pt_ev->TStart_Acc && t <= pt_ev->TStop_Acc) {
            //printf("Acc+Cool\n");
            Cp = Cfp(xm_p[0], pt_ev) / dxm_p[0];
            wm_p = Bfp(xm_p[0], pt_ev, pt_spec) / (Cfp(xm_p[0], pt_ev)) * dxm_p[0];

            //W^{+}_{m+1/2} W^{-}_{m+1/2}
            Wm(wm_p, &WP_pm, &WM_pm);

            //-------A[0]---------
            A[0] = 0;

            //-------C[0]---------
            C[0] = -1 * K * Cp*WP_pm;

            //-------B[0]----------
            B[0] = 1 + K * Cp * WM_pm + pt_ev->deltat / (pt_ev->T_esc_acc[0]);

            //-------R[0]----------
            R[0] = pt_ev->deltat * pt_ev->Q_inj[0] * pt_ev->T_inj_profile[T] + N[0];
           
            //---------------------------------------
        }           
        //----ONLY COOLING---------------------
        else {
            //printf("Cool\n");
            A[0] = 0;
            B[0] = 1 + K * Cooling(xm_m[0], pt_ev, pt_spec) + pt_ev->deltat / (pt_ev->T_esc_acc[0]);
            C[0] = -1 * K * Cooling(xm_p[0], pt_ev, pt_spec);
            R[0] = pt_ev->deltat * pt_ev->Q_inj[0] * pt_ev->T_inj_profile[T]+ N[0];
        }
        
        for (Gamma = 1; Gamma < E_N_SIZE; Gamma++) {
            //--------COOLING + ACC--------------
            K = pt_ev->deltat / dxm[Gamma];
            if (pt_ev->T_acc_profile[T]>0 && E_acc<pt_ev->E_acc_max){
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
                B[Gamma] = 1 + K * (Cm * WP_mm + Cp * WM_pm) + pt_ev->deltat / (pt_ev->T_esc_acc[Gamma]);
                R[Gamma] = pt_ev->deltat * pt_ev->Q_inj[Gamma] * pt_ev->T_inj_profile[T] + N[Gamma];
            }
            //-----ONLY COOLING-----------------
            else {
                A[Gamma] = 0;
                B[Gamma] = 1 + pt_ev->deltat / (pt_ev->T_esc_acc[Gamma]) + K * Cooling(xm_m[Gamma], pt_ev, pt_spec);
                C[Gamma] = -1 * K * Cooling(xm_p[Gamma], pt_ev, pt_spec);
                R[Gamma] = pt_ev->deltat * pt_ev->Q_inj[Gamma] * pt_ev->T_inj_profile[T] + N[Gamma];
            }

        }

        //--------G[N_SIZE]----------------------
        //--------COOLING + ACC--------------
        K = pt_ev->deltat / (dxm_m[E_N_SIZE]);
        if (pt_ev->T_acc_profile[T]>0  && E_acc<pt_ev->E_acc_max){
        //if (t >= pt_ev->TStart_Acc && t <= pt_ev->TStop_Acc) {
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
            B[E_N_SIZE] = 1 + K * Cm * WP_mm + pt_ev->deltat / (pt_ev->T_esc_rad[E_N_SIZE]);
            //--------R[N_SIZE]------------
            R[E_N_SIZE] = pt_ev->deltat * pt_ev->Q_inj[E_N_SIZE] * pt_ev->T_inj_profile[T]+ N[E_N_SIZE];
        } else {
            A[E_N_SIZE] = 0;
            B[E_N_SIZE] = 1 + pt_ev->deltat / (pt_ev->T_esc_rad[E_N_SIZE]) + K * Cooling(xm_m[E_N_SIZE], pt_ev, pt_spec);
            C[E_N_SIZE] = 0;
            R[E_N_SIZE] = pt_ev->deltat * pt_ev->Q_inj[E_N_SIZE] * pt_ev->T_inj_profile[T]+ N[E_N_SIZE];
        }


        //----------------------------------------

        //--------------END Loop over Gamma---------------------------------------

       
        
        //--------------TRIDIAG SYSTEM SOLVING--------------------
        //printf("********prima di sys****************\n");
        for (TMP = 0; TMP < E_SIZE; TMP++) {
            N_escaped[TMP] = N[TMP]*(1-exp(-pt_ev->deltat/pt_ev->T_esc_rad[TMP]))*(Vol_acc/V_rad);
            N1[TMP] = N[TMP];
            //printf("N=%e T=%d G=%d\n",N1[TMP],T,TMP);
        }
        if (solve_sys1(A, B, C, R, N1, E_SIZE) > 0) {
            printf("errore nella soluzione del sistema condizione di positivita' non soddisfatta\n ");
            exit(1);
        }
        
        
        //printf("********dopo sys********************\n");
        for (TMP = 0; TMP < E_SIZE; TMP++) {
            //iterpolate Ne from N1
            
            N[TMP] = N1[TMP];
        }

         if (pt_ev->T_acc_profile[T]>0 && E_acc<pt_ev->E_acc_max){
            E_acc_post=eval_E_acc(pt_ev->gamma, N, E_SIZE, Vol_acc);
            delta_E_acc=E_acc_post-E_acc_pre;
        }else{
            delta_E_acc=0;
            E_acc_post=0;
        }
        E_acc+=delta_E_acc;

        for (TMP = 0; TMP < pt_spec->gamma_grid_size; TMP++) {
            //iterpolate Ne from N
            pt_spec->Ne[TMP]=N_distr_interp(E_SIZE,pt_spec->griglia_gamma_Ne_log[TMP],pt_ev->gamma,N_escaped);
        }   
        
        //-----Cool Escaped----------
        K = pt_ev->deltat / (dxm_p[0]);
        A[0] = 0;
        B[0] = 1 + K * Cooling(xm_m[0], pt_ev, pt_spec) + pt_ev->deltat / (pt_ev->T_esc_Coeff_rad);
        C[0] = -1 * K * Cooling(xm_p[0], pt_ev, pt_spec);
        R[0] = N_escaped[0];
        for (Gamma = 1; Gamma < E_N_SIZE; Gamma++) {
            K = pt_ev->deltat / (dxm_m[E_N_SIZE]);
            A[Gamma] = 0;
            B[Gamma] = 1 + pt_ev->deltat / (pt_ev->T_esc_Coeff_rad) + K * Cooling(xm_m[Gamma], pt_ev, pt_spec);
            C[Gamma] = -1 * K * Cooling(xm_p[Gamma], pt_ev, pt_spec);
            R[Gamma] = N_escaped[Gamma];
        }
        A[E_N_SIZE] = 0;
        B[E_N_SIZE] = 1 + pt_ev->deltat / (pt_ev->T_esc_Coeff_rad) + K * Cooling(xm_m[E_N_SIZE], pt_ev, pt_spec);
        C[E_N_SIZE] = 0;
        R[E_N_SIZE] =  N_escaped[E_N_SIZE];
        for (TMP = 0; TMP < E_SIZE; TMP++) {
            N1[TMP] = N_escaped[TMP];
            //printf("N=%e T=%d G=%d\n",N1[TMP],T,TMP);
        }
        if (solve_sys1(A, B, C, R, N1, E_SIZE) > 0) {
            printf("errore nella soluzione del sistema condizione di positivita' non soddisfatta\n ");
            exit(1);
        }
        for (TMP = 0; TMP < E_SIZE; TMP++) {
            //iterpolate Ne from N1
            
            N_escaped[TMP] = N1[TMP];
        }

        
        
        //------------- OUT FILE and SED Computations ----------------
        OUT_FILE=(double)T-COUNT_FILE;


        if ((OUT_FILE >= 0) || (T==pt_ev->T_SIZE-1)) {
            
            //printf("-> NUM_OUT=%d T_SIZE=%d T=%d\n",NUM_OUT,pt_ev->T_SIZE,T);
            printf("acc=%e time=%e, delta_E_acc=%e, E_acc=%e, pre=%e, post=%e\n",pt_ev->T_acc_profile[T],t,delta_E_acc,E_acc,E_acc_pre,E_acc_post);
            for (Gamma = 0; Gamma < E_SIZE; Gamma++) {
                
                if (NUM_OUT<pt_ev->NUM_SET){
                    TMP = Gamma + (E_SIZE * NUM_OUT);
                    pt_ev->N_gamma[TMP] = N[Gamma]*(Vol_acc/V_rad)+N_escaped[Gamma];
                    pt_ev->N_escaped_gamma[TMP] = N_escaped[Gamma];
                    pt_ev->N_time[NUM_OUT]=t;
                    
                }
                
                
            }
            
            if (pt_ev->LOG_SET >=1){
                COUNT_FILE*=STEP_FILE;
            } else {
                COUNT_FILE += STEP_FILE;
            }
            NUM_OUT++;
        }
        //---------------------------------------
    }
    //--------------END Loop over Time-----------------------------------------------


  
    //freeing local dynamic arrays
    //printf("freeing local dynamic arrays start \n");
    free(x);
    free(xm_p);
    free(xm_m);
    free(dxm_p);
    free(dxm_m);
    free(dxm);
    free(A);
    free(B);
    free(C);
    free(R);
    free(N1);
    free(N);
    free(N_escaped);
    //free(T_inj);
    //printf("freeing local dynamic arrays stop \n");
    return;
}

double eval_E_acc(double *gamma, double *N, unsigned int gamma_size, double Vol_acc){
    return MEC2 * trapzd_array_arbritary_grid(gamma, N, gamma_size)*Vol_acc;
}


void alloc_temp_ev_array(double ** pt,int size){
        //free if already allocated
        //printf("temp ev pre %p\n",*pt);
        //printf("alloc n\n");
        if (*pt){
            free(*pt);
            //printf("freeing\n");
        }

        *pt = calloc(size, sizeof (double));
        //printf("post %p\n",*pt);
    }


double IntegrandCooolingEquilibrium( struct blob *pt, double gamma_1){
    return N_distr_interp(pt->gamma_grid_size,gamma_1,pt->griglia_gamma_Ne_log,pt->Q_inj_e_second)*exp(pt->gamma_cooling_eq*(1/gamma_1-(1.0/pt->Gamma)));
}


double IntegrateCooolingEquilibrium( struct blob *pt, double gamma, double T_esc ){

    double (*pf_K1) (struct blob * pt, double x);
    double a,b,res,delta;
    unsigned int integ_size;
    pf_K1 = &IntegrandCooolingEquilibrium;
    pt->Gamma=gamma;
    a=gamma;
    b=pt->griglia_gamma_Ne_log[pt->gamma_grid_size-1];
    delta=pt->griglia_gamma_Ne_log[pt->gamma_grid_size-1]-gamma;
    integ_size=delta*1000/(pt->griglia_gamma_Ne_log[pt->gamma_grid_size-1]-pt->griglia_gamma_Ne_log[0]);
    if (integ_size<3){
        integ_size=3;
    }

    res=integrale_trap_log_struct(pf_K1,pt,a,b,integ_size);
    return res*pt->gamma_cooling_eq*T_esc/(gamma*gamma);
}


void CooolingEquilibrium(struct blob * pt, double T_esc){
    //using Eq. 2.26 in Inoue&Takahara
    //http://adsabs.harvard.edu/doi/10.1086/177270
    
    double  a;
    unsigned int ID;
    a=3.0*MEC2/(4.0*vluce_cm*(pt->UB)*SIGTH);
    pt->gamma_cooling_eq=a/T_esc;
    
    for (ID = 0; ID < pt->gamma_grid_size ; ID++){
        
        pt->Ne[ID]=IntegrateCooolingEquilibrium(pt,
                                                pt->griglia_gamma_Ne_log[ID], 
                                                T_esc);
    }


}