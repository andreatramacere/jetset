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

void Init_temp_evolution(struct blob *pt_spec_rad, struct blob *pt_spec_acc, struct temp_ev *pt_ev, double luminosity_distance)
{
    //int grid_bounded_to_gamma;

    //unsigned int  gamma_grid_size;

    unsigned int i, TMP;

    double delta_t_eq_t_D, delta_t_eq_t_A, delta_t_eq_t_D_old, delta_t_eq_t_A_old;
    double delta_t_eq_t_DA, delta_t_eq_t_DA_old;


    double delta_log;
    double log_a, log_b;


    pt_ev->t_unit_rad = (pt_ev->R_rad_start / vluce_cm);
    pt_ev->t_unit_acc = (pt_ev->Delta_R_acc / vluce_cm);

    pt_ev->deltat = pt_ev->duration / (double)pt_ev->T_SIZE;

    //----- Diffusion coeff ----
    pt_ev->t_DA0 = pt_ev->t_D0 * 0.5;
    pt_ev->Diff_Coeff = 1.0 / (pt_ev->t_D0);
   
    //----- Acc coeff --------
    pt_ev->Acc_Coeff = 1.0 / (pt_ev->t_A0);
    

    pt_ev->Gamma_Max_Turb_L_max = Larmor_radius_to_gamma(pt_ev->Lambda_max_Turb, pt_spec_acc->B, pt_spec_acc->sin_psi);
    pt_ev->Gamma_Max_Turb_L_coher = Larmor_radius_to_gamma(pt_ev->Lambda_max_Turb * pt_ev->Lambda_choer_Turb_factor, pt_spec_acc->B, pt_spec_acc->sin_psi);
    pt_ev->Diff_coeff_CD = pt_ev->Diff_Coeff * pow((pt_ev->Gamma_Max_Turb_L_coher), pt_ev->Diff_Index);
    pt_ev->Diff_coeff_CA = pt_ev->Diff_Coeff * 2 * pow((pt_ev->Gamma_Max_Turb_L_coher), pt_ev->Diff_Index);

    //----- Esc coeff --------
    pt_ev->T_esc_Coeff_acc = pt_ev->T_esc_Coeff_R_by_c_acc * pt_ev->t_unit_acc;
    pt_ev->T_esc_Coeff_rad = pt_ev->T_esc_Coeff_R_by_c_rad * pt_ev->t_unit_rad;

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
        pt_ev->t_Sync_cool[i] = Sync_tcool(pt_spec_acc, pt_ev->g[i]);
        pt_ev->t_D[i] = (pt_ev->g[i] * pt_ev->g[i]) / f_Dp(pt_ev->g[i], pt_ev);
        pt_ev->t_DA[i] = 0.5 * pt_ev->t_D[i];
        pt_ev->t_A[i] = pt_ev->g[i] / f_A(pt_ev->g[i], pt_ev);
        pt_ev->t_Esc_acc[i] = f_Tesc(pt_ev->g[i], pt_ev->T_esc_Coeff_acc, pt_ev->Esc_Index_acc);
        pt_ev->t_Esc_rad[i] = f_Tesc(pt_ev->g[i], pt_ev->T_esc_Coeff_rad, pt_ev->Esc_Index_rad);
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
        pt_ev->gamma_eq_t_A = pt_spec_acc->griglia_gamma_Ne_log[0];
    }
    
    if (delta_t_eq_t_D_old >= 0 && pt_ev->gamma_eq_t_D < 0)
    {
        pt_ev->gamma_eq_t_D = pt_ev->Gamma_Max_Turb_L_max;
    }
    else if (delta_t_eq_t_D_old <= 0 && pt_ev->gamma_eq_t_D < 0)
    {
        pt_ev->gamma_eq_t_D = pt_spec_acc->griglia_gamma_Ne_log[0];
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
        pt_ev->gamma_eq_t_DA = pt_spec_acc->griglia_gamma_Ne_log[0];
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
    expansion_profile_pre_run(pt_spec_rad,pt_ev);
 }

void expansion_profile_pre_run(struct blob *pt_spec, struct temp_ev *pt_ev){
    unsigned int i;
    double deltat,t;
    //t_spec->B=pt_ev->B_rad;
    //pt_spec->R=pt_ev->R_rad_start;
   
    deltat = pt_ev->duration / (double)static_ev_arr_grid_size;
    pt_spec->beta_Gamma=eval_beta_gamma(pt_spec->BulkFactor);
    t=0;
    //pt_ev->theta_exp_rad=pt_ev->R_jet_exp*Deg_to_Rad;
    //pt_ev->theta_exp_AR=tan(pt_ev->theta_exp_rad);
    for (i = 0; i <= static_ev_arr_grid_size - 1; i++){
        pt_ev->time_blob_exp[i]=t;
        pt_ev->R_H_t_pre[i]=eval_R_H_jet_t( pt_spec,  pt_ev, t);
        if (pt_ev->do_Expansion ==1){                
                pt_ev->R_t_pre[i]=eval_R_jet_t(pt_spec,  pt_ev, pt_ev->time_blob_exp[i]);
                pt_ev->B_t_pre[i]=eval_B_jet_t(pt_spec,  pt_ev, pt_ev->R_t_pre[i], pt_ev->time_blob_exp[i]);
                //printf("R_H_jet_exp=%e R_H=%e R=%e \n",pt_ev->R_H_jet_exp,pt_ev->R_H_t_pre[i],pt_ev->R_t_pre[i]);
            }else{        
            pt_ev->R_t_pre[i]=pt_ev->R_rad_start;
            pt_ev->B_t_pre[i]=pt_ev->B_rad;
            
        }
    t=t+deltat;
    }
}

void Run_temp_evolution(struct blob *pt_spec_rad, struct blob *pt_spec_acc, struct temp_ev *pt_ev, int only_injection, int do_injection) {

    // if luminosity_distance is negative is evaluated internally
    // otherwise the passed value is used

    unsigned int i, E_SIZE, E_N_SIZE, Gamma, T, TMP,NUM_OUT;
    //double Q_scalig_factor;

    double STEP_FILE, COUNT_FILE, OUT_FILE;
    double *x, *N_swap, *N_acc, *N_rad, *N_escaped;
    double  t;
    //double g, t_D, t_DA, t_A, t_Sync_cool;
    double *xm_p, *xm_m;
    double *dxm_p, *dxm_m, *dxm;
    //double K, Cp, Cm, wm_p, wm_m;
    //double WP_pm, WM_mm, WP_mm, WM_pm;
    double *A, *B, *C, *R;
    double E_acc_pre,E_acc_post,delta_E_acc,E_acc;
    double Vol_acc, Vol_rad;
    double exp_factor;
    double t_ad,t_ad_esc,t_esc;
    

    
    //--------------ALLOCATION OF DYNAMICAL ARRAYS-----------------------
    //--------------SIZE DEFINITION------------------
    
    //---Useful size definition for 0..N-1----------
    E_SIZE = pt_ev->gamma_grid_size;
    E_N_SIZE = pt_ev->gamma_grid_size - 1;
    //------ALLOCATION-------------------------------
    
    //printf("alloc run1 start \n");
    //alloc_temp_ev_array(&(pt_ev->N_gamma), pt_ev->NUM_SET * E_SIZE);
    alloc_temp_ev_array(&(pt_ev->N_acc_gamma), pt_ev->NUM_SET * E_SIZE);
    alloc_temp_ev_array(&(pt_ev->N_rad_gamma), pt_ev->NUM_SET * E_SIZE);
    alloc_temp_ev_array(&(pt_ev->N_time), pt_ev->NUM_SET);
    alloc_temp_ev_array(&(pt_ev->T_esc_acc),E_SIZE);
    alloc_temp_ev_array(&(pt_ev->T_esc_rad),E_SIZE);
    //printf("alloc run1 stop \n");

    
    x = (double *)calloc(E_SIZE, sizeof(double));
    //--------Gamma energia normalizzata-------------
    N_swap = (double *) calloc(E_SIZE, sizeof (double));
    N_acc = (double *) calloc(E_SIZE, sizeof (double));
    N_rad = (double *) calloc(E_SIZE, sizeof (double));
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



    t = 0;
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        t = t + pt_ev->deltat;
        //T_inj[T] = Inj_temp_prof(t, pt_ev);
        //printf("t=%e,T_inj=%e\n",t,T_inj[T]);
    }
    
    for (TMP = 0; TMP < E_SIZE; TMP++) {
        N_rad[TMP]=0;
        pt_ev->T_esc_acc[TMP] = f_Tesc(x[TMP], pt_ev->T_esc_Coeff_acc, pt_ev->Esc_Index_acc);
        pt_ev->T_esc_rad[TMP] = f_Tesc(x[TMP], pt_ev->T_esc_Coeff_rad, pt_ev->Esc_Index_rad);
        N_escaped[TMP]=0;
        N_acc[TMP]=0;
    }

    //if (only_injection<10){
    //    for (TMP = 0; TMP < E_SIZE; TMP++) {
    //            //iterpolate N from Ne
    //            N_acc[TMP] =0;                
    //    }
    //}else{
    

    
    t = 0;
    pt_ev->t=t;
    NUM_OUT=0;
    delta_E_acc=0;
    E_acc_pre=0;
    E_acc_post=0;
    E_acc=0;
    Vol_acc=pi*(pt_ev->R_rad_start*pt_ev->R_rad_start*pt_ev->Delta_R_acc);
    Vol_rad=four_by_three_pi*(pt_ev->R_rad_start*pt_ev->R_rad_start*pt_ev->R_rad_start);
    pt_spec_rad->B=pt_ev->B_rad;
    pt_spec_rad->R=pt_ev->R_rad_start;
    pt_spec_acc->B=pt_ev->B_acc;
    InitRadiative(pt_spec_acc);
    InitRadiative(pt_spec_rad);
    pt_ev->R_H_jet_t=pt_ev->R_H_rad_start;
    pt_ev->R_jet_t=pt_ev->R_rad_start;
    exp_factor=1.0;
    t_ad_esc=-1.0;
    //pt_ev->theta_exp_rad=pt_ev->R_jet_exp*Deg_to_Rad;
    //pt_ev->theta_exp_AR=tan(pt_ev->theta_exp_rad);
    
    if (only_injection<1){
        for (TMP = 0; TMP < E_SIZE; TMP++) {
               //iterpolate N from Ne
               N_rad[TMP] =N_distr_interp(pt_spec_rad->gamma_grid_size,pt_ev->gamma[TMP],pt_spec_rad->griglia_gamma_Ne_log,pt_spec_rad->Ne);
               //*Vol_rad;
       }
    }
    

    //-----Injection Term-------
    for (TMP = 0; TMP < pt_ev->gamma_grid_size; TMP++)
    {   
        //iterpolate Q_inj from Q_inj_jetset
        pt_ev->Q_inj[TMP] =  N_distr_interp(pt_ev->Q_inj_jetset_gamma_grid_size,pt_ev->gamma[TMP],pt_ev->gamma_inj_jetset,pt_ev->Q_inj_jetset);
        //*Vol_acc;
    }
    
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        
        
        pt_ev->T_COUNTER=T;
        pt_ev->t=t;
        E_acc_pre=eval_E_acc(pt_ev->gamma, N_acc, E_SIZE, Vol_acc);
        
        // Evolve ACC Region
        if (do_injection>0){
            
            if (do_injection==1){
                time_evolve_emitters(pt_spec_acc,pt_ev,1,t,T,E_SIZE,E_N_SIZE,E_acc,pt_ev->T_esc_acc,N_escaped,N_acc,N_swap,A,B,C,R,x,xm_p,xm_m,dxm_p,dxm_m,dxm);
                //--------------UPDATE ACC ENERGY--------------------
                if (pt_ev->T_acc_profile[T]>0 && E_acc<pt_ev->E_acc_max){
                    E_acc_post=eval_E_acc(pt_ev->gamma, N_acc, E_SIZE, Vol_acc);
                    delta_E_acc=E_acc_post-E_acc_pre;
                }else{
                    delta_E_acc=0;
                    E_acc_post=0;
                }
                E_acc+=delta_E_acc;
                //--- Inj from ACC To Radiative
                for (TMP = 0; TMP < E_SIZE; TMP++) { 
                    N_escaped[TMP] = N_acc[TMP]*(1-exp(-pt_ev->deltat/pt_ev->T_esc_acc[TMP]))*Vol_acc/pt_spec_rad->Vol_sphere;
                }
            }else{
                for (TMP = 0; TMP < E_SIZE; TMP++) { 
                    N_escaped[TMP] = pt_ev->deltat * pt_ev->Q_inj[TMP] * pt_ev->T_inj_profile[T]*pt_spec_acc->Vol_sphere/pt_spec_rad->Vol_sphere;
                    N_acc[TMP]=N_escaped[TMP];
                }
            }
         }  
        
        
        
        //--------------UPDATE T_ACC PROFILE--------------------
        //if energy provided to electron exceeds limit, then switch off acceleration
        if (pt_ev->T_acc_profile[T]>0 && E_acc>=pt_ev->E_acc_max){
            pt_ev->T_acc_profile[T]=0;
        }
        E_acc+=delta_E_acc;
        // Evolve Rad Region
        if (pt_ev->do_Expansion ==1 ){
            
            exp_factor=update_jet_expansion(pt_spec_rad,pt_ev,t);
            Vol_rad=four_by_three_pi*(pt_ev->R_jet_t*pt_ev->R_jet_t*pt_ev->R_jet_t);     
            t_ad=Adiabatic_Cooling_time(pt_ev, pt_spec_rad, pt_ev->R_jet_t)/(3);
            pt_ev->T_esc_Coeff_rad = pt_ev->T_esc_Coeff_R_by_c_rad * pt_ev->R_jet_t/vluce_cm;;
            for (TMP = 0; TMP < E_SIZE; TMP++) {
                t_esc=f_Tesc(x[TMP], pt_ev->T_esc_Coeff_rad, pt_ev->Esc_Index_rad);
                pt_ev->T_esc_rad[TMP] = t_esc;
                if (t>=pt_ev->t_jet_exp){
                    t_ad_esc =(t_ad*t_esc)/(t_esc+t_ad);
                    pt_ev->T_esc_rad[TMP] = t_ad_esc;
                    //N_rad[TMP] = N_rad[TMP]*exp_factor;
                }else{
                   
                    pt_ev->T_esc_rad[TMP] = f_Tesc(x[TMP], pt_ev->T_esc_Coeff_rad, pt_ev->Esc_Index_rad);
                }
                
            }
        
        }else{
            pt_ev->T_esc_Coeff_rad = pt_ev->T_esc_Coeff_R_by_c_rad * pt_ev->R_rad_start/vluce_cm;
            for (TMP = 0; TMP < E_SIZE; TMP++) {
                pt_ev->T_esc_rad[TMP] = f_Tesc(x[TMP], pt_ev->T_esc_Coeff_rad, pt_ev->Esc_Index_rad);
            }
        }
        time_evolve_emitters(pt_spec_rad,pt_ev,2,t,T,E_SIZE,E_N_SIZE,E_acc,pt_ev->T_esc_rad,N_escaped,N_rad,N_swap,A,B,C,R,x,xm_p,xm_m,dxm_p,dxm_m,dxm);

        //------------- OUT FILE and SED Computations ----------------
        OUT_FILE=(double)T-COUNT_FILE;
        if ((OUT_FILE >= 0) || (T==0) || (T==pt_ev->T_SIZE-1)) {
        //if ((OUT_FILE >= 0) || (T==pt_ev->T_SIZE-1)) {
            
            //printf("-> NUM_OUT=%d T_SIZE=%d T=%d\n",NUM_OUT,pt_ev->T_SIZE,T);
            //printf(" time=%e, exp_factor=%e, T_esc_Coeff_rad=%e  T_ad=%e  T_ad_esc=%e\n",t,exp_factor,pt_ev->T_esc_Coeff_rad,t_ad,t_ad_esc);
            for (Gamma = 0; Gamma < E_SIZE; Gamma++) {
                
                if (NUM_OUT<pt_ev->NUM_SET){
                    TMP = Gamma + (E_SIZE * NUM_OUT);
                    //pt_ev->N_gamma[TMP] = N_acc[Gamma]*(Vol_acc/Vol_rad)+N_rad[Gamma];
                    pt_ev->N_rad_gamma[TMP] = N_rad[Gamma];
                    pt_ev->N_acc_gamma[TMP] = N_acc[Gamma];
                    
                    
                }
                
                
            }
            if (NUM_OUT<pt_ev->NUM_SET){
                pt_ev->N_time[NUM_OUT]=t;
                //printf("NUM_SET=%d NUM_OUT=%d t=%e T=%d T_SIZE=%d\n",pt_ev->NUM_SET,NUM_OUT,t,T,pt_ev->T_SIZE);
            }
            //if (T>0){
            if (pt_ev->LOG_SET >=1){
                COUNT_FILE*=STEP_FILE;
            } else {
                COUNT_FILE += STEP_FILE;
            }
            //}
            NUM_OUT++;
        }
        //---------------------------------------
        t = t + pt_ev->deltat;

        //if the last element OUT has been skipped, we set time negative to remove from the python wrapper
        if (NUM_OUT == pt_ev->NUM_SET-1){
           pt_ev->N_time[pt_ev->NUM_SET-1]=-1.0;
           //printf("NUM_SET=%d NUM_OUT=%d t=%e T=%d T_SIZE=%d\n",pt_ev->NUM_SET,NUM_OUT,t,T,pt_ev->T_SIZE); 
        }
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
    free(N_swap);
    free(N_escaped);
    free(N_acc);
    free(N_rad);
    //free(T_inj);
    //printf("freeing local dynamic arrays stop \n");
    return;
}

double time_blob_to_obs(double time_blob, struct blob *pt_spec){
    return time_blob*(1+pt_spec->z_cosm)*pt_spec->BulkFactor;
}

double time_obs_to_blob(double time_obs, struct blob *pt_spec){
    return time_obs/((1+pt_spec->z_cosm)*pt_spec->BulkFactor);
}

//double time_blob_to_RH(struct temp_ev *pt, struct blob *pt_spec,double R_H_jet_t){
 //   //t_obs=(pt->R_H_jet_t-pt->R_H_jet)/pt_spec->beta_Gamma*vluce_cm
//    double t_obs;
//    t_obs=(R_H_jet_t-pt->R_H_jet_exp)/(pt_spec->beta_Gamma*vluce_cm);
//    return time_obs_to_blob(t_obs,pt_spec);
//    //return (pt->R_H_jet_t-pt->R_H_jet)/(pt_spec->beta_Gamma*vluce_cm*(1+pt_spec->z_cosm)*pt_spec->BulkFactor);
//}

double eval_R_H_jet_t(struct blob *pt_spec, struct temp_ev *pt_ev, double time_blob){
    //position R_H measured in obs frame at time_blob
    //printf("t=%e, t_obs=%e, R_H_rad_start=%e, time_blob_to_obs=%e, betac=%e\n",time_blob,time_blob_to_obs(time_blob, pt_spec), pt_ev->R_H_rad_start,time_blob_to_obs(time_blob, pt_spec),pt_spec->beta_Gamma);
    return pt_ev->R_H_rad_start + (pt_spec->beta_Gamma*vluce_cm*time_blob_to_obs(time_blob, pt_spec));
}


double eval_R_jet_t(struct blob *pt_spec, struct temp_ev *pt_ev, double time_blob){
    //double R_H_t;
    //R_H_t=eval_R_H_jet_t(pt_spec,pt_ev,time_blob);
    if (time_blob<pt_ev->t_jet_exp){
        return pt_ev->R_rad_start;
    }else{
        return pt_ev->R_rad_start+(pt_ev->v_exp_by_c*vluce_cm*(time_blob - pt_ev->t_jet_exp));
    }
}

double eval_B_jet_t(struct blob *pt_spec, struct temp_ev *pt_ev,double R_jet_t, double time_blob){
    //printf(" R_H_rad_start=%e, R_H_jet_t=%e, B=%e\n", pt_ev->R_H_rad_start,R_H_jet_t,pt_ev->B_rad*pow(pt_ev->R_H_rad_start/R_H_jet_t, pt_ev->m_B));
    //double H0,H_t;
    if (time_blob<pt_ev->t_jet_exp){
        return pt_ev->B_rad;
    }else{
        //H0=pt_ev->R_rad_start/pt_ev->theta_exp_AR;
        //H_t=R_H_jet_t-pt_ev->R_H_jet_exp+H0;
        return pt_ev->B_rad*pow(pt_ev->R_rad_start/R_jet_t, pt_ev->m_B);
    }
}      

double update_jet_expansion(struct blob *pt_spec, struct temp_ev *pt_ev, double t){
    double R_jet_old,exp_factor;
    R_jet_old=pt_ev->R_jet_t;
    pt_spec->beta_Gamma=eval_beta_gamma(pt_spec->BulkFactor);

    pt_ev->R_H_jet_t=eval_R_H_jet_t( pt_spec,  pt_ev, t);
    pt_ev->R_jet_t=eval_R_jet_t(pt_spec,  pt_ev,t);
    pt_ev->B_t=eval_B_jet_t(pt_spec,  pt_ev,pt_ev->R_jet_t,t);
    
    pt_spec->B=pt_ev->B_t;
    pt_spec->R=pt_ev->R_jet_t;
    pt_spec->R_H=pt_ev->R_H_jet_t;
    InitRadiative(pt_spec);
    //if (pt_ev->R_H_jet_t>pt_ev->R_H_jet_exp){ 
    //    T_adiab=3*Adiabatic_Cooling_time(pt_ev,pt_spec,pt_ev->R_jet_t);
    //}else{
    //    T_adiab=-1;
    //}
    exp_factor= eval_N_expansiont_factor(R_jet_old,pt_ev->R_jet_t);
    return exp_factor;
 
}

double eval_N_expansiont_factor(double R_jet_old, double R_jet_new){
    return pow(R_jet_old/R_jet_new,3);

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
    double Uph;
    Uph=0;
    a=3.0*MEC2/(4.0*vluce_cm*(pt->UB + Uph)*SIGTH);
    pt->gamma_cooling_eq=a/T_esc;
    
    for (ID = 0; ID < pt->gamma_grid_size ; ID++){
        
        pt->Ne[ID]=IntegrateCooolingEquilibrium(pt,
                                                pt->griglia_gamma_Ne_log[ID], 
                                                T_esc);
    }


}