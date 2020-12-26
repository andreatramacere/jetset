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

void Init_temp_evolution(struct blob *pt_spec, struct temp_ev *pt_ev, double luminosity_distance)
{
    int grid_bounded_to_gamma;

    unsigned int  gamma_grid_size;

    unsigned int i, TMP;

    double delta_t_eq_t_D, delta_t_eq_t_A, delta_t_eq_t_D_old, delta_t_eq_t_A_old;

    double delta_log;
    double log_a, log_b;

    //FILE *fp;
    //char stringa[static_file_name_max_legth];

   

    build_Ne(pt_spec);

    pt_ev->t_unit = (pt_spec->R / vluce_cm);

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
    pt_ev->T_esc_Coeff = pt_ev->T_esc_Coeff_R_by_c * pt_ev->t_unit;


    pt_ev->gamma_eq_t_D = -1.0;
    pt_ev->gamma_eq_t_A = -1.0;
    
    delta_t_eq_t_A_old = 1.0;
    delta_t_eq_t_D_old = 1.0;

    log_a = log10(pt_ev->gmin_griglia);
    log_b = log10(pt_ev->gmax_griglia);
    delta_log = (log_b - log_a) / ((double)static_ev_arr_grid_size - 1);

    for (i = 0; i <= static_ev_arr_grid_size - 1; i++)
    {
        pt_ev->g[i] = pow(10, (log_a + delta_log * (double)(i)));
        pt_ev->t_Sync_cool[i] = Sync_tcool(pt_spec, pt_ev->g[i]);
        pt_ev->t_D[i] = (pt_ev->g[i] * pt_ev->g[i]) / f_Dp(pt_ev->g[i], pt_ev);
        pt_ev->t_DA[i] = 0.5 * pt_ev->t_D[i];
        pt_ev->t_A[i] = pt_ev->g[i] / f_Acc(pt_ev->g[i], pt_ev);
        pt_ev->t_Esc[i] = f_Tesc(pt_ev->g[i], pt_ev);
        delta_t_eq_t_D = pt_ev->t_Sync_cool[i] - pt_ev->t_D[i];
        delta_t_eq_t_A = pt_ev->t_Sync_cool[i] - pt_ev->t_A[i];

        if (i > 0 && pt_ev->g[i] < pt_ev->Gamma_Max_Turb_L_max)
        {
            if (delta_t_eq_t_D <= 0 && delta_t_eq_t_D_old >= 0)
            {
                pt_ev->gamma_eq_t_D = pt_ev->g[i];
            }

            if (delta_t_eq_t_A <= 0 && delta_t_eq_t_A_old >= 0)
            {
                pt_ev->gamma_eq_t_A = pt_ev->g[i];
            }
            delta_t_eq_t_A_old = delta_t_eq_t_A;
            delta_t_eq_t_D_old = delta_t_eq_t_D;

            //fprintf(fp, "%e %e %e %e %e\n", log10(g), log10(t_Sync_cool[]), log10(t_D), log10(t_DA), log10(t_A));
        }
    }
    //fclose(fp);

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

    

    //-----Injection and Escape Terms-------
    Init(pt_spec, luminosity_distance);
    EvalU_e(pt_spec);
    pt_ev->Q_scalig_factor = pt_ev->L_inj * pt_ev->deltat / pt_spec->E_tot_e;
    
    

    grid_bounded_to_gamma = pt_spec->grid_bounded_to_gamma;
    pt_spec->grid_bounded_to_gamma = 0;

    pt_spec->gmin_griglia = pt_ev->gmin_griglia;
    pt_spec->gmax_griglia = pt_ev->gmax_griglia;
    gamma_grid_size=pt_spec->gamma_grid_size;

    pt_spec->gamma_grid_size = pt_ev->gamma_grid_size;

    Init(pt_spec, luminosity_distance);

    pt_ev->Q_inj = (double *)calloc(pt_ev->gamma_grid_size, sizeof(double));
    pt_ev->gamma = (double *)calloc(pt_ev->gamma_grid_size, sizeof(double));
    
    for (i = 0; i < pt_ev->gamma_grid_size; i++)
    {
        pt_ev->gamma[i] = pt_spec->griglia_gamma_Ne_log[i];
    }

    for (TMP = 0; TMP < pt_ev->gamma_grid_size; TMP++)
    {
        //The injection term is the N_electron
        //passed from the SED computation at first step
        pt_ev->Q_inj[TMP] = pt_spec->Ne[TMP] * pt_ev->Q_scalig_factor;
    }

    pt_spec->grid_bounded_to_gamma=grid_bounded_to_gamma;
    pt_spec->gamma_grid_size = gamma_grid_size;
    
}

void Run_temp_evolution(struct blob *pt_spec, struct temp_ev *pt_ev) {

        // if luminosity_distance is negative is evaluated internally
        // otherwise the passed value is used

        unsigned int i, E_SIZE, E_N_SIZE, Gamma, T, TMP,NUM_OUT;
        double Q_scalig_factor;

        double STEP_FILE, COUNT_FILE, OUT_FILE;
        double *x, *N1;
        double  t;
        //double g, t_D, t_DA, t_A, t_Sync_cool;
        double *xm_p, *xm_m;
        double *dxm_p, *dxm_m, *dxm;
        double K, Cp, Cm, wm_p, wm_m;
        double WP_pm, WM_mm, WP_mm, WM_pm;
        double *A, *B, *C, *R, *T_inj;
        int grid_bounded_to_gamma;
        unsigned int gamma_grid_size;

        
    

        //----SCELTA DELLA VARIABILE ADIMENSIONALE-------------------------
        //    E=mec2*gamma, mec2 rest electron energy
        //    beta=v/c
        //    E_kin=E-mec2=mec2(gamma-1)
        //    x=E_kin/mec2=(g-1)=beta*gamma
        //------------------------------------------------------------------



    if (pt_spec->do_Sync > 0)
        {
            pt_ev->do_Sync_cooling = 1;
        }

    if (pt_ev->do_Compton_cooling>0){
        if (pt_spec->do_SSC>0){
            pt_ev->do_SSC_cooling=1;
        }

        if (pt_spec->do_EC_Disk>0){
            pt_ev->do_EC_cooling_Disk=1;
        }

        if (pt_spec->do_EC_BLR>0){
            pt_ev->do_EC_cooling_BLR=1;
        }

        if (pt_spec->do_EC_DT>0){
            pt_ev->do_EC_cooling_DT=1;
        }

        if (pt_spec->do_EC_Star>0){
            pt_ev->do_EC_cooling_Star=1;
        }

        if (pt_spec->do_EC_CMB>0){
            pt_ev->do_EC_cooling_CMB=1;
        }
    }

    //For temp_ev we need to have boundaris of gamma grid
    //unbound from gmin and gmax
    //and to change the grid buoundaries and size
    grid_bounded_to_gamma = pt_spec->grid_bounded_to_gamma;
    pt_spec->grid_bounded_to_gamma = 0;

    pt_spec->gmin_griglia = pt_ev->gmin_griglia;
    pt_spec->gmax_griglia = pt_ev->gmax_griglia;
    gamma_grid_size=pt_spec->gamma_grid_size;

    pt_spec->gamma_grid_size = pt_ev->gamma_grid_size;

    //--------------ALLOCATION OF DYNAMICAL ARRAYS-----------------------
    //--------------SIZE DEFINITION------------------
    //T_SIZE=10000;
    //SIZE=1000;
    //---Useful size definition for 0..N-1----------
    E_SIZE = pt_ev->gamma_grid_size;
    E_N_SIZE = pt_ev->gamma_grid_size - 1;
    //------ALLOCATION-------------------------------
    pt_ev->N_gamma = (double *)calloc(pt_ev->NUM_SET * E_SIZE, sizeof(double));
    pt_ev->N_time = (double *)calloc(pt_ev->T_SIZE, sizeof(double));
    x = (double *)calloc(E_SIZE, sizeof(double));
    //--------Gamma energia normalizzata-------------
    N1 = (double *) calloc(E_SIZE, sizeof (double));
    
    A = (double *) calloc(E_SIZE, sizeof (double));
    B = (double *) calloc(E_SIZE, sizeof (double));
    C = (double *) calloc(E_SIZE, sizeof (double));
    R = (double *) calloc(E_SIZE, sizeof (double));
    
    T_inj = (double *) calloc(pt_ev->T_SIZE, sizeof (double));

    pt_ev->T_esc = (double *) calloc(E_SIZE, sizeof (double));

    xm_p = (double *) calloc(E_SIZE, sizeof (double));
    xm_m = (double *) calloc(E_SIZE, sizeof (double));
    dxm_p = (double *) calloc(E_SIZE, sizeof (double));
    dxm_m = (double *) calloc(E_SIZE, sizeof (double));
    dxm = (double *) calloc(E_SIZE, sizeof (double));
    //--------------------------------------------------------------------



    //----------------Griglia X--------------------------------
    //x=gamma-1
    //griglia(x,1e-2,5e8,SIZE);
    for (i = 0; i <= E_SIZE - 1; i++) {
        x[i] = pt_spec->griglia_gamma_Ne_log[i] - 1;
        //printf("Gamma[%d]=%e\n",i,x[i]);
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
    //printf("x grid size=%d\n", E_SIZE);
    //printf("x (x=gamma-1) grid boundaries min=%e max=%e\n", x[0], x[E_SIZE - 1]);
    //---------------------------------------------------------
    
    //------SET TO WRITE ON OUTPUT -------------
    if (pt_ev->LOG_SET >=1){
        STEP_FILE=log10((double) pt_ev->T_SIZE)/(double) pt_ev->NUM_SET;
        //printf("STEP_FILE=%f\n", STEP_FILE);
        STEP_FILE =  pow(10,STEP_FILE);
        //("STEP_FILE=%f\n", STEP_FILE);
    } else{
        STEP_FILE = (double)pt_ev->T_SIZE / (double)pt_ev->NUM_SET;
    }
    COUNT_FILE = STEP_FILE;
    OUT_FILE=-1.0;
    //------------------------------------------


    

    t = 0;
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        t = t + pt_ev->deltat;
        //printf("T=%d\n",T);
        T_inj[T] = Inj_temp_prof(t, pt_ev);
        //fprintf(fp, "%e %e \n",t,T_inj[T]);
    }
    //fclose(fp);
    


    for (TMP = 0; TMP < E_SIZE; TMP++) {

        pt_ev->T_esc[TMP] = f_Tesc(x[TMP], pt_ev);
    }

  
    //--------------Loop over Gamma------------------------------------------
    t = 0;
    NUM_OUT=0;
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        //printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>  temporal evolution step=%d,%e,%e <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n", T, OUT_FILE,COUNT_FILE);
        if (pt_ev->do_Compton_cooling == 1) {
            printf(">>>>>>>>>>>>>>>> Eval Sync \n");
            //Run_SED(pt_spec, pt_ev);
            spettro_sincrotrone(1, pt_spec);
        }
        t = t + pt_ev->deltat;
        //printf("T=%d t=%e\n",T,t);

        //-------G[0]---------------------------

        //--ACC+COOLIG--------------------------
        K = pt_ev->deltat / (dxm_p[0]);
        if (t >= pt_ev->TStart_Acc && t <= pt_ev->TStop_Acc) {
            //printf("Acc+Cool\n");
            Cp = Cfp(xm_p[0], pt_ev) / dxm_p[0];
            wm_p = Bfp(xm_p[0], pt_ev, pt_spec) / (Cfp(xm_p[0], pt_ev)) * dxm_p[0];

            //W^{+}_{m+1/2} W^{-}_{m+1/2}
            Wm(wm_p, &WP_pm, &WM_pm);


            //W^{+}_{m+1/2}
            //WP_pm=Wm_p*exp(+wm_p*0.5);

            //W^{-}_{m+1/2}
            //WM_pm=Wm_p*exp(-wm_p*0.5);

            //-------A[0]---------
            A[0] = 0;

            //-------C[0]---------
            C[0] = -1 * K * Cp*WP_pm;

            //-------B[0]----------
            //B[0]=1+deltat/(pt_ev->T_esc[0]);
            B[0] = 1 + K * Cp * WM_pm + pt_ev->deltat / (pt_ev->T_esc[0]);

            //-------R[0]----------
            R[0] = pt_ev->deltat * pt_ev->Q_inj[0] * T_inj[T] + pt_spec->Ne[0];
            //printf("K=%e wm_p=%e Wm_p=%e WP_pm=%e\n",K,wm_p,Wm_p,WP_pm);
            //printf("A=%e B=%e C=%e R=%e\n",A[0],B[0],C[0],R[0]);
            //---------------------------------------
        }            //----ONLY COOLING---------------------
        else {
            //printf("Cool\n");
            A[0] = 0;
            B[0] = 1 + K * Cooling(xm_m[0], pt_ev, pt_spec) + pt_ev->deltat / (pt_ev->T_esc[0]);
            C[0] = -1 * K * Cooling(xm_p[0], pt_ev, pt_spec);
            R[0] = pt_ev->deltat * pt_ev->Q_inj[0] * T_inj[T] + pt_spec->Ne[0];
        }
        for (Gamma = 1; Gamma < E_N_SIZE; Gamma++) {
            //--------COOLING + ACC--------------
            K = pt_ev->deltat / dxm[Gamma];
            if (t >= pt_ev->TStart_Acc && t <= pt_ev->TStop_Acc) {
                //C_{m-1/2}
                Cm = Cfp(xm_m[Gamma], pt_ev) / dxm_m[Gamma];

                //C_{m+1/2]
                Cp = Cfp(xm_p[Gamma], pt_ev) / dxm_p[Gamma];

                //w_{m+1/2}
                wm_p = Bfp(xm_p[Gamma], pt_ev, pt_spec) / (Cfp(xm_p[Gamma], pt_ev)) * dxm_p[Gamma];

                //w_{m-1/2}
                wm_m = Bfp(xm_m[Gamma], pt_ev, pt_spec) / (Cfp(xm_m[Gamma], pt_ev)) * dxm_m[Gamma];

                //printf("gm_p=%e, gm_m=%e\n",gm_p[G],gm_m[G]);
                //printf("wm_p=%e, wm_m=%e\n",wm_p,wm_m);

                //      W^{+}_{m+1/2} W^{-}^{m+1/2}
                Wm(wm_p, &WP_pm, &WM_pm);

                //      W^{+}_{m-1/2} W^{-}^{m-1/2}
                Wm(wm_m, &WP_mm, &WM_mm);

                //W^{+}_{m+1/2}
                //WP_pm=Wm_p*exp(wm_p*0.5);

                //W^{+}_{m-1/2}
                //WP_mm=Wm_m*exp(wm_m*0.5);

                //W^{-}_{m+1/2}
                //WM_pm=Wm_p*exp(-wm_p*0.5);

                //W^{-}_{m-1/2}
                //WM_mm=Wm_m*exp(-wm_m*0.5);


                A[Gamma] = -1 * (K * Cm * WM_mm);
                C[Gamma] = -1 * (K * Cp * WP_pm);
                B[Gamma] = 1 + K * (Cm * WP_mm + Cp * WM_pm) + pt_ev->deltat / (pt_ev->T_esc[Gamma]);
                R[Gamma] = pt_ev->deltat * pt_ev->Q_inj[Gamma] * T_inj[T] + pt_spec->Ne[Gamma];
            }//-----ONLY COOLING-----------------
            else {
                A[Gamma] = 0;
                B[Gamma] = 1 + pt_ev->deltat / (pt_ev->T_esc[Gamma]) + K * Cooling(xm_m[Gamma], pt_ev, pt_spec);
                C[Gamma] = -1 * K * Cooling(xm_p[Gamma], pt_ev, pt_spec);
                R[Gamma] = pt_ev->deltat * pt_ev->Q_inj[Gamma] * T_inj[T] + pt_spec->Ne[Gamma];
            }

            //printf("A=%e B=%e C=%e R=%e\n",-A[G],B[G],-C[G],R[G]);
        }

        //--------G[N_SIZE]----------------------
        //--------COOLING + ACC--------------
        K = pt_ev->deltat / (dxm_m[E_N_SIZE]);
        if (t >= pt_ev->TStart_Acc && t <= pt_ev->TStop_Acc) {
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
            B[E_N_SIZE] = 1 + K * Cm * WP_mm + pt_ev->deltat / (pt_ev->T_esc[E_N_SIZE]);
            //--------R[N_SIZE]------------
            R[E_N_SIZE] = pt_ev->deltat * pt_ev->Q_inj[E_N_SIZE] * T_inj[T] + pt_spec->Ne[E_N_SIZE];
            //printf("A=%e B=%e C=%e R=%e\n",A[N_SIZE],B[N_SIZE],C[N_SIZE],R[N_SIZE]);
        } else {
            A[E_N_SIZE] = 0;
            B[E_N_SIZE] = 1 + pt_ev->deltat / (pt_ev->T_esc[E_N_SIZE]) + K * Cooling(xm_m[E_N_SIZE], pt_ev, pt_spec);
            C[E_N_SIZE] = 0;
            R[E_N_SIZE] = pt_ev->deltat * pt_ev->Q_inj[E_N_SIZE] * T_inj[T] + pt_spec->Ne[E_N_SIZE];
        }


        //----------------------------------------

        //--------------END Loop over Gamma---------------------------------------


        //--------------TRIDIAG SYSTEM SOLVING--------------------
        //printf("********prima di sys****************\n");
        for (TMP = 0; TMP < E_SIZE; TMP++) {
            N1[TMP] = pt_spec->Ne[TMP];
            //printf("N=%e T=%d G=%d\n",N1[TMP],T,TMP);
        }
        if (solve_sys1(A, B, C, R, N1, E_SIZE) > 0) {
            printf("errore nella soluzione del sistema condizione di positivita' non soddisfatta\n ");
            exit(1);
        }
        //printf("********dopo sys********************\n");
        for (TMP = 0; TMP < E_SIZE; TMP++) {
            pt_spec->Ne[TMP] = N1[TMP];
            //printf("N=%e T=%d G=%d\n",N1[TMP],T,TMP);
        }
        //--------------------------------------------------------


        //------------- OUT FILE and SED Computations ----------------
        OUT_FILE=(double)T-COUNT_FILE;


        if ((OUT_FILE >= 0) || (T==pt_ev->T_SIZE-1)) {
            //printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>  temporal evolution step=%d, t_sim frac=%e, t_sim=%e <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n", T, t/pt_ev->duration,t);
            //strcpy(old, pt_ev->STEM);
            //sprintf(name, "%s-TSTEP=%6.6d", pt_ev->STEM, T);
            //strcpy(pt_ev->STEM, name);
            //strcpy(pt_spec->STEM, name);
            //sprintf(name1, "distr-e-evol.dat");

            //sprintf(name_n_evol, "%s%s-distr-e-evol.dat", pt_ev->path, pt_ev->STEM);

            //fp1 = fopen(name_n_evol, "w");
            //if (fp1 == NULL)
            //{
            //    printf("warning non riesco ad aprire %s\n", fp);
            //    exit(1);
            //}

            //EvalU_e(pt_spec);
            //Run_SED(pt_spec);

            //Scrivi_N_file(pt_spec, name1, pt_spec->griglia_gamma_Ne_log, pt_spec->Ne);
            //strcpy(pt_ev->STEM, old);

            printf("-> NUM_OUT=%d T_SIZE=%d T=%d\n",NUM_OUT,pt_ev->T_SIZE,T);

            for (Gamma = 0; Gamma < E_SIZE; Gamma++) {
                
                if (NUM_OUT<pt_ev->NUM_SET){
                    TMP = Gamma + (E_SIZE * NUM_OUT);
                    pt_ev->N_gamma[TMP] = N1[Gamma];
                    pt_ev->N_time[NUM_OUT]=t;
                    
                }
                
                //if (pt_spec->Ne[Gamma] != 0) {
                    //fprintf(fp, "%e %e %e %e", t, log10(x[Gamma] + 1), log10(pt_spec->Ne[Gamma]), log10(pt_spec->Ne[Gamma]*(x[Gamma] + 1)*(x[Gamma] + 1)*(x[Gamma] + 1)));
                    //fprintf(fp1, " %e %e %e\n", log10(x[Gamma] + 1), log10(pt_spec->Ne[Gamma]), log10(pt_spec->Ne[Gamma] * (x[Gamma] + 1) * (x[Gamma] + 1) * (x[Gamma] + 1)));
                    //fprintf(fp,"%e %e",x[G]+1,pt->Ne[G]);
                    //fprintf(fp,"\n");
                    //fprintf(fp, "\n");
                //}
            }
            //fprintf(fp, "&&\n");
            ////fprintf(fp1,"&&\n");
            //if (T > 1) {
            //    fprintf(fp, "&&\n");
            //    //fprintf(fp1,"&&\n");
            //}
            //fclose(fp1);
            
            if (pt_ev->LOG_SET >=1){
                COUNT_FILE*=STEP_FILE;
            } else {
                COUNT_FILE += STEP_FILE;
            }
            NUM_OUT++;
        }
        //---------------------------------------
        //printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>> ------------------------  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
    }
    //--------------END Loop over Time-----------------------------------------------


    //fclose(fp);

    pt_spec->grid_bounded_to_gamma=grid_bounded_to_gamma;
    pt_spec->gamma_grid_size = gamma_grid_size;
    
    free(xm_p);
    free(xm_m);
    free(dxm_p);
    free(dxm);
    free(A);
    free(B);
    free(C);
    free(R);
    

    free(T_inj);
    return;
}

void free_tempe_ev(struct temp_ev *pt_ev)
{
    free(pt_ev->Q_inj);
    free(pt_ev->gamma);
    free(pt_ev->N_gamma);
    free(pt_ev->N_time);
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
    
    double  delta_t_cool, delta_t_cool_old, t_cool_sync;
    unsigned int ID;
    
    pt->gamma_cooling_eq=(10/Sync_cool(pt,10))/T_esc;
    
    for (ID = 0; ID < pt->gamma_grid_size ; ID++){
        
        pt->Ne[ID]=IntegrateCooolingEquilibrium(pt,
                                                pt->griglia_gamma_Ne_log[ID], 
                                                T_esc);
    }


}