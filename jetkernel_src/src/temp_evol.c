//=================================================================
//
//  FUNZIONI PER P-P Process
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

void temp_evolution(struct spettro *pt_spec, struct temp_ev *pt_ev) {
    unsigned long i, E_SIZE, E_N_SIZE, Gamma, T, TMP;
    double Q_scalig_factor;

    double STEP_FILE,COUNT_FILE,OUT_FILE;
    double *x, *N1;
    double duration, deltat, t;
    double g,t_D,t_DA,t_A,t_Sync_cool;
    double *xm_p, *xm_m;
    double *dxm_p, *dxm_m, *dxm;
    double K, Cp, Cm, wm_p, wm_m;
    double WP_pm, WM_mm, WP_mm, WM_pm;
    double *A, *B, *C, *R, *Q, *T_inj;
    //double t_cool_Sync, B2;
    //double t_acc_sist, t_acc_stoc;
    //double c2;
    char stringa[512];
    //double dato;
    //int dato_i;
    char name[512], old[512], name1[512];
    FILE *fp;





    //----SCELTA DELLA VARIABILE ADIMENSIONALE-------------------------
    //    E=mec2*gamma, mec2 rest electron energy
    //    beta=v/c
    //    E_kin=E-mec2=mec2(gamma-1)
    //    x=E_kin/mec2=(g-1)=beta*gamma
    //------------------------------------------------------------------


    /*
      fscanf(fp,"%s %d",stringa,&dato_i);
      SIZE=dato_i;
      printf("E_SIZE=%d\n",SIZE);
     */

    /*
      fscanf(fp,"%s %lf",stringa,&dato);
      pt_ev->gmin_inj=dato;
      printf("gmin_inj=%e\n",pt_ev->gmin_inj);


      fscanf(fp,"%s %lf",stringa,&dato);
      pt_ev->gmax_inj=dato;
      printf("gmax_inj=%e\n",pt_ev->gmax_inj);

      fscanf(fp,"%s %lf",stringa,&dato);
      pt_ev->alfa_inj=dato;
      printf("alfa_inj=%e\n",pt_ev->alfa_inj);

      fscanf(fp,"%s %lf",stringa,&dato);
      pt_ev->R=dato;
      printf("R=%e\n",pt_ev->R);

      fscanf(fp,"%s %lf",stringa,&dato);
      pt_ev->B=dato;
      printf("B=%e\n",pt_ev->B);
     */



    if (pt_spec->do_Sync>0){
    	pt_ev->do_Sync_cooling=1;
    }

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

    //--------------ALLOCATION OF DYNAMICAL ARRAYS-----------------------
    //--------------SIZE DEFINITION------------------
    //T_SIZE=10000;
    //SIZE=1000;
    //---Useful size definition for 0..N-1----------
    E_SIZE = pt_spec->gamma_grid_size;
    E_N_SIZE = pt_spec->gamma_grid_size - 1;
    //------ALLOCATION-------------------------------
    x = (double *) calloc(E_SIZE, sizeof (double));
    //--------Gamma energia normalizzata-------------
    N1 = (double *) calloc(E_SIZE, sizeof (double));
    //Pt->N
    //ESIZE=pt->gamma_grid_size
    //N =(double *)calloc(E_SIZE,sizeof(double));
    A = (double *) calloc(E_SIZE, sizeof (double));
    B = (double *) calloc(E_SIZE, sizeof (double));
    C = (double *) calloc(E_SIZE, sizeof (double));
    R = (double *) calloc(E_SIZE, sizeof (double));
    Q = (double *) calloc(E_SIZE, sizeof (double));
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
    //---------------------------------------------------------


    //--------Coeff and Phys Inizializations & Times Iinizialization--------

    //---Size and Mag. filed ----------
    //(1/8PI) assorbito in
    //Sync cooling cost vedi stoc_acc.h
    //B2=pt_ev->B*pt_ev->B;
    //printf("Sync cooling cost=%e\n",);

    //---------------------------------


    printf("\n\n");


    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>  temporal evolution: Parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    //---Duration Times definition ----
    pt_ev->t_unit = (pt_spec->R / vluce_cm);

    duration =  pt_ev->duration;

    //pt_ev->TStart_Inj = pt_ev->t_unit *pt_ev->TStart_Inj;
    //pt_ev->TStop_Inj = pt_ev->t_unit * pt_ev->TStop_Inj;


    //pt_ev->TStart_Acc = pt_ev->t_unit * pt_ev->TStart_Acc;
    //pt_ev->TStop_Acc = pt_ev->t_unit * pt_ev->TStop_Acc;

    deltat = duration / (double) pt_ev->T_SIZE;

    //----- Diffusion coeff ----
    pt_ev->t_DA0=pt_ev->t_D0*0.5;
    pt_ev->Diff_Coeff = 1.0 / (pt_ev->t_D0);
    printf("diff coeff %e\n", pt_ev->Diff_Coeff);
    printf("diff index %e\n", pt_ev->Diff_Index);

    //----- Acc coeff --------
    pt_ev->Acc_Coeff = 1.0 / (pt_ev->t_A0);
    printf("acc coeff %e\n", pt_ev->Acc_Coeff);

    printf("Lambda max Turbulence %e\n", pt_ev->Lambda_max_Turb);
    printf("Lambda coherence Turbulence %e\n", pt_ev->Lambda_max_Turb*pt_ev->Lambda_choer_Turb_factor);
    if (pt_ev->Lambda_max_Turb>pt_spec->R){
    	printf("Warning Lambda_max_Turb=%e cm> R =%e cm",pt_ev->Lambda_max_Turb,pt_spec->R);
    }

    pt_ev->Gamma_Max_Turb_L_max= Larmor_radius_to_gamma(pt_ev->Lambda_max_Turb,pt_spec->B, pt_spec->sin_psi);
    pt_ev->Gamma_Max_Turb_L_coher= Larmor_radius_to_gamma(pt_ev->Lambda_max_Turb*pt_ev->Lambda_choer_Turb_factor,pt_spec->B, pt_spec->sin_psi);
    pt_ev->Diff_coeff_CD=pt_ev->Diff_Coeff*pow((pt_ev->Gamma_Max_Turb_L_coher),pt_ev->Diff_Index);
    pt_ev->Diff_coeff_CA=pt_ev->Diff_Coeff*2*pow((pt_ev->Gamma_Max_Turb_L_coher),pt_ev->Diff_Index);

    printf("gamma max Turbulence Lambda max %e\n", pt_ev->Gamma_Max_Turb_L_max);
    printf("gamma max Turbulence Coherence length %e\n", pt_ev->Gamma_Max_Turb_L_coher);

    //---- x grids -------------
	printf("x grid size=%d\n",E_SIZE);
	printf("x (x=gamma-1) grid boundaries min=%e max=%e\n",x[0],x[E_SIZE-1]);

    //----- Esc coeff --------
    pt_ev->T_esc_Coeff *= pt_ev->t_unit;
    printf("T_esc Coeff %e (s) %e(R/c)\n", pt_ev->T_esc_Coeff,pt_ev->T_esc_Coeff/ pt_ev->t_unit);




    printf("T_SIZE=%d\n", pt_ev->T_SIZE);
    printf("(R/c)=%e(s) deltat=%e(s) \n", pt_ev->t_unit, deltat);
    printf("durata=%e(s), %e(R/c)\n", duration, duration / pt_ev->t_unit);
    //------SET TO WRITE ON OUTPUT -------------
    STEP_FILE=log10((double) pt_ev->T_SIZE)/(double) pt_ev->NUM_SET;
    printf("STEP_FILE=%f\n", STEP_FILE);
    STEP_FILE =  pow(10,STEP_FILE);
    printf("STEP_FILE=%f\n", STEP_FILE);
    COUNT_FILE = STEP_FILE;
    OUT_FILE=-1.0;
    //------------------------------------------


    printf("Tstart_Inj=%e(s) %e(R/c)\n", pt_ev->TStart_Inj, pt_ev->TStart_Inj / pt_ev->t_unit);
    printf("Tstop_Inj=%e(s) %e(R/c)\n", pt_ev->TStop_Inj, pt_ev->TStop_Inj / pt_ev->t_unit);
    printf("Tstart_Acc=%e(s) %e(R/c)\n", pt_ev->TStart_Acc, pt_ev->TStart_Acc / pt_ev->t_unit);
    printf("Tstop_Acc=%e(s) %e(R/c)\n", pt_ev->TStop_Acc, pt_ev->TStop_Acc / pt_ev->t_unit);
    //----------------------------------

    printf("Tesc=%e\n", pt_ev->T_esc_Coeff);
    printf("T_Cooling=  T_cooling(g=1)   =%e(s) %e(R/c)\n", Sync_tcool(pt_spec, 1.0), Sync_tcool(pt_spec, 1.0) / pt_ev->t_unit);
    printf("T_A0  =1/ACC_COEFF      =%e(s) %e(R/c)\n", pt_ev->t_A0, pt_ev->t_A0 / pt_ev->t_unit);
    printf("T_D0  =1/DIFF_COEFF     =%e(s) %e(R/c)\n", pt_ev->t_D0, pt_ev->t_D0/ pt_ev->t_unit);
    printf("T_DA0=1/(2DIFF_COEFF)  =%e(s) %e(R/c)\n", pt_ev->t_DA0, (pt_ev->t_DA0 / pt_ev->t_unit));
    printf("\n");



    printf("\n");

    t = 0;
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        t = t + deltat;
        //printf("T=%d\n",T);
        T_inj[T] = Inj_temp_prof(t, pt_ev);
    }

    sprintf(stringa, "%s%s-t_cool_t_acc.dat", pt_spec->path, pt_spec->STEM);



    fp = fopen(stringa, "w");
    if (fp == NULL) {
    	printf("warning non riesco ad aprire %s\n", fp);
        exit(1);
    }
    fprintf(fp, "#gamma t_Sync_cool t_D t_DA t_A\n");

    double delta_t_eq_t_D,delta_t_eq_t_A,delta_t_eq_t_D_old,delta_t_eq_t_A_old;
    pt_ev->gamma_eq_t_D=-1.0;
    pt_ev->gamma_eq_t_A=-1.0;
    for (i = 0; i <= E_SIZE - 1; i++) {
    	g=pt_spec->griglia_gamma_Ne_log[i];
    	t_Sync_cool=Sync_tcool(pt_spec, g);
    	t_D=(g*g)/f_Dp(g,pt_ev);
    	t_DA = 0.5 * t_D;
    	t_A=g/f_A(g,pt_ev);
    	delta_t_eq_t_D=t_Sync_cool-t_DA;
    	delta_t_eq_t_A=t_Sync_cool-t_A;


    	if (i>0 && g<pt_ev->Gamma_Max_Turb_L_max ){
    		if (delta_t_eq_t_D<=0 && delta_t_eq_t_D_old>=0){
    			pt_ev->gamma_eq_t_D=g;
    		}


			if (delta_t_eq_t_A<=0 && delta_t_eq_t_A>=0){
				pt_ev->gamma_eq_t_A=g;
			}
			delta_t_eq_t_A_old=delta_t_eq_t_A;
			delta_t_eq_t_D_old=delta_t_eq_t_D;

    	fprintf(fp, "%e %e %e %e %e\n",log10(g),log10(t_Sync_cool),log10(t_D),log10(t_DA),log10(t_A));
    	}
    }
    fclose(fp);

    // if gamma_eq still negative, then t_acc never crossed t_cool

    if (delta_t_eq_t_A_old>=0 && pt_ev->gamma_eq_t_A<0){
    	pt_ev->gamma_eq_t_A=pt_ev->Gamma_Max_Turb_L_max;
    }
    else if(delta_t_eq_t_A_old<=0 && pt_ev->gamma_eq_t_A<0){
    	pt_ev->gamma_eq_t_A=pt_spec->griglia_gamma_Ne_log[0];
    }
    if (delta_t_eq_t_D_old>=0 && pt_ev->gamma_eq_t_D<0){
        	pt_ev->gamma_eq_t_D=pt_ev->Gamma_Max_Turb_L_max;
	}
    else if(delta_t_eq_t_D_old<=0 && pt_ev->gamma_eq_t_D<0){
		pt_ev->gamma_eq_t_D=pt_spec->griglia_gamma_Ne_log[0];
	}


    if (pt_ev->gamma_eq_t_A>pt_ev->Gamma_Max_Turb_L_max){
    	pt_ev->gamma_eq_t_A=pt_ev->Gamma_Max_Turb_L_max;
    }
    if (pt_ev->gamma_eq_t_D>pt_ev->Gamma_Max_Turb_L_max){
    	pt_ev->gamma_eq_t_D=pt_ev->Gamma_Max_Turb_L_max;
    }

    printf("gamma_eq_Cool=T_cooling0/T_esc0=%e\n",
    		Sync_tcool(pt_spec, 1.0) / pt_ev->T_esc_Coeff);
    printf("gamma_eq_Sys= %e\n",pt_ev->gamma_eq_t_A);
    printf("gamma_eq_Diff= %e\n",pt_ev->gamma_eq_t_D);


    //-----Injection and Escape Terms-------
    Init(pt_spec);
    EvalU_e(pt_spec);
    Q_scalig_factor = pt_ev->L_inj *deltat / (pt_spec->U_e * pt_spec->Vol_sphere);
    printf("old N_e=%e \n",pt_spec->N);
    pt_spec->N= pt_spec->N*Q_scalig_factor;
    Init(pt_spec);
    printf("new N_e=%e,  multiplied by the Q_scalig_factor=%e, to match L_inj \n",pt_spec->N,Q_scalig_factor);
    for (TMP = 0; TMP < E_SIZE; TMP++) {
        //The injection term is the N_electron
        //passed from the SED computation at first step
        Q[TMP] = pt_spec->Ne[TMP] ;
        pt_ev->T_esc[TMP] = f_Tesc(x[TMP], pt_ev);
    }

    printf("L_inj=%e\n",pt_ev->L_inj);
    pt_spec->N_tot_e_Sferic = pt_spec->Vol_sphere * pt_spec->N*Q_scalig_factor;
    printf("Total number of electrons injected at each step   =%e\n", pt_spec->N_tot_e_Sferic);
    EvalU_e(pt_spec);
    printf("U_e   blob rest frame =%e erg/cm^3\n", pt_spec->U_e);
    printf("U_e/U_b =%e\n", pt_spec->U_e / pt_spec->UB);
    printf("E_tot (electron)  blob rest frame =%e erg \n", pt_spec->E_tot_e);
    printf("E_tot (electron)/delta_T (must be equal to L_inj)  =%e erg/s\n", pt_spec->E_tot_e/deltat);

    FindNe_NpGp(pt_spec);
    printf("Gamma_p of N(gamma)*gamma^2 = %e\n", pt_spec->Gamma_p2);
    printf("Gamma_p of N(gamma)*gamma^3 = %e\n", pt_spec->Gamma_p3);
    printf("Peak of  N(gamma)*gamma^2 = %e\n", pt_spec->Np2);
    printf("Peak of  N(gamma)*gamma^3 = %e\n", pt_spec-> Np3);



    strcpy(old, pt_spec->STEM);
    T=0;
    sprintf(name, "%s-TSTEP=%6.6d", pt_spec->STEM, T);
    strcpy(pt_spec->STEM, name);
    sprintf(name1, "distr-e-inj.dat");
    Scrivi_N_file(pt_spec, name1, pt_spec->griglia_gamma_Ne_log, pt_spec->Ne);
    strcpy(pt_spec->STEM, old);

    //---------------------------------------------------------------------


    printf("\n\n");

    //------------------------Loop over Time-----------------------------------------

    T = 0;



    //--------FILES IN USCITA-----------
    sprintf(stringa, "%s%s-distr-e-evol.dat", pt_spec->path, pt_spec->STEM);



    fp = fopen(stringa, "w");
    if (fp == NULL) {
        printf("warning non riesco ad aprire %s\n", fp);
        exit(1);
    }
    distr_e_header(fp);






    //--------------Loop over Gamma------------------------------------------
    t = 0;
    for (T = 0; T < pt_ev->T_SIZE; T++) {
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>  temporal evolution step=%d,%e,%e <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n", T, OUT_FILE,COUNT_FILE);
        if (pt_spec->attesa_compton_cooling == 1) {
            printf(">>>>>>>>>>>>>>>> Eval Sync \n");
            //Run_SED(pt_spec, pt_ev);
            spettro_sincrotrone(1, pt_spec);
        }
        t = t + deltat;
        //printf("T=%d t=%e\n",T,t);

        //-------G[0]---------------------------

        //--ACC+COOLIG--------------------------
        K = deltat / (dxm_p[0]);
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
            B[0] = 1 + K * Cp * WM_pm + deltat / (pt_ev->T_esc[0]);

            //-------R[0]----------
            R[0] = deltat * Q[0] * T_inj[T] + pt_spec->Ne[0];
            //printf("K=%e wm_p=%e Wm_p=%e WP_pm=%e\n",K,wm_p,Wm_p,WP_pm);
            //printf("A=%e B=%e C=%e R=%e\n",A[0],B[0],C[0],R[0]);
            //---------------------------------------
        }            //----ONLY COOLING---------------------
        else {
            //printf("Cool\n");
            A[0] = 0;
            B[0] = 1 + K * Cooling(xm_m[0], pt_ev, pt_spec) + deltat / (pt_ev->T_esc[0]);
            C[0] = -1 * K * Cooling(xm_p[0], pt_ev, pt_spec);
            R[0] = deltat * Q[0] * T_inj[T] + pt_spec->Ne[0];
        }
        for (Gamma = 1; Gamma < E_N_SIZE; Gamma++) {
            //--------COOLING + ACC--------------
            K = deltat / dxm[Gamma];
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
                B[Gamma] = 1 + K * (Cm * WP_mm + Cp * WM_pm) + deltat / (pt_ev->T_esc[Gamma]);
                R[Gamma] = deltat * Q[Gamma] * T_inj[T] + pt_spec->Ne[Gamma];
            }//-----ONLY COOLING-----------------
            else {
                A[Gamma] = 0;
                B[Gamma] = 1 + deltat / (pt_ev->T_esc[Gamma]) + K * Cooling(xm_m[Gamma], pt_ev, pt_spec);
                C[Gamma] = -1 * K * Cooling(xm_p[Gamma], pt_ev, pt_spec);
                R[Gamma] = deltat * Q[Gamma] * T_inj[T] + pt_spec->Ne[Gamma];
            }

            //printf("A=%e B=%e C=%e R=%e\n",-A[G],B[G],-C[G],R[G]);
        }

        //--------G[N_SIZE]----------------------
        //--------COOLING + ACC--------------
        K = deltat / (dxm_m[E_N_SIZE]);
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
            B[E_N_SIZE] = 1 + K * Cm * WP_mm + deltat / (pt_ev->T_esc[E_N_SIZE]);
            //--------R[N_SIZE]------------
            R[E_N_SIZE] = deltat * Q[E_N_SIZE] * T_inj[T] + pt_spec->Ne[E_N_SIZE];
            //printf("A=%e B=%e C=%e R=%e\n",A[N_SIZE],B[N_SIZE],C[N_SIZE],R[N_SIZE]);
        } else {
            A[E_N_SIZE] = 0;
            B[E_N_SIZE] = 1 + deltat / (pt_ev->T_esc[E_N_SIZE]) + K * Cooling(xm_m[E_N_SIZE], pt_ev, pt_spec);
            C[E_N_SIZE] = 0;
            R[E_N_SIZE] = deltat * Q[E_N_SIZE] * T_inj[T] + pt_spec->Ne[E_N_SIZE];
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
            strcpy(old, pt_spec->STEM);
            sprintf(name, "%s-TSTEP=%6.6d", pt_spec->STEM, T);
            strcpy(pt_spec->STEM, name);
            sprintf(name1, "distr-e-evol.dat");
            EvalU_e(pt_spec);
            Run_SED(pt_spec);
            Scrivi_N_file(pt_spec, name1, pt_spec->griglia_gamma_Ne_log, pt_spec->Ne);
            strcpy(pt_spec->STEM, old);

            for (Gamma = 0; Gamma < E_SIZE; Gamma++) {
                if (pt_spec->Ne[Gamma] != 0) {
                    fprintf(fp, "%e %e %e %e", t, log10(x[Gamma] + 1), log10(pt_spec->Ne[Gamma]), log10(pt_spec->Ne[Gamma]*(x[Gamma] + 1)*(x[Gamma] + 1)*(x[Gamma] + 1)));
                    //fprintf(fp,"%e %e",x[G]+1,pt->Ne[G]);
                    //fprintf(fp,"\n");
                    fprintf(fp, "\n");
                }
            }
            fprintf(fp, "&&\n");
            //fprintf(fp1,"&&\n");
            if (T > 1) {
                fprintf(fp, "&&\n");
                //fprintf(fp1,"&&\n");
            }

            COUNT_FILE*=STEP_FILE;

        }
        //---------------------------------------
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>> ------------------------  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
    }
    //--------------END Loop over Time-----------------------------------------------


    fclose(fp);

    return;
}
