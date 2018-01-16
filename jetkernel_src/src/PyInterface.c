#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"

/**
 * \file Blazar_SED.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief lettura input file
 * chimata sotto funzioni
 *
 */

//=========================================================================================

void show_blob(struct spettro pt ) {
    printf("verbose=%d\n", pt.verbose);
    printf("path=%s\n", pt.path);
    printf("STEM=%s\n", pt.STEM);
    printf("do_Sync=%d\n", pt.do_Sync);
    printf("do_SSC=%d\n", pt.do_SSC);
    printf("MODE=%s\n", pt.MODE);
    printf("PARTICLE=%s\n", pt.PARTICLE);
    printf("nu_seed_size=%d\n", pt.nu_seed_size);
    printf("nu_IC_size=%d\n", pt.nu_IC_size);
    printf("nu_start_Sync=%e\n", pt.nu_start_Sync);
    printf("nu_stop_Sync=%e\n", pt.nu_stop_Sync);
    printf("nu_start_SSC=%e\n", pt.nu_start_SSC);
    printf("nu_stop_SSC=%e\n", pt.nu_stop_SSC);
    printf("nu_sum_size=%d\n", pt.nu_sum_size);
    printf("nu_start_grid=%e\n", pt.nu_start_grid);
    printf("nu_stop_grid=%e\n", pt.nu_stop_grid);
    printf("B=%e\n", pt.B);
    printf("R=%e\n", pt.R);
    printf("BulkFactor=%e\n", pt.BulkFactor);
    printf("theta=%e\n", pt.theta);
    printf("z_cosm=%e\n", pt.z_cosm);
    printf("NH_pp=%e\n", pt.NH_pp);
    printf("N=%e\n", pt.N);
    printf("Norm_distr=%d\n", pt.Norm_distr);
    printf("DISTR=%s\n", pt.DISTR);
    printf("p=%e\n", pt.p);
    printf("p_1=%e\n", pt.p_1);
    printf("gamma_break=%e\n", pt.gamma_break);
    printf("gamma_cut=%e\n", pt.gamma_cut);
    printf("spit_index=%e\n", pt.spit_index);
    printf("spit_ratio=%e\n", pt.spit_ratio);
    printf("spit_cut=%e\n", pt.spit_cut);
    printf("spit_cut1=%e\n", pt.spit_cut1);
    printf("spit_cut2=%e\n", pt.spit_cut2);
    printf("r=%e\n", pt.r);
    printf("s=%e\n", pt.s);
    printf("gamma0_log_parab(lp,lppl)=%e\n", pt.gamma0_log_parab);
    printf("gammap_log_parab(lpep)=%e\n", pt.gammap_log_parab);
    printf("gmin=%e\n", pt.gmin);
    printf("gmax=%e\n", pt.gmax);
    printf("do_EC_Disk=%d\n", pt.do_EC_Disk);
    printf("do_EC_BLR=%d\n", pt.do_EC_BLR);
    printf("do_EC_DT=%d\n", pt.do_EC_DT);
    printf("disk type =%s\n", pt.disk_type);
    printf("nu_start_EC_BLR\n", pt.nu_start_EC_BLR);
    printf("nu_stop_EC_BLR\n", pt.nu_stop_EC_BLR);
    printf("Lum Diks %e\n", pt.L_disk);
    printf("tau BLR %e\n", pt.tau_BLR);
    printf("R_inner_Sw %e (Rs)\n", pt.R_inner_Sw);
    printf("R_ext_Sw %e (Rs)\n", pt.R_ext_Sw);
    printf("accr eff %e \n", pt.accr_eff);
    printf("T disk max (T max for MultiBB) %e\n", pt.T_disk_max);
    printf("dist disk BLR (cm))%e\n", pt.R_BLR_in);
    printf("test array=%e\n", pt.nuF_nu_SSC_obs[0]);
    printf("nu_start_EC_DT\n", pt.nu_start_EC_DT);
    printf("nu_stop_EC_DT\n", pt.nu_stop_EC_DT);
    printf("T_DT (T Dusty Torus) %e\n", pt.T_DT);
    printf("dist disk DT (cm))%e\n", pt.R_DT);
    printf("tau DT %e\n", pt.tau_DT);
}

void show_temp_ev(  struct temp_ev pt_ev){
    printf("do_Sync_cooling=%d\n", pt_ev.do_Sync_cooling);
    printf("do_SSC_cooling=%d\n", pt_ev.do_SSC_cooling);
    printf("do_EC_cooling_Disk=%d\n", pt_ev.do_EC_cooling_Disk);
    printf("do_EC_cooling_BLR=%d\n", pt_ev.do_EC_cooling_BLR);
    printf("do_EC_cooling_DT=%d\n", pt_ev.do_EC_cooling_DT);
    printf("do_EC_cooling_Star=%d\n", pt_ev.do_EC_cooling_Star);
    printf("do_EC_cooling_CMB=%d\n", pt_ev.do_EC_cooling_CMB);



    printf("L_inj %e (erg/s)\n", pt_ev.L_inj );
    printf("Diff_coeff %e (1/s), t_D=1/Diff_coeff %e (s)\n", pt_ev.Diff_Coeff, pt_ev.t_D0 );
    printf("Acc_coeff %e (1/s),  t_A=1/Acc_coeff %e (s)\n",  pt_ev.Acc_Coeff, pt_ev.t_A0 );
    printf("T_esc_Coeff %e (R/c)\n", pt_ev.T_esc_Coeff );
    printf("Esc_index %e\n", pt_ev.Esc_Index );
    printf("T_start_Acc %e (s)\n", pt_ev.TStart_Acc );
    printf("T_stop_Acc  %e (s)\n", pt_ev.TStop_Acc );
    printf("T_start_inj %e (s)\n", pt_ev.TStart_Inj );
    printf("T_stop_inj  %e (s)\n", pt_ev.TStop_Inj );
    printf("Num. out file  %d \n", pt_ev.NUM_SET) ;
    printf("T size %d \n", pt_ev.T_SIZE);
}
//=========================================================================================


//=========================================================================================
struct temp_ev MakeTempEv() {
    struct temp_ev ev_root;
    //ev_root.t_unit=; //unit time in light crossing time
    //double t_acc;

    ev_root.do_Sync_cooling = 0;
    ev_root.do_SSC_cooling = 0;
    ev_root.do_EC_cooling_Disk = 0;
    ev_root.do_EC_cooling_BLR = 0;
    ev_root.do_EC_cooling_DT = 0;
    ev_root.do_EC_cooling_Star = 0;
    ev_root.do_EC_cooling_CMB  = 0;





    ev_root.L_inj=1e39;
    //double *T_esc;
    ev_root.t_D0=1.0E4;
    ev_root.t_DA0=ev_root.t_D0*0.5;
    ev_root.t_A0=1.0E3;
    ev_root.Diff_Coeff=1.0/ev_root.t_D0;
    ev_root.Acc_Coeff=1.0/ev_root.t_A0;
    ev_root.Diff_Index=2.0;
    ev_root.Acc_Index=1.0;
    ev_root.T_esc_Coeff=2.0;
    ev_root.Esc_Index=0.0;
    ev_root.TStart_Inj=0.0;
    ev_root.TStop_Inj=3e4;
    ev_root.TStart_Acc=0.0;
    ev_root.TStop_Acc=3e4;
    ev_root.NUM_SET=50;
    ev_root.T_SIZE=1000;
    ev_root.duration=3e4;

    ev_root.Lambda_max_Turb= 1e15;
    ev_root.Lambda_choer_Turb_factor=0.1;
    ev_root.Gamma_Max_Turb_L_max=Larmor_radius_to_gamma(ev_root.Lambda_max_Turb,0.1, 1.0);
    ev_root.Gamma_Max_Turb_L_coher=Larmor_radius_to_gamma(ev_root.Lambda_max_Turb*ev_root.Lambda_choer_Turb_factor,0.1, 1.0);
    return ev_root;
}



struct spettro MakeBlob() {

    struct spettro spettro_root;
    spettro_root.spec_array_size=static_spec_arr_size;
    spettro_root.BESSEL_TABLE_DONE=0;
    spettro_root.verbose = 0;
    sprintf(spettro_root.path, "./");
    sprintf(spettro_root.STEM, "TEST");

    sprintf(spettro_root.PARTICLE, "leptons");
    spettro_root.do_Sync = 1;
    spettro_root.Sync_kernel=1;
    spettro_root.do_SSC = 1;
    spettro_root.do_IC=1;

    sprintf(spettro_root.MODE, "fast");

    spettro_root.nu_seed_size = 200;
    spettro_root.nu_IC_size = 100;
    spettro_root.gamma_grid_size = 1000;
    spettro_root.nu_start_Sync = 1e8;
    spettro_root.nu_stop_Sync = 1e20;
    spettro_root.nu_start_SSC = 1e16;
    spettro_root.nu_stop_SSC = 1e27;
    spettro_root.nu_sum_size = 200;
    spettro_root.nu_start_grid = 1e8;
    spettro_root.nu_stop_grid = 1e27;
    spettro_root.emiss_lim=1.0E-120;
    spettro_root.B = 0.1;
    spettro_root.sin_psi = 1.0;
    spettro_root.R = 1e15;
    sprintf(spettro_root.BEAMING_EXPR, "delta");
    spettro_root.BulkFactor = 10;
    spettro_root.beta_Gamma=eval_beta_gamma(spettro_root.BulkFactor);
    spettro_root.theta = 3.5;
    spettro_root.beam_obj=10.0;
    spettro_root.z_cosm = 0.1;
    spettro_root.NH_pp = 0.1;
    spettro_root.N = 10;
    spettro_root.Norm_distr = 1;
    spettro_root.Norm_distr_L_e_Sync=-1.0;
    spettro_root.Distr_e_done = 0;
    sprintf(spettro_root.DISTR, "lp");
    spettro_root.p = 2.0;
    spettro_root.p_1 = 3.0;
    spettro_root.gamma_break = 1.0e4;
    spettro_root.gamma_cut = 1.e4;
    spettro_root.spit_index = 2.4;
    spettro_root.spit_ratio = 1.0e12;
    spettro_root.spit_cut = 1.0e4;
    spettro_root.spit_cut1 = 1.0e6;
    spettro_root.spit_cut2 = 3.0e6;
    spettro_root.r = 0.4;
    spettro_root.s = 2.0;
    spettro_root.gamma0_log_parab = 1.0e4;
    spettro_root.gammap_log_parab = 1.0e4;
    spettro_root.gmin = 1.0e1;
    spettro_root.gmax = 1.0e5;
    spettro_root.gmin_griglia = -1.0;
    spettro_root.gmax_griglia = -1.0;
    
    
    spettro_root.do_EC_Disk = 0;
    spettro_root.do_EC_BLR = 0;
    spettro_root.do_EC_DT = 0;
    spettro_root.do_EC_CMB=0;
    spettro_root.do_EC_CMB_stat=0;
    spettro_root.do_EC_Star=0;

    spettro_root.nu_planck_min_factor=1E-4;
    spettro_root.nu_planck_max_factor=1E2;
    spettro_root.mono_planck_min_factor=0.5;
    spettro_root.mono_planck_max_factor=2.0;
    sprintf(spettro_root.disk_type, "BB");
    spettro_root.R_H=1E17;
    spettro_root.nu_start_EC_Disk = 1e13;
    spettro_root.nu_stop_EC_Disk = 1e26;
    spettro_root.nu_start_EC_BLR = 1e13;
    spettro_root.nu_stop_EC_BLR = 1e26;
    spettro_root.nu_start_EC_DT = 1e13;
    spettro_root.nu_stop_EC_BLR = 1e26;
    spettro_root.nu_start_EC_CMB = 1e13;
    spettro_root.nu_stop_EC_CMB = 1e30;


    spettro_root.L_disk = 1e47;
    spettro_root.tau_BLR = 1e-1;
    spettro_root.R_inner_Sw = 3.0;
    spettro_root.R_ext_Sw = 500;
    spettro_root.T_disk_max = 1e5;
    spettro_root.T_CMB_0=2.725;
    spettro_root.accr_eff = 0.1;
    spettro_root.R_BLR_in = 1e18;
    spettro_root.R_BLR_out=spettro_root.R_BLR_in*2;
    spettro_root.T_DT = 100;
    spettro_root.R_DT = 5.0e18;
    spettro_root.tau_DT = 1e-1;
    double test[1000];
    return spettro_root;
}
//=========================================================================================




//=========================================================================================

void Init(struct spettro *pt_base) {
    //struct spettro *pt_base;
    double (*pf) (struct spettro *, double);
    double test, prova;
    unsigned long i;
    char *SYSPATH;

    SYSPATH = getenv("BLAZARSED");
    sprintf(pt_base->SYSPATH,'%s', SYSPATH);

    if (pt_base->verbose) {
        printf("SYSPATH =%s\n", pt_base->SYSPATH);
        printf("STEM=%s\n", pt_base->STEM);
        printf("PATH =%s\n", pt_base->path);
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>> Satic Case Initilization <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    }
    //======================================
    // Arrays SetUp
    //======================================
    if (pt_base->nu_seed_size>=pt_base->spec_array_size){
    	pt_base->nu_seed_size=pt_base->spec_array_size-1;
        printf("nu_seed_size  was gt spec_array size \n ");
        printf("now set  to spec_array size %d \n ",pt_base->spec_array_size-1);

    }
    if (pt_base->nu_IC_size>=pt_base->spec_array_size){
    	pt_base->nu_IC_size=pt_base->spec_array_size-1;
    	printf("nu_IC_size  was gt spec_array size \n ");
    	printf("now set  to spec_array size %d\n ",pt_base->spec_array_size-1);

    }
    for (i = 0; i < pt_base->spec_array_size; i++) {
        pt_base->q_comp[i] = 0.0;
        pt_base->j_Sync[i] = 0.0;
        pt_base->alfa_Sync[i] = 0.0;
        pt_base->I_nu_Sync[i] = 0.0;
    }
    //printf("path=%s", spettro_root.path);

    //set file number counter
    pt_base->OUT_FILE = 1;



    //==========================================
    //Numerical Integration precision Setup
    //==========================================
    if (strcmp(pt_base->MODE, "accurate") == 0) {
        pt_base->gamma_grid_size = 10000;
        if (pt_base->verbose) {
            printf("gamma mesh set to value=%d for accurate integration \n",pt_base->gamma_grid_size);
        }
    }
    else if (strcmp(pt_base->MODE, "fast") == 0) {
        pt_base->gamma_grid_size = 1000;
        if (pt_base->verbose) {
            printf("gamma mesh set to value=%d for fast integration, \n",pt_base->gamma_grid_size);
        }
    }
    else  if (strcmp(pt_base->MODE, "custom") == 0) {
    	if (pt_base->verbose) {
    		printf("gamma mesh set to custom value=%d  \n",pt_base->gamma_grid_size);
    	}
    }
    else {
    	if (pt_base->verbose) {
			printf("MODE set to wrong value: %s, allowed= accurate,fast,custom",pt_base->MODE);
			exit(1);
      	}
    }

	if (fmod((double) pt_base->gamma_grid_size, 2.0) == 0) {
		pt_base->gamma_grid_size++;
		if (pt_base->verbose) {
			printf("!!!!!! gamma_grid_size has to be odd\n");
			printf("!!!!!!! pt->gamma_grid_size=%d\n", pt_base->gamma_grid_size);
		}
	}

    




    //=========================================
    // check on gamma grid
    //=========================================
    // gamma min griglia
    if (pt_base->gmin_griglia<0.0 || pt_base->gmin < pt_base->gmin_griglia ){
		if(pt_base->gmin>2.0){
			pt_base->gmin_griglia=pt_base->gmin/2.0;
		}
		else{
		   pt_base->gmin_griglia=1.0;
		}
    }


    if (pt_base->gmax_griglia<0.0 || pt_base->gmax > pt_base->gmax_griglia ){
    	pt_base->gmax_griglia=pt_base->gmax;
    }

    if (pt_base->gmin < pt_base->gmin_griglia ) {
        printf("gmin < gmin_griglia, it must be the oppsosite");
        exit(1);
    } 
    if (pt_base->gmax > pt_base->gmax_griglia ) {
        printf("gmax > gmax_griglia, it must be the oppsosite");
        exit(1);
    }

    //========================================================
    // Geometry Setup
    //========================================================
    pt_base->Vol_sphere = V_sphere(pt_base->R);
    pt_base->Surf_sphere = S_sphere(pt_base->R);
    SetBeaming(pt_base);
    pt_base->beta_Gamma=eval_beta_gamma(pt_base->BulkFactor);


    
    //========================================================
    // Synchrotron Parameter Initialization
    //========================================================
    pt_base->nu_B = (q_esu * pt_base->B) / (2 * pi * me_g * vluce_cm);
    pt_base->UB = pow(pt_base->B, 2.0) / (8.0 * pi); /*dens. ener. B */
    
    if (pt_base->verbose>0) {
        printf("gmin %e   gmax %e \n", pt_base->gmin, pt_base->gmax);
        printf("UB=%e \n", pt_base->UB);
        printf("nu_B_non_rel=%e \n", pt_base->nu_B);
        printf("beaming factor =%e\n", pt_base->beam_obj);
    }
    
    //COSTANTI PER ALFA=FIXED E KERNEL DELTA O KERNEL 2
    pt_base->C1_Sync_K53 = pow(3, 0.5) * pow(q_esu, 3.0) * pt_base->sin_psi;
    pt_base->C1_Sync_K53 *= pt_base->B / (MEC2) * one_by_four_pi;
    pt_base->C2_Sync_K53 = 2.0/(3*pt_base->nu_B);

    pt_base->C1_Sync_K_AVE= 4*pi*pow(3, 0.5) * pow(q_esu, 2.0)*pt_base->nu_B/(vluce_cm) * one_by_four_pi;
    pt_base->C2_Sync_K_AVE=1.0/(3*pt_base->nu_B);

    pt_base->C3_Sync_K53 = -1.0 * pow(3, 0.5) * pow(q_esu, 3.0) / (8 * pi * MEC2 * me_g);

    pt_base->COST_Sync_COOLING = SIGTH * vluce_cm/MEC2;


    //==================================
    //  Bessel Fucntion Setup
    //==================================
    //exit(1);
    if (pt_base->BESSEL_TABLE_DONE == 0){
    	tabella_Bessel(pt_base);
    }
    
    //========================================================
    // Compton Parameter Initialization
    //========================================================
    pt_base->COST_IC_K1 = 3.0 * SIGTH * vluce_cm / 4.0;
        pt_base->COST_IC_COOLING = (4.0/3.0) * SIGTH * vluce_cm*HPLANCK/MEC2;


    //========================================================
    // pp parameters initialization
    //========================================================
    pt_base->set_pp_racc_elec = 0;
    pt_base->set_pp_racc_gamma = 0;
    pt_base->pp_racc_elec = 1.0;
    pt_base->pp_racc_gamma = 1.0;


    //========================================================
    // Cosmological Setup
    //========================================================
    // Calcolo approssimato della distanza
    //pt_base->dist = Distanza_Lum_analyt(pt_base->z_cosm);
    //printf("Distanza Approssimativa=%e \n", pt_base->dist);
    // Calcolo rigoroso della distanza
    //pf = &distanza_z;
    //pt_base->dist = (vluce_cm * 1.0e-5 / H_0)*(1.0 + pt_base->z_cosm) * integrale_simp_struct(pf, pt_base, 0, pt_base->z_cosm, 10000);
    //pt_base->dist *= parsec * 1.0e6 * 1.0e2;
    pt_base->dist=dist_lum_cm(pt_base->z_cosm);
    
    if (pt_base->verbose) {
        printf("Distanza rigorosa=%e in Mpc \n", pt_base->dist/(1.0e6*1.0e2));
        printf("Distanza rigorosa=%e in cm \n", pt_base->dist);
    }

    //============================================
    // Numerical Tets
    //=============================================
    //pf = &test_int;
    //printf("************* test numerici **************** \n");
    //test=integrale_simp_log_struct(pf, pt_base, 10, 1e20, 10000);
    //printf("test log=%e\n", test);
    //test=integrale_simp_struct(pf, pt_base, 10, 1e20, 10000);
    //printf("test lin=%e\n", test);
    //printf("ris atteso=%e\n", (pow(1e20, 3)/3)-pow(10, 3)/3);



    //============================================
    //      Distre & Norm
    //============================================
    //Set Distre_e_done=0 significa che la distre nn e' stata
    //ancora calcolata


    pt_base->Distr_e_done = 0;
    pt_base->Distr_e_pp_done = 0;

    if (pt_base->verbose) {     
        printf("******************************  Geometry  *********************************\n");
        printf("Volume for Spherical Geom.=%e\n", pt_base->Vol_sphere);
    }
    
    if (strcmp(pt_base->PARTICLE, "leptons") == 0) {
        
        Genera_Ne(pt_base);
        pt_base->N_tot_e_Sferic = pt_base->Vol_sphere * pt_base->N;
        FindNe_NpGp(pt_base);
        EvalU_e(pt_base);
        
        if (pt_base->verbose) {     
            printf("********************       Leptonic Scenario       ********************\n");
            printf("type of distr=%d\n", pt_base->TIPO_DISTR);
            printf("*******  Leptonic Energetic   **********\n");
            printf("N_e=%e Ne/Ne_0=%e\n", pt_base->N, pt_base->N / pt_base->N_0e);
            printf("Total number of electrons    =%e\n", pt_base->N_tot_e_Sferic);
            printf("Gamma_p of N(gamma)*gamma^2 = %e\n", pt_base->Gamma_p2);
            printf("Gamma_p of N(gamma)*gamma^3 = %e\n", pt_base->Gamma_p3);
            printf("Peak of  N(gamma)*gamma^2 = %e\n", pt_base->Np2);
            printf("Peak of  N(gamma)*gamma^3 = %e\n", pt_base-> Np3);
            printf("U_e   blob rest frame =%e erg/cm^3\n", pt_base->U_e);
            printf("U_e/U_b =%e\n", pt_base->U_e / pt_base->UB);
            printf("E_tot (electron)  blob rest frame =%e erg     \n", pt_base->E_tot_e);
        }

    } else if (strcmp(pt_base->PARTICLE, "hadrons") == 0) {
        
        Genera_Np_Ne_pp(pt_base);        
        pt_base->N_tot_p_Sferic = pt_base->Vol_sphere * pt_base->N;             
        EvalU_p(pt_base);             
        pt_base->N_tot_e_Sferic = pt_base->Vol_sphere * pt_base->N_e_pp;
        EvalU_e(pt_base);
        FindNe_NpGp(pt_base);
             
        if (pt_base->verbose) {

            printf("***********************       Hadronic Scenario           ********************\n");
            printf("****** Generate Np and Ne form secondaries **************\n");
            printf("******************      Hadronic Energetic    **********\n");
            printf("N_p=%e N_p/N_0p=%e\n", pt_base->N, pt_base->N / pt_base->N_0p);
            printf("Total number of p    =%e\n", pt_base->N_tot_p_Sferic);
            printf("U_p   blob rest frame =%e erg/cm^3\n", pt_base->U_p);
            printf("U_p/U_b =%e\n", pt_base->U_p / pt_base->UB);
            printf("E_tot (protons)  blob rest frame =%e erg     \n", pt_base->E_tot_p);
            printf("*******  Scondaries Leptonic Energetic  ****************\n");
            printf("In this case N_0e=1, leptons come from pp, N_e_pp=%e\n", pt_base->N_e_pp);
            printf("Total number of secondary electrons    =%e\n", pt_base->N_tot_e_Sferic);
            printf("U_e   blob rest frame =%e erg/cm^3\n", pt_base->U_e);
            printf("U_e/U_b =%e\n", pt_base->U_e / pt_base->UB);
            printf("E_tot (electron)  blob rest frame =%e erg     \n", pt_base->E_tot_e);
            printf("Gamma_p of N(gamma)*gamma^2 = %e\n", pt_base->Gamma_p2);
            printf("Gamma_p of N(gamma)*gamma^3 = %e\n", pt_base->Gamma_p3);
            printf("Peak of  N(gamma)*gamma^2 = %e\n", pt_base->Np2);
            printf("Peak of  N(gamma)*gamma^3 = %e\n", pt_base-> Np3);
        }
    }

}
void Run_SED(struct spettro *pt_base){
	unsigned long i;
    if (pt_base->verbose) {
        printf("STEM=%s\n", pt_base->STEM);
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>> RUN      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    }
    //==================================================
    // Evaluate hadronic pp Spectrum
    //==================================================
    if (strcmp(pt_base->PARTICLE, "hadrons") == 0) {

        spettro_pp(1, pt_base);
    }


    //graph_integranda_j_nu_Sync(pt_base);
    //==================================================
    // Evaluate Synchrotron Spectrum
    //==================================================
    if (pt_base->do_Sync != 0) {
        spettro_sincrotrone(1, pt_base);
    }


    //==================================================
    // Evaluate SSC Spectrum
    //==================================================
    if (pt_base->do_SSC && pt_base->do_IC) {
        spettro_compton(1, pt_base);
    }


    //==================================================
    // Evaluate EC Spectrum
    //==================================================
	if (pt_base->do_IC) {
		if (pt_base->do_EC_Disk == 1 || pt_base->do_EC_BLR == 1
				|| pt_base->do_EC_DT == 1  || pt_base->do_EC_Star == 1
				|| pt_base->do_EC_CMB == 1 || pt_base->do_EC_CMB_stat) {
			spectra_External_Fields(1, pt_base);
			if (pt_base->do_EC_Star == 1) {
				if (pt_base->verbose) {
					printf("************* Disk ****************\n");
				}
				pt_base->EC = 4;
				spettro_EC(1, pt_base);
			}
			if (pt_base->do_EC_Disk == 1) {
				if (pt_base->verbose) {
					printf("************* Disk ****************\n");
				}
				pt_base->EC = 1;
				spettro_EC(1, pt_base);
			}
			if (pt_base->do_EC_BLR == 1) {
				if (pt_base->verbose) {
					printf("************* BLR ****************\n");
				}
				pt_base->EC = 2;
				spettro_EC(1, pt_base);
			}
			if (pt_base->do_EC_DT == 1) {
				if (pt_base->verbose) {
					printf("************* DT ****************\n");
				}
				pt_base->EC = 3;
				spettro_EC(1, pt_base);
			}
			if (pt_base->do_EC_CMB == 1) {
				if (pt_base->verbose) {
					printf("************* CMB ****************\n");
				}
				pt_base->EC = 5;
				spettro_EC(1, pt_base);
			}
            if (pt_base->do_EC_CMB_stat == 1) {
                if (pt_base->verbose) {
                    printf("************* CMB stat ****************\n");
                }
                printf("************* CMB stat ****************\n");
                pt_base->EC = 6;
                spettro_EC(1, pt_base);
            }
		}
	}
    //==================================================
    //Sum Up all the Spectral Components
    //==================================================
    spettro_somma_Sync_ic(1, pt_base);

    //==================================================
    // Energetic
    //==================================================
    //printf("Energetic computation (output to file)\n");


    //EnergeticOutput(pt_base);
    //CoolingRates(pt_base);



}

//=========================================================================================

//==================================================
//Funtions To access Ne and Spectral components form Python
//==================================================


double get_freq_array(double * arr, struct spettro * pt, unsigned long id){
	if ((id >=0) && (id <= pt->spec_array_size)){
		return arr[id];
	}
	else{
		printf("exceeded array size\n");
		exit(0);
	}
}


double get_elec_array(double * arr, struct spettro *pt, unsigned long id){
	if ((id>=0) && (id<=pt->gamma_grid_size)){
		return arr[id];
	}
	else{
		printf("exceeded array size\n");
		exit(0);
	}
}


//=========================================================================================



void SetBeaming(struct spettro *pt){

	if (strcmp(pt->BEAMING_EXPR, "delta") == 0) {
	        pt->beam_obj = pt->beam_obj;
    }

	else if (strcmp(pt->BEAMING_EXPR, "bulk_theta") == 0) {
		pt->beam_obj = get_beaming(pt->BulkFactor,pt->theta);
	}

	else {
		printf("BEAMING_EXPR variable set to wrong value, posible delta or bulk_theta \n");
		exit(0);
	}

	if (pt->verbose) {
	     printf("beaming set to  %e\n",pt->beam_obj);
	}
}

//=========================================================================================

void SetDistr(struct spettro *pt) {
    //8 is for secondary e- coming from pp

    /*** Associo ad ogni distribuzione di elettroni ***/
    if (strcmp(pt->DISTR, "pl") == 0) {
        pt->TIPO_DISTR = 0;
       
    }
    
    if (strcmp(pt->DISTR, "plc") == 0) {
        pt->TIPO_DISTR = 1;
    }
    
    if (strcmp(pt->DISTR, "bkn") == 0) {
        pt->TIPO_DISTR = 2;
    }
    
    if (strcmp(pt->DISTR, "lp") == 0) {
        pt->TIPO_DISTR = 3;
    }

    if (strcmp(pt->DISTR, "lpep") == 0) {
        pt->TIPO_DISTR = 4;
    }
    
    if (strcmp(pt->DISTR, "lppl") == 0) {
        pt->TIPO_DISTR = 5;
    }
    
    if (strcmp(pt->DISTR, "file") == 0) {
        pt->TIPO_DISTR = 6;
    }
    
    if (strcmp(pt->DISTR, "spitkov") == 0) {
        pt->TIPO_DISTR = 7;
    }
        
    if (pt->verbose) {
     printf("tipo di distribuzione %d\n",pt->TIPO_DISTR);
    }
}
//=========================================================================================
