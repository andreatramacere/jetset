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

void show_blob(struct blob pt ) {
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
    printf("nu_grid_size=%d\n", pt.nu_grid_size);
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
    printf("spit_temp=%e\n", pt.spit_temp);
    printf("spit_gamma_th=%e\n", pt.spit_gamma_th);
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
    printf("nu_start_EC_BLR %e\n", pt.nu_start_EC_BLR);
    printf("nu_stop_EC_BLR %e\n", pt.nu_stop_EC_BLR);
    printf("Lum Diks %e\n", pt.L_Disk);
    printf("tau BLR %e\n", pt.tau_BLR);
    printf("R_inner_Sw %e (Rs)\n", pt.R_inner_Sw);
    printf("R_ext_Sw %e (Rs)\n", pt.R_ext_Sw);
    printf("accr eff %e \n", pt.accr_eff);
    printf("T disk max (T max for MultiBB) %e\n", pt.T_Disk);
    printf("dist disk BLR (cm))%e\n", pt.R_BLR_in);
    printf("test array=%e\n", pt.nuF_nu_SSC_obs[0]);
    printf("nu_start_EC_DT %e\n", pt.nu_start_EC_DT);
    printf("nu_stop_EC_DT %e\n", pt.nu_stop_EC_DT);
    printf("T_DT (T Dusty Torus) %e\n", pt.T_DT);
    printf("dist disk DT (cm))%e\n", pt.R_DT);
    printf("tau DT %e\n", pt.tau_DT);
}

void show_temp_ev(  struct temp_ev pt_ev){
    printf("do_Sync_cooling=%d\n", pt_ev.do_Sync_cooling);
    printf("do_Compton_cooling=%d\n", pt_ev.do_Compton_cooling);
    


    printf("L_inj %e (erg/s)\n", pt_ev.L_inj );
    printf("Diff_coeff %e (1/s), t_D=1/Diff_coeff %e (s)\n", pt_ev.Diff_Coeff, pt_ev.t_D0 );
    printf("Acc_coeff %e (1/s),  t_A=1/Acc_coeff %e (s)\n",  pt_ev.Acc_Coeff, pt_ev.t_A0 );
    //printf("T_esc_Coeff %e (R/c)\n", pt_ev.T_esc_Coeff_R_by_c_acc );
    printf("Esc_index %e\n", pt_ev.Esc_Index_acc );
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

    ev_root.do_Sync_cooling = 1;
    ev_root.do_Compton_cooling = 0;
    ev_root.do_Expansion = 0;
    ev_root.do_Adiabatic_cooling = 1;
    ev_root.T_COUNTER=0;



    ev_root.L_inj=1e39;
    //double *T_esc;
    ev_root.t_D0=1.0E4;
    ev_root.t_DA0=ev_root.t_D0*0.5;
    ev_root.t_A0=1.0E3;
    ev_root.Diff_Coeff=1.0/ev_root.t_D0;
    ev_root.Acc_Coeff=1.0/ev_root.t_A0;
    ev_root.Diff_Index=2.0;
    ev_root.Acc_Index=1.0;
    ev_root.Esc_Index_acc=0;
    ev_root.Esc_Index_rad=0.0;
    ev_root.m_B=1.0;
    ev_root.B_rad=1.0;
    ev_root.B_acc=1.0;
    ev_root.B_t=1.0;
    //ev_root.m_R=1.0;
    ev_root.T_esc_Coeff_R_by_c_acc=2.0;
    ev_root.T_esc_Coeff_R_by_c_rad=2.0;
    
    ev_root.TStart_Inj=0.0;
    ev_root.TStop_Inj=3e4;
    ev_root.TStart_Acc=0.0;
    ev_root.TStop_Acc=3e4;
    ev_root.Inj_temp_slope=0.0;
    ev_root.NUM_SET=50;
    ev_root.T_SIZE=1000;
    ev_root.duration=3e4;
    ev_root.E_acc_max=1E200;
    ev_root.Delta_R_acc=1E13;
    //ev_root.R_jet=1E13;
    ev_root.v_exp_by_c=1;
    //ev_root.R_jet_exp=1E13;
    ev_root.t_jet_exp=1E5;
    ev_root.R_jet_t=1E16;
    ev_root.R_H_jet_t=1E17;
    ev_root.R_H_rad_start=1E17;
    ev_root.R_rad_start=1E16;;
    ev_root.gmin_griglia = 1.0e1;
    ev_root.gmax_griglia = 1.0e8;
    ev_root.gamma_grid_size =1E4;
    ev_root.Q_inj_jetset_gamma_grid_size=1E2;

    ev_root.Lambda_max_Turb = 1e15;
    ev_root.Lambda_choer_Turb_factor=0.1;
    ev_root.Gamma_Max_Turb_L_max=Larmor_radius_to_gamma(ev_root.Lambda_max_Turb,0.1, 1.0);
    ev_root.Gamma_Max_Turb_L_coher=Larmor_radius_to_gamma(ev_root.Lambda_max_Turb*ev_root.Lambda_choer_Turb_factor,0.1, 1.0);
    ev_root.LOG_SET = 1;
    ev_root.Q_inj=NULL;
    ev_root.gamma=NULL;
    ev_root.Q_inj_jetset=NULL;
    ev_root.gamma_inj_jetset=NULL;
    //ev_root.N_gamma=NULL;
    
    ev_root.N_rad_gamma=NULL;
    ev_root.N_acc_gamma=NULL;
    ev_root.N_time=NULL;
    ev_root.T_esc_acc=NULL;
    ev_root.T_esc_rad=NULL;
    //ev_root.T_esc_ad_rad=NULL;
    ev_root.T_inj_profile=NULL;
    ev_root.T_acc_profile=NULL;
    return ev_root;
}



struct blob MakeBlob() {

    struct blob spettro_root;
    spettro_root.spec_array_size=static_spec_arr_size;

    spettro_root.WRITE_TO_FILE=0;
    spettro_root.BESSEL_TABLE_DONE=0;
    spettro_root.verbose = 0;
    sprintf(spettro_root.path, "./");
    sprintf(spettro_root.STEM, "TEST");

    sprintf(spettro_root.PARTICLE, "electrons");
    spettro_root.do_Sync = 1;
    spettro_root.Sync_kernel=1;
    spettro_root.do_SSC = 1;
    spettro_root.do_IC=1;
    spettro_root.do_pp_gamma=0;
    spettro_root.do_bremss_ep=0;
    spettro_root.set_pp_racc_elec = 0;
    spettro_root.set_pp_racc_gamma = 0;
    spettro_root.set_pp_racc_nu_mu = 0;
    spettro_root.pp_racc_elec = 1.0;
    spettro_root.pp_racc_gamma = 1.0;
    spettro_root.pp_racc_nu_mu = 1.0;
    spettro_root.E_th_pp_delta_approx=0.1;
    spettro_root.E_pp_x_delta_approx=0.001;
    
    spettro_root.IC_adaptive_e_binning =0;
    spettro_root.do_IC_down_scattering =0;
    sprintf(spettro_root.MODE, "fast");
    //GRID SIZE FOR SEED
    spettro_root.nu_seed_size = 200;
    //GRID SIZE FOR IC
    spettro_root.nu_IC_size = 100;
    spettro_root.gamma_grid_size = 1000;
    spettro_root.gamma_custom_grid_size=1000;
    spettro_root.nu_start_Sync = 1e8;
    spettro_root.nu_stop_Sync = 1e20;
    spettro_root.nu_start_SSC = 1e16;
    spettro_root.nu_stop_SSC = 1e27;
    //GRID SIZE FOR INTERP
    spettro_root.nu_grid_size = 200;
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
    spettro_root.NH_cold_to_rel_e = 0.1;
    spettro_root.N = 10;
    spettro_root.Norm_distr = 1;
    //spettro_root.Norm_distr_L_e_Sync=-1.0;
    spettro_root.Distr_e_done = 0;
    sprintf(spettro_root.DISTR, "lp");
    spettro_root.grid_bounded_to_gamma=1;
    spettro_root.p = 2.0;
    spettro_root.p_1 = 3.0;
    spettro_root.gamma_break = 1.0e4;
    spettro_root.gamma_cut = 1.e4;
    spettro_root.spit_index = 2.4;
    spettro_root.spit_temp = 1.0e3;
    spettro_root.spit_gamma_th = 1.0e4;
    spettro_root.r = 0.4;
    spettro_root.s = 2.0;
    spettro_root.gamma0_log_parab = 1.0e4;
    spettro_root.gammap_log_parab = 1.0e4;
    spettro_root.gamma_inj = 1.0e3;
    spettro_root.gmin = 1.0e1;
    spettro_root.gmax = 1.0e5;
    spettro_root.gmin_secondaries=spettro_root.gmin;
    spettro_root.gmax_secondaries=spettro_root.gmax*mp_by_me;
    spettro_root.gmin_griglia = -1.0;
    spettro_root.gmax_griglia = -1.0;
    spettro_root.gamma_pile_up=1E5;
    spettro_root.gamma_pile_up_cut=5E5;
    spettro_root.alpha_pile_up=1;
    spettro_root.gamma_cooling_eq=0;

    spettro_root.EC_stat=0; 
    spettro_root.EC_stat_orig=0;
    spettro_root.EC_factor=1.0;
    spettro_root.do_EC_Disk = 0;
    spettro_root.do_EC_BLR = 0;
    spettro_root.do_EC_DT = 0;
    spettro_root.do_EC_CMB=0;

    spettro_root.do_EC_Star=0;
    spettro_root.do_Disk=0;
    spettro_root.do_DT=0;

    spettro_root.nu_planck_min_factor=1E-4;
    spettro_root.nu_planck_max_factor=1E2;
    spettro_root.mono_planck_min_factor=0.5;
    spettro_root.mono_planck_max_factor=2.0;
    sprintf(spettro_root.disk_type, "BB");
    spettro_root.R_H=1E17;
    spettro_root.R_H_orig=spettro_root.R_H;
    spettro_root.R_H_scale_factor=20.0;
    spettro_root.R_ext_factor=1.0;
    spettro_root.EC_theta_lim=5.0;
    spettro_root.M_BH = 1E9;

    spettro_root.theta_n_int=50;
    spettro_root.l_n_int=50;
    spettro_root.nu_start_EC_Disk = 1e13;
    spettro_root.nu_stop_EC_Disk = 1e26;
    spettro_root.nu_start_EC_BLR = 1e13;
    spettro_root.nu_stop_EC_BLR = 1e26;
    spettro_root.nu_start_EC_DT = 1e13;
    spettro_root.nu_stop_EC_BLR = 1e26;
    spettro_root.nu_start_EC_CMB = 1e13;
    spettro_root.nu_stop_EC_CMB = 1e30;


    spettro_root.L_Disk = 1e47;
    spettro_root.tau_BLR = 1e-1;
    spettro_root.R_inner_Sw = 3.0;
    spettro_root.R_ext_Sw = 500;
    spettro_root.T_Disk = 1e5;
    spettro_root.T_CMB_0=2.725;
    spettro_root.accr_eff = 0.08;
    spettro_root.R_BLR_in = 1e18;
    spettro_root.R_BLR_out=spettro_root.R_BLR_in*2;
    spettro_root.T_DT = 100;
    spettro_root.R_DT = 5.0e18;
    spettro_root.tau_DT = 1e-1;
    spettro_root.R_Star = 1e10;
    spettro_root.T_Star_max =1E5;

    spettro_root.gam=NULL;
    spettro_root.Q_inj_e_second=NULL;

    spettro_root.Ne=NULL;
    spettro_root.Ne_custom=NULL;
    spettro_root.Ne_IC=NULL;
    spettro_root.Ne_jetset=NULL;
    spettro_root.Ne_stat=NULL;
    
    spettro_root.griglia_gamma_Ne_log=NULL;
    spettro_root.gamma_e_custom=NULL;
    spettro_root.griglia_gamma_Ne_log_stat=NULL;
    spettro_root.griglia_gamma_Ne_log_IC=NULL;
    spettro_root.griglia_gamma_jetset_Ne_log=NULL;

    spettro_root.Np=NULL;
    spettro_root.Np_jetset=NULL;
    spettro_root.Np_custom=NULL;
     
    spettro_root.griglia_gamma_Np_log=NULL;
    spettro_root.griglia_gamma_jetset_Np_log=NULL;
    spettro_root.gamma_p_custom=NULL;
    spettro_root.Integrand_over_gamma_grid=NULL;
    
    return spettro_root;
}


//void MakeNe(struct spettro *pt_base){
//   build_Ne(pt_base);
//}

//=========================================================================================
void set_seed_freq_start(struct blob *pt_base){
    pt_base->nu_start_Sync = 1e6;
    pt_base->nu_stop_Sync = 1e20;
    pt_base->nu_start_SSC = 1e14;
    pt_base->nu_stop_SSC = 1e30;


    pt_base->nu_start_EC_Disk = 1e13;
    pt_base->nu_stop_EC_Disk = 1e30;
    pt_base->nu_start_EC_BLR = 1e13;
    pt_base->nu_stop_EC_BLR = 1e30;
    pt_base->nu_start_EC_DT = 1e13;
    pt_base->nu_start_EC_CMB = 1e13;
    pt_base->nu_stop_EC_CMB = 1e30;
}


//=========================================================================================
void InitRadiative(struct blob *pt_base){
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
        printf("Bessel Functions\n");
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
    pt_base->set_pp_racc_nu_mu = 0;
    pt_base->pp_racc_elec = 1.0;
    pt_base->pp_racc_gamma = 1.0;
    pt_base->pp_racc_nu_mu = 1.0;

    //========================================================

}

void Init(struct blob *pt_base, double luminosity_distance) {
    // if luminosity_distance is negative is evaluated internally
    // otherwise the passed value is used

    //struct spettro *pt_base;
    //double (*pf) (struct spettro *, double);
    
    unsigned int i;
    //char * ENV;
    pt_base->SYSPATH=getenv("BLAZARSED");
    set_seed_freq_start(pt_base);

    //pt_base->emiss_lim=1.0E-120;

    //sprintf(ENV,'%s',getenv("BLAZARSED"));
    //return;
    //printf("CIAO =%s\n",pt_base->SYSPATH);

    //sprintf(pt_base->SYSPATH,'%s', ENV);
    //return;
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
    	if (pt_base->verbose){
            printf("!!! Warning nu_seed_size  was gt spec_array size \n");
            printf("now set  to spec_array size %d \n ",pt_base->spec_array_size-1);
        }
    }
    if (pt_base->nu_IC_size>=pt_base->spec_array_size){
    	pt_base->nu_IC_size=pt_base->spec_array_size-1;
    	if (pt_base->verbose){
            printf("!!! Warning  nu_IC_size  was gt spec_array size \n ");
            printf("now set  to spec_array size %d\n ",pt_base->spec_array_size-1);
        }
    }
    for (i = 0; i < static_spec_arr_size; i++) {
        pt_base->q_comp[i] = 0.0;
        pt_base->j_Sync[i] = 0.0;
        pt_base->j_comp[i] = 0.0;
        pt_base->j_EC[i] = 0.0;
        pt_base->alfa_Sync[i] = 0.0;
        pt_base->I_nu_Sync[i] = 0.0;
        pt_base->j_pp_neutrino_tot[i]=0.0;
        pt_base->j_pp_neutrino_mu[i]=0.0;
        pt_base->j_pp_neutrino_e[i]=0.0;
        pt_base->j_pp_gamma[i]=0.0;
        pt_base->j_bremss_ep[i]=0.0;
    }

    for (i = 0; i < static_spec_arr_size; i++){
        pt_base->nuF_nu_Sync_obs[i]=0.0;
        pt_base->nuF_nu_SSC_obs[i]=0.0;
        pt_base->nuF_nu_EC_Disk_obs[i]=0;
        pt_base->nuF_nu_EC_BLR_obs[i]=0;
        pt_base->nuF_nu_EC_DT_obs[i]=0;
        pt_base->nuF_nu_EC_Star_obs[i]=0;
        pt_base->nuF_nu_EC_CMB_obs[i]=0;
        pt_base->nuF_nu_Disk_obs[i]=0;
        pt_base->nuF_nu_DT_obs[i]=0;
        pt_base->nuF_nu_Star_obs[i]=0;
        pt_base->nuFnu_pp_gamma_obs[i]=0;
        pt_base->nuFnu_pp_neutrino_tot_obs[i]=0;
        pt_base->nuFnu_pp_neutrino_mu_obs[i]=0;
        pt_base->nuFnu_pp_neutrino_e_obs[i]=0;
        pt_base->nuFnu_bremss_ep_obs[i]=0;
    }


    //set file number counter
    pt_base->OUT_FILE = 1;


    InitRadiative(pt_base);

    if (luminosity_distance<0){

        pt_base->dist = dist_lum_cm(pt_base->z_cosm);
    }
    else{
        pt_base->dist = luminosity_distance;
    }

    if (pt_base->verbose) {
        printf("Distanza rigorosa=%e in Mpc \n", pt_base->dist/(1.0e6*1.0e2));
        printf("Distanza rigorosa=%e in cm \n", pt_base->dist);
    }

    pt_base->Distr_e_done = 0;
    pt_base->Distr_e_pp_done = 0;

    if (pt_base->verbose) {     
        printf("******************************  Geometry  *********************************\n");
        printf("Volume for Spherical Geom.=%e\n", pt_base->Vol_sphere);
    }
    
    if (strcmp(pt_base->PARTICLE, "electrons") == 0) {
        InitNe(pt_base);
        pt_base->N_tot_e_Sferic = pt_base->Vol_sphere * pt_base->N_e;
        FindNe_NpGp(pt_base);
        EvalU_e(pt_base);
        
        if (pt_base->verbose) {     
            printf("********************       Leptonic Scenario       ********************\n");
            printf("type of distr=%d\n", pt_base->TIPO_DISTR);
            printf("*******  Leptonic Energetic   **********\n");
            printf("N_e=%e Ne/Ne_0=%e\n", pt_base->N_e, pt_base->N / pt_base->N_0e);
            printf("Total number of electrons    =%e\n", pt_base->N_tot_e_Sferic);
            printf("Gamma_p of N(gamma)*gamma^2 = %e\n", pt_base->Gamma_p2);
            printf("Gamma_p of N(gamma)*gamma^3 = %e\n", pt_base->Gamma_p3);
            printf("Peak of  N(gamma)*gamma^2 = %e\n", pt_base->Np2);
            printf("Peak of  N(gamma)*gamma^3 = %e\n", pt_base-> Np3);
            printf("U_e   blob rest frame =%e erg/cm^3\n", pt_base->U_e);
            printf("U_e/U_b =%e\n", pt_base->U_e / pt_base->UB);
            printf("E_tot (electron)  blob rest frame =%e erg     \n", pt_base->E_tot_e);
            printf("************************************************************************\n");
        }

    } else if (strcmp(pt_base->PARTICLE, "protons") == 0) {
        Init_Np_Ne_pp(pt_base);        
        pt_base->N_tot_p_Sferic = pt_base->Vol_sphere * pt_base->N_p;             
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
 
void Run_SED(struct blob *pt_base){
	//unsigned int i;
    if (pt_base->verbose) {
        printf("STEM=%s\n", pt_base->STEM);
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>> RUN      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    }
    //==================================================
    // Evaluate hadronic pp Spectrum
    //==================================================
    if ((strcmp(pt_base->PARTICLE, "protons") == 0) && pt_base->do_pp_gamma) {

        spettro_pp_gamma(1, pt_base);
        spettro_pp_neutrino(1,pt_base);
    }

    if ((strcmp(pt_base->PARTICLE, "protons") == 0) && pt_base->do_pp_neutrino) {
        spettro_pp_neutrino(1,pt_base);
    }

    //==================================================
    // Evaluate Bremms sp Spectrum
    //==================================================
    if (pt_base->do_bremss_ep) {
        spettro_bremss_ep(1, pt_base);
    }


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
		if (pt_base->do_EC_Disk == 1 || pt_base->do_EC_BLR == 1 || pt_base->do_EC_DT == 1  || pt_base->do_EC_Star == 1 || pt_base->do_EC_CMB == 1 || pt_base->do_Disk==1 || pt_base->do_DT==1 || pt_base->do_Star==1) 
        {
                
                spectra_External_Fields(1, pt_base);
                if (pt_base->do_EC_Star == 1) {
                    //if (pt_base->verbose) {
                    //    printf("************* Disk ****************\n");
                    //}
                    pt_base->EC = 4;
                    spettro_EC(1, pt_base);
                }
                if (pt_base->do_EC_Disk == 1 || pt_base->do_Disk==1) {
                    //if (pt_base->verbose) {
                    //    printf("************* Disk ****************\n");
                   // }
                    pt_base->EC = 1;
                    spettro_EC(1, pt_base);
                }
                if (pt_base->do_EC_BLR == 1) {
                    //if (pt_base->verbose) {
                    //    printf("************* BLR ****************\n");
                    //}
                    pt_base->EC = 2;
                    spettro_EC(1, pt_base);
                }
                if (pt_base->do_EC_DT == 1) {
                    //if (pt_base->verbose) {
                    //    printf("************* DT ****************\n");
                    // }
                    pt_base->EC = 3;
                    spettro_EC(1, pt_base);
                }
                if (pt_base->do_EC_CMB == 1) {
                    //if (pt_base->verbose) {
                    //    printf("************* CMB ****************\n");
                   // }
                    pt_base->EC = 5;
                    spettro_EC(1, pt_base);
                }
                //if (pt_base->do_EC_CMB_stat == 1) {
                //    if (pt_base->verbose) {
                //       printf("************* CMB stat ****************\n");
                //   }
                //       printf("************* CMB stat ****************\n");
                //    pt_base->EC = 6;
                //    spettro_EC(1, pt_base);
                //}
            }
        //printf("=>done\n");
    }
    //==================================================
    //Sum Up all the Spectral Components
    //==================================================
    common_grid_spectra(1, pt_base);

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
double get_array(double * arr, unsigned int id, unsigned int size){
	if ((id >=0) && (id <= size)){
		return arr[id];
	}
	else{
        printf("exceeded array size in get_spectral_array\n");
        exit(0);
	}
}



double get_spectral_array(double * arr, struct blob * pt, unsigned int id){
	if ((id >=0) && (id <= pt->nu_grid_size)){
		return arr[id];
	}
	else{
        printf("exceeded array size in get_spectral_array\n");
        exit(0);
	}
}


double get_elec_array(double * arr, struct blob *pt, unsigned int id){
	if ((id>=0) && (id<=pt->gamma_grid_size)){
		return arr[id];
	}
	else{
        printf("exceeded array size in get_elec_array\n");
        exit(0);
	}
}

double get_Q_inj_array(double *arr, struct temp_ev *pt_ev, unsigned int id)
{
    if ((id >= 0) && (id <= pt_ev->gamma_grid_size))
    {
        return arr[id];
    }
    else
    {
        printf("exceeded array size in get_Q_inj_array\n");
        exit(0);
    }
}

double get_temp_ev_N_gamma_array(double *arr, struct temp_ev *pt_ev, unsigned int row, unsigned int col)
{
    if (((col >= 0) && (col <= pt_ev->gamma_grid_size)) && ((row >= 0) && (row <= pt_ev->NUM_SET)))
    {
        return arr[row * pt_ev->gamma_grid_size + col];
    }
    else
    {
        printf("exceeded array size in get_temp_ev_N_gamma_array\n");
        exit(0);
    }
}

double get_temp_ev_N_time_array(double *arr, struct temp_ev *pt_ev, unsigned int id)
{
    if ((id >= 0) && (id <= pt_ev->NUM_SET))
    {
            return arr[id];
        }
        else
        {
            printf("exceeded array size in get_temp_ev_N_time_array\n");
            exit(0);
        }
    }

double get_temp_ev_gamma_array(double *arr, struct temp_ev *pt_ev, unsigned int id)
    {
        if ((id >= 0) && (id <= pt_ev->gamma_grid_size))
        {
            return arr[id];
        }
        else
        {
            printf("exceeded array size in get_temp_ev_gamma_array\n");
            exit(0);
        }
    }

void set_elec_array(double * arr,struct blob *pt, double val, unsigned int id){
    if ((id>=0) && (id<=pt->gamma_grid_size)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_elec_array\n");
            exit(0);
        }
}

void set_q_inj_user_array(double * arr,struct temp_ev *pt, double val, unsigned int id){
    if ((id>=0) && (id<=pt->Q_inj_jetset_gamma_grid_size)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_elec_array\n");
            exit(0);
        }
}

void set_temp_ev_Time_array(double * arr,struct temp_ev *pt, double val, unsigned int id){
    if ((id>=0) && (id<=pt->T_SIZE)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_T_inj_profile\n");
            exit(0);
        }
}


void set_elec_custom_array(double * arr, struct blob *pt,double val, unsigned int id){
    if ((id>=0) && (id<=pt->gamma_custom_grid_size)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_elec_custom_array\n");
            exit(0);
        }
}

void set_bessel_table(double *arr, struct blob *pt, double val, unsigned int id)
{
    if ((id >= 0) && (id <= static_bess_table_size))
    {
        arr[id] = val;
    }
    else
    {
        printf("exceeded array size in set_bessel_table\n");
        exit(0);
    }
}

double get_temp_ev_array_static(double *arr, unsigned int id){
    if ((id >= 0) && (id <= static_ev_arr_grid_size))
    {
        return arr[id];
    }
    else
    {
        printf("exceeded array size get_temp_ev_array_static\n");
        exit(0);
    }
}
//=========================================================================================
void SetBeaming(struct blob *pt){

	if (strcmp(pt->BEAMING_EXPR, "delta") == 0) {
	        pt->beam_obj = pt->beam_obj;
            pt->BulkFactor = pt->beam_obj;
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



