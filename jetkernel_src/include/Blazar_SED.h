

#ifndef SED_H
#define SED_H
#define G 6.673e-11 //SI
#define G_cgs 6.673e-8 //cgs
#define parsec    3.0856775807e16 /*parsec espresso in metri*/
#define H_0                     71.0 /*Km*sec^-1*Mpc^-1*/
#define q_0        0.55
#define Omega_lambda        0.73
#define Omega_matter        0.27
#define Omega_k             0.0
#define vluce_cm 2.99792458e+10 // cm/s
#define vluce_m  299792458      /* veocita' della luce in m/s */
#define vluce_km 299792.458     /* veocita' della luce in km/s */
#define HPLANCK		6.626075540e-27 // erg*s (h)
#define HPLANCK_TeV     4.135672e-27 // TeV*s (h)
#define HTPLANCK	1.0545726663e-27 // erg*s (htagl.=h/2pi)
#define q_esu		4.803206815e-10 // esu
#define erg_to_eV       6.24151e+11//
#define ev_to_erg       1.602176e-12//
#define Tev_to_erg       1.602176//
#define erg_to_TeV      0.6241512//
#define bn_to_cm2       1E-24//
#define SIGTH		6.652461618e-25 // cm^2 (sez. d'urto di thomson)
#define E_th_pp         1.2205E-3// TeV (th energy for pp Eq. 79 in Kelner et al. 2006, astro-ph.06066058v1, PHYSICAL REVIEW D 74, 034018 (2006))
#define Kpi            0.17//constant in Eq. 77 in Kelner et al. 2006, astro-ph.06066058v1, PHYSICAL REVIEW D 74, 034018 (2006))
#define Kpp_e          0.095//fraction of e- energy expressed as a fraction of Kpi
#define MEC2        8.187111e-07// erg  (me*c^2)
#define MPC2        1.5032764261e-3// erg  (mp*c^2)
#define MEC2_TeV    0.000000510998910// TeV  (me*c^2)
#define MPC2_TeV    0.000938272013// TeV  (mp*c^2)
#define MPIC2_TeV   0.0001349766//TeV  (mpi*c^2)
#define one_by_MEC2     1.221432045436937E6 // erg  (me*c^2)
#define me_g		9.109389754e-28 //massa di e- in  gr
#define mp_by_me 1836.15 // m_protne / m_elettrone
#define me 0.51/* massa a riposo di e- in MeV/c^2 */
#define m_sun 1.988992E33 //solar mass in gr
#define e_raggio 2.817940285e-13 /* raggio class. e- cm */
#define K_boltz 1.3806503e-16 /*costante di Boltzman espressa in erg/K */
#define sigma_steph_boltz 5.670400e-5 /* sigma Stephan_Boltzmann in erg*s^-1*cm^-2*K^-4*/
#define q 1.602+1e-19 /* carica dell 'elettrone in coulomb */
#define pi 3.1415926535897932384626433832795028841971693993751 /* pi */
#define four_pi  12.566370614359172953850573533118 // 4*pi
#define Deg_to_Rad 1.745329251994330e-2// deg_to_rad
#define four_by_three_pi 4.1887902047863909846168578443727 // (4/3)*pi
#define one_by_four_pi  0.079577471545947667884441881686257 // 1/(4*pi)
#define four_by_three 1.3333333333333333333333333333333333333333333 /* 4/3 */
#define three_by_four 0.75 /* 4/3 */
#define fake_nu_err 1.0
#define fake_flux_err 1.0
#define static_spec_arr_size 1000 /* num elementi dei vettori */
#define static_bess_table_size 1000 /* num elementi tabelle di Bessel */
#define Bessel_MAX 500.0
#define ELEMENTI_GAMMA_Contr_Comp 100
#define static_file_name_max_legth 512
#define LIM_LOSS_KN 1.0
#define min(a,b) (a<b) ? a:b;
#define max(a,b) (a>b) ? a:b;

/**
 * \file Blazar_SED.h
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief Header file principale
 *
 */

/************ ENV VARIBLE *************************/

/********************************     STRUTTURA BASE    ************************************/
struct spettro {
    int verbose;
    int BESSEL_TABLE_DONE;

    int CICCIO;

    char * SYSPATH;
    char STEM[256];
    char path[256];
    char DISTR[16];
    char disk_type[16];
    char MODE[16];
    char PARTICLE[16];
    int OUT_FILE;
    int START_FILE;
    int Num_file;
    int SSC, EC, TOT;
    int WRITE_TO_FILE;

    int do_Sync, do_SSC,do_IC,Sync_kernel;
    //int attesa_Sync_cooling, attesa_compton_cooling;

    unsigned long spec_array_size;

    int disk;
    double emiss_lim;
    //--- EMITTING SOURCE PARAMETERS
    double R; /* raggio blob sferica o raggio shell cilindrica */
    double B;
    double beam_obj;
    char BEAMING_EXPR[16];
    double theta; //viewing angle
    double BulkFactor;
    double beta_Gamma;
    double dist;
    double z_cosm;
    double Vol_sphere;
    double Surf_sphere;

    //----- Summed Spectra-----//
    //double nuF_nu_Sum_obs[static_spec_arr_size];


    //---- somma observer frame
    unsigned long nu_grid_size ;
    double nu_start_grid;
    double nu_stop_grid;

    double nu_grid[static_spec_arr_size];
    double nuFnu_sum_grid[static_spec_arr_size];

    double nuFnu_Sync_grid[static_spec_arr_size];
    double nuFnu_SSC_grid[static_spec_arr_size];
    double nuFnu_Disk_grid[static_spec_arr_size];
    double nuFnu_DT_grid[static_spec_arr_size];
    double nuFnu_Star_grid[static_spec_arr_size];
    double nuFnu_EC_CMB_grid[static_spec_arr_size];
    //double nuFnu_EC_CMB_stat_grid[static_spec_arr_size];
    double nuFnu_EC_BLR_grid[static_spec_arr_size];
    double nuFnu_EC_DT_grid[static_spec_arr_size];
    double nuFnu_EC_Disk_grid[static_spec_arr_size];
    double nuFnu_EC_Star_grid[static_spec_arr_size];

    //-----------Sync --------------//
    //--- CONST
    double C1_Sync_K53, C2_Sync_K53, C3_Sync_K53;
    double C1_Sync_K_AVE, C2_Sync_K_AVE;
    double COST_Sync_COOLING;

    double sin_psi; /*pitch angle*/
    double UB;
    double nu_B;

    //--- FREQ BOUNDARIES
    double nu_start_Sync;
    double nu_stop_Sync;
    double nu_start_Sync_obs;
    double nu_stop_Sync_obs;
    double nu_stop_Sync_ssc;
    unsigned long NU_INT_STOP_Sync_SSC;

    //--- FREQ/FLUX scalars
    double nu_peak_Sync_blob;
	double nuLnu_peak_Sync_blob;
	double nu_peak_Sync_src;
	double nuLnu_peak_Sync_src;
	double nu_peak_Sync_obs;
	double nuFnu_peak_Sync_obs;

	//--- FREQ/FLUX array
    double j_Sync[static_spec_arr_size];
    double alfa_Sync[static_spec_arr_size];
    double I_nu_Sync[static_spec_arr_size];
    double nu_Sync[static_spec_arr_size];
    double nu_Sync_obs[static_spec_arr_size];
    double n_Sync[static_spec_arr_size];
    double nuF_nu_Sync_obs[static_spec_arr_size];

    //--- Tabelle Bessel
    double F_Sync_x[static_bess_table_size];
    double F_ave_Sync_x[static_bess_table_size];
    double F_Sync_y[static_bess_table_size];
    double F_ave_Sync_y[static_bess_table_size];
    double log_F_Sync_x[static_bess_table_size];
    double log_F_Sync_y[static_bess_table_size];
    double log_F_ave_Sync_x[static_bess_table_size];
    double log_F_ave_Sync_y[static_bess_table_size];
    double t_Bessel_min, t_Bessel_max;
    double x_Bessel_min, x_Bessel_max;
    double x_ave_Bessel_min, x_ave_Bessel_max;
    double log_x_Bessel_min, log_x_Bessel_max;
    double log_x_ave_Bessel_min, log_x_ave_Bessel_max;

    //--------------------------------//

    //-----------pp-gamma-emission---//
    //--- CONST
    double NH_pp;

    //--- FREQ BOUNDARIES
    double nu_stop_pp_pred,nu_stop_pp;
    double nu_start_pp;
    unsigned long NU_INT_STOP_PP;

    //--- FREQ/FLUX array
    double j_pp[static_spec_arr_size];
    double nu_pp[static_spec_arr_size];
    double nuF_nu_pp_obs[static_spec_arr_size];

    //
    int set_pp_racc_gamma, set_pp_racc_elec;
    double pp_racc_gamma, pp_racc_elec;

    //--- FREQ/FLUX scalars
    double nu_peak_PP_blob;
    double nuLnu_peak_PP_blob;
    double nu_peak_PP_src;
    double nuLnu_peak_PP_src;
    double nu_peak_PP_obs;
    double nuFnu_peak_PP_obs;
    //--------------------------------//



    //-----------SSC-IC--------------//
    //--- CONST
    int ord_comp;
    double COST_IC_K1,COST_IC_COOLING ;

    //--- FREQ BOUNDARIES


    double nu_start_SSC;
    double nu_stop_SSC;
    double nu_start_SSC_obs;
    double nu_stop_SSC_obs;

    unsigned long NU_INT_STOP_COMPTON_SSC;

    //--- IC Kernel computation
    double Gamma;
    double nu; /*freq spettro sinc */
    double nu_1; /*freq spettro comp */
    double nu_compton_0; /* freq campo fot Sync per spettro IC */
    double * nu_seed;
    double * n_seed;

    //--- FREQ/FLUX array
    double q_comp[static_spec_arr_size];
    double j_comp[static_spec_arr_size];
    double j_EC[static_spec_arr_size];
    double nu_SSC[static_spec_arr_size];
    double nu_SSC_obs[static_spec_arr_size];
    double nuF_nu_SSC_obs[static_spec_arr_size];

    //--- FREQ/FLUX scalars
    double nu_peak_SSC_blob;
    double nuLnu_peak_SSC_blob;
    double nu_peak_SSC_src;
    double nuLnu_peak_SSC_src;
    double nu_peak_SSC_obs;
    double nuFnu_peak_SSC_obs;


    //-----------EC--------------//
    

    //Const
    int do_EC_Disk,do_EC_BLR,do_EC_DT,do_EC_Star,do_EC_CMB,EC_stat;
    int do_Disk,do_DT;
    double nu_planck_min_factor;
    double nu_planck_max_factor;
    double mono_planck_min_factor;
    double mono_planck_max_factor;
    //double EC_field_interp_factor;
    unsigned long theta_n_int;
    unsigned long l_n_int;

    double nu_blob_RF;
    double nu_disk_RF;
    double L_nu_disk_RF;

    //dist BLOB/DISK or STAR
    double R_H;

    //--- STAR
    //-PARAMTERS
    double L_Star;
    double T_Star_max;
    double R_Star;
    double Star_psi_1, Star_psi_2;
    //
    double Star_surface;
    double Star_mu_1, Star_mu_2;
    //-FREQ BOUNDARIES
	double nu_start_EC_Star;
	double nu_stop_EC_Star;
	double nu_start_EC_Star_obs;
	double nu_stop_EC_Star_obs;
	double nu_start_Star_obs;
	double nu_stop_Star_obs;
	double nu_start_Star;
	double nu_stop_Star;
    double nu_start_Star_DRF;
    double nu_stop_Star_DRF;

    unsigned long NU_INT_MAX_Star;
	unsigned long NU_INT_STOP_EC_Star;
	//-FREQ/FLUX arrays
	double I_nu_Star[static_spec_arr_size];
	double J_nu_Star_disk_RF[static_spec_arr_size];
	double I_nu_Star_disk_RF[static_spec_arr_size];
	double nu_Star[static_spec_arr_size];
	double nu_Star_obs[static_spec_arr_size];
	double nu_Star_disk_RF[static_spec_arr_size];
	double nuF_nu_Star_obs[static_spec_arr_size];
	double nu_EC_Star[static_spec_arr_size];
	double nu_EC_Star_obs[static_spec_arr_size];
	double nuF_nu_EC_Star_obs[static_spec_arr_size];
	double n_Star[static_spec_arr_size];
    double n_Star_DRF[static_spec_arr_size];

    //--- CMB
	//-PARAMTERS
	double T_CMB_0;
	double CMB_mu_1,CMB_mu_2;
	//-FREQ BOUNDARIES
	double nu_start_CMB, nu_stop_CMB;
	double nu_start_EC_CMB;
	double nu_stop_EC_CMB;
	double nu_start_EC_CMB_obs;
	double nu_stop_EC_CMB_obs;
    double nu_start_CMB_DRF;
    double nu_stop_CMB_DRF;
    unsigned long NU_INT_MAX_CMB,NU_INT_STOP_EC_CMB;
    //-FREQ/FLUX arrays
	double I_nu_CMB[static_spec_arr_size];
	double I_nu_CMB_disk_RF[static_spec_arr_size];
	double nu_CMB[static_spec_arr_size];
	double nu_CMB_disk_RF[static_spec_arr_size];

	double nu_EC_CMB[static_spec_arr_size];
	double nu_EC_CMB_obs[static_spec_arr_size];
	double nuF_nu_EC_CMB_obs[static_spec_arr_size];
    double n_CMB[static_spec_arr_size];
    double n_CMB_DRF[static_spec_arr_size];

    //TODO REMOVE UNSUED FROM THE CODE
    //--- CMB stat
    //-FREQ BOUNDARIES
    //double nu_start_CMB_stat, nu_stop_CMB_stat;
    //double nu_start_EC_CMB_stat;
    //double nu_stop_EC_CMB_stat;
    //double nu_start_EC_CMB_stat_obs;
    //double nu_stop_EC_CMB_stat_obs;
    //unsigned long NU_INT_MAX_CMB_stat,NU_INT_STOP_EC_CMB_stat;
    //-FREQ/FLUX arrays
    //double I_nu_CMB_stat[static_spec_arr_size];
    //double I_nu_CMB_disk_RF_stat[static_spec_arr_size];
    //double nu_CMB_stat[static_spec_arr_size];
    //double nu_CMB_disk_RF_stat[static_spec_arr_size];

    //double nu_EC_CMB_stat[static_spec_arr_size];
    //double nu_EC_CMB_stat_obs[static_spec_arr_size];
    //double nuF_nu_EC_CMB_stat_obs[static_spec_arr_size];
    //double n_CMB_stat[static_spec_arr_size];


    //--- DISK
    //-PARAMTERS
    double L_Disk;
    double L_Disk_radiative;
    double T_Disk;
    double accr_eff;
    double M_BH,accr_rate,L_Edd,accr_Edd;
    double R_inner_Sw, R_ext_Sw;
    double R_Disk_interp;
    //
    double R_Sw;
    double T_disk_max_4;
    double R_inner, R_ext;
    double Cost_disk_Mulit_BB;
    //double Cost_Norm_disk_Mulit_BB;
    double Disk_surface;
    double Disk_mu_1,Disk_mu_2;
    double Disk_geom_factor;
    //-FREQ BOUNDARIES
    double nu_disk_Multi_BB;
    double nu_start_EC_Disk;
    double nu_stop_EC_Disk;
    double nu_start_EC_Disk_obs;
    double nu_stop_EC_Disk_obs;
    double nu_start_Disk;
    double nu_stop_Disk;
    double nu_start_Disk_obs;
    double nu_stop_Disk_obs;
    double nu_start_Disk_DRF;
    double nu_stop_Disk_DRF;
    unsigned long NU_INT_MAX_Disk;
    unsigned long NU_INT_STOP_EC_Disk;

    //-FREQ/FLUX arrays
    double L_nu_Disk_disk_RF[static_spec_arr_size];
    double I_nu_Disk[static_spec_arr_size];
    //double J_nu_Disk_disk_RF[static_spec_arr_size];
    double I_nu_Disk_disk_RF[static_spec_arr_size];
    double nu_Disk[static_spec_arr_size];
    double nu_Disk_obs[static_spec_arr_size];
    double nu_Disk_disk_RF[static_spec_arr_size];
    double nuF_nu_Disk_obs[static_spec_arr_size];
    double nu_EC_Disk[static_spec_arr_size];
    double nu_EC_Disk_obs[static_spec_arr_size];
    double nuF_nu_EC_Disk_obs[static_spec_arr_size];
    double n_Disk[static_spec_arr_size];
    double n_Disk_DRF[static_spec_arr_size];

    //--- FREQ/FLUX scalars

    //--- BLR
    //-PARAMETERS
    double tau_BLR;
    double R_BLR_in;
    double R_BLR_out;
    double R_BLR_interp_val;
    double R_BLR_interp_start;

    //
    double BLR_Volume;
    double n0_BLR;
    double mu_j;
    double BLR_mu_1, BLR_mu_2;
    //double BLR_mu_r_J_1,BLR_mu_r_J_2;
    //double BLR_mu_r_in;
    double BLR_inner_Surface;
    double Delta_R_BLR;
    //double BLR_geom_factor;
    //-FREQ BOUNDARIES
    double nu_stop_BLR;
    double nu_start_BLR;
    double nu_stop_BLR_disk_RF;
    double nu_start_BLR_disk_RF;
    double nu_stop_EC_BLR;
    double nu_start_EC_BLR;
    double nu_stop_EC_BLR_obs;
    double nu_start_EC_BLR_obs;
    unsigned long NU_INT_MAX_BLR;
    unsigned long NU_INT_STOP_EC_BLR;
    
    //-FREQ/FLUX arrays
    double I_nu_BLR[static_spec_arr_size];
    double Lnu_BLR_disk_RF[static_spec_arr_size];
    double nu_BLR[static_spec_arr_size];
    double I_nu_BLR_disk_RF[static_spec_arr_size];
    double nuF_nu_EC_BLR_obs[static_spec_arr_size];
    double nu_EC_BLR[static_spec_arr_size];
    double nu_EC_BLR_obs[static_spec_arr_size];
    double nu_BLR_disk_RF[static_spec_arr_size];
    double n_BLR[static_spec_arr_size];
    double n_BLR_DRF[static_spec_arr_size];

    //--- DT
    //-PARAMETERS
    double T_DT;
    double tau_DT;
    double L_DT;
    double R_DT;
    double R_DT_interp_val ;
    double R_DT_interp_start;

    //
    double DT_mu_1,DT_mu_2;
    double DT_mu_r_J_1;
    double DT_mu_r_J_2;
    double DT_Volume;
    //-FREQ BOUNDARIES
    double nu_stop_DT;
    double nu_start_DT;
    double nu_start_DT_DRF;
    double nu_stop_DT_DRF;
    double nu_start_DT_obs;
    double nu_stop_DT_obs;
    double nu_stop_EC_DT;
    double nu_start_EC_DT;
    double nu_stop_EC_DT_obs;
    double nu_start_EC_DT_obs;
    unsigned long NU_INT_MAX_DT;
    unsigned long NU_INT_STOP_EC_DT;

    //-FREQ/FLUX arrays
    double I_nu_DT[static_spec_arr_size];
    double I_nu_DT_disk_RF[static_spec_arr_size];
    double nu_DT_obs[static_spec_arr_size];
    double nu_DT[static_spec_arr_size];
    double nu_DT_disk_RF[static_spec_arr_size];
    double nuF_nu_EC_DT_obs[static_spec_arr_size];
    double n_DT[static_spec_arr_size];
    double n_DT_DRF[static_spec_arr_size];
    double L_nu_DT_disk_RF[static_spec_arr_size];
    double nuF_nu_DT_obs[static_spec_arr_size];
    double nu_EC_DT[static_spec_arr_size];
    double nu_EC_DT_obs[static_spec_arr_size];


    //--- FREQ/FLUX scalars
    double nu_peak_EC_Disk_blob;
	double nuLnu_peak_EC_Disk_blob;
	double nu_peak_EC_Disk_src;
	double nuLnu_peak_EC_Disk_src;
	double nu_peak_EC_Disk_obs;
	double nuFnu_peak_EC_Disk_obs;
    double nu_peak_EC_BLR_blob;
	double nuLnu_peak_EC_BLR_blob;
	double nu_peak_EC_BLR_src;
	double nuLnu_peak_EC_BLR_src;
	double nu_peak_EC_BLR_obs;
	double nuFnu_peak_EC_BLR_obs;
	double nu_peak_EC_DT_blob;
	double nuLnu_peak_EC_DT_blob;
	double nu_peak_EC_DT_src;
	double nuLnu_peak_EC_DT_src;
	double nu_peak_EC_DT_obs;
	double nuFnu_peak_EC_DT_obs;

    double beaming_EC;

    


    //----------- INTEGRATION MESH--------------//
    unsigned long nu_seed_size;
    unsigned long nu_IC_size;

    //unsigned long mesh_intComp;
    //unsigned long mesh_intComp1;


    //----------- PARTICLE DISTRIBUTION --------------//
    int Norm_distr;
    //double Norm_distr_L_e_Sync;
    int Distr_e_done;
    int Distr_p_done;
    int Distr_e_pp_done;

    int TIPO_DISTR;
    double *Ne_custom;
    double *gamma_e_custom;
    unsigned long gamma_custom_grid_size;
    double *gam;
    double *Ne;
    double *Ne_IC;
    double *Ne_stat;
    double *Np;
    unsigned long gamma_grid_size;
    double * griglia_gamma_Ne_log;
    double * griglia_gamma_Ne_log_IC;
    double * griglia_gamma_Ne_log_stat;
    double * griglia_gamma_Np_log;
   
    //double *griglia_gamma_log_IC;
    //double *N_IC;

    unsigned long i_griglia_gamma;
    double N_tot_e_Sferic;
    double N_tot_p_Sferic;
    double N,N_e_pp;
    double N_0,N_0p,N_0e; /* costante di normalizzazione per distrib elettr staz */
    double gmin;
    double gmax;
    int grid_bounded_to_gamma;
    //unsigned long pt_griglia_max;
    double gmin_griglia;
    double gmax_griglia;
    double U_e, E_tot_e;
    double U_p, E_tot_p;
    double Gamma_p2; //gamma peak of N(gamma)*gamma^2
    double Gamma_p3; //gamma peak of N(gamma)*gamma^3
    double Np2, Np3; //peak of N(gamma)*gamma^2 and  N(gamma)*gamma^3

    //PL,PLC,BKN
    double p, p_1;
    double gamma_break;
    double gamma_cut;

    //LP+LPPL
    double s;
    double r;
    double gamma0_log_parab;
    //LPEP
    double gammap_log_parab;
    //LPEP PILEUP
    double gamma_inj;
    //Spit
    double spit_index,spit_temp,spit_gamma_th;

    //LPPL Pile-Up
    double gamma_pile_up,gamma_pile_up_cut,alpha_pile_up;
    double ratio_pile_up;

};

//===================================================================================

struct jet_energetic{
    double U_e,U_p,U_B;
    double U_Synch, U_Synch_DRF;
    double U_Disk, U_BLR, U_DT, U_CMB;
    double U_Disk_DRF, U_BLR_DRF, U_DT_DRF, U_CMB_DRF;
    double L_Sync_rf, L_SSC_rf, L_EC_Disk_rf,L_EC_BLR_rf, L_EC_DT_rf,L_EC_CMB_rf, L_PP_rf;
    double jet_L_Sync,jet_L_SSC, jet_L_EC_Disk, jet_L_EC_BLR, jet_L_EC_DT,jet_L_EC_CMB,jet_L_PP;
    double jet_L_rad,jet_L_kin, jet_L_tot, jet_L_e, jet_L_B, jet_L_p;
};



//=========================TEMPORAL EVOLUTION ================================

struct temp_ev{
    char STEM[256];
    char path[256];

	int do_EC_cooling_Disk;
	int do_EC_cooling_BLR;
	int do_EC_cooling_Star;
	int do_EC_cooling_CMB;
	int do_EC_cooling_DT;


	int do_Sync_cooling;
	int do_SSC_cooling;
    int do_Compton_cooling;



	double t_unit; //unit time in light crossing time
	//double t_acc;
	double L_inj;
	double *T_esc;
	double Diff_Coeff;
	double Diff_coeff_CD;
	double Diff_coeff_CA;
	double Acc_Coeff;
	double Diff_Index;
	double Acc_Index;
	double T_esc_Coeff;
	double Esc_Index;
	double Lambda_max_Turb;
	double Lambda_choer_Turb_factor;
	double Gamma_Max_Turb_L_max;
	double Gamma_Max_Turb_L_coher;
	double gamma_eq_t_D,gamma_eq_t_A;
	//double alfa_inj;
	double TStart_Inj;
	double TStop_Inj;
	double TStart_Acc;
	double TStop_Acc;
	double Inj_temp_slope;
	unsigned long NUM_SET;
	unsigned long T_SIZE;
	double duration,t_D0,t_DA0,t_A0;

};

//void griglia (double *, double,double,int);
double Cfp(double,struct temp_ev *);
double Bfp(double x, struct temp_ev *pt, struct spettro *pt_spec);
double f_Dp(double gamma,struct temp_ev *pt);
double f_DAp(double gamma,struct temp_ev *pt);
double f_A(double gamma,struct temp_ev *pt);


double f_Acc(double x, struct temp_ev *);
double f_Tesc(double,struct temp_ev *);
double Inj_temp_prof (double t,struct temp_ev *);
//double f_Inj(double x,struct temp_ev *);
void Wm (double, double* ,double *);
double Cooling(double,struct temp_ev *,struct spettro *);
int solve_sys1(double VX1[],double VX2[],double VX3[],double SX[],double u[],unsigned long size);
//===================================================================================






// File Interface
// void FileInput(int argc, char **argv, struct spettro *pt_spec,struct temp_ev *pt_temporal);


//=======================================================================================
/********************************     PyInterface    ************************************/
// PyInterface
struct spettro MakeBlob();
void MakeNe(struct spettro *pt_base);
struct temp_ev MakeTempEv();
void Init(struct spettro *pt);
void InitNe(struct spettro *pt);
//void build_photons(struct spettro *pt_base);
//void alloc_photons(double ** pt,int size);
void set_seed_freq_start(struct spettro *pt_base);
void Run_SED(struct spettro *pt_base);
void Run_temp_evolution(struct spettro *pt_spec, struct temp_ev *pt_ev);
double get_spectral_array(double * arr, struct spettro *pt, unsigned long id);
double get_elec_array(double * arr, struct spettro *pt, unsigned long id);
double set_elec_custom_array(double * arr, struct spettro *pt,double val, unsigned long id);
//===================================================================================







//===================================================================================
//ENERGETIC
double GetU_e(struct spettro *pt);
void EvalU_e(struct spettro *pt);
double GetE_tot(struct spettro *pt);
void CoolingRates(struct spettro * pt, struct temp_ev *pt_ev);
//===================================================================================



//===================================================================================
/************************************ FUNZIONI Distr N *****************************/
//void Genera_griglia_gamma_e_log(struct spettro *pt, double *griglia);
void alloc_N_distr(double ** pt,int size);

void Fill_N(struct spettro *pt, double *griglia_gamma_N_log, double *N);
void build_Ne_custom(struct spettro *pt,  unsigned int size);
void build_Ne(struct spettro *pt);
void Fill_Ne_IC(struct spettro *pt, double gmin, int stat_frame);
void Genera_Np_Ne_pp(struct spettro *pt);
double N_distr(struct spettro *, double Gamma);
double N_distr_U_e(struct spettro *, double Gamma);
double N_distr_U_p(struct spettro *, double Gamma);
double N_distr_integranda(struct spettro *, double Gamma);
void Scrivi_N_file(struct spettro *pt, char *name, double *g, double *N);
void FindNe_NpGp(struct spettro *pt);
double N_distr_interp(unsigned long size, double Gamma, double *griglia_gamma, double *N);
double Find_gmax(struct spettro *pt, double *g, double *N);
double pl_func(double Gamma,double p);
double plc_func(double Gamma,double gamma_cut,double p);
double bkn_func(double Gamma,double gamma_break,double p, double p_1);
//double pile_up_ratio(double Gamma,double sigma,double gamma_inj,double gamma_eq);
double bkn_pile_up_func(double Gamma,double gamma_inj, double p, double p_1,double gamma_eq, double gamma_cut ,double alpha);
double lp_func(double Gamma,double gamma0,double r, double s);
double lp_ep_func(double Gamma,double gamma_p,double r);
double lppl_func(double Gamma,double gamma0, double r, double s);
double pile_up_func(double Gamma, double gamma_pile_up_cut, double alpha_pile_up );
double lppl_pile_up_func(double Gamma,double gamma0, double gamma_inj, double r, double s,double gamma_eq,double gamma_cut, double alpha);
double spit_func(double Gamma,double gamma_th,double temp, double index);
//===================================================================================





//===================================================================================
/************************************ FUNZIONI VARIE *******************************/
void messaggio_errore();
//void manpage();
void flux_header(FILE *fp);
void flux_DISK_header(FILE *fp);
void distr_e_header(FILE *fp);
void somma_header(FILE *fp);
//===================================================================================



//===================================================================================
/************************************ SUMMED SPECTRA *******************************/
void interpola_somma(struct spettro *pt_j, double nu, unsigned long i);
void spettro_somma_Sync_ic(int num_file, struct spettro *);
//===================================================================================




//===================================================================================
/************************************ ENERGETIC *******************************/
double Power_Sync_Electron(struct spettro *pt);
double Uph_Sync(struct spettro *pt);
double Power_Sync_Electron_Integ(struct spettro *pt_N, double Gamma);

double Lum_Sync_at_nu (struct spettro *pt, double nu);
double	Lum_SSC_at_nu (struct spettro *pt , double nu_1);

void FindEpSp(double * nu_blob, double * nuFnu_obs, unsigned long NU_INT_MAX, struct spettro * pt,
        double * nu_peak_obs,
        double * nu_peak_src,
        double * nu_peak_blob,
        double * nuFnu_peak_obs,
        double * nuLnu_peak_src,
        double * nuLnu_peak_blob);


struct jet_energetic EnergeticOutput(struct spettro *pt, int write_file);

double I_nu_to_Uph(double * nu, double * I_nu, unsigned long NU_INT_STOP);
double PowerPhotons_blob_rest_frame(struct spettro *pt, double * nu, double * nu_Fnu, unsigned long NU_INT_STOP);
double PowerPhotons_disk_rest_frame(struct spettro *pt, double *nu, double *nu_Fnu, unsigned long NU_INT_STOP);

//===================================================================================





//======================================================================================
/************************************ GEOMETRY/REL.SPEC. *******************************/
double V_sphere(double R);
double S_sphere(double R);

double get_beaming(double BulkFactor, double theta);
void   SetBeaming(struct spettro *pt);
void SetDistr(struct spettro *pt);

double eval_beta_gamma(double  gamma);
double Larmor_radius(double gamma,double B, double sin_alpha);
double Larmor_radius_to_gamma(double Larmor_radius,double B, double sin_alpha);
//===================================================================================



//===========================================================================================
/**********************              FUNZIONI COSMOLOGICHE                *******************/
double Distanza_Lum_analyt(double);
double eta_Distanza_Lum_analyt(double, double);
double distanza_z( double);
double dist_lum_cm(double );
//===========================================================================================



//===========================================================================================
/********************************* FREQ/FLUX TRANSFORMATIONS JET/OBSERVER *******************/
double L_nu_Disk_to_F_nu(double L_nu, double z, double dist);
double L_nu_blob_to_F_nu(double L_nu_Sync, double beam_obj, double z, double dist);
double L_nu_src_to_F_nu(double L_nu_Sync, double beam_obj, double z, double dist);
double nuFnu_obs_to_nuLnu_src(double nuFnu_obs, double beam_obj, double z, double dist);
double nuFnu_obs_to_nuLnu_blob(double nuFnu_obs, double beam_obj, double z, double dist);
double I_nu_to_L_nu_blob(double I_nu, double R);
double I_nu_to_L_nu_src(double I_nu, double R, double beam_obj);
double j_nu_to_L_nu_blob(double j_nu, double R);
double j_nu_to_L_nu_src(double j_nu, double R, double beam_obj);
double nu_blob_to_nu_obs(double nu, double delta, double z);
double nu_blob_to_nu_src(double nu, double delta, double z);
double nu_disk_to_nu_obs_disk(double nu, double z);
double nu_obs_disk_to_nu_disk(double nu, double z);
double nu_obs_to_nu_blob(double nu, double delta, double z);
double nu_obs_to_nu_src(double nu, double z);
//===========================================================================================




//===========================================================================================
/****************************** FUNZIONI SPETTRO SINCROTRONE ********************************/
void spettro_sincrotrone(int Num_file, struct spettro *);
double F_K_53(struct spettro *, double x);
double F_K_ave(struct spettro *, double x);
double F_int_fix(struct spettro *,  unsigned long  ID);
double F_int_ave(struct spettro *, unsigned long  ID);
double Sync_self_abs_int(struct spettro *, unsigned long  ID);
double j_nu_Sync(struct spettro *);
double solve_S_nu_Sync(struct spettro * pt, unsigned long  NU_INT);
double alfa_nu_Sync(struct spettro *);
double integrale_Sync(double (*pf) (struct spettro *, unsigned long  ID), struct spettro * pt );
double Sync_tcool(struct spettro * , double g);
double Sync_cool(struct spettro * , double g);
//===========================================================================================





//===========================================================================================
/****************** FUNZION pp-gamma-emssion ************************************************/
void spettro_pp(int Num_file, struct spettro *pt);
double F_gamma(double x, double Ep_TeV);
double F_electrons(double x, double Ep_TeV);
double sigma_pp_inel(double Ep_TeV);
double rate_gamma_pp(struct spettro *pt);
double rate_electrons_pp(struct spettro *pt, double Gamma);
double pp_gamma_kernel(double gamma_p, double nu_pp,struct spettro * pt);
double pp_gamma_kernel_delta(struct spettro *pt, double Epi_TeV);
double pp_electron_kernel_delta(double Ee_TeV, struct spettro *pt);
double pp_electrons_kernel(double gamma_p, double Gamma_e, struct spettro *pt);
double integrale_pp_gamma_rate(double (*pf_pp_gamma_kernel) (double gamma_p, double nu_pp,struct spettro *pt), struct spettro * pt, unsigned long i_start);
double integrale_pp_electrons_rate(double (*pf_pp_gamma_kernel) (double gamma_p, double Gamma_e,struct spettro *pt), struct spettro * pt, unsigned long i_start);
//===========================================================================================




//========================================================================================
/************************************ PHOTON GRID FUNCTIONS *****************************/
void build_log_grid(double nu_start, double nu_stop, unsigned long SIZE, double * nu_grid);
unsigned long x_to_grid_index(double * nu_grid, double nu, unsigned long SIZE);
//===========================================================================================





//=================================================================================
/************************************ FUNCTIONS  IC  *****************************/
double f_compton_K1(struct spettro *, double Gamma);
double rate_compton_GR(struct spettro *);
double integrale_IC(double (*pf) (struct spettro *, double x), struct spettro * pt, double a, double b,int stat_frame);
double integrale_IC_cooling(struct spettro * pt, double a, double b, double gamma);
double compton_cooling(struct spettro *pt_spec, struct temp_ev *pt_ev, double gamma);
double f_compton_cooling(double b);
double I_nu_to_n(double I_nu, double nu);
//==================================================================================




//=============================================================================
/********************* FUNZIONI SPETTRO COMPTON *******************************/
void spettro_compton(int num_file, struct spettro *);
//=============================================================================



//===================================================================================
/********************* FUNZIONI EXTERNAL COMPOTON************************************/
void spectra_External_Fields(int Num_file, struct spettro *pt_d);
void spettro_EC(int num_file, struct spettro *);

/***  DISK PLANCK FUNCTIONS  *********/
double eval_T_disk(struct spettro *pt, double R);
double eval_nu_peak_planck(double T);
double f_planck_Multi_T(struct spettro *pt, double R, double nu);
double f_planck_Multi_T_norm(struct spettro *pt, double R,double nu);
double integrand_f_planck_Multi_T(struct spettro *pt,double R);
double f_planck(double temperatura, double nu);
double f_planck_norm(double temperatura, double nu);

/***  DISK ACCRETION POWER FUNCTIONS  **********/
double eval_R_Sw(double M_BH);
double eval_M_BH(double R_S);
double eval_accr_rate(double L_disk,double accr_eff);
double eval_L_Edd(double M_BH);
double eval_accr_Edd(double L_Edd, double accr_eff);

/***  EXTERNAL COMPOTON DISK/JET TRANSFORMATIONS  *******/
double eval_nu_min_blob_RF(struct spettro *pt, double mu1, double mu2, double nu_disk_RF );
double eval_nu_max_blob_RF(struct spettro *pt, double mu1, double mu2, double nu_disk_RF );
double I_nu_disk_RF_to_blob_RF(double I_nu_diks_RF, double nu_disk_RF, double nu_blob_RF, double beta, double Gamma);
double nu_blob_RF_to_nu_disk_RF(double nu_blob_RF, double Gamma, double beta, double mu_disk_RF);


/***  SPECTRAL/GEOMETRIC FUNCTIONS EC Star ****/
void Build_I_nu_Star(struct spettro *pt_d);
double eval_I_nu_Star_disk_RF(struct spettro *pt,double nu_Star_disk_RF);
double eval_J_nu_Star_disk_RF(struct spettro *pt, double I_nu_Star_disk_RF);
double eval_I_nu_Star_blob_RF(struct spettro *pt, double nu_blob_RF);
double integrand_I_nu_Star_blob_RF(struct spettro *pt, double mu);
double eval_Star_L_nu(struct spettro *pt, double nu_Star_disk_RF);
void set_Star_geometry(struct spettro *pt);

/***  SPECTRAL/GEOMETRIC FUNCTIONS EC DISK ****/
void set_Disk(struct spettro *pt);
void set_Disk_angles(struct spettro *pt);
void Build_I_nu_Disk(struct spettro *pt_d);
double Disk_Spectrum(struct spettro *pt, double nu_Disk_disk_RF);

//double eval_I_nu_Disk_disk_RF(struct spettro *pt,double Lnu_Disk_disk_RF);
//double eval_J_nu_Disk_disk_RF(struct spettro *pt, double nu_Disk_disk_RF, double I_nu_Disk_disk_RF);
//double eval_J_nu_Disk_disk_RF_single_BB(struct spettro *pt, double I_nu_Disk_disk_RF);
//double eval_J_nu_Disk_disk_RF_multi_BB(struct spettro *pt,  double  mu);

double eval_I_nu_Disk_disk_RF(struct spettro *pt, double nu_disk_RF);
double eval_I_nu_Disk_blob_RF(struct spettro *pt, double nu_disk_RF);
double eval_I_nu_theta_Disk(struct spettro *pt, double mu);
double integrand_I_nu_Disk_blob_RF(struct spettro *pt, double mu);
double integrand_I_nu_Disk_disk_RF(struct spettro *pt, double mu );
double eval_Disk_L_nu(struct spettro *pt, double nu_Disk_disk_RF);
void set_Disk_geometry(struct spettro *pt);
double eval_nu_peak_Disk(double T);

/***  SPECTRAL/GEOMETRIC FUNCTIONS EC CMB ****/
void Build_I_nu_CMB(struct spettro *pt_d);
//void Build_I_nu_CMB_stat(struct spettro *pt_d);
double eval_I_nu_CMB_disk_RF(double T_CMB,double nu_CMB_disk_RF);
double eval_I_nu_CMB_blob_RF(struct spettro *pt, double nu_blob_RF);
double integrand_I_nu_CMB_blob_RF(struct spettro *pt, double mu);
double eval_T_CMB_z(double z, double T_CMB_0);


/*** SPECTRAL/GEOMETRIC FUNCTIONS EC BLR  ****/
void Build_I_nu_BLR(struct spettro *pt);
double j_nu_BLR_integrand(struct spettro *pt, double l);
double eval_I_nu_theta_BLR(struct spettro *pt, double mu);
double integrand_I_nu_BLR_disk_RF(struct spettro *, double theta);
double eval_I_nu_BLR_disk_RF(struct spettro *pt);
double integrand_I_nu_BLR_blob_RF(struct spettro *pt, double theta);
double eval_I_nu_BLR_blob_RF(struct spettro *pt);
double eval_theta_max_BLR(struct spettro *pt);
void eval_l_values_BLR(struct spettro *pt, double mu, double l[]);
double eval_Lnu_BLR_disk_RF(struct spettro *pt, double nu_disk_RF);
//double integrand_J_nu_BLR_disk_RF(struct spettro *pt, double mu);
//double integrand_I_nu_BLR_blob_RF(struct spettro *pt, double mu);
void set_BLR_geometry(struct spettro *pt);
//double l_BLR(struct spettro *pt, double mu);

/***  FUNCTIONS Seed Photons EC DT  ***/
void Build_I_nu_DT(struct spettro *pt);
//double j_nu_DT_integrand(struct spettro *pt, double l);
double eval_I_nu_theta_DT(struct spettro *pt, double mu);
double integrand_I_nu_DT_disk_RF(struct spettro *, double theta);
double eval_I_nu_DT_disk_RF(struct spettro *pt);
double eval_I_nu_DT_blob_RF(struct spettro *pt);
double integrand_I_nu_DT_blob_RF(struct spettro *pt, double theta);
double eval_DT_L_nu(struct spettro *pt, double DT_disk_RF);
double eval_theta_max_DT(struct spettro *pt);
double eval_l_DT(struct spettro *pt, double mu);
double eval_circle_secant(double z,double R,double mu);
//===========================================================================================





//===========================================================================================
/*********************        FUNZIONI MATEMATICHE             *****************************/
double st_gamma(double x);
void beschb(double x, double *gam1, double *gam2, double *gampl,
        double *gammi);
double chebev(float a, float b, float c[], int m, float x);
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);
double bessel_K_53(struct spettro *, double x);
double bessel_K_pitch_ave(struct spettro *pt, double x);
void tabella_Bessel(struct spettro *);
double derivata(double (*pf) (struct spettro *, double x), struct spettro *pt, double x);
double integrale_trap_log_struct(double (*pf) (struct spettro *, double x),
        struct spettro * pt, double a, double b, unsigned long intervalli);
double integrale_simp_struct(double (*pf) (struct spettro *, double x),
        struct spettro * pt, double a, double b, unsigned long intervalli);
double integrale_simp(double (*pf) ( double x), double a, double b, unsigned long n_intervalli);
double trapzd_array_linear_grid(double *x, double *y, unsigned long SIZE);
double trapzd_array_arbritary_grid(double *x, double *y, unsigned long SIZE);
double test_int(struct spettro *, double x);
double test_int1(double x);
double log_lin_interp(double nu,  double * nu_grid, double nu_min, double nu_max, double *  flux_grid , unsigned long SIZE, double emiss_lim);
double log_log_interp(double log_x,  double * log_x_grid, double log_x_min, double log_x_max, double *  log_y_grid , unsigned long SIZE, double emiss_lim);

//===========================================================================================







//===========================================================================================
/**********************         FUNZIONI DI TEST                     ************************/
double test_lunghezza_vettore(double mesh, double a, double b, int Max_elem);
double test_solid_anlge(struct spettro *pt, double mu);

#endif




