

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
#define E_th_pp         0.2797E-3// Th energy pi production from pp (Kafexhiu, E., Aharonian, F., Taylor, A. M., & Vila, G. S. 2014, Physical Review D, 90, 123014, doi: 10.1103/PhysRevD.90.123014)
#define Kpi            0.17//constant in Eq. 77 in Kelner et al. 2006, astro-ph.06066058v1, PHYSICAL REVIEW D 74, 034018 (2006))
#define MEC2        8.187111e-07// erg  (me*c^2)
#define MPC2        1.5032764261e-3// erg  (mp*c^2)
#define MEMUC2_TeV    0.00010565// TeV  (me*c^2)
#define MEC2_MeV    0.510998910// TeV  (me*c^2)
#define MEC2_TeV    0.000000510998910// TeV  (me*c^2)
#define MPC2_TeV    0.000938272013// TeV  (mp*c^2)
#define MPI0C2_TeV   0.0001349766//TeV  (mpi*c^2)
#define MPICC2_TeV   0.00013957018//TeV  (mpi*c^2)
#define one_by_MEC2     1.221432045436937E6 // erg  (me*c^2)
#define brem_ee_th  2// MeV threshold for rel/non rel Bremsstrahlung
#define me_g		9.109389754e-28 //massa di e- in  gr
#define mp_by_me 1836.15 // m_protne / m_elettrone
#define me 0.51/* massa a riposo di e- in MeV/c^2 */
#define m_sun 1.988992E33 //solar mass in gr
#define one_by_alpha 137.035999206 //Fine-structure constant
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
#define static_spec_arr_grid_size 10000
#define static_ev_arr_grid_size 1000
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
// struct spectrum {
//     double nu_min, nu_max;
//     unsigned int NU_INT_MAX;
//     double j_nu[static_spec_arr_size];
//     double I_nu[static_spec_arr_size];
//     double n_nu[static_spec_arr_size];
//     double nuFnu_obs[static_spec_arr_size];
//     double nu[static_spec_arr_size];
//     double nu_obs[static_spec_arr_size];
//     double nuFnu_grid[static_spec_arr_grid_size];
// };

// struct spectrum_external{
//     double nu_min, nu_max;
//     unsigned int NU_INT_MAX;
//     double j_nu[static_spec_arr_size];
//     double I_nu[static_spec_arr_size];
//     double I_nu_DRF[static_spec_arr_size];
//     double n_nu[static_spec_arr_size];
//     double n_nu_DRF[static_spec_arr_size];
//     double nuFnu_obs[static_spec_arr_size];
//     double nu[static_spec_arr_size];
//     double nu_obs[static_spec_arr_size];
//     double nuFnu_grid[static_spec_arr_grid_size];
// };


struct blob {
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

    unsigned int spec_array_size;

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
    unsigned int nu_grid_size ;
    double nu_start_grid;
    double nu_stop_grid;

    double nu_grid[static_spec_arr_grid_size];
    double nuFnu_sum_grid[static_spec_arr_grid_size];

    double nuFnu_Sync_grid[static_spec_arr_grid_size];
    double nuFnu_SSC_grid[static_spec_arr_grid_size];
    double nuFnu_Disk_grid[static_spec_arr_grid_size];
    double nuFnu_DT_grid[static_spec_arr_grid_size];
    double nuFnu_Star_grid[static_spec_arr_grid_size];
    double nuFnu_EC_CMB_grid[static_spec_arr_grid_size];
    double nuFnu_pp_gamma_grid[static_spec_arr_grid_size];
    double nuFnu_pp_neutrino_tot_grid[static_spec_arr_grid_size];
    double nuFnu_pp_neutrino_mu_grid[static_spec_arr_grid_size];
    double nuFnu_pp_neutrino_e_grid[static_spec_arr_grid_size];
    double nuFnu_bremss_ep_grid[static_spec_arr_grid_size];
    double nuFnu_EC_BLR_grid[static_spec_arr_grid_size];
    double nuFnu_EC_DT_grid[static_spec_arr_grid_size];
    double nuFnu_EC_Disk_grid[static_spec_arr_grid_size];
    double nuFnu_EC_Star_grid[static_spec_arr_grid_size];

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
    unsigned int NU_INT_STOP_Sync_SSC;

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
    double MPI_kernel_delta;
    double MPI_kernel_delta_Emin;

    //--- FREQ BOUNDARIES
    double nu_stop_pp_gamma_pred,nu_stop_pp_gamma;
    double nu_start_pp_gamma;
    double nu_start_pp_gamma_obs;
    double nu_stop_pp_gamma_obs;
    unsigned int NU_INT_STOP_PP_GAMMA;

    double nu_stop_pp_neutrino_pred,nu_stop_pp_neutrino;
    double nu_start_pp_neutrino;
    double nu_start_pp_neutrino_obs;
    double nu_stop_pp_neutrino_obs;
    unsigned int NU_INT_STOP_PP_NUETRINO;
                 

    //
    int do_pp_gamma;
    int do_pp_neutrino;
    int set_pp_racc_gamma, set_pp_racc_elec, set_pp_racc_nu_mu;
    double pp_racc_gamma, pp_racc_elec, pp_racc_nu_mu;
    double E_th_pp_delta_approx,E_pp_x_delta_approx;
    double E_out_e_TeV_pp;

    //--- FREQ/FLUX scalars
    double nu_peak_PP_gamma_blob;
    double nuLnu_peak_PP_gamma_blob;
    double nu_peak_PP_gamma_src;
    double nuLnu_peak_PP_gamma_src;
    double nu_peak_PP_gamma_obs;
    double nuFnu_peak_PP_gamma_obs;

    double nu_peak_PP_neutrino_blob;
    double nuLnu_peak_PP_neutrino_blob;
    double nu_peak_PP_neutrino_src;
    double nuLnu_peak_PP_neutrino_src;
    double nu_peak_PP_neutrino_obs;
    double nuFnu_peak_PP_neutrino_obs;



    //--- FREQ/FLUX array
    
    double j_pp_gamma[static_spec_arr_size];
    double nu_pp_gamma[static_spec_arr_size];
    double nu_pp_gamma_obs[static_spec_arr_size];
    double nuFnu_pp_gamma_obs[static_spec_arr_size];

    double j_pp_neutrino_tot[static_spec_arr_size];
    double j_pp_neutrino_mu[static_spec_arr_size];
    double j_pp_neutrino_e[static_spec_arr_size];

    double nu_pp_neutrino_tot[static_spec_arr_size];
    double nu_pp_neutrino_mu[static_spec_arr_size];
    double nu_pp_neutrino_e[static_spec_arr_size];

    double nu_pp_neutrino_mu_obs[static_spec_arr_size];
    double nu_pp_neutrino_e_obs[static_spec_arr_size];
    double nu_pp_neutrino_tot_obs[static_spec_arr_size];

    double nuFnu_pp_neutrino_tot_obs[static_spec_arr_size];
    double nuFnu_pp_neutrino_mu_obs[static_spec_arr_size];
    double nuFnu_pp_neutrino_e_obs[static_spec_arr_size];

    //--------------------------------//

    //-----------pp-bremss_ep-emission---//
    //--- CONST


    //--- FREQ BOUNDARIES
    double nu_stop_bremss_ep_pred,nu_stop_bremss_ep;
    double nu_start_bremss_ep;
    double nu_start_bremss_ep_obs;
    double nu_stop_bremss_ep_obs;    
    unsigned int NU_INT_STOP_BREMSS_EP;

    //
    int do_bremss_ep;
   
    //--- FREQ/FLUX scalars
    double nu_peak_bremss_ep_blob;
    double nuLnu_peak_bremss_ep_blob;
    double nu_peak_bremss_ep_src;
    double nuLnu_peak_bremss_ep_src;
    double nu_peak_bremss_ep_obs;
    double nuFnu_peak_bremss_ep_obs;

    //--- FREQ/FLUX array
    double j_bremss_ep[static_spec_arr_size];
    double nu_bremss_ep[static_spec_arr_size];
    double nu_bremss_ep_obs[static_spec_arr_size];
    double nuFnu_bremss_ep_obs[static_spec_arr_size];
    //--------------------------------//



    //-----------SSC-IC--------------//
    //--- CONST
    int ord_comp;
    int IC_adaptive_e_binning;
    int do_IC_down_scattering; 
    double COST_IC_K1,COST_IC_COOLING ;

    //--- FREQ BOUNDARIES


    double nu_start_SSC;
    double nu_stop_SSC;
    double nu_start_SSC_obs;
    double nu_stop_SSC_obs;

    unsigned int NU_INT_STOP_COMPTON_SSC;

    //--- IC Kernel computation
    double Gamma;
    double nu; /*freq spettro sinc */
    double nu_1; /*freq spettro comp */
    double nu_compton_0; /* freq campo fot Sync per spettro IC */
    double * nu_seed;
    double * n_seed;
    double g_min_IC;

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
    int do_EC_Disk,do_EC_BLR,do_EC_DT,do_EC_Star,do_EC_CMB,EC_stat,EC_stat_orig;
    int do_Disk,do_DT,do_Star;
    double nu_planck_min_factor;
    double nu_planck_max_factor;
    double mono_planck_min_factor;
    double mono_planck_max_factor;
    //double EC_field_interp_factor;
    unsigned int theta_n_int;
    unsigned int l_n_int;

    double nu_blob_RF;
    double nu_disk_RF;
    double L_nu_disk_RF;

    //dist BLOB/DISK or STAR
    double R_H;
    double R_H_orig;
    double EC_factor;
    double R_H_scale_factor;
    double R_ext_factor;
    double EC_theta_lim;

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

    unsigned int NU_INT_MAX_Star;
	unsigned int NU_INT_STOP_EC_Star;
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
    unsigned int NU_INT_MAX_CMB,NU_INT_STOP_EC_CMB;
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
    //unsigned int NU_INT_MAX_CMB_stat,NU_INT_STOP_EC_CMB_stat;
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
    unsigned int NU_INT_MAX_Disk;
    unsigned int NU_INT_STOP_EC_Disk;

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
    unsigned int NU_INT_MAX_BLR;
    unsigned int NU_INT_STOP_EC_BLR;
    
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
    unsigned int NU_INT_MAX_DT;
    unsigned int NU_INT_STOP_EC_DT;

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
    unsigned int nu_seed_size;
    unsigned int nu_IC_size;
    
    //unsigned int mesh_intComp;
    //unsigned int mesh_intComp1;

    //----------- PARTICLE DISTRIBUTION --------------//
    int Norm_distr;
    //double Norm_distr_L_e_Sync;
    int Distr_e_done;
    int Distr_p_done;
    int Distr_e_pp_done;

    int TIPO_DISTR;
    double *Ne_custom;
    double *gamma_e_custom;
    double *Np_custom;
    double *gamma_p_custom;
    unsigned int gamma_custom_grid_size;
    double *gam;
    double *Ne;
    double *Ne_jetset;
    double *Ne_IC;
    double *Ne_stat;
    double *Np;
    double *Np_jetset;
    double *Q_inj_e_second;
    double *Integrand_over_gamma_grid;
    
    

    unsigned int gamma_grid_size;
    double * griglia_gamma_Ne_log;
    double * griglia_gamma_jetset_Ne_log;
    
    double * griglia_gamma_Ne_log_IC;
    double * griglia_gamma_Ne_log_stat;

    double * griglia_gamma_Np_log;
    double * griglia_gamma_jetset_Np_log;

    double T_esc_e_second;
    //double *griglia_gamma_log_IC;
    //double *N_IC;

    unsigned int i_griglia_gamma;
    double N_tot_e_Sferic;
    double N_tot_p_Sferic;
    double N;
    double N_e_pp, N_p, N_e;
    double NH_cold_to_rel_e;
    double N_0,N_0p,N_0e; /* costante di normalizzazione per distrib elettr staz */
    double gmin;
    double gmax;
    double gmin_secondaries;
    double gmax_secondaries;
    double gamma_cooling_eq;
    int grid_bounded_to_gamma;
    //unsigned int pt_griglia_max;
    double gmin_griglia;
    double gmax_griglia;
    double gmin_griglia_secondaries;
    double gmax_griglia_secondaries;
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
    double U_e, U_p_cold,U_B;
    double U_p, U_p_target;
    double U_Synch, U_Synch_DRF;
    double U_Disk, U_BLR, U_DT, U_CMB;
    double U_Disk_DRF, U_BLR_DRF, U_DT_DRF, U_CMB_DRF;
    double L_Sync_rf, L_SSC_rf, L_EC_Disk_rf,L_EC_BLR_rf, L_EC_DT_rf,L_EC_CMB_rf, L_pp_gamma_rf;
    double jet_L_Sync,jet_L_SSC, jet_L_EC_Disk, jet_L_EC_BLR, jet_L_EC_DT,jet_L_EC_CMB,jet_L_pp_gamma;
    double jet_L_rad,jet_L_kin, jet_L_tot, jet_L_e, jet_L_B, jet_L_p_cold, jet_L_p;
};



//=========================TEMPORAL EVOLUTION ================================

struct temp_ev{
    char STEM[256];
    char path[256];

	
    unsigned int T_COUNTER;

	int do_Sync_cooling;
    int do_Compton_cooling;
    int do_Expansion;
    int do_Adiabatic_cooling;
    //double Q_scaling_factor;

    double t_unit_rad, t_unit_acc; //unit time in light crossing time
    //double t_acc;

    //double *N_gamma;
    double *N_rad_gamma;
    double *N_acc_gamma;
    double *gamma;
    double *gamma_inj_jetset;
    double *N_time;
    double *Q_inj;
    double *Q_inj_jetset;
    double gmin_griglia, gmax_griglia; 
    unsigned int gamma_grid_size;
    unsigned int Q_inj_jetset_gamma_grid_size;
    double g[static_ev_arr_grid_size];
    double t_D[static_ev_arr_grid_size];
    double t_DA[static_ev_arr_grid_size];
    double t_A[static_ev_arr_grid_size];
    double t_Sync_cool[static_ev_arr_grid_size];
    double t_Esc_rad[static_ev_arr_grid_size];
    double t_Esc_acc[static_ev_arr_grid_size];
    //dev start
    double R_t_pre[static_ev_arr_grid_size];
    double B_t_pre[static_ev_arr_grid_size];
    double R_H_t_pre[static_ev_arr_grid_size];
    //double N_exp_fatcor[static_ev_arr_grid_size];
    double time_blob_exp[static_ev_arr_grid_size];
    //dev stop
    double L_inj;
	double *T_esc_acc, *T_esc_rad;
    //double *T_esc_ad_rad;
    double *T_inj_profile;
    double *T_acc_profile;
	double Diff_Coeff;
	double Diff_coeff_CD;
	double Diff_coeff_CA;
    //double m_R;
    double t;
	double Acc_Coeff;
	double Diff_Index;
	double Acc_Index;
    double m_B;
	double T_esc_Coeff_acc;
    double T_esc_Coeff_R_by_c_acc;
    double T_esc_Coeff_rad;
    double T_esc_Coeff_R_by_c_rad;
	double Esc_Index_acc;
    double Esc_Index_rad;
	double Lambda_max_Turb;
	double Lambda_choer_Turb_factor;
	double Gamma_Max_Turb_L_max;
	double Gamma_Max_Turb_L_coher;
	double gamma_eq_t_D, gamma_eq_t_DA ,gamma_eq_t_A;
    double E_acc_max;
    double Delta_R_acc;
    double R_rad_start;
    double R_H_rad_start;
    //double R_jet;
    //double R_jet_exp;
    //double theta_exp_rad;
    //double theta_exp_AR;
    
    double v_exp_by_c;
    double t_jet_exp;
    double R_jet_t;
    double R_H_jet_t;
    double B_acc, B_rad;
    double B_t;

    

	double TStart_Inj;
	double TStop_Inj;
	double TStart_Acc;
	double TStop_Acc;
	double Inj_temp_slope;
	unsigned int NUM_SET;
    unsigned int LOG_SET;
	unsigned int T_SIZE;
    //unsigned int T_EVALUATED;
	double duration,t_D0,t_DA0,t_A0;
    double deltat;

};

//void griglia (double *, double,double,int);
double Cfp(double,struct temp_ev *);
double Bfp(double x, struct temp_ev *pt, struct blob *pt_spec);
double f_Dp(double gamma,struct temp_ev *pt);
double f_DAp(double gamma,struct temp_ev *pt);
double f_A(double gamma,struct temp_ev *pt);


double f_Acc(double x, struct temp_ev *);
double f_Tesc(double, double coeff, double index);
double Inj_temp_prof (double t,struct temp_ev *);
void Init_Q_inj(struct temp_ev *pt_ev );
//double f_Inj(double x,struct temp_ev *);
void Wm (double, double* ,double *);
double Cooling(double,struct temp_ev *,struct blob *);
double Adiabatic_Cooling_time(struct temp_ev *pt, struct blob *pt_spec, double R_jet_t);
int solve_sys1(double VX1[],double VX2[],double VX3[],double SX[],double u[],unsigned int size);
//void free_tempe_ev(struct temp_ev *pt_ev);
void alloc_temp_ev_array(double ** pt,int size);
void CooolingEquilibrium(struct blob * pt, double T_esc);
double IntegrateCooolingEquilibrium( struct blob *pt,double gamma, double T_esc );
double IntegrandCooolingEquilibrium( struct blob *pt, double gamma_1);
double update_jet_expansion(struct blob *pt_spec, struct temp_ev *pt_ev, double t);
double eval_R_H_jet_t(struct blob *pt_spec, struct temp_ev *pt_ev, double time);
double eval_R_jet_t(struct blob *pt_spec, struct temp_ev *pt_ev, double time_blob);
double eval_B_jet_t(struct blob *pt_spec, struct temp_ev *pt_ev,double R_jet_t, double time_blob);
//double time_blob_to_RH(struct temp_ev *pt, struct blob *pt_spec,double R_H_jet_t);
double time_blob_to_obs(double time_blob, struct blob *pt_spec);
double time_obs_to_blob(double time_blob, struct blob *pt_spec);
double eval_N_expansiont_factor(double R_jet_old, double R_jet_new);
void expansion_profile_pre_run(struct blob *pt_spec, struct temp_ev *pt_ev);

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
                        double *dxm);
//===================================================================================

// File Interface
// void FileInput(int argc, char **argv, struct spettro *pt_spec,struct temp_ev *pt_temporal);

//=======================================================================================
/********************************     PyInterface    ************************************/
// PyInterface
struct blob MakeBlob( void );
//void MakeNe(struct spettro *pt_base);
struct temp_ev MakeTempEv( void);
void Init(struct blob *pt, double luminosity_distance);
void InitRadiative(struct blob *pt_base);

//void build_photons(struct spettro *pt_base);
//void alloc_photons(double ** pt,int size);
void set_seed_freq_start(struct blob *pt_base);
void Run_SED(struct blob *pt_base);
void Run_temp_evolution(struct blob *pt_spec_rad, struct blob *pt_spec_acc, struct temp_ev *pt_ev, int only_injection, int do_injection);
void Init_temp_evolution(struct blob *pt_spec_rad, struct blob *pt_spec_acc, struct temp_ev *pt_ev, double luminosity_distance);

double get_spectral_array(double * arr, struct blob *pt, unsigned int id);
double get_array(double * arr, unsigned int id, unsigned int size);
double get_elec_array(double * arr, struct blob *pt, unsigned int id);
double get_temp_ev_N_gamma_array(double *arr, struct temp_ev *pt_ev, unsigned int row, unsigned int col);
double get_temp_ev_N_time_array(double *arr, struct temp_ev *pt_ev, unsigned int id);
double get_temp_ev_gamma_array(double *arr, struct temp_ev *pt_ev, unsigned int id);
double get_Q_inj_array(double *arr, struct temp_ev *pt_ev, unsigned int id);
double get_temp_ev_array_static(double *arr, unsigned int id);

void set_q_inj_user_array(double * arr,struct temp_ev *pt, double val, unsigned int id);
void set_temp_ev_Time_array(double * arr,struct temp_ev *pt, double val, unsigned int id);
void set_elec_custom_array(double * arr, struct blob *pt,double val, unsigned int id);
void set_elec_array(double * arr,struct blob *pt, double val, unsigned int id);
void set_bessel_table(double *arr, struct blob *pt, double val, unsigned int id);

//===================================================================================







//===================================================================================
//ENERGETIC
double GetU_e(struct blob *pt);
double GetE_tot(struct blob *pt);
//void CoolingRates(struct blob * pt, struct temp_ev *pt_ev);
//===================================================================================



//===================================================================================
/************************************ FUNZIONI Distr N *****************************/
void alloc_N_distr(double ** pt,int size);
void setNgrid(struct blob *pt); 
void Fill_N(struct blob *pt, double *griglia_gamma_N_log, double *N);
void build_Ne_custom(struct blob *pt,  unsigned int size);
void build_Np_custom(struct blob *pt,  unsigned int size);

void build_Ne(struct blob *pt);
void build_Np(struct blob *pt);
void build_Ne_jetset(struct blob *pt);
void build_Np_jetset(struct blob *pt);
void Fill_Ne_IC(struct blob *pt, double gmin, int stat_frame);
void InitNe(struct blob *pt);
void Init_Np_Ne_pp(struct blob *pt);
double N_distr(struct blob *, double Gamma);
double N_distr_U_e(struct blob *, double Gamma);
double N_distr_U_p(struct blob *, double Gamma);
double N_distr_integranda(struct blob *, double Gamma);
//void Scrivi_N_file(struct blob *pt, char *name, double *g, double *N);
void FindNe_NpGp(struct blob *pt);
double N_distr_interp(unsigned int size, double Gamma, double *griglia_gamma, double *N);
double Find_gmax(struct blob *pt, double *g, double *N);
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
double N_tot(struct blob *pt, double (*pf_distr)(struct blob *, double x));
//===================================================================================





//===================================================================================
/************************************ FUNZIONI VARIE *******************************/
void messaggio_errore( void);
//void manpage();
void flux_header(FILE *fp);
void flux_DISK_header(FILE *fp);
void distr_e_header(FILE *fp);
void somma_header(FILE *fp);
//===================================================================================



//===================================================================================
/************************************ SUMMED SPECTRA *******************************/
void interpola_somma(struct blob *pt_j, double nu, unsigned int i);
void common_grid_spectra(int num_file, struct blob *);
//===================================================================================




//===================================================================================
/************************************ ENERGETIC *******************************/
double Power_Sync_Electron(struct blob *pt);
double Uph_Sync(struct blob *pt);
double Power_Sync_Electron_Integ(struct blob *pt_N, double Gamma);

double Lum_Sync_at_nu (struct blob *pt, double nu);
double	Lum_SSC_at_nu (struct blob *pt , double nu_1);
double eval_E_acc(double *gamma, double *N, unsigned int gamma_size, double Vol_acc);

void FindEpSp(double * nu_blob, double * nuFnu_obs, unsigned int NU_INT_MAX, struct blob * pt,
        double * nu_peak_obs,
        double * nu_peak_src,
        double * nu_peak_blob,
        double * nuFnu_peak_obs,
        double * nuLnu_peak_src,
        double * nuLnu_peak_blob);


struct jet_energetic EnergeticOutput(struct blob *pt);

void EvalU_e(struct blob *pt);
void EvalU_p(struct blob *pt); 
double I_nu_to_Uph(double * nu, double * I_nu, unsigned int NU_INT_STOP);
double PowerPhotons_blob_rest_frame(struct blob *pt, double *nu, double *nu_Fnu, unsigned int NU_INT_STOP);
double PowerPhotons_disk_rest_frame(struct blob *pt, double *nu, double *nu_Fnu, unsigned int NU_INT_STOP);

//===================================================================================





//======================================================================================
/************************************ GEOMETRY/REL.SPEC. *******************************/
double V_sphere(double R);
double S_sphere(double R);

double get_beaming(double BulkFactor, double theta);
void   SetBeaming(struct blob *pt);
void SetDistr(struct blob *pt);

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
double j_nu_src_to_L_nu_src(double j_nu, double R, double beam_obj);
double nu_blob_to_nu_obs(double nu, double delta, double z);
double nu_blob_to_nu_src(double nu, double delta, double z);
double nu_disk_to_nu_obs_disk(double nu, double z);
double nu_obs_disk_to_nu_disk(double nu, double z);
double nu_obs_to_nu_blob(double nu, double delta, double z);
double nu_obs_to_nu_src(double nu, double z);
//===========================================================================================




//===========================================================================================
/****************************** FUNZIONI SPETTRO SINCROTRONE ********************************/
void spettro_sincrotrone(int Num_file, struct blob *);
double F_K_53(struct blob *, double x);
double F_K_ave(struct blob *, double x);
double F_int_fix(struct blob *, unsigned int ID);
double F_int_ave(struct blob *, unsigned int  ID);
double Sync_self_abs_int(struct blob *, unsigned int ID);
double j_nu_Sync(struct blob *);
double solve_S_nu_Sync(struct blob *pt, unsigned int NU_INT);
double eval_S_nu_Sync(struct blob *pt, double j_nu, double alpha_nu);
double alfa_nu_Sync(struct blob *);
double integrale_Sync(double (*pf)(struct blob *, unsigned int ID), struct blob *pt);
double Sync_tcool(struct blob * , double g);
double Sync_cool(struct blob * , double g);
//===========================================================================================





//===========================================================================================
/****************** FUNZION pp-gamma-e-mu ************************************************/

void spettro_pp_gamma(int Num_file, struct blob *pt);
void spettro_pp_neutrino(int Num_file, struct blob *pt);

double F_gamma(double x, double Ep_TeV);
double F_electrons(double x, double Ep_TeV);
double F_neutrino_mu_1(double x, double Ep_TeV);

double sigma_pp_inel(double Ep_TeV);
double  check_pp_kernel(double res,struct blob *pt,double E_p_TeV, double x );


double rate_gamma_pp(struct blob *pt);
double rate_electrons_pp(struct blob *pt, double Gamma);
double rate_neutrino_mu_1_pp(struct blob *pt, double nu_nu_mu);

double pp_gamma_kernel(double gamma_p, double E_out_TeV,struct blob * pt);
double pp_gamma_kernel_delta(struct blob *pt, double E_pi);

double pp_electron_kernel_delta(struct blob *pt,double E_pi);
double pp_electrons_kernel(double gamma_p, double E_out_TeV, struct blob *pt);

double pp_neturino_mu_1_kernel(double gamma_p, double E_out_TeV, struct blob *pt);
double pp_neutrino_mu_1_kernel_delta(struct blob *pt, double E_pi);

double integrale_pp_second_high_en_rate(double (*pf_pp_kernel) (double gamma_p, double E_out_TeV, struct blob *pt), double E_out_TeV, struct blob * pt, unsigned int i_start);
double integrale_pp_second_low_en_rate(double (*pf_pp_delta_kernel) ( struct blob *pt,double E),
                              double (*E_min_pi) (double gamma_p, struct blob *pt),
                              double (*E_max_pi) (struct blob *pt),  
                              double E_out_TeV,
                              struct blob * pt);


unsigned int E_min_p_grid_even(struct blob *pt, double * gamma_p_grid, double E_start_TeV, unsigned int i_start, unsigned int gamma_p_grid_size );

double f_mu_2_pp(double x, double r);
double g_mu_pp(double x, double r);
double h_mu_1_pp(double x, double r);
double h_mu_2_pp(double x, double r);


double E_min_e_pp(double E_e, struct blob *pt);
double E_max_e_pp(struct blob *pt);
double E_min_neutrino_mu_1_pp(double E_mu, struct blob *pt);
double E_max_neutrino_mu_1_pp(struct blob *pt);
double E_min_gamma_pp(double E_gamma, struct blob *pt);
double E_max_gamma_pp(struct blob *pt);
//===========================================================================================

//===========================================================================================
/****************** FUNZION bremsstrahlung ************************************************/
void spettro_bremss_ep(int Num_file, struct blob *pt);
double j_nu_bremss_ee(struct blob * pt, double nu_out);
double j_nu_bremss_ep(struct blob * pt, double nu_out);
double b_ep_sigma(double gamma_e, double epsilon_gamma);
double b_ee_sigma(double gamma_e, double epsilon_gamma);
double b_ee_sigma_rel(double gamma_e, double epsilon_gamma);
double b_ee_sigma_non_rel(double gamma_e, double epsilon_gamma);
double ee_brem_F(double gamma_e, double x);
double b_ee_A_term(double gamma_e, double epsilon_gamma);
double bremss_sigma_1(double gamma_e, double epsilon_gamma);
double bremss_sigma_2(double gamma_e, double epsilon_gamma);
//===========================================================================================



//========================================================================================
/************************************ PHOTON GRID FUNCTIONS *****************************/
void build_log_grid(double nu_start, double nu_stop, unsigned int SIZE, double *nu_grid);
int x_to_grid_index(double *nu_grid, double nu, unsigned int SIZE);
//===========================================================================================





//=================================================================================
/************************************ FUNCTIONS  IC  *****************************/
double f_compton_K1(struct blob *, double Gamma);
void set_N_distr_for_Compton(struct blob *, double nu_in, double nu_out, int stat_frame);
double rate_compton_GR(struct blob *);
double integrale_IC(struct blob * pt, double a, double b,int stat_frame);
double integrale_IC_cooling(struct blob * pt, double a, double b, double gamma);
double compton_cooling(struct blob *pt_spec, struct temp_ev *pt_ev, double gamma);
double f_compton_cooling(double b);
double I_nu_to_n(double I_nu, double nu);
//==================================================================================




//=============================================================================
/********************* FUNZIONI SPETTRO COMPTON *******************************/
void spettro_compton(int num_file, struct blob *);
//=============================================================================



//===================================================================================
/********************* FUNZIONI EXTERNAL COMPOTON************************************/
void spectra_External_Fields(int Num_file, struct blob *pt_d);
void spettro_EC(int num_file, struct blob *);

void set_EC_stat_pre(struct blob *pt, double R_lim);
void set_EC_stat_post(struct blob *pt);

/***  DISK PLANCK FUNCTIONS  *********/
double eval_T_disk(struct blob *pt, double R);
double eval_nu_peak_planck(double T);
double f_planck_Multi_T(struct blob *pt, double R, double nu);
double f_planck_Multi_T_norm(struct blob *pt, double R,double nu);
double integrand_f_planck_Multi_T(struct blob *pt,double R);
double f_planck(double temperatura, double nu);
double f_planck_norm(double temperatura, double nu);

/***  DISK ACCRETION POWER FUNCTIONS  **********/
double eval_R_Sw(double M_BH);
double eval_M_BH(double R_S);
double eval_accr_rate(double L_disk,double accr_eff);
double eval_L_Edd(double M_BH);
double eval_accr_Edd(double L_Edd, double accr_eff);

/***  EXTERNAL COMPOTON DISK/JET TRANSFORMATIONS  *******/
double eval_nu_min_blob_RF(struct blob *pt, double mu1, double mu2, double nu_disk_RF );
double eval_nu_max_blob_RF(struct blob *pt, double mu1, double mu2, double nu_disk_RF );
double I_nu_disk_RF_to_blob_RF(double I_nu_diks_RF, double nu_disk_RF, double nu_blob_RF, double beta, double Gamma);
double nu_blob_RF_to_nu_disk_RF(double nu_blob_RF, double Gamma, double beta, double mu_disk_RF);


/***  SPECTRAL/GEOMETRIC FUNCTIONS EC Star ****/
void Build_I_nu_Star(struct blob *pt_d);
double eval_I_nu_Star_disk_RF(struct blob *pt,double nu_Star_disk_RF);
//double eval_J_nu_Star_disk_RF(struct spettro *pt, double I_nu_Star_disk_RF);
double eval_I_nu_Star_blob_RF(struct blob *pt, double nu_blob_RF);
double integrand_I_nu_Star_blob_RF(struct blob *pt, double mu);
double eval_Star_L_nu(struct blob *pt, double nu_Star_disk_RF);
void set_Star_geometry(struct blob *pt);

/***  SPECTRAL/GEOMETRIC FUNCTIONS EC DISK ****/
void set_Disk(struct blob *pt);
void set_Disk_angles(struct blob *pt);
void Build_I_nu_Disk(struct blob *pt_d);
double Disk_Spectrum(struct blob *pt, double nu_Disk_disk_RF);

//double eval_I_nu_Disk_disk_RF(struct spettro *pt,double Lnu_Disk_disk_RF);
//double eval_J_nu_Disk_disk_RF(struct spettro *pt, double nu_Disk_disk_RF, double I_nu_Disk_disk_RF);
//double eval_J_nu_Disk_disk_RF_single_BB(struct spettro *pt, double I_nu_Disk_disk_RF);
//double eval_J_nu_Disk_disk_RF_multi_BB(struct spettro *pt,  double  mu);

double eval_I_nu_Disk_disk_RF(struct blob *pt, double nu_disk_RF);
double eval_I_nu_Disk_blob_RF(struct blob *pt, double nu_disk_RF);
double eval_I_nu_theta_Disk(struct blob *pt, double mu);
double integrand_I_nu_Disk_blob_RF(struct blob *pt, double mu);
double integrand_I_nu_Disk_disk_RF(struct blob *pt, double mu);
double eval_Disk_L_nu(struct blob *pt, double nu_Disk_disk_RF);
void set_Disk_geometry(struct blob *pt);
double eval_nu_peak_Disk(double T);

/***  SPECTRAL/GEOMETRIC FUNCTIONS EC CMB ****/
void Build_I_nu_CMB(struct blob *pt_d);
//void Build_I_nu_CMB_stat(struct spettro *pt_d);
double eval_I_nu_CMB_disk_RF(double T_CMB,double nu_CMB_disk_RF);
double eval_I_nu_CMB_blob_RF(struct blob *pt, double nu_blob_RF);
double integrand_I_nu_CMB_blob_RF(struct blob *pt, double mu);
double eval_T_CMB_z(double z, double T_CMB_0);


/*** SPECTRAL/GEOMETRIC FUNCTIONS EC BLR  ****/
void Build_I_nu_BLR(struct blob *pt);
double j_nu_BLR_integrand(struct blob *pt, double l);
double eval_I_nu_theta_BLR(struct blob *pt, double mu);
double integrand_I_nu_BLR_disk_RF(struct blob *, double theta);
double eval_I_nu_BLR_disk_RF(struct blob *pt);
double integrand_I_nu_BLR_blob_RF(struct blob *pt, double theta);
double eval_I_nu_BLR_blob_RF(struct blob *pt);
double eval_theta_max_BLR(struct blob *pt);
void eval_l_values_BLR(struct blob *pt, double mu, double l[]);
double eval_Lnu_BLR_disk_RF(struct blob *pt, double nu_disk_RF);
//double integrand_J_nu_BLR_disk_RF(struct spettro *pt, double mu);
//double integrand_I_nu_BLR_blob_RF(struct spettro *pt, double mu);
void set_BLR_geometry(struct blob *pt);
//double l_BLR(struct spettro *pt, double mu);

/***  FUNCTIONS Seed Photons EC DT  ***/
void Build_I_nu_DT(struct blob *pt);
double j_nu_DT_integrand(struct blob *pt, double l);
double eval_I_nu_theta_DT(struct blob *pt, double mu, double theta);
double integrand_I_nu_DT_disk_RF(struct blob *, double theta);
double eval_I_nu_DT_disk_RF(struct blob *pt);
double eval_I_nu_DT_blob_RF(struct blob *pt);
double integrand_I_nu_DT_blob_RF(struct blob *pt, double theta);
double eval_DT_L_nu(struct blob *pt, double DT_disk_RF);
double eval_theta_max_DT(struct blob *pt);
double eval_l_DT(struct blob *pt, double mu);
double eval_circle_secant(double z,double R,double mu);
//===========================================================================================





//===========================================================================================
/*********************        FUNZIONI MATEMATICHE             *****************************/
double st_gamma(double x);
void beschb(double x, double *gam1, double *gam2, double *gampl,
        double *gammi);
double chebev(float a, float b, float c[], int m, float x);
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);
double bessel_K_53(struct blob *, double x);
double bessel_K_pitch_ave(struct blob *pt, double x);
void tabella_Bessel(struct blob *);
double derivata(double (*pf) (struct blob *, double x), struct blob *pt, double x);
double theta_heaviside(double x);
double integrale_trap_log_struct(double (*pf)(struct blob *, double x),
                                 struct blob *pt, double a, double b, unsigned int intervalli);
double integrale_simp_struct(double (*pf)(struct blob *, double x),
                             struct blob *pt, double a, double b, unsigned int intervalli);
double integrale_simp(double (*pf)(double x), double a, double b, unsigned int n_intervalli);
double integr_simp_grid_equilog(double * x, double *y, unsigned int size);
double trapzd_array_linear_grid(double *x, double *y, unsigned int SIZE);
double trapzd_array_arbritary_grid(double *x, double *y, unsigned int SIZE);
double test_int(struct blob *, double x);
double test_int1(double x);
double log_lin_interp(double nu, double *nu_grid, double nu_min, double nu_max, double *flux_grid, unsigned int SIZE, double emiss_lim);
double log_quad_interp(double x,  double * x_grid, double x_min, double x_max, double *  y_grid , unsigned int SIZE, double emiss_lim);
double log_log_interp(double log_x, double *log_x_grid, double log_x_min, double log_x_max, double *log_y_grid, unsigned int SIZE, double emiss_lim);
//===========================================================================================







//===========================================================================================
/**********************         FUNZIONI DI TEST                     ************************/
double test_lunghezza_vettore(double mesh, double a, double b, int Max_elem);
double test_solid_anlge(struct blob *pt, double mu);
#endif




