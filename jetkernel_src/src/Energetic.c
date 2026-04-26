//=========================================================================================
// Functions to Evaluate Energetic quantities
//=========================================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file Energetic.C
 * \author Andrea Tramacere
 * \date 19-09-2004
 * \brief funzioni per la
 * distribuzione energetica
 * degli elettroni sia nel caso
 * stazionario che per ET
 *
 */
//=========================================================================================
// Eval N(gamma) peaks
//=========================================================================================

void FindNe_NpGp(struct blob *pt) {
    unsigned int i;
    double N2, N3;

    if (pt->emitters.Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
        //Genera_Ne(pt);
    }

    pt->emitters.Gamma_p2 = 0;
    pt->emitters.Gamma_p3 = 0;
    pt->emitters.Np2 = 0;
    pt->emitters.Np3 = 0;
    for (i = 0; i < pt->emitters.gamma_grid_size; i++) {
        N2 = pt->emitters.Ne[i] * pt->emitters.griglia_gamma_Ne_log[i] * pt->emitters.griglia_gamma_Ne_log[i];
        N3 = N2 * pt->emitters.griglia_gamma_Ne_log[i];
        if (N2 > pt->emitters.Np2) {
            pt->emitters.Np2 = N2;
            pt->emitters.Gamma_p2 = pt->emitters.griglia_gamma_Ne_log[i];
        }
        if (N3 > pt->emitters.Np3) {
            pt->emitters.Np3 = N3;
            pt->emitters.Gamma_p3 = pt->emitters.griglia_gamma_Ne_log[i];
        }
    }
}


//=========================================================================================
// Ue and Up Functions
//=========================================================================================

void EvalU_e(struct blob *pt) {
    double (*pf_norm) (struct blob *, double x);


    if (pt->emitters.Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
        //Genera_Ne(pt);
    }
    //printf("N_0=%e\n", pt->emitters.N_0);
    //printf("N=%e N/N_0=%e\n", pt->emitters.N, pt->emitters.N / pt->emitters.N_0);
    pf_norm = &N_distr_U_e;
    pt->emitters.U_e = MEC2 * integrale_trap_log_struct(pf_norm,
            pt,
            pt->emitters.gmin_griglia,
            pt->emitters.gmax_griglia,
            10000);
    pt->emitters.E_tot_e = pt->emitters.U_e * pt->core.Vol_region;
}

void EvalU_p(struct blob *pt) {
    double (*pf_norm) (struct blob *, double x);


    if (pt->emitters.Distr_p_done == 0) {
        printf("No proton distribution calculated \n ");
        exit(0);
        //Genera_Ne(pt);
    }
    //printf("N_0=%e\n", pt->emitters.N_0);
    //printf("N=%e N/N_0=%e\n", pt->emitters.N, pt->emitters.N / pt->emitters.N_0);
    pf_norm = &N_distr_U_p;
    pt->emitters.U_p = MPC2 * integrale_trap_log_struct(pf_norm,
            pt,
            pt->emitters.gmin,
            pt->emitters.gmax,
            10000);
    pt->emitters.E_tot_p = pt->emitters.U_p * pt->core.Vol_region;
}

double GetU_e(struct blob *pt) {
    return pt->emitters.U_e;
}

double GetE_tot(struct blob *pt) {
    return pt->emitters.E_tot_e;
}

//================================================
//N(gamma) Integrand
//=================================================

double N_distr_U_e(struct blob *pt_N, double Gamma) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * Distribuzioni energetiche degli elettroni nel caso statico
     * per calcolare Ue
     *U_e =mec^2*Integ_g1^g2*N(gamma)*gamma
     *mec^2 e' fouri dall'integranda per velocizzare
     *l'integrale
     */
    //return N_distr(pt_N, Gamma) * Gamma;
    //!!!!!! ricordati di che si puo' usare N_distr
    // quando non usi i leptoni secondari
    return N_distr_interp(pt_N->emitters.gamma_grid_size, Gamma, pt_N->emitters.griglia_gamma_Ne_log, pt_N->emitters.Ne) * Gamma;
}

double N_distr_U_p(struct blob *pt_N, double Gamma) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * Distribuzioni energetiche degli elettroni nel caso statico
     * per calcolare Up
     *U_p =mpc^2*Integ_g1^g2*N(gamma)*gamma
     *mpc^2 e' fouri dall'integranda per velocizzare
     *l'integrale
     */
    //return N_distr(pt_N, Gamma) * Gamma;
    //!!!!!! ricordati di che si puo' usare N_distr
    // quando non usi i leptoni secondari
    return N_distr_interp(pt_N->emitters.gamma_grid_size, Gamma, pt_N->emitters.griglia_gamma_Np_log, pt_N->emitters.Np) * Gamma;
}




//========================================================================================




//=========================================================================================
// Find EsSp
//=========================================================================================

void FindEpSp(double * nu_blob, double * nuFnu_obs, unsigned int NU_INT_MAX, struct blob * pt,
        double * nu_peak_obs,
        double * nu_peak_src,
        double * nu_peak_blob,
        double * nuFnu_peak_obs,
        double * nuLnu_peak_src,
        double * nuLnu_peak_blob) {
	unsigned int i;
    
    *nu_peak_obs=nu_blob_to_nu_obs(nu_blob[0], pt->core.beam_obj, pt->core.z_cosm);
    *nu_peak_blob=nu_blob[0];
    *nuFnu_peak_obs = nuFnu_obs[0];
    
    for (i = 0; i <= NU_INT_MAX; i++) {
	//printf ("%e %e\n",nu_blob[i],nuFnu_obs[i]);
        if (nuFnu_obs[i] > *nuFnu_peak_obs) {
            *nuFnu_peak_obs = nuFnu_obs[i];
            *nu_peak_obs = nu_blob_to_nu_obs(nu_blob[i], pt->core.beam_obj, pt->core.z_cosm);
	    *nu_peak_blob=nu_blob[i];
	    //printf("%e %e %e\n",nu_peak_blob,nu_peak_obs,nuFnu_peak_obs);
        }
    }
    
    *nuLnu_peak_src = nuFnu_obs_to_nuLnu_src(*nuFnu_peak_obs, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);
    *nuLnu_peak_blob = nuFnu_obs_to_nuLnu_blob(*nuFnu_peak_obs, pt->core.beam_obj, pt->core.z_cosm, pt->core.dist);

    
    *nu_peak_src = nu_blob_to_nu_src(*nu_peak_blob, pt->core.beam_obj, pt->core.z_cosm);
}
//=========================================================================================


//=========================================================================================
//Total Power integrating  Lnu_blob
//=========================================================================================
//Function to Integrate the Total Power of emitted photons in the blob rest frame

double PowerPhotons_disk_rest_frame(struct blob *pt, double *nu_blob, double *nuFnu, unsigned int NU_INT_STOP)
{
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * Distribuzioni energetiche degli elettroni nel caso statico
     * per calcolare Ue
     */

    double Ptot, P1, P2, nu1, nu2;
    unsigned int i;

    Ptot = 0;
    nu1 = nu_blob[0];
    P1 = nuFnu_obs_to_nuLnu_src(nuFnu[0], pt->core.beam_obj, pt->core.z_cosm, pt->core.dist) / nu1;

    for (i = 1; i <= NU_INT_STOP; i++)
    {
        nu2 = nu_blob[i];
        P2 = nuFnu_obs_to_nuLnu_src(nuFnu[i], pt->core.beam_obj, pt->core.z_cosm, pt->core.dist) / nu2;
        Ptot += (P1 + P2) * (nu2 - nu1);
        nu1 = nu2;
        P1 = P2;
        //printf("%e %e %d\n",nu2, Ptot, i);
    }
    return Ptot * 0.5;
}

double PowerPhotons_blob_rest_frame(struct blob *pt, double *nu_blob, double *nuFnu, unsigned int NU_INT_STOP)
{
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     */

    double Ptot, P1, P2, nu1, nu2;
    unsigned int i;

    Ptot = 0;
    nu1 = nu_blob[0];
    P1 = nuFnu_obs_to_nuLnu_blob(nuFnu[0], pt->core.beam_obj, pt->core.z_cosm, pt->core.dist) / nu1;


    for (i = 1; i <= NU_INT_STOP; i++) {
        nu2 = nu_blob[i];
        P2 = nuFnu_obs_to_nuLnu_blob(nuFnu[i], pt->core.beam_obj, pt->core.z_cosm, pt->core.dist) / nu2;
        Ptot += (P1 + P2)*(nu2 - nu1);
        nu1 = nu2;
        P1 = P2;
        //printf("%e %e %d\n",nu2, Ptot, i);
    }
    return Ptot * 0.5;
}

//====================================================
//IC Lum at a given freq in the blob rest frame
//====================================================
double	Lum_SSC_at_nu (struct blob *pt , double nu_1) {
	double j_comp,q_comp,nuL_nu_comp;
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * calcola la Luminosita' emessa per SSC
     * alla frequenza nu
     */

    //printf("Eval Total Sync Power emitted by electrons\n");
    if (pt->emitters.Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }


    q_comp=rate_compton_GR(pt,nu_1);
    j_comp=q_comp*HPLANCK*nu_1;
    nuL_nu_comp=nu_1*j_nu_to_L_nu_blob(j_comp, pt->core.Vol_region); /*erg s^-1  Hz^-1 */

    return nuL_nu_comp;
}



//====================================================
//S Lum at a given freq in the blob rest frame
//====================================================


double	Lum_Sync_at_nu (struct blob *pt , double nu) {
    double j_nu, alpha_nu, S_nu, nuL_nu_Sync;
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * calcola la Luminosita' emessa per sincrotrone
     * alla frequenza nu
     */

    //printf("Eval Total Sync Power emitted by electrons\n");
    if (pt->emitters.Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }

    j_nu= j_nu_Sync(pt,nu);
    alpha_nu = alfa_nu_Sync(pt,nu);
    S_nu = eval_S_nu_Sync(pt, j_nu, alpha_nu);

    nuL_nu_Sync = I_nu_to_L_nu_blob(S_nu, pt->core.Surf_region)*nu; /*erg s^-1  Hz^-1 */

    return nuL_nu_Sync;
}


double Uph_Sync(struct blob *pt) {
	return I_nu_to_Uph(pt->Sync.spec.nu, pt->Sync.spec.I_nu, pt->Sync.NU_INT_STOP_Sync_SSC);
}



double Power_Sync_Electron(struct blob *pt) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * Distribuzioni energetiche degli elettroni nel caso statico
     * per calcolare Ue
     */

    //printf("Eval Total Sync Power emitted by electrons\n");
    if (pt->emitters.Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }




    double (*pf_N) (struct blob *, double x);
    double a;
    pf_N = &Power_Sync_Electron_Integ;
    a = integrale_trap_log_struct(pf_N,
            pt,
            pt->emitters.gmin_griglia,
            pt->emitters.gmax_griglia,
            10000);
    //printf("%e %e %e \n",a,,pt->core.Vol_region);

        return a * pt->Sync.UB * SIGTH * (four_by_three) * vluce_cm * pt->core.Vol_region *
                pt->Sync.sin_psi * pt->Sync.sin_psi;
}

//==============================
//Power Sync Integrand
//==============================

double Power_Sync_Electron_Integ(struct blob *pt_N, double Gamma) {
    return N_distr_interp(pt_N->emitters.gamma_grid_size, Gamma, pt_N->emitters.griglia_gamma_Ne_log, pt_N->emitters.Ne)
            * Gamma * Gamma;
    //(1.0 - (1.0 / (Gamma * Gamma)));
    //return N_distr(pt_N, Gamma) * Gamma * Gamma * (1.0 - (1.0 / (Gamma * Gamma)));
}
//=========================================================================================


//=========================================================================================
// Uph
//=========================================================================================
//Function to get the U_ph of a given photon field from I_nu 
//U_ph= 4pi* Inu*dnu

double I_nu_to_Uph(double * nu, double * I_nu, unsigned int NU_INT_STOP) {
    double Uph, n_nu1, n_nu2, nu1, nu2;
    unsigned int i;
    Uph = 0;
    nu1 = nu[0];
    n_nu1 = I_nu_to_n(I_nu[0], nu1);
    for (i = 1; i <= NU_INT_STOP; i++) {
        nu2 = nu[i];
        n_nu2 = I_nu_to_n(I_nu[i], nu[i]);
        Uph += (n_nu1 * nu1 + n_nu2 * nu2)*(nu2 - nu1);
        nu1 = nu2;
        n_nu1 = n_nu2;
        //printf("%e %d %e %e\n", Uph, i,nu2,nu1);
    }
    //0.5 from Trapezoidal rule
    return Uph * 0.5 * HPLANCK* four_pi;
}
//=========================================================================================




//=========================================================================================
// Energetic output
//=========================================================================================

struct jet_energetic EnergeticOutput(struct blob * pt) {
    double lum_factor,lum_factor_rad;
    //double L_rad, L_Sync, L_SSC, L_EC_Disk,L_EC_BLR, L_EC_DT, L_PP;
    //double L_kin, L_tot, L_e, L_B, L_p;
    struct jet_energetic energetic;
    //char f_Energetic[static_file_name_max_legth];
    //FILE *fp_Energetic;
    EvalU_e(pt);
    if (strcmp(pt->core.PARTICLE, "protons") == 0){
        EvalU_p(pt);
     }
    //lum_factor and  lum_factor_rad consistent with Eq. 3 and 4, Ghisellini 2010, doi:10.1111/j.1365-2966.2009.15898.x
    lum_factor_rad =0.25 *eval_beta_gamma(pt->core.BulkFactor) * pt->core.BulkFactor * pt->core.BulkFactor ;
    lum_factor = pi * pt->core.R * pt->core.R * vluce_cm * eval_beta_gamma(pt->core.BulkFactor) * pt->core.BulkFactor * pt->core.BulkFactor ;
    energetic.U_B= pt->Sync.UB;
    energetic.U_e= pt->emitters.U_e;
    energetic.jet_L_rad=0.;

    energetic.U_Synch = PowerPhotons_blob_rest_frame (pt, pt->Sync.spec.nu, pt->Sync.spec.nuFnu_obs, pt->Sync.NU_INT_STOP_Sync_SSC)/(4*pi*pt->core.R*pt->core.R*vluce_cm);
    energetic.U_BLR = I_nu_to_Uph(pt->BLR.spec.nu, pt->BLR.spec.I_nu, pt->BLR.spec.NU_INT_MAX);
    energetic.U_DT=I_nu_to_Uph(pt->DT.spec.nu, pt->DT.spec.I_nu, pt->DT.spec.NU_INT_MAX);
    energetic.U_CMB = I_nu_to_Uph(pt->CMB.spec.nu, pt->CMB.spec.I_nu, pt->CMB.spec.NU_INT_MAX);
    energetic.U_Disk = I_nu_to_Uph(pt->Disk.spec.nu, pt->Disk.spec.I_nu, pt->Disk.spec.NU_INT_MAX);
    energetic.U_Star = I_nu_to_Uph(pt->Star.spec.nu, pt->Star.spec.I_nu, pt->Star.spec.NU_INT_MAX);

    energetic.U_seed_tot = energetic.U_Synch+ energetic.U_BLR + energetic.U_DT + energetic.U_Disk + energetic.U_Star;

    energetic.U_Synch_DRF = energetic.U_Synch*(pt->core.beam_obj*pt->core.beam_obj*pt->core.beam_obj*pt->core.beam_obj);
    energetic.U_BLR_DRF = I_nu_to_Uph(pt->BLR.spec.nu_DRF, pt->BLR.spec.I_nu_DRF, pt->BLR.spec.NU_INT_MAX);
    energetic.U_DT_DRF = I_nu_to_Uph(pt->DT.spec.nu_DRF, pt->DT.spec.I_nu_DRF, pt->DT.spec.NU_INT_MAX);
    energetic.U_CMB_DRF = I_nu_to_Uph(pt->CMB.spec.nu_DRF, pt->CMB.spec.I_nu_DRF, pt->CMB.spec.NU_INT_MAX);
    energetic.U_Disk_DRF = I_nu_to_Uph(pt->Disk.spec.nu_DRF, pt->Disk.spec.I_nu_DRF, pt->Disk.spec.NU_INT_MAX);
    energetic.U_Star_DRF =  I_nu_to_Uph(pt->Star.spec.nu_DRF, pt->Star.spec.I_nu_DRF, pt->Star.spec.NU_INT_MAX);

    energetic.L_Sync_rf = PowerPhotons_blob_rest_frame (pt, pt->Sync.spec.nu, pt->Sync.spec.nuFnu_obs, pt->Sync.NU_INT_STOP_Sync_SSC);
    //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
    //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
    energetic.jet_L_Sync = energetic.L_Sync_rf * lum_factor_rad;
    energetic.jet_L_rad = +energetic.jet_L_Sync;
    
    
    if (pt->core.do_SSC) 
    {
        energetic.L_SSC_rf = PowerPhotons_blob_rest_frame(pt, pt->SSC.spec.nu, pt->SSC.spec.nuFnu_obs, pt->SSC.NU_INT_STOP_COMPTON_SSC);
        //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
        //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
        energetic.jet_L_SSC = energetic.L_SSC_rf * lum_factor_rad;
        energetic.jet_L_rad += energetic.jet_L_SSC;
    }   
    else
    {
        energetic.L_SSC_rf=0;
        energetic.jet_L_SSC=0;
    }

    if (strcmp(pt->core.PARTICLE, "protons") == 0) {
        energetic.U_p_target = pt->PP_gamma.NH_pp  * MPC2;
        energetic.U_p = pt->emitters.U_p;
        energetic.L_pp_gamma_rf = PowerPhotons_blob_rest_frame(pt, pt->PP_gamma.spec.nu, pt->PP_gamma.spec.nuFnu_obs, pt->PP_gamma.NU_INT_STOP_PP_GAMMA);
        //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
        //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
        energetic.jet_L_pp_gamma = energetic.L_pp_gamma_rf* lum_factor_rad;
        energetic.jet_L_rad += energetic.jet_L_pp_gamma;
        energetic.U_p_cold = 0.;
    }
    else
    {
        energetic.U_p_cold = pt->emitters.N * pt->emitters.NH_cold_to_rel_e * MPC2;
        energetic.U_p = 0.;
        energetic.U_p_target = 0.;
        energetic.L_pp_gamma_rf=0.;
        energetic.jet_L_pp_gamma=0.;
    }

    if (pt->core.do_EC_Disk == 1 ) {
        energetic.L_EC_Disk_rf = PowerPhotons_blob_rest_frame(pt, pt->Disk.ec.spec.nu, pt->Disk.ec.spec.nuFnu_obs, pt->Disk.ec.NU_INT_STOP);
        //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
        //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
        energetic.jet_L_EC_Disk = energetic.L_EC_Disk_rf * lum_factor_rad;
        energetic.jet_L_rad += energetic.jet_L_EC_Disk;
    }
    else
    {
        energetic.jet_L_EC_Disk = 0;
        energetic.L_EC_Disk_rf = 0;
    }

    if (pt->core.do_EC_Disk == 1 || pt->core.do_EC_BLR == 1)
    {
        energetic.L_EC_BLR_rf = PowerPhotons_blob_rest_frame(pt, pt->BLR.ec.spec.nu, pt->BLR.ec.spec.nuFnu_obs, pt->BLR.ec.NU_INT_STOP);
        //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
        //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
        energetic.jet_L_EC_BLR = energetic.L_EC_BLR_rf * lum_factor_rad;
        energetic.jet_L_rad += energetic.jet_L_EC_BLR;
    }
    else
    {
        energetic.L_EC_BLR_rf=0;
        energetic.jet_L_EC_BLR=0;
    }
    
    if (pt->core.do_EC_DT == 1) {
        energetic.L_EC_DT_rf = PowerPhotons_blob_rest_frame(pt, pt->DT.ec.spec.nu, pt->DT.ec.spec.nuFnu_obs, pt->DT.ec.NU_INT_STOP);
        //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
        //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
        energetic.jet_L_EC_DT = energetic.L_EC_DT_rf * lum_factor_rad;
        energetic.jet_L_rad += energetic.jet_L_EC_DT;
    }
    else
    {
        energetic.jet_L_EC_DT = 0;
        energetic.L_EC_DT_rf = 0;
    }
    
    if (pt->core.do_EC_CMB == 1)
    {
        energetic.L_EC_CMB_rf = PowerPhotons_blob_rest_frame(pt, pt->CMB.ec.spec.nu, pt->CMB.ec.spec.nuFnu_obs, pt->CMB.ec.NU_INT_STOP);
        //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
        //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
        energetic.jet_L_EC_CMB = energetic.L_EC_CMB_rf * lum_factor_rad;
        energetic.jet_L_rad += energetic.jet_L_EC_CMB;
    }
    else
    {
        energetic.jet_L_EC_CMB = 0;
        energetic.L_EC_CMB_rf = 0;
    }

    if (pt->core.do_EC_Star == 1)
    {
        energetic.L_EC_Star_rf = PowerPhotons_blob_rest_frame(pt, pt->Star.ec.spec.nu, pt->Star.ec.spec.nuFnu_obs, pt->Star.ec.NU_INT_STOP);
        //NOTE: PowerPhotons_blob_rest_frame*lum_factor_rad already takes into account
        //NOTE: U=L/(4 pi R^2 c) and pi R^2 U, the R^2 and pi cancel out  
        energetic.jet_L_EC_Star = energetic.L_EC_Star_rf * lum_factor_rad;
        energetic.jet_L_rad += energetic.jet_L_EC_Star;
    }
    else
    {
        energetic.L_EC_Star_rf = 0;
        energetic.jet_L_EC_Star = 0;
    }
    
    energetic.jet_L_e = pt->emitters.U_e * lum_factor;
    energetic.jet_L_p = lum_factor * energetic.U_p;
    energetic.jet_L_p_cold = lum_factor * energetic.U_p_cold;
    energetic.jet_L_B = pt->Sync.UB * lum_factor;
    energetic.jet_L_kin = energetic.jet_L_e + energetic.jet_L_p_cold + energetic.jet_L_p;
    energetic.jet_L_tot = energetic.jet_L_kin + energetic.jet_L_rad + energetic.jet_L_B ;

   
    return energetic;
}
