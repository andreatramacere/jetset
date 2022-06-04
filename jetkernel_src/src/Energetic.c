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

    if (pt->Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
        //Genera_Ne(pt);
    }

    pt->Gamma_p2 = 0;
    pt->Gamma_p3 = 0;
    pt->Np2 = 0;
    pt->Np3 = 0;
    for (i = 0; i < pt->gamma_grid_size; i++) {
        N2 = pt->Ne[i] * pt->griglia_gamma_Ne_log[i] * pt->griglia_gamma_Ne_log[i];
        N3 = N2 * pt->griglia_gamma_Ne_log[i];
        if (N2 > pt->Np2) {
            pt->Np2 = N2;
            pt->Gamma_p2 = pt->griglia_gamma_Ne_log[i];
        }
        if (N3 > pt->Np3) {
            pt->Np3 = N3;
            pt->Gamma_p3 = pt->griglia_gamma_Ne_log[i];
        }
    }
}


//=========================================================================================
// Ue and Up Functions
//=========================================================================================

void EvalU_e(struct blob *pt) {
    double (*pf_norm) (struct blob *, double x);


    if (pt->Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
        //Genera_Ne(pt);
    }
    //printf("N_0=%e\n", pt->N_0);
    //printf("N=%e N/N_0=%e\n", pt->N, pt->N / pt->N_0);
    pf_norm = &N_distr_U_e;
    pt-> U_e = MEC2 * integrale_trap_log_struct(pf_norm,
            pt,
            pt->gmin_griglia,
            pt->gmax_griglia,
            10000);
    pt->E_tot_e = pt->U_e * pt->Vol_sphere;
}

void EvalU_p(struct blob *pt) {
    double (*pf_norm) (struct blob *, double x);


    if (pt->Distr_p_done == 0) {
        printf("No proton distribution calculated \n ");
        exit(0);
        //Genera_Ne(pt);
    }
    //printf("N_0=%e\n", pt->N_0);
    //printf("N=%e N/N_0=%e\n", pt->N, pt->N / pt->N_0);
    pf_norm = &N_distr_U_p;
    pt-> U_p = MPC2 * integrale_trap_log_struct(pf_norm,
            pt,
            pt->gmin,
            pt->gmax,
            10000);
    pt->E_tot_p = pt->U_p * pt->Vol_sphere;
}

double GetU_e(struct blob *pt) {
    return pt->U_e;
}

double GetE_tot(struct blob *pt) {
    return pt->E_tot_e;
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
    return N_distr_interp(pt_N->gamma_grid_size, Gamma, pt_N->griglia_gamma_Ne_log, pt_N->Ne) * Gamma;
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
    return N_distr_interp(pt_N->gamma_grid_size, Gamma, pt_N->griglia_gamma_Np_log, pt_N->Np) * Gamma;
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
    
    *nu_peak_obs=nu_blob_to_nu_obs(nu_blob[0], pt->beam_obj, pt->z_cosm);
    *nu_peak_blob=nu_blob[0];
    *nuFnu_peak_obs = nuFnu_obs[0];
    
    for (i = 0; i <= NU_INT_MAX; i++) {
	//printf ("%e %e\n",nu_blob[i],nuFnu_obs[i]);
        if (nuFnu_obs[i] > *nuFnu_peak_obs) {
            *nuFnu_peak_obs = nuFnu_obs[i];
            *nu_peak_obs = nu_blob_to_nu_obs(nu_blob[i], pt->beam_obj, pt->z_cosm);
	    *nu_peak_blob=nu_blob[i];
	    //printf("%e %e %e\n",nu_peak_blob,nu_peak_obs,nuFnu_peak_obs);
        }
    }
    
    *nuLnu_peak_src = nuFnu_obs_to_nuLnu_src(*nuFnu_peak_obs, pt->beam_obj, pt->z_cosm, pt->dist);
    *nuLnu_peak_blob = nuFnu_obs_to_nuLnu_blob(*nuFnu_peak_obs, pt->beam_obj, pt->z_cosm, pt->dist);

    
    *nu_peak_src = nu_blob_to_nu_src(*nu_peak_blob, pt->beam_obj, pt->z_cosm);
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
    P1 = nuFnu_obs_to_nuLnu_src(nuFnu[0], pt->beam_obj, pt->z_cosm, pt->dist) / nu1;

    for (i = 1; i <= NU_INT_STOP; i++)
    {
        nu2 = nu_blob[i];
        P2 = nuFnu_obs_to_nuLnu_src(nuFnu[i], pt->beam_obj, pt->z_cosm, pt->dist) / nu2;
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
     * Distribuzioni energetiche degli elettroni nel caso statico
     * per calcolare Ue
     */

    double Ptot, P1, P2, nu1, nu2;
    unsigned int i;

    Ptot = 0;
    nu1 = nu_blob[0];
    P1 = nuFnu_obs_to_nuLnu_blob(nuFnu[0], pt->beam_obj, pt-> z_cosm, pt->dist) / nu1;


    for (i = 1; i <= NU_INT_STOP; i++) {
        nu2 = nu_blob[i];
        P2 = nuFnu_obs_to_nuLnu_blob(nuFnu[i], pt->beam_obj, pt-> z_cosm, pt->dist) / nu2;
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
    if (pt->Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }


    pt->nu_1=nu_1;
    q_comp=rate_compton_GR(pt);
    j_comp=q_comp*HPLANCK*nu_1;
    nuL_nu_comp=nu_1*j_nu_to_L_nu_blob(j_comp, pt->Vol_sphere); /*erg s^-1  Hz^-1 */

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
    if (pt->Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }

    pt->nu=nu;
    j_nu= j_nu_Sync(pt);
    alpha_nu = alfa_nu_Sync(pt);
    S_nu = eval_S_nu_Sync(pt, j_nu, alpha_nu);

    nuL_nu_Sync = I_nu_to_L_nu_blob(S_nu, pt->Surf_sphere)*nu; /*erg s^-1  Hz^-1 */

    return nuL_nu_Sync;
}


double Uph_Sync(struct blob *pt) {
	return I_nu_to_Uph(pt->nu_Sync, pt->I_nu_Sync, pt->NU_INT_STOP_Sync_SSC);
}



double Power_Sync_Electron(struct blob *pt) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * Distribuzioni energetiche degli elettroni nel caso statico
     * per calcolare Ue
     */

    //printf("Eval Total Sync Power emitted by electrons\n");
    if (pt->Distr_e_done == 0) {
        printf("No electron distribution calculated \n ");
        exit(0);
    }




    double (*pf_N) (struct blob *, double x);
    double a;
    pf_N = &Power_Sync_Electron_Integ;
    a = integrale_trap_log_struct(pf_N,
            pt,
            pt->gmin_griglia,
            pt->gmax_griglia,
            10000);
    //printf("%e %e %e \n",a,,pt->Vol_sphere);

        return a * pt->UB * SIGTH * (four_by_three) * vluce_cm * pt->Vol_sphere *
                pt->sin_psi * pt->sin_psi;
}

//==============================
//Power Sync Integrand
//==============================

double Power_Sync_Electron_Integ(struct blob *pt_N, double Gamma) {
    return N_distr_interp(pt_N->gamma_grid_size, Gamma, pt_N->griglia_gamma_Ne_log, pt_N->Ne)
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
    double lum_factor;
    //double L_rad, L_Sync, L_SSC, L_EC_Disk,L_EC_BLR, L_EC_DT, L_PP;
    //double L_kin, L_tot, L_e, L_B, L_p;
    struct jet_energetic energetic;
    //char f_Energetic[static_file_name_max_legth];
    //FILE *fp_Energetic;

    energetic.U_B= pt->UB;
    energetic.U_e= pt->U_e;
    energetic.jet_L_rad=0.;

    energetic.U_Synch = Uph_Sync(pt);
    energetic.U_BLR = I_nu_to_Uph(pt->nu_BLR, pt->I_nu_BLR, pt->NU_INT_MAX_BLR);
    energetic.U_DT=I_nu_to_Uph(pt->nu_DT, pt->I_nu_DT, pt->NU_INT_MAX_DT);
    energetic.U_CMB = I_nu_to_Uph(pt->nu_CMB, pt->I_nu_CMB, pt->NU_INT_MAX_CMB);
    energetic.U_Disk = I_nu_to_Uph(pt->nu_Disk, pt->I_nu_Disk, pt->NU_INT_MAX_Disk);

    energetic.U_Synch_DRF = Uph_Sync(pt) * pt->beam_obj * pt->beam_obj * pt->beam_obj * pt->beam_obj;
    energetic.U_BLR_DRF = I_nu_to_Uph(pt->nu_BLR_disk_RF, pt->I_nu_BLR_disk_RF, pt->NU_INT_MAX_BLR);
    energetic.U_DT_DRF = I_nu_to_Uph(pt->nu_DT_disk_RF, pt->I_nu_DT_disk_RF, pt->NU_INT_MAX_DT);
    energetic.U_CMB_DRF = I_nu_to_Uph(pt->nu_CMB_disk_RF, pt->I_nu_CMB_disk_RF, pt->NU_INT_MAX_CMB);
    energetic.U_Disk_DRF = I_nu_to_Uph(pt->nu_Disk_disk_RF, pt->I_nu_Disk_disk_RF, pt->NU_INT_MAX_Disk);

    energetic.L_Sync_rf = PowerPhotons_blob_rest_frame (pt, pt->nu_Sync, pt->nuF_nu_Sync_obs, pt->NU_INT_STOP_Sync_SSC);
    energetic.jet_L_Sync = energetic.L_Sync_rf * 0.25 * pt->BulkFactor * pt->BulkFactor;
    energetic.jet_L_rad = +energetic.jet_L_Sync;
    
    if (pt->do_SSC) 
    {
        energetic.L_SSC_rf = PowerPhotons_blob_rest_frame(pt, pt->nu_SSC, pt->nuF_nu_SSC_obs, pt->NU_INT_STOP_COMPTON_SSC);
        energetic.jet_L_SSC = energetic.L_SSC_rf * 0.25 * pt->BulkFactor * pt->BulkFactor;
        energetic.jet_L_rad += energetic.jet_L_SSC;
    }   
    else
    {
        energetic.L_SSC_rf=0;
        energetic.jet_L_SSC=0;
    }

    if (strcmp(pt->PARTICLE, "protons") == 0) {
        energetic.U_p_target = pt->NH_pp  * MPC2;
        energetic.U_p = pt->U_p;
        energetic.L_pp_gamma_rf = PowerPhotons_blob_rest_frame(pt, pt->nu_pp_gamma, pt->nuFnu_pp_gamma_obs, pt->NU_INT_STOP_PP_GAMMA);
        energetic.jet_L_pp_gamma = energetic.L_pp_gamma_rf* 0.25 * pt->BulkFactor * pt->BulkFactor;
        energetic.jet_L_rad += energetic.jet_L_pp_gamma;
        energetic.U_p_cold = 0.;
    }
    else
    {
        energetic.U_p_cold = pt->N * pt->NH_cold_to_rel_e * MPC2;
        energetic.U_p = 0.;
        energetic.U_p_target = 0.;
        energetic.L_pp_gamma_rf=0.;
        energetic.jet_L_pp_gamma=0.;
    }

    if (pt->do_EC_Disk == 1 ) {
        energetic.L_EC_Disk_rf = PowerPhotons_blob_rest_frame(pt, pt->nu_EC_Disk, pt->nuF_nu_EC_Disk_obs, pt->NU_INT_STOP_EC_Disk);
        energetic.jet_L_EC_Disk = energetic.L_EC_Disk_rf * 0.25 * pt->BulkFactor * pt->BulkFactor;
        energetic.jet_L_rad += energetic.jet_L_EC_Disk;
    }
    else
    {
        energetic.jet_L_EC_Disk = 0;
        energetic.L_EC_Disk_rf = 0;
    }

    if (pt->do_EC_Disk == 1 || pt->do_EC_BLR == 1)
    {
        energetic.L_EC_BLR_rf = PowerPhotons_blob_rest_frame(pt, pt->nu_EC_BLR, pt->nuF_nu_EC_BLR_obs, pt->NU_INT_STOP_EC_BLR);
        energetic.jet_L_EC_BLR = energetic.L_EC_BLR_rf * 0.25 * pt->BulkFactor * pt->BulkFactor;
        energetic.jet_L_rad += energetic.jet_L_EC_BLR;
    }
    else
    {
        energetic.L_EC_BLR_rf=0;
        energetic.jet_L_EC_BLR=0;
    }
    
    if (pt->do_EC_DT == 1) {
        energetic.L_EC_DT_rf = PowerPhotons_blob_rest_frame(pt, pt->nu_EC_DT, pt->nuF_nu_EC_DT_obs, pt->NU_INT_STOP_EC_DT);
        energetic.jet_L_EC_DT = energetic.L_EC_DT_rf * 0.25 * pt->BulkFactor * pt->BulkFactor;
        energetic.jet_L_rad += energetic.jet_L_EC_DT;
    }
    else
    {
        energetic.jet_L_EC_DT = 0;
        energetic.L_EC_DT_rf = 0;
    }
    
    if (pt->do_EC_CMB == 1)
    {
        energetic.L_EC_CMB_rf = PowerPhotons_blob_rest_frame(pt, pt->nu_EC_CMB, pt->nuF_nu_EC_CMB_obs, pt->NU_INT_STOP_EC_CMB);
        energetic.jet_L_EC_CMB = energetic.L_EC_CMB_rf * 0.25 * pt->BulkFactor * pt->BulkFactor;
        energetic.jet_L_rad += energetic.jet_L_EC_CMB;
    }
    else
    {
        energetic.jet_L_EC_CMB = 0;
        energetic.L_EC_CMB_rf = 0;
    }
    lum_factor = pi * pt->R * pt->R * vluce_cm * eval_beta_gamma(pt->BulkFactor) * pt->BulkFactor * pt->BulkFactor ;
    energetic.jet_L_e = pt->U_e * lum_factor;
    energetic.jet_L_p = lum_factor * energetic.U_p;
    energetic.jet_L_p_cold = lum_factor * energetic.U_p_cold;
    energetic.jet_L_B = pt->UB * lum_factor;
    energetic.jet_L_kin = energetic.jet_L_e + energetic.jet_L_p_cold + energetic.jet_L_p;
    energetic.jet_L_tot = energetic.jet_L_kin + energetic.jet_L_rad + energetic.jet_L_B ;

   
    return energetic;
}

