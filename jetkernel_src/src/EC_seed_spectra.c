//===============================================================
//
//                                FUNZIONI SPETTRO DISCO
//===============================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"
/**
 * \file spettro_disco.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief FUNZIONI SPETTRO DISCO
 *
 */

static size_t angle_dep_flat_index(unsigned int nu_id, unsigned int angle_id, unsigned int angle_n_int)
{
	return ((size_t)nu_id) * ((size_t)angle_n_int) + ((size_t)angle_id);
}

void reset_external_spectrum_angle_dep(struct spectrum_external *spec)
{
	if (spec == NULL) {
		return;
	}
	spec->angle_n_int = 0U;
	spec->angle_nu_size = 0U;
	spec->mu = NULL;
	spec->theta = NULL;
	spec->I_nu_theta = NULL;
	spec->I_nu_theta_DRF = NULL;
	spec->n_nu_theta = NULL;
	spec->n_nu_theta_DRF = NULL;
}

void free_external_spectrum_angle_dep(struct spectrum_external *spec)
{
	if (spec == NULL) {
		return;
	}

	if (spec->mu != NULL) {
		free(spec->mu);
		spec->mu = NULL;
	}
	if (spec->theta != NULL) {
		free(spec->theta);
		spec->theta = NULL;
	}
	if (spec->I_nu_theta != NULL) {
		free(spec->I_nu_theta);
		spec->I_nu_theta = NULL;
	}
	if (spec->I_nu_theta_DRF != NULL) {
		free(spec->I_nu_theta_DRF);
		spec->I_nu_theta_DRF = NULL;
	}
	if (spec->n_nu_theta != NULL) {
		free(spec->n_nu_theta);
		spec->n_nu_theta = NULL;
	}
	if (spec->n_nu_theta_DRF != NULL) {
		free(spec->n_nu_theta_DRF);
		spec->n_nu_theta_DRF = NULL;
	}
	spec->angle_n_int = 0U;
	spec->angle_nu_size = 0U;
}

int ensure_external_spectrum_angle_dep(struct spectrum_external *spec, unsigned int nu_size, unsigned int angle_n_int)
{
	size_t angle_size;
	size_t flat_size;
	int alloc_needed;

	if (spec == NULL) {
		return -1;
	}

	if ((nu_size == 0U) || (angle_n_int == 0U)) {
		free_external_spectrum_angle_dep(spec);
		return 0;
	}

	alloc_needed = 0;
	if ((spec->angle_n_int != angle_n_int) ||
		(spec->angle_nu_size != nu_size) ||
		(spec->mu == NULL) ||
		(spec->theta == NULL) ||
		(spec->I_nu_theta == NULL) ||
		(spec->I_nu_theta_DRF == NULL) ||
		(spec->n_nu_theta == NULL) ||
		(spec->n_nu_theta_DRF == NULL)) {
		alloc_needed = 1;
	}

	if (alloc_needed) {
		free_external_spectrum_angle_dep(spec);

		spec->mu = (double *)calloc((size_t)angle_n_int, sizeof(double));
		spec->theta = (double *)calloc((size_t)angle_n_int, sizeof(double));
		spec->I_nu_theta = (double *)calloc(((size_t)nu_size) * ((size_t)angle_n_int), sizeof(double));
		spec->I_nu_theta_DRF = (double *)calloc(((size_t)nu_size) * ((size_t)angle_n_int), sizeof(double));
		spec->n_nu_theta = (double *)calloc(((size_t)nu_size) * ((size_t)angle_n_int), sizeof(double));
		spec->n_nu_theta_DRF = (double *)calloc(((size_t)nu_size) * ((size_t)angle_n_int), sizeof(double));

		if ((spec->mu == NULL) ||
			(spec->theta == NULL) ||
			(spec->I_nu_theta == NULL) ||
			(spec->I_nu_theta_DRF == NULL) ||
			(spec->n_nu_theta == NULL) ||
			(spec->n_nu_theta_DRF == NULL)) {
			free_external_spectrum_angle_dep(spec);
			return -1;
		}

		spec->angle_n_int = angle_n_int;
		spec->angle_nu_size = nu_size;
	}

	angle_size = (size_t)spec->angle_n_int;
	flat_size = ((size_t)spec->angle_n_int) * ((size_t)spec->angle_nu_size);
	memset(spec->mu, 0, angle_size * sizeof(double));
	memset(spec->theta, 0, angle_size * sizeof(double));
	memset(spec->I_nu_theta, 0, flat_size * sizeof(double));
	memset(spec->I_nu_theta_DRF, 0, flat_size * sizeof(double));
	memset(spec->n_nu_theta, 0, flat_size * sizeof(double));
	memset(spec->n_nu_theta_DRF, 0, flat_size * sizeof(double));

	return 0;
}

void reset_blob_external_spectra_angle_dep(struct blob *pt)
{
	if (pt == NULL) {
		return;
	}
	reset_external_spectrum_angle_dep(&(pt->Disk.spec));
	reset_external_spectrum_angle_dep(&(pt->BLR.spec));
	reset_external_spectrum_angle_dep(&(pt->DT.spec));
	reset_external_spectrum_angle_dep(&(pt->Corona.spec));
	reset_external_spectrum_angle_dep(&(pt->Star.spec));
	reset_external_spectrum_angle_dep(&(pt->CMB.spec));
}

void free_blob_external_spectra_angle_dep(struct blob *pt)
{
	if (pt == NULL) {
		return;
	}
	free_external_spectrum_angle_dep(&(pt->Disk.spec));
	free_external_spectrum_angle_dep(&(pt->BLR.spec));
	free_external_spectrum_angle_dep(&(pt->DT.spec));
	free_external_spectrum_angle_dep(&(pt->Corona.spec));
	free_external_spectrum_angle_dep(&(pt->Star.spec));
	free_external_spectrum_angle_dep(&(pt->CMB.spec));
}

//===============================================================
// Evaluation of external radiative fields
//===============================================================

void spectra_External_Fields(int Num_file, struct blob *pt, int set_EC){

    //==================================================================
	//if (pt->core.verbose){
	if (pt->core.verbose > 0)
	{
		printf("**********************   Eval. seed photon fields for  EC       *******************************\n");
	}


    // ====================================================
    // approx gamma con beaming factor, impilcit teta circa= 1/gamma
    // Set nu start EC
	// not used in photon field computation
	// used only in analytic approx on screen
    //=====================================================
	pt->core.beaming_EC = pt->core.BulkFactor;

	//printf("spectra_External_Fields 1  R_H_orig=%e, R_H=%e\n", pt->core.R_H_orig, pt->core.R_H);
	if (pt->core.do_EC_Star==1 || pt->core.do_Star==1){
		//if (set_EC==1){
		//	set_EC_stat_pre(pt, -1);
		//}
    	Build_I_nu_Star(pt);
		//if (set_EC == 1)
		//{
		//	set_EC_stat_post(pt);
		//}
	}
	if (pt->core.do_EC_Disk == 1 || pt->core.do_EC_BLR == 1 || pt->core.do_Disk == 1 || pt->core.do_EC_DT == 1 || pt->core.do_DT ==1)
	{
		//if (set_EC == 1){
		//	set_EC_stat_pre(pt, pt->Disk.R_ext);
		//}
		Build_I_nu_Disk(pt);
		//if (set_EC == 1)
		//{
		//	set_EC_stat_post(pt);
		//}
	}
    if (pt->core.do_EC_BLR==1){
		//if (set_EC == 1)
		//{
		//	set_EC_stat_pre(pt, pt->BLR.R_BLR_out);
		//}
		Build_I_nu_BLR(pt);
		//if (set_EC == 1)
		//{
		//	set_EC_stat_post(pt);
		//}
	}
	if (pt->core.do_EC_DT==1 || pt->core.do_DT==1){
		//printf("EC_stat=%d, R_H=%e\n",pt->core.EC_stat,pt->core.R_H);
		//if (set_EC == 1)
		//{
		//	set_EC_stat_pre(pt, pt->DT.R_DT);
		//}
		Build_I_nu_DT(pt);
		//if (set_EC == 1)
		//{
		//	set_EC_stat_post(pt);
		//}
	}
	if (pt->core.do_EC_Corona==1 || pt->core.do_Corona==1){
		Build_I_nu_Corona(pt);
	}
	if (pt->core.do_EC_CMB==1){
		//if (set_EC == 1)
		//{
		//	set_EC_stat_pre(pt, -1);
		//}
		Build_I_nu_CMB(pt);
		//if (set_EC == 1)
		//{
		//	set_EC_stat_post(pt);
		//}
	}
	//if (pt->do_EC_CMB_stat==1){
    //	Build_I_nu_CMB_stat(pt);
    //}
	//printf("spectra_External_Fields 2  R_H_orig=%e, R_H=%e\n", pt->core.R_H_orig, pt->core.R_H);
	if (pt->core.verbose > 1)
	{
		printf("#-> ********************************\n\n");
	}
}
//=========================================================================================


//=========================================================================================
void Build_I_nu_Star(struct blob *pt){
	//char f_SED_star[static_file_name_max_legth];
	//FILE *fp_SED_star;
	double nu_peak_BB,nu_obs;
	unsigned int NU_INT,NU_INT_MAX;
	unsigned int ANGLE_INT, ANGLE_INT_MAX;
	int have_angle_storage;
	size_t angle_idx;
	double d_mu, mu_grid;
	double nu_start_disk_RF;
	double nu_stop_disk_RF;
	double nuL_nu_disk,F_nu_disk_obs;

	/*
	sprintf(f_SED_star, "%s%s-SED-star.dat",pt->core.path, pt->core.STEM);

	if (pt->core.WRITE_TO_FILE==1){
		fp_SED_star = fopen(f_SED_star, "w");
		if (fp_SED_star == NULL) {
			printf("unable to open %s\n ", fp_SED_star);
			exit(1);
		}
		flux_DISK_header(fp_SED_star);
	}
	*/

	set_Star_geometry(pt);

	//pt->Star_mu_1=0;
	//pt->Star_mu_2=1;
	
   
	nu_peak_BB=eval_nu_peak_Disk(pt->Star.T_Star);

	nu_start_disk_RF = nu_peak_BB*pt->core.nu_planck_min_factor;
	nu_stop_disk_RF  = nu_peak_BB*pt->core.nu_planck_max_factor;

	pt->Star.spec.nu_min = eval_nu_min_blob_RF(pt,pt->Star.mu_star, pt->Star.mu_star, nu_start_disk_RF);
	pt->Star.spec.nu_max  = eval_nu_max_blob_RF(pt,pt->Star.mu_star, pt->Star.mu_star, nu_stop_disk_RF);

	pt->Star.spec.nu_min_DRF = nu_start_disk_RF;
	pt->Star.spec.nu_max_DRF = nu_stop_disk_RF;


	NU_INT_MAX=pt->core.nu_seed_size-1;
	pt->Star.spec.NU_INT_MAX = NU_INT_MAX;
	have_angle_storage = 0;
	if (ensure_external_spectrum_angle_dep(&(pt->Star.spec), pt->core.nu_seed_size, pt->core.theta_n_int) == 0 &&
		pt->Star.spec.angle_n_int > 0U) {
		have_angle_storage = 1;
		ANGLE_INT_MAX = pt->Star.spec.angle_n_int - 1U;
		d_mu = (ANGLE_INT_MAX > 0U) ? (2.0 / (double)ANGLE_INT_MAX) : 0.0;
		for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
			mu_grid = -1.0 + d_mu * (double)ANGLE_INT;
			if (mu_grid > 1.0) {
				mu_grid = 1.0;
			}
			pt->Star.spec.mu[ANGLE_INT] = mu_grid;
			pt->Star.spec.theta[ANGLE_INT] = acos(mu_grid);
		}
	}


	pt->Star.spec.nu_min_obs=nu_disk_to_nu_obs_disk(nu_start_disk_RF , pt->core.z_cosm);
	pt->Star.spec.nu_max_obs=nu_disk_to_nu_obs_disk(nu_stop_disk_RF, pt->core.z_cosm);

	if (pt->core.verbose)
	{
		printf("-----------  Building I_nu Star     ----------- \n");

		printf("nu_start_Star=%e  nu_stop_Star=%e \n",
			   pt->Star.spec.nu_min,
			   pt->Star.spec.nu_max);

		printf("nu_start_Star_disk_RF=%e  nu_stop_Star_disk_RF=%e \n",
			   nu_start_disk_RF,
			   nu_stop_disk_RF);

		printf("nu_start_Star_obs=%e  nu_stop_Star_obs=%e \n",
			   pt->Star.spec.nu_min_obs,
			   pt->Star.spec.nu_max_obs);
		
	}
	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->core.nu_seed_size, pt->Star.spec.nu_DRF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->Star.spec.I_nu_DRF[NU_INT]=eval_I_nu_Star_disk_RF(pt, pt->Star.spec.nu_DRF[NU_INT]);
		//pt->Star.spec.J_nu_DRF[NU_INT]=eval_J_nu_Star_disk_RF(pt, pt->Star.spec.I_nu_DRF[NU_INT]);
	}



	build_log_grid( pt->Star.spec.nu_min,  pt->Star.spec.nu_max, pt->core.nu_seed_size, pt->Star.spec.nu);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		nu_obs = nu_disk_to_nu_obs_disk(pt->Star.spec.nu_DRF[NU_INT],pt->core.z_cosm);
		pt->Star.spec.nu_obs[NU_INT]=nu_obs;
		pt->Star.spec.I_nu[NU_INT]=eval_I_nu_Star_blob_RF(pt,pt->Star.spec.nu[NU_INT]);
		pt->Star.spec.n_nu[NU_INT] =I_nu_to_n(pt->Star.spec.I_nu[NU_INT], pt->Star.spec.nu[NU_INT]);
		//EC with n(gamma) transf
		pt->Star.spec.n_nu_DRF[NU_INT] = I_nu_to_n(pt->Star.spec.I_nu_DRF[NU_INT], pt->Star.spec.nu_DRF[NU_INT]);
		
		if (pt->Star.spec.I_nu[NU_INT]>pt->core.emiss_lim){
			pt->Star.spec.nu_max = pt->Star.spec.nu[NU_INT];
			pt->Star.spec.NU_INT_MAX = NU_INT;
		}
		else{
			pt->Star.spec.I_nu[NU_INT]=pt->core.emiss_lim;
			pt->Star.spec.n_nu[NU_INT] =I_nu_to_n(pt->Star.spec.I_nu[NU_INT], pt->Star.spec.nu[NU_INT]);

		}
		if (have_angle_storage) {
			for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
				angle_idx = angle_dep_flat_index(NU_INT, ANGLE_INT, pt->Star.spec.angle_n_int);
				pt->Star.spec.I_nu_theta_DRF[angle_idx] = pt->Star.spec.I_nu_DRF[NU_INT];
				pt->Star.spec.I_nu_theta[angle_idx] = pt->Star.spec.I_nu[NU_INT];
				pt->Star.spec.n_nu_theta_DRF[angle_idx] = pt->Star.spec.n_nu_DRF[NU_INT];
				pt->Star.spec.n_nu_theta[angle_idx] = pt->Star.spec.n_nu[NU_INT];
			}
		}

		nuL_nu_disk = eval_Star_L_nu(pt,pt->Star.spec.nu_DRF[NU_INT]) * pt->Star.spec.nu_DRF[NU_INT];
		F_nu_disk_obs= L_nu_Disk_to_F_nu(nuL_nu_disk / pt->Star.spec.nu_DRF[NU_INT], pt->core.z_cosm, pt->core.dist);
		pt->Star.spec.nuFnu_obs[NU_INT] = F_nu_disk_obs*nu_obs;
		if (pt->core.verbose > 1)
		{
			printf(" nu_Star_disk_RF=%e, nuF_nu_Star_obs=%e, nu_Star=%e, , I_nu_Star=%e,  nuL_nu_disk=%e, Star surface=%e nu_Star_obs=%e\n",
				   pt->Star.spec.nu_DRF[NU_INT],
				   pt->Star.spec.nuFnu_obs[NU_INT],
				   pt->Star.spec.nu[NU_INT],
				   pt->Star.spec.I_nu_DRF[NU_INT],
				   nuL_nu_disk,
				   pt->Star.Star_surface,
				   pt->Star.spec.nu_obs[NU_INT]);
		}

		/*
		if (pt->core.WRITE_TO_FILE==1){
			fprintf(fp_SED_star, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
				log10(nu_obs),
				log10(nu_obs * F_nu_disk_obs),
				nu_obs,
				nu_obs*F_nu_disk_obs,
				pt->Star.spec.nu_DRF[NU_INT],
				nuL_nu_disk);
		}
		*/

	}
	
	/*
	if (pt->core.WRITE_TO_FILE == 1)
	{
		fclose(fp_SED_star);
	}
	*/
}


//========================
// Star Spectral Functions
//========================

double eval_I_nu_Star_disk_RF(struct blob *pt,double nu_Star_disk_RF){
	return eval_Star_L_nu(pt,nu_Star_disk_RF)/(16*pi*pi*pt->Star.R_H_Star*pt->Star.R_H_Star);
}

// double integrand_I_nu_Star_blob_RF(struct blob *pt, double mu){
// 	int i;
// 	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->core.nu_blob_RF,pt->core.BulkFactor,pt->core.beta_Gamma,mu);

// 	i=x_to_grid_index( pt->Star.spec.nu_DRF,nu_disk_RF,pt->core.nu_seed_size);
// 	if (i>0){
// 		return pt->Star.spec.I_nu_DRF[i]*pt->core.BulkFactor*(1-pt->core.beta_Gamma*mu);
// 	}
// 	else{
// 		return 0;
// 	}
// }

double eval_I_nu_Star_blob_RF(struct blob *pt, double nu_blob_RF){
	int i;
	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(nu_blob_RF,pt->core.BulkFactor,pt->core.beta_Gamma,pt->Star.mu_star);
	i=x_to_grid_index( pt->Star.spec.nu_DRF,nu_disk_RF,pt->core.nu_seed_size);
	if (i>0){
		return pt->Star.spec.I_nu_DRF[i]*pt->core.BulkFactor*(1-pt->core.beta_Gamma*pt->Star.mu_star);
	}
	else{
		return 0;
	}
}

double eval_Star_L_nu(struct blob *pt, double nu_Star_disk_RF){
	return  pi*pt->Star.Star_surface *f_planck(pt->Star.T_Star, nu_Star_disk_RF);
}

double eval_Star_L(struct blob *pt, double T_Star){
	return  sigma_steph_boltz *T_Star*T_Star*T_Star*T_Star*pt->Star.Star_surface;
}


//========================
// Star Geometrical Functions
//========================

void set_Star_geometry(struct blob *pt){
	//double theta_c;
	pt->Star.theta_c_Star=asin(pt->Star.theta_Star/pt->Star.R_H_Star);
	
	// b=sqrt(pt->core.R_H*pt->core.R_H - pt->Star.R_Star*pt->Star.R_Star);	
	// mu1=b/pt->core.R_H;

	// pt->Star_mu_1=min(mu1,mu2);
	// pt->Star_mu_2=max(mu1,mu2);

	// if (pt->core.verbose){
	//printf("theta_c_Star=%20.20e\n",pt->Star.theta_c_Star);
	// 
	pt->Star.mu_star = cos(pt->Star.theta_Star*M_PI / 180.0);
	pt->Star.R_Star=sqrt(pt->Star.L_Star/(4*pi*pt->Star.T_Star*pt->Star.T_Star*pt->Star.T_Star*pt->Star.T_Star*sigma_steph_boltz));
	pt->Star.Star_surface=4*pi*pt->Star.R_Star*pt->Star.R_Star;
}
//=========================================================================================


//=========================================================================================
void Build_I_nu_CMB(struct blob *pt){
	double T_CMB_z;
	double nu_peak_CMB_z;
	unsigned int NU_INT,NU_INT_MAX;
	unsigned int ANGLE_INT, ANGLE_INT_MAX;
	int have_angle_storage;
	int i_mu;
	size_t angle_idx;
	double d_mu, mu_grid, nu_disk_RF, I_blob_theta;
	double nu_start_disk_RF;
	double nu_stop_disk_RF;

	pt->CMB.CMB_mu_1=-1.0;
	pt->CMB.CMB_mu_2=1.0;

	T_CMB_z=eval_T_CMB_z(pt->core.z_cosm,pt->CMB.T_CMB_0);
	//T_CMB_0=pt->CMB.T_CMB_0;

	nu_peak_CMB_z=eval_nu_peak_planck(T_CMB_z);
	//nu_peak_CMB_0=eval_nu_peak_planck(T_CMB_0);

	nu_start_disk_RF = nu_peak_CMB_z*pt->core.nu_planck_min_factor;
	nu_stop_disk_RF  = nu_peak_CMB_z*pt->core.nu_planck_max_factor;

	pt->CMB.spec.nu_min = eval_nu_min_blob_RF(pt,-1, 1, nu_start_disk_RF);
	pt->CMB.spec.nu_max  = eval_nu_max_blob_RF(pt,-1, 1, nu_stop_disk_RF);

	pt->CMB.spec.nu_min_DRF = nu_start_disk_RF;
	pt->CMB.spec.nu_max_DRF = nu_stop_disk_RF;
	//pt->nu_start_CMB_obs=nu_peak_CMB_0*pt->core.nu_planck_min_factor;
	//pt->nu_stop_CMB_obs=nu_peak_CMB_0*pt->core.nu_planck_max_factor;

	NU_INT_MAX=pt->core.nu_seed_size-1;
	pt->CMB.spec.NU_INT_MAX = NU_INT_MAX;
	have_angle_storage = 0;
	if (ensure_external_spectrum_angle_dep(&(pt->CMB.spec), pt->core.nu_seed_size, pt->core.theta_n_int) == 0 &&
		pt->CMB.spec.angle_n_int > 0U) {
		have_angle_storage = 1;
		ANGLE_INT_MAX = pt->CMB.spec.angle_n_int - 1U;
		d_mu = (ANGLE_INT_MAX > 0U) ? (2.0 / (double)ANGLE_INT_MAX) : 0.0;
		for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
			mu_grid = -1.0 + d_mu * (double)ANGLE_INT;
			if (mu_grid > 1.0) {
				mu_grid = 1.0;
			}
			pt->CMB.spec.mu[ANGLE_INT] = mu_grid;
			pt->CMB.spec.theta[ANGLE_INT] = acos(mu_grid);
		}
	}

	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->core.nu_seed_size, pt->CMB.spec.nu_DRF);
	
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
			pt->CMB.spec.I_nu_DRF[NU_INT]=eval_I_nu_CMB_disk_RF(T_CMB_z, pt->CMB.spec.nu_DRF[NU_INT]);
	}
	build_log_grid( pt->CMB.spec.nu_min,  pt->CMB.spec.nu_max, pt->core.nu_seed_size, pt->CMB.spec.nu);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->CMB.spec.I_nu[NU_INT]=eval_I_nu_CMB_blob_RF(pt,pt->CMB.spec.nu[NU_INT]);
		pt->CMB.spec.n_nu[NU_INT] =I_nu_to_n(pt->CMB.spec.I_nu[NU_INT], pt->CMB.spec.nu[NU_INT]);
		//EC with n(gamma) transf
		pt->CMB.spec.n_nu_DRF[NU_INT] = I_nu_to_n(pt->CMB.spec.I_nu_DRF[NU_INT], pt->CMB.spec.nu_DRF[NU_INT]);
		if (have_angle_storage) {
			for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
				angle_idx = angle_dep_flat_index(NU_INT, ANGLE_INT, pt->CMB.spec.angle_n_int);
				mu_grid = pt->CMB.spec.mu[ANGLE_INT];
				pt->CMB.spec.I_nu_theta_DRF[angle_idx] = pt->CMB.spec.I_nu_DRF[NU_INT];
				pt->CMB.spec.n_nu_theta_DRF[angle_idx] = pt->CMB.spec.n_nu_DRF[NU_INT];

				nu_disk_RF = nu_blob_RF_to_nu_disk_RF(pt->CMB.spec.nu[NU_INT],
				                                      pt->core.BulkFactor,
				                                      pt->core.beta_Gamma,
				                                      mu_grid);
				i_mu = x_to_grid_index(pt->CMB.spec.nu_DRF, nu_disk_RF, pt->core.nu_seed_size);
				if (i_mu > 0) {
					I_blob_theta = pt->CMB.spec.I_nu_DRF[i_mu] * pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu_grid);
				}
				else{
					I_blob_theta = 0.0;
				}

				pt->CMB.spec.I_nu_theta[angle_idx] = I_blob_theta;
				pt->CMB.spec.n_nu_theta[angle_idx] = I_nu_to_n(I_blob_theta, pt->CMB.spec.nu[NU_INT]);
			}
		}
	}
	
}



double eval_T_CMB_z(double z, double T_CMB_0){
		return T_CMB_0*(1+z);
}

double eval_I_nu_CMB_disk_RF(double T_CMB,double nu_CMB_disk_RF){
	return f_planck(T_CMB, nu_CMB_disk_RF);
}


double eval_I_nu_CMB_blob_RF(struct blob *pt, double nu_blob_RF){


	pt->core.nu_blob_RF=nu_blob_RF;
	double (*pf) (struct blob *, double x);
	pf = &integrand_I_nu_CMB_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return 0.5 * integrale_simp_struct(pf, pt, pt->CMB.CMB_mu_1, pt->CMB.CMB_mu_2, pt->core.theta_n_int);
}

double integrand_I_nu_CMB_blob_RF(struct blob *pt, double mu){
	int i=0;
 	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->core.nu_blob_RF,pt->core.BulkFactor,pt->core.beta_Gamma,mu);
	i=x_to_grid_index( pt->CMB.spec.nu_DRF,nu_disk_RF,pt->core.nu_seed_size);
	if (i>0){
		return pt->CMB.spec.I_nu_DRF[i]*pt->core.BulkFactor*(1-pt->core.beta_Gamma*mu);
	}
	else{
		return 0;
	}
}

//=========================================================================================







//=========================================================================================
void Build_I_nu_Disk(struct blob *pt){

	//char f_SED_disk[static_file_name_max_legth];
	//FILE *fp_SED_disk;
	double nu_peak_BB,nu_obs;
	unsigned int NU_INT,NU_INT_MAX;
	unsigned int ANGLE_INT, ANGLE_INT_MAX;
	int have_angle_storage;
	size_t angle_idx;
	double d_mu, mu_grid;
	double R_H_orig_angle, R_H_eval_angle, c_angle;
	double I_theta_DRF, I_theta_blob;
	double nu_start_disk_RF;
	double nu_stop_disk_RF;
	double nuL_nu_disk,F_nu_disk_obs;
	//printf("=> Ciccio 1\n");
	if (pt->core.verbose){
		printf("-----------  Building I_nu disk     ----------- \n");
	}

	/*
	if (pt->core.WRITE_TO_FILE==1){
		sprintf(f_SED_disk, "%s%s-SED-disk.dat",pt->core.path, pt->core.STEM);

		fp_SED_disk = fopen(f_SED_disk, "w");
		if (fp_SED_disk == NULL) {
			printf("unable to open %s\n ", f_SED_disk);
			exit(1);
		}
		flux_DISK_header(fp_SED_disk);
	}
	*/
	set_Disk(pt);
	set_Disk_geometry(pt);
	set_Disk_angles(pt);
	if (pt->core.disk == 1)
	{
		nu_peak_BB=eval_nu_peak_Disk(pt->Disk.T_Disk);
		nu_start_disk_RF = nu_peak_BB*pt->core.nu_planck_min_factor;
		nu_stop_disk_RF  = nu_peak_BB*pt->core.nu_planck_max_factor;
	}
	else if (pt->core.disk == 2)
	{
		nu_peak_BB=eval_nu_peak_Disk(pt->Disk.T_Disk);
		nu_start_disk_RF = nu_peak_BB*pt->core.nu_planck_min_factor;
		nu_stop_disk_RF  = nu_peak_BB*pt->core.nu_planck_max_factor;
		//double (*pf) (struct spettro *, double x);
		//pf = &Disk_Spectrum;
		//pt->Cost_Norm_disk_Mulit_BB= 1.0/
		//		integrale_simp_struct(pf, pt,nu_start_disk_RF, nu_stop_disk_RF, pt->core.theta_n_int);
		//printf( "%e\n",pt->Cost_Norm_disk_Mulit_BB);
	 }
	 else if (pt->core.disk == 3)
	 {
		 nu_peak_BB = eval_nu_peak_Disk(pt->Disk.T_Disk);
		 nu_start_disk_RF = nu_peak_BB * pt->core.mono_planck_min_factor;
		 nu_stop_disk_RF = nu_peak_BB * pt->core.mono_planck_max_factor;
	}
	else{
		printf("wrong disk type, option BB, MultiBB, Mono \n ");
		exit(1);
	}

	//if (pt->corona==1){
	//	pt->nu_stop_disk_RF=1E21;
	//}

	pt->Disk.spec.nu_min = eval_nu_min_blob_RF(pt, pt->Disk.Disk_mu_1, pt->Disk.Disk_mu_2, nu_start_disk_RF);
	pt->Disk.spec.nu_max = eval_nu_max_blob_RF(pt,pt->Disk.Disk_mu_1, pt->Disk.Disk_mu_2, nu_stop_disk_RF);

	pt->Disk.spec.nu_min_DRF = nu_start_disk_RF;
	pt->Disk.spec.nu_max_DRF = nu_stop_disk_RF;

		if (pt->core.verbose)
	{
		printf("nu_start_Disk=%e  nu_stop_Disk=%e \n",
			pt->Disk.spec.nu_min,
			pt->Disk.spec.nu_max);

		printf("nu_start_Disk_disk_RF=%e  nu_stop_Disk_disk_RF=%e \n",
			nu_start_disk_RF,
			nu_stop_disk_RF);
	}


	NU_INT_MAX=pt->core.nu_seed_size-1;
	pt->Disk.spec.NU_INT_MAX = NU_INT_MAX;
	have_angle_storage = 0;
	R_H_orig_angle = pt->core.R_H;
	R_H_eval_angle = R_H_orig_angle;
	c_angle = 1.0;
	if (ensure_external_spectrum_angle_dep(&(pt->Disk.spec), pt->core.nu_seed_size, pt->core.theta_n_int) == 0 &&
		pt->Disk.spec.angle_n_int > 0U) {
		have_angle_storage = 1;
		if (R_H_eval_angle > pt->Disk.R_Disk_interp) {
			R_H_eval_angle = pt->Disk.R_Disk_interp;
			c_angle = (pt->Disk.R_Disk_interp / R_H_orig_angle) * (pt->Disk.R_Disk_interp / R_H_orig_angle);
		}
		pt->core.R_H = R_H_eval_angle;
		set_Disk_angles(pt);
		ANGLE_INT_MAX = pt->Disk.spec.angle_n_int - 1U;
		d_mu = (ANGLE_INT_MAX > 0U) ? ((pt->Disk.Disk_mu_2 - pt->Disk.Disk_mu_1) / (double)ANGLE_INT_MAX) : 0.0;
		for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
			mu_grid = pt->Disk.Disk_mu_1 + d_mu * (double)ANGLE_INT;
			if (mu_grid < -1.0) {
				mu_grid = -1.0;
			}
			else if (mu_grid > 1.0) {
				mu_grid = 1.0;
			}
			pt->Disk.spec.mu[ANGLE_INT] = mu_grid;
			pt->Disk.spec.theta[ANGLE_INT] = acos(mu_grid);
		}
		pt->core.R_H = R_H_orig_angle;
		set_Disk_angles(pt);
	}


	pt->Disk.spec.nu_min_obs=nu_disk_to_nu_obs_disk(nu_start_disk_RF , pt->core.z_cosm);
	pt->Disk.spec.nu_max_obs=nu_disk_to_nu_obs_disk(nu_stop_disk_RF, pt->core.z_cosm);

	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->core.nu_seed_size, pt->Disk.spec.nu_DRF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->Disk.spec.L_nu_DRF[NU_INT] = eval_Disk_L_nu(pt, pt->Disk.spec.nu_DRF[NU_INT]);
		
	}
	for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++)
	{
		pt->Disk.spec.I_nu_DRF[NU_INT] = eval_I_nu_Disk_disk_RF(pt, pt->Disk.spec.nu_DRF[NU_INT]);		
	}

	build_log_grid( pt->Disk.spec.nu_min,  pt->Disk.spec.nu_max, pt->core.nu_seed_size, pt->Disk.spec.nu);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
 		nu_obs = nu_disk_to_nu_obs_disk(pt->Disk.spec.nu_DRF[NU_INT],pt->core.z_cosm);
		pt->Disk.spec.nu_obs[NU_INT]=nu_obs;
		pt->Disk.spec.I_nu[NU_INT] = eval_I_nu_Disk_blob_RF(pt, pt->Disk.spec.nu_DRF[NU_INT]);
		pt->Disk.spec.n_nu[NU_INT] =I_nu_to_n(pt->Disk.spec.I_nu[NU_INT], pt->Disk.spec.nu[NU_INT]);
		//EC with n(gamma) transf
		pt->Disk.spec.n_nu_DRF[NU_INT] = I_nu_to_n(pt->Disk.spec.I_nu_DRF[NU_INT], pt->Disk.spec.nu_DRF[NU_INT]);

		if (pt->core.verbose>1){
			printf(" nu_Disk_disk_RF=%e, I_nu_Disk_disk_RF=%e, nu_Disk=%e, , I_nu_Disk=%e\n",
				   pt->Disk.spec.nu_DRF[NU_INT],
				   pt->Disk.spec.I_nu_DRF[NU_INT],
				   pt->Disk.spec.nu[NU_INT],
				   pt->Disk.spec.I_nu[NU_INT]);
		}

		if (pt->Disk.spec.I_nu[NU_INT]>pt->core.emiss_lim){
			pt->Disk.spec.nu_max = pt->Disk.spec.nu[NU_INT];
			pt->Disk.spec.NU_INT_MAX = NU_INT;
		}
		else{
			pt->Disk.spec.I_nu[NU_INT]=pt->core.emiss_lim;
			pt->Disk.spec.n_nu[NU_INT] =I_nu_to_n(pt->Disk.spec.I_nu[NU_INT], pt->Disk.spec.nu[NU_INT]);

		}

		nuL_nu_disk = pt->Disk.spec.L_nu_DRF[NU_INT] * pt->Disk.spec.nu_DRF[NU_INT];
		F_nu_disk_obs= L_nu_Disk_to_F_nu(nuL_nu_disk / pt->Disk.spec.nu_DRF[NU_INT], pt->core.z_cosm, pt->core.dist);
		pt->Disk.spec.nuFnu_obs[NU_INT] = F_nu_disk_obs*nu_obs;
		/*
		if (pt->core.WRITE_TO_FILE==1){
			fprintf(fp_SED_disk, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
				log10(nu_obs),
				log10(nu_obs * F_nu_disk_obs),
				nu_obs,
				nu_obs*F_nu_disk_obs,
				pt->Disk.spec.nu_DRF[NU_INT],
				nuL_nu_disk);
		}
		*/


	}
	if (have_angle_storage) {
		pt->core.R_H = R_H_eval_angle;
		set_Disk_angles(pt);
		for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++) {
			pt->core.nu_disk_RF = pt->Disk.spec.nu_DRF[NU_INT];
			for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
				mu_grid = pt->Disk.spec.mu[ANGLE_INT];
				I_theta_DRF = c_angle * eval_I_nu_theta_Disk(pt, mu_grid);
				I_theta_blob = I_theta_DRF * pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu_grid);
				angle_idx = angle_dep_flat_index(NU_INT, ANGLE_INT, pt->Disk.spec.angle_n_int);
				pt->Disk.spec.I_nu_theta_DRF[angle_idx] = I_theta_DRF;
				pt->Disk.spec.I_nu_theta[angle_idx] = I_theta_blob;
				pt->Disk.spec.n_nu_theta_DRF[angle_idx] = I_nu_to_n(I_theta_DRF, pt->Disk.spec.nu_DRF[NU_INT]);
				pt->Disk.spec.n_nu_theta[angle_idx] = I_nu_to_n(I_theta_blob, pt->Disk.spec.nu[NU_INT]);
			}
		}
		pt->core.R_H = R_H_orig_angle;
		set_Disk_angles(pt);
	}
	pt->Disk.L_Disk_radiative = PowerPhotons_disk_rest_frame(pt, pt->Disk.spec.nu_DRF, pt->Disk.spec.nuFnu_obs, pt->Disk.spec.NU_INT_MAX);
	
	/*
	if (pt->core.WRITE_TO_FILE==1){
		fclose(fp_SED_disk);
	}
	*/
}



void set_Disk(struct blob *pt){
	double  nu_peak_BB;
	if (strcmp(pt->core.disk_type, "BB") == 0)
	{
		pt->core.disk = 1;
	}
	else if (strcmp(pt->core.disk_type, "MultiBB") == 0)
	{
		pt->core.disk = 2;
	}
	else if (strcmp(pt->core.disk_type, "Mono") == 0)
	{
		pt->core.disk = 3;
	}
	else
	{
		printf("wrong disk type, option BB, MultiBB, Mono \n ");
		exit(1);
	}

	pt->Disk.R_Sw=eval_R_Sw(pt->Disk.M_BH);
	//R_inner
	pt->Disk.R_inner = pt->Disk.R_inner_Sw * pt->Disk.R_Sw;
	//R_ext
	pt->Disk.R_ext = pt->Disk.R_ext_Sw * pt->Disk.R_Sw;
	pt->Disk.R_Disk_interp =   pt->Disk.R_ext*50.0;
	pt->Disk.L_Edd = eval_L_Edd(pt->Disk.M_BH);
	pt->Disk.accr_rate = eval_accr_rate(pt->Disk.L_Disk, pt->Disk.accr_eff);
	pt->Disk.accr_Edd = eval_accr_Edd(pt->Disk.L_Edd, pt->Disk.accr_eff);
	//as in Ghisellini 2009, but it is equivalent to the one in Eq. 5.43 in the Frank, King & Raine Book
	//but we use 8p in place of 16pi, because L_Disk is the L of uno disk surface
	//this must be evaluated before eval_T_disk
	pt->Disk.Cost_disk_Mulit_BB = pt->Disk.R_inner * pt->Disk.L_Disk / (8 * pi * sigma_steph_boltz * pt->Disk.accr_eff);

	if (pt->core.disk == 2)
	//multi BB
	{
		pt->Disk.T_Disk = eval_T_disk(pt, (49. / 36.) * pt->Disk.R_inner);
		//printf("Cost_disk_Mulit_BB = %e \n", pt->Disk.Cost_disk_Mulit_BB);
	}
	
	nu_peak_BB = eval_nu_peak_Disk(pt->Disk.T_Disk);

	if (pt->core.verbose){
		printf("T_max = %e (K)\n",pt->Disk.T_Disk);
		// energy corresponding to Tmax
		printf("E_max = %e (eV)\n",pt->Disk.T_Disk*K_boltz*erg_to_eV);
		// frequency corresponding to Tmax
		printf("nu_max = %e (Hz)\n",pt->Disk.T_Disk*K_boltz/HPLANCK);
		//Peak of the BB spectrum
		printf("nu_peak  = %e (Hz)\n",nu_peak_BB);
		printf("schwarzschild radius=%e\n", pt->Disk.R_Sw);
		printf("R_ext =%e (cm)\n", pt->Disk.R_ext);
		printf("R_inner =%e (cm)\n", pt->Disk.R_inner);

		printf("Black hole mass = %e (m_sun)\n", pt->Disk.M_BH);

		printf("Accr. rate = %e (g/s)\n", pt->Disk.accr_rate);
		printf("Accr. rate = %e (M_sun/year)\n", pt->Disk.accr_rate * 86400. * 365. / m_sun);
		printf("L_Edd = %e (erg/s)\n", pt->Disk.L_Edd);
		printf("L_Disk = %e (erg/s)\n", pt->Disk.L_Disk);
		printf("L_diks/L_edd = %e\n", pt->Disk.L_Disk / pt->Disk.L_Edd);

		printf("Accr_Edd = %e (g/s)\n", pt->Disk.accr_Edd);
		printf("Accr_Edd = %e (M_sun/year)\n", pt->Disk.accr_Edd * 86400. * 365. / m_sun);
	}
}

//========================
// Disk Spectral Functions
//========================


double Disk_Spectrum(struct blob *pt, double nu_Disk_disk_RF){
	double I;
	double (*pf)(struct blob *, double x);
	I=0;
	if (pt->core.disk == 1) {
		// in this case we use a normalized planck function
		I= f_planck_norm(pt->Disk.T_Disk, nu_Disk_disk_RF);
	}
	else if (pt->core.disk == 2) {
		//in this case we acutally integrate every annluar BB along the disk
		
		pf = &integrand_f_planck_Multi_T;
		pt->Disk.nu_disk_Multi_BB = nu_Disk_disk_RF;
		//printf("=> pt->Disk.nu_disk_Multi_BB %e\n",pt->Disk.nu_disk_Multi_BB);
		//pi is the angular part for a disk face
		I= pi *integrale_trap_log_struct(pf, pt, pt->Disk.R_inner * 1.01, pt->Disk.R_ext, 100);
	}
	else if (pt->core.disk==3){
		I= eval_nu_peak_Disk(pt->Disk.T_Disk)*(pt->core.mono_planck_max_factor-pt->core.mono_planck_min_factor);
	}
	return I*cos(pt->core.theta * Deg_to_Rad);
}

double eval_I_nu_theta_Disk(struct blob *pt, double mu)
{
	//double (*pf)(struct spettro *, double x);
	//unsigned int i;
	double  I,R,R_D;
	//pf = &j_nu_BLR_integrand;
	//pt->BLR.mu_j = mu;
	I=0;
    if (pt->core.disk == 1) {
		// in this case we use a normalized planck function
		I = f_planck_norm(pt->Disk.T_Disk, pt->core.nu_disk_RF)*pt->Disk.L_Disk * pt->Disk.Disk_geom_factor;
	}
	else if (pt->core.disk == 2) {
		//in this case we acutally integrate every annluar BB along the disk
		
	
		R=pt->core.R_H/mu;
		R_D = sqrt(R * R - pt->core.R_H * pt->core.R_H);
		I = f_planck_Multi_T(pt, R_D, pt->core.nu_disk_RF)/pi;
	}
	else if (pt->core.disk==3){
		I= eval_nu_peak_Disk(pt->Disk.T_Disk)*(pt->core.mono_planck_max_factor-pt->core.mono_planck_min_factor);
	}

	
	return I;
}

double integrand_I_nu_Disk_blob_RF(struct blob *pt, double mu)
{
	//double psi, sin_theta;
	//sin_theta=sqrt(1.0 - mu*mu);
	double f;
	f=  (pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu));
	//f=1/( (pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu)) * (pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu)) );
	return 2 * pi  * eval_I_nu_theta_Disk(pt, mu) *f;
}

double integrand_I_nu_Disk_disk_RF(struct blob *pt, double mu)
{
	//double psi, sin_theta;
	//sin_theta = sqrt(1.0 - mu * mu);
	//printf("=> %e %e\n", sin_theta, eval_I_nu_theta_Disk(pt, mu));
	return 2 * pi  * eval_I_nu_theta_Disk(pt, mu);
}

double eval_I_nu_Disk_blob_RF(struct blob *pt, double nu_disk_RF)
{
	double (*pf)(struct blob *, double x);
	double I,c,R_H_orig;
	//unsigned int i;
	pt->core.nu_disk_RF = nu_disk_RF;
	pf = &integrand_I_nu_Disk_blob_RF;

	c = 1.0;
	R_H_orig = pt->core.R_H;
	if (pt->core.R_H > pt->Disk.R_Disk_interp)
	{

		pt->core.R_H = pt->Disk.R_Disk_interp;
		c = (pt->Disk.R_Disk_interp / R_H_orig) * (pt->Disk.R_Disk_interp / R_H_orig);
	}

	set_Disk_angles(pt);
	I = integrale_simp_struct(pf, pt, pt->Disk.Disk_mu_1, pt->Disk.Disk_mu_2, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	set_Disk_angles(pt);
	return I * one_by_four_pi * c;

}

double eval_I_nu_Disk_disk_RF(struct blob *pt, double nu_disk_RF)
{
	double (*pf)(struct blob *, double x);
	double  I, R_H_orig, c;
	//unsigned int i;
	pt->core.nu_disk_RF = nu_disk_RF;
	pf = &integrand_I_nu_Disk_disk_RF;

	c = 1.0;
	R_H_orig = pt->core.R_H;
	if (pt->core.R_H > pt->Disk.R_Disk_interp)
	{

		pt->core.R_H = pt->Disk.R_Disk_interp;
		c = (pt->Disk.R_Disk_interp / R_H_orig) * (pt->Disk.R_Disk_interp / R_H_orig);
	}
	set_Disk_angles(pt);
	I = integrale_simp_struct(pf, pt, pt->Disk.Disk_mu_1, pt->Disk.Disk_mu_2, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	set_Disk_angles(pt);
	//printf("=> R_DT_interp=%e R_H=%e Disk_mu_1=%e Disk_mu_2=%e i=%e \n", pt->R_DT_interp, pt->core.R_H, pt->Disk.Disk_mu_1, pt->Disk.Disk_mu_2, I);
	return I * one_by_four_pi * c;
}

double eval_Disk_L_nu(struct blob *pt, double nu_Disk_disk_RF)
{
	if (pt->core.disk == 2) {
		//in this case no multiplication by L_Disk, because we acutally integrate every annluar BB along the disk
		//printf("=> %e\n", `(pt, nu_Disk_disk_RF));
		//printf("=> nu_Disk_disk_RF %e\n", nu_Disk_disk_RF);
		return  Disk_Spectrum(pt, nu_Disk_disk_RF);
	}
	else{
		return  pt->Disk.L_Disk *Disk_Spectrum(pt, nu_Disk_disk_RF);
	}
}

double eval_nu_peak_Disk(double T){
	return eval_nu_peak_planck(T);
}



//========================
// Disk Geometrical Functions
//========================

void set_Disk_angles(struct blob *pt)
{
	double mu1, mu2;
	mu1 = pt->core.R_H / sqrt(pt->core.R_H * pt->core.R_H + pt->Disk.R_inner * pt->Disk.R_inner);
	mu2 = pt->core.R_H / sqrt(pt->core.R_H * pt->core.R_H + pt->Disk.R_ext * pt->Disk.R_ext);
	//mu1=1.0/sqrt(1+((pt->Disk.R_inner*pt->Disk.R_inner)/(pt->core.R_H*pt->core.R_H)));
	//mu2 = 1.0 / sqrt(1 + ((pt->Disk.R_ext * pt->Disk.R_ext) / (pt->core.R_H * pt->core.R_H)));
	pt->Disk.Disk_mu_1 = min(mu1, mu2);
	pt->Disk.Disk_mu_2 = max(mu1, mu2);
}

void set_Disk_geometry(struct blob *pt){

	pt->Disk.Disk_surface=pi*((pt->Disk.R_ext * pt->Disk.R_ext) - (pt->Disk.R_inner*pt->Disk.R_inner) );
	pt->Disk.Disk_geom_factor = (1.0) / (four_pi * pt->core.R_H * pt->core.R_H * (pt->Disk.Disk_surface / (pt->core.R_H * pt->core.R_H)));
}




//=========================================================================================
void Build_I_nu_BLR(struct blob *pt){

	//-------------------------------------
	// we follow the method in Donea&Protheroe https://arxiv.org/abs/astro-ph/0202068v1
	//-------------------------------------
	//double nu_stop_BLR_blob_RF;
	unsigned int NU_INT,NU_INT_MAX;
	unsigned int ANGLE_INT, ANGLE_INT_MAX;
	int have_angle_storage;
	size_t angle_idx;
	double d_theta, theta_grid, mu_grid;
	double R_H_orig_angle, R_H_eval_angle, c_angle;
	double theta_max_angle;
	double I_theta_DRF, I_theta_blob, geom_theta;
	double I_nu_theta_disk_RF,I_nu_theta_blob_RF;
	//char f_BLR_disk[static_file_name_max_legth];
	//FILE *fp_BLR_disk;

	/*
	if (pt->core.WRITE_TO_FILE==1){
		sprintf(f_BLR_disk, "%s%s-I_nu_BLR.dat",pt->core.path, pt->core.STEM);

		fp_BLR_disk = fopen(f_BLR_disk, "w");
		if (fp_BLR_disk == NULL) {
			printf("unable to open %s\n ", fp_BLR_disk);
			exit(1);
		}
	}
	*/
	//flux_DISK_header(fp_BLR_disk);
	if (pt->core.verbose){

		printf("-----------  Building I_nu BLR     ----------- \n");
	}
	set_BLR_geometry(pt);
	//printf("=>R_H=%e BLR_mu_1=%e BLR_mu_2=%e\n",pt->core.R_H,pt->BLR.BLR_mu_1,pt->BLR.BLR_mu_2);

	pt->BLR.BLR_mu_1 = 1.0;
	pt->BLR.BLR_mu_2 = cos(eval_theta_max_BLR(pt));
	
	//if (pt->BLR.tau_BLR>0.9){
	//	printf ("!!! Waring, the fraction of L_Disk reaching DT is (1-tau_BLR)\n");
	//	printf ("!!! if tau_BLR=1.0 no DT photons will be generated\n");

	//}

	pt->BLR.spec.nu_min_DRF=pt->Disk.spec.nu_DRF[0];
	pt->BLR.spec.nu_max_DRF=pt->Disk.spec.nu_DRF[pt->Disk.spec.NU_INT_MAX];

	pt->BLR.spec.nu_min = eval_nu_max_blob_RF(pt, pt->BLR.BLR_mu_1, pt->BLR.BLR_mu_2, pt->BLR.spec.nu_min_DRF);
	pt->BLR.spec.nu_max  = eval_nu_max_blob_RF(pt,pt->BLR.BLR_mu_1, pt->BLR.BLR_mu_2, pt->BLR.spec.nu_max_DRF);

	pt->BLR.R_BLR_interp_val = pt->BLR.R_BLR_out * 50.0;
	pt->BLR.R_BLR_interp_start = pt->BLR.R_BLR_out * 50.0;
	//printf("=>R_H=%e BLR_mu_1=%e BLR_mu_2=%e nu1=%e nu2=%e  nu1 d=%e nu2 d=%e\n", pt->core.R_H, pt->BLR.BLR_mu_1, pt->BLR.BLR_mu_2, pt->BLR.spec.nu_min, pt->BLR.spec.nu_max, pt->BLR.spec.nu_min_DRF, pt->BLR.spec.nu_max_DRF);
	pt->BLR.n0_BLR = pt->BLR.tau_BLR / (SIGTH * (pt->BLR.R_BLR_out - pt->BLR.R_BLR_in));
	

	if (pt->core.verbose)
	{
		printf("BLR_mu_1=%e BLR_mu_2=%e\n", pt->BLR.BLR_mu_1, pt->BLR.BLR_mu_2);

		printf("n0_BLR=%e \n", pt->BLR.n0_BLR);

		printf("nu_start_BLR_disk_RF=%e  nu_stop_BLR_disk_RF=%e \n",
				   pt->BLR.spec.nu_min_DRF,
				   pt->BLR.spec.nu_max_DRF);

		printf("nu_start_BLR=%e  nu_stop_BLR=%e \n",
					pt->BLR.spec.nu_min,
					pt->BLR.spec.nu_max);
	}


	NU_INT_MAX = pt->core.nu_seed_size-1;
	pt->BLR.spec.NU_INT_MAX=NU_INT_MAX;
	have_angle_storage = 0;
	R_H_orig_angle = pt->core.R_H;
	R_H_eval_angle = R_H_orig_angle;
	c_angle = 1.0;
	if (ensure_external_spectrum_angle_dep(&(pt->BLR.spec), pt->core.nu_seed_size, pt->core.theta_n_int) == 0 &&
		pt->BLR.spec.angle_n_int > 0U) {
		have_angle_storage = 1;
		if (R_H_eval_angle > pt->BLR.R_BLR_interp_start) {
			R_H_eval_angle = pt->BLR.R_BLR_interp_val;
			c_angle = (pt->BLR.R_BLR_interp_val / R_H_orig_angle) * (pt->BLR.R_BLR_interp_val / R_H_orig_angle);
		}
		pt->core.R_H = R_H_eval_angle;
		theta_max_angle = eval_theta_max_BLR(pt);
		ANGLE_INT_MAX = pt->BLR.spec.angle_n_int - 1U;
		d_theta = (ANGLE_INT_MAX > 0U) ? (theta_max_angle / (double)ANGLE_INT_MAX) : 0.0;
		for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
			theta_grid = d_theta * (double)ANGLE_INT;
			mu_grid = cos(theta_grid);
			pt->BLR.spec.theta[ANGLE_INT] = theta_grid;
			pt->BLR.spec.mu[ANGLE_INT] = mu_grid;
		}
		pt->core.R_H = R_H_orig_angle;
	}
	//This is evaluating the angular pattern
	//It does not depends on frequency, because each region
	//is emitting the same spectrum
	I_nu_theta_disk_RF = eval_I_nu_BLR_disk_RF(pt);
	I_nu_theta_blob_RF = eval_I_nu_BLR_blob_RF(pt);

	build_log_grid( pt->BLR.spec.nu_min_DRF,  pt->BLR.spec.nu_max_DRF, pt->core.nu_seed_size, pt->BLR.spec.nu_DRF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		//pt->BLR.spec.L_nu_DRF[NU_INT] = eval_Lnu_BLR_disk_RF(pt, pt->BLR.spec.nu_DRF[NU_INT]);
		pt->BLR.spec.L_nu_DRF[NU_INT] = eval_Lnu_BLR_disk_RF(pt,pt->Disk.spec.L_nu_DRF[NU_INT]);
	}
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->BLR.spec.I_nu_DRF[NU_INT] = I_nu_theta_disk_RF * pt->BLR.spec.L_nu_DRF[NU_INT];
	}


	build_log_grid( pt->BLR.spec.nu_min,  pt->BLR.spec.nu_max, pt->core.nu_seed_size, pt->BLR.spec.nu);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		//we have to pass nu_BLR_disk_RF, because we integrate
		//the I' expressed in terms of I
		pt->BLR.spec.I_nu[NU_INT] = I_nu_theta_blob_RF * pt->BLR.spec.L_nu_DRF[NU_INT];
		pt->BLR.spec.n_nu[NU_INT] =I_nu_to_n(pt->BLR.spec.I_nu[NU_INT], pt->BLR.spec.nu[NU_INT]);
		//EC with n(gamma) transf
		pt->BLR.spec.n_nu_DRF[NU_INT] = I_nu_to_n(pt->BLR.spec.I_nu_DRF[NU_INT], pt->BLR.spec.nu_DRF[NU_INT]);

		if (pt->BLR.spec.I_nu[NU_INT]>pt->core.emiss_lim){
			pt->BLR.spec.nu_max = pt->BLR.spec.nu[NU_INT];
			pt->BLR.spec.NU_INT_MAX = NU_INT;
		}
		else{
			pt->BLR.spec.I_nu[NU_INT]=pt->core.emiss_lim;
			pt->BLR.spec.n_nu[NU_INT] =I_nu_to_n(pt->BLR.spec.I_nu[NU_INT], pt->BLR.spec.nu[NU_INT]);
		}

		if (pt->core.verbose>1){
			printf(" nu_BLR_disk_RF=%e, I_nu_BLR_disk_RF=%e, nu_BLR=%e, , I_nu_BLR=%e\n",
					pt->BLR.spec.nu_DRF[NU_INT],
					pt->BLR.spec.I_nu_DRF[NU_INT],
					pt->BLR.spec.nu[NU_INT],
					pt->BLR.spec.I_nu[NU_INT]);
		}
		/*
		if (pt->core.WRITE_TO_FILE==1){

			fprintf(fp_BLR_disk, "%4.4e\t %4.4e\t %4.4e\t %4.4e \n",
				log10(pt->BLR.spec.nu_DRF[NU_INT]),
				log10(pt->BLR.spec.I_nu_DRF[NU_INT]),
				log10(pt->BLR.spec.nu[NU_INT]),
				log10(pt->BLR.spec.I_nu[NU_INT]));
		}
		*/

	}
	if (have_angle_storage) {
		pt->core.R_H = R_H_eval_angle;
		for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++) {
			for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
				mu_grid = pt->BLR.spec.mu[ANGLE_INT];
				geom_theta = eval_I_nu_theta_BLR(pt, mu_grid);
				I_theta_DRF = c_angle * geom_theta * pt->BLR.spec.L_nu_DRF[NU_INT];
				I_theta_blob = I_theta_DRF * pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu_grid);
				angle_idx = angle_dep_flat_index(NU_INT, ANGLE_INT, pt->BLR.spec.angle_n_int);
				pt->BLR.spec.I_nu_theta_DRF[angle_idx] = I_theta_DRF;
				pt->BLR.spec.I_nu_theta[angle_idx] = I_theta_blob;
				pt->BLR.spec.n_nu_theta_DRF[angle_idx] = I_nu_to_n(I_theta_DRF, pt->BLR.spec.nu_DRF[NU_INT]);
				pt->BLR.spec.n_nu_theta[angle_idx] = I_nu_to_n(I_theta_blob, pt->BLR.spec.nu[NU_INT]);
			}
		}
		pt->core.R_H = R_H_orig_angle;
	}
	/*
	if (pt->core.WRITE_TO_FILE == 1)
	{
		fclose(fp_BLR_disk);
	}
	*/
}



//========================
// BLR Spectral Functions
//========================

//-------------------------------------
// we follow the method in Donea&Protheroe https://arxiv.org/abs/astro-ph/0202068v1
//-------------------------------------

double j_nu_BLR_integrand(struct blob *pt, double l)
{
	//unsigned int i;
	double L, r2;

	//i = x_to_grid_index(pt->BLR.spec.nu_DRF, pt->core.nu_disk_RF, pt->core.nu_seed_size);
	
	r2 = (pt->core.R_H * pt->core.R_H) - 2.0 * pt->core.R_H * l * pt->BLR.mu_j + l * l;
	
	
	//L = eval_Disk_L_nu(pt, pt->core.nu_disk_RF) * pt->BLR.n0_BLR * SIGTH;
	if ((r2 > (pt->BLR.R_BLR_out * pt->BLR.R_BLR_out)) || (r2 < (pt->BLR.R_BLR_in * pt->BLR.R_BLR_in)))
	{
		L=0.0;
	}
	else{
		L =1.0/ (four_pi * four_pi * r2);
	}
	return L;
}

double eval_I_nu_theta_BLR(struct blob *pt, double mu)
{
	double (*pf)(struct blob *, double x);
	//unsigned int i;
	double l_values[3], I;
	
	pf = &j_nu_BLR_integrand;
	pt->BLR.mu_j=mu;
	
	eval_l_values_BLR(pt, mu, l_values);
	if(pt->core.R_H<pt->BLR.R_BLR_out){
		I = integrale_simp_struct(pf, pt, 0, l_values[0], pt->core.l_n_int)+ integrale_simp_struct(pf, pt, l_values[1], l_values[2], pt->core.l_n_int);
	}
	else{
		I = integrale_simp_struct(pf, pt, 0, l_values[2], pt->core.l_n_int);
	}
	//printf("mu=%e, l0=%e, l1=%e, l2=%e, delta=%e\n", mu, l_values[0], l_values[1], l_values[2], l_values[2]- l_values[1]);
	return I;
}

double integrand_I_nu_BLR_blob_RF(struct blob *pt, double theta)
{
	//double psi
	//double mu,mu1,c;
	//mu = cos(theta);
	//mu1 = (pt->core.beta_Gamma - mu )/(pt->core.beta_Gamma*mu - 1.0 );
	//c=(pt->core.BulkFactor * pt->core.BulkFactor * pt->core.BulkFactor );
	double f;
	//c = c * (1.0 + pt->core.BulkFactor * mu + 1.0) * (1.0 + pt->core.BulkFactor * mu + 1.0) * (1.0 + pt->core.BulkFactor * mu + 1.0);
	f=pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * cos(theta));
	//f=1/( (pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * cos(theta))) * (pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * cos(theta))));
	return 2 * pi * sin(theta) * eval_I_nu_theta_BLR(pt, cos(theta)) *f;
}

double integrand_I_nu_BLR_disk_RF(struct blob * pt, double theta)
{
	//double psi;
	return 2 * pi * sin(theta) * eval_I_nu_theta_BLR( pt,  cos(theta));
}

double eval_I_nu_BLR_disk_RF(struct blob *pt)
{
	double (*pf)(struct blob *, double x);
	double theta_min, theta_max, I, R_H_orig,c;

	//pt->core.nu_disk_RF=nu_disk_RF;
	pf = &integrand_I_nu_BLR_disk_RF;
	c=1.0;
	R_H_orig = pt->core.R_H;
	if (pt->core.R_H > pt->BLR.R_BLR_interp_start)
	{
		
		pt->core.R_H = pt->BLR.R_BLR_interp_val;
		c = (pt->BLR.R_BLR_interp_val / R_H_orig) * (pt->BLR.R_BLR_interp_val / R_H_orig);
		//printf("=>R_H=%e R_H_orig=%e  pt->R_BLR_interp=%e\n",pt->core.R_H,R_H_orig,pt->R_BLR_interp);
	}
	theta_min=0.0;
	theta_max = eval_theta_max_BLR(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	//printf("=>R_H=%e R_BLR_inter=%e I=%e %e %e c=%e\n ",pt->core.R_H,pt->R_BLR_interp, I, theta_min, theta_max,c);
	return I*one_by_four_pi*c;
}


double eval_I_nu_BLR_blob_RF(struct blob *pt)
{
	double (*pf)(struct blob *, double x);
	double theta_min, theta_max, I, R_H_orig,c;
	// we use directly nu_disk_RF
	// because we integrate the I' expressed as I
	//pt->core.nu_disk_RF = nu_disk_RF;
	pf = &integrand_I_nu_BLR_blob_RF;
	
	c=1.0;
	R_H_orig = pt->core.R_H;
	if (pt->core.R_H > pt->BLR.R_BLR_interp_start)
	{

		pt->core.R_H = pt->BLR.R_BLR_interp_val;
		c = (pt->BLR.R_BLR_interp_val / R_H_orig) * (pt->BLR.R_BLR_interp_val / R_H_orig);
		//printf("=>R_H=%e R_H_orig=%e  pt->R_BLR_interp=%e\n",pt->core.R_H,R_H_orig,pt->R_BLR_interp);
	}
	theta_min = 0.0;
	theta_max = eval_theta_max_BLR(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	//printf("=>BLR R_H=%e R_B=%e I=%e %e %e c=%e\n ", pt->core.R_H, pt->BLR.R_BLR_out, I, theta_min, theta_max, c);
	return I*one_by_four_pi*c;
}

//double eval_Lnu_BLR_disk_RF(struct blob *pt, double nu_disk_RF)
double eval_Lnu_BLR_disk_RF(struct blob *pt, double Disk_L_nu)
{
	return Disk_L_nu* pt->BLR.n0_BLR *SIGTH;
	//return eval_Disk_L_nu(pt, nu_disk_RF) * pt->BLR.n0_BLR *SIGTH;
}




//========================
// BLR Geometrical Functions
//========================

double eval_theta_max_BLR(struct blob *pt)
{
	double theta_max;
	if (pt->core.R_H > pt->BLR.R_BLR_out)
	{
		theta_max = asin(pt->BLR.R_BLR_out / pt->core.R_H);
		//theta_max = 2 * (pi * 0.5 - acos(pt->BLR.R_BLR_out / pt->core.R_H))	;
	}
	else
	{
		theta_max = pi;
	}

	
	return theta_max;
	
}

void eval_l_values_BLR(struct blob *pt, double mu, double l[])
{
	double s;

		
	s = mu * mu + (pt->BLR.R_BLR_in / pt->core.R_H) * (pt->BLR.R_BLR_in / pt->core.R_H) - 1.0;
	if (s < 0.0){
		l[0] = 0.0;
		l[1] = 0.0;
	}else
	{
		l[1] = pt->core.R_H * mu + pt->core.R_H * sqrt(s);
		l[0]= pt->core.R_H * mu - pt->core.R_H * sqrt(s);
	}
	if (l[1] < 0.0){
		l[1] = 0.;
	}

	if (l[0] < 0.0){
		l[0] = 0.;
	}

	s = mu * mu + (pt->BLR.R_BLR_out / pt->core.R_H) * (pt->BLR.R_BLR_out / pt->core.R_H) - 1.0;
	if (s < 0.0){
		l[2] = 0;
	}
	else{

		l[2] = pt->core.R_H * mu + pt->core.R_H * sqrt(s);
		
		if (l[2] < 0.0){
			l[2] = 0.;
		}
	}
}

void set_BLR_geometry(struct blob *pt)
{

	pt->BLR.BLR_Volume = (4. / 3.) * pi * ((pt->BLR.R_BLR_out * pt->BLR.R_BLR_out * pt->BLR.R_BLR_out) - (pt->BLR.R_BLR_in * pt->BLR.R_BLR_in * pt->BLR.R_BLR_in));
	pt->BLR.BLR_inner_Surface = 4 * pi * (pt->BLR.R_BLR_in * pt->BLR.R_BLR_in);
	pt->BLR.Delta_R_BLR = pt->BLR.R_BLR_out - pt->BLR.R_BLR_in;

	

	/*
	if (pt->core.R_H < pt->BLR.R_BLR_in)
	{
		pt->BLR.BLR_mu_1 = -1.0;
		pt->BLR.BLR_mu_2 = pt->core.R_H / sqrt(pt->core.R_H * pt->core.R_H + pt->Disk.R_ext * pt->Disk.R_ext);
		//pt->BLR_mu_r_J_2 = pt->core.R_H / sqrt(pt->core.R_H * pt->core.R_H + pt->Disk.R_ext * pt->Disk.R_ext);
		//pt->BLR_mu_r_J_1 = -1.0;
		//pt->BLR_geom_factor=(1.0)/(4*pi*4*pi);
		//pt->BLR_geom_factor*=1.0/(pt->BLR.R_BLR_in*pt->BLR.R_BLR_in);
	}
	else if (pt->core.R_H >= pt->BLR.R_BLR_in && pt->core.R_H < pt->BLR.R_BLR_out)
	{
		pt->BLR.BLR_mu_1 = -1.0;
		pt->BLR.BLR_mu_2 = pt->core.R_H / sqrt(pt->core.R_H * pt->core.R_H + pt->Disk.R_ext * pt->Disk.R_ext);
		//pt->BLR_mu_r_in = sqrt(pt->core.R_H * pt->core.R_H - pt->BLR.R_BLR_in * pt->BLR.R_BLR_in) / pt->core.R_H;
		//pt->BLR_mu_r_J_2 = pt->core.R_H / sqrt(pt->core.R_H * pt->core.R_H + pt->Disk.R_ext * pt->Disk.R_ext);
		//pt->BLR_mu_r_J_1 = -1.0;
	}
	else
	{
		pt->BLR.BLR_mu_2 = pt->core.R_H / sqrt(pt->core.R_H * pt->core.R_H + pt->Disk.R_ext * pt->Disk.R_ext);
		pt->BLR.BLR_mu_1 = sqrt(pt->core.R_H * pt->core.R_H - pt->BLR.R_BLR_out * pt->BLR.R_BLR_out) / pt->core.R_H;
		//pt->BLR_mu_r_in = sqrt(pt->core.R_H * pt->core.R_H - pt->BLR.R_BLR_in * pt->BLR.R_BLR_in) / pt->core.R_H;
		//pt->BLR_mu_r_J_1 = pt->BLR.BLR_mu_1;
		//pt->BLR_mu_r_J_2 = pt->BLR.BLR_mu_2;
		//pt->BLR_geom_factor=(1.0)/(4*pi*4*pi);
		//pt->BLR_geom_factor*=1.0/(pt->BLR.R_BLR_in*pt->BLR.R_BLR_in);
	}
	*/
}

//=========================================================================================

void Build_I_nu_DT(struct blob *pt){
	//FILE *fp_SED_DT;
	//char f_SED_DT[static_file_name_max_legth];
	unsigned int NU_INT,NU_INT_MAX;
	unsigned int ANGLE_INT, ANGLE_INT_MAX;
	int have_angle_storage;
	size_t angle_idx;
	double d_theta, theta_grid, mu_grid;
	double R_H_orig_angle, R_H_eval_angle, c_angle;
	double theta_max_angle;
	double I_theta_DRF, I_theta_blob, geom_theta;
	double I_nu_theta_disk_RF, I_nu_theta_blob_RF;
	double nu_peak_DT_disk_RF;
	double nu_start_DT_disk_RF,nu_stop_DT_disk_RF;
	double nu_obs;
	double nuL_nu_DT,F_nu_DT_obs;

	/*
	if (pt->core.WRITE_TO_FILE==1){
		sprintf(f_SED_DT, "%s%s-SED-DT.dat",
					pt->core.path, pt->core.STEM);

		fp_SED_DT = fopen(f_SED_DT, "w");
		if (fp_SED_DT == NULL) {
			printf("unable to open %s\n ", f_SED_DT);
			exit(1);
		}
		flux_DISK_header(fp_SED_DT);
	}
	*/

	if (pt->core.verbose){

		printf("-----------  Building I_nu DT     ----------- \n");
	}


	//if (pt->BLR.tau_BLR>0.9){
	//		printf ("!!! Waring, the fraction of L_Disk reaching DT is (1-tau_BLR)\n");
	//		printf ("!!! if tau_BLR=1.0 no DT photons will be generated\n");

	//}

	nu_peak_DT_disk_RF=eval_nu_peak_planck(pt->DT.T_DT);

	nu_start_DT_disk_RF=nu_peak_DT_disk_RF*pt->core.nu_planck_min_factor;
	nu_stop_DT_disk_RF=nu_peak_DT_disk_RF*pt->core.nu_planck_max_factor;


	pt->DT.DT_mu_1 = 1.0;
	pt->DT.DT_mu_2 = cos(eval_theta_max_DT(pt));

	pt->DT.spec.nu_min = eval_nu_max_blob_RF(pt, pt->DT.DT_mu_1, pt->DT.DT_mu_2, nu_start_DT_disk_RF);
	pt->DT.spec.nu_max  = eval_nu_max_blob_RF(pt,pt->DT.DT_mu_1, pt->DT.DT_mu_2, nu_stop_DT_disk_RF);

	pt->DT.spec.nu_min_DRF=nu_start_DT_disk_RF;
	pt->DT.spec.nu_max_DRF=nu_stop_DT_disk_RF;

	pt->DT.spec.nu_min_obs=nu_disk_to_nu_obs_disk(nu_start_DT_disk_RF , pt->core.z_cosm);
	pt->DT.spec.nu_max_obs=nu_disk_to_nu_obs_disk(nu_stop_DT_disk_RF, pt->core.z_cosm);

	if (pt->core.verbose){
		printf("nu_start_DT (blob frame) =%e \n",
					pt->DT.spec.nu_min);
		printf("nu_stop_DT (blob frame) =%e \n",
						pt->DT.spec.nu_max);
		printf("nu_start_DT (disk frame) =%e \n",
				nu_start_DT_disk_RF);
		printf("nu_stop_DT (disk frame) =%e \n",
				nu_stop_DT_disk_RF);
	}

	NU_INT_MAX = pt->core.nu_seed_size-1;
	pt->DT.spec.NU_INT_MAX = NU_INT_MAX;
	have_angle_storage = 0;
	R_H_orig_angle = pt->core.R_H;
	R_H_eval_angle = R_H_orig_angle;
	c_angle = 1.0;
	if (ensure_external_spectrum_angle_dep(&(pt->DT.spec), pt->core.nu_seed_size, pt->core.theta_n_int) == 0 &&
		pt->DT.spec.angle_n_int > 0U) {
		have_angle_storage = 1;
		if (R_H_eval_angle > (pt->DT.R_DT * 50.0)) {
			R_H_eval_angle = pt->DT.R_DT * 50.0;
			c_angle = (R_H_eval_angle / R_H_orig_angle) * (R_H_eval_angle / R_H_orig_angle);
		}
		pt->core.R_H = R_H_eval_angle;
		theta_max_angle = eval_theta_max_DT(pt);
		ANGLE_INT_MAX = pt->DT.spec.angle_n_int - 1U;
		d_theta = (ANGLE_INT_MAX > 0U) ? (theta_max_angle / (double)ANGLE_INT_MAX) : 0.0;
		for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
			theta_grid = d_theta * (double)ANGLE_INT;
			mu_grid = cos(theta_grid);
			pt->DT.spec.theta[ANGLE_INT] = theta_grid;
			pt->DT.spec.mu[ANGLE_INT] = mu_grid;
		}
		pt->core.R_H = R_H_orig_angle;
	}

	pt->DT.R_DT_interp_val = pt->DT.R_DT* 50.0;
	pt->DT.R_DT_interp_start = pt->DT.R_DT * 50.0;

	pt->DT.DT_Volume=(4./3.)*pi*pt->DT.R_DT*pt->DT.R_DT*pt->DT.R_DT;

	//This is evaluating the angular pattern
	//It does not depends on frequency, because each region
	//is emitting the same spectrum
	I_nu_theta_disk_RF = eval_I_nu_DT_disk_RF(pt);
	I_nu_theta_blob_RF = eval_I_nu_DT_blob_RF(pt);

	build_log_grid( nu_start_DT_disk_RF,  nu_stop_DT_disk_RF, pt->core.nu_seed_size, pt->DT.spec.nu_DRF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->DT.spec.L_nu_DRF[NU_INT] = eval_DT_L_nu(pt, pt->DT.spec.nu_DRF[NU_INT]);
	}
	for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++){
		pt->DT.spec.I_nu_DRF[NU_INT] = I_nu_theta_disk_RF * pt->DT.spec.L_nu_DRF[NU_INT];
	}


	build_log_grid( pt->DT.spec.nu_min,  pt->DT.spec.nu_max, pt->core.nu_seed_size, pt->DT.spec.nu);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {


		nu_obs = nu_disk_to_nu_obs_disk(pt->DT.spec.nu_DRF[NU_INT], pt->core.z_cosm);

		pt->DT.spec.nu_obs[NU_INT]=nu_obs;

		pt->DT.spec.I_nu[NU_INT] = I_nu_theta_blob_RF * pt->DT.spec.L_nu_DRF[NU_INT];
		pt->DT.spec.n_nu[NU_INT] =I_nu_to_n(pt->DT.spec.I_nu[NU_INT], pt->DT.spec.nu[NU_INT]);
		//EC with n(gamma) transf
		pt->DT.spec.n_nu_DRF[NU_INT] = I_nu_to_n(pt->DT.spec.I_nu_DRF[NU_INT], pt->DT.spec.nu_DRF[NU_INT]);

		if (pt->DT.spec.I_nu[NU_INT]>pt->core.emiss_lim){
			pt->DT.spec.nu_max = pt->DT.spec.nu[NU_INT];
			pt->DT.spec.NU_INT_MAX = NU_INT;
		}
		else{
			pt->DT.spec.I_nu[NU_INT]=pt->core.emiss_lim;
			pt->DT.spec.n_nu[NU_INT] =I_nu_to_n(pt->DT.spec.I_nu[NU_INT], pt->DT.spec.nu[NU_INT]);
		}

		nuL_nu_DT = pt->DT.spec.L_nu_DRF[NU_INT] * pt->DT.spec.nu_DRF[NU_INT];
		//if (pt->BLR.tau_BLR<1.0){
		//	nuL_nu_DT = pt->DT.spec.L_nu_DRF[NU_INT] * pt->DT.spec.nu_DRF[NU_INT];
		//}
		//else {
		//	nuL_nu_DT = 0;
		//}
		F_nu_DT_obs = L_nu_Disk_to_F_nu(nuL_nu_DT / pt->DT.spec.nu_DRF[NU_INT], pt->core.z_cosm, pt->core.dist);

		pt->DT.spec.nuFnu_obs[NU_INT] = F_nu_DT_obs*nu_obs;
		if (pt->core.verbose>1){
			printf(" nu_DT_disk_RF=%e, I_nu_DT_disk_RF=%e, nu_DT=%e, I_nu_DT=%e\n",
				pt->DT.spec.nu_DRF[NU_INT],
				pt->DT.spec.I_nu_DRF[NU_INT],
				pt->DT.spec.nu[NU_INT],
				pt->DT.spec.I_nu[NU_INT]);
		}
		/*
		if (pt->core.WRITE_TO_FILE==1){
			fprintf(fp_SED_DT, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
									log10(nu_obs),
									log10(nu_obs * F_nu_DT_obs),
									nu_obs,
									nu_obs*F_nu_DT_obs,
									pt->DT.spec.nu_DRF[NU_INT],
									nuL_nu_DT);
		}
		*/
	}
	/*
	if (pt->core.WRITE_TO_FILE==1){
		fclose(fp_SED_DT);
	}
	*/
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->DT.spec.n_nu[NU_INT] =I_nu_to_n(pt->DT.spec.I_nu[NU_INT], pt->DT.spec.nu[NU_INT]);
	}
	if (have_angle_storage) {
		pt->core.R_H = R_H_eval_angle;
		for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++) {
			for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
				mu_grid = pt->DT.spec.mu[ANGLE_INT];
				theta_grid = pt->DT.spec.theta[ANGLE_INT];
				geom_theta = eval_I_nu_theta_DT(pt, mu_grid, theta_grid);
				I_theta_DRF = c_angle * geom_theta * pt->DT.spec.L_nu_DRF[NU_INT];
				I_theta_blob = I_theta_DRF * pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu_grid);
				angle_idx = angle_dep_flat_index(NU_INT, ANGLE_INT, pt->DT.spec.angle_n_int);
				pt->DT.spec.I_nu_theta_DRF[angle_idx] = I_theta_DRF;
				pt->DT.spec.I_nu_theta[angle_idx] = I_theta_blob;
				pt->DT.spec.n_nu_theta_DRF[angle_idx] = I_nu_to_n(I_theta_DRF, pt->DT.spec.nu_DRF[NU_INT]);
				pt->DT.spec.n_nu_theta[angle_idx] = I_nu_to_n(I_theta_blob, pt->DT.spec.nu[NU_INT]);
			}
		}
		pt->core.R_H = R_H_orig_angle;
	}
}

//========================
// Torus Spectral Functions
//========================

double j_nu_DT_integrand(struct blob *pt, double l)
{
	//unsigned int i;
	double L, r2;

	//i = x_to_grid_index(pt->BLR.spec.nu_DRF, pt->core.nu_disk_RF, pt->core.nu_seed_size);
	
	r2 = (pt->core.R_H * pt->core.R_H) - 2.0 * pt->DT.R_DT * l * pt->BLR.mu_j + l * l;
	
	
	//L = eval_Disk_L_nu(pt, pt->core.nu_disk_RF) * pt->BLR.n0_BLR * SIGTH;
	if  (r2 > (pt->DT.R_DT * pt->DT.R_DT) )
	{
		L=0.0;
	}
	else{
		L =1.0/ (four_pi * four_pi * r2*pt->DT.R_DT);
	}
	return L;
}


double eval_I_nu_theta_DT(struct blob *pt, double mu, double theta)
{
	//double (*pf)(struct spettro *, double x);
	//unsigned int i;
	double l, I,cos_theta_norm,alpha;
	//double I;
	
	
	//i = x_to_grid_index(pt->DT.spec.nu_DRF, pt->core.nu_disk_RF, pt->core.nu_seed_size);
	if (pt->core.R_H < pt->DT.R_DT)
	{
		I = 1.0 / (4 * pi * 4 * pi * pt->DT.R_DT * pt->DT.R_DT);
	}	
	else
	{
		l = eval_l_DT(pt, mu);
		//pf = &j_nu_DT_integrand;
		//pt->BLR.mu_j = mu;

		
		alpha = acos(l * sin(theta) / pt->DT.R_DT);

		cos_theta_norm = cos(pi - (alpha + 0.5 * pi - theta ));

		I = cos_theta_norm	 / ((4 * pi * pi * pt->core.R_H * pt->core.R_H) * ((pt->DT.R_DT / pt->core.R_H) * (pt->DT.R_DT / pt->core.R_H)));

		//I = integrale_simp_struct(pf, pt, 0, l, pt->core.l_n_int);

	}
			
	//I = integrale_simp_struct(pf, pt, 0, l, pt->core.l_n_int);
	return I;
}

double integrand_I_nu_DT_blob_RF(struct blob *pt, double theta)
{
	//double psi;
	double f;
	//f=1/( (pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * cos(theta))) * (pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * cos(theta))));
	f=pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * cos(theta));
	return 2 * pi * sin(theta) * eval_I_nu_theta_DT(pt, cos(theta), theta)*f;
}

double integrand_I_nu_DT_disk_RF(struct blob *pt, double theta)
{
	//double psi;
	return 2 * pi * sin(theta) * eval_I_nu_theta_DT(pt, cos(theta), theta);
}

double eval_I_nu_DT_disk_RF(struct blob *pt )
{
	double (*pf)(struct blob *, double x);
	double theta_min, theta_max, I, R_H_orig, c;

	
	//now integrating only over angles
	//not needed anymore
	//pt->core.nu_disk_RF = nu_disk_RF;

	pf = &integrand_I_nu_DT_disk_RF;

	
	c = 1.0;
	R_H_orig = pt->core.R_H;
	if (pt->core.R_H > pt->DT.R_DT_interp_start)
	{

		pt->core.R_H = pt->DT.R_DT_interp_val;
		c = (pt->DT.R_DT_interp_val / R_H_orig) * (pt->DT.R_DT_interp_val / R_H_orig);
	}
	theta_min = 0.0;
	theta_max = eval_theta_max_DT(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	return I * one_by_four_pi * c;
}

double eval_I_nu_DT_blob_RF(struct blob *pt )
{
	double (*pf)(struct blob *, double x);
	double theta_min, theta_max, I, R_H_orig, c;

	//now integrating only over angles
	//not needed anymore
	//pt->core.nu_disk_RF = nu_disk_RF;

	pf = &integrand_I_nu_DT_blob_RF;

	c=1.0;
	R_H_orig = pt->core.R_H;
	if (pt->core.R_H > pt->DT.R_DT_interp_start)
	{

		pt->core.R_H = pt->DT.R_DT_interp_val;
		c = (pt->DT.R_DT_interp_val / R_H_orig) * (pt->DT.R_DT_interp_val / R_H_orig);
	}
	theta_min = 0.0;
	theta_max = eval_theta_max_DT(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	//printf("=>DT  R_H=%e R_D=%e I=%e %e %e c=%e\n ", pt->core.R_H, pt->DT.R_DT, I, theta_min, theta_max, c);
	return I * one_by_four_pi * c;
}

double eval_DT_L_nu(struct blob *pt, double DT_disk_RF)
{
	return pt->Disk.L_Disk_radiative * pt->DT.tau_DT * f_planck_norm(pt->DT.T_DT, DT_disk_RF);
}

//========================
// Torus Geometrical Functions
//========================
double eval_theta_max_DT(struct blob *pt)
{
	double theta_max;

	if (pt->core.R_H >= pt->DT.R_DT){
		theta_max= asin(pt->DT.R_DT / pt->core.R_H);
	}
	else{
		theta_max=pi;
	}

	return theta_max;
}

double eval_l_DT(struct blob *pt, double mu)
{
	double s,l;

	s = mu * mu + (pt->DT.R_DT / pt->core.R_H) * (pt->DT.R_DT / pt->core.R_H) - 1.0;
	if (s < 0.0){
		l=0.0;
	}
	else{
		l= pt->core.R_H * mu - pt->core.R_H * sqrt(s);
	}
	if (l < 0.0){
		l = 0.;
	}
	return l;
}

//=========================================================================================

static double integrand_f_nu_Corona_norm(struct blob *pt, double nu_Corona_disk_RF)
{
	return f_nu_Corona(pt, nu_Corona_disk_RF);
}

static double eval_dist_blob_corona(struct blob *pt)
{
	return fabs(pt->core.R_H - pt->Corona.R_H_Corona);
}

void Build_I_nu_Corona(struct blob *pt)
{
	unsigned int NU_INT, NU_INT_MAX;
	unsigned int ANGLE_INT, ANGLE_INT_MAX;
	int have_angle_storage;
	size_t angle_idx;
	double nu_start_Corona_disk_RF, nu_stop_Corona_disk_RF;
	double nu_ref_low, nu_ref_high;
	double nu_obs;
	double nuL_nu_Corona, F_nu_Corona_obs;
	double norm_int;
	double d_mu, mu_grid;
	double R_H_orig_angle, R_H_eval_angle, c_angle, R_H_test_angle;
	double dist_blob_corona;
	double I_theta_DRF, I_theta_blob;
	double (*pf)(struct blob *, double x);

	if (pt->core.verbose){
		printf("-----------  Building I_nu Corona     ----------- \n");
	}

	set_Corona_geometry(pt);
	set_Corona_angles(pt);

	nu_ref_high = pt->Corona.nu_cut_Corona;
	if (nu_ref_high <= 0.0){
		nu_ref_high = 1.0;
	}
	nu_ref_low = nu_ref_high;
	if (pt->Corona.nu_cut_low_Corona > 0.0){
		nu_ref_low = min(nu_ref_low, pt->Corona.nu_cut_low_Corona);
		nu_ref_high = max(nu_ref_high, pt->Corona.nu_cut_low_Corona);
	}
	nu_start_Corona_disk_RF = nu_ref_low * pt->core.nu_planck_min_factor;
	nu_stop_Corona_disk_RF = nu_ref_high * pt->core.nu_planck_max_factor;
	if (nu_start_Corona_disk_RF <= 0.0){
		nu_start_Corona_disk_RF = 1.0;
	}
	if (nu_stop_Corona_disk_RF <= nu_start_Corona_disk_RF){
		nu_stop_Corona_disk_RF = nu_start_Corona_disk_RF * 10.0;
	}

	pt->Corona.spec.nu_min = eval_nu_min_blob_RF(pt, pt->Corona.Corona_mu_1, pt->Corona.Corona_mu_2, nu_start_Corona_disk_RF);
	pt->Corona.spec.nu_max = eval_nu_max_blob_RF(pt, pt->Corona.Corona_mu_1, pt->Corona.Corona_mu_2, nu_stop_Corona_disk_RF);
	pt->Corona.spec.nu_min_DRF = nu_start_Corona_disk_RF;
	pt->Corona.spec.nu_max_DRF = nu_stop_Corona_disk_RF;
	pt->Corona.spec.nu_min_obs = nu_disk_to_nu_obs_disk(nu_start_Corona_disk_RF, pt->core.z_cosm);
	pt->Corona.spec.nu_max_obs = nu_disk_to_nu_obs_disk(nu_stop_Corona_disk_RF, pt->core.z_cosm);

	NU_INT_MAX = pt->core.nu_seed_size - 1;
	pt->Corona.spec.NU_INT_MAX = NU_INT_MAX;
	pt->Corona.R_Corona_interp_val = pt->Corona.R_Corona * 50.0;
	pt->Corona.R_Corona_interp_start = pt->Corona.R_Corona * 50.0;
	have_angle_storage = 0;
	R_H_orig_angle = pt->core.R_H;
	R_H_eval_angle = R_H_orig_angle;
	c_angle = 1.0;
	if (ensure_external_spectrum_angle_dep(&(pt->Corona.spec), pt->core.nu_seed_size, pt->core.theta_n_int) == 0 &&
		pt->Corona.spec.angle_n_int > 0U) {
		have_angle_storage = 1;
		dist_blob_corona = fabs(R_H_orig_angle - pt->Corona.R_H_Corona);
		if ((dist_blob_corona > pt->Corona.R_Corona_interp_start) && (dist_blob_corona > 0.0)) {
			if (R_H_orig_angle >= pt->Corona.R_H_Corona) {
				R_H_test_angle = pt->Corona.R_H_Corona + pt->Corona.R_Corona_interp_val;
			}
			else{
				R_H_test_angle = pt->Corona.R_H_Corona - pt->Corona.R_Corona_interp_val;
			}
			R_H_eval_angle = max(R_H_test_angle, 0.0);
			c_angle = (pt->Corona.R_Corona_interp_val / dist_blob_corona) * (pt->Corona.R_Corona_interp_val / dist_blob_corona);
		}

		pt->core.R_H = R_H_eval_angle;
		set_Corona_angles(pt);
		ANGLE_INT_MAX = pt->Corona.spec.angle_n_int - 1U;
		d_mu = (ANGLE_INT_MAX > 0U) ? ((pt->Corona.Corona_mu_2 - pt->Corona.Corona_mu_1) / (double)ANGLE_INT_MAX) : 0.0;
		for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
			mu_grid = pt->Corona.Corona_mu_1 + d_mu * (double)ANGLE_INT;
			if (mu_grid < -1.0) {
				mu_grid = -1.0;
			}
			else if (mu_grid > 1.0) {
				mu_grid = 1.0;
			}
			pt->Corona.spec.mu[ANGLE_INT] = mu_grid;
			pt->Corona.spec.theta[ANGLE_INT] = acos(mu_grid);
		}
		pt->core.R_H = R_H_orig_angle;
		set_Corona_angles(pt);
	}

	pf = &integrand_f_nu_Corona_norm;
	norm_int = integrale_trap_log_struct(pf, pt, nu_start_Corona_disk_RF, nu_stop_Corona_disk_RF, 400);
	if (norm_int > 0.0){
		pt->Corona.f_Corona_norm = 1.0 / norm_int;
	}
	else{
		pt->Corona.f_Corona_norm = 0.0;
	}

	build_log_grid(nu_start_Corona_disk_RF, nu_stop_Corona_disk_RF, pt->core.nu_seed_size, pt->Corona.spec.nu_DRF);
	for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++) {
		pt->Corona.spec.L_nu_DRF[NU_INT] = eval_Corona_L_nu(pt, pt->Corona.spec.nu_DRF[NU_INT]);
		pt->Corona.spec.I_nu_DRF[NU_INT] = eval_I_nu_Corona_disk_RF(pt, pt->Corona.spec.nu_DRF[NU_INT]);
	}

	build_log_grid(pt->Corona.spec.nu_min, pt->Corona.spec.nu_max, pt->core.nu_seed_size, pt->Corona.spec.nu);
	for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++) {
		nu_obs = nu_disk_to_nu_obs_disk(pt->Corona.spec.nu_DRF[NU_INT], pt->core.z_cosm);
		pt->Corona.spec.nu_obs[NU_INT] = nu_obs;
		pt->Corona.spec.I_nu[NU_INT] = eval_I_nu_Corona_blob_RF(pt, pt->Corona.spec.nu_DRF[NU_INT]);
		pt->Corona.spec.n_nu[NU_INT] = I_nu_to_n(pt->Corona.spec.I_nu[NU_INT], pt->Corona.spec.nu[NU_INT]);
		pt->Corona.spec.n_nu_DRF[NU_INT] = I_nu_to_n(pt->Corona.spec.I_nu_DRF[NU_INT], pt->Corona.spec.nu_DRF[NU_INT]);

		if (pt->Corona.spec.I_nu[NU_INT] > pt->core.emiss_lim){
			pt->Corona.spec.nu_max = pt->Corona.spec.nu[NU_INT];
			pt->Corona.spec.NU_INT_MAX = NU_INT;
		}
		else{
			pt->Corona.spec.I_nu[NU_INT] = pt->core.emiss_lim;
			pt->Corona.spec.n_nu[NU_INT] = I_nu_to_n(pt->Corona.spec.I_nu[NU_INT], pt->Corona.spec.nu[NU_INT]);
		}

		nuL_nu_Corona = pt->Corona.spec.L_nu_DRF[NU_INT] * pt->Corona.spec.nu_DRF[NU_INT];
		F_nu_Corona_obs = L_nu_Disk_to_F_nu(nuL_nu_Corona / pt->Corona.spec.nu_DRF[NU_INT], pt->core.z_cosm, pt->core.dist);
		pt->Corona.spec.nuFnu_obs[NU_INT] = F_nu_Corona_obs * nu_obs;
	}
	if (have_angle_storage) {
		pt->core.R_H = R_H_eval_angle;
		set_Corona_angles(pt);
		for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++) {
			pt->core.nu_disk_RF = pt->Corona.spec.nu_DRF[NU_INT];
			for (ANGLE_INT = 0; ANGLE_INT <= ANGLE_INT_MAX; ANGLE_INT++) {
				mu_grid = pt->Corona.spec.mu[ANGLE_INT];
				I_theta_DRF = c_angle * eval_I_nu_theta_Corona(pt, mu_grid);
				I_theta_blob = I_theta_DRF * pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu_grid);
				angle_idx = angle_dep_flat_index(NU_INT, ANGLE_INT, pt->Corona.spec.angle_n_int);
				pt->Corona.spec.I_nu_theta_DRF[angle_idx] = I_theta_DRF;
				pt->Corona.spec.I_nu_theta[angle_idx] = I_theta_blob;
				pt->Corona.spec.n_nu_theta_DRF[angle_idx] = I_nu_to_n(I_theta_DRF, pt->Corona.spec.nu_DRF[NU_INT]);
				pt->Corona.spec.n_nu_theta[angle_idx] = I_nu_to_n(I_theta_blob, pt->Corona.spec.nu[NU_INT]);
			}
		}
		pt->core.R_H = R_H_orig_angle;
		set_Corona_angles(pt);
	}
}

double f_nu_Corona(struct blob *pt, double nu_Corona_disk_RF)
{
	double f;
	if (nu_Corona_disk_RF <= 0.0 || pt->Corona.nu_cut_Corona <= 0.0){
		return 0.0;
	}
	f = pow(nu_Corona_disk_RF, -pt->Corona.alpha_Corona) * exp(-nu_Corona_disk_RF / pt->Corona.nu_cut_Corona);
	if (pt->Corona.nu_cut_low_Corona > 0.0){
		f *= exp(-pt->Corona.nu_cut_low_Corona / nu_Corona_disk_RF);
	}
	return f;
}

double eval_Corona_L_nu(struct blob *pt, double nu_Corona_disk_RF)
{
	return pt->Corona.L_Corona * pt->Corona.f_Corona_norm * f_nu_Corona(pt, nu_Corona_disk_RF);
}

double eval_I_nu_theta_Corona(struct blob *pt, double mu)
{
	(void)mu;
	return eval_Corona_L_nu(pt, pt->core.nu_disk_RF) * pt->Corona.Corona_geom_factor;
}

double integrand_I_nu_Corona_blob_RF(struct blob *pt, double mu)
{
	double f;
	f = pt->core.BulkFactor * (1.0 - pt->core.beta_Gamma * mu);
	return 2 * pi * eval_I_nu_theta_Corona(pt, mu) * f;
}

double integrand_I_nu_Corona_disk_RF(struct blob *pt, double mu)
{
	return 2 * pi * eval_I_nu_theta_Corona(pt, mu);
}

double eval_I_nu_Corona_disk_RF(struct blob *pt, double nu_Corona_disk_RF)
{
	double (*pf)(struct blob *, double x);
	double I, R_H_orig, dist_blob_corona, c, R_H_test;

	pt->core.nu_disk_RF = nu_Corona_disk_RF;
	pf = &integrand_I_nu_Corona_disk_RF;

	c = 1.0;
	R_H_orig = pt->core.R_H;
	dist_blob_corona = eval_dist_blob_corona(pt);
	if (dist_blob_corona > pt->Corona.R_Corona_interp_start && dist_blob_corona > 0.0)
	{
		if (pt->core.R_H >= pt->Corona.R_H_Corona){
			R_H_test = pt->Corona.R_H_Corona + pt->Corona.R_Corona_interp_val;
		}
		else{
			R_H_test = pt->Corona.R_H_Corona - pt->Corona.R_Corona_interp_val;
		}
		pt->core.R_H = max(R_H_test, 0.0);
		c = (pt->Corona.R_Corona_interp_val / dist_blob_corona) * (pt->Corona.R_Corona_interp_val / dist_blob_corona);
	}

	set_Corona_angles(pt);
	I = integrale_simp_struct(pf, pt, pt->Corona.Corona_mu_1, pt->Corona.Corona_mu_2, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	set_Corona_angles(pt);

	return I * one_by_four_pi * c;
}

double eval_I_nu_Corona_blob_RF(struct blob *pt, double nu_Corona_disk_RF)
{
	double (*pf)(struct blob *, double x);
	double I, R_H_orig, dist_blob_corona, c, R_H_test;

	pt->core.nu_disk_RF = nu_Corona_disk_RF;
	pf = &integrand_I_nu_Corona_blob_RF;

	c = 1.0;
	R_H_orig = pt->core.R_H;
	dist_blob_corona = eval_dist_blob_corona(pt);
	if (dist_blob_corona > pt->Corona.R_Corona_interp_start && dist_blob_corona > 0.0)
	{
		if (pt->core.R_H >= pt->Corona.R_H_Corona){
			R_H_test = pt->Corona.R_H_Corona + pt->Corona.R_Corona_interp_val;
		}
		else{
			R_H_test = pt->Corona.R_H_Corona - pt->Corona.R_Corona_interp_val;
		}
		pt->core.R_H = max(R_H_test, 0.0);
		c = (pt->Corona.R_Corona_interp_val / dist_blob_corona) * (pt->Corona.R_Corona_interp_val / dist_blob_corona);
	}

	set_Corona_angles(pt);
	I = integrale_simp_struct(pf, pt, pt->Corona.Corona_mu_1, pt->Corona.Corona_mu_2, pt->core.theta_n_int);
	pt->core.R_H = R_H_orig;
	set_Corona_angles(pt);

	return I * one_by_four_pi * c;
}

void set_Corona_angles(struct blob *pt)
{
	double mu1, mu2, denom, dist_blob_corona;
	mu1 = 1.0;
	dist_blob_corona = eval_dist_blob_corona(pt);
	denom = sqrt(dist_blob_corona * dist_blob_corona + pt->Corona.R_Corona * pt->Corona.R_Corona);
	if (denom > 0.0){
		mu2 = dist_blob_corona / denom;
	}
	else{
		mu2 = 0.0;
	}
	pt->Corona.Corona_mu_1 = min(mu1, mu2);
	pt->Corona.Corona_mu_2 = max(mu1, mu2);
}

void set_Corona_geometry(struct blob *pt)
{
	pt->Corona.Corona_surface = pi * pt->Corona.R_Corona * pt->Corona.R_Corona;
	if (pt->Corona.Corona_surface > 0.0){
		pt->Corona.Corona_geom_factor = 1.0 / (four_pi * pt->Corona.Corona_surface);
	}
	else{
		pt->Corona.Corona_geom_factor = 0.0;
	}
}

//=========================================================================================








//========================
// Accretion Power Physical Functions
//========================
//double eval_R_Sw(double L_Disk, double accr_eff, double T_disk_max_4){
//	double a;
//	a= pow(0.140836,4) * L_Disk / (pi * sigma_steph_boltz * accr_eff * T_disk_max_4);
//	return pow(a, 0.5);
//}

double eval_R_Sw(double M_BH)
{
	return 2 * M_BH*m_sun * G_cgs / ( vluce_cm*vluce_cm);
}

double eval_M_BH(double R_Sw){
	return R_Sw * vluce_cm * vluce_cm / (2 * G_cgs);
}


double eval_accr_rate(double L_Disk,double accr_eff){
 return L_Disk/(vluce_cm * vluce_cm*accr_eff);
}


double eval_L_Edd(double M_BH){
	return 1.3E38*M_BH;
}

double eval_accr_Edd(double L_Edd, double accr_eff){
	return L_Edd/(vluce_cm*vluce_cm*accr_eff);
}
//=================================================================


//========================
// Planckian Physical Functions
//========================

double eval_nu_peak_planck(double T){
	//Peak of the BB spectrum
	//wien law -> http://en.wikipedia.org/wiki/Wien%27s_displacement_law
	//the 1.39 is to move to get peak fo nu*F(nu)
	return 1.39*5.879e10*T;
}

double eval_T_disk(struct blob *pt, double R)
{
	double  T_disco_r;
	T_disco_r = pt->Disk.Cost_disk_Mulit_BB/(R*R*R) * (1 - pow((pt->Disk.R_inner / R), 0.5));
	T_disco_r = pow(T_disco_r, 0.25);
	//printf("=> T_disco_r %e\n",T_disco_r);
	return T_disco_r;
}

double f_planck_Multi_T(struct blob *pt, double R ,double nu) {
	if (R>pt->Disk.R_ext || R<pt->Disk.R_inner) {
		return 0.0;
	}
	return f_planck(eval_T_disk(pt, R), nu);
}

double f_planck_Multi_T_norm(struct blob *pt, double R, double nu) {
	double T;
	T = eval_T_disk(pt, R);
	return f_planck_norm(T,nu);
}

double integrand_f_planck_Multi_T(struct blob *pt, double R){
	//printf("=> %e %e %e\n", f_planck_Multi_T(pt, R, pt->Disk.nu_disk_Multi_BB), R, pt->Disk.nu_disk_Multi_BB);
	return 2*pi*f_planck_Multi_T(pt,R,pt->Disk.nu_disk_Multi_BB)*R;
}

double f_planck(double T, double nu) {
    double a;
    a = 2 * HPLANCK * pow(nu, 3) / pow(vluce_cm, 2);
    a *= 1.0 / (exp((HPLANCK * nu) / (K_boltz * T)) - 1);
    return a;
}

double f_planck_norm(double T, double nu) {
    return f_planck(T, nu) / ((sigma_steph_boltz / pi) * T * T * T * T);
}



//=================================================================


//========================
//GENERIC SPECTRAL TRANSFORAMTION FUNCTIONS
//========================
double eval_nu_min_blob_RF(struct blob *pt, double mu1, double mu2, double nu_disk_RF ){
	double a,nu_1,nu_2;
	nu_1=nu_disk_RF*pt->core.BulkFactor*(1-pt->core.beta_Gamma*mu1);
	nu_2=nu_disk_RF*pt->core.BulkFactor*(1-pt->core.beta_Gamma*mu2);
	a=  min(nu_1,nu_2);
	return a;
}



double eval_nu_max_blob_RF(struct blob *pt, double mu1, double mu2, double nu_disk_RF ){
	double a,nu_1,nu_2;
	nu_1=nu_disk_RF*pt->core.BulkFactor*(1-pt->core.beta_Gamma*mu1);
	nu_2=nu_disk_RF*pt->core.BulkFactor*(1-pt->core.beta_Gamma*mu2);
	a=  max(nu_1,nu_2);
	return a;
}

double nu_blob_RF_to_nu_disk_RF(double nu_blob_RF, double Gamma, double beta, double mu_disk_RF){
	return nu_blob_RF/(Gamma*(1-beta*mu_disk_RF));

}

double I_nu_disk_RF_to_blob_RF(double I_nu_diks_RF, double nu_disk_RF, double nu_blob_RF, double beta, double Gamma){
	return I_nu_diks_RF/(beta*Gamma)*(nu_blob_RF*nu_blob_RF*nu_blob_RF)/(nu_disk_RF*nu_disk_RF*nu_disk_RF);
}

//========================
//GENERIC GEOMETRICAL TRANSFORAMTION FUNCTIONS
//========================


double eval_circle_secant(double z, double R, double mu)
{
	double x1, x2, y1, y2, m, b, c, a;
	m = tan(acos(mu));
	b = -2 * z;
	c = z * z - R * R;
	a = 1 + m * m;
	if ((b * b - 4 * a * c) > 0)
	{
		x1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
		x2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
		//printf("%e %e\n",x1,x2);
		y1 = m * x1;
		y2 = m * x2;
		return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
	}
	else
		return 0;
}
