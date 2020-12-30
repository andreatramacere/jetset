/***************************************************************************/
/*                   COMMON GRID SPECTRA                                   */
/***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"

/**
 * \file grid_spectra.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief spectra over common grid
 *
 */

void common_grid_spectra(int Num_file, struct blob * pt) {
	double nu_obs, nu_min, nu_max;
	//char somma_obs_log_log[static_file_name_max_legth];
	//char somma_obs[static_file_name_max_legth],somma_obs_src[static_file_name_max_legth];
	double log_nu_start,k;
	unsigned int I_MAX, i;
	//FILE *fp, *fpll, *fpll_src;
	//const char *s;


	//somma_log_log_src_header(fpll_src);
	
	//frequeze osservate a terra
	nu_min = pt->nu_start_grid;
	nu_max = pt->nu_stop_grid;

	I_MAX = pt->nu_grid_size;

	k = (log10(nu_max) - log10(nu_min));
	log_nu_start = log10(nu_min);

	

	for (i = 0; i < I_MAX; i++) {

		nu_obs = pow(10, log_nu_start + k * (double) i / (double) I_MAX);


		//nu_obs=pt->beam_obj*nu/(1+pt->z_cosm);
		//printf("nu=%e nu_obs=%e, i=%d, i_max=%d\n",nu,nu_obs,i,I_MAX);
		interpola_somma(pt, nu_obs,i);
		//pt->nuF_nu_Sum_obs[i]= pt->nuFnu_somma_grid;
		pt->nu_grid[i] = nu_obs;

		//if(nuF_nu_obs>1.e-60){
		//printf("nuF_nu_obs=%e\n********************\n",nuF_nu_obs);
		
		if (pt->nuFnu_sum_grid[i] == 0)
		{
			pt->nuFnu_sum_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_Sync_grid[i] == 0)
		{
			pt->nuFnu_Sync_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_SSC_grid[i] == 0)
		{
			pt->nuFnu_SSC_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_Disk_grid[i] == 0)
		{
			pt->nuFnu_Disk_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_DT_grid[i] == 0)
		{
			pt->nuFnu_DT_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_Star_grid[i] == 0)
		{
			pt->nuFnu_Star_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_Disk_grid[i] == 0)
		{
			pt->nuFnu_EC_Disk_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_BLR_grid[i] == 0)
		{
			pt->nuFnu_EC_BLR_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_DT_grid[i] == 0)
		{
			pt->nuFnu_EC_DT_grid[i] = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_Star_grid[i] == 0)
		{
			pt->nuFnu_EC_Star_grid[i] = pt->emiss_lim;
		}
		if (pt->nuFnu_EC_CMB_grid[i] == 0)
		{
			pt->nuFnu_EC_CMB_grid[i] = pt->emiss_lim;
		}
		if (pt->nuFnu_bremss_ep_grid[i] == 0)
		{
			pt->nuFnu_bremss_ep_grid[i] = pt->emiss_lim;
		}
		if (pt->nuFnu_pp_gamma_grid[i] == 0)
		{
			pt->nuFnu_pp_gamma_grid[i] = pt->emiss_lim;
		}
		if (pt->nuFnu_pp_neutrino_tot_obs[i] == 0)
		{
			pt->nuFnu_pp_neutrino_tot_obs[i] = pt->emiss_lim;
		}
		if (pt->nuFnu_pp_neutrino_mu_obs[i] == 0)
		{
			pt->nuFnu_pp_neutrino_mu_obs[i] = pt->emiss_lim;
		}
		if (pt->nuFnu_pp_neutrino_e_obs[i] == 0)
		{
			pt->nuFnu_pp_neutrino_e_obs[i] = pt->emiss_lim;
		}
		
	}
	
	return;
}

void interpola_somma(struct blob *pt_j, double nu_obs, unsigned int i)
{
	double interp_flux;

	pt_j->nuFnu_sum_grid[i] = 0;

	//Sync
	if (pt_j->do_Sync >= 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_Sync_obs,  pt_j->nu_start_Sync_obs,pt_j->nu_stop_Sync_obs, pt_j->nuF_nu_Sync_obs , pt_j->nu_seed_size, pt_j->emiss_lim);
		//printf("Sync interp_flux=%e\n",interp_flux);
		//printf("Sync interp_flux=%e %e %e %e  %lu \n",interp_flux,nu_obs,  pt_j->nu_start_Sync_obs,pt_j->nu_stop_Sync_obs, pt_j->nu_seed_size);
		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_Sync_grid[i] =  interp_flux;
		}
		else {
			pt_j->nuFnu_Sync_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_Sync_grid[i];
	}

	//SSC
	if (pt_j->do_SSC) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_SSC_obs,  pt_j->nu_start_SSC_obs,pt_j->nu_stop_SSC_obs, pt_j->nuF_nu_SSC_obs , pt_j->nu_IC_size, pt_j->emiss_lim);
		//printf("SSC interp_flux=%e %e %e  %lu \n",interp_flux,  pt_j->nu_start_SSC_obs,pt_j->nu_stop_SSC_obs, pt_j->nu_IC_size);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_SSC_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_SSC_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_SSC_grid[i];
	}

	//EC Disk
	if (pt_j->do_EC_Disk == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_Disk_obs,  pt_j->nu_start_EC_Disk_obs,pt_j->nu_stop_EC_Disk_obs, pt_j->nuF_nu_EC_Disk_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_Disk_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_EC_Disk_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_EC_Disk_grid[i];
	}

	//EC BLR
	if (pt_j->do_EC_BLR == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_BLR_obs,  pt_j->nu_start_EC_BLR_obs,pt_j->nu_stop_EC_BLR_obs, pt_j->nuF_nu_EC_BLR_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_BLR_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_EC_BLR_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_EC_BLR_grid[i];
	}


	//EC DT
	if (pt_j->do_EC_DT == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_DT_obs,  pt_j->nu_start_EC_DT_obs,pt_j->nu_stop_EC_DT_obs, pt_j->nuF_nu_EC_DT_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_DT_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_EC_DT_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_EC_DT_grid[i];
	}


	//EC Star
	if (pt_j->do_EC_Star == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_Star_obs,  pt_j->nu_start_EC_Star_obs,pt_j->nu_stop_EC_Star_obs, pt_j->nuF_nu_EC_Star_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_Star_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_EC_Star_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_EC_Star_grid[i];
	}

	//EC CMB
	if (pt_j->do_EC_CMB == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_CMB_obs,  pt_j->nu_start_EC_CMB_obs,pt_j->nu_stop_EC_CMB_obs, pt_j->nuF_nu_EC_CMB_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_CMB_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_EC_CMB_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_EC_CMB_grid[i];
	}

	//nuFnu_pp_gamma_grid
	if (pt_j->do_pp_gamma == 1)
	{
		interp_flux = log_lin_interp(nu_obs, pt_j->nu_pp_gamma_obs, pt_j->nu_start_pp_gamma_obs, pt_j->nu_stop_pp_gamma_obs, pt_j->nuFnu_pp_gamma_obs, pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_pp_gamma_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_pp_gamma_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_pp_gamma_grid[i];
	} 

	//nuFnu_bress_ep_grid
	if (pt_j->do_bremss_ep == 1)
	{
		interp_flux = log_lin_interp(nu_obs, pt_j->nu_bremss_ep_obs, pt_j->nu_start_bremss_ep_obs, pt_j->nu_stop_bremss_ep_obs, pt_j->nuFnu_bremss_ep_obs, pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_bremss_ep_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_bremss_ep_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_bremss_ep_grid[i];
	} 

	//Disk
	if (pt_j->do_EC_Disk==1 || pt_j->do_EC_BLR==1 || pt_j->do_Disk==1) {

		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_Disk_obs,  pt_j->nu_start_Disk_obs,pt_j->nu_stop_Disk_obs, pt_j->nuF_nu_Disk_obs , pt_j->nu_seed_size, pt_j->emiss_lim);
		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_Disk_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_Disk_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_Disk_grid[i];
	}


	//Dusty Torus
	if (pt_j->do_EC_DT==1 || pt_j->do_DT==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_DT_obs,  pt_j->nu_start_DT_obs,pt_j->nu_stop_DT_obs, pt_j->nuF_nu_DT_obs , pt_j->nu_seed_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_DT_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_DT_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_DT_grid[i];
	}

	//Star
	if (pt_j->do_Star==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_Star_obs,  pt_j->nu_start_Star_obs,pt_j->nu_stop_Star_obs, pt_j->nuF_nu_Star_obs , pt_j->nu_seed_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_Star_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_Star_grid[i] = 0;
		}
		pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_Star_grid[i];
		printf("=> pt->nuFnu_Star_grid %e\n", pt_j->nuFnu_Star_grid[i]);
	}

	//Neutrino
	if (pt_j->do_pp_neutrino==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_pp_neutrino_e_obs,  pt_j->nu_start_pp_neutrino_obs,pt_j->nu_stop_pp_neutrino_obs, pt_j->nuFnu_pp_neutrino_e_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_pp_neutrino_e_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_pp_neutrino_e_grid[i] = 0;
		}
		
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_pp_neutrino_mu_obs,  pt_j->nu_start_pp_neutrino_obs,pt_j->nu_stop_pp_neutrino_obs, pt_j->nuFnu_pp_neutrino_mu_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_pp_neutrino_mu_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_pp_neutrino_mu_grid[i] = 0;
		}

		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_pp_neutrino_tot_obs,  pt_j->nu_start_pp_neutrino_obs,pt_j->nu_stop_pp_neutrino_obs, pt_j->nuFnu_pp_neutrino_tot_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_pp_neutrino_tot_grid[i] = interp_flux;
		}
		else {
			pt_j->nuFnu_pp_neutrino_tot_grid[i] = 0;
		}
		
		//NETURINO NOT SUMMED TO PHOTONS!!
		//pt_j->nuFnu_sum_grid[i] += pt_j->nuFnu_Star_grid[i];
		//printf("=> pt->nuFnu_Star_grid %e\n", pt_j->nuFnu_Star_grid[i]);
	}



	return;
}



