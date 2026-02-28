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
	nu_min = pt->core.nu_start_grid;
	nu_max = pt->core.nu_stop_grid;

	I_MAX = pt->core.nu_grid_size;

	k = (log10(nu_max) - log10(nu_min));
	log_nu_start = log10(nu_min);

	

	for (i = 0; i < I_MAX; i++) {

		nu_obs = pow(10, log_nu_start + k * (double) i / (double) I_MAX);


		//nu_obs=pt->core.beam_obj*nu/(1+pt->core.z_cosm);
		//printf("nu=%e nu_obs=%e, i=%d, i_max=%d\n",nu,nu_obs,i,I_MAX);
		interpola_somma(pt, nu_obs,i);
		//pt->nuF_nu_Sum_obs[i]= pt->nuFnu_somma_grid;
		pt->core.nu_grid[i] = nu_obs;

		//if(nuF_nu_obs>1.e-60){
		//printf("nuF_nu_obs=%e\n********************\n",nuF_nu_obs);
		
		if (pt->core.nuFnu_sum_grid[i] == 0)
		{
			pt->core.nuFnu_sum_grid[i] = pt->core.emiss_lim;
		}

		if (pt->Sync.spec.nuFnu_grid[i] == 0)
		{
			pt->Sync.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->SSC.spec.nuFnu_grid[i] == 0)
		{
			pt->SSC.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->Disk.spec.nuFnu_grid[i] == 0)
		{
			pt->Disk.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->DT.spec.nuFnu_grid[i] == 0)
		{
			pt->DT.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->Star.spec.nuFnu_grid[i] == 0)
		{
			pt->Star.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->Disk.ec.spec.nuFnu_grid[i] == 0)
		{
			pt->Disk.ec.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->BLR.ec.spec.nuFnu_grid[i] == 0)
		{
			pt->BLR.ec.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->DT.ec.spec.nuFnu_grid[i] == 0)
		{
			pt->DT.ec.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}

		if (pt->Star.ec.spec.nuFnu_grid[i] == 0)
		{
			pt->Star.ec.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}
		if (pt->CMB.ec.spec.nuFnu_grid[i] == 0)
		{
			pt->CMB.ec.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}
		if (pt->Bremss_ep.spec.nuFnu_grid[i] == 0)
		{
			pt->Bremss_ep.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}
		if (pt->PP_gamma.spec.nuFnu_grid[i] == 0)
		{
			pt->PP_gamma.spec.nuFnu_grid[i] = pt->core.emiss_lim;
		}
		if (pt->PP_neutrino.spec_tot.nuFnu_obs[i] == 0)
		{
			pt->PP_neutrino.spec_tot.nuFnu_obs[i] = pt->core.emiss_lim;
		}
		if (pt->PP_neutrino.spec_mu.nuFnu_obs[i] == 0)
		{
			pt->PP_neutrino.spec_mu.nuFnu_obs[i] = pt->core.emiss_lim;
		}
		if (pt->PP_neutrino.spec_e.nuFnu_obs[i] == 0)
		{
			pt->PP_neutrino.spec_e.nuFnu_obs[i] = pt->core.emiss_lim;
		}
		
	}
	
	return;
}

void interpola_somma(struct blob *pt_j, double nu_obs, unsigned int i)
{
	double interp_flux;

	pt_j->core.nuFnu_sum_grid[i] = 0;

	//Sync
	if (pt_j->core.do_Sync >= 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->Sync.spec.nu_obs,  pt_j->Sync.spec.nu_min_obs,pt_j->Sync.spec.nu_max_obs, pt_j->Sync.spec.nuFnu_obs , pt_j->core.nu_seed_size, pt_j->core.emiss_lim);
		//printf("Sync interp_flux=%e\n",interp_flux);
		//printf("Sync interp_flux=%e %e %e %e  %lu \n",interp_flux,nu_obs,  pt_j->Sync.spec.nu_min_obs,pt_j->Sync.spec.nu_max_obs, pt_j->core.nu_seed_size);
		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->Sync.spec.nuFnu_grid[i] =  interp_flux;
		}
		else {
			pt_j->Sync.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->Sync.spec.nuFnu_grid[i];
	}

	//SSC
	if (pt_j->core.do_SSC) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->SSC.spec.nu_obs,  pt_j->SSC.spec.nu_min_obs,pt_j->SSC.spec.nu_max_obs, pt_j->SSC.spec.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);
		//printf("SSC interp_flux=%e %e %e  %lu \n",interp_flux,  pt_j->SSC.spec.nu_min_obs,pt_j->SSC.spec.nu_max_obs, pt_j->core.nu_IC_size);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->SSC.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->SSC.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->SSC.spec.nuFnu_grid[i];
	}

	//EC Disk
	if (pt_j->core.do_EC_Disk == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->Disk.ec.spec.nu_obs,  pt_j->Disk.ec.spec.nu_min_obs,pt_j->Disk.ec.spec.nu_max_obs, pt_j->Disk.ec.spec.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->Disk.ec.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->Disk.ec.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->Disk.ec.spec.nuFnu_grid[i];
	}

	//EC BLR
	if (pt_j->core.do_EC_BLR == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->BLR.ec.spec.nu_obs,  pt_j->BLR.ec.spec.nu_min_obs,pt_j->BLR.ec.spec.nu_max_obs, pt_j->BLR.ec.spec.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->BLR.ec.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->BLR.ec.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->BLR.ec.spec.nuFnu_grid[i];
	}


	//EC DT
	if (pt_j->core.do_EC_DT == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->DT.ec.spec.nu_obs,  pt_j->DT.ec.spec.nu_min_obs,pt_j->DT.ec.spec.nu_max_obs, pt_j->DT.ec.spec.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->DT.ec.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->DT.ec.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->DT.ec.spec.nuFnu_grid[i];
	}


	//EC Star
	if (pt_j->core.do_EC_Star == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->Star.ec.spec.nu_obs,  pt_j->Star.ec.spec.nu_min_obs,pt_j->Star.ec.spec.nu_max_obs, pt_j->Star.ec.spec.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->Star.ec.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->Star.ec.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->Star.ec.spec.nuFnu_grid[i];
	}

	//EC CMB
	if (pt_j->core.do_EC_CMB == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->CMB.ec.spec.nu_obs,  pt_j->CMB.ec.spec.nu_min_obs,pt_j->CMB.ec.spec.nu_max_obs, pt_j->CMB.ec.spec.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->CMB.ec.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->CMB.ec.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->CMB.ec.spec.nuFnu_grid[i];
	}

	//nuFnu_pp_gamma_grid
	if (pt_j->PP_gamma.do_pp_gamma == 1)
	{
		interp_flux = log_lin_interp(nu_obs, pt_j->PP_gamma.spec.nu_obs, pt_j->PP_gamma.spec.nu_min_obs, pt_j->PP_gamma.spec.nu_max_obs, pt_j->PP_gamma.spec.nuFnu_obs, pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->PP_gamma.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->PP_gamma.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->PP_gamma.spec.nuFnu_grid[i];
	} 

	//nuFnu_bress_ep_grid
	if (pt_j->Bremss_ep.do_bremss_ep == 1)
	{
		interp_flux = log_lin_interp(nu_obs, pt_j->Bremss_ep.spec.nu_obs, pt_j->Bremss_ep.spec.nu_min_obs, pt_j->Bremss_ep.spec.nu_max_obs, pt_j->Bremss_ep.spec.nuFnu_obs, pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->Bremss_ep.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->Bremss_ep.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->Bremss_ep.spec.nuFnu_grid[i];
	} 

	//Disk
	if (pt_j->core.do_EC_Disk==1 || pt_j->core.do_EC_BLR==1 || pt_j->core.do_Disk==1) {

		interp_flux=log_lin_interp( nu_obs,  pt_j->Disk.spec.nu_obs,  pt_j->Disk.spec.nu_min_obs,pt_j->Disk.spec.nu_max_obs, pt_j->Disk.spec.nuFnu_obs , pt_j->core.nu_seed_size, pt_j->core.emiss_lim);
		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->Disk.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->Disk.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->Disk.spec.nuFnu_grid[i];
	}


	//Dusty Torus
	if (pt_j->core.do_EC_DT==1 || pt_j->core.do_DT==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->DT.spec.nu_obs,  pt_j->DT.spec.nu_min_obs,pt_j->DT.spec.nu_max_obs, pt_j->DT.spec.nuFnu_obs , pt_j->core.nu_seed_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->DT.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->DT.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->DT.spec.nuFnu_grid[i];
	}

	//Star
	if (pt_j->core.do_Star==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->Star.spec.nu_obs,  pt_j->Star.spec.nu_min_obs,pt_j->Star.spec.nu_max_obs, pt_j->Star.spec.nuFnu_obs , pt_j->core.nu_seed_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->Star.spec.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->Star.spec.nuFnu_grid[i] = 0;
		}
		pt_j->core.nuFnu_sum_grid[i] += pt_j->Star.spec.nuFnu_grid[i];
		//printf("=> pt->Star.spec.nuFnu_grid %e\n", pt_j->Star.spec.nuFnu_grid[i]);
	}

	//Neutrino
	if (pt_j->PP_neutrino.do_pp_neutrino==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->PP_neutrino.spec_e.nu_obs,  pt_j->PP_neutrino.spec_tot.nu_min_obs,pt_j->PP_neutrino.spec_tot.nu_max_obs, pt_j->PP_neutrino.spec_e.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->PP_neutrino.spec_e.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->PP_neutrino.spec_e.nuFnu_grid[i] = 0;
		}
		
		interp_flux=log_lin_interp( nu_obs,  pt_j->PP_neutrino.spec_mu.nu_obs,  pt_j->PP_neutrino.spec_tot.nu_min_obs,pt_j->PP_neutrino.spec_tot.nu_max_obs, pt_j->PP_neutrino.spec_mu.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->PP_neutrino.spec_mu.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->PP_neutrino.spec_mu.nuFnu_grid[i] = 0;
		}

		interp_flux=log_lin_interp( nu_obs,  pt_j->PP_neutrino.spec_tot.nu_obs,  pt_j->PP_neutrino.spec_tot.nu_min_obs,pt_j->PP_neutrino.spec_tot.nu_max_obs, pt_j->PP_neutrino.spec_tot.nuFnu_obs , pt_j->core.nu_IC_size, pt_j->core.emiss_lim);

		if (interp_flux > pt_j->core.emiss_lim) {
			pt_j->PP_neutrino.spec_tot.nuFnu_grid[i] = interp_flux;
		}
		else {
			pt_j->PP_neutrino.spec_tot.nuFnu_grid[i] = 0;
		}
		
		//NETURINO NOT SUMMED TO PHOTONS!!
		//pt_j->core.nuFnu_sum_grid[i] += pt_j->Star.spec.nuFnu_grid[i];
		//printf("=> pt->Star.spec.nuFnu_grid %e\n", pt_j->Star.spec.nuFnu_grid[i]);
	}



	return;
}



