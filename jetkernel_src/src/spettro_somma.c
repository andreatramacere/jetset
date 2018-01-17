/***************************************************************************/
/*                   SOMMA SPETTRO IC*Sync                                 */
/***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"

/**
 * \file spettro_Compton.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief Somma IC+Sync
 *
 */

void spettro_somma_Sync_ic(int Num_file, struct spettro * pt) {
	double nu_obs,nu_src, nu, nu_min, nu_max;
	char somma_obs_log_log[static_file_name_max_legth];
	char somma_obs[static_file_name_max_legth],somma_obs_src[static_file_name_max_legth];
	double log_nu, log_nu_start, nuF_nu_obs, k;
	unsigned long I_MAX, i;
	FILE *fp, *fpll, *fpll_src;
	const char *s;
	//printf("**********************  CALCOLO DELLO SPETTRO SOMMA  *******************************\n");

	sprintf(somma_obs, "%s%s-somma.dat",
			pt->path, pt->STEM);

	fp = fopen(somma_obs, "w");
	if (fp == NULL) {
		printf("non posso aprire %s\n ", somma_obs);
		perror(s);
		exit(1);
		;
	}


	sprintf(somma_obs_log_log, "%s%s-somma-log-log.dat",
			pt->path, pt->STEM);

	fpll = fopen(somma_obs_log_log, "w");
	if (fpll == NULL) {
		printf("non posso aprire %s\n ", somma_obs_log_log);
		perror(s);
		exit(1);
	}



	sprintf(somma_obs_src, "%s%s-somma-log-log-src.dat",
			pt->path, pt->STEM);
	fpll_src = fopen(somma_obs_src, "w");
	if (fpll_src == NULL) {
		printf("non posso aprire %s\n ", somma_obs_src);
		perror(s);
		exit(1);
	}

	//frequeze osservate a terra
	nu_min = pt->nu_start_grid;
	nu_max = pt->nu_stop_grid;

	I_MAX = pt->nu_sum_size;

	k = (log10(nu_max) - log10(nu_min));
	log_nu_start = log10(nu_min);

	somma_header(fp);
	somma_log_log_header(fpll);
	somma_log_log_src_header(fpll_src);

	for (i = 0; i < I_MAX; i++) {

		nu_obs = pow(10, log_nu_start + k * (double) i / (double) I_MAX);


		//nu_obs=pt->beam_obj*nu/(1+pt->z_cosm);
		//printf("nu=%e nu_obs=%e now call\n",nu,nu_obs);
		interpola_somma(pt, nu_obs);
		pt->nuF_nu_Sum_obs[i]= pt->nuFnu_somma_grid;
		pt->nu_Sum_obs[i]= nu_obs;

		//if(nuF_nu_obs>1.e-60){
		//printf("nuF_nu_obs=%e\n********************\n",nuF_nu_obs);
		fprintf(fp, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\n",
				nu_obs,
				pt->nuFnu_somma_grid,
				pt->nuFnu_Sync_grid,
				pt->nuFnu_comp_grid,
				pt->nuFnu_Disk_grid,
				pt->nuFnu_DT_grid,
				pt->nuFnu_Star_grid,
				pt->nuFnu_EC_Disk_grid,
				pt->nuFnu_EC_BLR_grid,
				pt->nuFnu_EC_DT_grid,
				pt->nuFnu_EC_Star_grid,
				pt->nuFnu_EC_CMB_grid,
				pt->nuFnu_EC_CMB_stat_grid);



		if (pt->nuFnu_somma_grid == 0) {
			pt->nuFnu_somma_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_Sync_grid == 0) {
			pt->nuFnu_Sync_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_comp_grid == 0) {
			pt->nuFnu_comp_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_Disk_grid == 0) {
			pt->nuFnu_Disk_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_DT_grid == 0) {
			pt->nuFnu_DT_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_Star_grid == 0) {
			pt->nuFnu_Star_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_Disk_grid == 0) {
			pt->nuFnu_EC_Disk_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_BLR_grid == 0) {
			pt->nuFnu_EC_BLR_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_DT_grid == 0) {
			pt->nuFnu_EC_DT_grid = pt->emiss_lim;
		}

		if (pt->nuFnu_EC_Star_grid == 0) {
			pt->nuFnu_EC_Star_grid = pt->emiss_lim;
		}
		if (pt->nuFnu_EC_CMB_grid == 0) {
			pt->nuFnu_EC_CMB_grid = pt->emiss_lim;
		}
		if (pt->nuFnu_EC_CMB_stat_grid == 0) {
			pt->nuFnu_EC_CMB_stat_grid = pt->emiss_lim;
		}

		fprintf(fpll, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\n",
				log10(nu_obs),
				log10(pt->nuFnu_somma_grid),
				log10(pt->nuFnu_Sync_grid),
				log10(pt->nuFnu_comp_grid),
				log10(pt->nuFnu_Disk_grid),
				log10(pt->nuFnu_DT_grid),
				log10(pt->nuFnu_Star_grid),
				log10(pt->nuFnu_EC_Disk_grid),
				log10(pt->nuFnu_EC_BLR_grid),
				log10(pt->nuFnu_EC_DT_grid),
				log10(pt->nuFnu_EC_Star_grid),
				log10(pt->nuFnu_EC_CMB_grid),
				log10(pt->nuFnu_EC_CMB_stat_grid)
		);
		//fprintf(fp1, "%e %e \n",
		//        log10(nu_obs),
		//        log10(nuF_nu_obs));

		//}

		fprintf(fpll_src, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\t %4.4e\n",
				log10(nu_obs_to_nu_src(nu_obs, pt->z_cosm)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_somma_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_Sync_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_comp_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_Disk_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_DT_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_Star_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_EC_Disk_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_EC_BLR_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_EC_DT_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_EC_Star_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_EC_CMB_grid,pt->beam_obj, pt->z_cosm, pt->dist)),
				log10(nuFnu_obs_to_nuLnu_src(pt->nuFnu_EC_CMB_stat_grid,pt->beam_obj, pt->z_cosm, pt->dist))
		);

	}
	fclose(fp);
	fclose(fpll_src);
	fclose(fpll);
	return;
}

void interpola_somma(struct spettro *pt_j, double nu_obs) {
	double interp_flux;



	pt_j->nuFnu_somma_grid = 0;


	//Sync
	if (pt_j->do_Sync >= 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_Sync_obs,  pt_j->nu_start_Sync_obs,pt_j->nu_stop_Sync_obs, pt_j->nuF_nu_Sync_obs , pt_j->nu_seed_size, pt_j->emiss_lim);
		//printf("Sync interp_flux=%e\n",interp_flux);
		//printf("Sync interp_flux=%e %e %e %e  %lu \n",interp_flux,nu_obs,  pt_j->nu_start_Sync_obs,pt_j->nu_stop_Sync_obs, pt_j->nu_seed_size);
		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_Sync_grid =  interp_flux;
		}
		else {
			pt_j->nuFnu_Sync_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_Sync_grid;
	}

	//SSC
	if (pt_j->do_SSC) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_SSC_obs,  pt_j->nu_start_SSC_obs,pt_j->nu_stop_SSC_obs, pt_j->nuF_nu_SSC_obs , pt_j->nu_IC_size, pt_j->emiss_lim);
		//printf("SSC interp_flux=%e %e %e  %lu \n",interp_flux,  pt_j->nu_start_SSC_obs,pt_j->nu_stop_SSC_obs, pt_j->nu_IC_size);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_comp_grid =   interp_flux;
		}
		else {
			pt_j->nuFnu_comp_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_comp_grid;
	}

	//EC Disk
	if (pt_j->do_EC_Disk == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_Disk_obs,  pt_j->nu_start_EC_Disk_obs,pt_j->nu_stop_EC_Disk_obs, pt_j->nuF_nu_EC_Disk_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_Disk_grid =  interp_flux;
		}
		else {
			pt_j->nuFnu_EC_Disk_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_EC_Disk_grid;
	}

	//EC BLR
	if (pt_j->do_EC_BLR == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_BLR_obs,  pt_j->nu_start_EC_BLR_obs,pt_j->nu_stop_EC_BLR_obs, pt_j->nuF_nu_EC_BLR_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_BLR_grid =  interp_flux;
		}
		else {
			pt_j->nuFnu_EC_BLR_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_EC_BLR_grid;
	}


	//EC DT
	if (pt_j->do_EC_DT == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_DT_obs,  pt_j->nu_start_EC_DT_obs,pt_j->nu_stop_EC_DT_obs, pt_j->nuF_nu_EC_DT_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_DT_grid =   interp_flux;
		}
		else {
			pt_j->nuFnu_EC_DT_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_EC_DT_grid;
	}


	//EC Star
	if (pt_j->do_EC_Star == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_Star_obs,  pt_j->nu_start_EC_Star_obs,pt_j->nu_stop_EC_Star_obs, pt_j->nuF_nu_EC_Star_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_Star_grid =   interp_flux;
		}
		else {
			pt_j->nuFnu_EC_Star_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_EC_Star_grid;
	}

	//EC CMB
	if (pt_j->do_EC_CMB == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_CMB_obs,  pt_j->nu_start_EC_CMB_obs,pt_j->nu_stop_EC_CMB_obs, pt_j->nuF_nu_EC_CMB_obs , pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_CMB_grid =   interp_flux;
		}
		else {
			pt_j->nuFnu_EC_CMB_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_EC_CMB_grid;
	}

	//EC CMB 
	if (pt_j->do_EC_CMB_stat == 1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_EC_CMB_stat_obs,  pt_j->nu_start_EC_CMB_stat_obs,pt_j->nu_stop_EC_CMB_stat_obs, pt_j->nuF_nu_EC_CMB_stat_obs, pt_j->nu_IC_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_EC_CMB_stat_grid =   interp_flux;
		}
		else {
			pt_j->nuFnu_EC_CMB_stat_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_EC_CMB_stat_grid;
	}

	//Disk
	if (pt_j->do_EC_Disk==1 || pt_j->do_EC_BLR==1) {

		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_Disk_obs,  pt_j->nu_start_Disk_obs,pt_j->nu_stop_Disk_obs, pt_j->nuF_nu_Disk_obs , pt_j->nu_seed_size, pt_j->emiss_lim);
		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_Disk_grid =   interp_flux;
		}
		else {
			pt_j->nuFnu_Disk_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_Disk_grid;
	}


	//Dusty Torus
	if (pt_j->do_EC_DT==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_DT_obs,  pt_j->nu_start_DT_obs,pt_j->nu_stop_DT_obs, pt_j->nuF_nu_DT_obs , pt_j->nu_seed_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_DT_grid =  interp_flux;
		}
		else {
			pt_j->nuFnu_DT_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_DT_grid;
	}

	//Star
	if (pt_j->do_EC_Star==1) {
		interp_flux=log_lin_interp( nu_obs,  pt_j->nu_Star_obs,  pt_j->nu_start_Star_obs,pt_j->nu_stop_Star_obs, pt_j->nuF_nu_Star_obs , pt_j->nu_seed_size, pt_j->emiss_lim);

		if (interp_flux > pt_j->emiss_lim) {
			pt_j->nuFnu_Star_grid =  interp_flux;
		}
		else {
			pt_j->nuFnu_Star_grid = 0;

		}
		pt_j->nuFnu_somma_grid += pt_j->nuFnu_Star_grid;
	}

	return;
}



