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

//===============================================================
// Evaluation of external radiative fields
//===============================================================

void spectra_External_Fields(int Num_file, struct blob *pt) {

    //==================================================================
	//if (pt->verbose){
	if (pt->verbose > 0)
	{
		printf("**********************   Eval. seed photon fields for  EC       *******************************\n");
	}


    // ====================================================
    // approx gamma con beaming factor, impilcit teta circa= 1/gamma
    // Set nu start EC
	// not used in photon field computation
	// used only in analytic approx on screen
    //=====================================================
	pt->beaming_EC = pt->BulkFactor;


	if (pt->do_EC_Star==1 || pt->do_Star==1){
		set_EC_stat_pre(pt, pt->R_Star);
    	Build_I_nu_Star(pt);
    }
	if (pt->do_EC_Disk == 1 || pt->do_EC_BLR == 1 || pt->do_Disk == 1 || pt->do_EC_DT == 1 || pt->do_DT ==1)
	{
		set_EC_stat_pre(pt, pt->R_ext);
		Build_I_nu_Disk(pt);
		set_EC_stat_post(pt);
    }
    if (pt->do_EC_BLR==1){
		set_EC_stat_pre(pt, pt->R_BLR_out);
		Build_I_nu_BLR(pt);
		set_EC_stat_post(pt);
	}
    if (pt->do_EC_DT==1 || pt->do_DT==1){
		//printf("EC_stat=%d, R_H=%e\n",pt->EC_stat,pt->R_H);
		set_EC_stat_pre(pt, pt->R_DT);
    	Build_I_nu_DT(pt);
		//printf("EC_stat=%d, R_H=%e\n",pt->EC_stat,pt->R_H);
		set_EC_stat_post(pt);
    }
    if (pt->do_EC_CMB==1){
		set_EC_stat_pre(pt, -1);
    	Build_I_nu_CMB(pt);
		set_EC_stat_post(pt);
    }

    //if (pt->do_EC_CMB_stat==1){
    //	Build_I_nu_CMB_stat(pt);
    //}
	if (pt->verbose > 1)
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
	double nu_start_disk_RF;
	double nu_stop_disk_RF;
	double nuL_nu_disk,F_nu_disk_obs;

	/*
	sprintf(f_SED_star, "%s%s-SED-star.dat",pt->path, pt->STEM);

	if (pt->WRITE_TO_FILE==1){
		fp_SED_star = fopen(f_SED_star, "w");
		if (fp_SED_star == NULL) {
			printf("unable to open %s\n ", fp_SED_star);
			exit(1);
		}
		flux_DISK_header(fp_SED_star);
	}
	*/

	//set_Star_geometry(pt);

	pt->Star_mu_1=0;
	pt->Star_mu_2=1;
	pt->Star_surface = 4 * pi * pt->R_Star * pt->R_Star;
   
	nu_peak_BB=eval_nu_peak_Disk(pt->T_Star_max);

	nu_start_disk_RF = nu_peak_BB*pt->nu_planck_min_factor;
	nu_stop_disk_RF  = nu_peak_BB*pt->nu_planck_max_factor;

	pt->nu_start_Star = eval_nu_min_blob_RF(pt,pt->Star_mu_1, pt->Star_mu_2, nu_start_disk_RF);
	pt->nu_stop_Star  = eval_nu_max_blob_RF(pt,pt->Star_mu_1, pt->Star_mu_2, nu_stop_disk_RF);

	pt->nu_start_Star_DRF = nu_start_disk_RF;
	pt->nu_stop_Star_DRF = nu_stop_disk_RF;


	NU_INT_MAX=pt->nu_seed_size-1;
	pt->NU_INT_MAX_Star = NU_INT_MAX;


	pt->nu_start_Star_obs=nu_disk_to_nu_obs_disk(nu_start_disk_RF , pt->z_cosm);
	pt->nu_stop_Star_obs=nu_disk_to_nu_obs_disk(nu_stop_disk_RF, pt->z_cosm);

	if (pt->verbose)
	{
		printf("-----------  Building I_nu Star     ----------- \n");

		printf("nu_start_Star=%e  nu_stop_Star=%e \n",
			   pt->nu_start_Star,
			   pt->nu_stop_Star);

		printf("nu_start_Star_disk_RF=%e  nu_stop_Star_disk_RF=%e \n",
			   nu_start_disk_RF,
			   nu_stop_disk_RF);

		printf("nu_start_Star_obs=%e  nu_stop_Star_obs=%e \n",
			   pt->nu_start_Star_obs,
			   pt->nu_stop_Star_obs);
	}

	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->nu_seed_size, pt->nu_Star_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->I_nu_Star_disk_RF[NU_INT]=eval_I_nu_Star_disk_RF(pt, pt->nu_Star_disk_RF[NU_INT]);
		//pt->J_nu_Star_disk_RF[NU_INT]=eval_J_nu_Star_disk_RF(pt, pt->I_nu_Star_disk_RF[NU_INT]);

	}



	build_log_grid( pt->nu_start_Star,  pt->nu_stop_Star, pt->nu_seed_size, pt->nu_Star);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		nu_obs = nu_disk_to_nu_obs_disk(pt->nu_Star_disk_RF[NU_INT],pt->z_cosm);
		pt->nu_Star_obs[NU_INT]=nu_obs;

		pt->I_nu_Star[NU_INT]=eval_I_nu_Star_blob_RF(pt,pt->nu_Star[NU_INT]);
		pt->n_Star[NU_INT] =I_nu_to_n(pt->I_nu_Star[NU_INT], pt->nu_Star[NU_INT]);
		//EC with n(gamma) transf
		pt->n_Star_DRF[NU_INT] = I_nu_to_n(pt->J_nu_Star_disk_RF[NU_INT], pt->nu_Star_disk_RF[NU_INT]);
		

		if (pt->I_nu_Star[NU_INT]>pt->emiss_lim){
			pt->nu_stop_Star = pt->nu_Star[NU_INT];
			pt->NU_INT_MAX_Star = NU_INT;
		}
		else{
			pt->I_nu_Star[NU_INT]=pt->emiss_lim;
			pt->n_Star[NU_INT] =I_nu_to_n(pt->I_nu_Star[NU_INT], pt->nu_Star[NU_INT]);

		}

		nuL_nu_disk = eval_Star_L_nu(pt,pt->nu_Star_disk_RF[NU_INT]) * pt->nu_Star_disk_RF[NU_INT];
		F_nu_disk_obs= L_nu_Disk_to_F_nu(nuL_nu_disk / pt->nu_Star_disk_RF[NU_INT], pt-> z_cosm, pt-> dist);
		pt->nuF_nu_Star_obs[NU_INT] = F_nu_disk_obs*nu_obs;
		if (pt->verbose > 1)
		{
			printf(" nu_Star_disk_RF=%e, nuF_nu_Star_obs=%e, nu_Star=%e, , I_nu_Star=%e,  nuL_nu_disk=%e, Star surface=%e nu_Star_obs=%e\n",
				   pt->nu_Star_disk_RF[NU_INT],
				   pt->nuF_nu_Star_obs[NU_INT],
				   pt->nu_Star[NU_INT],
				   pt->I_nu_Star_disk_RF[NU_INT],
				   nuL_nu_disk,
				   pt->Star_surface,
				   pt->nu_Star_obs[NU_INT]);
		}

		/*
		if (pt->WRITE_TO_FILE==1){
			fprintf(fp_SED_star, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
				log10(nu_obs),
				log10(nu_obs * F_nu_disk_obs),
				nu_obs,
				nu_obs*F_nu_disk_obs,
				pt->nu_Star_disk_RF[NU_INT],
				nuL_nu_disk);
		}
		*/

	}
	/*
	if (pt->WRITE_TO_FILE == 1)
	{
		fclose(fp_SED_star);
	}
	*/
}


//========================
// Star Spectral Functions
//========================

double eval_I_nu_Star_disk_RF(struct blob *pt,double nu_Star_disk_RF){
	return f_planck(pt->T_Star_max, nu_Star_disk_RF);
	
}

//TODO!! THIS MUST BE IMPROVED FOR PROPER ANGULAR INTEGRATION
//double eval_J_nu_Star_disk_RF(struct spettro *pt, double I_nu_Star_disk_RF){
//	return I_nu_Star_disk_RF*(2.0*pi/(4*pi))*(pt->Star_mu_2- pt->Star_mu_1);
//}


double integrand_I_nu_Star_blob_RF(struct blob *pt, double mu){
	int i;
	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->nu_blob_RF,pt->BulkFactor,pt->beta_Gamma,mu);

	i=x_to_grid_index( pt->nu_Star_disk_RF,nu_disk_RF,pt->nu_seed_size);
	if (i>0){
		return pt->J_nu_Star_disk_RF[i]*pt->BulkFactor*(1-pt->beta_Gamma*mu);
	}
	else{
		return 0;
	}
}



double eval_I_nu_Star_blob_RF(struct blob *pt, double nu_blob_RF){
	pt->nu_blob_RF=nu_blob_RF;
	double (*pf) (struct blob *, double x);
	pf = &integrand_I_nu_Star_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return  0.5*integrale_simp_struct(pf, pt, pt->Star_mu_1, pt->Star_mu_2,pt->theta_n_int);
}

double eval_Star_L_nu(struct blob *pt, double nu_Star_disk_RF){
	return  pt->Star_surface *eval_I_nu_Star_disk_RF(pt, nu_Star_disk_RF)*pi;
}


//========================
// Star Geometrical Functions
//========================

void set_Star_geometry(struct blob *pt){
	double mu1,mu2,mu,mu_star,c,b,a,psi_star;

	mu=1.0;
	b=pt->R_H-pt->R_Star*cos(pt->Star_psi_1* pi/ 180.0);
	a=pt->R_Star*sin(pt->Star_psi_1* pi/ 180.0);
	c=sqrt(a*a+b*b);
	printf("%e %e %e\n",a,b,c);
	mu_star=b/c;

	//max cos=> min angle (0-pi)
	mu1=min(mu,mu_star);

	if (pt->verbose){
		printf("mu1=%e mu_star=%e  mu=%e \n",mu1,mu_star,mu);
	}


	b=sqrt(pt->R_H*pt->R_H - pt->R_Star*pt->R_Star  );
	mu_star=b/pt->R_H;
	psi_star = asin(mu_star) * 180.0 / pi;


	mu2=min(pt->Star_psi_2,psi_star);
	mu2=cos(mu2* pi/ 180.0);


	if (psi_star<=pt->Star_psi_2){
		c=pt->R_H;
		a=pt->R_Star;
		b=sqrt(c*c - a*a);
		mu2=b/c;
	}else{
		b=pt->R_H - pt->R_Star*cos(pt->Star_psi_2* pi/ 180.0);
		a=pt->R_Star*sin(pt->Star_psi_2* pi/ 180.0);
		c=sqrt(b*b + a*a);
		mu2=b/c;
	}
	//printf("mu2=%e mu_star=%e  psi_star=%e Star_psi_2=%e \n",mu2,mu_star,psi_star,pt->Star_psi_2);

	pt->Star_mu_1=min(mu1,mu2);
	pt->Star_mu_2=max(mu1,mu2);


	if (pt->verbose){
		printf("mu1=%e mu2=%e psi_star=%e ,mu=%e \n",pt->Star_mu_1,pt->Star_mu_2,psi_star,mu);
	}

	pt->Star_surface=4*pi*pt->R_Star*pt->R_Star;

}
//=========================================================================================


//=========================================================================================
void Build_I_nu_CMB(struct blob *pt){
	double T_CMB_z;
	double nu_peak_CMB_z;
	unsigned int NU_INT,NU_INT_MAX;
	double nu_start_disk_RF;
	double nu_stop_disk_RF;

	pt->CMB_mu_1=-1.0;
	pt->CMB_mu_2=1.0;

	T_CMB_z=eval_T_CMB_z(pt->z_cosm,pt->T_CMB_0);
	//T_CMB_0=pt->T_CMB_0;

	nu_peak_CMB_z=eval_nu_peak_planck(T_CMB_z);
	//nu_peak_CMB_0=eval_nu_peak_planck(T_CMB_0);

	nu_start_disk_RF = nu_peak_CMB_z*pt->nu_planck_min_factor;
	nu_stop_disk_RF  = nu_peak_CMB_z*pt->nu_planck_max_factor;

	pt->nu_start_CMB = eval_nu_min_blob_RF(pt,-1, 1, nu_start_disk_RF);
	pt->nu_stop_CMB  = eval_nu_max_blob_RF(pt,-1, 1, nu_stop_disk_RF);

	pt->nu_start_CMB_DRF = nu_start_disk_RF;
	pt->nu_stop_CMB_DRF = nu_stop_disk_RF;
	//pt->nu_start_CMB_obs=nu_peak_CMB_0*pt->nu_planck_min_factor;
	//pt->nu_stop_CMB_obs=nu_peak_CMB_0*pt->nu_planck_max_factor;

	NU_INT_MAX=pt->nu_seed_size-1;
	pt->NU_INT_MAX_CMB = NU_INT_MAX;

	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->nu_seed_size, pt->nu_CMB_disk_RF);
	
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
			pt->I_nu_CMB_disk_RF[NU_INT]=eval_I_nu_CMB_disk_RF(T_CMB_z, pt->nu_CMB_disk_RF[NU_INT]);
	}
	build_log_grid( pt->nu_start_CMB,  pt->nu_stop_CMB, pt->nu_seed_size, pt->nu_CMB);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->I_nu_CMB[NU_INT]=eval_I_nu_CMB_blob_RF(pt,pt->nu_CMB[NU_INT]);
		pt->n_CMB[NU_INT] =I_nu_to_n(pt->I_nu_CMB[NU_INT], pt->nu_CMB[NU_INT]);
		//EC with n(gamma) transf
		pt->n_CMB_DRF[NU_INT] = I_nu_to_n(pt->I_nu_CMB_disk_RF[NU_INT], pt->nu_CMB_disk_RF[NU_INT]);
	}
	
}



double eval_T_CMB_z(double z, double T_CMB_0){
		return T_CMB_0*(1+z);
}

double eval_I_nu_CMB_disk_RF(double T_CMB,double nu_CMB_disk_RF){
	return f_planck(T_CMB, nu_CMB_disk_RF);
}


double eval_I_nu_CMB_blob_RF(struct blob *pt, double nu_blob_RF){


	pt->nu_blob_RF=nu_blob_RF;
	double (*pf) (struct blob *, double x);
	pf = &integrand_I_nu_CMB_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return 0.5 * integrale_simp_struct(pf, pt, pt->CMB_mu_1, pt->CMB_mu_2, pt->theta_n_int);
}

double integrand_I_nu_CMB_blob_RF(struct blob *pt, double mu){
	int i=0;
 	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->nu_blob_RF,pt->BulkFactor,pt->beta_Gamma,mu);
	i=x_to_grid_index( pt->nu_CMB_disk_RF,nu_disk_RF,pt->nu_seed_size);
	if (i>0){
		return pt->I_nu_CMB_disk_RF[i]*pt->BulkFactor*(1-pt->beta_Gamma*mu);
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
	double nu_start_disk_RF;
	double nu_stop_disk_RF;
	double nuL_nu_disk,F_nu_disk_obs;
	//printf("=> Ciccio 1\n");
	if (pt->verbose){
		printf("-----------  Building I_nu disk     ----------- \n");
	}

	/*
	if (pt->WRITE_TO_FILE==1){
		sprintf(f_SED_disk, "%s%s-SED-disk.dat",pt->path, pt->STEM);

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
	if (pt->disk == 1)
	{
		nu_peak_BB=eval_nu_peak_Disk(pt->T_Disk);
		nu_start_disk_RF = nu_peak_BB*pt->nu_planck_min_factor;
		nu_stop_disk_RF  = nu_peak_BB*pt->nu_planck_max_factor;
	}
	else if (pt->disk == 2)
	{
		nu_peak_BB=eval_nu_peak_Disk(pt->T_Disk);
		nu_start_disk_RF = nu_peak_BB*pt->nu_planck_min_factor;
		nu_stop_disk_RF  = nu_peak_BB*pt->nu_planck_max_factor;
		//double (*pf) (struct spettro *, double x);
		//pf = &Disk_Spectrum;
		//pt->Cost_Norm_disk_Mulit_BB= 1.0/
		//		integrale_simp_struct(pf, pt,nu_start_disk_RF, nu_stop_disk_RF, pt->theta_n_int);
		//printf( "%e\n",pt->Cost_Norm_disk_Mulit_BB);
	 }
	 else if (pt->disk == 3)
	 {
		 nu_peak_BB = eval_nu_peak_Disk(pt->T_Disk);
		 nu_start_disk_RF = nu_peak_BB * pt->mono_planck_min_factor;
		 nu_stop_disk_RF = nu_peak_BB * pt->mono_planck_max_factor;
	}
	else{
		printf("wrong disk type, option BB, MultiBB, Mono \n ");
		exit(1);
	}

	pt->nu_start_Disk = eval_nu_min_blob_RF(pt, pt->Disk_mu_1, pt->Disk_mu_2, nu_start_disk_RF);
	pt->nu_stop_Disk = eval_nu_max_blob_RF(pt,pt->Disk_mu_1, pt->Disk_mu_2, nu_stop_disk_RF);

	pt->nu_start_Disk_DRF = nu_start_disk_RF;
	pt->nu_stop_Disk_DRF = nu_stop_disk_RF;

		if (pt->verbose)
	{
		printf("nu_start_Disk=%e  nu_stop_Disk=%e \n",
			pt->nu_start_Disk,
			pt->nu_stop_Disk);

		printf("nu_start_Disk_disk_RF=%e  nu_stop_Disk_disk_RF=%e \n",
			nu_start_disk_RF,
			nu_stop_disk_RF);
	}


	NU_INT_MAX=pt->nu_seed_size-1;
	pt->NU_INT_MAX_Disk = NU_INT_MAX;


	pt->nu_start_Disk_obs=nu_disk_to_nu_obs_disk(nu_start_disk_RF , pt->z_cosm);
	pt->nu_stop_Disk_obs=nu_disk_to_nu_obs_disk(nu_stop_disk_RF, pt->z_cosm);

	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->nu_seed_size, pt->nu_Disk_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->L_nu_Disk_disk_RF[NU_INT] = eval_Disk_L_nu(pt, pt->nu_Disk_disk_RF[NU_INT]);
		
	}
	for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++)
	{
		pt->I_nu_Disk_disk_RF[NU_INT] = eval_I_nu_Disk_disk_RF(pt, pt->nu_Disk_disk_RF[NU_INT]);		
	}

	build_log_grid( pt->nu_start_Disk,  pt->nu_stop_Disk, pt->nu_seed_size, pt->nu_Disk);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
 		nu_obs = nu_disk_to_nu_obs_disk(pt->nu_Disk_disk_RF[NU_INT],pt->z_cosm);
		pt->nu_Disk_obs[NU_INT]=nu_obs;
		pt->I_nu_Disk[NU_INT] = eval_I_nu_Disk_blob_RF(pt, pt->nu_Disk_disk_RF[NU_INT]);
		pt->n_Disk[NU_INT] =I_nu_to_n(pt->I_nu_Disk[NU_INT], pt->nu_Disk[NU_INT]);
		//EC with n(gamma) transf
		pt->n_Disk_DRF[NU_INT] = I_nu_to_n(pt->I_nu_Disk_disk_RF[NU_INT], pt->nu_Disk_disk_RF[NU_INT]);

		if (pt->verbose>1){
			printf(" nu_Disk_disk_RF=%e, I_nu_Disk_disk_RF=%e, nu_Disk=%e, , I_nu_Disk=%e\n",
				   pt->nu_Disk_disk_RF[NU_INT],
				   pt->I_nu_Disk_disk_RF[NU_INT],
				   pt->nu_Disk[NU_INT],
				   pt->I_nu_Disk[NU_INT]);
		}

		if (pt->I_nu_Disk[NU_INT]>pt->emiss_lim){
			pt->nu_stop_Disk = pt->nu_Disk[NU_INT];
			pt->NU_INT_MAX_Disk = NU_INT;
		}
		else{
			pt->I_nu_Disk[NU_INT]=pt->emiss_lim;
			pt->n_Disk[NU_INT] =I_nu_to_n(pt->I_nu_Disk[NU_INT], pt->nu_Disk[NU_INT]);

		}

		nuL_nu_disk = pt->L_nu_Disk_disk_RF[NU_INT] * pt->nu_Disk_disk_RF[NU_INT];
		F_nu_disk_obs= L_nu_Disk_to_F_nu(nuL_nu_disk / pt->nu_Disk_disk_RF[NU_INT], pt-> z_cosm, pt-> dist);
		pt->nuF_nu_Disk_obs[NU_INT] = F_nu_disk_obs*nu_obs;
		/*
		if (pt->WRITE_TO_FILE==1){
			fprintf(fp_SED_disk, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
				log10(nu_obs),
				log10(nu_obs * F_nu_disk_obs),
				nu_obs,
				nu_obs*F_nu_disk_obs,
				pt->nu_Disk_disk_RF[NU_INT],
				nuL_nu_disk);
		}
		*/


	}
	pt->L_Disk_radiative = PowerPhotons_disk_rest_frame(pt, pt->nu_Disk_disk_RF, pt->nuF_nu_Disk_obs, pt->NU_INT_MAX_Disk);
	
	/*
	if (pt->WRITE_TO_FILE==1){
		fclose(fp_SED_disk);
	}
	*/
}



void set_Disk(struct blob *pt){
	double  nu_peak_BB;
	if (strcmp(pt->disk_type, "BB") == 0)
	{
		pt->disk = 1;
	}
	else if (strcmp(pt->disk_type, "MultiBB") == 0)
	{
		pt->disk = 2;
	}
	else if (strcmp(pt->disk_type, "Mono") == 0)
	{
		pt->disk = 3;
	}
	else
	{
		printf("wrong disk type, option BB, MultiBB, Mono \n ");
		exit(1);
	}

	pt->R_Sw=eval_R_Sw(pt->M_BH);
	//R_inner
	pt->R_inner = pt->R_inner_Sw * pt->R_Sw;
	//R_ext
	pt->R_ext = pt->R_ext_Sw * pt->R_Sw;
	pt->R_Disk_interp =   pt->R_ext*50.0;
	pt->L_Edd = eval_L_Edd(pt->M_BH);
	pt->accr_rate = eval_accr_rate(pt->L_Disk, pt->accr_eff);
	pt->accr_Edd = eval_accr_Edd(pt->L_Edd, pt->accr_eff);
	//as in Ghisellini 2009, but it is equivalent to the one in Eq. 5.43 in the Frank, King & Raine Book
	//but we use 8p in place of 16pi, because L_Disk is the L of uno disk surface
	//this must be evaluated before eval_T_disk
	pt->Cost_disk_Mulit_BB = pt->R_inner * pt->L_Disk / (8 * pi * sigma_steph_boltz * pt->accr_eff);

	if (pt->disk == 2)
	//multi BB
	{
		pt->T_Disk = eval_T_disk(pt, (49. / 36.) * pt->R_inner);
		//printf("Cost_disk_Mulit_BB = %e \n", pt->Cost_disk_Mulit_BB);
	}
	
	nu_peak_BB = eval_nu_peak_Disk(pt->T_Disk);

	if (pt->verbose){
		printf("T_max = %e (K)\n",pt->T_Disk);
		// energy corresponding to Tmax
		printf("E_max = %e (eV)\n",pt->T_Disk*K_boltz*erg_to_eV);
		// frequency corresponding to Tmax
		printf("nu_max = %e (Hz)\n",pt->T_Disk*K_boltz/HPLANCK);
		//Peak of the BB spectrum
		printf("nu_peak  = %e (Hz)\n",nu_peak_BB);
		printf("schwarzschild radius=%e\n", pt->R_Sw);
		printf("R_ext =%e (cm)\n", pt->R_ext);
		printf("R_inner =%e (cm)\n", pt->R_inner);

		printf("Black hole mass = %e (m_sun)\n", pt->M_BH);

		printf("Accr. rate = %e (g/s)\n", pt->accr_rate);
		printf("Accr. rate = %e (M_sun/year)\n", pt->accr_rate * 86400. * 365. / m_sun);
		printf("L_Edd = %e (erg/s)\n", pt->L_Edd);
		printf("L_Disk = %e (erg/s)\n", pt->L_Disk);
		printf("L_diks/L_edd = %e\n", pt->L_Disk / pt->L_Edd);

		printf("Accr_Edd = %e (g/s)\n", pt->accr_Edd);
		printf("Accr_Edd = %e (M_sun/year)\n", pt->accr_Edd * 86400. * 365. / m_sun);
	}
}

//========================
// Disk Spectral Functions
//========================


double Disk_Spectrum(struct blob *pt, double nu_Disk_disk_RF){
	double I;
	double (*pf)(struct blob *, double x);
	I=0;
	if (pt->disk == 1) {
		// in this case we use a normalized planck function
		I= f_planck_norm(pt->T_Disk, nu_Disk_disk_RF);
	}
	else if (pt->disk == 2) {
		//in this case we acutally integrate every annluar BB along the disk
		
		pf = &integrand_f_planck_Multi_T;
		pt->nu_disk_Multi_BB = nu_Disk_disk_RF;
		//printf("=> pt->nu_disk_Multi_BB %e\n",pt->nu_disk_Multi_BB);
		//pi is the angular part for a disk face
		I= pi *integrale_trap_log_struct(pf, pt, pt->R_inner * 1.01, pt->R_ext, 100);
	}
	else if (pt->disk==3){
		I= eval_nu_peak_Disk(pt->T_Disk)*(pt->mono_planck_max_factor-pt->mono_planck_min_factor);
	}
	return I;
}

double eval_I_nu_theta_Disk(struct blob *pt, double mu)
{
	//double (*pf)(struct spettro *, double x);
	//unsigned int i;
	double  I,R,R_D;
	//pf = &j_nu_BLR_integrand;
	//pt->mu_j = mu;
	I=0;
    if (pt->disk == 1) {
		// in this case we use a normalized planck function
		I = f_planck_norm(pt->T_Disk, pt->nu_disk_RF)*pt->L_Disk * pt->Disk_geom_factor;
	}
	else if (pt->disk == 2) {
		//in this case we acutally integrate every annluar BB along the disk
		
	
		R=pt->R_H/mu;
		R_D = sqrt(R * R - pt->R_H * pt->R_H);
		I = f_planck_Multi_T(pt, R_D, pt->nu_disk_RF)/pi;
	}
	else if (pt->disk==3){
		I= eval_nu_peak_Disk(pt->T_Disk)*(pt->mono_planck_max_factor-pt->mono_planck_min_factor);
	}

	
	return I;
}

double integrand_I_nu_Disk_blob_RF(struct blob *pt, double mu)
{
	//double psi, sin_theta;
	//sin_theta=sqrt(1.0 - mu*mu);
	return 2 * pi  * eval_I_nu_theta_Disk(pt, mu) * pt->BulkFactor * (1.0 - pt->beta_Gamma * mu);
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
	pt->nu_disk_RF = nu_disk_RF;
	pf = &integrand_I_nu_Disk_blob_RF;

	c = 1.0;
	R_H_orig = pt->R_H;
	if (pt->R_H > pt->R_Disk_interp)
	{

		pt->R_H = pt->R_Disk_interp;
		c = (pt->R_Disk_interp / R_H_orig) * (pt->R_Disk_interp / R_H_orig);
	}

	set_Disk_angles(pt);
	I = integrale_simp_struct(pf, pt, pt->Disk_mu_1, pt->Disk_mu_2, pt->theta_n_int);
	pt->R_H = R_H_orig;
	set_Disk_angles(pt);
	return I * one_by_four_pi * c;

}

double eval_I_nu_Disk_disk_RF(struct blob *pt, double nu_disk_RF)
{
	double (*pf)(struct blob *, double x);
	double  I, R_H_orig, c;
	//unsigned int i;
	pt->nu_disk_RF = nu_disk_RF;
	pf = &integrand_I_nu_Disk_disk_RF;

	c = 1.0;
	R_H_orig = pt->R_H;
	if (pt->R_H > pt->R_Disk_interp)
	{

		pt->R_H = pt->R_Disk_interp;
		c = (pt->R_Disk_interp / R_H_orig) * (pt->R_Disk_interp / R_H_orig);
	}
	set_Disk_angles(pt);
	I = integrale_simp_struct(pf, pt, pt->Disk_mu_1, pt->Disk_mu_2, pt->theta_n_int);
	pt->R_H = R_H_orig;
	set_Disk_angles(pt);
	//printf("=> R_DT_interp=%e R_H=%e Disk_mu_1=%e Disk_mu_2=%e i=%e \n", pt->R_DT_interp, pt->R_H, pt->Disk_mu_1, pt->Disk_mu_2, I);
	return I * one_by_four_pi * c;
}

double eval_Disk_L_nu(struct blob *pt, double nu_Disk_disk_RF)
{
	if (pt->disk == 2) {
		//in this case no multiplication by L_Disk, because we acutally integrate every annluar BB along the disk
		//printf("=> %e\n", Disk_Spectrum(pt, nu_Disk_disk_RF));
		//printf("=> nu_Disk_disk_RF %e\n", nu_Disk_disk_RF);
		return  Disk_Spectrum(pt, nu_Disk_disk_RF);
	}
	else{
		return  pt->L_Disk *Disk_Spectrum(pt, nu_Disk_disk_RF);
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
	mu1 = pt->R_H / sqrt(pt->R_H * pt->R_H + pt->R_inner * pt->R_inner);
	mu2 = pt->R_H / sqrt(pt->R_H * pt->R_H + pt->R_ext * pt->R_ext);
	//mu1=1.0/sqrt(1+((pt->R_inner*pt->R_inner)/(pt->R_H*pt->R_H)));
	//mu2 = 1.0 / sqrt(1 + ((pt->R_ext * pt->R_ext) / (pt->R_H * pt->R_H)));
	pt->Disk_mu_1 = min(mu1, mu2);
	pt->Disk_mu_2 = max(mu1, mu2);
}

void set_Disk_geometry(struct blob *pt){

	pt->Disk_surface=pi*((pt->R_ext * pt->R_ext) - (pt->R_inner*pt->R_inner) );
	pt->Disk_geom_factor = (1.0) / (four_pi * pt->R_H * pt->R_H * (pt->Disk_surface / (pt->R_H * pt->R_H)));
}




//=========================================================================================
void Build_I_nu_BLR(struct blob *pt){

	//-------------------------------------
	// we follow the method in Donea&Protheroe https://arxiv.org/abs/astro-ph/0202068v1
	//-------------------------------------
	//double nu_stop_BLR_blob_RF;
	unsigned int NU_INT,NU_INT_MAX;
	double I_nu_theta_disk_RF,I_nu_theta_blob_RF;
	//char f_BLR_disk[static_file_name_max_legth];
	//FILE *fp_BLR_disk;

	/*
	if (pt->WRITE_TO_FILE==1){
		sprintf(f_BLR_disk, "%s%s-I_nu_BLR.dat",pt->path, pt->STEM);

		fp_BLR_disk = fopen(f_BLR_disk, "w");
		if (fp_BLR_disk == NULL) {
			printf("unable to open %s\n ", fp_BLR_disk);
			exit(1);
		}
	}
	*/
	//flux_DISK_header(fp_BLR_disk);
	if (pt->verbose){

		printf("-----------  Building I_nu BLR     ----------- \n");
	}
	set_BLR_geometry(pt);
	//printf("=>R_H=%e BLR_mu_1=%e BLR_mu_2=%e\n",pt->R_H,pt->BLR_mu_1,pt->BLR_mu_2);

	pt->BLR_mu_1 = 1.0;
	pt->BLR_mu_2 = cos(eval_theta_max_BLR(pt));
	
	//if (pt->tau_BLR>0.9){
	//	printf ("!!! Waring, the fraction of L_Disk reaching DT is (1-tau_BLR)\n");
	//	printf ("!!! if tau_BLR=1.0 no DT photons will be generated\n");

	//}

	pt->nu_start_BLR_disk_RF=pt->nu_Disk_disk_RF[0];
	pt->nu_stop_BLR_disk_RF=pt->nu_Disk_disk_RF[pt->NU_INT_MAX_Disk];

	pt->nu_start_BLR = eval_nu_max_blob_RF(pt, pt->BLR_mu_1, pt->BLR_mu_2, pt->nu_start_BLR_disk_RF);
	pt->nu_stop_BLR  = eval_nu_max_blob_RF(pt,pt->BLR_mu_1, pt->BLR_mu_2, pt->nu_stop_BLR_disk_RF);

	pt->R_BLR_interp_val = pt->R_BLR_out * 50.0;
	pt->R_BLR_interp_start = pt->R_BLR_out * 50.0;
	//printf("=>R_H=%e BLR_mu_1=%e BLR_mu_2=%e nu1=%e nu2=%e  nu1 d=%e nu2 d=%e\n", pt->R_H, pt->BLR_mu_1, pt->BLR_mu_2, pt->nu_start_BLR, pt->nu_stop_BLR, pt->nu_start_BLR_disk_RF, pt->nu_stop_BLR_disk_RF);
	pt->n0_BLR = pt->tau_BLR / (SIGTH * (pt->R_BLR_out - pt->R_BLR_in));
	

	if (pt->verbose)
	{
		printf("BLR_mu_1=%e BLR_mu_2=%e\n", pt->BLR_mu_1, pt->BLR_mu_2);

		printf("n0_BLR=%e \n", pt->n0_BLR);

		printf("nu_start_BLR_disk_RF=%e  nu_stop_BLR_disk_RF=%e \n",
				   pt->nu_start_BLR_disk_RF,
				   pt->nu_stop_BLR_disk_RF);

		printf("nu_start_BLR=%e  nu_stop_BLR=%e \n",
					pt->nu_start_BLR,
					pt->nu_stop_BLR);
	}


	NU_INT_MAX = pt->nu_seed_size-1;
	pt->NU_INT_MAX_BLR=NU_INT_MAX;
	//This is evaluating the angular pattern
	//It does not depends on frequency, because each region
	//is emitting the same spectrum
	I_nu_theta_disk_RF = eval_I_nu_BLR_disk_RF(pt);
	I_nu_theta_blob_RF = eval_I_nu_BLR_blob_RF(pt);

	build_log_grid( pt->nu_start_BLR_disk_RF,  pt->nu_stop_BLR_disk_RF, pt->nu_seed_size, pt->nu_BLR_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->Lnu_BLR_disk_RF[NU_INT] = eval_Lnu_BLR_disk_RF(pt, pt->nu_BLR_disk_RF[NU_INT]);
	}
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->I_nu_BLR_disk_RF[NU_INT] = I_nu_theta_disk_RF * pt->Lnu_BLR_disk_RF[NU_INT];
	}


	build_log_grid( pt->nu_start_BLR,  pt->nu_stop_BLR, pt->nu_seed_size, pt->nu_BLR);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		//we have to pass nu_BLR_disk_RF, because we integrate
		//the I' expressed in terms of I
		pt->I_nu_BLR[NU_INT] = I_nu_theta_blob_RF * pt->Lnu_BLR_disk_RF[NU_INT];
		pt->n_BLR[NU_INT] =I_nu_to_n(pt->I_nu_BLR[NU_INT], pt->nu_BLR[NU_INT]);
		//EC with n(gamma) transf
		pt->n_BLR_DRF[NU_INT] = I_nu_to_n(pt->I_nu_BLR_disk_RF[NU_INT], pt->nu_BLR_disk_RF[NU_INT]);

		if (pt->I_nu_BLR[NU_INT]>pt->emiss_lim){
			pt->nu_stop_BLR = pt->nu_BLR[NU_INT];
			pt->NU_INT_MAX_BLR = NU_INT;
		}
		else{
			pt->I_nu_BLR[NU_INT]=pt->emiss_lim;
			pt->n_BLR[NU_INT] =I_nu_to_n(pt->I_nu_BLR[NU_INT], pt->nu_BLR[NU_INT]);
		}

		if (pt->verbose>1){
			printf(" nu_BLR_disk_RF=%e, I_nu_BLR_disk_RF=%e, nu_BLR=%e, , I_nu_BLR=%e\n",
					pt->nu_BLR_disk_RF[NU_INT],
					pt->I_nu_BLR_disk_RF[NU_INT],
					pt->nu_BLR[NU_INT],
					pt->I_nu_BLR[NU_INT]);
		}
		/*
		if (pt->WRITE_TO_FILE==1){

			fprintf(fp_BLR_disk, "%4.4e\t %4.4e\t %4.4e\t %4.4e \n",
				log10(pt->nu_BLR_disk_RF[NU_INT]),
				log10(pt->I_nu_BLR_disk_RF[NU_INT]),
				log10(pt->nu_BLR[NU_INT]),
				log10(pt->I_nu_BLR[NU_INT]));
		}
		*/

	}
	/*
	if (pt->WRITE_TO_FILE == 1)
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

	//i = x_to_grid_index(pt->nu_BLR_disk_RF, pt->nu_disk_RF, pt->nu_seed_size);
	
	r2 = (pt->R_H * pt->R_H) - 2.0 * pt->R_H * l * pt->mu_j + l * l;
	
	
	//L = eval_Disk_L_nu(pt, pt->nu_disk_RF) * pt->n0_BLR * SIGTH;
	if ((r2 > (pt->R_BLR_out * pt->R_BLR_out)) || (r2 < (pt->R_BLR_in * pt->R_BLR_in)))
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
	pt->mu_j=mu;
	
	eval_l_values_BLR(pt, mu, l_values);
	if(pt->R_H<pt->R_BLR_out){
		I = integrale_simp_struct(pf, pt, 0, l_values[0], pt->l_n_int)+ integrale_simp_struct(pf, pt, l_values[1], l_values[2], pt->l_n_int);
	}
	else{
		I = integrale_simp_struct(pf, pt, 0, l_values[2], pt->l_n_int);
	}
	//printf("mu=%e, l0=%e, l1=%e, l2=%e, delta=%e\n", mu, l_values[0], l_values[1], l_values[2], l_values[2]- l_values[1]);
	return I;
}

double integrand_I_nu_BLR_blob_RF(struct blob *pt, double theta)
{
	//double psi
	//double mu,mu1,c;
	//mu = cos(theta);
	//mu1 = (pt->beta_Gamma - mu )/(pt->beta_Gamma*mu - 1.0 );
	//c=(pt->BulkFactor * pt->BulkFactor * pt->BulkFactor );
	//c = c * (1.0 + pt->BulkFactor * mu + 1.0) * (1.0 + pt->BulkFactor * mu + 1.0) * (1.0 + pt->BulkFactor * mu + 1.0);

	return 2 * pi * sin(theta) * eval_I_nu_theta_BLR(pt, cos(theta)) * pt->BulkFactor * (1.0 - pt->beta_Gamma * cos(theta));
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

	//pt->nu_disk_RF=nu_disk_RF;
	pf = &integrand_I_nu_BLR_disk_RF;
	c=1.0;
	R_H_orig = pt->R_H;
	if (pt->R_H > pt->R_BLR_interp_start)
	{
		
		pt->R_H = pt->R_BLR_interp_val;
		c = (pt->R_BLR_interp_val / R_H_orig) * (pt->R_BLR_interp_val / R_H_orig);
		//printf("=>R_H=%e R_H_orig=%e  pt->R_BLR_interp=%e\n",pt->R_H,R_H_orig,pt->R_BLR_interp);
	}
	theta_min=0.0;
	theta_max = eval_theta_max_BLR(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->theta_n_int);
	pt->R_H = R_H_orig;
	//printf("=>R_H=%e R_BLR_inter=%e I=%e %e %e c=%e\n ",pt->R_H,pt->R_BLR_interp, I, theta_min, theta_max,c);
	return I*one_by_four_pi*c;
}


double eval_I_nu_BLR_blob_RF(struct blob *pt)
{
	double (*pf)(struct blob *, double x);
	double theta_min, theta_max, I, R_H_orig,c;
	// we use directly nu_disk_RF
	// because we integrate the I' expressed as I
	//pt->nu_disk_RF = nu_disk_RF;
	pf = &integrand_I_nu_BLR_blob_RF;
	
	c=1.0;
	R_H_orig = pt->R_H;
	if (pt->R_H > pt->R_BLR_interp_start)
	{

		pt->R_H = pt->R_BLR_interp_val;
		c = (pt->R_BLR_interp_val / R_H_orig) * (pt->R_BLR_interp_val / R_H_orig);
		//printf("=>R_H=%e R_H_orig=%e  pt->R_BLR_interp=%e\n",pt->R_H,R_H_orig,pt->R_BLR_interp);
	}
	theta_min = 0.0;
	theta_max = eval_theta_max_BLR(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->theta_n_int);
	pt->R_H = R_H_orig;
	//printf("=>BLR R_H=%e R_B=%e I=%e %e %e c=%e\n ", pt->R_H, pt->R_BLR_out, I, theta_min, theta_max, c);
	return I*one_by_four_pi*c;
}

double eval_Lnu_BLR_disk_RF(struct blob *pt, double nu_disk_RF)
{
	return eval_Disk_L_nu(pt, nu_disk_RF) * pt->n0_BLR *SIGTH;
}




//========================
// BLR Geometrical Functions
//========================

double eval_theta_max_BLR(struct blob *pt)
{
	double theta_max;
	if (pt->R_H > pt->R_BLR_out)
	{
		theta_max = asin(pt->R_BLR_out / pt->R_H);
		//theta_max = 2 * (pi * 0.5 - acos(pt->R_BLR_out / pt->R_H))	;
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

		
	s = mu * mu + (pt->R_BLR_in / pt->R_H) * (pt->R_BLR_in / pt->R_H) - 1.0;
	if (s < 0.0){
		l[0] = 0.0;
		l[1] = 0.0;
	}else
	{
		l[1] = pt->R_H * mu + pt->R_H * sqrt(s);
		l[0]= pt->R_H * mu - pt->R_H * sqrt(s);
	}
	if (l[1] < 0.0){
		l[1] = 0.;
	}

	if (l[0] < 0.0){
		l[0] = 0.;
	}

	s = mu * mu + (pt->R_BLR_out / pt->R_H) * (pt->R_BLR_out / pt->R_H) - 1.0;
	if (s < 0.0){
		l[2] = 0;
	}
	else{

		l[2] = pt->R_H * mu + pt->R_H * sqrt(s);
		
		if (l[2] < 0.0){
			l[2] = 0.;
		}
	}
}

void set_BLR_geometry(struct blob *pt)
{

	pt->BLR_Volume = (4. / 3.) * pi * ((pt->R_BLR_out * pt->R_BLR_out * pt->R_BLR_out) - (pt->R_BLR_in * pt->R_BLR_in * pt->R_BLR_in));
	pt->BLR_inner_Surface = 4 * pi * (pt->R_BLR_in * pt->R_BLR_in);
	pt->Delta_R_BLR = pt->R_BLR_out - pt->R_BLR_in;

	

	/*
	if (pt->R_H < pt->R_BLR_in)
	{
		pt->BLR_mu_1 = -1.0;
		pt->BLR_mu_2 = pt->R_H / sqrt(pt->R_H * pt->R_H + pt->R_ext * pt->R_ext);
		//pt->BLR_mu_r_J_2 = pt->R_H / sqrt(pt->R_H * pt->R_H + pt->R_ext * pt->R_ext);
		//pt->BLR_mu_r_J_1 = -1.0;
		//pt->BLR_geom_factor=(1.0)/(4*pi*4*pi);
		//pt->BLR_geom_factor*=1.0/(pt->R_BLR_in*pt->R_BLR_in);
	}
	else if (pt->R_H >= pt->R_BLR_in && pt->R_H < pt->R_BLR_out)
	{
		pt->BLR_mu_1 = -1.0;
		pt->BLR_mu_2 = pt->R_H / sqrt(pt->R_H * pt->R_H + pt->R_ext * pt->R_ext);
		//pt->BLR_mu_r_in = sqrt(pt->R_H * pt->R_H - pt->R_BLR_in * pt->R_BLR_in) / pt->R_H;
		//pt->BLR_mu_r_J_2 = pt->R_H / sqrt(pt->R_H * pt->R_H + pt->R_ext * pt->R_ext);
		//pt->BLR_mu_r_J_1 = -1.0;
	}
	else
	{
		pt->BLR_mu_2 = pt->R_H / sqrt(pt->R_H * pt->R_H + pt->R_ext * pt->R_ext);
		pt->BLR_mu_1 = sqrt(pt->R_H * pt->R_H - pt->R_BLR_out * pt->R_BLR_out) / pt->R_H;
		//pt->BLR_mu_r_in = sqrt(pt->R_H * pt->R_H - pt->R_BLR_in * pt->R_BLR_in) / pt->R_H;
		//pt->BLR_mu_r_J_1 = pt->BLR_mu_1;
		//pt->BLR_mu_r_J_2 = pt->BLR_mu_2;
		//pt->BLR_geom_factor=(1.0)/(4*pi*4*pi);
		//pt->BLR_geom_factor*=1.0/(pt->R_BLR_in*pt->R_BLR_in);
	}
	*/
}

//=========================================================================================

void Build_I_nu_DT(struct blob *pt){
	//FILE *fp_SED_DT;
	//char f_SED_DT[static_file_name_max_legth];
	unsigned int NU_INT,NU_INT_MAX;
	double I_nu_theta_disk_RF, I_nu_theta_blob_RF;
	double nu_peak_DT_disk_RF;
	double nu_start_DT_disk_RF,nu_stop_DT_disk_RF;
	double nu_obs;
	double nuL_nu_DT,F_nu_DT_obs;

	/*
	if (pt->WRITE_TO_FILE==1){
		sprintf(f_SED_DT, "%s%s-SED-DT.dat",
					pt->path, pt->STEM);

		fp_SED_DT = fopen(f_SED_DT, "w");
		if (fp_SED_DT == NULL) {
			printf("unable to open %s\n ", f_SED_DT);
			exit(1);
		}
		flux_DISK_header(fp_SED_DT);
	}
	*/

	if (pt->verbose){

		printf("-----------  Building I_nu DT     ----------- \n");
	}


	//if (pt->tau_BLR>0.9){
	//		printf ("!!! Waring, the fraction of L_Disk reaching DT is (1-tau_BLR)\n");
	//		printf ("!!! if tau_BLR=1.0 no DT photons will be generated\n");

	//}

	nu_peak_DT_disk_RF=eval_nu_peak_planck(pt->T_DT);

	nu_start_DT_disk_RF=nu_peak_DT_disk_RF*pt->nu_planck_min_factor;
	nu_stop_DT_disk_RF=nu_peak_DT_disk_RF*pt->nu_planck_max_factor;


	pt->DT_mu_1 = 1.0;
	pt->DT_mu_2 = cos(eval_theta_max_DT(pt));

	pt->nu_start_DT = eval_nu_max_blob_RF(pt, pt->DT_mu_1, pt->DT_mu_2, nu_start_DT_disk_RF);
	pt->nu_stop_DT  = eval_nu_max_blob_RF(pt,pt->DT_mu_1, pt->DT_mu_2, nu_stop_DT_disk_RF);

	pt->nu_start_DT_DRF=nu_start_DT_disk_RF;
	pt->nu_stop_DT_DRF=nu_stop_DT_disk_RF;

	pt->nu_start_DT_obs=nu_disk_to_nu_obs_disk(nu_start_DT_disk_RF , pt->z_cosm);
	pt->nu_stop_DT_obs=nu_disk_to_nu_obs_disk(nu_stop_DT_disk_RF, pt->z_cosm);

	if (pt->verbose){
		printf("nu_start_DT (blob frame) =%e \n",
					pt->nu_start_DT);
		printf("nu_stop_DT (blob frame) =%e \n",
						pt->nu_stop_DT);
		printf("nu_start_DT (disk frame) =%e \n",
				nu_start_DT_disk_RF);
		printf("nu_stop_DT (disk frame) =%e \n",
				nu_stop_DT_disk_RF);
	}

	NU_INT_MAX = pt->nu_seed_size-1;
	pt->NU_INT_MAX_DT = NU_INT_MAX;

	pt->R_DT_interp_val = pt->R_DT* 50.0;
	pt->R_DT_interp_start = pt->R_DT * 50.0;

	pt->DT_Volume=(4./3.)*pi*pt->R_DT*pt->R_DT*pt->R_DT;

	//This is evaluating the angular pattern
	//It does not depends on frequency, because each region
	//is emitting the same spectrum
	I_nu_theta_disk_RF = eval_I_nu_DT_disk_RF(pt);
	I_nu_theta_blob_RF = eval_I_nu_DT_blob_RF(pt);

	build_log_grid( nu_start_DT_disk_RF,  nu_stop_DT_disk_RF, pt->nu_seed_size, pt->nu_DT_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->L_nu_DT_disk_RF[NU_INT] = eval_DT_L_nu(pt, pt->nu_DT_disk_RF[NU_INT]);
	}
	for (NU_INT = 0; NU_INT <= NU_INT_MAX; NU_INT++){
		pt->I_nu_DT_disk_RF[NU_INT] = I_nu_theta_disk_RF * pt->L_nu_DT_disk_RF[NU_INT];
	}


	build_log_grid( pt->nu_start_DT,  pt->nu_stop_DT, pt->nu_seed_size, pt->nu_DT);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {


		nu_obs = nu_disk_to_nu_obs_disk(pt->nu_DT_disk_RF[NU_INT], pt->z_cosm);

		pt->nu_DT_obs[NU_INT]=nu_obs;

		pt->I_nu_DT[NU_INT] = I_nu_theta_blob_RF * pt->L_nu_DT_disk_RF[NU_INT];
		pt->n_DT[NU_INT] =I_nu_to_n(pt->I_nu_DT[NU_INT], pt->nu_DT[NU_INT]);
		//EC with n(gamma) transf
		pt->n_DT_DRF[NU_INT] = I_nu_to_n(pt->I_nu_DT_disk_RF[NU_INT], pt->nu_DT_disk_RF[NU_INT]);

		if (pt->I_nu_DT[NU_INT]>pt->emiss_lim){
			pt->nu_stop_DT = pt->nu_DT[NU_INT];
			pt->NU_INT_MAX_DT = NU_INT;
		}
		else{
			pt->I_nu_DT[NU_INT]=pt->emiss_lim;
			pt->n_DT[NU_INT] =I_nu_to_n(pt->I_nu_DT[NU_INT], pt->nu_DT[NU_INT]);
		}

		nuL_nu_DT = pt->L_nu_DT_disk_RF[NU_INT] * pt->nu_DT_disk_RF[NU_INT];
		//if (pt->tau_BLR<1.0){
		//	nuL_nu_DT = pt->L_nu_DT_disk_RF[NU_INT] * pt->nu_DT_disk_RF[NU_INT];
		//}
		//else {
		//	nuL_nu_DT = 0;
		//}
		F_nu_DT_obs = L_nu_Disk_to_F_nu(nuL_nu_DT / pt->nu_DT_disk_RF[NU_INT], pt-> z_cosm, pt-> dist);

		pt->nuF_nu_DT_obs[NU_INT] = F_nu_DT_obs*nu_obs;
		if (pt->verbose>1){
			printf(" nu_DT_disk_RF=%e, I_nu_DT_disk_RF=%e, nu_DT=%e, I_nu_DT=%e\n",
				pt->nu_DT_disk_RF[NU_INT],
				pt->I_nu_DT_disk_RF[NU_INT],
				pt->nu_DT[NU_INT],
				pt->I_nu_DT[NU_INT]);
		}
		/*
		if (pt->WRITE_TO_FILE==1){
			fprintf(fp_SED_DT, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
									log10(nu_obs),
									log10(nu_obs * F_nu_DT_obs),
									nu_obs,
									nu_obs*F_nu_DT_obs,
									pt->nu_DT_disk_RF[NU_INT],
									nuL_nu_DT);
		}
		*/
	}
	/*
	if (pt->WRITE_TO_FILE==1){
		fclose(fp_SED_DT);
	}
	*/
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->n_DT[NU_INT] =I_nu_to_n(pt->I_nu_DT[NU_INT], pt->nu_DT[NU_INT]);
	}
}

//========================
// Torus Spectral Functions
//========================

double j_nu_DT_integrand(struct blob *pt, double l)
{
	//unsigned int i;
	double L, r2;

	//i = x_to_grid_index(pt->nu_BLR_disk_RF, pt->nu_disk_RF, pt->nu_seed_size);
	
	r2 = (pt->R_H * pt->R_H) - 2.0 * pt->R_DT * l * pt->mu_j + l * l;
	
	
	//L = eval_Disk_L_nu(pt, pt->nu_disk_RF) * pt->n0_BLR * SIGTH;
	if  (r2 > (pt->R_DT * pt->R_DT) )
	{
		L=0.0;
	}
	else{
		L =1.0/ (four_pi * four_pi * r2*pt->R_DT);
	}
	return L;
}


double eval_I_nu_theta_DT(struct blob *pt, double mu, double theta)
{
	//double (*pf)(struct spettro *, double x);
	//unsigned int i;
	double l, I,cos_theta_norm,alpha;
	//double I;
	
	
	//i = x_to_grid_index(pt->nu_DT_disk_RF, pt->nu_disk_RF, pt->nu_seed_size);
	if (pt->R_H < pt->R_DT)
	{
		I = 1.0 / (4 * pi * 4 * pi * pt->R_DT * pt->R_DT);
	}	
	else
	{
		l = eval_l_DT(pt, mu);
		//pf = &j_nu_DT_integrand;
		//pt->mu_j = mu;

		
		alpha = acos(l * sin(theta) / pt->R_DT);

		cos_theta_norm = cos(pi - (alpha + 0.5 * pi - theta ));

		I = cos_theta_norm	 / ((4 * pi * pi * pt->R_H * pt->R_H) * ((pt->R_DT / pt->R_H) * (pt->R_DT / pt->R_H)));

		//I = integrale_simp_struct(pf, pt, 0, l, pt->l_n_int);

	}
			
	//I = integrale_simp_struct(pf, pt, 0, l, pt->l_n_int);
	return I;
}

double integrand_I_nu_DT_blob_RF(struct blob *pt, double theta)
{
	//double psi;
	return 2 * pi * sin(theta) * eval_I_nu_theta_DT(pt, cos(theta), theta) * pt->BulkFactor * (1.0 - pt->beta_Gamma * cos(theta));
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
	//pt->nu_disk_RF = nu_disk_RF;

	pf = &integrand_I_nu_DT_disk_RF;

	
	c = 1.0;
	R_H_orig = pt->R_H;
	if (pt->R_H > pt->R_DT_interp_start)
	{

		pt->R_H = pt->R_DT_interp_val;
		c = (pt->R_DT_interp_val / R_H_orig) * (pt->R_DT_interp_val / R_H_orig);
	}
	theta_min = 0.0;
	theta_max = eval_theta_max_DT(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->theta_n_int);
	pt->R_H = R_H_orig;
	return I * one_by_four_pi * c;
}

double eval_I_nu_DT_blob_RF(struct blob *pt )
{
	double (*pf)(struct blob *, double x);
	double theta_min, theta_max, I, R_H_orig, c;

	//now integrating only over angles
	//not needed anymore
	//pt->nu_disk_RF = nu_disk_RF;

	pf = &integrand_I_nu_DT_blob_RF;

	c=1.0;
	R_H_orig = pt->R_H;
	if (pt->R_H > pt->R_DT_interp_start)
	{

		pt->R_H = pt->R_DT_interp_val;
		c = (pt->R_DT_interp_val / R_H_orig) * (pt->R_DT_interp_val / R_H_orig);
	}
	theta_min = 0.0;
	theta_max = eval_theta_max_DT(pt);

	I = integrale_simp_struct(pf, pt, theta_min, theta_max, pt->theta_n_int);
	pt->R_H = R_H_orig;
	//printf("=>DT  R_H=%e R_D=%e I=%e %e %e c=%e\n ", pt->R_H, pt->R_DT, I, theta_min, theta_max, c);
	return I * one_by_four_pi * c;
}

double eval_DT_L_nu(struct blob *pt, double DT_disk_RF)
{
	//unsigned int i;
	//i = x_to_grid_index(pt->nu_Disk_disk_RF, DT_disk_RF, pt->nu_seed_size);
	return pt->L_Disk_radiative * pt->tau_DT * f_planck_norm(pt->T_DT, DT_disk_RF);
}

//========================
// Torus Geometrical Functions
//========================
double eval_theta_max_DT(struct blob *pt)
{
	double theta_max;

	if (pt->R_H >= pt->R_DT){
		theta_max= asin(pt->R_DT / pt->R_H);
	}
	else{
		theta_max=pi;
	}

	return theta_max;
}

double eval_l_DT(struct blob *pt, double mu)
{
	double s,l;

	s = mu * mu + (pt->R_DT / pt->R_H) * (pt->R_DT / pt->R_H) - 1.0;
	if (s < 0.0){
		l=0.0;
	}
	else{
		l= pt->R_H * mu - pt->R_H * sqrt(s);
	}
	if (l < 0.0){
		l = 0.;
	}
	return l;
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
	return 4*pi*G_cgs*M_BH*m_sun*vluce_cm/0.3;
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
	T_disco_r = pt->Cost_disk_Mulit_BB/(R*R*R) * (1 - pow((pt->R_inner / R), 0.5));
	T_disco_r = pow(T_disco_r, 0.25);
	//printf("=> T_disco_r %e\n",T_disco_r);
	return T_disco_r;
}

double f_planck_Multi_T(struct blob *pt, double R ,double nu) {
	if (R>pt->R_ext || R<pt->R_inner) {
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
	//printf("=> %e %e %e\n", f_planck_Multi_T(pt, R, pt->nu_disk_Multi_BB), R, pt->nu_disk_Multi_BB);
	return 2*pi*f_planck_Multi_T(pt,R,pt->nu_disk_Multi_BB)*R;
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
	nu_1=nu_disk_RF*pt->BulkFactor*(1-pt->beta_Gamma*mu1);
	nu_2=nu_disk_RF*pt->BulkFactor*(1-pt->beta_Gamma*mu2);
	a=  min(nu_1,nu_2);
	return a;
}



double eval_nu_max_blob_RF(struct blob *pt, double mu1, double mu2, double nu_disk_RF ){
	double a,nu_1,nu_2;
	nu_1=nu_disk_RF*pt->BulkFactor*(1-pt->beta_Gamma*mu1);
	nu_2=nu_disk_RF*pt->BulkFactor*(1-pt->beta_Gamma*mu2);
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
