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

void spectra_External_Fields(int Num_file, struct spettro *pt) {

    //==================================================================
	if (pt->verbose){
		printf("-----------  Eval. seed photon fields for  EC  ----------- \n");
	}


    // ====================================================
    // approx gamma con beaming factor, impilcit teta circa= 1/gamma
    // Set nu start EC
	// not used in photon field computation
	// used only in analytic approx on screen
    //=====================================================
	pt->beaming_EC = pt->BulkFactor;



    //T
    
    if (pt->do_EC_Star==1){
    	Build_I_nu_Star(pt);
    }
    if (pt->do_EC_Disk==1 || pt->do_EC_BLR==1){
    	Build_I_nu_Disk(pt);
    }
    if (pt->do_EC_BLR==1){
		Build_I_nu_BLR(pt);
	}
    if (pt->do_EC_DT==1){
    	Build_I_nu_DT(pt);
    }
    if (pt->do_EC_CMB==1){
    	Build_I_nu_CMB(pt);
    }

    if (pt->do_EC_CMB_stat==1){
    	Build_I_nu_CMB_stat(pt);
    }

}
//=========================================================================================


//=========================================================================================
void Build_I_nu_Star(struct spettro *pt){
	char f_SED_star[512];
	FILE *fp_SED_star;
	double nu_peak_BB,nu_obs;
	unsigned long NU_INT,NU_INT_MAX;
	double nu_start_disk_RF,nu_start_blob_RF;
	double nu_stop_disk_RF,nu_stop_blob_RF;
	double nuL_nu_disk,F_nu_disk_obs;

	sprintf(f_SED_star, "%s%s-SED-star.dat",pt->path, pt->STEM);

	fp_SED_star = fopen(f_SED_star, "w");
	if (fp_SED_star == NULL) {
		printf("unable to open %s\n ", fp_SED_star);
		exit(1);
	}
	flux_DISK_header(fp_SED_star);



	set_Star_geometry(pt);

	nu_peak_BB=eval_nu_peak_Disk(pt->T_Star_max);

	nu_start_disk_RF = nu_peak_BB*pt->nu_planck_min_factor;
	nu_stop_disk_RF  = nu_peak_BB*pt->nu_planck_max_factor;

	pt->nu_start_Star = eval_nu_min_blob_RF(pt,pt->Star_mu_1, pt->Star_mu_2, nu_start_disk_RF);
	pt->nu_stop_Star  = eval_nu_max_blob_RF(pt,pt->Star_mu_1, pt->Star_mu_2, nu_stop_disk_RF);

	if (pt->verbose){
		printf("-----------  Building I_nu Star     ----------- \n");

		printf("nu_start_Star=%e  nu_stop_Star=%e \n",
			pt->nu_start_Star,
			pt->nu_stop_Star);

		printf("nu_start_Star_disk_RF=%e  nu_stop_Star_disk_RF=%e \n",
			nu_start_disk_RF,
			nu_stop_disk_RF);
	}


	NU_INT_MAX=pt->nu_seed_size-1;
	pt->NU_INT_MAX_Star = NU_INT_MAX;


	pt->nu_start_Star_obs=nu_disk_to_nu_obs_disk(nu_start_disk_RF , pt->z_cosm);
	pt->nu_stop_Star_obs=nu_disk_to_nu_obs_disk(nu_start_disk_RF, pt->z_cosm);

	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->nu_seed_size, pt->nu_Star_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->I_nu_Star_disk_RF[NU_INT]=eval_I_nu_Star_disk_RF(pt, pt->nu_Star_disk_RF[NU_INT]);
		pt->J_nu_Star_disk_RF[NU_INT]=eval_J_nu_Star_disk_RF(pt, pt->I_nu_Star_disk_RF[NU_INT]);

	}



	build_log_grid( pt->nu_start_Star,  pt->nu_stop_Star, pt->nu_seed_size, pt->nu_Star);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		nu_obs = nu_disk_to_nu_obs_disk(pt->nu_Star_disk_RF[NU_INT],pt->z_cosm);
		pt->nu_Star_obs[NU_INT]=nu_obs;

		pt->I_nu_Star[NU_INT]=eval_I_nu_Star_blob_RF(pt,pt->nu_Star[NU_INT]);
		pt->n_Star[NU_INT] =I_nu_to_n(pt->I_nu_Star[NU_INT], pt->nu_Star[NU_INT]);
		if (pt->verbose>1){
			printf(" nu_Star_disk_RF=%e, I_nu_Star_disk_RF=%e, nu_Star=%e, , I_nu_Star=%e\n",
					pt->nu_Star_disk_RF[NU_INT],
					pt->J_nu_Star_disk_RF[NU_INT],
					pt->nu_Star[NU_INT],
					pt->I_nu_Star[NU_INT]);
		}

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

		fprintf(fp_SED_star, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
			log10(nu_obs),
			log10(nu_obs * F_nu_disk_obs),
			nu_obs,
			nu_obs*F_nu_disk_obs,
			pt->nu_Star_disk_RF[NU_INT],
			nuL_nu_disk);



	}

	fclose(fp_SED_star);

}


//========================
// Star Spectral Functions
//========================

double eval_I_nu_Star_disk_RF(struct spettro *pt,double nu_Star_disk_RF){
	return f_planck(pt->T_Star_max, nu_Star_disk_RF);;
}

double eval_J_nu_Star_disk_RF(struct spettro *pt, double I_nu_Star_disk_RF){
	return I_nu_Star_disk_RF*(2.0*pi/(4*pi))*(pt->Star_mu_2- pt->Star_mu_1);
}


double integrand_I_nu_Star_blob_RF(struct spettro *pt, double mu){
	unsigned long i;
	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->nu_blob_RF,pt->BulkFactor,pt->beta_Gamma,mu);

	i=x_to_grid_index( pt->nu_Star_disk_RF,nu_disk_RF,pt->nu_seed_size);
	if (i>0){
		return pt->J_nu_Star_disk_RF[i]*pt->BulkFactor*(1-pt->beta_Gamma*mu);
	}
	else{
		return 0;
	}
}



double eval_I_nu_Star_blob_RF(struct spettro *pt, double nu_blob_RF){
	unsigned long i;
	pt->nu_blob_RF=nu_blob_RF;
	double (*pf) (struct spettro *, double x);
	pf = &integrand_I_nu_Star_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return  0.5*integrale_simp_struct(pf, pt, pt->Star_mu_1, pt->Star_mu_2,100.0);
}

double eval_Star_L_nu(struct spettro *pt, double nu_Star_disk_RF){
	return  pt->Star_surface *eval_I_nu_Star_disk_RF(pt, nu_Star_disk_RF)*pi;
}


//========================
// Star Geometrical Functions
//========================

void set_Star_geometry(struct spettro *pt){
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
void Build_I_nu_CMB(struct spettro *pt){
	double T_CMB_z,T_CMB_0;
	double nu_peak_CMB_z,nu_peak_CMB_0;
	unsigned long NU_INT,NU_INT_MAX;
	double nu_start_disk_RF;
	double nu_stop_disk_RF;

	pt->CMB_mu_1=-1.0;
	pt->CMB_mu_2=1.0;

	T_CMB_z=eval_T_CMB_z(pt->z_cosm,pt->T_CMB_0);
	T_CMB_0=pt->T_CMB_0;

	nu_peak_CMB_z=eval_nu_peak_planck(T_CMB_z);
	nu_peak_CMB_0=eval_nu_peak_planck(T_CMB_0);

	nu_start_disk_RF = nu_peak_CMB_z*pt->nu_planck_min_factor;
	nu_stop_disk_RF  = nu_peak_CMB_z*pt->nu_planck_max_factor;

	pt->nu_start_CMB = eval_nu_min_blob_RF(pt,-1, 1, nu_start_disk_RF);
	pt->nu_stop_CMB  = eval_nu_max_blob_RF(pt,-1, 1, nu_stop_disk_RF);

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

	}

}




void Build_I_nu_CMB_stat(struct spettro *pt){
	double T_CMB_z,T_CMB_0;
	double nu_peak_CMB_z,nu_peak_CMB_0;
	unsigned long NU_INT,NU_INT_MAX;
	double nu_start_disk_RF;
	double nu_stop_disk_RF;

	pt->CMB_mu_1=-1.0;
	pt->CMB_mu_2=1.0;

	T_CMB_z=eval_T_CMB_z(pt->z_cosm,pt->T_CMB_0);
	T_CMB_0=pt->T_CMB_0;

	nu_peak_CMB_z=eval_nu_peak_planck(T_CMB_z);
	nu_peak_CMB_0=eval_nu_peak_planck(T_CMB_0);

	nu_start_disk_RF = nu_peak_CMB_z*pt->nu_planck_min_factor;
	nu_stop_disk_RF  = nu_peak_CMB_z*pt->nu_planck_max_factor;

	pt->nu_start_CMB_stat = nu_start_disk_RF;
	pt->nu_stop_CMB_stat  = nu_stop_disk_RF;

    //printf("AZZo nu_star_CMB=%e    nu_stop_CMB=%e\n",
    //                pt->nu_start_CMB_stat,
    //                pt->nu_stop_CMB_stat);

	//pt->nu_start_CMB_obs=nu_peak_CMB_0*pt->nu_planck_min_factor;
	//pt->nu_stop_CMB_obs=nu_peak_CMB_0*pt->nu_planck_max_factor;

	NU_INT_MAX=pt->nu_seed_size-1;
	pt->NU_INT_MAX_CMB = NU_INT_MAX;

	build_log_grid( nu_start_disk_RF,  nu_stop_disk_RF, pt->nu_seed_size, pt->nu_CMB_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
			//printf("AZZo nu_CMB_disk_RF=%e\n",pt->nu_CMB_disk_RF[NU_INT]);
			pt->I_nu_CMB_disk_RF[NU_INT]=eval_I_nu_CMB_disk_RF(T_CMB_z, pt->nu_CMB_disk_RF[NU_INT]);
	}

	build_log_grid( pt->nu_start_CMB_stat,  pt->nu_stop_CMB_stat, pt->nu_seed_size, pt->nu_CMB_stat);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		//printf("AZZo nu_CMB_stat=%e\n",pt->nu_CMB_stat[NU_INT]);
		pt->I_nu_CMB_stat[NU_INT]=pt->I_nu_CMB_disk_RF[NU_INT];
		pt->n_CMB_stat[NU_INT] =I_nu_to_n(pt->I_nu_CMB_stat[NU_INT], pt->nu_CMB_stat[NU_INT]);

	}

}



double eval_T_CMB_z(double z, double T_CMB_0){
		return T_CMB_0*(1+z);
}

double eval_I_nu_CMB_disk_RF(double T_CMB,double nu_CMB_disk_RF){
	return f_planck(T_CMB, nu_CMB_disk_RF);
}


double eval_I_nu_CMB_blob_RF(struct spettro *pt, double nu_blob_RF){
	unsigned long i;
	pt->nu_blob_RF=nu_blob_RF;
	double (*pf) (struct spettro *, double x);
	pf = &integrand_I_nu_CMB_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return  0.5*integrale_simp_struct(pf, pt, pt->CMB_mu_1, pt->CMB_mu_2,100.0);
}

double integrand_I_nu_CMB_blob_RF(struct spettro *pt, double mu){
	unsigned long i=0;
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
void Build_I_nu_Disk(struct spettro *pt){

	char f_SED_disk[512];
	FILE *fp_SED_disk;
	double nu_peak_BB,nu_obs;
	unsigned long NU_INT,NU_INT_MAX;
	double nu_start_disk_RF,nu_start_blob_RF;
	double nu_stop_disk_RF,nu_stop_blob_RF;
	double nuL_nu_disk,F_nu_disk_obs;

	if (pt->verbose){
		printf("-----------  Building I_nu disk     ----------- \n");
	}

	sprintf(f_SED_disk, "%s%s-SED-disk.dat",pt->path, pt->STEM);

	fp_SED_disk = fopen(f_SED_disk, "w");
	if (fp_SED_disk == NULL) {
		printf("unable to open %s\n ", f_SED_disk);
		exit(1);
	}
	flux_DISK_header(fp_SED_disk);


	set_Disk(pt);
	set_Disk_geometry(pt);

	if (strcmp(pt->disk_type, "BB") == 0) {
		pt->disk = 1;
		nu_peak_BB=eval_nu_peak_Disk(pt->T_disk_max);
		nu_start_disk_RF = nu_peak_BB*pt->nu_planck_min_factor;
		nu_stop_disk_RF  = nu_peak_BB*pt->nu_planck_max_factor;
	}
	else if (strcmp(pt->disk_type, "MultiBB") == 0) {
		pt->disk = 2;
		nu_peak_BB=eval_nu_peak_Disk(pt->T_disk_max);
		nu_start_disk_RF = nu_peak_BB*pt->nu_planck_min_factor;
		nu_stop_disk_RF  = nu_peak_BB*pt->nu_planck_max_factor;
		double (*pf) (struct spettro *, double x);
		pf = &Disk_Spectrum;

		pt->Cost_Norm_disk_Mulit_BB= 1.0/
				integrale_simp_struct(pf, pt,nu_start_disk_RF, nu_stop_disk_RF, 1000.0);
		printf( "%e\n",pt->Cost_Norm_disk_Mulit_BB);
	 }
	else if (strcmp(pt->disk_type, "Mono") == 0) {
			pt->disk = 3;
			nu_peak_BB=eval_nu_peak_Disk(pt->T_disk_max);
			nu_start_disk_RF = nu_peak_BB*pt->mono_planck_min_factor;
			nu_stop_disk_RF  = nu_peak_BB*pt->mono_planck_max_factor;
	}
	else{
		printf("wrong disk type, option BB, MultiBB, Mono \n ");
		exit(1);
	}

	pt->nu_start_Disk = eval_nu_min_blob_RF(pt,pt->Disk_mu_1, pt->Disk_mu_2, nu_start_disk_RF);
	pt->nu_stop_Disk  = eval_nu_max_blob_RF(pt,pt->Disk_mu_1, pt->Disk_mu_2, nu_stop_disk_RF);

	if (pt->verbose){
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
		pt->I_nu_Disk_disk_RF[NU_INT]=eval_I_nu_Disk_disk_RF(pt, pt->nu_Disk_disk_RF[NU_INT]);
		pt->J_nu_Disk_disk_RF[NU_INT]=eval_J_nu_Disk_disk_RF(pt, pt->I_nu_Disk_disk_RF[NU_INT]);
	}


	build_log_grid( pt->nu_start_Disk,  pt->nu_stop_Disk, pt->nu_seed_size, pt->nu_Disk);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
 		nu_obs = nu_disk_to_nu_obs_disk(pt->nu_Disk_disk_RF[NU_INT],pt->z_cosm);
		pt->nu_Disk_obs[NU_INT]=nu_obs;
 		pt->I_nu_Disk[NU_INT]=eval_I_nu_Disk_blob_RF(pt,pt->nu_Disk[NU_INT]);
 		pt->n_Disk[NU_INT] =I_nu_to_n(pt->I_nu_Disk[NU_INT], pt->nu_Disk[NU_INT]);

		if (pt->verbose>1){
			printf(" nu_Disk_disk_RF=%e, I_nu_Disk_disk_RF=%e, nu_Disk=%e, , I_nu_Disk=%e\n",
					pt->nu_Disk_disk_RF[NU_INT],
					pt->J_nu_Disk_disk_RF[NU_INT],
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

		nuL_nu_disk = eval_Disk_L_nu(pt,pt->nu_Disk_disk_RF[NU_INT]) * pt->nu_Disk_disk_RF[NU_INT];
		F_nu_disk_obs= L_nu_Disk_to_F_nu(nuL_nu_disk / pt->nu_Disk_disk_RF[NU_INT], pt-> z_cosm, pt-> dist);
		pt->nuF_nu_Disk_obs[NU_INT] = F_nu_disk_obs*nu_obs;

		fprintf(fp_SED_disk, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
			log10(nu_obs),
			log10(nu_obs * F_nu_disk_obs),
			nu_obs,
			nu_obs*F_nu_disk_obs,
			pt->nu_Disk_disk_RF[NU_INT],
			nuL_nu_disk);



	}

	fclose(fp_SED_disk);
}




void set_Disk(struct spettro *pt){
	double  nu_peak_BB;



	double M_BH,accr_rate,L_Edd,accr_Edd;

	pt->T_disk_max_4 = pt->T_disk_max * pt->T_disk_max * pt->T_disk_max * pt->T_disk_max;

	pt->R_Sw = eval_R_Sw(pt->L_disk,pt->accr_eff, pt->T_disk_max_4);

	//R_inner
	pt->R_inner = pt->R_inner_Sw * pt->R_Sw;
	//R_ext
	pt->R_ext = pt->R_ext_Sw * pt->R_Sw;
	pt->Cost_disk_Mulit_BB = 3*pt->R_Sw * pt->L_disk / (16 * pi * sigma_steph_boltz * pt->accr_eff);
	if (pt->verbose){
		printf("T_max = %e (K)\n",pt->T_disk_max);
		// energy corresponding to Tmax
		printf("E_max = %e (eV)\n",pt->T_disk_max*K_boltz*erg_to_eV);
		// frequency corresponding to Tmax
		printf("nu_max = %e (Hz)\n",pt->T_disk_max*K_boltz/HPLANCK);
		//Peak of the BB spectrum
		nu_peak_BB=eval_nu_peak_Disk(pt->T_disk_max);
		printf("nu_peak for a single BB = %e (Hz)\n",nu_peak_BB);
		printf("schwarzschild radius=%e\n", pt->R_Sw);
		printf("R_ext =%e (cm)\n", pt->R_ext);
		printf("R_inner =%e (cm)\n", pt->R_inner);
		M_BH=eval_M_BH(pt->R_Sw);
		printf("Black hole mass = %e (m_sun)\n",M_BH/m_sun );
		accr_rate=eval_accr_rate(pt->L_disk,pt->accr_eff);
		printf("Accr. rate = %e (g/s)\n",accr_rate);
		printf("Accr. rate = %e (M_sun/year)\n",accr_rate*86400.*365./m_sun);
		L_Edd=eval_L_Edd(M_BH);
		printf("L_Edd = %e (erg/s)\n",L_Edd);
		printf("L_disk = %e (erg/s)\n", pt->L_disk);
		printf("L_diks/L_edd = %e\n",pt->L_disk/L_Edd);
		accr_Edd=eval_accr_Edd(L_Edd,pt->accr_eff);


		printf("Accr_Edd = %e (g/s)\n",accr_Edd);
		printf("Accr_Edd = %e (M_sun/year)\n",accr_Edd*86400.*365./m_sun);
	}
}

//========================
// Disk Spectral Functions
//========================


double Disk_Spectrum(struct spettro *pt, double nu_Disk_disk_RF){

	if (pt->disk == 1) {
		return f_planck_norm(pt->T_disk_max, nu_Disk_disk_RF);
	}
	else if (pt->disk == 2) {
		double (*pf) (struct spettro *, double x);
		pf = &integrand_f_planck_Multi_T;
		pt->nu_disk_Multi_BB = nu_Disk_disk_RF;
		return integrale_trap_log_struct(pf, pt, pt->R_inner * 1.01, pt->R_ext, 100.0);
	}
	else if (pt->disk==3){
		return eval_nu_peak_Disk(pt->T_disk_max)*(pt->mono_planck_max_factor-pt->mono_planck_min_factor);
	}
}

double eval_I_nu_Disk_disk_RF(struct spettro *pt,double nu_Disk_disk_RF){
	double I;

	if (pt->disk == 2) {
		I = pt->L_disk *Disk_Spectrum(pt, nu_Disk_disk_RF)*pt->Disk_geom_factor*pt->Cost_Norm_disk_Mulit_BB;
	}
	else{
		I = pt->L_disk *Disk_Spectrum(pt, nu_Disk_disk_RF)*pt->Disk_geom_factor;
	}
	return I;
}

double eval_J_nu_Disk_disk_RF(struct spettro *pt, double I_nu_Disk_disk_RF){
	return I_nu_Disk_disk_RF*(2.0*pi/(4*pi))*(pt->Disk_mu_2- pt->Disk_mu_1);
}

double integrand_I_nu_Disk_blob_RF(struct spettro *pt, double mu){
	unsigned long i=0;
 	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->nu_blob_RF,pt->BulkFactor,pt->beta_Gamma,mu);
	i=x_to_grid_index( pt->nu_Disk_disk_RF,nu_disk_RF,pt->nu_seed_size);
 	if (i>0){
		return pt->J_nu_Disk_disk_RF[i]*pt->BulkFactor*(1-pt->beta_Gamma*mu);
	}
	else{
		return 0;
	}
}



double eval_I_nu_Disk_blob_RF(struct spettro *pt, double nu_blob_RF){
	unsigned long i;
	pt->nu_blob_RF=nu_blob_RF;
	double (*pf) (struct spettro *, double x);
	pf = &integrand_I_nu_Disk_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return  0.5*integrale_simp_struct(pf, pt, pt->Disk_mu_1, pt->Disk_mu_2,100.0);
}


double eval_Disk_L_nu(struct spettro *pt, double nu_Disk_disk_RF){
	if (pt->disk == 2) {
		return  pt->L_disk *Disk_Spectrum(pt, nu_Disk_disk_RF)*pt->Cost_Norm_disk_Mulit_BB;
	}
	else{
		return  pt->L_disk *Disk_Spectrum(pt, nu_Disk_disk_RF);
	}
}

double eval_nu_peak_Disk(double T){
	return eval_nu_peak_planck(T);
}



//========================
// Disk Geometrical Functions
//========================

void set_Disk_geometry(struct spettro *pt){
	double mu1,mu2;
	mu1=pt->R_H/sqrt(pt->R_H*pt->R_H+pt->R_inner*pt->R_inner);
	mu2=pt->R_H/sqrt(pt->R_H*pt->R_H+pt->R_ext*pt->R_ext);

	pt->Disk_mu_1=min(mu1,mu2);
	pt->Disk_mu_2=max(mu1,mu2);

	pt->Disk_surface=2*pi*((pt->R_ext * pt->R_ext) - (pt->R_inner*pt->R_inner) );
	pt->Disk_geom_factor=(1.0)/(pi*pt->Disk_surface);

}




//=========================================================================================
void Build_I_nu_BLR(struct spettro *pt){
	double nu_start_BLR_blob_RF,nu_stop_BLR_blob_RF;
	unsigned long NU_INT,NU_INT_MAX;
	char f_BLR_disk[512];
	FILE *fp_BLR_disk;

	sprintf(f_BLR_disk, "%s%s-I_nu_BLR.dat",pt->path, pt->STEM);

	fp_BLR_disk = fopen(f_BLR_disk, "w");
	if (fp_BLR_disk == NULL) {
		printf("unable to open %s\n ", fp_BLR_disk);
		exit(1);
	}
	//flux_DISK_header(fp_BLR_disk);
	if (pt->verbose){

		printf("-----------  Building I_nu BLR     ----------- \n");
	}
	set_BLR_geometry(pt);

	if (pt->tau_BLR>0.9){
		printf ("!!! Waring, the fraction of L_disk reaching DT is (1-tau_BLR)\n");
		printf ("!!! if tau_BLR=1.0 no DT photons will be generated\n");

	}

	pt->nu_start_BLR_disk_RF=pt->nu_Disk_disk_RF[0];
	pt->nu_stop_BLR_disk_RF=pt->nu_Disk_disk_RF[pt->NU_INT_MAX_Disk];

	pt->nu_start_BLR = eval_nu_min_blob_RF(pt,pt->BLR_mu_1, pt->BLR_mu_2, pt->nu_start_BLR_disk_RF);
	pt->nu_stop_BLR  = eval_nu_max_blob_RF(pt,pt->BLR_mu_1, pt->BLR_mu_2, pt->nu_stop_BLR_disk_RF);

	if (pt->verbose){
		printf("nu_start_BLR_disk_RF=%e  nu_stop_BLR_disk_RF=%e \n",
					pt->nu_start_BLR_disk_RF,
					pt->nu_stop_BLR_disk_RF);

		printf("nu_start_BLR=%e  nu_stop_BLR=%e \n",
					pt->nu_start_BLR,
					pt->nu_stop_BLR);
	}


	NU_INT_MAX = pt->nu_seed_size-1;
	pt->NU_INT_MAX_BLR=NU_INT_MAX;

	build_log_grid( pt->nu_start_BLR_disk_RF,  pt->nu_stop_BLR_disk_RF, pt->nu_seed_size, pt->nu_BLR_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
			pt->J_nu_BLR_disk_RF[NU_INT]=eval_J_nu_BLR_disk_RF(pt,pt->I_nu_Disk_disk_RF[NU_INT]);
	}


	build_log_grid( pt->nu_start_BLR,  pt->nu_stop_BLR, pt->nu_seed_size, pt->nu_BLR);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->I_nu_BLR[NU_INT]=eval_I_nu_BLR_blob_RF(pt,pt->nu_BLR[NU_INT]);
		pt->n_BLR[NU_INT] =I_nu_to_n(pt->I_nu_BLR[NU_INT], pt->nu_BLR[NU_INT]);
		//pt->nu_stop_BLR = pt->nu_BLR[NU_INT];
		//pt->NU_INT_MAX_BLR = NU_INT;

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
					pt->J_nu_BLR_disk_RF[NU_INT],
					pt->nu_BLR[NU_INT],
					pt->I_nu_BLR[NU_INT]);
		}
		fprintf(fp_BLR_disk, "%4.4e\t %4.4e\t %4.4e\t %4.4e \n",
			log10(pt->nu_BLR_disk_RF[NU_INT]),
			log10(pt->J_nu_BLR_disk_RF[NU_INT]),
			log10(pt->nu_BLR[NU_INT]),
			log10(pt->I_nu_BLR[NU_INT]));



	}

	fclose(fp_BLR_disk);
}






//========================
// BLR Spectral Functions
//========================
double eval_J_nu_BLR_disk_RF(struct spettro *pt,double I_nu_Disk_disk_RF){
	double (*pf) (struct spettro *, double x);
	pf = &integrand_J_nu_BLR_disk_RF;
	pt->I_nu_disk_RF=I_nu_Disk_disk_RF;
	return (2*pi/(4*pi))*integrale_simp_struct(pf, pt, pt->BLR_mu_r_J_1, pt->BLR_mu_r_J_2, 100.0);
}

double integrand_J_nu_BLR_disk_RF(struct spettro *pt,double mu){
	//printf("mu=%e, l=%e\n",mu,l_BLR(pt,mu));
	return one_by_four_pi*pt->I_nu_disk_RF*pi*pt->Disk_surface*pt->tau_BLR* l_BLR(pt,mu)/(pt->Delta_R_BLR*pt->BLR_inner_Surface);
}

double integrand_I_nu_BLR_blob_RF(struct spettro *pt, double mu){
	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->nu_blob_RF,pt->BulkFactor,pt->beta_Gamma,mu);
	unsigned long i;
	i=x_to_grid_index( pt->nu_BLR_disk_RF,nu_disk_RF,pt->nu_seed_size);

	if (i>0){

		return pt->J_nu_BLR_disk_RF[i]*pt->BulkFactor*(1-pt->beta_Gamma*mu);
	}
	else{
		return 0;
	}
}

double eval_I_nu_BLR_blob_RF(struct spettro *pt, double nu_blob_RF){
	unsigned long i;
	pt->nu_blob_RF=nu_blob_RF;
	double (*pf) (struct spettro *, double x);
	pf = &integrand_I_nu_BLR_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return  0.5*integrale_simp_struct(pf, pt, pt->BLR_mu_1, pt->BLR_mu_2,100.0);
}



//========================
// BLR Geometrical Functions
//========================
double l_BLR(struct spettro *pt, double mu){
	double S1,S2,c,l;
	l=0;
	if (pt->R_H<pt->R_BLR_in){
		S1=eval_circle_secant(pt->R_H,pt->R_BLR_in,mu);
		S2=eval_circle_secant(pt->R_H,pt->R_BLR_out,mu);
		l= (S2-S1)*0.5;
	}
	else if (pt->R_H>=pt->R_BLR_in && pt->R_H<=pt->R_BLR_out){
		if (mu>pt->BLR_mu_r_in && mu<=1.0){
			S1=eval_circle_secant(pt->R_H,pt->R_BLR_in,mu);
			S2=eval_circle_secant(pt->R_H,pt->R_BLR_out,mu);
			c=eval_circle_secant_ratio(pt->R_H,pt->R_BLR_out,mu);
			//printf("%e %e %e\n",S1,S2,S2-S1-(pt->R_BLR_out-pt->R_H)/mu);
			l= (S2-S1-(S2*(c+1/c)));
		}
		else if (mu>0 && mu<=pt->BLR_mu_r_in){
			S2=eval_circle_secant(pt->R_H,pt->R_BLR_out,mu);
			//printf("%e %e\n",S1,S2-S1-(pt->R_BLR_out-pt->R_H)/mu);
			c=eval_circle_secant_ratio(pt->R_H,pt->R_BLR_out,mu);

			l= (S2 -(S2*(c+1/c)));
		}
		else if(mu==0){
			l= sqrt(pt->R_BLR_out*pt->R_BLR_out - pt->R_H*pt->R_H)*0.5;
		}
		else{
			S2=eval_circle_secant(pt->R_H,pt->R_BLR_out,mu);
			//printf("%e %e\n",S1,S2-S1-(pt->R_BLR_out-pt->R_H)/mu);
			c=eval_circle_secant_ratio(pt->R_H,pt->R_BLR_out,mu);
			l= (S2*(c+1/c));

		}
	}
	else{
		if (mu>pt->BLR_mu_r_in && mu<=1.0){
			S1=eval_circle_secant(pt->R_H,pt->R_BLR_in,mu);
			S2=eval_circle_secant(pt->R_H,pt->R_BLR_out,mu);
			//printf("%e \n",S2-S1);
			l= (S2-S1);
		}
		else{
			l= eval_circle_secant(pt->R_H,pt->R_BLR_out,mu);
		}

	}

	return l;
}




void set_BLR_geometry(struct spettro *pt){

	pt->BLR_Volume=(4./3.)*pi*((pt->R_BLR_out*pt->R_BLR_out*pt->R_BLR_out)-(pt->R_BLR_in*pt->R_BLR_in*pt->R_BLR_in));
	pt->BLR_inner_Surface=4*pi*(pt->R_BLR_in*pt->R_BLR_in);
	pt->Delta_R_BLR=pt->R_BLR_out-pt->R_BLR_in;

	if (pt->R_H<pt->R_BLR_in){
		pt->BLR_mu_1=-1.0;
		pt->BLR_mu_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);
		pt->BLR_mu_r_J_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);
		pt->BLR_mu_r_J_1=-1.0;
		//pt->BLR_geom_factor=(1.0)/(4*pi*4*pi);
		//pt->BLR_geom_factor*=1.0/(pt->R_BLR_in*pt->R_BLR_in);

	}
	else if (pt->R_H>=pt->R_BLR_in && pt->R_H<pt->R_BLR_out){
		pt->BLR_mu_1=-1.0;
		pt->BLR_mu_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);
		pt->BLR_mu_r_in=sqrt(pt->R_H*pt->R_H - pt->R_BLR_in*pt->R_BLR_in)/pt->R_H;
		pt->BLR_mu_r_J_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);
		pt->BLR_mu_r_J_1=-1.0;
	}
	else{
		pt->BLR_mu_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);
		pt->BLR_mu_1=sqrt(pt->R_H*pt->R_H - pt->R_BLR_out*pt->R_BLR_out)/pt->R_H;
		pt->BLR_mu_r_in=sqrt(pt->R_H*pt->R_H - pt->R_BLR_in*pt->R_BLR_in)/pt->R_H;
		pt->BLR_mu_r_J_1=pt->BLR_mu_1;
		pt->BLR_mu_r_J_2=pt->BLR_mu_2;
		//pt->BLR_geom_factor=(1.0)/(4*pi*4*pi);
		//pt->BLR_geom_factor*=1.0/(pt->R_BLR_in*pt->R_BLR_in);
	}
}

//=========================================================================================

void Build_I_nu_DT(struct spettro *pt){
	FILE *fp_SED_DT;
	char f_SED_DT[512];
	unsigned long NU_INT,NU_INT_MAX;
	double nu_peak_DT_disk_RF,nu_torus_disk_RF;
	double nu_start_DT_disk_RF,nu_stop_DT_disk_RF;
	double nu_obs;
	double nuL_nu_DT,F_nu_DT_obs;

	sprintf(f_SED_DT, "%s%s-SED-DT.dat",
	            pt->path, pt->STEM);

	fp_SED_DT = fopen(f_SED_DT, "w");
	if (fp_SED_DT == NULL) {
		printf("unable to open %s\n ", f_SED_DT);
		exit(1);
	}
	flux_DISK_header(fp_SED_DT);

	if (pt->verbose){

		printf("-----------  Building I_nu DT     ----------- \n");
	}

	set_DT_geometry(pt);

	if (pt->tau_BLR>0.9){
			printf ("!!! Waring, the fraction of L_disk reaching DT is (1-tau_BLR)\n");
			printf ("!!! if tau_BLR=1.0 no DT photons will be generated\n");

	}

	nu_peak_DT_disk_RF=eval_nu_peak_planck(pt->T_DT);

	nu_start_DT_disk_RF=nu_peak_DT_disk_RF*pt->nu_planck_min_factor;
	nu_stop_DT_disk_RF=nu_peak_DT_disk_RF*pt->nu_planck_max_factor;

	pt->nu_start_DT = eval_nu_min_blob_RF(pt,pt->DT_mu_1, pt->DT_mu_2, nu_start_DT_disk_RF);
	pt->nu_stop_DT  = eval_nu_max_blob_RF(pt,pt->DT_mu_1, pt->DT_mu_2, nu_stop_DT_disk_RF);

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

	build_log_grid( nu_start_DT_disk_RF,  nu_stop_DT_disk_RF, pt->nu_seed_size, pt->nu_DT_disk_RF);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
			pt->J_nu_DT_disk_RF[NU_INT] = eval_J_nu_DT_disk_RF( pt,pt->nu_DT_disk_RF[NU_INT]);
	}


	build_log_grid( pt->nu_start_DT,  pt->nu_stop_DT, pt->nu_seed_size, pt->nu_DT);
	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {


		nu_obs = nu_disk_to_nu_obs_disk(pt->nu_DT_disk_RF[NU_INT], pt->z_cosm);

		pt->nu_DT_obs[NU_INT]=nu_obs;

		pt->I_nu_DT[NU_INT]=eval_I_nu_DT_blob_RF( pt, pt->nu_DT[NU_INT]);
		pt->n_DT[NU_INT] =I_nu_to_n(pt->I_nu_DT[NU_INT], pt->nu_DT[NU_INT]);


		if (pt->I_nu_DT[NU_INT]>pt->emiss_lim){
			pt->nu_stop_DT = pt->nu_DT[NU_INT];
			pt->NU_INT_MAX_DT = NU_INT;
		}
		else{
			pt->I_nu_DT[NU_INT]=pt->emiss_lim;
			pt->n_DT[NU_INT] =I_nu_to_n(pt->I_nu_DT[NU_INT], pt->nu_DT[NU_INT]);
		}


		if (pt->tau_BLR<1.0){
			nuL_nu_DT = eval_DT_L_nu(pt,pt->nu_DT_disk_RF[NU_INT]) * pt->nu_DT_disk_RF[NU_INT];
		}
		else {
			nuL_nu_DT = 0;
		}
		F_nu_DT_obs = L_nu_Disk_to_F_nu(nuL_nu_DT / pt->nu_DT_disk_RF[NU_INT], pt-> z_cosm, pt-> dist);

		pt->nuF_nu_DT_obs[NU_INT] = F_nu_DT_obs*nu_obs;
		if (pt->verbose>1){
			printf(" nu_DT_disk_RF=%e, I_nu_DT_disk_RF=%e, nu_DT=%e, I_nu_DT=%e\n",
				pt->nu_DT_disk_RF[NU_INT],
				pt->J_nu_DT_disk_RF[NU_INT],
				pt->nu_DT[NU_INT],
				pt->I_nu_DT[NU_INT]);
		}
		fprintf(fp_SED_DT, "%4.4e\t %4.4e\t %4.4e\t %4.4e\t%4.4e\t%4.4e \n",
								log10(nu_obs),
								log10(nu_obs * F_nu_DT_obs),
								nu_obs,
								nu_obs*F_nu_DT_obs,
								pt->nu_DT_disk_RF[NU_INT],
								nuL_nu_DT);
	}

	fclose(fp_SED_DT);


	for (NU_INT = 0; NU_INT<= NU_INT_MAX; NU_INT++) {
		pt->n_DT[NU_INT] =I_nu_to_n(pt->I_nu_DT[NU_INT], pt->nu_DT[NU_INT]);
	}
}

//========================
// Torus Spectral Functions
//========================
double eval_DT_L_nu(struct spettro *pt, double nu_torus){
	return  pt->L_disk*pt->tau_DT*f_planck_norm(pt->T_DT, nu_torus);
}

double eval_J_nu_DT_disk_RF(struct spettro *pt,double nu_torus_disk_RF){
	double (*pf) (struct spettro *, double x);
		pf = &integrand_J_nu_DT_disk_RF;
		pt->nu_disk_RF=nu_torus_disk_RF;
		return (2*pi/(4*pi))*integrale_simp_struct(pf, pt, pt->DT_mu_r_J_1, pt->DT_mu_r_J_2, 100.0);
}




double integrand_J_nu_DT_disk_RF(struct spettro *pt,double mu){
	if (pt->R_H<pt->R_DT){
		return one_by_four_pi*pt->L_disk * pt->tau_DT*(1-pt->tau_BLR)*
				f_planck_norm(pt->T_DT, pt->nu_disk_RF)/
				pt->DT_Volume* l_TORUS(pt,  mu);
	}
	else {
		return one_by_four_pi*pt->L_disk * pt->tau_DT*
				f_planck_norm(pt->T_DT, pt->nu_disk_RF)/
						pt->DT_Volume* l_TORUS(pt,  mu);
	}

	return 0;
}




double integrand_I_nu_DT_blob_RF(struct spettro *pt, double mu){
	double nu_disk_RF=nu_blob_RF_to_nu_disk_RF(pt->nu_blob_RF,pt->BulkFactor,pt->beta_Gamma,mu);
	unsigned long i;
	i=x_to_grid_index( pt->nu_DT_disk_RF,nu_disk_RF,pt->nu_seed_size);
	if (i>0){
		return pt->J_nu_DT_disk_RF[i]*pt->BulkFactor*(1-pt->beta_Gamma*mu);
	}
	else {
		return 0;
	}
}

double eval_I_nu_DT_blob_RF(struct spettro *pt, double nu_blob_RF){
	unsigned long i;
	pt->nu_blob_RF=nu_blob_RF;
	double (*pf) (struct spettro *, double x);
	pf = &integrand_I_nu_DT_blob_RF;
	//0.5 comes from 2pi/(4pi)
	return  0.5*integrale_simp_struct(pf, pt, pt->DT_mu_1, pt->DT_mu_2,100.0);
}




//========================
// Torus Geometrical Functions
//========================
double l_TORUS(struct spettro *pt, double mu){
	double l;
	l=0;
	if (pt->R_H<pt->R_DT){
		l=eval_circle_secant(pt->R_H,pt->R_DT,mu)*0.5;
	}
	else{
		l=eval_circle_secant(pt->R_H,pt->R_DT,mu);
	}

	return l;
}


void set_DT_geometry(struct spettro * pt){
	if (pt->R_H<pt->R_DT){
		pt->DT_mu_1=-1.0;
		pt->DT_mu_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);

		pt->DT_Volume=(4./3.)*pi*pt->R_DT*pt->R_DT*pt->R_DT;
		pt->DT_mu_r_J_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);
		pt->DT_mu_r_J_1=-1.0;
	}
	else{
		pt->DT_mu_1=sqrt(pt->R_H*pt->R_H - pt->R_DT*pt->R_DT)/pt->R_H;
		pt->DT_mu_2=1.0;

		pt->DT_mu_r_J_2=pt->R_H/sqrt(pt->R_H*pt->R_H + pt->R_ext*pt->R_ext);
		pt->DT_mu_r_J_1=sqrt(pt->R_H*pt->R_H - pt->R_DT*pt->R_DT)/pt->R_H;

		pt->DT_Volume=(4./3.)*pi*pt->R_DT*pt->R_DT*pt->R_DT;
	}
}


//=========================================================================================








//========================
// Accretion Power Physical Functions
//========================
double eval_R_Sw(double L_disk, double accr_eff, double T_disk_max_4){
	double a;
	a= pow(0.140836,4) * L_disk / (pi * sigma_steph_boltz * accr_eff * T_disk_max_4);
	return pow(a, 0.5);
}

double eval_M_BH(double R_Sw){
	return R_Sw * vluce_cm * vluce_cm / (2 * G_cgs);
}


double eval_accr_rate(double L_disk,double accr_eff){
 return L_disk/(vluce_cm * vluce_cm*accr_eff);
}


double eval_L_Edd(double M_BH){
	return 4*pi*G_cgs*M_BH*vluce_cm/0.3;
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


double T_disk(struct spettro *pt, double R){
	double a, T_disco_r;
	T_disco_r = pt->Cost_disk_Mulit_BB/(R*R*R) * (1 - pow((pt->R_inner / R), 0.5));
	T_disco_r = pow(T_disco_r, 0.25);
	return T_disco_r;
}

double f_planck_Multi_T(struct spettro *pt, double R ,double nu) {
    return f_planck(T_disk(pt,R) , nu);
}

double f_planck_Multi_T_norm(struct spettro *pt, double R, double nu) {
	double T;
	T=T_disk(pt,R);
    return f_planck_norm(T,nu);
}

double integrand_f_planck_Multi_T(struct spettro *pt, double R){
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
double eval_nu_min_blob_RF(struct spettro *pt, double mu1, double mu2, double nu_disk_RF ){
	double a,nu_1,nu_2;
	nu_1=nu_disk_RF*pt->BulkFactor*(1-pt->beta_Gamma*mu1);
	nu_2=nu_disk_RF*pt->BulkFactor*(1-pt->beta_Gamma*mu2);
	a=  min(nu_1,nu_2);
	return a;
}



double eval_nu_max_blob_RF(struct spettro *pt, double mu1, double mu2, double nu_disk_RF ){
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

double eval_circle_secant_ratio(double z,double R,double mu){
	double x2,x1,l,m,b,c,a;
	m=tan(acos(mu));
	b=-2*z;
	c=z*z-R*R;
	a=1+m*m;
	if ((b*b-4*a*c)>0){
		x1=(-b +sqrt(b*b-4*a*c))/(2*a);
		x2=(-b -sqrt(b*b-4*a*c))/(2*a);
		//printf("%e %e\n",x1,x2);
		return fabs(x1)/fabs(x2);
	}
	else return 0;
}


double eval_circle_secant(double z,double R,double mu){
	double x1,x2,y1,y2,l,m,b,c,a;
	m=tan(acos(mu));
	b=-2*z;
	c=z*z-R*R;
	a=1+m*m;
	if ((b*b-4*a*c)>0){
		x1=(-b +sqrt(b*b-4*a*c))/(2*a);
		x2=(-b -sqrt(b*b-4*a*c))/(2*a);
		//printf("%e %e\n",x1,x2);
		y1=m*x1;
		y2=m*x2;
		return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}
	else return 0;
}


//=========================================================================================






