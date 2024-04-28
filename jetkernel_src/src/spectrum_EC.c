//=========================================================================================
//                   CALCOLO DELLO SPETTRO EC
//=========================================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file spettro_Compton.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief CALCOLO DELLO SPETTRO EC
 *
 */


double f_psi_EC_ring(double R_ext,double R_H, double mu_s,double beaming,double phi){
 	double x2, cos_psi, mu_star;
 	x2 =  (R_ext*R_ext)+(R_H*R_H);
 	mu_star = R_H/sqrt(x2);
	cos_psi=(mu_s*mu_star)+(sqrt((1-(mu_s*mu_s)))*sqrt(1-(mu_star*mu_star)))*cos(phi);
 	return ((1 - cos_psi) * (1 - cos_psi)) * pow(beaming,6) / x2;
	}


double f_psi_EC_sphere(double R_ext,double R_H, double mu_s, double mu_re, double beaming,double phi){
	double x2, cos_psi, mu_star;
	x2= (R_ext*R_ext)+(R_H*R_H) -2*R_H*R_ext*mu_re;
	mu_star=sqrt(1- ((R_ext*R_ext/x2) *(1-mu_re*mu_re)) );
	//printf("x2=%e mu_re=%e mu_star=%e\n", x2,mu_re,mu_star);
	if (mu_star>1){
		return 0.;
	}
	cos_psi=(mu_s*mu_star) +( sqrt(1-(mu_star*mu_star))*sqrt(1-(mu_s*mu_s))*cos(phi));
	return ((1 - cos_psi) * (1 - cos_psi)) * pow(beaming,6)/x2;
}

double beaming_pattern_EC(double theta_s, double R_ext, double R_H, double Gamma, int geom){
	// geom = 0 sphere
	// geom = 1 ring
	
	unsigned  int_size = 100;
	double phi[100], mu_re[100], y[100], z[100];
	unsigned int i,j;
	double delta_phi, beaming, mu_s,mu_re_max,mu_re_min,delta_mu_re,bp;

	beaming = get_beaming(Gamma, theta_s);

	mu_re_min=-1;
	if (R_H > R_ext){
         mu_re_max= R_ext/ R_H;
	
    }else{
		mu_re_max = 1;
	}	

	delta_phi = 2 * pi / (int_size-1);
	delta_mu_re  = (mu_re_max-mu_re_min)/ (int_size-1);
	
	mu_s = cos(Deg_to_Rad * theta_s);
	for (i = 0; i < int_size; i++)
	{
		phi[i] = 0 + (i * delta_phi);
		mu_re[i] = mu_re_min + (i*delta_mu_re);

	}
	
	if (geom==0){
		for (i = 0; i < int_size; i++){
			for (j = 0; j < int_size; j++){
			
		 	z[j]= f_psi_EC_sphere(R_ext, R_H, mu_s, mu_re[j],beaming, phi[i]);
	
			}	
			y[i]= trapzd_array_linear_grid(mu_re, z, int_size);
		}
		bp= trapzd_array_linear_grid(phi, y, int_size);
	
	}else if(geom==1){
		for (i = 0; i < int_size; i++){
			y[i]=  f_psi_EC_ring(R_ext, R_H, mu_s, beaming, phi[i]);
			}
		bp=trapzd_array_linear_grid(phi, y, int_size);
	}else{
		printf("wrong geometry for beaming pattern \n ");
        exit(0);	
	}
	return bp;
}

double scaling_function_EC(double theta_s, double R_ext, double R_H_in, double R_H_orig, double Gamma){
	double y_theta, y_theta_0;

	y_theta = beaming_pattern_EC(theta_s, R_ext, R_H_orig, Gamma,0);
	y_theta_0 = beaming_pattern_EC(theta_s, R_ext, R_H_in, Gamma,0);
	return y_theta / y_theta_0;
}

void update_EC_for_bp(struct blob *pt, double nuFnu_obs_ref, double R_ext_emit, unsigned int SIZE, double *nuFnu_obs, double *nu_obs)
{
	double s_bp, s_actual, nuFnu_obs_max;
	unsigned int I_MAX, NU_INT;

	s_bp = scaling_function_EC(pt->theta, R_ext_emit, 0, pt->R_H_orig, pt->BulkFactor);

	I_MAX = pt->nu_IC_size - 1;	
	nuFnu_obs_max = nuFnu_obs[0];
	for (NU_INT = 0; NU_INT < I_MAX; NU_INT++)
	{
		if (nuFnu_obs[NU_INT] > nuFnu_obs_max)
		{
			nuFnu_obs_max = nuFnu_obs[NU_INT];
		}
	}
	s_actual = nuFnu_obs_max / nuFnu_obs_ref;
	for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++)
	{
		if (nuFnu_obs[NU_INT] > pt->emiss_lim)
		{
			nuFnu_obs[NU_INT] = nuFnu_obs[NU_INT] * (s_bp / s_actual);
			nu_obs[NU_INT] = nu_obs[NU_INT] * pow((s_bp /(pt->BulkFactor* s_actual)),0.25);
			if (nuFnu_obs[NU_INT] <= pt->emiss_lim)
			{
					nuFnu_obs[NU_INT] = pt->emiss_lim;
			}
			
		}
	}

}

double get_EC_reference(struct blob *pt, double *nuFnu_obs)
{
	double nuFnu_obs_ref;
	unsigned int I_MAX, NU_INT;
	

	I_MAX = pt->nu_IC_size - 1;

	nuFnu_obs_ref = nuFnu_obs[0];
	for (NU_INT = 0; NU_INT < I_MAX; NU_INT++)
	{
		if (nuFnu_obs[NU_INT] > nuFnu_obs_ref)
		{
			nuFnu_obs_ref = nuFnu_obs[NU_INT];
		}
	}
	return nuFnu_obs_ref;
}

int set_condition_EC_correction(struct blob *pt,double R_ext_emit)
{
	int do_EC_correction =0;
	if ((pt->R_H > (R_ext_emit * pt->R_ext_emit_factor)) && (pt->EC_stat == 1) && R_ext_emit > 0){
		do_EC_correction =1;
	}
	//printf("do_EC_correction=%d \n",do_EC_correction);
	return do_EC_correction;
}

void set_EC_stat_pre(struct blob *pt, double R_ext_emit)
{
	
	//printf("set_EC_stat_pre 1 R_ext_emit =%e, R_H_orig=%e, R_H=%e\n", R_ext_emit, pt->R_H_orig, pt->R_H);

	if (set_condition_EC_correction(pt, R_ext_emit) > 0 && R_ext_emit > 0 && pt->EC_stat==1)
	{
		pt->R_H_scale_factor = pt->BulkFactor / get_beaming( pt->BulkFactor,pt->theta);
		if ((pt->R_H / (R_ext_emit * pt->R_ext_emit_factor)) > pt->R_H_scale_factor)
		{
			pt->EC_stat = 0;
		}
		//pt->R_H = R_ext_emit * pt->R_H_scale_factor;
	}
}

void set_EC_stat_post(struct blob *pt)
{
	pt->EC_stat = pt->EC_stat_orig;
	pt->R_H = pt->R_H_orig;
	pt->R_H_scale_factor=1.0;
}


void spettro_EC(int Num_file, struct blob *pt) {
    double L_nu_EC, F_nu_EC_obs,nu_peak;
    double  gmax,numax_KN,numax_TH;
    double j_nu_disk;
    double * freq_array, *freq_array_obs;
    double * nuFnu_obs_array;
    double * nu_start_EC, * nu_stop_EC, * nu_start_EC_obs, * nu_stop_EC_obs, nu_seed_max;
    unsigned int * NU_INT_STOP_EC;
    unsigned int NU_INT, I_MAX, stop;
	double R_ext_emit;
	void *(*eval_j_ptr)(void * args);
	//============================================================
    //         inizio  loop sulle freq per spettro  compton
    //============================================================
    pt->TOT = 0;

    pt->ord_comp = 1;
    pt->SSC = 0;
    stop = 0;

    gmax=Find_gmax(pt,pt->Ne,pt->griglia_gamma_Ne_log);



    if (pt->verbose>0) {
    	printf("**********************  EC spectrum        *******************************\n");
    	printf("gmax from Ne>0 = %e\n", gmax);
        printf("-------------------------------------------------------------------\n");

    }

	//TODO check nu_seed_max for EC_stat=1
	if (pt->EC == 1) {

    	freq_array_obs=pt->nu_EC_Disk_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_Disk_obs;
    	freq_array=pt->nu_EC_Disk;
    	nu_seed_max =  pt->nu_stop_Disk;
    	nu_start_EC = &(pt->nu_start_EC_Disk);
    	nu_stop_EC = &(pt->nu_stop_EC_Disk);
		nu_start_EC_obs = &(pt->nu_start_EC_Disk_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_Disk_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_Disk);
		R_ext_emit= pt->R_ext;
		if (pt->verbose>0) {
    		printf("nu_star_Disk=%e    nu_stop_Disk=%e\n",
    			pt->nu_start_Disk,
    			pt->nu_stop_Disk);
    		printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
    		printf("-----------------------------------------------------------------\n");
    	}
    }else if (pt->EC == 2) {
    	freq_array_obs=pt->nu_EC_BLR_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_BLR_obs;
    	freq_array=pt->nu_EC_BLR;
    	nu_seed_max =  pt->nu_stop_BLR;
    	nu_start_EC = &(pt->nu_start_EC_BLR);
    	nu_stop_EC = &(pt->nu_stop_EC_BLR);
    	nu_start_EC_obs = &(pt->nu_start_EC_BLR_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_BLR_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_BLR);
		R_ext_emit = pt->R_BLR_out;
		if (pt->verbose>0) {
			printf("nu_star_BLR=%e    nu_stop_BLR=%e\n",
					pt->nu_start_BLR,
					pt->nu_stop_BLR);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}
    }
    else if (pt->EC == 3) {
    	freq_array_obs=pt->nu_EC_DT_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_DT_obs;
    	freq_array=pt->nu_EC_DT;
    	nu_seed_max =  pt->nu_stop_DT;
    	nu_start_EC = &(pt->nu_start_EC_DT);
    	nu_stop_EC = &(pt->nu_stop_EC_DT);
    	nu_start_EC_obs = &(pt->nu_start_EC_DT_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_DT_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_DT);
		R_ext_emit = pt->R_DT;
		//printf("pre R_H=%e R_DT=%e EC_Stat=%d\n",pt->R_H,pt->R_DT,pt->EC_stat);
    	if (pt->verbose>0) {
			printf("nu_star_DT=%e    nu_stop_DT=%e\n",
					pt->nu_start_DT,
					pt->nu_stop_DT);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}

    }
    else if (pt->EC == 4) {
    	freq_array_obs=pt->nu_EC_Star_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_Star_obs;
    	freq_array=pt->nu_EC_Star;
    	nu_seed_max =  pt->nu_stop_Star;
    	nu_start_EC = &(pt->nu_start_EC_Star);
    	nu_stop_EC = &(pt->nu_stop_EC_Star);
    	nu_start_EC_obs = &(pt->nu_start_EC_Star_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_Star_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_Star);
		R_ext_emit = -1;
		if (pt->verbose>0) {
			printf("nu_start_Star=%e    nu_stop_Star=%e\n",
					pt->nu_start_Star,
					pt->nu_stop_Star);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}

    }
    else if (pt->EC == 5) {
    	freq_array_obs=pt->nu_EC_CMB_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_CMB_obs;
    	freq_array=pt->nu_EC_CMB;
    	nu_seed_max =  pt->nu_stop_CMB;
    	nu_start_EC = &(pt->nu_start_EC_CMB);
    	nu_stop_EC = &(pt->nu_stop_EC_CMB);
    	nu_start_EC_obs = &(pt->nu_start_EC_CMB_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_CMB_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_CMB);
		R_ext_emit = -1;
		
		if (pt->verbose>0) {
    		printf("nu_start_CMB=%e    nu_stop_CMB=%e\n",
    				pt->nu_start_CMB,
    				pt->nu_stop_CMB);
    		printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
    		printf("-----------------------------------------------------------------\n");
    	}

    }else{
		printf("wrong EC \n ");
        exit(0);
		
	}
	//printf("spettro_EC 1 R_H=%e c=%d \n", pt->R_H, set_condition_EC_correction(pt, pt->R_DT));
	set_EC_stat_pre(pt, R_ext_emit);
	//printf("spettro_EC 2 R_H=%e c=%d \n", pt->R_H, set_condition_EC_correction(pt, pt->R_DT));

	gmax=Find_gmax(pt,pt->Ne,pt->griglia_gamma_Ne_log);
	numax_KN = 1000 * gmax * MEC2 / HPLANCK;
	numax_TH = 1000 * (4.0 / 3.0) * pow(gmax, 2) * nu_seed_max;
	if (HPLANCK * nu_seed_max * gmax / MEC2 > 0.1) {
		*nu_stop_EC = numax_KN;
	} else {
		*nu_stop_EC = numax_TH;
	}


	if (HPLANCK * nu_seed_max * pt->Gamma_p3 / MEC2 > 0.1) {
		nu_peak= pt->Gamma_p3 * nu_seed_max;
	} else {
		nu_peak= (4.0 / 3.0) * pow(pt->Gamma_p3, 2) * nu_seed_max;
	}


   	//*nu_stop_EC = 100. * (4.0 / 3.0) * pow(gmax, 2) * nu_seed_max;

	if (pt->do_Sync>0){
		*nu_start_EC =  pt->nu_peak_Sync_blob;
	}else{
		*nu_start_EC =  nu_seed_max/10;
	}

   	*nu_stop_EC_obs = nu_blob_to_nu_obs(*nu_stop_EC, pt->beam_obj, pt->z_cosm);
   	*nu_start_EC_obs=nu_blob_to_nu_obs(*nu_start_EC, pt->beam_obj, pt->z_cosm);


   	build_log_grid(*nu_start_EC,  *nu_stop_EC, pt->nu_IC_size, freq_array);
   	build_log_grid(*nu_start_EC_obs,  *nu_stop_EC_obs, pt->nu_IC_size, freq_array_obs);

	


   	if (pt->verbose>0) {
		printf("nu_start_EC=%e nu_stop_EC=%e nu_peak(estim.)=%e\n",
				*nu_start_EC,
				*nu_stop_EC,
				nu_peak);
		printf("nu_start_SSC=%e nu_stop_SSC=%e\n",
				pt->nu_start_SSC,
				pt->nu_stop_SSC);
		printf("SSC=%d EC=%d TOT=%d\n",
				pt->SSC, pt->EC,
				pt->TOT);
		printf("emiss limit=%e\n", pt->emiss_lim);
	}

	I_MAX = pt->nu_IC_size -1;
	stop=0;

	eval_j_ptr = &eval_j_EC;
    threaded_j_evaluation(pt, eval_j_ptr, pt->j_EC, freq_array, *nu_start_EC, *nu_stop_EC,I_MAX,pt->N_THREADS);

	for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++) {
        nuFnu_obs_array[NU_INT]=pt->emiss_lim;
        if (pt->verbose>1) {
            printf("#-> nu_em=%e  nu_obs=%e  i=%d\n", freq_array[NU_INT], freq_array_obs[NU_INT], NU_INT);
        }
        if ((freq_array[NU_INT] >= *nu_start_EC) && (freq_array[NU_INT] <= *nu_stop_EC)) {
			if (!stop) {
				
				if (pt->verbose > 1) {
					printf("#-> q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,
							pt->q_comp[NU_INT], NU_INT, pt->j_EC[NU_INT],
							freq_array[NU_INT]);
				}

				//nu_src = nu_blob_to_nu_src(freq_array[NU_INT], pt->beam_obj,
				//		pt->z_cosm);
				if (pt->EC_stat == 1)
				{
					L_nu_EC = j_nu_src_to_L_nu_src(pt->j_EC[NU_INT], pt->Vol_region,
												   pt->beam_obj);
				}
				else{
					L_nu_EC = j_nu_to_L_nu_src(pt->j_EC[NU_INT], pt->Vol_region,
						pt->beam_obj);
				}
				//nuL_nu_EC = L_nu_EC * nu_src;
				F_nu_EC_obs = L_nu_src_to_F_nu(L_nu_EC, pt->beam_obj,
						pt->z_cosm, pt->dist);
                
				nuFnu_obs_array[NU_INT] = F_nu_EC_obs * freq_array_obs[NU_INT];
                

				if (pt->j_EC[NU_INT] < pt->emiss_lim) {
					//out=0;
					if (freq_array[NU_INT] > numax_TH) {
						stop = 1;
					}
					nuFnu_obs_array[NU_INT] = pt->emiss_lim;
					pt->j_EC[NU_INT] = pt->emiss_lim;
					pt->q_comp[NU_INT] = pt->emiss_lim;
				}
				//else{
				//	out=1;
				//}

				if (stop == 1 && freq_array[NU_INT] > numax_TH) {

					*nu_stop_EC_obs = freq_array_obs[NU_INT];
					*nu_stop_EC = freq_array[NU_INT];
					*NU_INT_STOP_EC = NU_INT;
					if (pt->verbose > 1) {
						printf("%e %d\n ", freq_array[NU_INT], NU_INT);
					}
				}
			}
            else{
				nuFnu_obs_array[NU_INT]=pt->emiss_lim;
				pt->j_EC[NU_INT] = pt->emiss_lim;
				pt->q_comp[NU_INT]=pt->emiss_lim;

			}
            //===========================================
            // FILES output nu dnu nuFnu dnuFnu
            //===========================================
            /*
			if (pt->WRITE_TO_FILE==1){
				if (!stop && out) {
					fprintf(fp_EC, "%4.4e\t%4.4e\t%4.4e\t %4.4e\t%4.4e\t%4.4e\n",
							log10(freq_array_obs[NU_INT]),
							log10(F_nu_EC_obs * freq_array_obs[NU_INT]),
							freq_array_obs[NU_INT],
							F_nu_EC_obs*freq_array_obs[NU_INT],
							nu_src,
							nuL_nu_EC);
				}
			}
			*/
            if (pt->verbose>1) {
                printf("#-> ********************************\n\n");
            }
            //==========================  END of Loop ove frequencies ====================================

        }
    }

    //Se ancora non ha trovato nu_stop
    if (!stop) {
    	*nu_stop_EC_obs = freq_array_obs[NU_INT-1];
        *nu_stop_EC = freq_array[NU_INT-1];
        *NU_INT_STOP_EC = NU_INT-1;
        if (pt->verbose > 1) {
            printf("%e %d\n ", freq_array[NU_INT-1], NU_INT-1);
        }
    	if (pt->verbose>0) {
			if (pt->EC == 1) {

				printf("nu_stop_EC_Disk=%e NU_INT_STOP_EC_Disk=%d\n", pt->nu_stop_EC_Disk, pt->NU_INT_STOP_EC_Disk);
			}
			if (pt->EC == 2) {

				printf("nu_stop_EC_BLR=%e NU_INT_STOP_EC_BLR=%d\n", pt->nu_stop_EC_BLR, pt->NU_INT_STOP_EC_BLR);
			}
			if (pt->EC == 3) {

				printf("nu_stop_EC_DT=%e NU_INT_STOP_EC_DT=%d\n", pt->nu_stop_EC_DT, pt->NU_INT_STOP_EC_DT);
			}
			if (pt->EC == 4) {

				printf("nu_stop_EC_Star=%e NU_INT_STOP_EC_Star=%d\n", pt->nu_stop_EC_Star, pt->NU_INT_STOP_EC_Star);
			}
			if (pt->EC == 5) {

				printf("nu_stop_EC_CMB=%e NU_INT_STOP_EC_CMB=%d\n", pt->nu_stop_EC_CMB, pt->NU_INT_STOP_EC_CMB);
			}
            //if (pt->EC == 6) {

            //    printf("nu_stop_EC_CMB_stat=%e NU_INT_STOP_EC_CMB_stat=%d\n", pt->nu_stop_EC_CMB_stat, pt->NU_INT_STOP_EC_CMB_stat);
            //}
    	}
    }
	//printf("spettro_EC 3 R_H=%e c=%d \n", pt->R_H, set_condition_EC_correction(pt, pt->R_DT));

	set_EC_stat_post(pt);
	//printf("spettro_EC 4 R_H=%e c=%d \n", pt->R_H, set_condition_EC_correction(pt, pt->R_DT));

	//spectra_External_Fields(1, pt,1);
	//printf("spettro_EC 5 R_H=%e c=%d \n", pt->R_H, set_condition_EC_correction(pt, pt->R_DT));

	//===========================================
	//    trova nu peak e Flux peak
	//===========================================
	if (pt->EC == 1)
	{
		FindEpSp(freq_array, nuFnu_obs_array, pt->NU_INT_STOP_EC_Disk, pt,
					&(pt->nu_peak_EC_Disk_obs),
					&(pt->nu_peak_EC_Disk_src),
					&(pt->nu_peak_EC_Disk_blob),
					&(pt->nuFnu_peak_EC_Disk_obs),
					&(pt->nuLnu_peak_EC_Disk_src),
					&(pt->nuLnu_peak_EC_Disk_blob));

		if (pt->verbose > 0)
		{
			printf("nu_stop_EC_Disk=%e NU_INT_STOP_EC_Disk=%d\n", pt->nu_stop_EC_Disk, pt->NU_INT_STOP_EC_Disk);
			printf("EC Disk ");
			printf("nu_EC_blob peak=%e\n", pt->nu_peak_EC_Disk_blob);
			printf("nu_EC_src  peak=%e\n", pt->nu_peak_EC_Disk_src);
			printf("nu_EC_obs  peak=%e\n", pt->nu_peak_EC_Disk_obs);

			printf("nuFnu EC  blob    peak=%e\n", pt->nuFnu_peak_EC_Disk_obs);
			printf("nuLnu EC  src     peak=%e\n", pt->nuLnu_peak_EC_Disk_src);
			printf("nuLnu EC  obs     peak=%e\n", pt->nuLnu_peak_EC_Disk_blob);
		}
	}

	if (pt->EC == 2)
	{
		FindEpSp(pt->nu_EC_BLR, nuFnu_obs_array, pt->NU_INT_STOP_EC_BLR, pt,
					&(pt->nu_peak_EC_BLR_obs),
					&(pt->nu_peak_EC_BLR_src),
					&(pt->nu_peak_EC_BLR_blob),
					&(pt->nuFnu_peak_EC_BLR_obs),
					&(pt->nuLnu_peak_EC_BLR_src),
					&(pt->nuLnu_peak_EC_BLR_blob));
		if (pt->verbose > 0)
		{
			printf("nu_stop_EC_BLR=%e NU_INT_STOP_EC_BLR=%d\n", pt->nu_stop_EC_BLR, pt->NU_INT_STOP_EC_BLR);
			printf("EC BLR ");
			printf("nu_EC_blob peak=%e\n", pt->nu_peak_EC_BLR_blob);
			printf("nu_EC_src  peak=%e\n", pt->nu_peak_EC_BLR_src);
			printf("nu_EC_obs  peak=%e\n", pt->nu_peak_EC_BLR_obs);

			printf("nuFnu EC  blob    peak=%e\n", pt->nuFnu_peak_EC_BLR_obs);
			printf("nuLnu EC  src     peak=%e\n", pt->nuLnu_peak_EC_BLR_src);
			printf("nuLnu EC  obs     peak=%e\n", pt->nuLnu_peak_EC_BLR_blob);
		}
	}

	if (pt->EC == 3)
	{
		FindEpSp(pt->nu_EC_DT, nuFnu_obs_array, pt->NU_INT_STOP_EC_DT, pt,
					&(pt->nu_peak_EC_DT_obs),
					&(pt->nu_peak_EC_DT_src),
					&(pt->nu_peak_EC_DT_blob),
					&(pt->nuFnu_peak_EC_DT_obs),
					&(pt->nuLnu_peak_EC_DT_src),
					&(pt->nuLnu_peak_EC_DT_blob));
		if (pt->verbose > 0)
		{
			printf("nu_stop_EC_DT=%e NU_INT_STOP_EC_DT=%d\n", pt->nu_stop_EC_DT, pt->NU_INT_STOP_EC_DT);
			printf("EC DT ");
			printf("nu_EC_blob peak=%e\n", pt->nu_peak_EC_DT_blob);
			printf("nu_EC_src  peak=%e\n", pt->nu_peak_EC_DT_src);
			printf("nu_EC_obs  peak=%e\n", pt->nu_peak_EC_DT_obs);

			printf("nuFnu EC  blob    peak=%e\n", pt->nuFnu_peak_EC_DT_obs);
			printf("nuLnu EC  src     peak=%e\n", pt->nuLnu_peak_EC_DT_src);
			printf("nuLnu EC  obs     peak=%e\n", pt->nuLnu_peak_EC_DT_blob);
		}
	}

	//printf("=>done\n");
	return;
}
//=========================================================================================




void  * eval_j_EC(void *data){
    unsigned int NU_INT;
    double nu_IC_out;
	struct j_args *thread_args = data;
    for (NU_INT = thread_args->NU_INT_START; NU_INT <= thread_args->NU_INT_STOP; NU_INT++) {
        nu_IC_out=thread_args->nu_array[NU_INT];
        thread_args->blob_pt->q_comp[NU_INT] = 0.;
        thread_args->blob_pt->j_EC[NU_INT] = 0.;
       
        if (thread_args->blob_pt->verbose > 1) {
                printf("#->1 in eval_j_EC   NU_INT=%d eval_j_EC  nu_1=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }

		thread_args->blob_pt->q_comp[NU_INT] = rate_compton_GR(thread_args->blob_pt,nu_IC_out);
		if (thread_args->blob_pt->EC_stat == 1){
			//in this case we have q_comp in the disk frame, so j_nu is in the disk rest frame
			//and we have to use also the scattered nu in the disk rest frame
			thread_args->blob_pt->j_EC[NU_INT]=thread_args->blob_pt->q_comp[NU_INT]*HPLANCK*thread_args->nu_array[NU_INT]*thread_args->blob_pt->beam_obj;
		}
		else{
			thread_args->blob_pt->j_EC[NU_INT] = thread_args->blob_pt->q_comp[NU_INT] * HPLANCK * thread_args->nu_array[NU_INT];
		}
        if (thread_args->blob_pt->verbose > 1) {
                printf("#->2 in  eval_j_EC NU_INT=%d q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,NU_INT,
                        thread_args->blob_pt->q_comp[NU_INT], NU_INT, thread_args->blob_pt->j_EC[NU_INT],
                        thread_args->nu_array[NU_INT]);
        }
    }
	return NULL; 
}


