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

static double eval_R_H_EC(struct blob *pt, double R_H)
{
	if (pt->core.EC == 6){
		return fabs(R_H - pt->Corona.R_H_Corona);
	}
	return R_H;
}

void update_EC_for_bp(struct blob *pt, double nuFnu_obs_ref, double R_ext_emit, unsigned int SIZE, double *nuFnu_obs, double *nu_obs)
{
	double s_bp, s_actual, nuFnu_obs_max, R_H_orig_eval;
	unsigned int I_MAX, NU_INT;
	(void)SIZE;

	R_H_orig_eval = eval_R_H_EC(pt, pt->core.R_H_orig);
	s_bp = scaling_function_EC(pt->core.theta, R_ext_emit, 0, R_H_orig_eval, pt->core.BulkFactor);

	I_MAX = pt->core.nu_IC_size - 1;	
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
		if (nuFnu_obs[NU_INT] > pt->core.emiss_lim)
		{
			nuFnu_obs[NU_INT] = nuFnu_obs[NU_INT] * (s_bp / s_actual);
			nu_obs[NU_INT] = nu_obs[NU_INT] * pow((s_bp /(pt->core.BulkFactor* s_actual)),0.25);
			if (nuFnu_obs[NU_INT] <= pt->core.emiss_lim)
			{
					nuFnu_obs[NU_INT] = pt->core.emiss_lim;
			}
			
		}
	}

}

double get_EC_reference(struct blob *pt, double *nuFnu_obs)
{
	double nuFnu_obs_ref;
	unsigned int I_MAX, NU_INT;
	

	I_MAX = pt->core.nu_IC_size - 1;

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
	double R_H_eval;
	int do_EC_correction =0;
	R_H_eval = eval_R_H_EC(pt, pt->core.R_H);
	if ((R_H_eval > (R_ext_emit * pt->core.R_ext_emit_factor)) && (pt->core.EC_stat == 1) && R_ext_emit > 0){
		do_EC_correction =1;
	}
	//printf("do_EC_correction=%d \n",do_EC_correction);
	return do_EC_correction;
}

void set_EC_stat_pre(struct blob *pt, double R_ext_emit)
{
	double R_H_eval;
	
	//printf("set_EC_stat_pre 1 R_ext_emit =%e, R_H_orig=%e, R_H=%e\n", R_ext_emit, pt->core.R_H_orig, pt->core.R_H);
	R_H_eval = eval_R_H_EC(pt, pt->core.R_H);

	if (set_condition_EC_correction(pt, R_ext_emit) > 0 && R_ext_emit > 0 && pt->core.EC_stat==1)
	{
		pt->core.R_H_scale_factor = pt->core.BulkFactor / get_beaming( pt->core.BulkFactor,pt->core.theta);
		if ((R_H_eval / (R_ext_emit * pt->core.R_ext_emit_factor)) > pt->core.R_H_scale_factor)
		{
			pt->core.EC_stat = 0;
		}
		//pt->core.R_H = R_ext_emit * pt->core.R_H_scale_factor;
	}
}

void set_EC_stat_post(struct blob *pt)
{
	pt->core.EC_stat = pt->core.EC_stat_orig;
	pt->core.R_H = pt->core.R_H_orig;
	pt->core.R_H_scale_factor=1.0;
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
	struct ec_comp *ec = NULL;
	//============================================================
    //         inizio  loop sulle freq per spettro  compton
    //============================================================
    pt->core.TOT = 0;

    pt->core.ord_comp = 1;
    pt->core.SSC = 0;
    stop = 0;

    gmax=Find_gmax(pt,pt->emitters.Ne,pt->emitters.griglia_gamma_Ne_log);



    if (pt->core.verbose>0) {
    	printf("**********************  EC spectrum        *******************************\n");
    	printf("gmax from Ne>0 = %e\n", gmax);
        printf("-------------------------------------------------------------------\n");

    }

	//TODO check nu_seed_max for EC_stat=1
	if (pt->core.EC == 1) {

    	ec = &pt->Disk.ec;
    	freq_array_obs=pt->Disk.ec.spec.nu_obs;
    	nuFnu_obs_array=pt->Disk.ec.spec.nuFnu_obs;
    	freq_array=pt->Disk.ec.spec.nu;
    	nu_seed_max =  pt->Disk.spec.nu_max;
    	nu_start_EC = &(pt->Disk.ec.spec.nu_min);
    	nu_stop_EC = &(pt->Disk.ec.spec.nu_max);
		nu_start_EC_obs = &(pt->Disk.ec.spec.nu_min_obs);
    	nu_stop_EC_obs = &(pt->Disk.ec.spec.nu_max_obs);
    	NU_INT_STOP_EC= &(pt->Disk.ec.NU_INT_STOP);
		R_ext_emit= pt->Disk.R_ext;
		if (pt->core.verbose>0) {
    		printf("nu_star_Disk=%e    nu_stop_Disk=%e\n",
    			pt->Disk.spec.nu_min,
    			pt->Disk.spec.nu_max);
    		printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
    		printf("-----------------------------------------------------------------\n");
    	}
    }else if (pt->core.EC == 2) {
    	ec = &pt->BLR.ec;
    	freq_array_obs=pt->BLR.ec.spec.nu_obs;
    	nuFnu_obs_array=pt->BLR.ec.spec.nuFnu_obs;
    	freq_array=pt->BLR.ec.spec.nu;
    	nu_seed_max =  pt->BLR.spec.nu_max;
    	nu_start_EC = &(pt->BLR.ec.spec.nu_min);
    	nu_stop_EC = &(pt->BLR.ec.spec.nu_max);
    	nu_start_EC_obs = &(pt->BLR.ec.spec.nu_min_obs);
    	nu_stop_EC_obs = &(pt->BLR.ec.spec.nu_max_obs);
    	NU_INT_STOP_EC= &(pt->BLR.ec.NU_INT_STOP);
		R_ext_emit = pt->BLR.R_BLR_out;
		if (pt->core.verbose>0) {
			printf("nu_star_BLR=%e    nu_stop_BLR=%e\n",
					pt->BLR.spec.nu_min,
					pt->BLR.spec.nu_max);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}
    }
    else if (pt->core.EC == 3) {
    	ec = &pt->DT.ec;
    	freq_array_obs=pt->DT.ec.spec.nu_obs;
    	nuFnu_obs_array=pt->DT.ec.spec.nuFnu_obs;
    	freq_array=pt->DT.ec.spec.nu;
    	nu_seed_max =  pt->DT.spec.nu_max;
    	nu_start_EC = &(pt->DT.ec.spec.nu_min);
    	nu_stop_EC = &(pt->DT.ec.spec.nu_max);
    	nu_start_EC_obs = &(pt->DT.ec.spec.nu_min_obs);
    	nu_stop_EC_obs = &(pt->DT.ec.spec.nu_max_obs);
    	NU_INT_STOP_EC= &(pt->DT.ec.NU_INT_STOP);
		R_ext_emit = pt->DT.R_DT;
		//printf("pre R_H=%e R_DT=%e EC_Stat=%d\n",pt->core.R_H,pt->DT.R_DT,pt->core.EC_stat);
    	if (pt->core.verbose>0) {
			printf("nu_star_DT=%e    nu_stop_DT=%e\n",
					pt->DT.spec.nu_min,
					pt->DT.spec.nu_max);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}

    }
    else if (pt->core.EC == 4) {
    	ec = &pt->Star.ec;
    	freq_array_obs=pt->Star.ec.spec.nu_obs;
    	nuFnu_obs_array=pt->Star.ec.spec.nuFnu_obs;
    	freq_array=pt->Star.ec.spec.nu;
    	nu_seed_max =  pt->Star.spec.nu_max;
    	nu_start_EC = &(pt->Star.ec.spec.nu_min);
    	nu_stop_EC = &(pt->Star.ec.spec.nu_max);
    	nu_start_EC_obs = &(pt->Star.ec.spec.nu_min_obs);
    	nu_stop_EC_obs = &(pt->Star.ec.spec.nu_max_obs);
    	NU_INT_STOP_EC= &(pt->Star.ec.NU_INT_STOP);
		R_ext_emit = -1;
		if (pt->core.verbose>0) {
			printf("nu_start_Star=%e    nu_stop_Star=%e\n",
					pt->Star.spec.nu_min,
					pt->Star.spec.nu_max);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}

    }
	    else if (pt->core.EC == 5) {
	    	ec = &pt->CMB.ec;
	    	freq_array_obs=pt->CMB.ec.spec.nu_obs;
    	nuFnu_obs_array=pt->CMB.ec.spec.nuFnu_obs;
    	freq_array=pt->CMB.ec.spec.nu;
    	nu_seed_max =  pt->CMB.spec.nu_max;
    	nu_start_EC = &(pt->CMB.ec.spec.nu_min);
    	nu_stop_EC = &(pt->CMB.ec.spec.nu_max);
    	nu_start_EC_obs = &(pt->CMB.ec.spec.nu_min_obs);
    	nu_stop_EC_obs = &(pt->CMB.ec.spec.nu_max_obs);
    	NU_INT_STOP_EC= &(pt->CMB.ec.NU_INT_STOP);
		R_ext_emit = -1;
		
		if (pt->core.verbose>0) {
    		printf("nu_start_CMB=%e    nu_stop_CMB=%e\n",
    				pt->CMB.spec.nu_min,
    				pt->CMB.spec.nu_max);
    		printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
	    		printf("-----------------------------------------------------------------\n");
	    	}

	    }
	    else if (pt->core.EC == 6) {
	    	ec = &pt->Corona.ec;
	    	freq_array_obs=pt->Corona.ec.spec.nu_obs;
	    	nuFnu_obs_array=pt->Corona.ec.spec.nuFnu_obs;
	    	freq_array=pt->Corona.ec.spec.nu;
	    	nu_seed_max =  pt->Corona.spec.nu_max;
	    	nu_start_EC = &(pt->Corona.ec.spec.nu_min);
	    	nu_stop_EC = &(pt->Corona.ec.spec.nu_max);
	    	nu_start_EC_obs = &(pt->Corona.ec.spec.nu_min_obs);
	    	nu_stop_EC_obs = &(pt->Corona.ec.spec.nu_max_obs);
	    	NU_INT_STOP_EC= &(pt->Corona.ec.NU_INT_STOP);
			R_ext_emit = pt->Corona.R_Corona;
			if (pt->core.verbose>0) {
	    		printf("nu_start_Corona=%e    nu_stop_Corona=%e\n",
	    				pt->Corona.spec.nu_min,
	    				pt->Corona.spec.nu_max);
	    		printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
	    		printf("-----------------------------------------------------------------\n");
	    	}

	    }else{
			printf("wrong EC \n ");
	        exit(0);
		
	}
	//printf("spettro_EC 1 R_H=%e c=%d \n", pt->core.R_H, set_condition_EC_correction(pt, pt->DT.R_DT));
	set_EC_stat_pre(pt, R_ext_emit);
	//printf("spettro_EC 2 R_H=%e c=%d \n", pt->core.R_H, set_condition_EC_correction(pt, pt->DT.R_DT));

	gmax=Find_gmax(pt,pt->emitters.Ne,pt->emitters.griglia_gamma_Ne_log);
	numax_KN = 1000 * gmax * MEC2 / HPLANCK;
	numax_TH = 1000 * (4.0 / 3.0) * pow(gmax, 2) * nu_seed_max;
	if (HPLANCK * nu_seed_max * gmax / MEC2 > 0.1) {
		*nu_stop_EC = numax_KN;
	} else {
		*nu_stop_EC = numax_TH;
	}


	if (HPLANCK * nu_seed_max * pt->emitters.Gamma_p3 / MEC2 > 0.1) {
		nu_peak= pt->emitters.Gamma_p3 * nu_seed_max;
	} else {
		nu_peak= (4.0 / 3.0) * pow(pt->emitters.Gamma_p3, 2) * nu_seed_max;
	}


   	//*nu_stop_EC = 100. * (4.0 / 3.0) * pow(gmax, 2) * nu_seed_max;

	if (pt->core.do_Sync>0){
		*nu_start_EC =  pt->Sync.spec.nu_peak_blob;
	}else{
		*nu_start_EC =  nu_seed_max/10;
	}

   	*nu_stop_EC_obs = nu_blob_to_nu_obs(*nu_stop_EC, pt->core.beam_obj, pt->core.z_cosm);
   	*nu_start_EC_obs=nu_blob_to_nu_obs(*nu_start_EC, pt->core.beam_obj, pt->core.z_cosm);


   	build_log_grid(*nu_start_EC,  *nu_stop_EC, pt->core.nu_IC_size, freq_array);
   	build_log_grid(*nu_start_EC_obs,  *nu_stop_EC_obs, pt->core.nu_IC_size, freq_array_obs);

	


   	if (pt->core.verbose>0) {
		printf("nu_start_EC=%e nu_stop_EC=%e nu_peak(estim.)=%e\n",
				*nu_start_EC,
				*nu_stop_EC,
				nu_peak);
		printf("nu_start_SSC=%e nu_stop_SSC=%e\n",
				pt->SSC.spec.nu_min,
				pt->SSC.spec.nu_max);
		printf("SSC=%d EC=%d TOT=%d\n",
				pt->core.SSC, pt->core.EC,
				pt->core.TOT);
		printf("emiss limit=%e\n", pt->core.emiss_lim);
	}

	I_MAX = pt->core.nu_IC_size -1;
	stop=0;

	eval_j_ptr = &eval_j_EC;
    threaded_j_evaluation(pt, eval_j_ptr, ec->spec.j_nu, freq_array, *nu_start_EC, *nu_stop_EC,I_MAX,pt->core.N_THREADS);

	for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++) {
        nuFnu_obs_array[NU_INT]=pt->core.emiss_lim;
        if (pt->core.verbose>1) {
            printf("#-> nu_em=%e  nu_obs=%e  i=%d\n", freq_array[NU_INT], freq_array_obs[NU_INT], NU_INT);
        }
        if ((freq_array[NU_INT] >= *nu_start_EC) && (freq_array[NU_INT] <= *nu_stop_EC)) {
			if (!stop) {
				
				if (pt->core.verbose > 1) {
					printf("#-> q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,
							pt->SSC.q_comp[NU_INT], NU_INT, ec->spec.j_nu[NU_INT],
							freq_array[NU_INT]);
				}

				//nu_src = nu_blob_to_nu_src(freq_array[NU_INT], pt->core.beam_obj,
				//		pt->core.z_cosm);
				if (pt->core.EC_stat == 1)
				{
					L_nu_EC = j_nu_src_to_L_nu_src(ec->spec.j_nu[NU_INT], pt->core.Vol_region,
												   pt->core.beam_obj);
				}
				else{
					L_nu_EC = j_nu_to_L_nu_src(ec->spec.j_nu[NU_INT], pt->core.Vol_region,
						pt->core.beam_obj);
				}
				//nuL_nu_EC = L_nu_EC * nu_src;
				F_nu_EC_obs = L_nu_src_to_F_nu(L_nu_EC, pt->core.beam_obj,
						pt->core.z_cosm, pt->core.dist);
                
				nuFnu_obs_array[NU_INT] = F_nu_EC_obs * freq_array_obs[NU_INT];
                

				if (ec->spec.j_nu[NU_INT] < pt->core.emiss_lim) {
					//out=0;
					if (freq_array[NU_INT] > numax_TH) {
						stop = 1;
					}
					nuFnu_obs_array[NU_INT] = pt->core.emiss_lim;
					ec->spec.j_nu[NU_INT] = pt->core.emiss_lim;
					pt->SSC.q_comp[NU_INT] = pt->core.emiss_lim;
				}
				//else{
				//	out=1;
				//}

				if (stop == 1 && freq_array[NU_INT] > numax_TH) {

					*nu_stop_EC_obs = freq_array_obs[NU_INT];
					*nu_stop_EC = freq_array[NU_INT];
					*NU_INT_STOP_EC = NU_INT;
					if (pt->core.verbose > 1) {
						printf("%e %d\n ", freq_array[NU_INT], NU_INT);
					}
				}
			}
            else{
				nuFnu_obs_array[NU_INT]=pt->core.emiss_lim;
				ec->spec.j_nu[NU_INT] = pt->core.emiss_lim;
				pt->SSC.q_comp[NU_INT]=pt->core.emiss_lim;

			}
            //===========================================
            // FILES output nu dnu nuFnu dnuFnu
            //===========================================
            /*
			if (pt->core.WRITE_TO_FILE==1){
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
            if (pt->core.verbose>1) {
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
        if (pt->core.verbose > 1) {
            printf("%e %d\n ", freq_array[NU_INT-1], NU_INT-1);
        }
    	if (pt->core.verbose>0) {
			if (pt->core.EC == 1) {

				printf("nu_stop_EC_Disk=%e NU_INT_STOP_EC_Disk=%d\n", pt->Disk.ec.spec.nu_max, pt->Disk.ec.NU_INT_STOP);
			}
			if (pt->core.EC == 2) {

				printf("nu_stop_EC_BLR=%e NU_INT_STOP_EC_BLR=%d\n", pt->BLR.ec.spec.nu_max, pt->BLR.ec.NU_INT_STOP);
			}
			if (pt->core.EC == 3) {

				printf("nu_stop_EC_DT=%e NU_INT_STOP_EC_DT=%d\n", pt->DT.ec.spec.nu_max, pt->DT.ec.NU_INT_STOP);
			}
			if (pt->core.EC == 4) {

				printf("nu_stop_EC_Star=%e NU_INT_STOP_EC_Star=%d\n", pt->Star.ec.spec.nu_max, pt->Star.ec.NU_INT_STOP);
			}
				if (pt->core.EC == 5) {

					printf("nu_stop_EC_CMB=%e NU_INT_STOP_EC_CMB=%d\n", pt->CMB.ec.spec.nu_max, pt->CMB.ec.NU_INT_STOP);
				}
				if (pt->core.EC == 6) {

					printf("nu_stop_EC_Corona=%e NU_INT_STOP_EC_Corona=%d\n", pt->Corona.ec.spec.nu_max, pt->Corona.ec.NU_INT_STOP);
				}
            //if (pt->core.EC == 6) {

            //    printf("nu_stop_EC_CMB_stat=%e NU_INT_STOP_EC_CMB_stat=%d\n", pt->nu_stop_EC_CMB_stat, pt->NU_INT_STOP_EC_CMB_stat);
            //}
    	}
    }
	//printf("spettro_EC 3 R_H=%e c=%d \n", pt->core.R_H, set_condition_EC_correction(pt, pt->DT.R_DT));

	set_EC_stat_post(pt);
	//printf("spettro_EC 4 R_H=%e c=%d \n", pt->core.R_H, set_condition_EC_correction(pt, pt->DT.R_DT));

	//spectra_External_Fields(1, pt,1);
	//printf("spettro_EC 5 R_H=%e c=%d \n", pt->core.R_H, set_condition_EC_correction(pt, pt->DT.R_DT));

	//===========================================
	//    trova nu peak e Flux peak
	//===========================================
	if (pt->core.EC == 1)
	{
		FindEpSp(freq_array, nuFnu_obs_array, pt->Disk.ec.NU_INT_STOP, pt,
					&(pt->Disk.ec.spec.nu_peak_obs),
					&(pt->Disk.ec.spec.nu_peak_src),
					&(pt->Disk.ec.spec.nu_peak_blob),
					&(pt->Disk.ec.spec.nuFnu_peak_obs),
					&(pt->Disk.ec.spec.nuLnu_peak_src),
					&(pt->Disk.ec.spec.nuLnu_peak_blob));

		if (pt->core.verbose > 0)
		{
			printf("nu_stop_EC_Disk=%e NU_INT_STOP_EC_Disk=%d\n", pt->Disk.ec.spec.nu_max, pt->Disk.ec.NU_INT_STOP);
			printf("EC Disk ");
			printf("nu_EC_blob peak=%e\n", pt->Disk.ec.spec.nu_peak_blob);
			printf("nu_EC_src  peak=%e\n", pt->Disk.ec.spec.nu_peak_src);
			printf("nu_EC_obs  peak=%e\n", pt->Disk.ec.spec.nu_peak_obs);

			printf("nuFnu EC  blob    peak=%e\n", pt->Disk.ec.spec.nuFnu_peak_obs);
			printf("nuLnu EC  src     peak=%e\n", pt->Disk.ec.spec.nuLnu_peak_src);
			printf("nuLnu EC  obs     peak=%e\n", pt->Disk.ec.spec.nuLnu_peak_blob);
		}
	}

	if (pt->core.EC == 2)
	{
		FindEpSp(pt->BLR.ec.spec.nu, nuFnu_obs_array, pt->BLR.ec.NU_INT_STOP, pt,
					&(pt->BLR.ec.spec.nu_peak_obs),
					&(pt->BLR.ec.spec.nu_peak_src),
					&(pt->BLR.ec.spec.nu_peak_blob),
					&(pt->BLR.ec.spec.nuFnu_peak_obs),
					&(pt->BLR.ec.spec.nuLnu_peak_src),
					&(pt->BLR.ec.spec.nuLnu_peak_blob));
		if (pt->core.verbose > 0)
		{
			printf("nu_stop_EC_BLR=%e NU_INT_STOP_EC_BLR=%d\n", pt->BLR.ec.spec.nu_max, pt->BLR.ec.NU_INT_STOP);
			printf("EC BLR ");
			printf("nu_EC_blob peak=%e\n", pt->BLR.ec.spec.nu_peak_blob);
			printf("nu_EC_src  peak=%e\n", pt->BLR.ec.spec.nu_peak_src);
			printf("nu_EC_obs  peak=%e\n", pt->BLR.ec.spec.nu_peak_obs);

			printf("nuFnu EC  blob    peak=%e\n", pt->BLR.ec.spec.nuFnu_peak_obs);
			printf("nuLnu EC  src     peak=%e\n", pt->BLR.ec.spec.nuLnu_peak_src);
			printf("nuLnu EC  obs     peak=%e\n", pt->BLR.ec.spec.nuLnu_peak_blob);
		}
	}

	if (pt->core.EC == 3)
	{
		FindEpSp(pt->DT.ec.spec.nu, nuFnu_obs_array, pt->DT.ec.NU_INT_STOP, pt,
					&(pt->DT.ec.spec.nu_peak_obs),
					&(pt->DT.ec.spec.nu_peak_src),
					&(pt->DT.ec.spec.nu_peak_blob),
					&(pt->DT.ec.spec.nuFnu_peak_obs),
					&(pt->DT.ec.spec.nuLnu_peak_src),
					&(pt->DT.ec.spec.nuLnu_peak_blob));
		if (pt->core.verbose > 0)
		{
			printf("nu_stop_EC_DT=%e NU_INT_STOP_EC_DT=%d\n", pt->DT.ec.spec.nu_max, pt->DT.ec.NU_INT_STOP);
			printf("EC DT ");
			printf("nu_EC_blob peak=%e\n", pt->DT.ec.spec.nu_peak_blob);
			printf("nu_EC_src  peak=%e\n", pt->DT.ec.spec.nu_peak_src);
			printf("nu_EC_obs  peak=%e\n", pt->DT.ec.spec.nu_peak_obs);

			printf("nuFnu EC  blob    peak=%e\n", pt->DT.ec.spec.nuFnu_peak_obs);
			printf("nuLnu EC  src     peak=%e\n", pt->DT.ec.spec.nuLnu_peak_src);
			printf("nuLnu EC  obs     peak=%e\n", pt->DT.ec.spec.nuLnu_peak_blob);
			}
		}

	if (pt->core.EC == 6)
	{
		FindEpSp(pt->Corona.ec.spec.nu, nuFnu_obs_array, pt->Corona.ec.NU_INT_STOP, pt,
					&(pt->Corona.ec.spec.nu_peak_obs),
					&(pt->Corona.ec.spec.nu_peak_src),
					&(pt->Corona.ec.spec.nu_peak_blob),
					&(pt->Corona.ec.spec.nuFnu_peak_obs),
					&(pt->Corona.ec.spec.nuLnu_peak_src),
					&(pt->Corona.ec.spec.nuLnu_peak_blob));
		if (pt->core.verbose > 0)
		{
			printf("nu_stop_EC_Corona=%e NU_INT_STOP_EC_Corona=%d\n", pt->Corona.ec.spec.nu_max, pt->Corona.ec.NU_INT_STOP);
			printf("EC Corona ");
			printf("nu_EC_blob peak=%e\n", pt->Corona.ec.spec.nu_peak_blob);
			printf("nu_EC_src  peak=%e\n", pt->Corona.ec.spec.nu_peak_src);
			printf("nu_EC_obs  peak=%e\n", pt->Corona.ec.spec.nu_peak_obs);

			printf("nuFnu EC  blob    peak=%e\n", pt->Corona.ec.spec.nuFnu_peak_obs);
			printf("nuLnu EC  src     peak=%e\n", pt->Corona.ec.spec.nuLnu_peak_src);
			printf("nuLnu EC  obs     peak=%e\n", pt->Corona.ec.spec.nuLnu_peak_blob);
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
        thread_args->blob_pt->SSC.q_comp[NU_INT] = 0.;
        thread_args->j_array[NU_INT] = 0.;
       
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->1 in eval_j_EC   NU_INT=%d eval_j_EC  nu_1=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }

		thread_args->blob_pt->SSC.q_comp[NU_INT] = rate_compton_GR(thread_args->blob_pt,nu_IC_out);
		if (thread_args->blob_pt->core.EC_stat == 1){
			//in this case we have q_comp in the disk frame, so j_nu is in the disk rest frame
			//and we have to use also the scattered nu in the disk rest frame
			thread_args->j_array[NU_INT]=thread_args->blob_pt->SSC.q_comp[NU_INT]*HPLANCK*thread_args->nu_array[NU_INT]*thread_args->blob_pt->core.beam_obj;
		}
		else{
			thread_args->j_array[NU_INT] = thread_args->blob_pt->SSC.q_comp[NU_INT] * HPLANCK * thread_args->nu_array[NU_INT];
		}
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->2 in  eval_j_EC NU_INT=%d q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,NU_INT,
                        thread_args->blob_pt->SSC.q_comp[NU_INT], NU_INT, thread_args->j_array[NU_INT],
                        thread_args->nu_array[NU_INT]);
        }
    }
	return NULL; 
}
