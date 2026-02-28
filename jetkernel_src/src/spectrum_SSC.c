//=========================================================================================
//                   CALCOLO DELLO SPETTRO COMPTON
//=========================================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"
#include <pthread.h>

/**
 * \file spettro_Compton.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief CALCOLO DELLO SPETTRO COMPTON
 *
 */


void spettro_compton(int Num_file, struct blob *pt){
    double nu_peak;
    double L_nu_SSC, F_nu_SSC_obs;
    double gmax,numax_KN,numax_TH,nu_min_TH_1,nu_min_TH_2;
    unsigned int NU_INT, I_MAX, stop;
    void *(*eval_j_ptr)(void * args);
    
    
    //============================================================
    //         inizio  loop sulle freq per spettro  compton
    //============================================================
    pt->core.TOT=0;
    pt->core.EC=0;
    pt->core.ord_comp=1;
    pt->core.SSC=pt->core.do_SSC;
    

    // massima e minima freq compton
    gmax=Find_gmax(pt,pt->emitters.Ne,pt->emitters.griglia_gamma_Ne_log);
    numax_KN=100000*gmax*MEC2/HPLANCK;
    numax_TH=100000*(4.0/3.0)*pow(gmax, 2)*pt->Sync.spec.nu_max;
    if (HPLANCK*pt->Sync.spec.nu_max*gmax/MEC2>0.1){
        pt->SSC.spec.nu_max=numax_KN;
    }
    else{
        pt->SSC.spec.nu_max=numax_TH;
    }
   

    if (HPLANCK *pt->Sync.spec.nu_peak_blob * pt->emitters.Gamma_p3 / MEC2 > 0.1) {
		nu_peak = pt->emitters.Gamma_p3 * pt->Sync.spec.nu_peak_blob;
	} else {
		nu_peak = (4.0 / 3.0) * pow(pt->emitters.Gamma_p3, 2) * pt->Sync.spec.nu_peak_blob;
	}

    nu_min_TH_1=pt->emitters.gmax*pt->emitters.gmax*pt->Sync.spec.nu_min;
    nu_min_TH_2=pt->emitters.gmin*pt->emitters.gmin*pt->Sync.spec.nu_min;
    pt->SSC.spec.nu_min=min(nu_min_TH_1,nu_min_TH_2);

    pt->SSC.spec.nu_min_obs =nu_blob_to_nu_obs(pt->SSC.spec.nu_min, pt->core.beam_obj, pt->core.z_cosm);
    pt->SSC.spec.nu_max_obs = nu_blob_to_nu_obs(pt->SSC.spec.nu_max, pt->core.beam_obj, pt->core.z_cosm);

    build_log_grid(pt->SSC.spec.nu_min,  pt->SSC.spec.nu_max, pt->core.nu_IC_size, pt->SSC.spec.nu);
    build_log_grid(pt->SSC.spec.nu_min_obs,  pt->SSC.spec.nu_max_obs, pt->core.nu_IC_size, pt->SSC.spec.nu_obs);


	I_MAX = pt->core.nu_IC_size-1;
	if (pt->core.verbose>0) {
		printf("**********************  SSC spectrum 1st Order   ****************************\n");
		printf("gmax from Ne>0 = %e\n", gmax);
		printf("nu_star_Sync_ssc=%e nu_stop_Sync_ssc=%e\n",
				pt->Sync.spec.nu_min,
				pt->Sync.nu_stop_Sync_ssc);
		printf("nu_stop_compton_TH=%e nu_stop_compton_KN=%e\n",
				numax_TH,
				numax_KN);
		printf("nu_start_comp=%e nu_stop_compton=%e nu_peak(estim.)=%e\n",
				pt->SSC.spec.nu_min,
				pt->SSC.spec.nu_max,
				nu_peak);
		printf("SSC=%d EC=%d TOT=%d\n",
				pt->core.SSC, pt->core.EC,
				pt->core.TOT);
		printf("Number of freq to eval=%d\n",I_MAX);
	}

	stop=0;
    eval_j_ptr = &eval_j_SSC;
    threaded_j_evaluation(pt, eval_j_ptr, pt->SSC.spec.j_nu,pt->SSC.spec.nu,pt->SSC.spec.nu_min, pt->SSC.spec.nu_max,I_MAX,pt->core.N_THREADS);
    for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++) {
        
        //pt->core.nu_1=pt->SSC.spec.nu[NU_INT];
        if(pt->core.verbose>1){
            printf("#-> nu_em=%e  nu_obs=%e  i=%d\n", pt->SSC.spec.nu[NU_INT], pt->SSC.spec.nu_obs[NU_INT], NU_INT);
        }
        if((pt->SSC.spec.nu[NU_INT]>=pt->SSC.spec.nu_min) &&(pt->SSC.spec.nu[NU_INT]<=pt->SSC.spec.nu_max)){
			if (!stop) {
				if (pt->core.verbose > 1) {
					printf("#-> q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,
							pt->SSC.q_comp[NU_INT], NU_INT, pt->SSC.spec.j_nu[NU_INT],
							pt->SSC.spec.nu[NU_INT]);
				}
				//nu_src = nu_blob_to_nu_src(pt->SSC.spec.nu[NU_INT], pt->core.beam_obj,
				//		pt->core.z_cosm);
				L_nu_SSC = j_nu_to_L_nu_src(pt->SSC.spec.j_nu[NU_INT], pt->core.Vol_region,
						pt->core.beam_obj);
				//nuL_nu_SSC = L_nu_SSC * nu_src;
				F_nu_SSC_obs = L_nu_src_to_F_nu(L_nu_SSC, pt->core.beam_obj,
						pt->core.z_cosm, pt->core.dist);
				pt->SSC.spec.nuFnu_obs[NU_INT] = F_nu_SSC_obs
						* pt->SSC.spec.nu_obs[NU_INT];

				if (pt->core.verbose > 1) {
					printf("nu_stop_comp_SSC=%e NU_INT=%d\n ", pt->SSC.spec.nu_max,
							NU_INT);
				}

			
                if (pt->SSC.spec.j_nu[NU_INT] < pt->core.emiss_lim) {
                    pt->SSC.spec.j_nu[NU_INT] = pt->core.emiss_lim;
                    pt->SSC.spec.nuFnu_obs[NU_INT] = pt->core.emiss_lim;
                }
				
				pt->SSC.NU_INT_STOP_COMPTON_SSC = NU_INT;
				
			}
            else{
				pt->SSC.spec.j_nu[NU_INT]=pt->core.emiss_lim;
				pt->SSC.q_comp[NU_INT]=pt->core.emiss_lim;
				pt->SSC.spec.nuFnu_obs[NU_INT]=pt->core.emiss_lim;
			 }
            
           
            if(pt->core.verbose>1){
                printf("#-> ********************************\n\n");
            }
            //==========================  END of Loop ove frequencies ====================================
        }
    }
    //Se ancora non ha trovato nu_stop
    //if (!stop){
    //    pt->SSC.spec.nu_max = pt->SSC.spec.nu[NU_INT-1];
    //    pt->SSC.spec.nu_max_obs = pt->SSC.spec.nu_obs[NU_INT-1];
    //    pt->SSC.NU_INT_STOP_COMPTON_SSC = NU_INT-1;
    //    if (pt->core.verbose > 1) {
    //        printf("%e %d\n ", pt->SSC.spec.nu[NU_INT-1], NU_INT-1);
    //    }
    
    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================
    FindEpSp(pt->SSC.spec.nu, pt->SSC.spec.nuFnu_obs,   pt->SSC.NU_INT_STOP_COMPTON_SSC, pt,
            &(pt->SSC.spec.nu_peak_obs),
            &(pt->SSC.spec.nu_peak_src),
            &(pt->SSC.spec.nu_peak_blob),
            &(pt->SSC.spec.nuFnu_peak_obs),
            &(pt->SSC.spec.nuLnu_peak_src),
            &(pt->SSC.spec.nuLnu_peak_blob));
    
    if (pt->core.verbose>0) {
        printf("nu_stop=%e NU_INT_STOP_COMPTON_SSC=%d\n", pt->SSC.spec.nu_max, pt->SSC.NU_INT_STOP_COMPTON_SSC);

        printf("nu_SSC_blob peak=%e\n", pt->SSC.spec.nu_peak_blob);
        printf("nu_SSC_src   peak=%e\n", pt->SSC.spec.nu_peak_src);
        printf("nu_SSC_obs  peak=%e\n", pt->SSC.spec.nu_peak_obs);
        
        printf("nuFnu SSC  blob    peak=%e\n", pt->SSC.spec.nuFnu_peak_obs);
        printf("nuLnu SSC  src      peak=%e\n", pt->SSC.spec.nuLnu_peak_src);
        printf("nuLnu SSC  obs     peak=%e\n", pt->SSC.spec.nuLnu_peak_blob);
    }
    
   
    return ;
}
//=========================================================================================


void  * eval_j_SSC(void *data){
    unsigned int NU_INT;
    double nu_IC_out;
    struct j_args *thread_args = data;
    for (NU_INT = thread_args->NU_INT_START; NU_INT <= thread_args->NU_INT_STOP; NU_INT++) {
        nu_IC_out=thread_args->nu_array[NU_INT];
        thread_args->blob_pt->SSC.q_comp[NU_INT] = 0.;
        thread_args->blob_pt->SSC.spec.j_nu[NU_INT] = 0.;
       
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->1 in eval_j_SSC   NU_INT=%d eval_j_SSC  nu_1=%e \n", NU_INT, thread_args->nu_array[NU_INT]);
        }
        
        thread_args->blob_pt->SSC.q_comp[NU_INT] = rate_compton_GR(thread_args->blob_pt,nu_IC_out);
        thread_args->blob_pt->SSC.spec.j_nu[NU_INT] = thread_args->blob_pt->SSC.q_comp[NU_INT] *HPLANCK * thread_args->nu_array[NU_INT];
        if (thread_args->blob_pt->core.verbose > 1) {
                printf("#->2 in  eval_j_SSC NU_INT=%d q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,NU_INT,
                        thread_args->blob_pt->SSC.q_comp[NU_INT], NU_INT, thread_args->blob_pt->SSC.spec.j_nu[NU_INT],
                        thread_args->nu_array[NU_INT]);
        }
    }
    //}
   return NULL; 
}