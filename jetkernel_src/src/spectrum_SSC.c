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
   
    
    
    //============================================================
    //         inizio  loop sulle freq per spettro  compton
    //============================================================
    pt->TOT=0;
    pt->EC=0;
    pt->ord_comp=1;
    pt->SSC=pt->do_SSC;
    

    // massima e minima freq compton
    gmax=Find_gmax(pt,pt->Ne,pt->griglia_gamma_Ne_log);
    numax_KN=100000*gmax*MEC2/HPLANCK;
    numax_TH=100000*(4.0/3.0)*pow(gmax, 2)*pt->nu_stop_Sync;
    if (HPLANCK*pt->nu_stop_Sync*gmax/MEC2>0.1){
        pt->nu_stop_SSC=numax_KN;
    }
    else{
        pt->nu_stop_SSC=numax_TH;
    }
   

    if (HPLANCK *pt->nu_peak_Sync_blob * pt->Gamma_p3 / MEC2 > 0.1) {
		nu_peak = pt->Gamma_p3 * pt->nu_peak_Sync_blob;
	} else {
		nu_peak = (4.0 / 3.0) * pow(pt->Gamma_p3, 2) * pt->nu_peak_Sync_blob;
	}

    nu_min_TH_1=pt->gmax*pt->gmax*pt->nu_start_Sync;
    nu_min_TH_2=pt->gmin*pt->gmin*pt->nu_start_Sync;
    pt->nu_start_SSC=min(nu_min_TH_1,nu_min_TH_2);

    pt->nu_start_SSC_obs =nu_blob_to_nu_obs(pt->nu_start_SSC, pt->beam_obj, pt->z_cosm);
    pt->nu_stop_SSC_obs = nu_blob_to_nu_obs(pt->nu_stop_SSC, pt->beam_obj, pt->z_cosm);

    build_log_grid(pt->nu_start_SSC,  pt->nu_stop_SSC, pt->nu_IC_size, pt->nu_SSC);
    build_log_grid(pt->nu_start_SSC_obs,  pt->nu_stop_SSC_obs, pt->nu_IC_size, pt->nu_SSC_obs);



    


	I_MAX = pt->nu_IC_size-1;
	if (pt->verbose>0) {
		printf("**********************  SSC spectrum 1st Order   ****************************\n");
		printf("gmax from Ne>0 = %e\n", gmax);
		printf("nu_star_Sync_ssc=%e nu_stop_Sync_ssc=%e\n",
				pt->nu_start_Sync,
				pt->nu_stop_Sync_ssc);
		printf("nu_stop_compton_TH=%e nu_stop_compton_KN=%e\n",
				numax_TH,
				numax_KN);
		printf("nu_start_comp=%e nu_stop_compton=%e nu_peak(estim.)=%e\n",
				pt->nu_start_SSC,
				pt->nu_stop_SSC,
				nu_peak);
		printf("SSC=%d EC=%d TOT=%d\n",
				pt->SSC, pt->EC,
				pt->TOT);
		printf("Number of freq to eval=%d\n",I_MAX);
	}

	stop=0;
    for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++) {
        

        pt->nu_1=pt->nu_SSC[NU_INT];
        //pt->j_comp[NU_INT]=pt->emiss_lim;
        //pt->nuF_nu_SSC_obs[NU_INT]=pt->emiss_lim;
        if(pt->verbose>1){
            printf("#-> nu_em=%e  nu_obs=%e  i=%d\n", pt->nu_SSC[NU_INT], pt->nu_SSC_obs[NU_INT], NU_INT);
        }
        if((pt->nu_SSC[NU_INT]>=pt->nu_start_SSC) &&(pt->nu_SSC[NU_INT]<=pt->nu_stop_SSC)){
			if (!stop) {
				pt->q_comp[NU_INT] = rate_compton_GR(pt);
				pt->j_comp[NU_INT] = pt->q_comp[NU_INT] *
				HPLANCK * pt->nu_SSC[NU_INT];
				if (pt->verbose > 1) {
					printf("#-> q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,
							pt->q_comp[NU_INT], NU_INT, pt->j_comp[NU_INT],
							pt->nu_SSC[NU_INT]);
				}
				//nu_src = nu_blob_to_nu_src(pt->nu_SSC[NU_INT], pt->beam_obj,
				//		pt->z_cosm);
				L_nu_SSC = j_nu_to_L_nu_src(pt->j_comp[NU_INT], pt->Vol_sphere,
						pt->beam_obj);
				//nuL_nu_SSC = L_nu_SSC * nu_src;
				F_nu_SSC_obs = L_nu_src_to_F_nu(L_nu_SSC, pt->beam_obj,
						pt->z_cosm, pt->dist);
				pt->nuF_nu_SSC_obs[NU_INT] = F_nu_SSC_obs
						* pt->nu_SSC_obs[NU_INT];

				if (pt->verbose > 1) {
					printf("nu_stop_comp_SSC=%e NU_INT=%d\n ", pt->nu_stop_SSC,
							NU_INT);
				}

				//if (pt->j_comp[NU_INT] < pt->emiss_lim) {
				//	out = 0;

				//	if (pt->nu_SSC[NU_INT] > numax_TH) {
				//		stop = 1;
				//	}
                if (pt->j_comp[NU_INT] < pt->emiss_lim) {
                    pt->j_comp[NU_INT] = pt->emiss_lim;
                    pt->nuF_nu_SSC_obs[NU_INT] = pt->emiss_lim;
                }
				//}
                //
				//else {
				//	out = 1;
				//}

				//if (stop == 1 && pt->nu_SSC[NU_INT] > numax_TH) {
				//	F_nu_SSC_obs = pt->emiss_lim;
				//pt->nu_stop_SSC = pt->nu_SSC[NU_INT];
				//pt->nu_stop_SSC_obs = pt->nu_SSC_obs[NU_INT];
				pt->NU_INT_STOP_COMPTON_SSC = NU_INT;
				//	if (pt->verbose > 1) {
				//		printf("%e %d\n ", pt->nu_SSC[NU_INT], NU_INT);
				//	}
				//}
			}
            else{
				pt->j_comp[NU_INT]=pt->emiss_lim;
				pt->q_comp[NU_INT]=pt->emiss_lim;
				pt->nuF_nu_SSC_obs[NU_INT]=pt->emiss_lim;
			 }
            
           
            if(pt->verbose>1){
                printf("#-> ********************************\n\n");
            }
            //==========================  END of Loop ove frequencies ====================================
        }
    }
    //Se ancora non ha trovato nu_stop
    //if (!stop){
    //    pt->nu_stop_SSC = pt->nu_SSC[NU_INT-1];
    //    pt->nu_stop_SSC_obs = pt->nu_SSC_obs[NU_INT-1];
    //    pt->NU_INT_STOP_COMPTON_SSC = NU_INT-1;
    //    if (pt->verbose > 1) {
    //        printf("%e %d\n ", pt->nu_SSC[NU_INT-1], NU_INT-1);
    //    }
    
    
    //===========================================
    //    trova nu peak e Flux peak
    //===========================================
    FindEpSp(pt->nu_SSC, pt->nuF_nu_SSC_obs,   pt->NU_INT_STOP_COMPTON_SSC, pt,
            &(pt->nu_peak_SSC_obs),
            &(pt->nu_peak_SSC_src),
            &(pt->nu_peak_SSC_blob),
            &(pt->nuFnu_peak_SSC_obs),
            &(pt->nuLnu_peak_SSC_src),
            &(pt->nuLnu_peak_SSC_blob));
    
    if (pt->verbose>0) {
        printf("nu_stop=%e NU_INT_STOP_COMPTON_SSC=%d\n", pt->nu_stop_SSC, pt->NU_INT_STOP_COMPTON_SSC);

        printf("nu_SSC_blob peak=%e\n", pt->nu_peak_SSC_blob);
        printf("nu_SSC_src   peak=%e\n", pt->nu_peak_SSC_src);
        printf("nu_SSC_obs  peak=%e\n", pt->nu_peak_SSC_obs);
        
        printf("nuFnu SSC  blob    peak=%e\n", pt->nuFnu_peak_SSC_obs);
        printf("nuLnu SSC  src      peak=%e\n", pt->nuLnu_peak_SSC_src);
        printf("nuLnu SSC  obs     peak=%e\n", pt->nuLnu_peak_SSC_blob);
    }
    
   
    return ;
}
//=========================================================================================







