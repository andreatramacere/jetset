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


void spettro_EC(int Num_file, struct spettro *pt) {
    double  nu_src;
    double L_nu_EC, nuL_nu_EC, F_nu_EC_obs,nu_peak;
    double  gmax,numax_KN,numax_TH;
    double j_nu_disk;
    double * freq_array, *freq_array_obs;
    double * nuFnu_obs_array;
    double * nu_start_EC, * nu_stop_EC, * nu_start_EC_obs, * nu_stop_EC_obs, nu_seed_max;
    unsigned long * NU_INT_STOP_EC;
    unsigned long l, NU_INT, I_MAX, stop,out;
    char f_EC[static_file_name_max_legth];
    FILE *fp_EC;

    //=================================
    // apre i files dove scrive i dati di SSC
    // HEADER FILES
    //=================================
    if (pt->EC == 1) {
    	sprintf(f_EC, "%s%s-EC-Diks.dat",
    	                        pt->path, pt->STEM);
    }

    if (pt->EC == 2) {
        sprintf(f_EC, "%s%s-EC-BLR.dat",
                pt->path, pt->STEM);
    }

    if (pt->EC == 3) {
        sprintf(f_EC, "%s%s-EC-DT.dat",
                pt->path, pt->STEM);
    }
    if (pt->EC == 4) {
    	sprintf(f_EC, "%s%s-EC-Star.dat",
                    pt->path, pt->STEM);
    }
    if (pt->EC == 5) {
        	sprintf(f_EC, "%s%s-EC-CMB.dat",
                        pt->path, pt->STEM);
        }
    if (pt->EC == 6) {
            sprintf(f_EC, "%s%s-EC-CMB-stat.dat",
                        pt->path, pt->STEM);
        }
    fp_EC = fopen(f_EC, "w");
    if (fp_EC == NULL) {
        printf("non posso aprire %s\n ", f_EC);
        exit(1);
    }

    flux_header(fp_EC);

    //==================================================================


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
    	if (pt->verbose>0) {
    		printf("nu_star_Disk=%e    nu_stop_Disk=%e\n",
    			pt->nu_start_Disk,
    			pt->nu_stop_Disk);
    		printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
    		printf("-----------------------------------------------------------------\n");
    	}
    }
    if (pt->EC == 2) {
    	freq_array_obs=pt->nu_EC_BLR_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_BLR_obs;
    	freq_array=pt->nu_EC_BLR;
    	nu_seed_max =  pt->nu_stop_BLR;
    	nu_start_EC = &(pt->nu_start_EC_BLR);
    	nu_stop_EC = &(pt->nu_stop_EC_BLR);
    	nu_start_EC_obs = &(pt->nu_start_EC_BLR_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_BLR_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_BLR);
    	if (pt->verbose>0) {
			printf("nu_star_BLR=%e    nu_stop_BLR=%e\n",
					pt->nu_start_BLR,
					pt->nu_stop_BLR);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}
    }
    if (pt->EC == 3) {
    	freq_array_obs=pt->nu_EC_DT_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_DT_obs;
    	freq_array=pt->nu_EC_DT;
    	nu_seed_max =  pt->nu_stop_DT;
    	nu_start_EC = &(pt->nu_start_EC_DT);
    	nu_stop_EC = &(pt->nu_stop_EC_DT);
    	nu_start_EC_obs = &(pt->nu_start_EC_DT_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_DT_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_DT);
    	if (pt->verbose>0) {
			printf("nu_star_DT=%e    nu_stop_DT=%e\n",
					pt->nu_start_DT,
					pt->nu_stop_DT);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}

    }
    if (pt->EC == 4) {
    	freq_array_obs=pt->nu_EC_Star_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_Star_obs;
    	freq_array=pt->nu_EC_Star;
    	nu_seed_max =  pt->nu_stop_Star;
    	nu_start_EC = &(pt->nu_start_EC_Star);
    	nu_stop_EC = &(pt->nu_stop_EC_Star);
    	nu_start_EC_obs = &(pt->nu_start_EC_Star_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_Star_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_Star);
    	if (pt->verbose>0) {
			printf("nu_start_Star=%e    nu_stop_Star=%e\n",
					pt->nu_start_Star,
					pt->nu_stop_Star);
			printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
			printf("-----------------------------------------------------------------\n");
    	}

    }
    if (pt->EC == 5) {
    	freq_array_obs=pt->nu_EC_CMB_obs;
    	nuFnu_obs_array=pt->nuF_nu_EC_CMB_obs;
    	freq_array=pt->nu_EC_CMB;
    	nu_seed_max =  pt->nu_stop_CMB;
    	nu_start_EC = &(pt->nu_start_EC_CMB);
    	nu_stop_EC = &(pt->nu_stop_EC_CMB);
    	nu_start_EC_obs = &(pt->nu_start_EC_CMB_obs);
    	nu_stop_EC_obs = &(pt->nu_stop_EC_CMB_obs);
    	NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_CMB);
    	if (pt->verbose>0) {
    		printf("nu_start_CMB=%e    nu_stop_CMB=%e\n",
    				pt->nu_start_CMB,
    				pt->nu_stop_CMB);
    		printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
    		printf("-----------------------------------------------------------------\n");
    	}

    }
    if (pt->EC == 6) {
        freq_array_obs=pt->nu_EC_CMB_stat_obs;
        nuFnu_obs_array=pt->nuF_nu_EC_CMB_stat_obs;
        freq_array=pt->nu_EC_CMB_stat;
        nu_seed_max =  pt->nu_stop_CMB_stat;
        nu_start_EC = &(pt->nu_start_EC_CMB_stat);
        nu_stop_EC = &(pt->nu_stop_EC_CMB_stat);
        nu_start_EC_obs = &(pt->nu_start_EC_CMB_stat_obs);
        nu_stop_EC_obs = &(pt->nu_stop_EC_CMB_stat_obs);
        NU_INT_STOP_EC= &(pt->NU_INT_STOP_EC_CMB_stat);
        if (pt->verbose>0) {
            printf("nu_start_CMB=%e    nu_stop_CMB=%e\n",
                    pt->nu_start_CMB_stat,
                    pt->nu_stop_CMB_stat);
            printf("these freq. are boosted from the DISK frame  into the BLOB frame\n");
            printf("-----------------------------------------------------------------\n");
        }


    }

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

    //pt->nu_stop_compton = *nu_start_EC_obs;
    //pt->nu_start_compton = *nu_stop_EC_obs;


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

   	}

	I_MAX = pt->nu_IC_size-1;
	stop=0;
	for (NU_INT = 0; NU_INT <= I_MAX; NU_INT++) {

        pt->nu_1 = freq_array[NU_INT];
        nuFnu_obs_array[NU_INT]=pt->emiss_lim;
        pt->j_EC[NU_INT] = pt->emiss_lim;

        if (pt->verbose>1) {
            printf("#-> nu_em=%e  nu_obs=%e  i=%d\n", freq_array[NU_INT], freq_array_obs[NU_INT], NU_INT);
        }
        if ((freq_array[NU_INT] >= *nu_start_EC) && (freq_array[NU_INT] <= *nu_stop_EC)) {
			if (!stop) {

				pt->q_comp[NU_INT] = rate_compton_GR(pt);
                if (pt->EC == 6){
                    //in this case we have q_comp in the disk frame, so j_nu is in the disk rest frame
                    //and have to get the out nu in the disk rest frame 
                    j_nu_disk=pt->q_comp[NU_INT]*HPLANCK*freq_array[NU_INT]*pt->beam_obj;
                    
                    //now we go back to the blob 
                    //this is the j_nu in the blob frame at nu_blob 
                    //evaluated from j_disk at nu_disk   
                    pt->j_EC[NU_INT]=j_nu_disk/(pt->beam_obj*pt->beam_obj);

                }
                else{

                    pt->j_EC[NU_INT] = pt->q_comp[NU_INT] *
                    HPLANCK * freq_array[NU_INT];
                }

				if (pt->verbose > 1) {
					printf("#-> q_comp[%d]=%e j[%d]=%e nu_1=%e \n", NU_INT,
							pt->q_comp[NU_INT], NU_INT, pt->j_EC[NU_INT],
							freq_array[NU_INT]);
				}

				nu_src = nu_blob_to_nu_src(freq_array[NU_INT], pt->beam_obj,
						pt->z_cosm);
				L_nu_EC = j_nu_to_L_nu_src(pt->j_EC[NU_INT], pt->Vol_sphere,
						pt->beam_obj);
				nuL_nu_EC = L_nu_EC * nu_src;
				F_nu_EC_obs = L_nu_src_to_F_nu(L_nu_EC, pt->beam_obj,
						pt->z_cosm, pt->dist);
                
				nuFnu_obs_array[NU_INT] = F_nu_EC_obs * freq_array_obs[NU_INT];
                


				if (pt->j_EC[NU_INT] < pt->emiss_lim) {
					out=0;
					if (freq_array[NU_INT] > nu_peak) {
						stop = 1;
					}
					nuFnu_obs_array[NU_INT] = pt->emiss_lim;
					pt->j_EC[NU_INT] = pt->emiss_lim;
					pt->q_comp[NU_INT] = pt->emiss_lim;
				}
				else{
					out=1;
				}

				if (stop == 1 && freq_array[NU_INT] > nu_peak) {

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
            if (!stop && out) {
                fprintf(fp_EC, "%4.4e\t%4.4e\t%4.4e\t %4.4e\t%4.4e\t%4.4e\n",
                        log10(freq_array_obs[NU_INT]),
                        log10(F_nu_EC_obs * freq_array_obs[NU_INT]),
                        freq_array_obs[NU_INT],
                        F_nu_EC_obs*freq_array_obs[NU_INT],
                        nu_src,
                        nuL_nu_EC);
            }
            if (pt->verbose>1) {
                printf("#-> ********************************\n\n");
            }
            //==========================  END of Loop ove frequencies ====================================

        }
    }


    //Se ancora non ha trovato nu_stop
    if (!stop) {
    	*nu_stop_EC_obs= freq_array_obs[NU_INT];
    	*nu_stop_EC = freq_array[NU_INT];
    	*NU_INT_STOP_EC = NU_INT;
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
            if (pt->EC == 6) {

                printf("nu_stop_EC_CMB_stat=%e NU_INT_STOP_EC_CMB_stat=%d\n", pt->nu_stop_EC_CMB_stat, pt->NU_INT_STOP_EC_CMB_stat);
            }
    	}
    }

    fclose(fp_EC);

    //===========================================
    //    trova nu peak e Flux peak
    //===========================================
    if (pt->EC == 1) {
    	FindEpSp(freq_array, nuFnu_obs_array, pt->NU_INT_STOP_EC_Disk, pt,
                    &(pt->nu_peak_EC_Disk_obs),
                    &(pt->nu_peak_EC_Disk_src),
                    &(pt->nu_peak_EC_Disk_blob),
                    &(pt->nuFnu_peak_EC_Disk_obs),
                    &(pt->nuLnu_peak_EC_Disk_src),
                    &(pt->nuLnu_peak_EC_Disk_blob));

    	if (pt->verbose>0) {
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

    if (pt->EC == 2) {
        FindEpSp(pt->nu_EC_BLR, nuFnu_obs_array, pt->NU_INT_STOP_EC_BLR, pt,
                &(pt->nu_peak_EC_BLR_obs),
                &(pt->nu_peak_EC_BLR_src),
                &(pt->nu_peak_EC_BLR_blob),
                &(pt->nuFnu_peak_EC_BLR_obs),
                &(pt->nuLnu_peak_EC_BLR_src),
                &(pt->nuLnu_peak_EC_BLR_blob));
        if (pt->verbose>0) {
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

    if (pt->EC == 3) {
        FindEpSp(pt->nu_EC_DT, nuFnu_obs_array, pt->NU_INT_STOP_EC_DT, pt,
                &(pt->nu_peak_EC_DT_obs),
                &(pt->nu_peak_EC_DT_src),
                &(pt->nu_peak_EC_DT_blob),
                &(pt->nuFnu_peak_EC_DT_obs),
                &(pt->nuLnu_peak_EC_DT_src),
                &(pt->nuLnu_peak_EC_DT_blob));
        if (pt->verbose>0) {
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

    return;
}
//=========================================================================================







