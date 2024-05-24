//=========================================================================================
//   FUNZIONI COMPTON
//
//=========================================================================================

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"
/**
 * \file funzioni_compton.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief funzioni per il Compton
 *
 */




//=========================================================================================
// Rate Compton
//=========================================================================================
double rate_compton_GR(struct blob *pt_GR, double nu_IC_out) {
    /**
     * \author Andrea Tramacere
     * \date 27-04-2004
     * \brief funzioni per il Compton
     *  CALCOLO DEL RATE COMPTON GENERALE metodo GRINDLAY 1985
     */
    double rate_comp=0;
    //double (*pf_K) (struct blob *, double x);
    double nu_IC_out_stat;

    //pf_K = &f_compton_K1;

    nu_IC_out_stat = nu_IC_out * pt_GR->beam_obj;

    if (pt_GR->verbose>1) {
        printf("GR\n");
        printf("#-> SSC=%d EC=%d\n", pt_GR->SSC, pt_GR->EC);
        printf("#-> gmin=%e gmax=%e\n", pt_GR->gmin, pt_GR->gmax);
    }

    //SSC
    if (nu_IC_out < pt_GR->nu_stop_SSC && pt_GR->ord_comp == 1) {
        if (pt_GR->SSC == 1 && pt_GR->EC == 0) {
            if (pt_GR->verbose>1) {
                printf("nu_start_Sync=%e\n", pt_GR->nu_start_Sync);
                printf("nu_stop_Sync_ssc=%e\n", pt_GR->nu_stop_Sync_ssc);
            }
            pt_GR->nu_seed = pt_GR->nu_Sync;
            pt_GR->n_seed = pt_GR->n_Sync;
            //pt_GR->griglia_gamma_log_IC=pt_GR->griglia_gamma_Ne_log;
            //pt_GR->N_IC=pt_GR->Ne;
            rate_comp = integrale_IC(pt_GR,
                    pt_GR->nu_start_Sync,
                    pt_GR->nu_stop_Sync_ssc,
                    0,
                    nu_IC_out);
        }
    }
    //EC Disk
    if (nu_IC_out < pt_GR->nu_stop_EC_Disk && pt_GR->ord_comp == 1) {
		if (pt_GR->SSC == 0 && pt_GR->EC == 1) {

			if (pt_GR->verbose>1) {
				printf("Disk\n");
				printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->nu_start_Disk);
                printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->nu_stop_Disk);
            }
            if (pt_GR->EC_stat == 0)
            {
                pt_GR->nu_seed = pt_GR->nu_Disk;
                pt_GR->n_seed = pt_GR->n_Disk;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_Disk,
                                         pt_GR->nu_stop_Disk,
                                         pt_GR->EC_stat,
                                         nu_IC_out);
            }
            else{
                pt_GR->nu_seed = pt_GR->nu_Disk_disk_RF;
                pt_GR->n_seed = pt_GR->n_Disk_DRF;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_Disk_DRF,
                                         pt_GR->nu_stop_Disk_DRF,
                                         pt_GR->EC_stat,
                                         nu_IC_out_stat);
            }
			
		}
    }
    //EC BLR
    if (nu_IC_out < pt_GR->nu_stop_EC_BLR && pt_GR->ord_comp == 1) {
    	if (pt_GR->SSC == 0 && pt_GR->EC == 2) {

            if (pt_GR->verbose>1) {
                printf("BLR\n");
                printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->nu_start_BLR);
                printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->nu_stop_BLR);
            }

            if (pt_GR->EC_stat == 0)
            {
                pt_GR->nu_seed = pt_GR->nu_BLR;
                pt_GR->n_seed = pt_GR->n_BLR;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_BLR,
                                         pt_GR->nu_stop_BLR,
                                         pt_GR->EC_stat,
                                         nu_IC_out);
            }
            else
            {             
                pt_GR->nu_seed = pt_GR->nu_BLR_disk_RF;
                pt_GR->n_seed = pt_GR->n_BLR_DRF;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_BLR_disk_RF,
                                         pt_GR->nu_stop_BLR_disk_RF,
                                         pt_GR->EC_stat,
                                         nu_IC_out_stat);

            }
        }
    }
    //EC DT
    if (nu_IC_out < pt_GR->nu_stop_EC_DT && pt_GR->ord_comp == 1) {
    	if (pt_GR->SSC == 0 && pt_GR->EC == 3) {
            if (pt_GR->verbose>1) {
                printf("DT\n");
                printf("(blob rest frame) nu_start_EC_seed DT=%e\n", pt_GR->nu_start_DT);
                printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->nu_stop_DT);
            }

            if (pt_GR->EC_stat == 0)
            {
                pt_GR->nu_seed = pt_GR->nu_DT;
                pt_GR->n_seed = pt_GR->n_DT;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_DT,
                                         pt_GR->nu_stop_DT,
                                         pt_GR->EC_stat,
                                         nu_IC_out);
            }
            else
            {
                pt_GR->nu_seed = pt_GR->nu_DT_disk_RF;
                pt_GR->n_seed = pt_GR->n_DT_DRF;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_DT,
                                         pt_GR->nu_stop_DT_DRF,
                                         pt_GR->EC_stat,
                                         nu_IC_out_stat);
            }
       }
    }

    //EC Star
    if (nu_IC_out < pt_GR->nu_stop_EC_Star && pt_GR->ord_comp == 1) {
    	if (pt_GR->SSC == 0 && pt_GR->EC == 4) {

		   if (pt_GR->verbose>1) {
			   printf("DT\n");
               printf("(blob rest frame) nu_start_EC_seed Star=%e\n", pt_GR->nu_start_Star);
               printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->nu_stop_Star);
           }
		   pt_GR->nu_seed = pt_GR->nu_Star;
		   pt_GR->n_seed = pt_GR->n_Star;
           if (pt_GR->EC_stat == 0)
           {
               pt_GR->nu_seed = pt_GR->nu_Star;
               pt_GR->n_seed = pt_GR->n_Star;
               rate_comp = integrale_IC(pt_GR,
                                        pt_GR->nu_start_Star,
                                        pt_GR->nu_stop_Star,
                                        pt_GR->EC_stat,
                                        nu_IC_out);
           }
           else
           {
               pt_GR->nu_seed = pt_GR->nu_Star_disk_RF;
               pt_GR->n_seed = pt_GR->n_Star_DRF;
               rate_comp = integrale_IC(pt_GR,
                                        pt_GR->nu_start_Star_DRF,
                                        pt_GR->nu_stop_Star_DRF,
                                        pt_GR->EC_stat,
                                        nu_IC_out_stat);
           }
       }
    }

    //EC CMB
    if (nu_IC_out < pt_GR->nu_stop_EC_CMB && pt_GR->ord_comp == 1) {
    	if (pt_GR->SSC == 0 && pt_GR->EC == 5) {

    		if (pt_GR->verbose>1) {
    			printf("CMB\n");
    			printf("nu_start_CMB_seed=%e\n", pt_GR->nu_start_CMB);
    			printf("nu_stop_CMB_seed=%e\n", pt_GR->nu_stop_CMB);
    		}
    		
            if (pt_GR->EC_stat == 0)
            {
                pt_GR->nu_seed = pt_GR->nu_CMB;
                pt_GR->n_seed = pt_GR->n_CMB;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_CMB,
                                         pt_GR->nu_stop_CMB,
                                         pt_GR->EC_stat,
                                         nu_IC_out);
            }
            else
            {
                pt_GR->nu_seed = pt_GR->nu_CMB_disk_RF;
                pt_GR->n_seed = pt_GR->n_CMB_DRF;
                rate_comp = integrale_IC(pt_GR,
                                         pt_GR->nu_start_CMB_DRF,
                                         pt_GR->nu_stop_CMB_DRF,
                                         pt_GR->EC_stat,
                                         nu_IC_out_stat);
            }
        }
    }

    

    return rate_comp;
}
//=========================================================================================


double f_compton_bulk(struct blob *pt_K1, double g, double nu_IC_out, double nu_IC_in_1, double nu_IC_in_2) {
    double cost, rate;
    rate=0;
    if (nu_IC_out >=  nu_IC_in_1 &&  nu_IC_out <nu_IC_in_2) {
        cost = pt_K1->COST_IC_K1/nu_IC_in_1;
       
        rate = cost;
    }

    return rate;
}

//=========================================================================================
// Function to evaluate the kernel of IC emission
// Band & Grindlay pg 138, 1985 ApJ 298
// nu'=nu_IC_in
// nu=nu_IC_out
//=========================================================================================
double f_compton_K1(struct blob *pt_K1, double g, double nu_IC_out, double nu_IC_in) {
    /**
     * \funzione f_compton_K1
     * \author Andrea Tramacere
     * \date 27-04-2004
     * \brief
     *
     * g = Gamma degli e-
     */
    double cost, rate,a, c, k, nu_1_min, nu_1_max, g2;
    double epsilon_0, epsilon_1,Gamma_e;
    //pt_K1->Gamma = g;
    g2 = g*g;
    epsilon_0 = HPLANCK * nu_IC_in*one_by_MEC2;
    epsilon_1 = HPLANCK * nu_IC_out*one_by_MEC2;
    nu_1_min = nu_IC_in/(4.0*g2);
    nu_1_max = 4.0 * nu_IC_in * g2 / (1.0 + 4.0*g*epsilon_0);
  

    //=================================================
    rate=0.0;

    if (nu_IC_out > nu_1_max || nu_IC_out < nu_1_min ) {
       rate=0.0;
    }
    if (nu_IC_out >=  nu_1_min &&  nu_IC_out <nu_IC_in) {

        if (pt_K1->do_IC_down_scattering==1){
        //------------------------------------------
        //Eq 8 Jones 1968 
        //Eq IV.I  Band & Grindlay 1985 ApJ 298
        //This is the down-scattering and is optional
        cost = pt_K1->COST_IC_K1 / (4.0*(g2*g2) * nu_IC_in);
        k=4.0*g2*nu_IC_out/nu_IC_in ;
        rate=k-1;
        rate *= cost;
        }else{
            rate=0;
        }
    }

    if (nu_IC_out >= nu_IC_in && nu_IC_out <= nu_1_max) {
        //-----------------------
        //Eq 44 Jones 1968
        //Eq IV.I  Band & Grindlay 1985 ApJ 298
        k=nu_IC_out / (nu_IC_in * 4.0 *( g2 - epsilon_1*g));
        //this condition is superfluous!!!
        //if (k>1.0/(4*g2) && k<=1){
            //printf("nu_1=%e nu_min=%e nu_max=%e gamma=%e, k2=%e 1/(4g^2)=%e\n",nu_IC_out,nu_1_min,nu_1_max,g,k,(1.0/(4*g2)));

        Gamma_e=4.0*g*epsilon_0;

        cost = pt_K1->COST_IC_K1 / ((g2) *nu_IC_in);

        a = 2.0 * k * log(k) ;

        a = a + (1+2*k)*(1-k);

        c = 0.5*(1-k)*(Gamma_e*k)*(Gamma_e*k)/(1+4.0*k*Gamma_e);

        rate = a+c;
        rate *= cost;
            //printf("2\n");
        //}
        //else{
        // rate=0;
         //printf("3\n");
        //}
    }
    return rate;
}
//=========================================================================================

void set_N_distr_for_Compton(struct blob * pt, double nu_in, double nu_out, int stat_frame, double * Ne_IC, double * griglia_gamma_Ne_log_IC)
{
    double epsilon_0, epsilon_1,g_min_IC;
    epsilon_0 = HPLANCK * nu_in * one_by_MEC2;
    epsilon_1 = HPLANCK * nu_out * one_by_MEC2;
    
    //Eq 7.111 Dermer&Menon 2009
    g_min_IC = 0.5 * epsilon_1 *(1 + sqrt(1.0 + (1.0 / (epsilon_1 * epsilon_0))));
    
    
    if (pt->EC_stat == 1)
    {
        g_min_IC = g_min_IC / pt->beam_obj;
    }
    if (g_min_IC > pt->gmin_griglia)
    {
        Fill_Ne_IC(pt, g_min_IC, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);
    }
    else
    {
        Fill_Ne_IC(pt, pt->gmin_griglia, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);
    }
}

//=========================================================================================
// IC INTEGRATION METHOD TRAPEZOIDAL/SIMPSON_GRID_EQUI_LOG
// a,b: boundaries for photon integration
// returns [emitted photons, cm-3, s-1, Hz-1, sterad-1]
// the [sterad-1] comes from n_seed
//=========================================================================================
double integrale_IC( struct blob * pt, double a, double b, int stat_frame, double nu_IC_out) {
    double integr_nu, nu_IC_in;
    
    unsigned int ID,ID_gamma;
    double *Integrand_over_gamma_grid, *Ne_IC, *griglia_gamma_Ne_log_IC, *integr_gamma;
    Integrand_over_gamma_grid = (double *) calloc(pt->gamma_grid_size, sizeof (double));
    griglia_gamma_Ne_log_IC =  (double *) calloc(pt->gamma_grid_size, sizeof (double));
    integr_gamma = (double *) calloc(pt->nu_seed_size, sizeof (double));
    Ne_IC = (double *) calloc(pt->gamma_grid_size, sizeof (double));
    double ic_kernel;
    integr_nu = 0.0;
   


    set_N_distr_for_Compton(pt, b, nu_IC_out, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);

    if (pt->verbose>1) {
        printf("***** Integrale  IC ******\n");
        printf("i=%d\n", ID);
        printf("nu=%e a=%e b=%e  g_min_grid=%e g_max_grid=%e\n", pt->nu_seed[ID], a, b, griglia_gamma_Ne_log_IC[0], griglia_gamma_Ne_log_IC[pt->gamma_grid_size - 1]);
    }
 
    for (ID=0; ID<pt->nu_seed_size; ID++){
        if (pt->nu_seed[ID] <= b && pt->nu_seed[ID] >= a){
            nu_IC_in= pt->nu_seed[ID];

            //Integration over electron Lorentz factor
            for (ID_gamma = 0; ID_gamma < pt->gamma_grid_size ; ID_gamma++){
                if (pt->bulk_compton == 0){
                    ic_kernel=f_compton_K1(pt, griglia_gamma_Ne_log_IC[ID_gamma], nu_IC_out, nu_IC_in);
                }else{
                    ic_kernel=f_compton_bulk(pt, griglia_gamma_Ne_log_IC[ID_gamma], nu_IC_out,  pt->nu_seed[ID], pt->nu_seed[ID+1]);
                }
                Integrand_over_gamma_grid[ID_gamma] =ic_kernel * Ne_IC[ID_gamma];
                
            }
            integr_gamma[ID]= pt->n_seed[ID]*integr_simp_grid_equilog(griglia_gamma_Ne_log_IC, Integrand_over_gamma_grid, pt->gamma_grid_size);

        }else{
            integr_gamma[ID]=0;
        }
    }
    integr_nu=trapzd_array_arbritary_grid( pt->nu_seed,integr_gamma, pt->nu_seed_size);

    //============================================================
    //0.75 fattore di correzione di GOULD
    //has been moved to spetto_sincrotrone.c
    //============================================================
    free(Integrand_over_gamma_grid);
    free(Ne_IC);
    free(griglia_gamma_Ne_log_IC);
    free(integr_gamma);
    return integr_nu;
}
//=========================================================================================





//=========================================================================================
// Cooling Compton
//=========================================================================================
double compton_cooling(struct blob *pt_spec, struct temp_ev *pt_ev, double gamma) {
    /**
     * \author Andrea Tramacere
     * \date 27-04-2004
     * \brief funzioni per il Compton
     *  CALCOLO DEL COMPTON Cooling sezione d'urto da Moderski et al. 2005 MNRAS 363
     */
    double comp_cooling;

    comp_cooling=0;
    
    if (pt_spec->verbose>1) {
        printf("GR\n");
        printf("#-> SSC=%d EC=%d\n", pt_spec->SSC, pt_spec->EC);
        printf("#-> gmin=%e gmax=%e\n", pt_spec->gmin, pt_spec->gmax);
    }

    //SSC
    if (pt_spec->do_Sync) {
        if (pt_spec->verbose>1) {
            printf("nu_start_Sync=%e\n", pt_spec->nu_start_Sync);
            printf("nu_stop_Sync_ssc=%e\n", pt_spec->nu_stop_Sync_ssc);

        }
         
        pt_spec->nu_seed = pt_spec->nu_Sync;
        pt_spec->n_seed = pt_spec->n_Sync;
        comp_cooling += integrale_IC_cooling(pt_spec,
                pt_spec->nu_start_Sync,
                pt_spec->nu_stop_Sync_ssc,
                gamma);
        //printf("evaluate IC cooling, gamma=%e cooling_rate=%e, Sync_cooling_rate_ratio=%e\n",gamma,comp_cooling,comp_cooling/Sync_cool(pt_spec->B,gamma));
    }

    //EC Disk

    if (pt_spec->do_EC_Disk == 1 ) {

        if (pt_spec->verbose>1) {
            printf("Disk\n");
            printf("nu_start_EC_seed=%e\n", pt_spec->nu_start_Disk);
            printf("nu_stop_EC_seed=%e\n", pt_spec->nu_stop_Disk);
        }
        pt_spec->nu_seed = pt_spec->nu_Disk;
        pt_spec->n_seed = pt_spec->n_Disk;
        comp_cooling += integrale_IC_cooling(pt_spec,
                pt_spec->nu_start_Disk,
                pt_spec->nu_stop_Disk,
                gamma);
        //printf("%e\n",rate_comp);
    }

    //EC BLR
    if (pt_spec->do_EC_BLR == 1 ) {

    	if (pt_spec->verbose>1) {
    		printf("BLR\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->nu_start_BLR);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->nu_stop_BLR);
    	}
    	pt_spec->nu_seed = pt_spec->nu_BLR;
    	pt_spec->n_seed = pt_spec->n_BLR;
    	comp_cooling += integrale_IC_cooling(pt_spec,
    			pt_spec->nu_start_BLR,
    			pt_spec->nu_stop_BLR,
    			gamma);
    	//printf("%e\n",rate_comp);
    }


    //EC DT
    if (pt_spec->do_EC_DT == 1 ) {

    	if (pt_spec->verbose>1) {
    		printf("DT\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->nu_start_DT);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->nu_stop_DT);
    	}
    	pt_spec->nu_seed = pt_spec->nu_DT;
    	pt_spec->n_seed = pt_spec->n_DT;
    	comp_cooling += integrale_IC_cooling(pt_spec,
    			pt_spec->nu_start_DT,
    			pt_spec->nu_stop_DT,
    			gamma);
    	//printf("%e\n",rate_comp);
    }

    //EC Star
    if (pt_spec->do_EC_Star == 1 ) {

    	if (pt_spec->verbose>1) {
    		printf("Star\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->nu_start_Star);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->nu_stop_Star);
    	}
    	pt_spec->nu_seed = pt_spec->nu_Star;
    	pt_spec->n_seed = pt_spec->n_Star;
    	comp_cooling += integrale_IC_cooling(pt_spec,
    			pt_spec->nu_start_Star,
    			pt_spec->nu_stop_Star,
    			gamma);
    	//printf("%e\n",rate_comp);
    }

    //EC CMB
    if (pt_spec->do_EC_CMB == 1 ) {

    	if (pt_spec->verbose>1) {
    		printf("CMB\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->nu_start_CMB);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->nu_stop_CMB);
    	}
    	pt_spec->nu_seed = pt_spec->nu_CMB;
    	pt_spec->n_seed = pt_spec->n_CMB;
    	comp_cooling += integrale_IC_cooling(pt_spec,
    			pt_spec->nu_start_CMB,
    			pt_spec->nu_stop_CMB,
    			gamma);
    	//printf("%e\n",rate_comp);
    }

    //printf("evaluate IC cooling, gamma=%e cooling_rate=%e\n",gamma,comp_cooling);
    return comp_cooling;
}
//=========================================================================================





//=========================================================================================
// INTEGRAZIONE DEL COMPTON COOLING CON METODO TRAPEZIO
//=========================================================================================
double integrale_IC_cooling(struct blob * pt, double a, double b, double gamma) {
    double nu1, nu2, integr_nu;
    double y_nu1, y_nu2;
    double delta_nu,b_kn;
    unsigned int i;

    i = 0;
    integr_nu=0;
    while (pt->nu_seed[i] < a) {
        //  printf("i=%d\n",i);
        i++;
    }

    if (pt->verbose>1) {
        printf("***** Integrale IC cooling ******\n");
        printf("i=%d\n", i);
        printf("nu=%e a=%e i=%d\n", pt->nu_seed[i], a, i);
    }

    nu1 = pt->nu_seed[i];
    b_kn=4*gamma*pt->nu_seed[i]*HPLANCK*one_by_MEC2;
    y_nu1 = pt->n_seed[i] * f_compton_cooling(b_kn)*nu1;
    

    while (pt->nu_seed[i + 1] <= b && pt->nu_seed[i + 1] >= a) {
       

        b_kn=4*gamma*pt->nu_seed[i+1]*HPLANCK*one_by_MEC2;
        //printf("b=%e nu=%e f_kn=%e\n",b_kn,pt->nu_seed[i+1],f_compton_cooling(b_kn));

        
        nu2=pt->nu_seed[i+1];
        y_nu2 = pt->n_seed[i + 1] * f_compton_cooling(b_kn)*nu2;


        delta_nu = nu2 - nu1;

        integr_nu += (y_nu2 + y_nu1) * delta_nu;
        nu1 = nu2;
        y_nu1 = y_nu2;
        i++;
    }
    integr_nu *= gamma * gamma * pt->COST_IC_COOLING;
    //printf("integr_nu=%e\n",integr_nu);
    //============================================================
    //lo 0.5 viene dalla regola del trapezio dell'integrale in nu
    //4PI comes from Uph  integrated over angles
    //============================================================
    return (integr_nu * 0.5)*4*pi;

}
//=========================================================================================






//=========================================================================================
// Kernel per il  Compton coolig, Moderski et al. 2005 MNRAS 363
// Eq. 3
// in my code b=4*gamma*pt->nu_seed[i+1]*HPLANCK*one_by_MEC2
// in the paper b=4*gamma*h*nu/mec^2
// I use the approximation in Eq. 3
// I use b<1000 that gives a better connection compared to the value of b<10000
// used in the paper
// f_KN=1/(1+b)^1.5 if b<1000
// else 9/(b^2)*(log(b)-11/6)
// I_nu_Sync=>I_nu_seed
//=========================================================================================
double f_compton_cooling(double b) {
    if (b < 1000.0) {
        return  1.0 / pow((1 + b), 1.5);
    }
    else {
        return 9.0 / (2.0 * b * b)*(log(b) - 11. / 6.);
    }
}
//=========================================================================================
