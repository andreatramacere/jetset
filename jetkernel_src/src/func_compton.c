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

    double * nu_seed;
    double * n_seed;
    unsigned int nu_seed_size;

    nu_seed_size=pt_GR->core.nu_seed_size;
    nu_IC_out_stat = nu_IC_out * pt_GR->core.beam_obj;

    if (pt_GR->core.verbose>1) {
        printf("GR\n");
        printf("#-> SSC=%d EC=%d\n", pt_GR->core.SSC, pt_GR->core.EC);
        printf("#-> gmin=%e gmax=%e\n", pt_GR->emitters.gmin, pt_GR->emitters.gmax);
    }

    //SSC
    if (nu_IC_out < pt_GR->SSC.spec.nu_max && pt_GR->core.ord_comp == 1) {
        if (pt_GR->core.SSC == 1 && pt_GR->core.EC == 0) {
            if (pt_GR->core.verbose>1) {
                printf("nu_start_Sync=%e\n", pt_GR->Sync.spec.nu_min);
                printf("nu_stop_Sync_ssc=%e\n", pt_GR->Sync.nu_stop_Sync_ssc);
            }
            nu_seed = pt_GR->Sync.spec.nu;
            n_seed = pt_GR->Sync.spec.n_nu;
            //pt_GR->griglia_gamma_log_IC=pt_GR->emitters.griglia_gamma_Ne_log;
            //pt_GR->N_IC=pt_GR->emitters.Ne;
            rate_comp = integrale_IC(pt_GR,
                    nu_seed,
                    n_seed,
                    nu_seed_size,
                    pt_GR->Sync.spec.nu_min,
                    pt_GR->Sync.nu_stop_Sync_ssc,
                    0,
                    nu_IC_out);
        }
    }
    //EC Disk
    if (nu_IC_out < pt_GR->Disk.ec.spec.nu_max && pt_GR->core.ord_comp == 1) {
		if (pt_GR->core.SSC == 0 && pt_GR->core.EC == 1) {

			if (pt_GR->core.verbose>1) {
				printf("Disk\n");
				printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->Disk.spec.nu_min);
                printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->Disk.spec.nu_max);
            }
            if (pt_GR->core.EC_stat == 0)
            {
                nu_seed = pt_GR->Disk.spec.nu;
                n_seed = pt_GR->Disk.spec.n_nu;
                rate_comp = integrale_IC(pt_GR,
                                        nu_seed,
                                        n_seed,
                                        nu_seed_size,
                                        pt_GR->Disk.spec.nu_min,
                                        pt_GR->Disk.spec.nu_max,
                                        pt_GR->core.EC_stat,
                                        nu_IC_out);
            }
            else{
                nu_seed = pt_GR->Disk.spec.nu_DRF;
                n_seed = pt_GR->Disk.spec.n_nu_DRF;
                rate_comp = integrale_IC(pt_GR,
                                        nu_seed,
                                        n_seed,
                                        nu_seed_size,
                                        pt_GR->Disk.spec.nu_min_DRF,
                                        pt_GR->Disk.spec.nu_max_DRF,
                                        pt_GR->core.EC_stat,
                                        nu_IC_out_stat);
            }
			
		}
    }
    //EC BLR
    if (nu_IC_out < pt_GR->BLR.ec.spec.nu_max && pt_GR->core.ord_comp == 1) {
    	if (pt_GR->core.SSC == 0 && pt_GR->core.EC == 2) {

            if (pt_GR->core.verbose>1) {
                printf("BLR\n");
                printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->BLR.spec.nu_min);
                printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->BLR.spec.nu_max);
            }

            if (pt_GR->core.EC_stat == 0)
            {
                nu_seed = pt_GR->BLR.spec.nu;
                n_seed = pt_GR->BLR.spec.n_nu;
                rate_comp = integrale_IC(pt_GR,
                                        nu_seed,
                                        n_seed,
                                        nu_seed_size,
                                        pt_GR->BLR.spec.nu_min,
                                        pt_GR->BLR.spec.nu_max,
                                        pt_GR->core.EC_stat,
                                        nu_IC_out);
            }
            else
            {             
                nu_seed = pt_GR->BLR.spec.nu_DRF;
                n_seed = pt_GR->BLR.spec.n_nu_DRF;
                rate_comp = integrale_IC(pt_GR,
                                         nu_seed,
                                         n_seed,
                                         nu_seed_size,
                                         pt_GR->BLR.spec.nu_min_DRF,
                                         pt_GR->BLR.spec.nu_max_DRF,
                                         pt_GR->core.EC_stat,
                                         nu_IC_out_stat);

            }
        }
    }
    //EC DT
    if (nu_IC_out < pt_GR->DT.ec.spec.nu_max && pt_GR->core.ord_comp == 1) {
    	if (pt_GR->core.SSC == 0 && pt_GR->core.EC == 3) {
            if (pt_GR->core.verbose>1) {
                printf("DT\n");
                printf("(blob rest frame) nu_start_EC_seed DT=%e\n", pt_GR->DT.spec.nu_min);
                printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->DT.spec.nu_max);
            }

            if (pt_GR->core.EC_stat == 0){
               
                nu_seed = pt_GR->DT.spec.nu;
                n_seed = pt_GR->DT.spec.n_nu;
                rate_comp = integrale_IC(pt_GR,
                                         nu_seed,
                                         n_seed,
                                         nu_seed_size,
                                         pt_GR->DT.spec.nu_min,
                                         pt_GR->DT.spec.nu_max,
                                         pt_GR->core.EC_stat,
                                         nu_IC_out);
            }
            else
            {
                nu_seed = pt_GR->DT.spec.nu_DRF;
                n_seed = pt_GR->DT.spec.n_nu_DRF;
                rate_comp = integrale_IC(pt_GR,
                                         nu_seed,
                                         n_seed,
                                         nu_seed_size,
                                         pt_GR->DT.spec.nu_min,
                                         pt_GR->DT.spec.nu_max_DRF,
                                         pt_GR->core.EC_stat,
                                         nu_IC_out_stat);
            }
       }
    }

    //EC Star
    if (nu_IC_out < pt_GR->Star.ec.spec.nu_max && pt_GR->core.ord_comp == 1) {
    	if (pt_GR->core.SSC == 0 && pt_GR->core.EC == 4) {

		   if (pt_GR->core.verbose>1) {
			   printf("DT\n");
               printf("(blob rest frame) nu_start_EC_seed Star=%e\n", pt_GR->Star.spec.nu_min);
               printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->Star.spec.nu_max);
           }
		   nu_seed = pt_GR->Star.spec.nu;
		   n_seed = pt_GR->Star.spec.n_nu;
           if (pt_GR->core.EC_stat == 0)
           {
               nu_seed = pt_GR->Star.spec.nu;
               n_seed = pt_GR->Star.spec.n_nu;
               rate_comp = integrale_IC(pt_GR,
                                        nu_seed,
                                        n_seed,
                                        nu_seed_size,
                                        pt_GR->Star.spec.nu_min,
                                        pt_GR->Star.spec.nu_max,
                                        pt_GR->core.EC_stat,
                                        nu_IC_out);
           }
           else
           {
               nu_seed = pt_GR->Star.spec.nu_DRF;
               n_seed = pt_GR->Star.spec.n_nu_DRF;
               rate_comp = integrale_IC(pt_GR,
                                        nu_seed,
                                        n_seed,
                                        nu_seed_size,
                                        pt_GR->Star.spec.nu_min_DRF,
                                        pt_GR->Star.spec.nu_max_DRF,
                                        pt_GR->core.EC_stat,
                                        nu_IC_out_stat);
           }
       }
    }

    //EC CMB
    if (nu_IC_out < pt_GR->CMB.ec.spec.nu_max && pt_GR->core.ord_comp == 1) {
    	if (pt_GR->core.SSC == 0 && pt_GR->core.EC == 5) {

    		if (pt_GR->core.verbose>1) {
    			printf("CMB\n");
    			printf("nu_start_CMB_seed=%e\n", pt_GR->CMB.spec.nu_min);
    			printf("nu_stop_CMB_seed=%e\n", pt_GR->CMB.spec.nu_max);
    		}
    		
            if (pt_GR->core.EC_stat == 0)
            {
                nu_seed = pt_GR->CMB.spec.nu;
                n_seed = pt_GR->CMB.spec.n_nu;
                rate_comp = integrale_IC(pt_GR,
                                         nu_seed,
                                         n_seed,
                                         nu_seed_size,
                                         pt_GR->CMB.spec.nu_min,
                                         pt_GR->CMB.spec.nu_max,
                                         pt_GR->core.EC_stat,
                                         nu_IC_out);
            }
            else
            {
                nu_seed = pt_GR->CMB.spec.nu_DRF;
                n_seed = pt_GR->CMB.spec.n_nu_DRF;
                rate_comp = integrale_IC(pt_GR,
                                         nu_seed,
                                         n_seed,
                                         nu_seed_size,
                                         pt_GR->CMB.spec.nu_min_DRF,
                                         pt_GR->CMB.spec.nu_max_DRF,
                                         pt_GR->core.EC_stat,
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
        cost = pt_K1->core.COST_IC_K1/nu_IC_in_1;
       
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

        if (pt_K1->core.do_IC_down_scattering==1){
        //------------------------------------------
        //Eq 8 Jones 1968 
        //Eq IV.I  Band & Grindlay 1985 ApJ 298
        //This is the down-scattering and is optional
        cost = pt_K1->core.COST_IC_K1 / (4.0*(g2*g2) * nu_IC_in);
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

        cost = pt_K1->core.COST_IC_K1 / ((g2) *nu_IC_in);

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
    
    
    if (pt->core.EC_stat == 1)
    {
        g_min_IC = g_min_IC / pt->core.beam_obj;
    }
    if (g_min_IC > pt->emitters.gmin_griglia)
    {
        Fill_Ne_IC(pt, g_min_IC, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);
    }
    else
    {
        Fill_Ne_IC(pt, pt->emitters.gmin_griglia, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);
    }
}

//=========================================================================================
// IC INTEGRATION METHOD TRAPEZOIDAL/SIMPSON_GRID_EQUI_LOG
// a,b: boundaries for photon integration
// returns [emitted photons, cm-3, s-1, Hz-1, sterad-1]
// the [sterad-1] comes from n_seed
//=========================================================================================
//double integrale_IC( struct blob * pt, double a, double b, int stat_frame, double nu_IC_out) 
double integrale_IC(struct blob *pt, const double *nu_seed, const double *n_seed, unsigned int nu_seed_size, double a, double b, int stat_frame, double nu_IC_out){
    double integr_nu, nu_IC_in;
    
    unsigned int ID,ID_gamma;
    double *Integrand_over_gamma_grid, *Ne_IC, *griglia_gamma_Ne_log_IC, *integr_gamma;
    Integrand_over_gamma_grid = (double *) calloc(pt->emitters.gamma_grid_size, sizeof (double));
    griglia_gamma_Ne_log_IC =  (double *) calloc(pt->emitters.gamma_grid_size, sizeof (double));
    integr_gamma = (double *) calloc(nu_seed_size, sizeof (double));
    Ne_IC = (double *) calloc(pt->emitters.gamma_grid_size, sizeof (double));
    double ic_kernel;
    integr_nu = 0.0;
   


    set_N_distr_for_Compton(pt, b, nu_IC_out, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);

    for (ID=0; ID<nu_seed_size; ID++){
        if (nu_seed[ID] <= b && nu_seed[ID] >= a){
            nu_IC_in= nu_seed[ID];

            //Integration over electron Lorentz factor
            for (ID_gamma = 0; ID_gamma < pt->emitters.gamma_grid_size ; ID_gamma++){
                if (pt->core.bulk_compton == 0){
                    ic_kernel=f_compton_K1(pt, griglia_gamma_Ne_log_IC[ID_gamma], nu_IC_out, nu_IC_in);
                }else{
                    if (ID<nu_seed_size-1){
                        ic_kernel=f_compton_bulk(pt, griglia_gamma_Ne_log_IC[ID_gamma], nu_IC_out,   nu_seed[ID],  nu_seed[ID+1]);
                    }else{
                        ic_kernel=f_compton_bulk(pt, griglia_gamma_Ne_log_IC[ID_gamma], nu_IC_out,   nu_seed[ID-1],  nu_seed[ID]);
                    }
                }    
                Integrand_over_gamma_grid[ID_gamma] =ic_kernel * Ne_IC[ID_gamma];
                
            }
            integr_gamma[ID]= n_seed[ID]*integr_simp_grid_equilog(griglia_gamma_Ne_log_IC, Integrand_over_gamma_grid, pt->emitters.gamma_grid_size);

        }else{
            integr_gamma[ID]=0;
        }
    }
    integr_nu=trapzd_array_arbritary_grid( nu_seed,integr_gamma, nu_seed_size);

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
    double * nu_seed;
    double * n_seed;

    unsigned int nu_seed_size;

    nu_seed_size=pt_spec->core.nu_seed_size;
    
    if (pt_spec->core.verbose>1) {
        printf("GR\n");
        printf("#-> SSC=%d EC=%d\n", pt_spec->core.SSC, pt_spec->core.EC);
        printf("#-> gmin=%e gmax=%e\n", pt_spec->emitters.gmin, pt_spec->emitters.gmax);
    }

    //SSC
    if (pt_spec->core.do_Sync) {
        if (pt_spec->core.verbose>1) {
            printf("nu_start_Sync=%e\n", pt_spec->Sync.spec.nu_min);
            printf("nu_stop_Sync_ssc=%e\n", pt_spec->Sync.nu_stop_Sync_ssc);

        }
         
        nu_seed = pt_spec->Sync.spec.nu;
        n_seed = pt_spec->Sync.spec.n_nu;
        comp_cooling += integrale_IC_cooling(pt_spec,
                                             nu_seed,
                                             n_seed,
                                             nu_seed_size,
                                             pt_spec->Sync.spec.nu_min,
                                             pt_spec->Sync.nu_stop_Sync_ssc,
                                             gamma);
        //printf("evaluate IC cooling, gamma=%e cooling_rate=%e, Sync_cooling_rate_ratio=%e\n",gamma,comp_cooling,comp_cooling/Sync_cool(pt_spec->core.B,gamma));
    }

    //EC Disk

    if (pt_spec->core.do_EC_Disk == 1 ) {

        if (pt_spec->core.verbose>1) {
            printf("Disk\n");
            printf("nu_start_EC_seed=%e\n", pt_spec->Disk.spec.nu_min);
            printf("nu_stop_EC_seed=%e\n", pt_spec->Disk.spec.nu_max);
        }
        nu_seed = pt_spec->Disk.spec.nu;
        n_seed = pt_spec->Disk.spec.n_nu;
        comp_cooling += integrale_IC_cooling(pt_spec,
                nu_seed,
                n_seed,
                nu_seed_size,
                pt_spec->Disk.spec.nu_min,
                pt_spec->Disk.spec.nu_max,
                gamma);
        //printf("%e\n",rate_comp);
    }

    //EC BLR
    if (pt_spec->core.do_EC_BLR == 1 ) {

    	if (pt_spec->core.verbose>1) {
    		printf("BLR\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->BLR.spec.nu_min);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->BLR.spec.nu_max);
    	}
    	nu_seed = pt_spec->BLR.spec.nu;
    	n_seed = pt_spec->BLR.spec.n_nu;
    	comp_cooling += integrale_IC_cooling(pt_spec,
                nu_seed,
                n_seed,
                nu_seed_size,
    			pt_spec->BLR.spec.nu_min,
    			pt_spec->BLR.spec.nu_max,
    			gamma);
    	//printf("%e\n",rate_comp);
    }


    //EC DT
    if (pt_spec->core.do_EC_DT == 1 ) {

    	if (pt_spec->core.verbose>1) {
    		printf("DT\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->DT.spec.nu_min);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->DT.spec.nu_max);
    	}
    	nu_seed = pt_spec->DT.spec.nu;
    	n_seed = pt_spec->DT.spec.n_nu;
    	comp_cooling += integrale_IC_cooling(pt_spec,
                nu_seed,
                n_seed,
                nu_seed_size,
    			pt_spec->DT.spec.nu_min,
    			pt_spec->DT.spec.nu_max,
    			gamma);
    	//printf("%e\n",rate_comp);
    }

    //EC Star
    if (pt_spec->core.do_EC_Star == 1 ) {

    	if (pt_spec->core.verbose>1) {
    		printf("Star\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->Star.spec.nu_min);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->Star.spec.nu_max);
    	}
    	nu_seed = pt_spec->Star.spec.nu;
    	n_seed = pt_spec->Star.spec.n_nu;
    	comp_cooling += integrale_IC_cooling(pt_spec,
                nu_seed,
                n_seed,
                nu_seed_size,
    			pt_spec->Star.spec.nu_min,
    			pt_spec->Star.spec.nu_max,
    			gamma);
    	//printf("%e\n",rate_comp);
    }

    //EC CMB
    if (pt_spec->core.do_EC_CMB == 1 ) {

    	if (pt_spec->core.verbose>1) {
    		printf("CMB\n");
    		printf("nu_start_EC_seed=%e\n", pt_spec->CMB.spec.nu_min);
    		printf("nu_stop_EC_seed=%e\n", pt_spec->CMB.spec.nu_max);
    	}
    	nu_seed = pt_spec->CMB.spec.nu;
    	n_seed = pt_spec->CMB.spec.n_nu;
    	comp_cooling += integrale_IC_cooling(pt_spec,
                nu_seed,
                n_seed,
                nu_seed_size,
    			pt_spec->CMB.spec.nu_min,
    			pt_spec->CMB.spec.nu_max,
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
//double integrale_IC_cooling(struct blob * pt, double a, double b, double gamma) 
double integrale_IC_cooling(struct blob *pt, const double *nu_seed, const double *n_seed, unsigned int nu_seed_size, double a, double b, double gamma) {
    double nu1, nu2, integr_nu;
    double y_nu1, y_nu2;
    double delta_nu,b_kn;
    unsigned int i;

    i = 0;
    integr_nu=0;
    while (i<nu_seed_size-1 && nu_seed[i] < a) {
        //  printf("i=%d\n",i);
        i++;
    }

    //if (pt->core.verbose>1) {
    //    printf("***** Integrale IC cooling ******\n");
    //    printf("i=%d\n", i);
    //    printf("nu=%e a=%e i=%d\n", nu_seed[i], a, i);
    //}

    nu1 = nu_seed[i];
    b_kn=4*gamma*nu_seed[i]*HPLANCK*one_by_MEC2;
    y_nu1 = n_seed[i] * f_compton_cooling(b_kn)*nu1;
    

    while ( i<nu_seed_size-1 && nu_seed[i + 1] <= b && nu_seed[i + 1] >= a) {
       

        b_kn=4*gamma*nu_seed[i+1]*HPLANCK*one_by_MEC2;
        //printf("b=%e nu=%e f_kn=%e\n",b_kn,pt->nu_seed[i+1],f_compton_cooling(b_kn));

        
        nu2=nu_seed[i+1];
        y_nu2 = n_seed[i + 1] * f_compton_cooling(b_kn)*nu2;


        delta_nu = nu2 - nu1;

        integr_nu += (y_nu2 + y_nu1) * delta_nu;
        nu1 = nu2;
        y_nu1 = y_nu2;
        i++;
    }
    integr_nu *= gamma * gamma * pt->core.COST_IC_COOLING;
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
