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
double rate_compton_GR(struct spettro *pt_GR) {
    /**
     * \author Andrea Tramacere
     * \date 27-04-2004
     * \brief funzioni per il Compton
     *  CALCOLO DEL RATE COMPTON GENERALE metodo GRINDLAY 1985
     */
    double rate_comp=0;
    double (*pf_K) (struct spettro *, double x);
    double nu_1_original;
    int i;
    pf_K = &f_compton_K1;



    if (pt_GR->verbose>1) {
        printf("GR\n");
        printf("#-> SSC=%d EC=%d\n", pt_GR->SSC, pt_GR->EC);
        printf("#-> gmin=%e gmax=%e\n", pt_GR->gmin, pt_GR->gmax);
    }

    //SSC
    if (pt_GR->nu_1 < pt_GR->nu_stop_SSC && pt_GR->ord_comp == 1) {
        if (pt_GR->SSC == 1 && pt_GR->EC == 0) {
            if (pt_GR->verbose>1) {
                printf("nu_start_Sync=%e\n", pt_GR->nu_start_Sync);
                printf("nu_stop_Sync_ssc=%e\n", pt_GR->nu_stop_Sync_ssc);
            }
            pt_GR->nu_seed = pt_GR->nu_Sync;
            pt_GR->n_seed = pt_GR->n_Sync;
            //pt_GR->griglia_gamma_log_IC=pt_GR->griglia_gamma_Ne_log;
            //pt_GR->N_IC=pt_GR->Ne;
            rate_comp = integrale_IC(pf_K,
                    pt_GR,
                    pt_GR->nu_start_Sync,
                    pt_GR->nu_stop_Sync_ssc,
                    0);
        }
    }
    //EC Diks
    if (pt_GR->nu_1 < pt_GR->nu_stop_EC_Disk && pt_GR->ord_comp == 1) {
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
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_start_Disk,
                                         pt_GR->nu_stop_Disk,
                                         pt_GR->EC_stat);
            }
            else{
                nu_1_original = pt_GR->nu_1;
                pt_GR->nu_1 = pt_GR->nu_1 * pt_GR->beam_obj;

                pt_GR->nu_seed = pt_GR->nu_Disk_disk_RF;
                pt_GR->n_seed = pt_GR->n_Disk_DRF;
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_start_Disk_DRF,
                                         pt_GR->nu_stop_Disk_DRF,
                                         pt_GR->EC_stat);

                pt_GR->nu_1 = nu_1_original;
            }
			
		}
    }
    //EC BLR
    if (pt_GR->nu_1 < pt_GR->nu_stop_EC_BLR && pt_GR->ord_comp == 1) {
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
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_start_BLR,
                                         pt_GR->nu_stop_BLR,
                                         pt_GR->EC_stat);
            }
            else
            {
                nu_1_original = pt_GR->nu_1;
                pt_GR->nu_1 = pt_GR->nu_1 * pt_GR->beam_obj;
                
                pt_GR->nu_seed = pt_GR->nu_BLR_disk_RF;
                pt_GR->n_seed = pt_GR->n_BLR_DRF;
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_start_BLR_disk_RF,
                                         pt_GR->nu_stop_BLR_disk_RF,
                                         pt_GR->EC_stat);
                pt_GR->nu_1 = nu_1_original;
            }
        }
    }
    //EC DT
    if (pt_GR->nu_1 < pt_GR->nu_stop_EC_DT && pt_GR->ord_comp == 1) {
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
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_start_DT,
                                         pt_GR->nu_stop_DT,
                                         pt_GR->EC_stat);
            }
            else
            {
                nu_1_original = pt_GR->nu_1;
                pt_GR->nu_1 = pt_GR->nu_1 * pt_GR->beam_obj;

                pt_GR->nu_seed = pt_GR->nu_DT_disk_RF;
                pt_GR->n_seed = pt_GR->n_DT_DRF;
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_stop_DT_DRF,
                                         pt_GR->nu_stop_DT_DRF,
                                         pt_GR->EC_stat);
                pt_GR->nu_1 = nu_1_original;
            }
       }
    }

    //EC Star
    if (pt_GR->nu_1 < pt_GR->nu_stop_EC_Star && pt_GR->ord_comp == 1) {
    	if (pt_GR->SSC == 0 && pt_GR->EC == 4) {

		   if (pt_GR->verbose>1) {
			   printf("DT\n");
               printf("(blob rest frame) nu_start_EC_seed Star=%e\n", pt_GR->nu_start_Star);
               printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->nu_stop_Star);
           }
		   pt_GR->nu_seed = pt_GR->nu_Star;
		   pt_GR->n_seed = pt_GR->n_Star;
           //pt_GR->griglia_gamma_log_IC=pt_GR->griglia_gamma_Ne_log;
           //pt_GR->N_IC=pt_GR->Ne;
		   //rate_comp = integrale_IC(pf_K,
		   //	   pt_GR,
		   //		   pt_GR->nu_start_Star,
		   //		   pt_GR->nu_stop_Star,
		   //		   0);
		   //printf("%e\n",rate_comp);
           if (pt_GR->EC_stat == 0)
           {
               pt_GR->nu_seed = pt_GR->nu_Star;
               pt_GR->n_seed = pt_GR->n_Star;
               rate_comp = integrale_IC(pf_K,
                                        pt_GR,
                                        pt_GR->nu_start_Star,
                                        pt_GR->nu_stop_Star,
                                        pt_GR->EC_stat);
           }
           else
           {
               nu_1_original = pt_GR->nu_1;
               pt_GR->nu_1 = pt_GR->nu_1 * pt_GR->beam_obj;

               pt_GR->nu_seed = pt_GR->nu_Star_disk_RF;
               pt_GR->n_seed = pt_GR->n_Star_DRF;
               rate_comp = integrale_IC(pf_K,
                                        pt_GR,
                                        pt_GR->nu_start_Star_DRF,
                                        pt_GR->nu_stop_Star_DRF,
                                        pt_GR->EC_stat);
               pt_GR->nu_1 = nu_1_original;
           }
       }
    }

    //EC CMB
    if (pt_GR->nu_1 < pt_GR->nu_stop_EC_CMB && pt_GR->ord_comp == 1) {
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
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_start_CMB,
                                         pt_GR->nu_stop_CMB,
                                         pt_GR->EC_stat);
            }
            else
            {
                nu_1_original = pt_GR->nu_1;
                pt_GR->nu_1 = pt_GR->nu_1 * pt_GR->beam_obj;
                
                pt_GR->nu_seed = pt_GR->nu_CMB_disk_RF;
                pt_GR->n_seed = pt_GR->n_CMB_DRF;
                rate_comp = integrale_IC(pf_K,
                                         pt_GR,
                                         pt_GR->nu_start_CMB_DRF,
                                         pt_GR->nu_stop_CMB_DRF,
                                         pt_GR->EC_stat);
                pt_GR->nu_1 = nu_1_original;
            }
        }
    }

    

 /*    //EC CMB stat
    if (pt_GR->nu_1 < pt_GR->nu_stop_EC_CMB && pt_GR->ord_comp == 1) {
        if (pt_GR->SSC == 0 && pt_GR->EC == 6) {

            if (pt_GR->verbose>1) {
                printf("CMB\n");
                printf("nu_start_CMB_seed=%e\n", pt_GR->nu_start_CMB_stat);
                printf("nu_stop_CMB_seed=%e\n", pt_GR->nu_stop_CMB_stat);
            }

            //sed photon field in the disk rest frame
            pt_GR->nu_seed = pt_GR->nu_CMB_stat;
            pt_GR->n_seed = pt_GR->n_CMB_stat;

            // nu out comp in disk rest frame
            // nu seed is alread in disk rest frame 

            nu_1_original=pt_GR->nu_1;
            pt_GR->nu_1=pt_GR->nu_1*pt_GR->beam_obj;


            //pt_GR->griglia_gamma_log_IC=pt_GR->griglia_gamma_Ne_log_stat;
            //pt_GR->N_IC=pt_GR->Ne_stat;
            rate_comp = integrale_IC(pf_K,
                    pt_GR,
                    pt_GR->nu_start_CMB_stat,
                    pt_GR->nu_stop_CMB_stat,
                    1);
            
           

            pt_GR->nu_1=nu_1_original;



            //printf("rate comp%e\n",rate_comp);
        }

    } */


    return rate_comp;
}
//=========================================================================================






//=========================================================================================
// Function to evaluate the kernel of IC emission
// Band & Grindlay pg 138, 1985 ApJ 298
// nu'=nu_compton_0
// nu=nu_1
//=========================================================================================
double f_compton_K1(struct spettro *pt_K1, double g) {
    /**
     * \funzione f_compton_K1
     * \author Andrea Tramacere
     * \date 27-04-2004
     * \brief
     *
     * g = Gamma degli e-
     */
    double cost, rate,a, c, k, nu_1_min, nu_1_max, Gamma, g2;
    double epsilon_0, epsilon_1,Gamma_e;
    pt_K1->Gamma = g;
    g2 = g*g;
    epsilon_0 = HPLANCK * pt_K1->nu_compton_0*one_by_MEC2;
    epsilon_1 = HPLANCK * pt_K1->nu_1*one_by_MEC2;
    nu_1_min = pt_K1->nu_compton_0/(4.0*g2);
    nu_1_max = 4.0 * pt_K1->nu_compton_0 * (g2) / (1.0 + 4.0*g*epsilon_0);


    //=================================================
    //Se la frequenza del fotone prodotto per IC
    //fuori dal range derivato dalla cinamatica
    //allora il contributo 0 ed esco

    rate=0.0;

    if (pt_K1->nu_1 > nu_1_max || pt_K1->nu_1 < nu_1_min ) {
       rate=0.0;
    }
    if (pt_K1->nu_1 >=  nu_1_min &&  pt_K1->nu_1 < pt_K1->nu_compton_0) {

        //------------------------------------------
        //Eq 8 Jones 1968
        //This should become optional
        //cost = pt_K1->COST_IC_K1 / (4.0*(g2*g2) * pt_K1->nu_compton_0);
        //k=4.0*g2*pt_K1->nu_1/pt_K1->nu_compton_0 ;
        //rate=k-1;
        //rate *= cost;

        //This I don't remember where it is coming from
        //rate = (k-1)*(1.0+2.0/q)-2.0*log(k);


        rate=0.;
        //printf("1\n");
        //printf("nu_1=%e nu_min=%e nu_max=%e gamma=%e, k1=%e\n",pt_K1->nu_1,nu_1_min,nu_1_max,g,k);
    }

    if (pt_K1->nu_1 >= pt_K1->nu_compton_0 && pt_K1->nu_1 <= nu_1_max) {


        k = pt_K1->nu_1 / (pt_K1->nu_compton_0 * 4.0 * ( g2 - g*epsilon_1));
        if (k>1.0/(4*g2) && k<=1){
            //printf("nu_1=%e nu_min=%e nu_max=%e gamma=%e, k2=%e 1/(4g^2)=%e\n",pt_K1->nu_1,nu_1_min,nu_1_max,g,k,(1.0/(4*g2)));

            Gamma_e=4.0*g*epsilon_0;

            cost = pt_K1->COST_IC_K1 / ((g2) * pt_K1->nu_compton_0);

            a = 2.0 * k * log(k) ;

            a = a + (1+2*k)*(1-k);


            c = 0.5*(1-k)*(Gamma_e*k)*(Gamma_e*k)/(1+4.0*k*Gamma_e);




            rate = a+c;
            rate *= cost;
            //printf("2\n");
        }
        else{
         rate=0;
         //printf("3\n");
        }
    }


    //==========================================
    //CHI is q in Eq 4.1
    //q=  nu_1/((4*g^2*nu_compton_0)*(1-h*nu_1/(h*mec^2*g))
    //CHI=nu_1/((nu_compton_0)*(4*g*g - 4*g*h*nu_1/(h*mec^2))
    //CHI=nu_1/(nu_compton_0*(4*g*g-k_epsilon_1))
    //==========================================
    //CHI = pt_K1->nu_1 / (pt_K1->nu_compton_0 * (4.0 * g2 - k_epsilon_1));
    //cost = pt_K1->COST_IC_K1 / ((g2) * pt_K1->nu_compton_0);

    //=========================================
    //bracket term in Eq. 4.1
    //braket term=
    //2*q*ln(q)+(1+2*q)*(1-q)+0.5*(1-q)*(k_epsilon_0*q)^2/(1+k_epsilon_0*q)
    //a=2*CHI*ln(CHI)+(1+2*CHI)*(1-CHI)=2*CHI*ln(CHI)+(1+CHI-2*CHI^2)
    //c=CHI*k_epsilon_0
    //b=c^2/(1+c)
    //
    //braket term=a+b*0.5*(1-CHI)
    //==========================================
    //This piece of Eq. 4.1 :
    //2.0*CHI*log(CHI)+(1.0+CHI-2.0*CHI*CHI);
    //can be sub by a 4th polinomial regression
    //a= 0.81847 - 2.5358 * CHI + 5.1355 * CHI*CHI - 5.0344 * CHI*CHI*CHI + 1.62 * CHI*CHI*CHI*CHI;
    
    //==========================================
    //a = 2.0 * CHI * log(CHI)+(1.0 + CHI - 2.0 * CHI * CHI);
    //c = (CHI * k_epsilon_0);
    //b = (c * c) / (1.0 + c);
    //a += b * 0.5 * (1.0 - CHI);
    //a *= cost;

    return rate;
}
//=========================================================================================

void set_N_distr_for_Compton(struct spettro * pt, double nu_in, double nu_out, int stat_frame)
{
    double epsilon_0, epsilon_1, g_min_BG;
    epsilon_0 = HPLANCK * nu_in * one_by_MEC2;
    epsilon_1 = HPLANCK * nu_out * one_by_MEC2;
    
    g_min_BG = 0.5 * epsilon_1 *(1 + sqrt(1.0 + (1.0 / (epsilon_1 * epsilon_0))));

    
    if (pt->EC_stat == 1)
    {
        g_min_BG = g_min_BG / pt->beam_obj;
    }
    if (g_min_BG > pt->gmin_griglia)
    {
        Fill_Ne_IC(pt, g_min_BG, stat_frame);
    }
    else
    {
        Fill_Ne_IC(pt, pt->gmin_griglia, stat_frame);
    }
}

//=========================================================================================
// INTEGRAZIONE SSC TRAPEZOIDALE/METODO DI SIMPSON E GRIGLIA EQUI-LOG
// TEST MODIFICO INTERFACCIA
// nu_Sync=>nu_seed
// I_nu_Sync=>I_nu_seed
// a,b: boundaries for photon integration
//=========================================================================================
double integrale_IC(double (*pf) (struct spettro *, double x), struct spettro * pt, double a, double b, int stat_frame) {
    double nu1, nu2, integr_gamma, integr_nu;
    double g3, g1, y_g1, y_g2, y_g3, y_nu1, y_nu2, g;
    double delta_g, delta_nu;
    unsigned long i;
    double (*pf_K1) (struct spettro *, double x);

    pf_K1 = &f_compton_K1;
    integr_gamma = 0.0;
    integr_nu = 0.0;
    i = 0;

   
    
    g = pt->gmin_griglia;

    //wit this choice gmin is set to its lowest possible
    //value for enay nu_seed in the range (a,b)
    set_N_distr_for_Compton(pt, b, pt->nu_1, stat_frame);

    while (pt->nu_seed[i] < a && i<pt->nu_seed_size) {
        i++;
    }

    if (pt->verbose>1) {
        printf("***** Integrale  IC ******\n");
        printf("i=%d\n", i);
        printf("nu=%e a=%e b=%e  g_min_grid=%e g_max_grid=%e\n", pt->nu_seed[i], a, b, pt->griglia_gamma_Ne_log_IC[0], pt->griglia_gamma_Ne_log_IC[pt->gamma_grid_size - 1]);
    }

    nu1=0.;
    y_nu1=0.;

    if (i<pt->nu_seed_size){
        nu1 = pt->nu_seed[i];
        y_nu1 = pt->n_seed[i];
    }


    while (pt->nu_seed[i + 1] <= b && pt->nu_seed[i + 1] >= a && i<pt->nu_seed_size-1) {


        integr_gamma = 0.0;

        g1 = pt->griglia_gamma_Ne_log_IC[0];
        pt->nu_compton_0 = pt->nu_seed[i];

        if (pt->adaptive_e_binning ==1){
            //wit this choice gmin is set to its lowest valeu
            //the actual nu_seed value
            set_N_distr_for_Compton(pt, pt->nu_compton_0, pt->nu_1, stat_frame);
        }
        
        y_g1 = f_compton_K1(pt, g1) * pt->Ne_IC[0];



        for (pt->i_griglia_gamma = 1; pt->i_griglia_gamma < pt->gamma_grid_size - 1; pt->i_griglia_gamma++) {



            y_g2 = f_compton_K1(pt, pt->griglia_gamma_Ne_log_IC[pt->i_griglia_gamma]) * pt->Ne_IC[pt->i_griglia_gamma];

            pt->i_griglia_gamma++;
            g3 = pt->griglia_gamma_Ne_log_IC[pt->i_griglia_gamma];
            y_g3 = f_compton_K1(pt, g3) * pt->Ne_IC[pt->i_griglia_gamma];

            delta_g = (g3 - g1);
            integr_gamma += (delta_g)*(y_g1 + 4.0 * y_g2 + y_g3);

            y_g1 = y_g3;
            g1 = g3;


        }


        nu2 = pt->nu_seed[i + 1];
        y_nu2 = pt->n_seed[i + 1];
        delta_nu = nu2 - nu1;

        integr_nu += (y_nu2 + y_nu1) * delta_nu * (integr_gamma);
        nu1 = nu2;
        y_nu1 = y_nu2;
        i++;

    }


    pt->gmin_griglia=g;



    //============================================================
    //0.75 fattore di correzione di GOULD
    //has been moved to spetto_sincrotrone.c

    //0.5/3.0 viene dall'int simpson in gamma
    //lo 0.5 viene dalla regola del trapezio dell'integrale in nu
    //============================================================
    return (integr_nu * 0.5)*(0.5 / 3.0);
}
//=========================================================================================





//=========================================================================================
// Cooling Compton
//=========================================================================================
double compton_cooling(struct spettro *pt_spec, struct temp_ev *pt_ev, double gamma) {
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
    if (pt_ev->do_SSC_cooling) {
        if (pt_spec->verbose>1) {
            printf("nu_start_Sync=%e\n", pt_spec->nu_start_Sync);
            printf("nu_stop_Sync_ssc=%e\n", pt_spec->nu_stop_Sync_ssc);

        }
         
        pt_spec->nu_seed = pt_spec->nu_Sync;
        pt_spec->n_seed = pt_spec->n_Sync;
        //0.75 is the Gould correction factor
        //has been moved to spetto_sincrotrone.c
        comp_cooling += integrale_IC_cooling(pt_spec,
                pt_spec->nu_start_Sync,
                pt_spec->nu_stop_Sync_ssc,
                gamma);
        //printf("evaluate IC cooling, gamma=%e cooling_rate=%e, Sync_cooling_rate_ratio=%e\n",gamma,comp_cooling,comp_cooling/Sync_cool(pt_spec->B,gamma));
    }

    //EC Disk
    if (pt_ev->do_EC_cooling_Disk == 1 ) {

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
    if (pt_ev->do_EC_cooling_BLR == 1 ) {

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
    if (pt_ev->do_EC_cooling_DT == 1 ) {

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
    if (pt_ev->do_EC_cooling_Star == 1 ) {

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
    if (pt_ev->do_EC_cooling_CMB == 1 ) {

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
double integrale_IC_cooling(struct spettro * pt, double a, double b, double gamma) {
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
