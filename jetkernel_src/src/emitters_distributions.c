//==============================================================================
//  FUNZIONI CHE COSTRUISCONO LA DISTRIBUZIONE ELETTRONICA
//==============================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file distrib_elettr.c
 * \author Andrea Tramacere
 * \date 19-09-2004
 * \brief funzioni per la
 * distribuzione energetica
 * degli elettroni sia nel caso
 * stazionario che per ET
 *
 */

//==============================================================
/// Genera la griglia su gamma che viene usata per integrare in gamma il
// il Sync ed l'IC. I punti con indice pari sono equispaziati nel logar
// tmo fra di loro, quelli con indice dispari sono la media fra il prece
// dente pari ed il successivo pari
//==============================================================

void Genera_griglia_gamma_N_log(struct blob *pt, double * griglia_gamma_N_log, double gmin_griglia, double gmax_griglia) {
	unsigned int i;
    double delta_log;
    double log_a, log_b;
    if (pt->verbose>1) {
        printf("Generete log gamma_grid for N \n");
        printf("size is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
    }
    
    log_a = log10(gmin_griglia);
    log_b = log10(gmax_griglia);
    delta_log = (log_b - log_a) / ((double) pt->gamma_grid_size - 1);
    //PUNTI CON INDICE PARI LOG
    for (i = 0; i < pt->gamma_grid_size; i += 2) {
        griglia_gamma_N_log[i] = pow(10, (log_a + delta_log * (double) (i)));
        //printf("i=%d griglia_gamma_Ne_log=%e\n",i,pt->griglia_gamma_Ne_log[i]);
    }
    //PUNTI CON INDICE DISPARI LIN
    for (i = 1; i < pt->gamma_grid_size; i += 2) {
        griglia_gamma_N_log[i] =
                (griglia_gamma_N_log[i - 1] + griglia_gamma_N_log[i + 1])*0.5;
        //printf("i=%d griglia_gamma_Ne_log=%e\n",i,pt->griglia_gamma_Ne_log[i]);
    }
}

void setNgrid(struct blob *pt)
{
    //==========================================
    //Numerical Integration precision Setup
    //==========================================

    double  *gmin, *gmax , *gmin_griglia, *gmax_griglia;
    unsigned int *gamma_grid_size;
    
    if (strcmp(pt->PARTICLE, "secondaries_el") == 0)
    {
      gamma_grid_size = &(pt->gamma_grid_size);
      gmax = &(pt->gmax_secondaries);
      gmin = &(pt->gmin_secondaries);
      gmax_griglia = &(pt->gmax_griglia_secondaries);
      gmin_griglia = &(pt->gmin_griglia_secondaries);
    }
    else{
        gamma_grid_size = &(pt->gamma_grid_size);
        gmax = &(pt->gmax);
        gmin = &(pt->gmin);
        gmax_griglia = &(pt->gmax_griglia);
        gmin_griglia = &(pt->gmin_griglia);

    }
    if (strcmp(pt->MODE, "accurate") == 0)
    {
        *gamma_grid_size = 10000;
        if (pt->verbose)
        {
            printf("gamma mesh set to value=%d for accurate integration \n", *gamma_grid_size);
        }
    }
    else if (strcmp(pt->MODE, "fast") == 0)
    {
        *gamma_grid_size = 1000;
        if (pt->verbose)
        {
            printf("gamma mesh set to value=%d for fast integration, \n", *gamma_grid_size);
        }
    }
    else if (strcmp(pt->MODE, "custom") == 0)
    {
        if (pt->verbose)
        {
            printf("gamma mesh set to custom value=%d  \n", *gamma_grid_size);
        }
    }
    else
    {
        if (pt->verbose)
        {
            printf("MODE set to wrong value: %s, allowed= accurate,fast,custom", pt->MODE);
            exit(1);
        }
    }

    if ( (int)(*gamma_grid_size)%2 == 0)
    {
        (*gamma_grid_size) ++;
        if (pt->verbose)
        {
            printf("!! gamma_grid_size has to be odd\n");
            printf("!! pt->gamma_grid_size=%d\n", (*gamma_grid_size));
        }
    }

    //=========================================
    // check on gamma grid
    //=========================================
    // gamma min griglia
    if (*gmin_griglia < 0.0 || *gmin < *gmin_griglia)
    {
        //		if(pt->gmin>2.0){
        //			pt->gmin_griglia=pt->gmin/2.0;
        //		}
        //		else{
        //		   pt->gmin_griglia=1.0;
        //		}
        *gmin_griglia = *gmin;
    }

    if (*gmax_griglia < 0.0 || *gmax > *gmax_griglia)
    {
        *gmax_griglia = *gmax;
    }

    if (*gmin < *gmin_griglia)
    {
        printf("gmin < gmin_griglia, it must be the oppsosite");
        exit(1);
    }
    if (*gmax > *gmax_griglia)
    {
        printf("gmax > gmax_griglia, it must be the oppsosite");
        exit(1);
    }

    if (pt->verbose > 1)
    {
        printf("Set array per Ne \n");
        printf("elements number is pt->gamma_grid_size=%d\n", *gamma_grid_size);
    }

    if (pt->grid_bounded_to_gamma == 1)
    {
        *gmax_griglia = *gmax;
        *gmin_griglia = *gmin;
        if (strcmp(pt->PARTICLE, "secondaries_el") == 0)
        {
            *gmin_griglia=1.0;
        }
    }
}

//========================================
// Genera la  N[i] per e-
//========================================

void build_Ne(struct blob *pt) {
   

    //printf("Set array per Ne %s \n",pt->DISTR);
    //printf("build_Ne Set array per Ne 1 %s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_Ne_log),pt->gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log, pt->gmin_griglia, pt->gmax_griglia);
    //printf("build_Ne Set array per Ne 2%s \n",pt->DISTR);
    alloc_N_distr(&(pt->Ne),pt->gamma_grid_size);

    //stationary frame
    //printf("stationary build_Ne Set array per Ne 3%s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_Ne_log_stat),pt->gamma_grid_size);
    //printf("stationarybuild_Ne Set array per Ne 4%s \n",pt->DISTR);
    alloc_N_distr(&(pt->Ne_stat),pt->gamma_grid_size);
    
    //if (pt->verbose>1) {
    //    printf("DONE \n");
    //}

    //N for IC and simpson equilog integration over N_grid
    //printf("N for IC and simpson build_Ne Set array per Ne 5%s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_Ne_log_IC),pt->gamma_grid_size);
    //printf("N for IC and simpson build_Ne Set array per Ne 6%s \n",pt->DISTR);
    alloc_N_distr(&(pt->Ne_IC),pt->gamma_grid_size);
    //printf("N for IC and simpson build_Ne Set array per Ne 7%s \n",pt->DISTR);
    alloc_N_distr(&(pt->Integrand_over_gamma_grid),pt->gamma_grid_size);

}

void build_Ne_secondaries(struct blob *pt) {
    //printf("Set array per Ne %s \n",pt->DISTR);
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_Ne_log),pt->gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log,pt->gmin_griglia_secondaries, pt->gmax_griglia_secondaries);
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->Ne),pt->gamma_grid_size);

    //stationary frame
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_Ne_log_stat),pt->gamma_grid_size);
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->Ne_stat),pt->gamma_grid_size);
    //if (pt->verbose>1) {
    //    printf("DONE \n");
    //}

    //N for IC and simpson equilog integration over N_grid
    alloc_N_distr(&(pt->griglia_gamma_Ne_log_IC),pt->gamma_grid_size);
    alloc_N_distr(&(pt->Ne_IC),pt->gamma_grid_size);
    alloc_N_distr(&(pt->Integrand_over_gamma_grid),pt->gamma_grid_size);

}

void build_Q_inj_e_second(struct blob *pt) {
    //printf("build_Q_inj_e_second Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->Q_inj_e_second),pt->gamma_grid_size);
}

void build_Np(struct blob *pt)
{
    //printf("build_Np Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_Np_log), pt->gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Np_log,pt->gmin_griglia, pt->gmax_griglia);
    //printf("build_Np Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->Np), pt->gamma_grid_size);
}

void build_Np_jetset(struct blob *pt) {
    //printf("build_Np_jetset Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_jetset_Np_log), pt->gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_jetset_Np_log,pt->gmin_griglia, pt->gmax_griglia);
    //printf("build_Np_jetset Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->Np_jetset), pt->gamma_grid_size);

}

void build_Ne_jetset(struct blob *pt) {
    //printf("build_Ne_jetset Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->griglia_gamma_jetset_Ne_log),pt->gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_jetset_Ne_log, pt->gmin_griglia, pt->gmax_griglia);
    //printf("build_Ne_jetset Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->Ne_jetset),pt->gamma_grid_size);
}


void Fill_Ne_IC(struct blob *pt, double g_min_IC, int stat_frame) {
    unsigned int i,i_start;
    //double gmin_grid;
    i_start=0;
    while (pt->griglia_gamma_Ne_log[i_start] < g_min_IC && i_start < pt->gamma_grid_size) {
        i_start++;
    }
    if (i_start % 2 != 0) {
        i_start = max(0,i_start-1);
    }
    g_min_IC=pt->griglia_gamma_Ne_log[i_start];
    
    if (pt->verbose>1) {
        printf("Set array per Ne IC\n");
        printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
    }
    //printf("Set array per Ne %s \n",pt->DISTR);
    
    if (strcmp(pt->PARTICLE, "protons") == 0) {
        
        if(pt->IC_adaptive_e_binning ==1){
            Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log_IC,g_min_IC, pt->gmax_griglia_secondaries);
        }else{
            Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log_IC,pt->gmin_griglia_secondaries, pt->gmax_griglia_secondaries);
        }
        
    }
    else{
        if(pt->IC_adaptive_e_binning ==1){
            //gmin_grid=max(pt->gmin_griglia,g_min_IC/10);
            Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log_IC,g_min_IC, pt->gmax_griglia);
        }else{
            Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log_IC,pt->gmin_griglia, pt->gmax_griglia);
        }
    }
    SetDistr(pt);
    for (i = 0; i < pt->gamma_grid_size; i++) {
        if(pt->IC_adaptive_e_binning ==1){
            if (pt->griglia_gamma_Ne_log_IC[i]>=g_min_IC){
             pt->Ne_IC[i] = N_distr_interp(pt->gamma_grid_size,
                                    pt->griglia_gamma_Ne_log_IC[i],
                                    pt->griglia_gamma_Ne_log,
                                    pt->Ne);
            }else{
                 pt->Ne_IC[i]=0;
                }
        }else{
            pt->Ne_IC[i] = pt->Ne[i]; 
        }                 
        if (stat_frame==1){
            //the delta^2 in Ne_stat is also correct because we use electron density
            //so the relativistic invariant is
            //N/(V*gamma^2)=N'/(V'gamma'2^)
            pt->Ne_IC[i]*=pt->beam_obj*pt->beam_obj;

            //This transformation is correct
            //the grid is shifted by a factor of delta, hence the integration
            //boundaries are properly updated but the value of N[i] is still the
            //value of N(gamma') as in the formula 6.133 in Dermer&Menon
            pt->griglia_gamma_Ne_log_IC[i]*=pt->beam_obj;
        }
    }
}



void build_Ne_custom(struct blob *pt,  unsigned int size) {
    pt->gamma_custom_grid_size=size;
    if (pt->verbose>1) {
        printf("Set array for Ne for from_array mode \n");
        printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
    }
    //printf("build_Ne_custom Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->gamma_e_custom),size);
    alloc_N_distr(&(pt->Ne_custom),size);

}

void build_Np_custom(struct blob *pt,  unsigned int size) {
    pt->gamma_custom_grid_size=size;
    if (pt->verbose>1) {
        printf("Set array for Np for from_array mode \n");
        printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
    }
    //printf("build_Np_custom Set array per Ne %s \n",pt->DISTR);
    alloc_N_distr(&(pt->gamma_p_custom),size);
    alloc_N_distr(&(pt->Np_custom),size);

}



void InitNe(struct blob *pt){
    //double (*pf_distr)(struct blob *, double x);
    //pf_distr = &N_distr_integranda;

    //printf("==> InitNe start\n");
    //printf("==> setNgrid\n");
    setNgrid(pt);
    //printf("==> build_Ne\n");
    build_Ne(pt);
    //printf("==> SetDistr\n");
    SetDistr(pt);
    //printf("==> SetDistr\n");
    Fill_N(pt, pt->griglia_gamma_Ne_log, pt->Ne);


	//This flag is set to 1 to know that
	pt->Distr_e_done = 1;

    pt->N_0e = pt->N_0;
    //printf("==> N_tot\n");
    pt->N_e  = N_tot(pt, N_distr_integranda);
    //printf("==> InitNe stop\n");
}


//========================================
// Genera la  N[i] per pp ed e- secondari
//========================================
void Init_Np_Ne_pp(struct blob *pt)
{
    //char *name;
    //double (*pf_distr) (struct blob *, double x);
    //pf_distr = &N_distr_integranda;

    pt->gmin_secondaries=pt->gmin;
    pt->gmax_secondaries=pt->gmax*mp_by_me;
    setNgrid(pt);
    build_Np(pt);
    SetDistr(pt);
    if (pt->verbose>1) {
        printf("********** protons ***********\n");
        printf("set array for Np\n");
        printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
        printf("TIPO_DISTR %d\n", pt->TIPO_DISTR);
    }
    Fill_N(pt, pt->griglia_gamma_Np_log, pt->Np);
    //
    //This flag si set to 1 to know that
    //N(gamma) has been properly initialized and filled
    //printf("-->\n" );
    //printf("--> N0e %e N0 %e N0p %e\n", pt->N_0e, pt->N_0, pt->N_0p);
    
    pt->Distr_p_done = 1;
    pt->N_0p = pt->N_0;
    
    //printf("--> N0e %e N0 %e N0p %e\n", pt->N_0e, pt->N_0, pt->N_0p);
    
    pt->N_p = N_tot(pt, N_distr_integranda);
    //name = "distr-p.dat";
    //Scrivi_N_file(pt, name, pt->griglia_gamma_Np_log, pt->Np);

    // Secondaries e- from pp

    //Set N to e- from pp
    sprintf(pt->PARTICLE, "secondaries_el");
    setNgrid(pt);
    build_Ne_secondaries(pt);
    build_Q_inj_e_second(pt);
    SetDistr(pt);
    Fill_N(pt, pt->griglia_gamma_Ne_log, pt->Q_inj_e_second);
    CooolingEquilibrium(pt,pt->T_esc_e_second);
    pt->Distr_e_done = 1;
    pt->N_0e = pt->N_0;
    

    
    //printf("--> N0e %e N0 %e N0p %e\n", pt->N_0e, pt->N_0, pt->N_0p);
    //printf("-->\n");
    pt->N_e_pp = N_tot(pt, N_distr_integranda);

    //if (pt->verbose > 1)
    //{
    //    printf("****** secondary leptons *****\n");
    //    printf("set array for secondary Ne\n");
    //    printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
    //    printf("N_e_pp =%e\n", pt->N_e_pp);
    //}

    //name = "distr-e-from-pp.dat";
    //Scrivi_N_file(pt, name, pt->griglia_gamma_Ne_log, pt->Ne);
    
    //set back pt->N_0 to the proton value and particle name
    pt->N_0 = pt->N_0p;
    sprintf(pt->PARTICLE, "protons");
    SetDistr(pt);
}





//========================================
// Trova il gmax da N[i]>o
//========================================

double Find_gmax(struct blob *pt, double *N, double *g) {
	unsigned int i;
    double gmax;
    gmax = g[0];

    for (i = 0; i < pt->gamma_grid_size; i++) {
        if (N[i] > 0) {
            gmax = g[i];
        }
    }
    return gmax;
}
//=====================================================







//========================================
// RIEMPIE IL VETTORE  N[i]
//========================================

void Fill_N(struct blob *pt, double * griglia_gamma_N_log, double * N) {
	unsigned int i;
    //integranda Disre e
    double (*pf_norm) (struct blob *, double x);

    pt->N_0 = 1.0;
    //=========================================
    // interpolate custom Ne/p
    //=========================================
    if (pt->TIPO_DISTR == 0)
    {   
        if (strcmp(pt->PARTICLE, "protons") == 0){
            for (i = 0; i < pt->gamma_grid_size; i++)
            
            {
                N[i] = N_distr_interp(pt->gamma_custom_grid_size,
                                    griglia_gamma_N_log[i],
                                    pt->gamma_p_custom,
                                    pt->Np_custom);
            }
        }else{
            for (i = 0; i < pt->gamma_grid_size; i++)
            
            {
                N[i] = N_distr_interp(pt->gamma_custom_grid_size,
                                    griglia_gamma_N_log[i],
                                    pt->gamma_e_custom,
                                    pt->Ne_custom);
            }
        }

    }
    else if (pt->TIPO_DISTR==10){
        if (strcmp(pt->PARTICLE, "protons") == 0){
            for (i = 0; i < pt->gamma_grid_size; i++)
            
            {
                N[i] = pt->Np_jetset[i];
            }
        }else{
            for (i = 0; i < pt->gamma_grid_size; i++)
            
            {
                N[i] = pt->Ne_jetset[i];
            }
        }

    }

    //=========================================
    // fill defined Ne/p
    //=========================================
    else if (pt->TIPO_DISTR != -1){

        //Normalization
        
        
        if (pt->Norm_distr == 1 && pt->TIPO_DISTR != -1)
        {
            pf_norm = &N_distr_integranda;
            pt->N_0 = integrale_trap_log_struct(pf_norm, pt, pt->gmin, pt->gmax, 10000);
        }

        for (i = 0; i < pt->gamma_grid_size; i++)
        {
            N[i] = N_distr(pt, griglia_gamma_N_log[i]);
        }
    }

    //if distr is e- from pp te
    //the distribution is filled with the injection
    //by the function N_distr
    else if (pt->TIPO_DISTR == -1){
        for (i = 0; i < pt->gamma_grid_size; i++)
        {
            N[i] = N_distr(pt, griglia_gamma_N_log[i]);
        }
    }
    else {
        printf("TIPO_DISTR set to wrong value: %d\n",pt->TIPO_DISTR);
        exit(1);
    }

    //pt->Distr_e_done = 1;

}   


double pl_func(double Gamma,double p){
    return   pow(Gamma, -p) ;
}

double plc_func(double Gamma,double gamma_cut,double p){
  return   pow(Gamma, -p) * exp(-Gamma / gamma_cut);
}

double bkn_func(double Gamma,double gamma_break,double p, double p_1){
     double a=0.;
     if (Gamma< gamma_break){
        a= pow(Gamma, -p);
      } else {
        a=  pow(gamma_break, -(p -p_1)) * pow(Gamma, -p_1);
      }

     return a;
}

double lp_func(double Gamma,double gamma0,double r, double s){
    return  pow((Gamma / (gamma0)),(-s -r*log10(Gamma /gamma0 )));
}


double lp_ep_func(double Gamma,double gamma_p,double r){
    return  pow(10, (-r * pow(log10(Gamma /gamma_p), 2)));
}

double lppl_func(double Gamma,double gamma0, double r, double s){
    double a=0.;


    if (Gamma < gamma0) {
        a= pow((Gamma / gamma0), (- s));
    }
    if (Gamma >= gamma0 ) {
        a= pow((Gamma / (gamma0)),
                (-s -r*log10(Gamma /gamma0 )));
    }

    //printf("Gamma=%e %e %e\n", Gamma,gamma0,a);
    return a;
    }


//double pile_up_ratio(double Gamma,double sigma,double gamma_eq){
//     return pow(gamma_eq*,-sigma-2.)*st_gamma(sigma-1.)/st_gamma(2.*sigma+2.);
//}

double pile_up_func(double Gamma, double gamma_eq, double alpha){
    return Gamma*Gamma*exp(-pow((Gamma/gamma_eq),alpha)/alpha);
}

double lppl_pile_up_func(double Gamma,double gamma0, double gamma_inj,double r, double s,double gamma_eq, double ratio_pile_up ,double alpha){
    double a,b,s1;
    //ratio_pile_up;
    b=0;
    a=0;
    
    //ratio_pile_up=lppl_func(gamma_cut,gamma0, r, s);
    //ratio_pile_up*=1.0/pile_up_func(gamma_cut,gamma_eq,alpha);
    //ratio_pile_up=pow(gamma_eq,-s-2.0)*st_gamma(s-1.0)/st_gamma(2*s+2.0);
    //printf("gamma_inj %e \n",gamma_inj);
    s1=s+0.5;
    if (Gamma< gamma_inj){
        b= pow(Gamma/gamma0,s1);
    } else {
        b= pow(gamma_inj/gamma0,s1+s)*lppl_func(Gamma,gamma0, r, s);
    }
    
    
    
    //a=  pow(gamma_inj,2.0*p_1+1)*(2.*p_1+1.);
    a= ratio_pile_up*pile_up_func(Gamma,gamma_eq,alpha);



    return (a+b);
}

double bkn_pile_up_func(double Gamma,double gamma_inj, double p, double p_1,double gamma_eq, double gamma_cut ,double alpha){
    double a,b,ratio_pile_up;
    b=0;
    a=0;
    //ratio_pile_up=bkn_func(gamma_cut,gamma_break, p, p_1);
    //ratio_pile_up*=1.0/pile_up_func(gamma_cut,gamma_eq,alpha);

    //ratio_pile_up=st_gamma(p_1-1.0)/st_gamma(2*p_1+2.0)*pow(gamma_eq,-p_1-2.0);
    //if (Gamma < gamma_cut*0.01) {

//        b=bkn_func(Gamma,gamma_break, p, p_1);
  //  }



    //if (Gamma >=gamma_cut*0.01 && Gamma<gamma_eq*10){

    if (Gamma<gamma_inj){

         b= 1.0/(2.0*p_1+1.0)*pow(gamma_inj,-p_1-2.0)*pow(Gamma,p_1+1);
    }
    else {
         b= 1.0/(2.0*p_1+1.0)*pow(gamma_inj,p_1-1.0)*pow(Gamma,-p_1)*exp(-(Gamma/gamma_cut));


         ratio_pile_up=st_gamma(p_1-1.)/st_gamma(2.*p_1+2.);
         ratio_pile_up*= pow(gamma_inj,p_1-1.0)*pow(gamma_eq,-p_1-2.0);
         a= ratio_pile_up*pile_up_func(Gamma,gamma_eq,alpha);
   }


    return (a+b);
}


double spit_func(double Gamma,double gamma_th,double temp, double index){
    double a,b,c,f=0.;
    if (Gamma < gamma_th) {
        a=(Gamma*Gamma)/(2.0*temp*temp*temp);
        b=exp(-(Gamma/temp));

        f= a*b;
    } else {
        a=( gamma_th* gamma_th)/(2.0* temp* temp* temp);
        b=exp(-( gamma_th/ temp));
        c= pow((Gamma/ gamma_th), -index);
        f= a*b*c;

    }

    return f;
}



//==============================================================
// N_distr
//==============================================================

double N_distr(struct blob *pt_N, double Gamma) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * questa funzione restituisce la distribuzione energetica                          \n
     * richiesti dalla funzione chiamante. Le funzioni chiamanti sono quelle per        \n
     * il calcolo degli spettri di sincrotrone ed di IC/EC                              \n
     *
     */

    double a ;


    a=0.;

    if (Gamma >= pt_N->gmin_secondaries && Gamma <= pt_N->gmax_secondaries && pt_N->TIPO_DISTR == -1) {
        //pt_N->Gamma = Gamma;
        a= vluce_cm * pt_N->NH_pp * MEC2_TeV * bn_to_cm2 * rate_electrons_pp(pt_N, Gamma);
    }else{

        a= N_distr_integranda(pt_N,Gamma)*pt_N->N/pt_N->N_0;

    }


    return a;


}

double N_tot(struct blob *pt, double (*pf_distr)(struct blob *, double x))
{
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * questa funzione restituisce il numero tototale di particelle \n
     *
     */

    double a;
    a = 0.;

    a= integrale_trap_log_struct(pf_distr,
                                pt,
                                pt->gmin,
                                pt->gmax,
                                10000);

    //if the distr is not secondaries or interpolated
    if (pt->TIPO_DISTR > 0)
    {
        a = a * pt->N / pt->N_0;
    }

    return a;
}

//==============================================================
//   funzione integranda per la distribuzione degli e-
//    per calcolare il coeff di norm
//==============================================================

double N_distr_integranda(struct blob *pt_N, double Gamma) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * funzione che resitituisce le integrande per il calcolo  \n
     * del coefficiente di normalizzazione delle distribuzioni \n
     * elettroniche statiche                                   \n
     */

    double a;
    a=0.;

    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax){

        //Secondaris e Distribution has not analytical expression
        //it is taken from the N array, throug log-lin interpolation
        if (  pt_N->TIPO_DISTR == -1) {
            a= N_distr_interp(pt_N->gamma_grid_size,
                              Gamma,
                              pt_N->griglia_gamma_Ne_log,
                              pt_N->Ne);
        }

        if (  pt_N->TIPO_DISTR == 0) {
            a= N_distr_interp(pt_N->gamma_custom_grid_size,
                                Gamma,
                                pt_N->gamma_e_custom,
                                pt_N->Ne_custom);
        }

        //PL

        if (  pt_N->TIPO_DISTR == 1) {
            a= pl_func(Gamma, pt_N->p) ;
        }

        //PL CON EXP CUTOFF
        if (  pt_N->TIPO_DISTR == 2) {
            a=  plc_func(Gamma, pt_N->gamma_cut, pt_N->p);
        }

        //PL BROCKEN
        if (  pt_N->TIPO_DISTR == 3) {
           a= bkn_func(Gamma, pt_N->gamma_break,pt_N->p,pt_N->p_1);
        }

         //LOG PARABOLA
        if (  pt_N->TIPO_DISTR == 4) {
            a= lp_func(Gamma,pt_N->gamma0_log_parab,pt_N->r,pt_N->s);

        }


        //LOG PARABOLA CON PICCO
        if (  pt_N->TIPO_DISTR == 5) {
            a= lp_ep_func(Gamma,pt_N->gammap_log_parab,pt_N->r);
        }

        //LOG PARABOLA CON PL
        if (  pt_N->TIPO_DISTR == 6) {
            a= lppl_func(Gamma,pt_N->gamma0_log_parab,pt_N->r,pt_N->s);

        }

        //Spit
        //double emin, norm, x, x_min;
        if (  pt_N->TIPO_DISTR == 7) {
            a= spit_func(Gamma,pt_N->spit_gamma_th,pt_N->spit_temp,pt_N->spit_index);
        }



        //LOG PARABOLA CON PL e PILE-UP
        if (  pt_N->TIPO_DISTR == 8) {
            a= lppl_pile_up_func( Gamma,pt_N->gamma0_log_parab,pt_N->gamma_inj,pt_N->r,pt_N->s,pt_N->gamma_pile_up, pt_N->ratio_pile_up ,pt_N->alpha_pile_up);
        }

         //PL BROCKEN-PILEUP
        if (  pt_N->TIPO_DISTR == 9) {
            a= bkn_pile_up_func(Gamma,pt_N->gamma_break,pt_N->p,pt_N->p_1,pt_N->gamma_pile_up, pt_N->gamma_pile_up_cut ,pt_N->alpha_pile_up);
        }




    }


    return a;
}

double N_distr_interp(unsigned int size, double Gamma, double *griglia_gamma, double *N) {
	//size: input grid size
    //Gamme: output  gamma
    //griglia_gamma: input gamma_grid
    //N: input N
    unsigned int i;
    double gamma_piu, gamma_meno, Npiu, Nmeno, g, a;
    i = 0;
    while (griglia_gamma[i] < Gamma && i < size) {
        i++;
    }
    //i--;
    //printf("G=%e G_file=%e\n",pt->griglia_gamma_Ne_log[i],G_File[count]);
    if (i > 0 && i < size && N[i] > 0 && N[i - 1] > 0) {
        gamma_piu = log10(griglia_gamma[i]);
        gamma_meno = log10(griglia_gamma[i - 1]);
        Npiu = log10(N[i]);
        Nmeno = log10(N[i - 1]);
        g = log10(Gamma);
        a = ((g - gamma_meno) / (gamma_piu - gamma_meno))*(Npiu - Nmeno);
        a += Nmeno;
        //printf("%d %e %e %e %e %e %e\n",i,gamma_piu,gamma_meno,N[i],Npiu,Nmeno,a);
        return pow(10, a);
    } else {
        return 0;
    }
}


void alloc_N_distr(double ** pt,int size){
        //printf("pre %p\n",*pt);
        //printf("alloc n\n");
        //if (*pt==NULL){
        //   printf("is  NULL\n");
        //}
        if (*pt){
            //printf("freeing\n");
            //printf("%e\n",pt[0]);
            free(*pt);
            //printf("free\n");
        }

        *pt = calloc(size, sizeof (double));
        //*pt= mallot(size * sizeof (double));
        //printf("post %p\n",*pt);

    }

//=========================================================================================

void SetDistr(struct blob *pt) {
    //-1 is for secondary e- coming from pp

    /*** Associo ad ogni distribuzione di elettroni ***/

    if (strcmp(pt->PARTICLE, "secondaries_el") == 0)
    {
        pt->TIPO_DISTR = -1;
    }
    else
    {
        if (strcmp(pt->DISTR, "from_array") == 0)
        {
            pt->TIPO_DISTR = 0;
        }

        if (strcmp(pt->DISTR, "pl") == 0)
        {
            pt->TIPO_DISTR = 1;
        }

        if (strcmp(pt->DISTR, "plc") == 0)
        {
            pt->TIPO_DISTR = 2;
        }

        if (strcmp(pt->DISTR, "bkn") == 0)
        {
            pt->TIPO_DISTR = 3;
        }

        if (strcmp(pt->DISTR, "lp") == 0)
        {
            pt->TIPO_DISTR = 4;
        }

        if (strcmp(pt->DISTR, "lpep") == 0)
        {
            pt->TIPO_DISTR = 5;
        }

        if (strcmp(pt->DISTR, "lppl") == 0)
        {
            pt->TIPO_DISTR = 6;
        }

        if (strcmp(pt->DISTR, "spitkov") == 0)
        {
            pt->TIPO_DISTR = 7;
        }

        if (strcmp(pt->DISTR, "lppl_pile_up") == 0)
        {
            pt->TIPO_DISTR = 8;
        }

        if (strcmp(pt->DISTR, "bkn_pile_up") == 0)
        {
            pt->TIPO_DISTR = 9;
        }

        if (strcmp(pt->DISTR, "jetset") == 0)
        {
            pt->TIPO_DISTR = 10;
        }

    }
    

   





    //if (pt->verbose) {
     //printf("tipo di distribuzione %d\n",pt->TIPO_DISTR);
    //}
}
//=========================================================================================
