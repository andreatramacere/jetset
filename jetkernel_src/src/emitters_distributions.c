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
    if (pt->core.verbose>1) {
        printf("Generete log gamma_grid for N \n");
        printf("size is pt->emitters.gamma_grid_size=%d\n", pt->emitters.gamma_grid_size);
    }
    
    log_a = log10(gmin_griglia);
    log_b = log10(gmax_griglia);
    delta_log = (log_b - log_a) / ((double) pt->emitters.gamma_grid_size - 1);
    //PUNTI CON INDICE PARI LOG
    for (i = 0; i < pt->emitters.gamma_grid_size; i += 2) {
        griglia_gamma_N_log[i] = pow(10, (log_a + delta_log * (double) (i)));
        //printf("i=%d griglia_gamma_Ne_log=%e\n",i,pt->emitters.griglia_gamma_Ne_log[i]);
    }
    //PUNTI CON INDICE DISPARI LIN
    for (i = 1; i < pt->emitters.gamma_grid_size; i += 2) {
        griglia_gamma_N_log[i] =
                (griglia_gamma_N_log[i - 1] + griglia_gamma_N_log[i + 1])*0.5;
        //printf("i=%d griglia_gamma_Ne_log=%e\n",i,pt->emitters.griglia_gamma_Ne_log[i]);
    }
}

void setNgrid(struct blob *pt)
{
    //==========================================
    //Numerical Integration precision Setup
    //==========================================

    double  *gmin, *gmax , *gmin_griglia, *gmax_griglia;
    unsigned int *gamma_grid_size;
    
    if (strcmp(pt->core.PARTICLE, "secondaries_el") == 0)
    {
      gamma_grid_size = &(pt->emitters.gamma_grid_size);
      gmax = &(pt->emitters.gmax_secondaries);
      gmin = &(pt->emitters.gmin_secondaries);
      gmax_griglia = &(pt->emitters.gmax_griglia_secondaries);
      gmin_griglia = &(pt->emitters.gmin_griglia_secondaries);
    }
    else{
        gamma_grid_size = &(pt->emitters.gamma_grid_size);
        gmax = &(pt->emitters.gmax);
        gmin = &(pt->emitters.gmin);
        gmax_griglia = &(pt->emitters.gmax_griglia);
        gmin_griglia = &(pt->emitters.gmin_griglia);

    }
    if (strcmp(pt->core.MODE, "accurate") == 0)
    {
        *gamma_grid_size = 10000;
        if (pt->core.verbose)
        {
            printf("gamma mesh set to value=%d for accurate integration \n", *gamma_grid_size);
        }
    }
    else if (strcmp(pt->core.MODE, "fast") == 0)
    {
        *gamma_grid_size = 1000;
        if (pt->core.verbose)
        {
            printf("gamma mesh set to value=%d for fast integration, \n", *gamma_grid_size);
        }
    }
    else if (strcmp(pt->core.MODE, "custom") == 0)
    {
        if (pt->core.verbose)
        {
            printf("gamma mesh set to custom value=%d  \n", *gamma_grid_size);
        }
    }
    else
    {
        if (pt->core.verbose)
        {
            printf("MODE set to wrong value: %s, allowed= accurate,fast,custom", pt->core.MODE);
            exit(1);
        }
    }

    if ( (int)(*gamma_grid_size)%2 == 0)
    {
        (*gamma_grid_size) ++;
        if (pt->core.verbose)
        {
            printf("!! gamma_grid_size has to be odd\n");
            printf("!! pt->emitters.gamma_grid_size=%d\n", (*gamma_grid_size));
        }
    }

    //=========================================
    // check on gamma grid
    //=========================================
    // gamma min griglia
    if (*gmin_griglia < 0.0 || *gmin < *gmin_griglia)
    {
        //		if(pt->emitters.gmin>2.0){
        //			pt->emitters.gmin_griglia=pt->emitters.gmin/2.0;
        //		}
        //		else{
        //		   pt->emitters.gmin_griglia=1.0;
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

    if (pt->core.verbose > 1)
    {
        printf("Set array per Ne \n");
        printf("elements number is pt->emitters.gamma_grid_size=%d\n", *gamma_grid_size);
    }

    if (pt->emitters.grid_bounded_to_gamma == 1)
    {
        *gmax_griglia = *gmax;
        *gmin_griglia = *gmin;
        if (strcmp(pt->core.PARTICLE, "secondaries_el") == 0)
        {
            *gmin_griglia=1.0;
        }
    }
}

//========================================
// Genera la  N[i] per e-
//========================================

void build_Ne(struct blob *pt) {
   

    //printf("Set array per Ne %s \n",pt->core.DISTR);
    //printf("build_Ne Set array per Ne 1 %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.griglia_gamma_Ne_log),pt->emitters.gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->emitters.griglia_gamma_Ne_log, pt->emitters.gmin_griglia, pt->emitters.gmax_griglia);
    //printf("build_Ne Set array per Ne 2%s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.Ne),pt->emitters.gamma_grid_size);

    //stationary frame
    //printf("stationary build_Ne Set array per Ne 3%s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.griglia_gamma_Ne_log_stat),pt->emitters.gamma_grid_size);
    //printf("stationarybuild_Ne Set array per Ne 4%s \n",pt->core.DISTR);
    //alloc_N_distr(&(pt->Ne_stat),pt->emitters.gamma_grid_size);
    
    //if (pt->core.verbose>1) {
    //    printf("DONE \n");
    //}

    //N for IC and simpson equilog integration over N_grid
    //printf("N for IC and simpson build_Ne Set array per Ne 5%s \n",pt->core.DISTR);
    //alloc_N_distr(&(pt->griglia_gamma_Ne_log_IC),pt->emitters.gamma_grid_size);
    //printf("N for IC and simpson build_Ne Set array per Ne 6%s \n",pt->core.DISTR);
    //alloc_N_distr(&(pt->Ne_IC),pt->emitters.gamma_grid_size);
    //printf("N for IC and simpson build_Ne Set array per Ne 7%s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.Integrand_over_gamma_grid),pt->emitters.gamma_grid_size);

}

void build_Ne_secondaries(struct blob *pt) {
    //printf("Set array per Ne %s \n",pt->core.DISTR);
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.griglia_gamma_Ne_log),pt->emitters.gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->emitters.griglia_gamma_Ne_log,pt->emitters.gmin_griglia_secondaries, pt->emitters.gmax_griglia_secondaries);
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.Ne),pt->emitters.gamma_grid_size);

    //stationary frame
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.griglia_gamma_Ne_log_stat),pt->emitters.gamma_grid_size);
    //printf("build_Ne_secondaries Set array per Ne %s \n",pt->core.DISTR);
    //alloc_N_distr(&(pt->Ne_stat),pt->emitters.gamma_grid_size);
    //if (pt->core.verbose>1) {
    //    printf("DONE \n");
    //}

    //N for IC and simpson equilog integration over N_grid
    //alloc_N_distr(&(pt->griglia_gamma_Ne_log_IC),pt->emitters.gamma_grid_size);
    //alloc_N_distr(&(pt->Ne_IC),pt->emitters.gamma_grid_size);
    alloc_N_distr(&(pt->emitters.Integrand_over_gamma_grid),pt->emitters.gamma_grid_size);

}

void build_Q_inj_e_second(struct blob *pt) {
    //printf("build_Q_inj_e_second Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.Q_inj_e_second),pt->emitters.gamma_grid_size);
}

void build_Np(struct blob *pt)
{
    //printf("build_Np Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.griglia_gamma_Np_log), pt->emitters.gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->emitters.griglia_gamma_Np_log,pt->emitters.gmin_griglia, pt->emitters.gmax_griglia);
    //printf("build_Np Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.Np), pt->emitters.gamma_grid_size);
}

void build_Np_jetset(struct blob *pt) {
    //printf("build_Np_jetset Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.griglia_gamma_jetset_Np_log), pt->emitters.gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->emitters.griglia_gamma_jetset_Np_log,pt->emitters.gmin_griglia, pt->emitters.gmax_griglia);
    //printf("build_Np_jetset Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.Np_jetset), pt->emitters.gamma_grid_size);
    alloc_N_distr(&(pt->emitters.Ne_jetset),pt->emitters.gamma_grid_size);
}

void build_Ne_jetset(struct blob *pt) {
    //printf("build_Ne_jetset Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.griglia_gamma_jetset_Ne_log),pt->emitters.gamma_grid_size);
    Genera_griglia_gamma_N_log(pt, pt->emitters.griglia_gamma_jetset_Ne_log, pt->emitters.gmin_griglia, pt->emitters.gmax_griglia);
    //printf("build_Ne_jetset Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.Ne_jetset),pt->emitters.gamma_grid_size);
}


void Fill_Ne_IC(struct blob *pt, double g_min_IC, int stat_frame, double * Ne_IC, double * griglia_gamma_Ne_log_IC) {
    unsigned int i,i_start;
    //double gmin_grid;
    i_start=0;
    while (pt->emitters.griglia_gamma_Ne_log[i_start] < g_min_IC && i_start < pt->emitters.gamma_grid_size) {
        i_start++;
    }
    if (i_start % 2 != 0) {
        i_start = max(0,i_start-1);
    }
    g_min_IC=pt->emitters.griglia_gamma_Ne_log[i_start];
    
    if (pt->core.verbose>1) {
        printf("Set array per Ne IC\n");
        printf("elements number is pt->emitters.gamma_grid_size=%d\n", pt->emitters.gamma_grid_size);
    }
    //printf("Set array per Ne %s \n",pt->core.DISTR);
    
    if (strcmp(pt->core.PARTICLE, "protons") == 0) {
        
        if(pt->core.IC_adaptive_e_binning ==1){
            Genera_griglia_gamma_N_log(pt, griglia_gamma_Ne_log_IC,g_min_IC, pt->emitters.gmax_griglia_secondaries);
        }else{
            Genera_griglia_gamma_N_log(pt, griglia_gamma_Ne_log_IC,pt->emitters.gmin_griglia_secondaries, pt->emitters.gmax_griglia_secondaries);
        }
        
    }
    else{
        if(pt->core.IC_adaptive_e_binning ==1){
            //gmin_grid=max(pt->emitters.gmin_griglia,g_min_IC/10);
            Genera_griglia_gamma_N_log(pt, griglia_gamma_Ne_log_IC,g_min_IC, pt->emitters.gmax_griglia);
        }else{
            Genera_griglia_gamma_N_log(pt, griglia_gamma_Ne_log_IC,pt->emitters.gmin_griglia, pt->emitters.gmax_griglia);
        }
    }
    SetDistr(pt);
    for (i = 0; i < pt->emitters.gamma_grid_size; i++) {
        if(pt->core.IC_adaptive_e_binning ==1){
            if (griglia_gamma_Ne_log_IC[i]>=g_min_IC){
                Ne_IC[i] = N_distr_interp(pt->emitters.gamma_grid_size,
                                    griglia_gamma_Ne_log_IC[i],
                                    pt->emitters.griglia_gamma_Ne_log,
                                    pt->emitters.Ne);
            }else{
                Ne_IC[i]=0;
                }
        }else{
            Ne_IC[i] = pt->emitters.Ne[i]; 
        }                 
        if (stat_frame==1){
            //the delta^2 in Ne_stat is also correct because we use electron density
            //so the relativistic invariant is
            //N/(V*gamma^2)=N'/(V'gamma'2^)
            Ne_IC[i]*=pt->core.beam_obj*pt->core.beam_obj;

            //This transformation is correct
            //the grid is shifted by a factor of delta, hence the integration
            //boundaries are properly updated but the value of N[i] is still the
            //value of N(gamma') as in the formula 6.133 in Dermer&Menon
            griglia_gamma_Ne_log_IC[i]*=pt->core.beam_obj;
        }
    }
}



void build_Ne_custom(struct blob *pt,  unsigned int size) {
    pt->emitters.gamma_custom_grid_size=size;
    if (pt->core.verbose>1) {
        printf("Set array for Ne for from_array mode \n");
        printf("elements number is pt->emitters.gamma_grid_size=%d\n", pt->emitters.gamma_grid_size);
    }
    //printf("build_Ne_custom Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.gamma_e_custom),size);
    alloc_N_distr(&(pt->emitters.Ne_custom),size);

}

void build_Np_custom(struct blob *pt,  unsigned int size) {
    pt->emitters.gamma_custom_grid_size=size;
    if (pt->core.verbose>1) {
        printf("Set array for Np for from_array mode \n");
        printf("elements number is pt->emitters.gamma_grid_size=%d\n", pt->emitters.gamma_grid_size);
    }
    //printf("build_Np_custom Set array per Ne %s \n",pt->core.DISTR);
    alloc_N_distr(&(pt->emitters.gamma_p_custom),size);
    alloc_N_distr(&(pt->emitters.Np_custom),size);

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
    Fill_N(pt, pt->emitters.griglia_gamma_Ne_log, pt->emitters.Ne);


	//This flag is set to 1 to know that
	pt->emitters.Distr_e_done = 1;

    pt->emitters.N_0e = pt->emitters.N_0;
    //printf("==> N_tot\n");
    pt->emitters.N_e  = N_tot(pt, N_distr_integranda);
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

    pt->emitters.gmin_secondaries=pt->emitters.gmin;
    pt->emitters.gmax_secondaries=pt->emitters.gmax*mp_by_me;
    setNgrid(pt);
    build_Np(pt);
    SetDistr(pt);
    if (pt->core.verbose>1) {
        printf("********** protons ***********\n");
        printf("set array for Np\n");
        printf("elements number is pt->emitters.gamma_grid_size=%d\n", pt->emitters.gamma_grid_size);
        printf("TIPO_DISTR %d\n", pt->emitters.TIPO_DISTR);
    }
    Fill_N(pt, pt->emitters.griglia_gamma_Np_log, pt->emitters.Np);
    //
    //This flag si set to 1 to know that
    //N(gamma) has been properly initialized and filled
    //printf("-->\n" );
    //printf("--> N0e %e N0 %e N0p %e\n", pt->emitters.N_0e, pt->emitters.N_0, pt->emitters.N_0p);
    
    pt->emitters.Distr_p_done = 1;
    pt->emitters.N_0p = pt->emitters.N_0;
    
    //printf("--> N0e %e N0 %e N0p %e\n", pt->emitters.N_0e, pt->emitters.N_0, pt->emitters.N_0p);
    
    pt->emitters.N_p = N_tot(pt, N_distr_integranda);
    //name = "distr-p.dat";
    //Scrivi_N_file(pt, name, pt->emitters.griglia_gamma_Np_log, pt->emitters.Np);

    // Secondaries e- from pp

    //Set N to e- from pp
    sprintf(pt->core.PARTICLE, "secondaries_el");
    setNgrid(pt);
    build_Ne_secondaries(pt);
    build_Q_inj_e_second(pt);
    SetDistr(pt);
    pt->PP_gamma.pp_racc_elec=rate_electrons_pp(pt, pt->emitters.griglia_gamma_Ne_log[0],1);
    Fill_N(pt, pt->emitters.griglia_gamma_Ne_log, pt->emitters.Q_inj_e_second);
    CoolingEquilibrium(pt,pt->emitters.T_esc_e_second);
    //Filling Ne_jetset with secondaries
    unsigned int i;
    for (i = 0; i < pt->emitters.gamma_grid_size; i++) {
        pt->emitters.Ne_jetset[i]=pt->emitters.Ne[i];
    }
    pt->emitters.Distr_e_done = 1;
    pt->emitters.N_0e = pt->emitters.N_0;
    

    //printf("--> N0e %e N0 %e N0p %e\n", pt->emitters.N_0e, pt->emitters.N_0, pt->emitters.N_0p);
    //printf("-->\n");
    pt->emitters.N_e_pp = N_tot(pt, N_distr_integranda);

    //if (pt->core.verbose > 1)
    //{
    //    printf("****** secondary leptons *****\n");
    //    printf("set array for secondary Ne\n");
    //    printf("elements number is pt->emitters.gamma_grid_size=%d\n", pt->emitters.gamma_grid_size);
    //    printf("N_e_pp =%e\n", pt->emitters.N_e_pp);
    //}

    //name = "distr-e-from-pp.dat";
    //Scrivi_N_file(pt, name, pt->emitters.griglia_gamma_Ne_log, pt->emitters.Ne);
    
    //set back pt->emitters.N_0 to the proton value and particle name
    pt->emitters.N_0 = pt->emitters.N_0p;
    sprintf(pt->core.PARTICLE, "protons");
    SetDistr(pt);
}





//========================================
// Trova il gmax da N[i]>o
//========================================

double Find_gmax(struct blob *pt, double *N, double *g) {
	unsigned int i;
    double gmax;
    gmax = g[0];

    for (i = 0; i < pt->emitters.gamma_grid_size; i++) {
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

    pt->emitters.N_0 = 1.0;
    //=========================================
    // interpolate custom Ne/p
    //=========================================
    if (pt->emitters.TIPO_DISTR == 0)
    {   
        if (strcmp(pt->core.PARTICLE, "protons") == 0){
            for (i = 0; i < pt->emitters.gamma_grid_size; i++)
            
            {
                N[i] = N_distr_interp(pt->emitters.gamma_custom_grid_size,
                                    griglia_gamma_N_log[i],
                                    pt->emitters.gamma_p_custom,
                                    pt->emitters.Np_custom);
            }
        }else{
            for (i = 0; i < pt->emitters.gamma_grid_size; i++)
            
            {
                N[i] = N_distr_interp(pt->emitters.gamma_custom_grid_size,
                                    griglia_gamma_N_log[i],
                                    pt->emitters.gamma_e_custom,
                                    pt->emitters.Ne_custom);
            }
        }

    }
    else if (pt->emitters.TIPO_DISTR==10){
        if (strcmp(pt->core.PARTICLE, "protons") == 0){
            for (i = 0; i < pt->emitters.gamma_grid_size; i++)
            
            {
                N[i] = pt->emitters.Np_jetset[i];
            }
        }else{
            for (i = 0; i < pt->emitters.gamma_grid_size; i++)
            
            {
                N[i] = pt->emitters.Ne_jetset[i];
            }
        }

    }

    //=========================================
    // fill defined Ne/p
    //=========================================
    else if (pt->emitters.TIPO_DISTR != -1){

        //Normalization
        
        
        if (pt->emitters.Norm_distr == 1 && pt->emitters.TIPO_DISTR != -1)
        {
            pf_norm = &N_distr_integranda;
            pt->emitters.N_0 = integrale_trap_log_struct(pf_norm, pt, pt->emitters.gmin, pt->emitters.gmax, 10000);
        }

        for (i = 0; i < pt->emitters.gamma_grid_size; i++)
        {
            N[i] = N_distr(pt, griglia_gamma_N_log[i]);
        }
    }

    //if distr is e- from pp te
    //the distribution is filled with the injection
    //by the function N_distr
    else if (pt->emitters.TIPO_DISTR == -1){
        for (i = 0; i < pt->emitters.gamma_grid_size; i++)
        {
            N[i] = N_distr(pt, griglia_gamma_N_log[i]);
        }
    }
    else {
        printf("TIPO_DISTR set to wrong value: %d\n",pt->emitters.TIPO_DISTR);
        exit(1);
    }

    //pt->emitters.Distr_e_done = 1;
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
    if (Gamma >= pt_N->emitters.gmin_secondaries && Gamma <= pt_N->emitters.gmax_secondaries && pt_N->emitters.TIPO_DISTR == -1) {
        
        a= vluce_cm * pt_N->PP_gamma.NH_pp * MEC2_TeV * bn_to_cm2 * rate_electrons_pp(pt_N, Gamma,-1);
    }else{

        a= N_distr_integranda(pt_N,Gamma)*pt_N->emitters.N/pt_N->emitters.N_0;

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
                                pt->emitters.gmin,
                                pt->emitters.gmax,
                                10000);

    //if the distr is not secondaries or interpolated
    if (pt->emitters.TIPO_DISTR > 0)
    {
        a = a * pt->emitters.N / pt->emitters.N_0;
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

    if (Gamma >= pt_N->emitters.gmin && Gamma <= pt_N->emitters.gmax){

        //Secondaris e Distribution has not analytical expression
        //it is taken from the N array, throug log-lin interpolation
        if (  pt_N->emitters.TIPO_DISTR == -1) {
            a= N_distr_interp(pt_N->emitters.gamma_grid_size,
                              Gamma,
                              pt_N->emitters.griglia_gamma_Ne_log,
                              pt_N->emitters.Ne);
        }

        if (  pt_N->emitters.TIPO_DISTR == 0) {
            a= N_distr_interp(pt_N->emitters.gamma_custom_grid_size,
                                Gamma,
                                pt_N->emitters.gamma_e_custom,
                                pt_N->emitters.Ne_custom);
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
    //printf("G=%e G_file=%e\n",pt->emitters.griglia_gamma_Ne_log[i],G_File[count]);
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

    if (strcmp(pt->core.PARTICLE, "secondaries_el") == 0)
    {
        pt->emitters.TIPO_DISTR = -1;
    }
    else
    {
        if (strcmp(pt->core.DISTR, "from_array") == 0)
        {
            pt->emitters.TIPO_DISTR = 0;
        }

        if (strcmp(pt->core.DISTR, "jetset") == 0)
        {
            pt->emitters.TIPO_DISTR = 10;
        }

    }
    //if (pt->core.verbose) {
     //printf("tipo di distribuzione %d\n",pt->emitters.TIPO_DISTR);
    //}
}
//=========================================================================================
