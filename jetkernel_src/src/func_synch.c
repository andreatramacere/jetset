//=================================================================
//
//  FUNZIONI SPETTRO SINCROTRONE
//=================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"

/**
 * \file funzioni_sincrotrone.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief spettro di sincrotrone
 *
 */





//=========================================================================================
// Sync F(X) log-log interpolation
//=========================================================================================
double F_K_53(struct blob * pt, double x){
    return log_log_interp(log10(x), pt->log_F_Sync_x, pt->log_x_Bessel_min, pt->log_x_Bessel_max, pt->log_F_Sync_y,static_bess_table_size,0  );
}
double F_K_23(struct blob * pt, double x){
    return log_log_interp(log10(x), pt->log_G_Sync_x, pt->log_x_Bessel_min, pt->log_x_Bessel_max, pt->log_G_Sync_y,static_bess_table_size,0  );
}

double F_K_ave(struct blob *pt, double x){
    return log_log_interp(log10(x), pt->log_F_ave_Sync_x, pt->log_x_ave_Bessel_min, pt->log_x_ave_Bessel_max, pt->log_F_ave_Sync_y,static_bess_table_size,0  );

}
//=========================================================================================


//=========================================================================================
// j_nu Sync integrands
//=========================================================================================
double F_int_fix(struct blob * pt,unsigned int  ID, double nu_sync){
    //PITCH ANGLE FIXED
    double a, y,g;
    g=pt->griglia_gamma_Ne_log[ID];
    y=(nu_sync/(g*g))*pt->C2_Sync_K53;
    a=F_K_53(pt, y);
    a*=pt->Ne[ID];
    return a;
}

double F_int_fix_parallel(struct blob * pt,unsigned int  ID, double nu_sync){
    //PITCH ANGLE FIXED
    double a, y,g;
    g=pt->griglia_gamma_Ne_log[ID];
    y=(nu_sync/(g*g))*pt->C2_Sync_K53;
    a=F_K_23(pt, y);
    a*=pt->Ne[ID];
    return a;
}


double F_int_ave(struct blob * pt,unsigned int  ID, double nu_sync){
    //PITCH ANGLE AVE
    //The Astrophysical Journal, 334:L5-L8,1988 November 1
    double a, y,g;
    g=pt->griglia_gamma_Ne_log[ID];
    y=nu_sync/(g*g)*pt->C2_Sync_K_AVE;
    a=F_K_ave(pt, y);
    a*=pt->Ne[ID];
    return a;
}
//=========================================================================================





//=========================================================================================
//    integrand for alfa_nu_Sync
//=========================================================================================
double Sync_self_abs_int(struct blob *pt,unsigned int  ID, double nu_sync){
    double a,g, y, x1, x2, y1, y2, delta;
    
    

    //----------derivative-----------------
    //if i=0 forward
    //if i>0 back
    //-----------------------------------
    if(ID==0){
        x1=pt->griglia_gamma_Ne_log[ID];
        x2=pt->griglia_gamma_Ne_log[ID+1];
        y1=pt->Ne[ID]/(x1*x1);
        y2=pt->Ne[ID+1]/(x2*x2);
    }
    else{
        x1=pt->griglia_gamma_Ne_log[ID-1];
        x2=pt->griglia_gamma_Ne_log[ID];
        y1=pt->Ne[ID-1]/(x1*x1);
        y2=pt->Ne[ID]/(x2*x2);
    }
    delta=x2-x1;
    a=(y2-y1)/delta;
    g=pt->griglia_gamma_Ne_log[ID];
    if (pt->Sync_kernel==0){
    	y=(nu_sync/(g*g))*pt->C2_Sync_K53;
    	a*=(g*g)*F_K_53(pt, y);
    }
    else{
    	y=nu_sync/(g*g)*pt->C2_Sync_K_AVE;
    	a*=(g*g)*F_K_ave(pt, y);
    }
    //This is fixing numerical instabilities 
    if( a>0.0){
        a=0.0;
    }
    return a;   
}
//=========================================================================================





//=========================================================================================
// Radiative transfer solution for sefl abs
// see Kataoka Thesis, page 299
// tau_nu in I_nu=alfa_nu*R, so we multiply by 0.5
double solve_S_nu_Sync(struct blob * pt, unsigned int  NU_INT){
	double S_nu;
    // pt->I_nu_Sync[NU_INT] = 0.0;


	// if (pt->do_Sync == 2) {
	// 	tau_nu = 2 * pt->R_sync_self_abs * pt->alfa_Sync[NU_INT];
	// 	if (tau_nu > 1e-4) {
	// 		pt->I_nu_Sync[NU_INT] =
	// 				(pt->j_Sync[NU_INT] / pt->alfa_Sync[NU_INT])*
	// 				(1 - exp(-tau_nu*0.5));

    //     } else {
	// 		pt->I_nu_Sync[NU_INT] =
	// 				(pt->j_Sync[NU_INT] / pt->alfa_Sync[NU_INT])*
	// 				( tau_nu*0.5 - (1.0 / 4.0) * tau_nu * tau_nu*0.5*0.5);
			
	// 	}
	// }

	// //==========================
	// //Radiative solution for no self abs
	// //limit of S_nu,alfa->0=(4/3)*R
	// if (pt->do_Sync == 1) {
	// 	pt->I_nu_Sync[NU_INT]=pt->j_Sync[NU_INT] * pt->R_sync;
	// }
	// if (pt->verbose>1) {
	// 	printf("#-> nu=%e j=%e alfa=%e tau_nu=%e  I_nu=%e\n", pt->nu_Sync[NU_INT], pt->j_Sync[NU_INT],
	// 			pt->alfa_Sync[NU_INT], tau_nu, pt->I_nu_Sync[NU_INT]);
	// }
    
    S_nu = eval_S_nu_Sync(pt, pt->j_Sync[NU_INT], pt->alfa_Sync[NU_INT]);
    pt->I_nu_Sync[NU_INT]=I_nu_to_L_nu_blob(S_nu,pt->Surf_region)/(16*pi*pt->R_sync_n_photons*pt->R_sync_n_photons); 
    return S_nu;
}

double eval_S_nu_Sync(struct blob *pt, double j_Sync, double alfa_Sync)
{
    double S_nu, tau_nu;
    S_nu=0;
    if (pt->do_Sync == 2)
    {
        tau_nu = 2 * pt->R_sync_self_abs *  alfa_Sync;
        if (tau_nu > 1e-4)
        {
          
            S_nu = ( j_Sync /  alfa_Sync) *
                (1.0 - (2 / (tau_nu * tau_nu)) * (1 - exp(-tau_nu) * (tau_nu + 1)));
        }
        else
        {
           
            S_nu = ( j_Sync /  alfa_Sync) *
                ((2.0 / 3.0) * tau_nu - (1.0 / 4.0) * tau_nu * tau_nu);
        }
    }

    //==========================
    
    if (pt->do_Sync == 1)
    {
       
        S_nu =  j_Sync * pt->R_sync;
    }
    
    return S_nu;
}



void set_R_Sync(struct blob * pt){
    double R_sync_Shell;
    if (strcmp(pt->GEOMETRY, "spherical") == 0) {
            pt->R_sync_self_abs = pt->R;
            //Radiative solution for no self abs
            //limit of S_nu,alfa->0=(4/3)*R
            pt->R_sync = pt->R* four_by_three;
            pt->n_sync_corr_factor=0.75;
            pt->R_sync_n_photons=pt->R;

        }

        else if (strcmp(pt->GEOMETRY, "spherical_shell") == 0) {    
            R_sync_Shell=(1 - (1-pt->h_sh)*(1-pt->h_sh)*(1-pt->h_sh))*pt->R_sh;
            pt->R_sync_self_abs = R_sync_Shell;
            pt->R_sync = R_sync_Shell*(four_by_three);
            pt->n_sync_corr_factor=1.0;
            //NOTE: This is a crude approximation, should be improved 
            pt->R_sync_n_photons=pt->R_sh;
        }
        else {
            printf("GEOMETRY variable set to wrong value, possible spherical or spherical_shell \n");
            exit(0);
        }

    }


//=========================================================================================






//=========================================================================================
//  Synchrotron emissivity j_nu_Sync
//=========================================================================================
double j_nu_Sync(struct blob * f, double nu_sync){
    double a;
    double (*pf_fint) (struct blob * ,unsigned int  ID, double nu_sync);
    /*** segli in base al kernel ***/
    if (f->Sync_kernel==0){
		pf_fint=&F_int_fix;
		a=integrale_Sync(pf_fint, f,  nu_sync);
		return a*f->C1_Sync_K53;
    }
    else {
    	pf_fint=&F_int_ave;
    	a=integrale_Sync(pf_fint, f, nu_sync);
    	return a*f->C1_Sync_K_AVE;
    }
}
//=========================================================================================

double eval_Sync_polarization(struct blob * f, double nu_sync){
    double p_num,p_den,pol;
    pol=0;
    double (*pf_fint) (struct blob * ,unsigned int  ID, double nu_sync);
    //EQ 6.37 R&L, for integration over N(gamma)
    //Since The integral is additive, 6.37 holds
    //p_num=P_ort-P_parallel=> integral [F(x) + G(X)] - integral [F(X) - G(X)] =  2*integral [G(X)]
    //p_den=P_ort+P_parallel=> integral [F(x) + G(X)] + integral [F(X) - G(X)] =  2*integral [F(X)]
    //pol=integral [G(X)]/integral [F(X)]
  
        pf_fint=&F_int_fix_parallel;
        p_num=integrale_Sync(pf_fint, f,  nu_sync);
        pf_fint=&F_int_fix;
        p_den=integrale_Sync(pf_fint, f,  nu_sync);
        pol=(p_num)/(p_den);
  
    return pol;
}





//=========================================================================================
// Synch self abs alfa_nu_Sync
//=========================================================================================
double alfa_nu_Sync(struct blob * f, double nu_sync){
    double a;
    double (*pf_fint1) (struct blob * , unsigned int ID, double nu_sync);
    pf_fint1=&Sync_self_abs_int;
    a=integrale_Sync(pf_fint1, f,nu_sync);
    return a*f->C3_Sync_K53*(f->B)/(nu_sync*nu_sync);
}
//=========================================================================================





//=========================================================================================
// Sync INTEGRATION WITH SIMPSON AND GRIGLIA EQUI-LOG
//=========================================================================================
double integrale_Sync(double (*pf) (struct blob *, unsigned int  ID, double nu_sync), struct blob * pt, double nu_sync ) {

    unsigned int  ID;
    double *Integrand_over_gamma_grid;
    double integral;
    integral =0;
    Integrand_over_gamma_grid = (double *) calloc(pt->gamma_grid_size, sizeof (double));
    //double test;
    for (ID = 0; ID < pt->gamma_grid_size ; ID++){
        Integrand_over_gamma_grid[ID] =pf(pt,ID, nu_sync);
    }
    integral= integr_simp_grid_equilog(pt->griglia_gamma_Ne_log, Integrand_over_gamma_grid, pt->gamma_grid_size);
    free(Integrand_over_gamma_grid);
    return integral;
}
//=========================================================================================




//=========================================================================================
double Sync_tcool(struct blob * pt, double g){
    	return g/Sync_cool(pt,g);
}
//=========================================================================================




//=========================================================================================
double Sync_cool(struct blob * pt, double g){
	double beta_gamma,c;
	beta_gamma=eval_beta_gamma(g);
    if (pt->Sync_kernel==0){
    	c=2.0*beta_gamma*beta_gamma*pt->UB*g*g*pt->sin_psi*pt->sin_psi;

    }
    else{
    	c=four_by_three*beta_gamma*beta_gamma*pt->UB*g*g;

    }

    return pt->COST_Sync_COOLING*c;
}
//=========================================================================================
