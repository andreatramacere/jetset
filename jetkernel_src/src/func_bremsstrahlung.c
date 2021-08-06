#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"
/**
 * \file bremsstrahlung.c
 * \author Andrea Tramacere
 * \date 22-12-2020
 * \brief bremsstrahlung emission
 *
 */

//=======================================================================
//   bremsstrahlung following Baring et al 1999
//   https://iopscience.iop.org/article/10.1086/306829/pdf
//=======================================================================

double j_nu_bremss_ee(struct blob * pt, double nu_out) {
    //double integr;
    unsigned int ID_gamma;
    for (ID_gamma = 0; ID_gamma < pt->gamma_grid_size ; ID_gamma++){
            pt->Integrand_over_gamma_grid[ID_gamma] =b_ee_sigma(pt->griglia_gamma_Ne_log[ID_gamma],nu_out*HPLANCK*one_by_MEC2)*pt->Ne[ID_gamma];
        }
    return vluce_cm*one_by_four_pi*one_by_MEC2*HPLANCK*HPLANCK*nu_out*integr_simp_grid_equilog(pt->griglia_gamma_Ne_log, pt->Integrand_over_gamma_grid, pt->gamma_grid_size);
}

double j_nu_bremss_ep(struct blob *pt, double nu_out) {
    //double integr;
    unsigned int ID_gamma;
    for (ID_gamma = 0; ID_gamma < pt->gamma_grid_size ; ID_gamma++){
            pt->Integrand_over_gamma_grid[ID_gamma] =b_ep_sigma(pt->griglia_gamma_Ne_log[ID_gamma],nu_out*HPLANCK*one_by_MEC2)*pt->Ne[ID_gamma];
        }
    return vluce_cm*one_by_four_pi*one_by_MEC2*HPLANCK*HPLANCK*nu_out*integr_simp_grid_equilog(pt->griglia_gamma_Ne_log, pt->Integrand_over_gamma_grid, pt->gamma_grid_size);
}

double b_ep_sigma(double gamma_e, double epsilon_gamma){
    ////ep valid only for Eout>10MeV
    //if  (epsilon_gamma>10/MEC2_MeV){
    return bremss_sigma_1(gamma_e,epsilon_gamma);
    //}
    //return 0.;
}

double b_ee_sigma(double gamma_e, double epsilon_gamma){
    //cm2
    if  (gamma_e*MEC2_MeV<=brem_ee_th){
        return b_ee_sigma_non_rel(gamma_e,epsilon_gamma);
    }else{
        return b_ee_sigma_rel(gamma_e,epsilon_gamma);
    } 
    return 0.;
}

double b_ee_sigma_rel(double gamma_e, double epsilon_gamma){
    //Eq A1
    //For gamma_e*MEC2_MeV>=2
    //https://iopscience.iop.org/article/10.1086/306829/pdf
    //double res;

    return b_ee_A_term(gamma_e,epsilon_gamma)*(bremss_sigma_1(gamma_e,epsilon_gamma)+bremss_sigma_2(gamma_e,epsilon_gamma));
}

double b_ee_sigma_non_rel(double gamma_e, double epsilon_gamma){
    //Eq A5
    //https://iopscience.iop.org/article/10.1086/306829/pdf
    //For gamma_e*MEC2_MeV<2
    double x,a,res;
    x=4*epsilon_gamma/(gamma_e*gamma_e -1);
    a=0.25*(gamma_e*gamma_e-1);
    res=0;
    if ((epsilon_gamma>0) && (epsilon_gamma<a)){
     res= 4.*e_raggio*e_raggio/(15*one_by_alpha*epsilon_gamma)*ee_brem_F(gamma_e,x);
    }
    return res;
}

double ee_brem_F(double gamma_e, double x){
    //Eq A6
    //https://iopscience.iop.org/article/10.1086/306829/pdf
    double beta_e,f_B,f_C,a1,a2,a3;
    beta_e= sqrt(1 - (gamma_e*gamma_e));
    
    f_B = 1 + 0.5*((gamma_e*gamma_e) - 1);
    f_C = 10 * x * gamma_e * beta_e * (2 + gamma_e * beta_e);
    f_C *= 1.0/(1 + (x*x) * ((gamma_e*gamma_e) - 1));

    a1 = (17 -( 3 * x*x / ((2 - x) *(2 - x))) - f_C) * sqrt(1 - x);
    a2 = 12 * (2 - x) - (7*x*x /(2 - x)) - (3 * x * x * x * x /((2 - x)*(2 - x)*(2 - x)));
    a3 = log((1 + sqrt(1 - x)) /sqrt(x));

    return f_B*a1 +(a2*a3);
}


double b_ee_A_term(double gamma_e, double epsilon_gamma){
    //Eq A4
    //https://iopscience.iop.org/article/10.1086/306829/pdf
    return 1 -((8/3)*pow((gamma_e-1),0.2)/(gamma_e+1)*pow((epsilon_gamma/gamma_e),1/3));
}

double bremss_sigma_1(double gamma_e, double epsilon_gamma){
    //Eq A2 
    //https://iopscience.iop.org/article/10.1086/306829/pdf
    double a,a1,a2,x,res;
    if((gamma_e-epsilon_gamma)<0){
        res =0.;
    
    }else{
      x=epsilon_gamma/gamma_e;
        a =4.*e_raggio*e_raggio/(one_by_alpha*epsilon_gamma);
        a1=1 + ((1./3.)-x)*(1-x);
        a2=log(2*gamma_e*((1/x) -1 ))-0.5;

        res= a*a1*a2;
    }
    
    return res;
}

double bremss_sigma_2(double gamma_e, double epsilon_gamma){
    //Eq A3 
    //https://iopscience.iop.org/article/10.1086/306829/pdf
    double res,a,a1,a2,a3,one_by_eps;
    a2=0;  
    one_by_eps=1.0/epsilon_gamma;

    if (epsilon_gamma<=0.5){
        a=16*(1-epsilon_gamma+(epsilon_gamma*epsilon_gamma))*log(gamma_e*one_by_eps);
        a1=(-1.0/(one_by_eps*one_by_eps));
        a2 += 3.0/(one_by_eps) -4 +(4*epsilon_gamma) -(8*epsilon_gamma*epsilon_gamma);
        a3 = -2*(1-2*epsilon_gamma)*log(1-(2*epsilon_gamma));
        a3 *= (0.5*one_by_eps*one_by_eps*one_by_eps) -0.5*one_by_eps +3.0*one_by_eps -2 +4*one_by_eps;
        res=a+a1+a2+a3; 
    }else{
        a=2*one_by_eps;
        a1=(4-one_by_eps+(0.25*one_by_eps))*log(2*gamma_e);
        a2=-2+ (2*one_by_eps) -(0.625*one_by_eps*one_by_eps);
        res = a*(a1+a2);
    }

    return res*e_raggio*e_raggio/(one_by_alpha*3*gamma_e);

}


