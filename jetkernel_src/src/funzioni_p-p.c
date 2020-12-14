//=================================================================
//
//  FUNZIONI PER P-P Process
//=================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file funzioni_pp.c
 * \author Andrea Tramacere
 * \date 04-05-2010
 * \brief spettri dal processo pp
 *
 */

//=========================================================================================
// PP inelastic cross section
//
//=========================================================================================

double sigma_pp_inel(double Ep_TeV) {
    //From Eq.79 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //L=ln(E_p/1TeV)
    //
    double a,b,c;
    double k, L;
    if (Ep_TeV<E_th_pp){
        c=0;
    }
    else{
        k = E_th_pp / Ep_TeV;
        L = log(1/k);
        a=30.7-0.96*L+0.18*L*L;
        b=(1.0 - pow(k,1.9));
        c= a*b*b*b;
    }
    return c;
}
//=========================================================================================

unsigned int E_min_p_grid_even(struct spettro *pt, double * gamma_p_grid, double E_start_TeV, unsigned int i_start, unsigned int gamma_p_grid_size ){
    //return the startin grid even index
    double  gamma_p_min;
    unsigned i_start_out;
    i_start_out=i_start;
    if (E_start_TeV<pt->E_th_pp_delta_approx){
        E_start_TeV=pt->E_th_pp_delta_approx;
    }
    gamma_p_min =  E_start_TeV/ MPC2_TeV;
    while ((gamma_p_grid[i_start_out] < gamma_p_min) && (i_start_out < gamma_p_grid_size)) {
        i_start_out++;
    }
    //!!!!!!!
    //i_start deve essere pari !!!
    //!!!!!!!
    if (i_start_out % 2 != 0) {
        i_start_out++;
    }

    if (i_start_out>= gamma_p_grid_size) {
        i_start_out= gamma_p_grid_size;
    }
    if (E_start_TeV)
    return i_start_out;
}

double  check_pp_kernel(double res,struct spettro *pt,double E_p_TeV, double x ){
    double res_out;
    res_out=res;
    if (E_p_TeV<pt->E_th_pp_delta_approx){
        res_out=0.0;

    }
    if (x<pt->E_pp_x_delta_approx){
        res_out=0.0;

    }
    return res_out;
}


//=========================================================================================
// PP -> e- production
//
//=========================================================================================
double rate_electrons_pp(struct spettro *pt, double Gamma_e) {
    //From Eq.71 and 78 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //e- from mu->e +mu_nu +mu_e 
    double  Ee_TeV, a1, a2, Emin_pi,Emax_pi;
    unsigned int i_start;
    double (*pf_K) (double gamma_p, double nu_pp, struct spettro * pt);
    double (*pf_K_delta) (struct spettro * pt, double x);
    double (*pf_E_min) (double gamma_p, struct spettro * pt);
    double (*pf_E_max) (struct spettro * pt);
    double res;

    pf_K= &pp_electrons_kernel;
    pf_K_delta = &pp_electron_kernel_delta;
    pf_E_min = &E_min_e_pp;
    pf_E_max = &E_max_e_pp;
    //!!!! MPICC2_TeV RO MPI0C2_TeV
    pt->MPI_kernel_delta=MPICC2_TeV;
    pt->MPI_kernel_delta_Emin=MPICC2_TeV;


    if (pt->set_pp_racc_elec == 0) {
        pt->set_pp_racc_elec = 1;
        Ee_TeV =  pt->E_th_pp_delta_approx;
        pt->E_out_e_TeV_pp=Ee_TeV;
        //printf("==>1 pp_racc_elec=%e\n",pt->pp_racc_elec);
        i_start = E_min_p_grid_even(pt,pt->griglia_gamma_Np_log,Ee_TeV, 0, pt->gamma_grid_size );
        //Eq. 71
        a2 = integrale_pp_second_high_en_rate(pf_K, Ee_TeV, pt, i_start);
        //Eq. 78
        a1= integrale_pp_second_low_en_rate(pf_K_delta,pf_E_min,pf_E_max,Ee_TeV,pt);
        //printf("==>a1=%e, pt->a2=%e\n",a1,a2);
        pt->pp_racc_elec = a2 / a1;
        //printf("==>2 pp_racc_elec=%e\n",pt->pp_racc_elec);
        //exit(0);
    }
    Ee_TeV = Gamma_e * MEC2_TeV;
    pt->E_out_e_TeV_pp=Ee_TeV;
    //Eq. 71
    if (Ee_TeV > pt->E_th_pp_delta_approx) {
        i_start = E_min_p_grid_even(pt,pt->griglia_gamma_Np_log,E_min_e_pp(Ee_TeV,pt), 0, pt->gamma_grid_size );
        //i_start=0;
        res= integrale_pp_second_high_en_rate(pf_K, Ee_TeV, pt, i_start);
    } else {        
        //Eq. 78
        res= integrale_pp_second_low_en_rate(pf_K_delta,pf_E_min,pf_E_max,Ee_TeV,pt);
        
    }
    //printf("==>Gamma_e=%e, pt->E_out_e_TeV_pp=%e pt->pp_racc_elec=%e\n",Gamma_e,pt->E_out_e_TeV_pp,pt->pp_racc_elec);
    return res;
}

double E_min_e_pp(double E_e, struct spettro *pt){
    double psida;
    psida=  1.0 ;
    return (E_e/psida)+(pt->MPI_kernel_delta_Emin * pt->MPI_kernel_delta_Emin)*psida/ (4 * E_e);
}
double E_max_e_pp(struct spettro *pt){
    return (pt->gmax*MPC2_TeV - MPC2_TeV);
}

double f_mu_2_pp(double x, double r){
    double a,b;
    a=g_mu_pp(x,r)*theta_heaviside(x-r);
    b=(h_mu_1_pp(x,r) + h_mu_2_pp(x,r))*theta_heaviside(r-x);
    //printf("==> g=%e,h1=%e,h2=%e (x-r)=%e, (r-x)=%e r=%e x=%e a=%e b=%e\n",g_mu_pp(x,r),h_mu_1_pp(x,r),h_mu_2_pp(x,r),(x-r),(r-x),r,x,a,b);
    return a+b;
}


double g_mu_pp(double x, double r){
    double a,b;
    a=(3.-2.*r)/(9.*(1. - r)*(1 - r));
    b=(9.* x*x -6.*log(x) -4.0*x*x*x -5.0 );    
    return a*b;
}

double h_mu_1_pp(double x, double r){
    double a,b;
    a=(3.-2.*r)/(9.*(1. - r)*(1 - r));
    b=(9.* r*r -6.*log(r) -4.0*r*r*r -5.0 );    
    return a*b;
}

double h_mu_2_pp(double x, double r){
    double a,b;
    a=(1. + 2.*r)*(r -x)/(9.*r*r);
    b=(9.*(r+x) -4.*((r*r) +r*x +(x*x)));  
    return a*b;
}

double pp_electron_kernel_delta(struct spettro *pt,double E_pi) {
    double qe, Ep0_TeV, gamma_p,res,f;
    Ep0_TeV = E_pi / (Kpi) + MPC2_TeV;
    gamma_p = Ep0_TeV / MPC2_TeV;
    
    qe = pt->pp_racc_elec / (Kpi) * sigma_pp_inel(Ep0_TeV)*
        N_distr_interp(pt->gamma_grid_size, gamma_p, pt->griglia_gamma_Np_log, pt->Np);
    f=f_mu_2_pp(( pt->E_out_e_TeV_pp/E_pi),0.573);
    res =2.0*qe*f/E_pi;
    //printf("==>gamma_p=%e f=%e q=%e E_pi=%e pt->pp_racc_elec=%e, sigma_pp_inel(Ep0_TeV)=%e N_distr=%e\n",gamma_p,f,qe,E_pi,pt->pp_racc_elec,sigma_pp_inel(Ep0_TeV),N_distr_interp(pt->gamma_grid_size, gamma_p, pt->griglia_gamma_Np_log, pt->Np));
    //return 2.0 * qe / sqrt(E_pi * E_pi - pt->MPI_kernel_delta * pt->MPI_kernel_delta);
    //printf("==>gamma_p=%e E_out_e_TeV=%e, qe=%e, res=%e  sigma_pp_inel(Ep0_TeV)=%e, N_distr_interp=%e\n",gamma_p,pt->E_out_e_TeV_pp,qe,res,sigma_pp_inel(Ep0_TeV), N_distr_interp(pt->gamma_grid_size, gamma_p, pt->griglia_gamma_Np_log, pt->Np));
    return res;
}

double pp_electrons_kernel(double gamma_p, double E_out_TeV, struct spettro *pt) {
    double Ep_TeV,x,res;
    Ep_TeV = gamma_p*MPC2_TeV;
    x=E_out_TeV / Ep_TeV;
    ///printf("Ee_TeV=%e gamma_e=%e, F_e=%e\n",Ep_TeV,gamma_p,(F_electrons((E_out_TeV / Ep_TeV), Ep_TeV) / pt->griglia_gamma_Np_log[pt->i_griglia_gamma]));
    res= sigma_pp_inel(Ep_TeV) *
            pt->Np[pt->i_griglia_gamma] *
            F_electrons((x), Ep_TeV) / pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
    return check_pp_kernel(res,pt,Ep_TeV,x);
}

double F_electrons(double x, double Ep_TeV) {
    //From Eq.62 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //x=E_e/E_pi
    //L=ln(E_p/1TeV)

    double B, beta, k, lx, xb, L, a,res;
    if (x<=0.001){
       res=0.; 
    }
    else{
        lx = log(x);
        L = log(Ep_TeV);

        B = 1 / (69.5 + 2.65 * L + 0.3 * L * L);
        beta = 1.0 / pow((0.201 + 0.062 * L + 0.00042 * L * L), 0.25);
        k = (0.279 + 0.141 * L + 0.0172 * L * L) / (0.3 + (2.3 + L)*(2.3 + L));
        xb = pow(x, beta);

        a = (1 + k * lx * lx);
        //printf("%e %e %e %e %e %e\n",B,beta,x,xb,L,lx);
        res= B * (a * a * a) / (x * (1 + 0.3 / xb)) * (-lx)*(lx)*(lx)*(lx)*(lx);
    }
    return res;
}


//=========================================================================================
// PP -> nu_mu production
//
//=========================================================================================

double rate_neutrino_mu_1_pp(struct spettro *pt, double nu_nu_mu) {
    //From Eq.71 and 78 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //For pi->mu+nu_mu
    double  Emu_TeV, a1, a2, Emin_pi,Emax_pi;
    unsigned int i_start;
    double (*pf_K) (double gamma_p, double nu_mu_nu, struct spettro * pt);
    double (*pf_K_delta) (struct spettro * pt, double x);
    double (*pf_E_min) (double gamma_p, struct spettro * pt);
    double (*pf_E_max) (struct spettro * pt);
    pf_K = &pp_neturino_mu_1_kernel;
    pf_K_delta = &pp_neutrino_mu_1_kernel_delta;
    pf_E_min = &E_min_neutrino_mu_1_pp;
    pf_E_max = &E_max_neutrino_mu_1_pp;

    pt->MPI_kernel_delta=MPICC2_TeV;
    pt->MPI_kernel_delta_Emin=MPICC2_TeV;
    if (pt->set_pp_racc_nu_mu == 0) {
        pt->set_pp_racc_nu_mu = 1;
        Emu_TeV =  pt->E_th_pp_delta_approx;
       
        i_start = E_min_p_grid_even(pt,pt->griglia_gamma_Np_log,Emu_TeV, 0, pt->gamma_grid_size );
        
        //Eq. 71
        a2 = integrale_pp_second_high_en_rate(pf_K,Emu_TeV, pt, i_start);
        
        //Eq. 78
        a1=integrale_pp_second_low_en_rate(pf_K_delta,pf_E_min,pf_E_max,Emu_TeV,pt);
        pt->pp_racc_nu_mu = a2 / a1;
 
    }


    Emu_TeV =nu_nu_mu*HPLANCK;
    
    if (Emu_TeV > pt->E_th_pp_delta_approx) {
        //Eq. 71
        i_start = E_min_p_grid_even(pt,pt->griglia_gamma_Np_log,Emu_TeV, 0, pt->gamma_grid_size );
        //i_start=0;
        return integrale_pp_second_high_en_rate(pf_K,Emu_TeV, pt, i_start);
        //printf("i_start=%d gamma_p_min=%e\n", i_start, gamma_p_min);
    } else {
        //Eq. 78
        return integrale_pp_second_low_en_rate(pf_K_delta,pf_E_min,pf_E_max,Emu_TeV,pt);
    }

    return 0.;

}

double E_min_neutrino_mu_1_pp(double E_mu, struct spettro * pt){
    double psida;
    psida=  1.0 - (MEMUC2_TeV*MEMUC2_TeV)/(pt->MPI_kernel_delta_Emin*pt->MPI_kernel_delta_Emin)*0.5;
    return (E_mu/psida)+(pt->MPI_kernel_delta_Emin * pt->MPI_kernel_delta_Emin) / (4 * E_mu)*psida;
}
double E_max_neutrino_mu_1_pp(struct spettro *pt){
    return (pt->gmax*MPC2_TeV - MPC2_TeV);
}


double pp_neutrino_mu_1_kernel_delta(struct spettro *pt,double E_pi ) {
    double q_nu_mu, Ep0_TeV, gamma_p;
    Ep0_TeV = E_pi / (Kpi) + MPC2_TeV;
    gamma_p = Ep0_TeV / MPC2_TeV;
   
    q_nu_mu = pt->pp_racc_nu_mu / (Kpi) * sigma_pp_inel(Ep0_TeV)*
            N_distr_interp(pt->gamma_grid_size, gamma_p, pt->griglia_gamma_Np_log, pt->Np);
    q_nu_mu = 2.0 * q_nu_mu / sqrt(E_pi * E_pi - pt->MPI_kernel_delta * pt->MPI_kernel_delta);

    return q_nu_mu;
}

double pp_neturino_mu_1_kernel(double gamma_p, double E_out_TeV, struct spettro *pt) {
    double Ep_TeV,res,x;
    Ep_TeV = gamma_p*MPC2_TeV;
    x=E_out_TeV / Ep_TeV;
    res= sigma_pp_inel(Ep_TeV) *
            pt->Np[pt->i_griglia_gamma] *
            F_neutrino_mu_1((E_out_TeV / Ep_TeV), Ep_TeV) / pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
     return check_pp_kernel(res,pt,Ep_TeV,x);
}


double F_neutrino_mu_1(double x, double Ep_TeV) {
    //From Eq.66 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //x=E_nu_mu/E_p
    //L=ln(E_p/1TeV)
    double B, beta, k, ly, y, yb, L;
    double a1, a2, a3, a4, a5, res;
    
    y=x/0.427;
    
    if (y>1. || x<0.001) {
        res =0.;
    }
    else{
    
        ly = log(y);
        L = log(Ep_TeV);

        B = 1.75 + 0.204*L + 0.010*L*L;
        beta = 1.0 /  (1.67 + 0.111 * L + 0.0038 * L * L) ;
        k = 1.07 -0.086*L + 0.002*L*L;
        yb = pow(y, beta);

        a1=B*ly/y;

        a2=(1-yb)/(1 + k*yb *(1-yb));
        
        a3=1/ly;

        a4=-(4*beta*yb)/(1-yb);

        a5=-(4*k*beta*yb*(1-2*yb))/(1+k*yb*(1-yb));

        res= a1*(a2*a2*a2*a2)*(a3+a4+a5);
    }
    
  

    return res;
}





//=========================================================================================
// PP -> gamma production
//
//=========================================================================================


double rate_gamma_pp(struct spettro *pt) {
    //From Eq.71 and 78 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //I change the integratio variable from Ep
    //to gamma_p that is more conveninet to the
    //code structure
    //
    double  Emin_pi,Emax_pi, E_gamma_TeV, a1, a2;
    double (*pf_K) (double gamma_p, double nu_pp, struct spettro * pt);
    double (*pf_K_delta) (struct spettro * pt, double x);
    double (*pf_E_min) (double gamma_p, struct spettro * pt);
    double (*pf_E_max) (struct spettro * pt);
    unsigned int i_start;

    pf_K = &pp_gamma_kernel;
    pf_K_delta = &pp_gamma_kernel_delta;
    pf_E_min = &E_min_gamma_pp;
    pf_E_max = &E_max_gamma_pp;
    pt->MPI_kernel_delta=MPI0C2_TeV;
    pt->MPI_kernel_delta_Emin=MPI0C2_TeV;

    //Here we find the n~ (reported in the paper)
    //to find the connection between the standard kernel and the delta-approx
    // At 100 GeV
    //set_pp_racc_gamma=0 means you have to evaluate the connection factor
    if (pt->set_pp_racc_gamma == 0) {
        pt->set_pp_racc_gamma = 1;
        E_gamma_TeV =  pt->E_th_pp_delta_approx;
        
        i_start = E_min_p_grid_even(pt,pt->griglia_gamma_Np_log,E_gamma_TeV, 0, pt->gamma_grid_size );
        //Eq. 71
        pf_K = &pp_gamma_kernel;
        a2 = integrale_pp_second_high_en_rate(pf_K, E_gamma_TeV, pt, i_start);

        //Eq. 78
        a1= integrale_pp_second_low_en_rate(pf_K_delta,pf_E_min,pf_E_max,E_gamma_TeV,pt);

        pt->pp_racc_gamma = a2 / a1;
       
    }

    E_gamma_TeV = pt->nu_1 * HPLANCK_TeV;

    //Eq. 71
    if (E_gamma_TeV > pt->E_th_pp_delta_approx) {
        pf_K = &pp_gamma_kernel;
        i_start = E_min_p_grid_even(pt,pt->griglia_gamma_Np_log,E_gamma_TeV, 0, pt->gamma_grid_size );
        //i_start=0;
        return integrale_pp_second_high_en_rate(pf_K, E_gamma_TeV, pt, i_start);

    } else {
        //Eq. 78
        return integrale_pp_second_low_en_rate(pf_K_delta,pf_E_min,pf_E_max,E_gamma_TeV,pt);

    }

    return 0;
}

double E_min_gamma_pp(double E_gamma, struct spettro *pt){
    return (E_gamma)+(pt->MPI_kernel_delta_Emin * pt->MPI_kernel_delta_Emin) / (4 * E_gamma);
}
double E_max_gamma_pp(struct spettro *pt){
    return (pt->gmax*MPC2_TeV - MPC2_TeV);
}

double pp_gamma_kernel_delta(struct spettro *pt, double E_pi) {
    double qpi,Ep0_TeV, gamma_p;
    Ep0_TeV = MPC2_TeV+ E_pi/Kpi;
    gamma_p=Ep0_TeV/MPC2_TeV;
  
    qpi = pt->pp_racc_gamma / (Kpi)*sigma_pp_inel(Ep0_TeV) * N_distr_interp(pt->gamma_grid_size, gamma_p, pt->griglia_gamma_Np_log, pt->Np);
    return 2.0 * qpi / sqrt(E_pi * E_pi - pt->MPI_kernel_delta * pt->MPI_kernel_delta);
    
}

double pp_gamma_kernel(double gamma_p, double E_out_TeV, struct spettro *pt) {
    double Ep_TeV,x,res;
    Ep_TeV = gamma_p*MPC2_TeV;
    x=E_out_TeV / Ep_TeV;  
    res= sigma_pp_inel(Ep_TeV) *
            pt->Np[pt->i_griglia_gamma] *
            F_gamma((E_out_TeV / Ep_TeV), Ep_TeV) / pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
    return check_pp_kernel(res,pt,Ep_TeV,x);
}


double F_gamma(double x, double Ep_TeV) {
    //From Eq.58 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //x=E_gamma/E_p
    //L=ln(E_p/1TeV)
    double B, beta, k, lx, xb, L, a, a1, a2,res;
    if ( x<0.001) {
        res =0.;
    }
    else{
        lx = log(x);
        L = log(Ep_TeV);

        B = 1.30 + 0.14 * L + 0.011 * L*L;
        beta = 1.0 / (1.79 + 0.11 * L + 0.008 * L * L);
        k = 1.0 / (0.801 + 0.049 * L + 0.014 * L * L);
        xb = pow(x, beta);

        a = (1 - xb) / (1 + k * xb * (1 - xb));
        a1 = (4 * beta * xb) / (1 - xb);
        a2 = (4 * k * beta * xb * (1 - 2 * xb)) / (1 + k * xb * (1 - xb));

        res= B * (lx / x)*(a * a * a * a)*((1 / lx) - a1 - a2);
    }
    return res;
}

//==================================================================
// INTEGRAZIONE PER PP->secondaries  CON METODO TRAPZ E GRIGLIA EQUI-LOG
//
//==================================================================
double integrale_pp_second_low_en_rate(double (*pf_pp_delta_kernel) (struct spettro *pt, double E),
                              double (*E_min_pi) (double gamma_p, struct spettro *pt),
                              double (*E_max_pi) (struct spettro *pt),  
                              double E_out_TeV,
                              struct spettro * pt) {
    //we split in two the integration gdrid
    //to have a better sampling at low E
    //and avoid the spike
    double Emin_pi,E_mid,Emax_pi, a1 ,a2;

    Emin_pi = E_min_pi(E_out_TeV,pt);
    E_mid=Emin_pi*1.1;
    Emax_pi = E_max_pi(pt);
    
    a1=0;
    a2=0;
    
    if (E_mid<=Emax_pi){
        a2= integrale_trap_log_struct(pf_pp_delta_kernel,
            pt,
            E_mid,
            Emax_pi,
            500);
    }
    else{
        E_mid=Emax_pi;
    }
    a1= integrale_trap_log_struct(pf_pp_delta_kernel,
        pt,
        Emin_pi,
        E_mid,
        1000);
    
    //printf("==>a1=%e, a2=%e, E1=%e, E2=%e, E3=%e \n",a1,a2,Emin_pi,E_mid,Emax_pi);
    return a1+a2;
}

                              


double integrale_pp_second_high_en_rate(double (*pf_pp_kernel) (double gamma_p, double E, struct spettro *pt),
                              double E_out_TeV,
                              struct spettro * pt, 
                              unsigned int i_start) {

    double integr, y1, y2, y3, x1, x3;
    double delta;
    //double nu_pp;

    integr = 0;

    //define i_start
    //this the i_th element to which corresponds
    //the gamma_pp tresh hold enegy
    //it is replacing E_gamma in Eq. 71

  
    //nu_pp = pt->nu_1;
    if (i_start<pt->gamma_grid_size -1){
        pt->i_griglia_gamma = i_start;
        x1 = pt->griglia_gamma_Np_log[i_start];
        y1 = pf_pp_kernel(x1, E_out_TeV, pt);
        //CON QUESTO LOOP L'EVENTUALE PUNTO IN ECCESSO
        //CON INDICE DISPARI VIENE SALTATO
        for (pt->i_griglia_gamma = i_start + 1; pt->i_griglia_gamma < pt->gamma_grid_size - 1; pt->i_griglia_gamma++) {
            //printf("i=%d\n",pt->i_griglia_gamma);

            y2 = pf_pp_kernel(pt->griglia_gamma_Np_log[pt->i_griglia_gamma], E_out_TeV, pt);
            pt->i_griglia_gamma++;
            x3 = pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
            y3 = pf_pp_kernel(pt->griglia_gamma_Np_log[pt->i_griglia_gamma], E_out_TeV, pt);

            //printf("i=%d x1=%e x2=%e  x3=%e ",
            //pt->i_griglia_gamma-1,x1,pt->griglia_gamma_Ne_log[pt->i_griglia_gamma-1],x3);

            //QUESTO DELTA RIMANE QUI
            //PERCHE' LA GRIGLIA NON E' EQUISPACED
            //NON PUO ANDARE FUORI DAL LOOP
            delta = (x3 - x1);
            integr += (y1 + 4.0 * y2 + y3) * delta;
            y1 = y3;
            x1 = x3;
            //printf("y1,2,3=%e,%e,%e\n",y1,y2,y3);
            //printf("delta=%e y1=%e y2=%e y3=%e\n",delta,y1,y2,y3);
            //printf("integr=%e\n",integr);
        }
    }
    if (pt->verbose) {
        printf("Integr=%e\n", integr);
    }
    //printf("i_start=%d integr=%e\n",i_start,integr);
    return integr * (0.5 / 3.0);
}
//==================================================================


