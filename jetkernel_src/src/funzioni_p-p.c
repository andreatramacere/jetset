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
    double k, L;
    k = E_th_pp / Ep_TeV;
    L = log(Ep_TeV);
    if (Ep_TeV < E_th_pp) {
        return 0;
    } else {
        return (34.3 + 1.88 * L + 0.25 * L * L)*(1 - k * k * k * k)*(1 - k * k * k * k);
        //return 34;
    }
}
//=========================================================================================




//=========================================================================================
// PP -> e- production
//
//=========================================================================================

double rate_electrons_pp(struct spettro *pt, double Gamma_e) {
    //From Eq.71 and 78 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //I change the integratio variable from Ep
    //to gamma_p that is more conveninet to the
    //code structure
    //
    double gamma_p_min, Ee_TeV, Gamma_1, a1, a2;
    unsigned int i_start;
    double (*pf_K) (double gamma_p, double nu_pp, struct spettro * pt);

    if (pt->set_pp_racc_elec == 0) {
        pt->set_pp_racc_elec = 1;
        Ee_TeV = 0.1;
        Gamma_1 = pt->Gamma;
        pt->Gamma = Ee_TeV / MEC2_TeV;
        gamma_p_min = Ee_TeV / MPC2_TeV;
        i_start = 0;
        while (pt->griglia_gamma_Np_log[i_start] < gamma_p_min) {
            i_start++;
        }
        //!!!!!!!
        //i_start deve essere pari !!!
        //!!!!!!!
        if (i_start % 2 != 0) {
            i_start++;
        }
        pf_K = &pp_electrons_kernel;
        a2 = integrale_pp_electrons_rate(pf_K, pt, i_start);

        a1 = pp_electron_kernel_delta(Ee_TeV, pt);
        pt->pp_racc_elec = a2 / a1;
        //printf("gamma_p_min=%e Ee_TeV=%e i_start=%d\n", gamma_p_min, Ee_TeV,i_start);
        //printf("e- >>>>>>>>>>%e %e %e \n", pt->pp_racc_elec, a1, a2);
        pt->Gamma = Gamma_1;
    }


    Ee_TeV = Gamma_e * MEC2_TeV;
    //Eq. 71
    if (Ee_TeV > 0.1) {
        gamma_p_min = Ee_TeV / MPC2_TeV;
        i_start = 0;
        while (pt->griglia_gamma_Np_log[i_start] < gamma_p_min) {
            i_start++;
        }
        //!!!!!!!
        //i_start deve essere pari !!!
        //!!!!!!!
        if (i_start % 2 != 0) {
            i_start++;
        }
        //printf("i_start=%d gamma_p_min=%e Ee=%e\n", i_start, gamma_p_min, Ee_TeV);
        pf_K = &pp_electrons_kernel;
        //printf("i_start=%d gamma_p_min=%e\n", i_start, gamma_p_min);
        return integrale_pp_electrons_rate(pf_K, pt, i_start);
    } else {

        //Eq. 77
        return pp_electron_kernel_delta(Ee_TeV, pt);
    }



}

double pp_electron_kernel_delta(double Ee_TeV, struct spettro *pt) {
    double qe, Ep_TeV, gamma_p;
    Ep_TeV = Ee_TeV / (Kpp_e) + MPC2_TeV;
    gamma_p = Ep_TeV / MPC2_TeV;
    qe = pt->pp_racc_elec / (Kpp_e) * sigma_pp_inel(Ep_TeV)*
            N_distr_interp(pt->gamma_grid_size, gamma_p, pt->griglia_gamma_Np_log, pt->Np);
    //        N_distr(pt,gamma_p);
    return qe;
}

double pp_electrons_kernel(double gamma_p, double Gamma_e, struct spettro *pt) {
    //printf("gamma_p=%e Ep_TeV=%e Eg_TeV=%e\n", gamma_p, gamma_p * MPC2_TeV, nu_pp * HPLANCK_TeV);
    double Ep_TeV;
    Ep_TeV = gamma_p*MPC2_TeV;
    return sigma_pp_inel(Ep_TeV) *
            pt->Np[pt->i_griglia_gamma] *
            F_electrons((Gamma_e * MEC2_TeV / Ep_TeV), Ep_TeV) / pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
}

double F_electrons(double x, double Ep_TeV) {
    //From Eq.62 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //x=E_gamma/E_p
    //L=ln(E_p/1TeV)
    double B, beta, k, lx, xb, L, a;

    lx = log(x);
    L = log(Ep_TeV);

    B = 1 / (69.5 + 2.65 * L + 0.3 * L * L);
    beta = 1.0 / pow((0.201 + 0.062 * L + 0.00042 * L * L), 0.25);
    k = (0.279 + 0.141 * L + 0.0172 * L * L) / (0.3 + (2.3 + L)*(2.3 + L));
    xb = pow(x, beta);

    a = (1 + k * lx * lx);
    //printf("%e %e %e %e %e %e\n",B,beta,x,xb,L,lx);
    return B * (a * a * a) / (x * (1 + 0.3 / xb)) * (-lx)*(lx)*(lx)*(lx)*(lx);
}

//==================================================================
// INTEGRAZIONE PER PP->electrons  CON METODO DI SIMPSON E GRIGLIA EQUI-LOG
//
//==================================================================

double integrale_pp_electrons_rate(double (*pf_pp_electrons_kernel) (double gamma_p, double Gamma_e, struct spettro *pt), struct spettro * pt, unsigned int i_start) {

    double integr, y1, y2, y3, x1, x3;
    double delta;
    double Gamma_e;

    integr = 0;

    //define i_start
    //this the i_th element to which corresponds
    //the gamma_pp tresh hold enegy
    //it is replacing E_gamma in Eq. 71


    Gamma_e = pt->Gamma;
    //printf("Gamma_e=%e\n",Gamma_e);
    pt->i_griglia_gamma = i_start;
    x1 = pt->griglia_gamma_Np_log[i_start];
    y1 = pf_pp_electrons_kernel(x1, Gamma_e, pt);

    //CON QUESTO LOOP L'EVENTUALE PUNTO IN ECCESSO
    //CON INDICE DISPARI VIENE SALTATO
    for (pt->i_griglia_gamma = i_start + 1; pt->i_griglia_gamma < pt->gamma_grid_size - 1; pt->i_griglia_gamma++) {
        //printf("i=%d\n",pt->i_griglia_gamma);

        y2 = pf_pp_electrons_kernel(pt->griglia_gamma_Np_log[pt->i_griglia_gamma], Gamma_e, pt);
        pt->i_griglia_gamma++;
        x3 = pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
        y3 = pf_pp_electrons_kernel(pt->griglia_gamma_Np_log[pt->i_griglia_gamma], Gamma_e, pt);

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
    if (pt->verbose) {
        printf("Integr=%e\n", integr);
    }
    return integr * (0.5 / 3.0);
}
//=========================================================================================



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
    double gamma_p_min, Emin_pi,Emax_pi, E_gamma, a1, a2, nu1;
    unsigned int i_start;
    double (*pf_K) (double gamma_p, double nu_pp, struct spettro * pt);
    double (*pf_K1) (struct spettro *, double x);

    if (pt->set_pp_racc_gamma == 0) {
        pt->set_pp_racc_gamma = 1;
        E_gamma = 0.1;
        nu1 = pt->nu_1;
        pt->nu_1 = E_gamma / HPLANCK_TeV;
        gamma_p_min = E_gamma / MPC2_TeV;
        i_start = 0;
        while (pt->griglia_gamma_Np_log[i_start] < gamma_p_min) {
            i_start++;
        }
        //!!!!!!!
        //i_start deve essere pari !!!
        //!!!!!!!
        if (i_start % 2 != 0) {
            i_start++;
        }
        pf_K = &pp_gamma_kernel;
        a2 = integrale_pp_gamma_rate(pf_K, pt, i_start);

     
        pf_K1 = &pp_gamma_kernel_delta;
        Emin_pi = (E_gamma)+(MPIC2_TeV * MPIC2_TeV) / (4 * E_gamma);
        Emax_pi = (pt->gmax*MPC2_TeV - MPC2_TeV)*Kpi;
        a1=integrale_trap_log_struct(pf_K1,
            pt,
            Emin_pi,
            Emax_pi,
            1000);
        //a1 = integrale_pp_gamma_rate(pf_K, pt, i_start);


        pt->pp_racc_gamma = a2 / a1;
        //printf("gamma_p_min=%e E_gamma=%e\n", gamma_p_min, E_gamma);
        //printf("gamma >>>>>>>>>>%e %e %e\n", pt->pp_racc_gamma, a1, a2);
        pt->nu_1 = nu1;
    }


    E_gamma = pt->nu_1 * HPLANCK_TeV;


    //Eq. 71
    if (E_gamma > 0.1) {
        gamma_p_min = E_gamma / MPC2_TeV;
        pf_K = &pp_gamma_kernel;
        i_start = 0;
        while (pt->griglia_gamma_Np_log[i_start] < gamma_p_min) {
            i_start++;
        }
        //!!!!!!!
        //i_start deve essere pari !!!
        //!!!!!!!
        if (i_start % 2 != 0) {
            i_start++;
        }
        //printf("i_start=%d gamma_p_min=%e E_gamma=%e\n", i_start, gamma_p_min, E_gamma);
        return integrale_pp_gamma_rate(pf_K, pt, i_start);

    } else {
        //Eq. 78
        //Emin pi

        Emin_pi = (E_gamma)+(MPIC2_TeV * MPIC2_TeV) / (4 * E_gamma);
        Emax_pi = (pt->gmax*MPC2_TeV - MPC2_TeV)*Kpi;
        //printf("gamma_p_min=%e E_gamma=%e\n",gamma_p_min,E_gamma);
        pf_K1 = &pp_gamma_kernel_delta;
        return integrale_trap_log_struct(pf_K1,
            pt,
            Emin_pi,
            Emax_pi,
            10000);
    }
}


double pp_gamma_kernel_delta(struct spettro *pt, double Epi_TeV) {
    double qpi,Ep_TeV, Gamma;
    Ep_TeV = Epi_TeV/Kpi +MPC2_TeV;
    Gamma=Ep_TeV/MPC2_TeV;
    qpi = pt->pp_racc_gamma / (Kpi)*sigma_pp_inel(Ep_TeV) * N_distr_interp(pt->gamma_grid_size, Gamma, pt->griglia_gamma_Np_log, pt->Np);
    return 2.0 * qpi / sqrt(Epi_TeV * Epi_TeV - MPIC2_TeV * MPIC2_TeV);
}

double pp_gamma_kernel(double gamma_p, double nu_pp, struct spettro *pt) {
    //printf("gamma_p=%e Ep_TeV=%e Eg_TeV=%e\n", gamma_p, gamma_p * MPC2_TeV, nu_pp * HPLANCK_TeV);
    double Ep_TeV;
    Ep_TeV = gamma_p*MPC2_TeV;
    return sigma_pp_inel(Ep_TeV) *
            pt->Np[pt->i_griglia_gamma] *
            F_gamma((nu_pp * HPLANCK_TeV / Ep_TeV), Ep_TeV) / pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
}

double F_gamma(double x, double Ep_TeV) {
    //From Eq.58 in Kelner et al. 2006
    //astro-ph.06066058v1
    //PHYSICAL REVIEW D 74, 034018 (2006)
    //x=E_gamma/E_p
    //L=ln(E_p/1TeV)
    double B, beta, k, lx, xb, L, a, a1, a2;

    lx = log(x);
    L = log(Ep_TeV);

    B = 1.30 + 0.14 * L + 0.011 * L*L;
    beta = 1.0 / (1.79 + 0.11 * L + 0.008 * L * L);
    k = 1.0 / (0.801 + 0.049 * L + 0.014 * L * L);
    xb = pow(x, beta);

    a = (1 - xb) / (1 + k * xb * (1 - xb));
    a1 = (4 * beta * xb) / (1 - xb);
    a2 = (4 * k * beta * xb * (1 - 2 * xb)) / (1 + k * xb * (1 - xb));

    return B * (lx / x)*(a * a * a * a)*((1 / lx) - a1 - a2);
    //return 3.06*exp(-9.47*pow(x,0.75));
}

//==================================================================
// INTEGRAZIONE PER PP->gamma  CON METODO DI SIMPSON E GRIGLIA EQUI-LOG
//
//==================================================================

double integrale_pp_gamma_rate(double (*pf_pp_gamma_kernel) (double gamma_p, double nu_pp, struct spettro *pt), struct spettro * pt, unsigned int i_start) {

    double integr, y1, y2, y3, x1, x3;
    double delta;
    double nu_pp;

    integr = 0;

    //define i_start
    //this the i_th element to which corresponds
    //the gamma_pp tresh hold enegy
    //it is replacing E_gamma in Eq. 71

  
    nu_pp = pt->nu_1;
    pt->i_griglia_gamma = i_start;
    x1 = pt->griglia_gamma_Np_log[i_start];
    y1 = pf_pp_gamma_kernel(x1, nu_pp, pt);
    //CON QUESTO LOOP L'EVENTUALE PUNTO IN ECCESSO
    //CON INDICE DISPARI VIENE SALTATO
    for (pt->i_griglia_gamma = i_start + 1; pt->i_griglia_gamma < pt->gamma_grid_size - 1; pt->i_griglia_gamma++) {
        //printf("i=%d\n",pt->i_griglia_gamma);

        y2 = pf_pp_gamma_kernel(pt->griglia_gamma_Np_log[pt->i_griglia_gamma], nu_pp, pt);
        pt->i_griglia_gamma++;
        x3 = pt->griglia_gamma_Np_log[pt->i_griglia_gamma];
        y3 = pf_pp_gamma_kernel(pt->griglia_gamma_Np_log[pt->i_griglia_gamma], nu_pp, pt);

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
    if (pt->verbose) {
        printf("Integr=%e\n", integr);
    }
    //printf("i_start=%d integr=%e\n",i_start,integr);
    return integr * (0.5 / 3.0);
}
//==================================================================


