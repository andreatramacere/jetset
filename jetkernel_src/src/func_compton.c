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

static unsigned int lower_bound_gamma_grid(const double *gamma_grid, unsigned int gamma_grid_size, double gmin) {
    unsigned int lo, hi, mid;

    lo = 0;
    hi = gamma_grid_size;
    while (lo < hi) {
        mid = lo + (hi - lo) / 2;
        if (gamma_grid[mid] < gmin) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo;
}

static size_t external_angle_flat_index(unsigned int nu_id, unsigned int angle_id, unsigned int angle_n_int)
{
    return ((size_t)nu_id) * ((size_t)angle_n_int) + ((size_t)angle_id);
}

static double clamp_to_interval(double x, double x_min, double x_max)
{
    if (x < x_min) {
        return x_min;
    }
    if (x > x_max) {
        return x_max;
    }
    return x;
}

static int use_ec_angle_dep_full(const struct blob *pt,
                                 const struct spectrum_external *spec,
                                 int use_drf,
                                 unsigned int nu_seed_size)
{
    const double *n_theta;

    if (pt == NULL || spec == NULL) {
        return 0;
    }
    if (pt->core.EC_kernel != EC_KERNEL_ANGLE_DEP_FULL) {
        return 0;
    }
    if (spec->mu == NULL || spec->angle_n_int < 2U || spec->angle_nu_size < nu_seed_size) {
        return 0;
    }

    n_theta = (use_drf != 0) ? spec->n_nu_theta_DRF : spec->n_nu_theta;
    if (n_theta == NULL) {
        return 0;
    }

    return 1;
}

static double eval_ec_dsigma_depsilon_s_angle_dep(double gamma,
                                                   double epsilon,
                                                   double epsilon_s,
                                                   double one_minus_cos_psi)
{
    double y, bar_epsilon, ratio, Xi, prefactor;
    double epsilon_min, epsilon_s_max;

    if (gamma <= 1.0 || epsilon <= 0.0 || epsilon_s <= 0.0 || epsilon_s >= gamma) {
        return 0.0;
    }
    if (one_minus_cos_psi <= 0.0) {
        return 0.0;
    }

    y = 1.0 - (epsilon_s / gamma);
    if (y <= 0.0) {
        return 0.0;
    }

    bar_epsilon = gamma * epsilon * one_minus_cos_psi;
    if (bar_epsilon <= 0.0) {
        return 0.0;
    }

    epsilon_min = epsilon_s / (2.0 * one_minus_cos_psi * gamma * (gamma - epsilon_s));
    if (epsilon < epsilon_min) {
        return 0.0;
    }

    epsilon_s_max = (2.0 * one_minus_cos_psi * epsilon * gamma * gamma) /
                    (1.0 + 2.0 * one_minus_cos_psi * epsilon * gamma);
    if (epsilon_s > epsilon_s_max) {
        return 0.0;
    }

    ratio = epsilon_s / (gamma * bar_epsilon * y);
    Xi = y + (1.0 / y) - (2.0 * ratio) + (ratio * ratio);
    if (Xi <= 0.0 || !isfinite(Xi)) {
        return 0.0;
    }

    prefactor = (3.0 * SIGTH) / (8.0 * gamma * bar_epsilon);
    return prefactor * Xi;
}

static double eval_ec_phi_integral_kernel(double gamma,
                                          double epsilon,
                                          double epsilon_s,
                                          const double *a_phi,
                                          unsigned int n_phi,
                                          double dphi)
{
    unsigned int i_phi;
    double one_minus_cos_psi, dsigma, phi_integral;

    phi_integral = 0.0;
    for (i_phi = 0; i_phi < n_phi; i_phi++) {
        one_minus_cos_psi = a_phi[i_phi];
        if (one_minus_cos_psi <= 0.0) {
            continue;
        }
        dsigma = eval_ec_dsigma_depsilon_s_angle_dep(gamma, epsilon, epsilon_s, one_minus_cos_psi);
        if (dsigma > 0.0) {
            phi_integral += one_minus_cos_psi * dsigma;
        }
    }

    return phi_integral * dphi;
}

static double integrale_IC_angle_dep_full(struct blob *pt,
                                          const struct spectrum_external *spec,
                                          int use_drf,
                                          const double *nu_seed,
                                          unsigned int nu_seed_size,
                                          double a,
                                          double b,
                                          int stat_frame,
                                          double nu_IC_out)
{
    unsigned int gamma_grid_size;
    unsigned int angle_n_int;
    unsigned int n_phi;
    unsigned int ID, ID_mu, ID_gamma, gamma_start;
    double dphi;
    double mu_s, mu_s_prime, sin_mu_s_prime;
    double epsilon_s;
    double epsilon_in;
    double gamma_min_kin;
    double mu, sin_mu, cos_psi, one_minus_cos_psi;
    double sin_arg;
    double mu_integral, dmu;
    double kernel_phi;
    double nu_integral;
    const double *n_theta;
    double *Integrand_over_gamma_grid, *Ne_IC, *griglia_gamma_Ne_log_IC, *integr_nu, *integr_mu, *cos_phi, *A_grid;
    size_t mu_idx;

    if (pt == NULL || spec == NULL || nu_seed == NULL) {
        return 0.0;
    }
    if (nu_seed_size == 0U || nu_seed_size > spec->angle_nu_size || spec->angle_n_int < 2U || spec->mu == NULL) {
        return 0.0;
    }

    n_theta = (use_drf != 0) ? spec->n_nu_theta_DRF : spec->n_nu_theta;
    if (n_theta == NULL) {
        return 0.0;
    }

    gamma_grid_size = pt->emitters.gamma_grid_size;
    angle_n_int = spec->angle_n_int;
    n_phi = pt->core.EC_angle_n_phi;
    if (n_phi < 4U) {
        n_phi = 4U;
    }
    dphi = (2.0 * pi) / (double)n_phi;

    Integrand_over_gamma_grid = (double *)calloc(gamma_grid_size, sizeof(double));
    griglia_gamma_Ne_log_IC = (double *)calloc(gamma_grid_size, sizeof(double));
    Ne_IC = (double *)calloc(gamma_grid_size, sizeof(double));
    integr_nu = (double *)calloc(nu_seed_size, sizeof(double));
    integr_mu = (double *)calloc(angle_n_int, sizeof(double));
    cos_phi = (double *)calloc(n_phi, sizeof(double));
    A_grid = (double *)calloc((size_t)angle_n_int * (size_t)n_phi, sizeof(double));

    if (Integrand_over_gamma_grid == NULL ||
        griglia_gamma_Ne_log_IC == NULL ||
        Ne_IC == NULL ||
        integr_nu == NULL ||
        integr_mu == NULL ||
        cos_phi == NULL ||
        A_grid == NULL) {
        free(Integrand_over_gamma_grid);
        free(griglia_gamma_Ne_log_IC);
        free(Ne_IC);
        free(integr_nu);
        free(integr_mu);
        free(cos_phi);
        free(A_grid);
        return 0.0;
    }

    for (ID = 0U; ID < n_phi; ID++) {
        cos_phi[ID] = cos(((double)ID + 0.5) * dphi);
    }

    mu_s = cos(pt->core.theta * Deg_to_Rad);
    mu_s_prime = (mu_s - pt->core.beta_Gamma) / (1.0 - pt->core.beta_Gamma * mu_s);
    mu_s_prime = clamp_to_interval(mu_s_prime, -1.0, 1.0);
    sin_arg = 1.0 - mu_s_prime * mu_s_prime;
    if (sin_arg < 0.0) {
        sin_arg = 0.0;
    }
    sin_mu_s_prime = sqrt(sin_arg);
    epsilon_s = HPLANCK * nu_IC_out * one_by_MEC2;
    if (epsilon_s <= 0.0) {
        free(Integrand_over_gamma_grid);
        free(griglia_gamma_Ne_log_IC);
        free(Ne_IC);
        free(integr_nu);
        free(integr_mu);
        free(cos_phi);
        free(A_grid);
        return 0.0;
    }

    for (ID_mu = 0U; ID_mu < angle_n_int; ID_mu++) {
        mu = clamp_to_interval(spec->mu[ID_mu], -1.0, 1.0);
        sin_arg = 1.0 - mu * mu;
        if (sin_arg < 0.0) {
            sin_arg = 0.0;
        }
        sin_mu = sqrt(sin_arg);
        for (ID = 0U; ID < n_phi; ID++) {
            cos_psi = mu * mu_s_prime + sin_mu * sin_mu_s_prime * cos_phi[ID];
            cos_psi = clamp_to_interval(cos_psi, -1.0, 1.0);
            one_minus_cos_psi = 1.0 - cos_psi;
            if (one_minus_cos_psi < 0.0) {
                one_minus_cos_psi = 0.0;
            }
            A_grid[external_angle_flat_index(ID_mu, ID, n_phi)] = one_minus_cos_psi;
        }
    }

    set_N_distr_for_Compton(pt, b, nu_IC_out, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);

    for (ID = 0U; ID < nu_seed_size; ID++) {
        if (nu_seed[ID] < a || nu_seed[ID] > b) {
            integr_nu[ID] = 0.0;
            continue;
        }

        epsilon_in = HPLANCK * nu_seed[ID] * one_by_MEC2;
        if (epsilon_in <= 0.0) {
            integr_nu[ID] = 0.0;
            continue;
        }

        gamma_min_kin = 0.5 * epsilon_s * (1.0 + sqrt(1.0 + (1.0 / (epsilon_s * epsilon_in))));
        gamma_start = lower_bound_gamma_grid(griglia_gamma_Ne_log_IC, gamma_grid_size, gamma_min_kin);
        for (ID_gamma = 0U; ID_gamma < gamma_start; ID_gamma++) {
            Integrand_over_gamma_grid[ID_gamma] = 0.0;
        }

        for (ID_mu = 0U; ID_mu < angle_n_int; ID_mu++) {
            mu_idx = external_angle_flat_index(ID, ID_mu, angle_n_int);
            if (n_theta[mu_idx] <= 0.0) {
                integr_mu[ID_mu] = 0.0;
                continue;
            }

            for (ID_gamma = gamma_start; ID_gamma < gamma_grid_size; ID_gamma++) {
                kernel_phi = eval_ec_phi_integral_kernel(griglia_gamma_Ne_log_IC[ID_gamma],
                                                         epsilon_in,
                                                         epsilon_s,
                                                         &(A_grid[external_angle_flat_index(ID_mu, 0U, n_phi)]),
                                                         n_phi,
                                                         dphi);
                Integrand_over_gamma_grid[ID_gamma] = Ne_IC[ID_gamma] * kernel_phi;
            }

            integr_mu[ID_mu] = n_theta[mu_idx] *
                               integr_simp_grid_equilog(griglia_gamma_Ne_log_IC,
                                                        Integrand_over_gamma_grid,
                                                        gamma_grid_size);
        }

        mu_integral = 0.0;
        for (ID_mu = 0U; ID_mu + 1U < angle_n_int; ID_mu++) {
            dmu = fabs(spec->mu[ID_mu + 1U] - spec->mu[ID_mu]);
            mu_integral += 0.5 * dmu * (integr_mu[ID_mu] + integr_mu[ID_mu + 1U]);
        }

        integr_nu[ID] = vluce_cm * mu_integral;
    }

    nu_integral = trapzd_array_arbritary_grid((double *)nu_seed, integr_nu, nu_seed_size);

    free(Integrand_over_gamma_grid);
    free(griglia_gamma_Ne_log_IC);
    free(Ne_IC);
    free(integr_nu);
    free(integr_mu);
    free(cos_phi);
    free(A_grid);

    return nu_integral;
}

static double eval_external_ec_rate(struct blob *pt,
                                    const struct spectrum_external *spec,
                                    unsigned int nu_seed_size,
                                    double nu_IC_out,
                                    double nu_IC_out_stat,
                                    const char *field_name)
{
    const double *nu_seed;
    const double *n_seed;
    double nu_min, nu_max, nu_out_eval;
    int use_drf;
    int do_angle_dep_full;

    if (pt == NULL || spec == NULL) {
        return 0.0;
    }

    use_drf = (pt->core.EC_stat != 0) ? 1 : 0;
    if (use_drf != 0) {
        nu_seed = spec->nu_DRF;
        n_seed = spec->n_nu_DRF;
        nu_min = spec->nu_min_DRF;
        nu_max = spec->nu_max_DRF;
        nu_out_eval = nu_IC_out_stat;
    } else {
        nu_seed = spec->nu;
        n_seed = spec->n_nu;
        nu_min = spec->nu_min;
        nu_max = spec->nu_max;
        nu_out_eval = nu_IC_out;
    }

    do_angle_dep_full = use_ec_angle_dep_full(pt, spec, use_drf, nu_seed_size);
    if (pt->core.verbose > 1) {
        printf("%s external field transformation=blob, EC angular kernel=%s, angular photon density=%s, axisymmetric field=true, density=n_nu(theta) per Hz per sr, n_mu=%u, n_phi=%u\n",
               field_name,
               do_angle_dep_full ? "angle_dep_full" : "isotropic",
               do_angle_dep_full ? "enabled" : "disabled",
               do_angle_dep_full ? spec->angle_n_int : 0U,
               pt->core.EC_angle_n_phi);
    }

    if (do_angle_dep_full) {
        return integrale_IC_angle_dep_full(pt,
                                           spec,
                                           use_drf,
                                           nu_seed,
                                           nu_seed_size,
                                           nu_min,
                                           nu_max,
                                           pt->core.EC_stat,
                                           nu_out_eval);
    }

    return integrale_IC(pt,
                        nu_seed,
                        n_seed,
                        nu_seed_size,
                        nu_min,
                        nu_max,
                        pt->core.EC_stat,
                        nu_out_eval);
}
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
    double rate_comp;
    double nu_IC_out_stat;
    const double *nu_seed;
    const double *n_seed;
    unsigned int nu_seed_size;

    rate_comp = 0.0;
    nu_seed_size = pt_GR->core.nu_seed_size;
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
    if (pt_GR->core.SSC == 0 && pt_GR->core.ord_comp == 1) {
        switch (pt_GR->core.EC) {
            case 1:
                if (nu_IC_out < pt_GR->Disk.ec.spec.nu_max) {
                    if (pt_GR->core.verbose > 1) {
                        printf("Disk\n");
                        printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->Disk.spec.nu_min);
                        printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->Disk.spec.nu_max);
                    }
                    rate_comp = eval_external_ec_rate(pt_GR, &(pt_GR->Disk.spec), nu_seed_size, nu_IC_out, nu_IC_out_stat, "Disk");
                }
                break;

            case 2:
                if (nu_IC_out < pt_GR->BLR.ec.spec.nu_max) {
                    if (pt_GR->core.verbose > 1) {
                        printf("BLR\n");
                        printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->BLR.spec.nu_min);
                        printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->BLR.spec.nu_max);
                    }
                    rate_comp = eval_external_ec_rate(pt_GR, &(pt_GR->BLR.spec), nu_seed_size, nu_IC_out, nu_IC_out_stat, "BLR");
                }
                break;

            case 3:
                if (nu_IC_out < pt_GR->DT.ec.spec.nu_max) {
                    if (pt_GR->core.verbose > 1) {
                        printf("DT\n");
                        printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->DT.spec.nu_min);
                        printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->DT.spec.nu_max);
                    }
                    rate_comp = eval_external_ec_rate(pt_GR, &(pt_GR->DT.spec), nu_seed_size, nu_IC_out, nu_IC_out_stat, "DT");
                }
                break;

            case 4:
                if (nu_IC_out < pt_GR->Star.ec.spec.nu_max) {
                    if (pt_GR->core.verbose > 1) {
                        printf("Star\n");
                        printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->Star.spec.nu_min);
                        printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->Star.spec.nu_max);
                    }
                    rate_comp = eval_external_ec_rate(pt_GR, &(pt_GR->Star.spec), nu_seed_size, nu_IC_out, nu_IC_out_stat, "Star");
                }
                break;

            case 5:
                if (nu_IC_out < pt_GR->CMB.ec.spec.nu_max) {
                    if (pt_GR->core.verbose > 1) {
                        printf("CMB\n");
                        printf("nu_start_CMB_seed=%e\n", pt_GR->CMB.spec.nu_min);
                        printf("nu_stop_CMB_seed=%e\n", pt_GR->CMB.spec.nu_max);
                    }
                    rate_comp = eval_external_ec_rate(pt_GR, &(pt_GR->CMB.spec), nu_seed_size, nu_IC_out, nu_IC_out_stat, "CMB");
                }
                break;

            case 6:
                if (nu_IC_out < pt_GR->Corona.ec.spec.nu_max) {
                    if (pt_GR->core.verbose > 1) {
                        printf("Corona\n");
                        printf("(blob rest frame) nu_start_EC_seed=%e\n", pt_GR->Corona.spec.nu_min);
                        printf("(blob rest frame) nu_stop_EC_seed=%e\n", pt_GR->Corona.spec.nu_max);
                    }
                    rate_comp = eval_external_ec_rate(pt_GR, &(pt_GR->Corona.spec), nu_seed_size, nu_IC_out, nu_IC_out_stat, "Corona");
                }
                break;

            default:
                break;
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
    unsigned int ID, ID_gamma, gamma_start, gamma_grid_size;
    double *Integrand_over_gamma_grid, *Ne_IC, *griglia_gamma_Ne_log_IC, *integr_gamma;
    double ic_kernel;

    gamma_grid_size = pt->emitters.gamma_grid_size;
    integr_nu = 0.0;
    Integrand_over_gamma_grid = (double *) calloc(gamma_grid_size, sizeof(double));
    griglia_gamma_Ne_log_IC = (double *) calloc(gamma_grid_size, sizeof(double));
    integr_gamma = (double *) calloc(nu_seed_size, sizeof(double));
    Ne_IC = (double *) calloc(gamma_grid_size, sizeof(double));

    if (Integrand_over_gamma_grid == NULL || griglia_gamma_Ne_log_IC == NULL || integr_gamma == NULL || Ne_IC == NULL) {
        free(Integrand_over_gamma_grid);
        free(griglia_gamma_Ne_log_IC);
        free(integr_gamma);
        free(Ne_IC);
        return 0.0;
    }

    set_N_distr_for_Compton(pt, b, nu_IC_out, stat_frame, Ne_IC, griglia_gamma_Ne_log_IC);

    for (ID = 0; ID < nu_seed_size; ID++) {
        if (nu_seed[ID] <= b && nu_seed[ID] >= a) {
            nu_IC_in = nu_seed[ID];

            if (pt->core.bulk_compton == 0) {
                if (pt->core.do_IC_down_scattering == 0 && nu_IC_in > nu_IC_out) {
                    integr_gamma[ID] = 0.0;
                    continue;
                }

                if (nu_IC_in > 0.0 && nu_IC_out > 0.0) {
                    double g_min_kin;
                    if (nu_IC_in > nu_IC_out) {
                        g_min_kin = sqrt(nu_IC_in / (4.0 * nu_IC_out));
                    } else {
                        double epsilon_0, epsilon_1;
                        epsilon_0 = HPLANCK * nu_IC_in * one_by_MEC2;
                        epsilon_1 = HPLANCK * nu_IC_out * one_by_MEC2;
                        g_min_kin = 0.5 * epsilon_1 * (1.0 + sqrt(1.0 + (1.0 / (epsilon_1 * epsilon_0))));
                    }
                    gamma_start = lower_bound_gamma_grid(griglia_gamma_Ne_log_IC, gamma_grid_size, g_min_kin);
                } else {
                    gamma_start = 0;
                }

                for (ID_gamma = 0; ID_gamma < gamma_start; ID_gamma++) {
                    Integrand_over_gamma_grid[ID_gamma] = 0.0;
                }
                for (ID_gamma = gamma_start; ID_gamma < gamma_grid_size; ID_gamma++) {
                    ic_kernel = f_compton_K1(pt, griglia_gamma_Ne_log_IC[ID_gamma], nu_IC_out, nu_IC_in);
                    Integrand_over_gamma_grid[ID_gamma] = ic_kernel * Ne_IC[ID_gamma];
                }
            } else {
                double nu_left, nu_right;
                if (ID < nu_seed_size - 1) {
                    nu_left = nu_seed[ID];
                    nu_right = nu_seed[ID + 1];
                } else {
                    nu_left = (ID == 0) ? nu_seed[ID] : nu_seed[ID - 1];
                    nu_right = nu_seed[ID];
                }

                ic_kernel = f_compton_bulk(pt, 1.0, nu_IC_out, nu_left, nu_right);
                if (ic_kernel == 0.0) {
                    integr_gamma[ID] = 0.0;
                    continue;
                }

                for (ID_gamma = 0; ID_gamma < gamma_grid_size; ID_gamma++) {
                    Integrand_over_gamma_grid[ID_gamma] = ic_kernel * Ne_IC[ID_gamma];
                }
            }

            integr_gamma[ID] = n_seed[ID] * integr_simp_grid_equilog(griglia_gamma_Ne_log_IC, Integrand_over_gamma_grid, gamma_grid_size);
        } else {
            integr_gamma[ID] = 0.0;
        }
    }
    integr_nu = trapzd_array_arbritary_grid((double *) nu_seed, integr_gamma, nu_seed_size);

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
    const double *nu_seed;
    const double *n_seed;

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

	    //EC Corona
	    if (pt_spec->core.do_EC_Corona == 1 ) {

	    	if (pt_spec->core.verbose>1) {
	    		printf("Corona\n");
	    		printf("nu_start_EC_seed=%e\n", pt_spec->Corona.spec.nu_min);
		    		printf("nu_stop_EC_seed=%e\n", pt_spec->Corona.spec.nu_max);
		    	}
		    	nu_seed = pt_spec->Corona.spec.nu;
		    	n_seed = pt_spec->Corona.spec.n_nu;
		    	comp_cooling += integrale_IC_cooling(pt_spec,
		                nu_seed,
		                n_seed,
	                nu_seed_size,
	    			pt_spec->Corona.spec.nu_min,
	    			pt_spec->Corona.spec.nu_max,
	    			gamma);
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
