#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Blazar_SED.h"

/*
 * Internal gamma-gamma absorption workflow
 * ---------------------------------------
 * 1) `reset_internal_abs_store`/`free_internal_abs_store` initialize and release
 *    per-component caches for BLR, DT, and Corona.
 * 2) `eval_internal_abs_tau` computes tau(nu) for one component:
 *    - resolve the target component and geometry,
 *    - sample the soft-photon seed field (`sample_seed_field`),
 *    - build integration grids in path length (R_H), angle (mu), and seed frequency,
 *    - integrate sigma_gg * n_soft over (nu_soft, mu, R_H) with trapezoids.
 * 3) Results are stored in `pt->core.internal_abs.<component>` (`nu_tau`, `tau`,
 *    flags, and numeric setup) and can be recomputed by
 *    `recompute_internal_absorption_tau`.
 * 4) `get_internal_abs_tau_at_nu` interpolates each enabled component in log-log
 *    space and returns the summed opacity at a requested observed frequency.
 *
 * Implementation note: during evaluation the solver temporarily changes
 * `pt->core.R_H` to sample geometry-dependent fields and always restores it.
 */

#define INTABS_MIN_Y 1.0e-200
#define INTABS_MIN_ONE_MINUS_MU 1.0e-20

typedef enum {
    INTABS_COMP_INVALID = 0,
    INTABS_COMP_BLR = 1,
    INTABS_COMP_DT = 2,
    INTABS_COMP_CORONA = 3
} intabs_comp_t;

/* Reset one component bookkeeping and cached tau pointers to a known empty state. */
static void init_internal_abs_component(struct internal_abs_component *comp) {
    if (comp == NULL) {
        return;
    }

    comp->is_enabled = 0;
    comp->is_valid = 0;
    comp->use_R_H_profile_extrapolation = 0;
    comp->use_sigma_gamma_gamma_fast = 0;
    comp->peak_mode = 0;
    comp->N_soft = 0;
    comp->N_hard = 0;
    comp->N_R_H = 0;
    comp->N_theta = 0;
    comp->tau_size = 0;
    comp->nu_min = 0.0;
    comp->nu_src_max = 0.0;
    comp->nu_tau = NULL;
    comp->tau = NULL;
}

/* Free one component tau cache (if allocated) and mark the cache as invalid. */
static void free_internal_abs_component(struct internal_abs_component *comp) {
    if (comp == NULL) {
        return;
    }

    if (comp->nu_tau != NULL) {
        free(comp->nu_tau);
        comp->nu_tau = NULL;
    }
    if (comp->tau != NULL) {
        free(comp->tau);
        comp->tau = NULL;
    }
    comp->tau_size = 0;
    comp->is_valid = 0;
}

/*
 * Reset all internal-absorption component slots on a blob.
 * This does not free memory; it zeroes the component metadata/handles.
 */
void reset_internal_abs_store(struct blob *pt) {
    if (pt == NULL) {
        return;
    }
    init_internal_abs_component(&(pt->core.internal_abs.BLR));
    init_internal_abs_component(&(pt->core.internal_abs.DT));
    init_internal_abs_component(&(pt->core.internal_abs.Corona));
}

/* Release all per-component tau caches stored in the blob. */
void free_internal_abs_store(struct blob *pt) {
    if (pt == NULL) {
        return;
    }
    free_internal_abs_component(&(pt->core.internal_abs.BLR));
    free_internal_abs_component(&(pt->core.internal_abs.DT));
    free_internal_abs_component(&(pt->core.internal_abs.Corona));
}

/* Map the user-facing seed field name to the internal component identifier. */
static intabs_comp_t parse_internal_abs_component(const char *seed_photons_name) {
    if (seed_photons_name == NULL) {
        return INTABS_COMP_INVALID;
    }
    if (strcmp(seed_photons_name, "BLR") == 0) {
        return INTABS_COMP_BLR;
    }
    if (strcmp(seed_photons_name, "DT") == 0) {
        return INTABS_COMP_DT;
    }
    if (strcmp(seed_photons_name, "Corona") == 0) {
        return INTABS_COMP_CORONA;
    }
    return INTABS_COMP_INVALID;
}

/* Return the selected internal-absorption component struct inside `pt`. */
static struct internal_abs_component *get_internal_abs_component_ptr(struct blob *pt, intabs_comp_t comp_id) {
    if (pt == NULL) {
        return NULL;
    }

    if (comp_id == INTABS_COMP_BLR) {
        return &(pt->core.internal_abs.BLR);
    }
    if (comp_id == INTABS_COMP_DT) {
        return &(pt->core.internal_abs.DT);
    }
    if (comp_id == INTABS_COMP_CORONA) {
        return &(pt->core.internal_abs.Corona);
    }

    return NULL;
}

/*
 * Ensure `comp->nu_tau` and `comp->tau` exist with the requested size.
 * Reallocates if size changed; zeroes existing arrays if size is unchanged.
 */
static int ensure_tau_arrays(struct internal_abs_component *comp, unsigned int tau_size) {
    if (comp == NULL) {
        return -1;
    }

    if (tau_size == 0) {
        free_internal_abs_component(comp);
        return 0;
    }

    if ((comp->nu_tau == NULL) || (comp->tau == NULL) || (comp->tau_size != tau_size)) {
        if (comp->nu_tau != NULL) {
            free(comp->nu_tau);
            comp->nu_tau = NULL;
        }
        if (comp->tau != NULL) {
            free(comp->tau);
            comp->tau = NULL;
        }

        comp->nu_tau = (double *)calloc((size_t)tau_size, sizeof(double));
        comp->tau = (double *)calloc((size_t)tau_size, sizeof(double));
        if ((comp->nu_tau == NULL) || (comp->tau == NULL)) {
            if (comp->nu_tau != NULL) {
                free(comp->nu_tau);
                comp->nu_tau = NULL;
            }
            if (comp->tau != NULL) {
                free(comp->tau);
                comp->tau = NULL;
            }
            comp->tau_size = 0;
            return -1;
        }
    } else {
        memset(comp->nu_tau, 0, ((size_t)tau_size) * sizeof(double));
        memset(comp->tau, 0, ((size_t)tau_size) * sizeof(double));
    }

    comp->tau_size = tau_size;
    return 0;
}

struct internal_abs_eval_workspace {
    double *nu_soft;
    double *n_soft;
    double *eps_soft;
    double *dnu_soft;
    double *mu_grid;
    double *one_minus_mu;
    double *dmu_grid;
    double *R_H_grid;
    double *d_rh_grid;
    double *nu_soft_ref;
    double *n_soft_ref;
    double *nu_soft_tmp;
    double *n_soft_tmp;
};

/* Initialize all workspace pointers to NULL so cleanup is always safe. */
static void init_internal_abs_eval_workspace(struct internal_abs_eval_workspace *ws) {
    if (ws == NULL) {
        return;
    }
    memset(ws, 0, sizeof(*ws));
}

/* Free every temporary array used by the internal-absorption integration. */
static void free_internal_abs_eval_workspace(struct internal_abs_eval_workspace *ws) {
    if (ws == NULL) {
        return;
    }
    if (ws->nu_soft != NULL) {
        free(ws->nu_soft);
        ws->nu_soft = NULL;
    }
    if (ws->n_soft != NULL) {
        free(ws->n_soft);
        ws->n_soft = NULL;
    }
    if (ws->eps_soft != NULL) {
        free(ws->eps_soft);
        ws->eps_soft = NULL;
    }
    if (ws->dnu_soft != NULL) {
        free(ws->dnu_soft);
        ws->dnu_soft = NULL;
    }
    if (ws->mu_grid != NULL) {
        free(ws->mu_grid);
        ws->mu_grid = NULL;
    }
    if (ws->one_minus_mu != NULL) {
        free(ws->one_minus_mu);
        ws->one_minus_mu = NULL;
    }
    if (ws->dmu_grid != NULL) {
        free(ws->dmu_grid);
        ws->dmu_grid = NULL;
    }
    if (ws->R_H_grid != NULL) {
        free(ws->R_H_grid);
        ws->R_H_grid = NULL;
    }
    if (ws->d_rh_grid != NULL) {
        free(ws->d_rh_grid);
        ws->d_rh_grid = NULL;
    }
    if (ws->nu_soft_ref != NULL) {
        free(ws->nu_soft_ref);
        ws->nu_soft_ref = NULL;
    }
    if (ws->n_soft_ref != NULL) {
        free(ws->n_soft_ref);
        ws->n_soft_ref = NULL;
    }
    if (ws->nu_soft_tmp != NULL) {
        free(ws->nu_soft_tmp);
        ws->nu_soft_tmp = NULL;
    }
    if (ws->n_soft_tmp != NULL) {
        free(ws->n_soft_tmp);
        ws->n_soft_tmp = NULL;
    }
}

/*
 * Shared solver exit path: restore `pt->core.R_H`, free temporaries,
 * and invalidate component cache on errors.
 */
static int finalize_internal_abs_eval(struct blob *pt,
                                      struct internal_abs_component *comp,
                                      double R_H_saved,
                                      int status,
                                      struct internal_abs_eval_workspace *ws) {
    if (pt != NULL) {
        pt->core.R_H = R_H_saved;
    }
    free_internal_abs_eval_workspace(ws);
    if ((status < 0) && (comp != NULL)) {
        comp->is_valid = 0;
    }
    return status;
}

/*
 * Build the target seed-photon intensity field and expose the sampled DRF
 * arrays plus active frequency bounds for the requested component.
 */
static void build_seed_spectrum(struct blob *pt,
                                intabs_comp_t comp_id,
                                double **nu_grid,
                                double **n_grid,
                                double *nu_start,
                                double *nu_stop) {
    Build_I_nu_Disk(pt);

    if (comp_id == INTABS_COMP_BLR) {
        Build_I_nu_BLR(pt);
        *nu_grid = pt->BLR.spec.nu_DRF;
        *n_grid = pt->BLR.spec.n_nu_DRF;
        *nu_start = pt->BLR.spec.nu_min_DRF;
        *nu_stop = pt->BLR.spec.nu_max_DRF;
    } else if (comp_id == INTABS_COMP_DT) {
        Build_I_nu_DT(pt);
        *nu_grid = pt->DT.spec.nu_DRF;
        *n_grid = pt->DT.spec.n_nu_DRF;
        *nu_start = pt->DT.spec.nu_min_DRF;
        *nu_stop = pt->DT.spec.nu_max_DRF;
    } else {
        Build_I_nu_Corona(pt);
        *nu_grid = pt->Corona.spec.nu_DRF;
        *n_grid = pt->Corona.spec.n_nu_DRF;
        *nu_start = pt->Corona.spec.nu_min_DRF;
        *nu_stop = pt->Corona.spec.nu_max_DRF;
    }
}

/*
 * Sample the current seed field on a compact positive-frequency grid.
 *
 * Modes:
 * - `peak != 0`: return a one-point representation around the peak frequency
 *   with a normalization derived from the integrated spectrum.
 * - `peak == 0`: filter very weak tails, then log-resample to `N_soft` points
 *   using log-log interpolation.
 */
static int sample_seed_field(struct blob *pt,
                             intabs_comp_t comp_id,
                             unsigned int N_soft,
                             int peak,
                             double *nu_out,
                             double *n_out) {
    unsigned int i;
    unsigned int size_grid;
    unsigned int n_valid;
    unsigned int n_sel;
    unsigned int i_max;
    unsigned int left;
    double y_max;
    double threshold;
    double integral;
    double dnu;
    double log_x_min;
    double log_x_max;
    double log_step;
    double lx;
    double log_x0;
    double log_x1;
    double log_y0;
    double log_y1;
    double t;
    double *nu_grid;
    double *n_grid;
    double nu_start;
    double nu_stop;
    double *x = NULL;
    double *y = NULL;
    double *x_f = NULL;
    double *y_f = NULL;

    if ((pt == NULL) || (nu_out == NULL) || (n_out == NULL) || (N_soft == 0)) {
        return -1;
    }

    build_seed_spectrum(pt, comp_id, &nu_grid, &n_grid, &nu_start, &nu_stop);

    size_grid = pt->core.nu_grid_size;
    if (size_grid < 2) {
        return -1;
    }

    n_valid = 0;
    for (i = 0; i < size_grid; ++i) {
        if ((nu_grid[i] >= nu_start) && (nu_grid[i] <= nu_stop) && (nu_grid[i] > 0.0) &&
            (n_grid[i] > 0.0) && isfinite(n_grid[i])) {
            n_valid += 1U;
        }
    }
    if (n_valid < 2) {
        return -1;
    }

    x = (double *)malloc(((size_t)n_valid) * sizeof(double));
    y = (double *)malloc(((size_t)n_valid) * sizeof(double));
    if ((x == NULL) || (y == NULL)) {
        if (x != NULL) {
            free(x);
        }
        if (y != NULL) {
            free(y);
        }
        return -1;
    }

    n_sel = 0;
    y_max = 0.0;
    for (i = 0; i < size_grid; ++i) {
        if ((nu_grid[i] >= nu_start) && (nu_grid[i] <= nu_stop) && (nu_grid[i] > 0.0) &&
            (n_grid[i] > 0.0) && isfinite(n_grid[i])) {
            x[n_sel] = nu_grid[i];
            y[n_sel] = n_grid[i];
            if (y[n_sel] > y_max) {
                y_max = y[n_sel];
            }
            n_sel += 1U;
        }
    }

    if ((n_sel < 2) || (y_max <= 0.0)) {
        free(x);
        free(y);
        return -1;
    }

    if (peak) {
        i_max = 0;
        for (i = 1; i < n_sel; ++i) {
            if (y[i] > y[i_max]) {
                i_max = i;
            }
        }

        integral = 0.0;
        for (i = 1; i < n_sel; ++i) {
            dnu = x[i] - x[i - 1];
            integral += 0.5 * (y[i - 1] + y[i]) * dnu;
        }

        nu_out[0] = x[i_max];
        if ((integral > 0.0) && (x[i_max] > 0.0) && (y[i_max] > 0.0)) {
            n_out[0] = integral / x[i_max];
        } else {
            n_out[0] = y[i_max];
        }

        free(x);
        free(y);
        return 0;
    }

    threshold = y_max / 100.0;
    n_valid = 0;
    for (i = 0; i < n_sel; ++i) {
        if (y[i] > threshold) {
            n_valid += 1U;
        }
    }

    if (n_valid >= 2) {
        x_f = (double *)malloc(((size_t)n_valid) * sizeof(double));
        y_f = (double *)malloc(((size_t)n_valid) * sizeof(double));
        if ((x_f == NULL) || (y_f == NULL)) {
            if (x_f != NULL) {
                free(x_f);
            }
            if (y_f != NULL) {
                free(y_f);
            }
            free(x);
            free(y);
            return -1;
        }

        n_valid = 0;
        for (i = 0; i < n_sel; ++i) {
            if (y[i] > threshold) {
                x_f[n_valid] = x[i];
                y_f[n_valid] = y[i];
                n_valid += 1U;
            }
        }

        free(x);
        free(y);
        x = x_f;
        y = y_f;
        n_sel = n_valid;
    }

    if (n_sel < 2) {
        free(x);
        free(y);
        return -1;
    }

    if (N_soft == 1) {
        nu_out[0] = x[0];
        n_out[0] = y[0];
        free(x);
        free(y);
        return 0;
    }

    log_x_min = log10(x[0]);
    log_x_max = log10(x[n_sel - 1]);
    if ((!isfinite(log_x_min)) || (!isfinite(log_x_max)) || (log_x_max <= log_x_min)) {
        free(x);
        free(y);
        return -1;
    }

    log_step = (log_x_max - log_x_min) / ((double)N_soft - 1.0);
    left = 0;
    for (i = 0; i < N_soft; ++i) {
        lx = log_x_min + log_step * (double)i;
        nu_out[i] = pow(10.0, lx);

        while (((left + 1U) < n_sel) && (nu_out[i] > x[left + 1U])) {
            left += 1U;
        }

        if (left >= (n_sel - 1U)) {
            n_out[i] = INTABS_MIN_Y;
            continue;
        }

        log_x0 = log10(x[left]);
        log_x1 = log10(x[left + 1U]);
        log_y0 = log10(y[left]);
        log_y1 = log10(y[left + 1U]);

        if (log_x1 == log_x0) {
            n_out[i] = y[left];
        } else {
            t = (lx - log_x0) / (log_x1 - log_x0);
            n_out[i] = pow(10.0, log_y0 + t * (log_y1 - log_y0));
        }

        if ((!isfinite(n_out[i])) || (n_out[i] <= 0.0)) {
            n_out[i] = INTABS_MIN_Y;
        }
    }

    free(x);
    free(y);
    return 0;
}

/* Pair-production cross section sigma_{gamma-gamma}(s) for center-of-mass energy `s`. */
static double sigma_gamma_gamma(double s) {
    double beta;
    double beta2;
    double beta4;
    double term;

    if (s < 1.0) {
        return 0.0;
    }

    beta = sqrt(1.0 - 1.0 / s);
    beta2 = beta * beta;
    beta4 = beta2 * beta2;
    term = (3.0 - beta4) * log((1.0 + beta) / (1.0 - beta)) - 2.0 * beta * (2.0 - beta2);

    return (0.75 * SIGTH * 0.5) * (1.0 - beta2) * term;
}

/*
 * Precompute sigma_{gamma-gamma}(s) on a compact transformed grid:
 * u = (s - 1) / (s + GG_A),  s = (1 + GG_A*u) / (1 - u),  u in [0, 1).
 */
void init_sigma_gamma_gamma_table(struct blob *pt) {
    int i;
    double u;
    double s;
    struct internal_abs_store *store;

    if (pt == NULL) {
        return;
    }

    store = &(pt->core.internal_abs);
    for (i = 0; i <= GG_NTAB; ++i) {
        u = (double)i / (double)GG_NTAB;

        if ((i == 0) || (i == GG_NTAB)) {
            store->gg_tab[i] = 0.0;
        } else {
            s = (1.0 + GG_A * u) / (1.0 - u);
            store->gg_tab[i] = sigma_gamma_gamma(s);
        }
    }
}

/* Fast linear-interpolated sigma_{gamma-gamma}(s) lookup on the precomputed table. */
static inline double sigma_gamma_gamma_fast(const struct internal_abs_store *store, double s) {
    double u;
    double x;
    double f;
    int i;

    if ((store == NULL) || (s < 1.0)) {
        return 0.0;
    }

    u = (s - 1.0) / (s + GG_A);
    if (u <= 0.0) {
        return 0.0;
    }
    if (u >= 1.0) {
        return 0.0;
    }

    x = u * GG_NTAB;
    i = (int)x;
    if (i >= GG_NTAB) {
        return 0.0;
    }

    f = x - (double)i;
    return store->gg_tab[i] * (1.0 - f) + store->gg_tab[i + 1] * f;
}

/* Replace a zero/undefined saved grid size with a conservative default value. */
static unsigned int sanitize_grid_size(unsigned int value, unsigned int fallback_value) {
    if (value == 0U) {
        return fallback_value;
    }
    return value;
}

/*
 * Interpolate one component tau table at `nu_obs` in log-log space.
 * Returns a tiny positive floor outside the low-energy side to avoid zeros.
 */
static double interp_tau_component(const struct internal_abs_component *comp, double nu_obs) {
    const double eps_tau = 1.0e-300;
    unsigned int ID;
    unsigned int i_max;
    double tau1;
    double tau2;
    double nu1;
    double nu2;
    double log_tau;
    double t;

    if ((comp == NULL) || (comp->is_enabled == 0) || (comp->is_valid == 0) || (comp->tau_size == 0U) ||
        (comp->nu_tau == NULL) || (comp->tau == NULL) || (nu_obs <= 0.0) || (!isfinite(nu_obs))) {
        return 0.0;
    }

    i_max = comp->tau_size - 1U;
    if (nu_obs < comp->nu_tau[0]) {
        return eps_tau;
    }
    if (nu_obs >= comp->nu_tau[i_max]) {
        return (comp->tau[i_max] > 0.0) ? comp->tau[i_max] : eps_tau;
    }

    ID = (unsigned int)x_to_grid_index(comp->nu_tau, nu_obs, comp->tau_size);
    if (ID >= i_max) {
        return eps_tau;
    }

    nu1 = comp->nu_tau[ID];
    nu2 = comp->nu_tau[ID + 1U];
    if ((nu1 <= 0.0) || (nu2 <= nu1)) {
        return eps_tau;
    }

    tau1 = (comp->tau[ID] > eps_tau) ? comp->tau[ID] : eps_tau;
    tau2 = (comp->tau[ID + 1U] > eps_tau) ? comp->tau[ID + 1U] : eps_tau;
    t = (log10(nu_obs) - log10(nu1)) / (log10(nu2) - log10(nu1));
    log_tau = log10(tau1) + t * (log10(tau2) - log10(tau1));

    return pow(10.0, log_tau);
}

/*
 * Core internal-absorption integration routine for one seed component.
 *
 * Inputs control the quadrature grids:
 * - `N_soft`: soft-photon frequency samples (or one sample in peak mode),
 * - `N_hard`: number of gamma-ray frequencies where tau is stored,
 * - `N_R_H`: samples along propagation distance,
 * - `N_theta`: angular samples in cos(theta)=mu.
 *
 * Integration flow:
 * 1) validate arguments, select component, and prepare work arrays;
 * 2) build a reference soft field and derive `nu_min` if needed;
 * 3) build gamma-ray grid `comp->nu_tau` from `nu_min` to `nu_src_max`;
 * 4) for each R_H sample, obtain seed spectra (resampled or extrapolated),
 *    convert to dimensionless energies, and build mu-grid geometry;
 * 5) for each gamma frequency, integrate trapezoidally over
 *    nu_soft -> mu -> R_H, using `sigma_gamma_gamma(s)`;
 * 6) store tau and cache setup into `pt->core.internal_abs.<component>`.
 *
 * During evaluation `pt->core.R_H` is temporarily changed for sampling and
 * restored on every return path by `finalize_internal_abs_eval`.
 */
int eval_internal_abs_tau(struct blob *pt,
                          const char *seed_photons_name,
                          double nu_min,
                          unsigned int N_soft,
                          unsigned int N_hard,
                          unsigned int N_R_H,
                          unsigned int N_theta,
                          int use_R_H_profile_extrapolation,
                          int peak,
                          double nu_src_max) {
    int status;
    intabs_comp_t comp_id;
    struct internal_abs_component *comp;
    double R_H_saved_input;
    double R_H_saved;
    double R_H_saved_eval;
    double R_H_sample;
    double distance_blob_from_seed_field_geom_center;
    double R_seed;
    double R_H_ref;
    double corona_side;
    double nu_min_eff;
    double nu_soft_max;
    double ratio;
    double mu_min;
    double mu_max;
    double one_minus_mu;
    double eps_gamma;
    double s_value;
    double integrand;
    double prev_integrand;
    double nu_integral;
    double prev_mu_integral;
    double mu_integral;
    double prev_rh_integral;
    double tau_gamma;
    double R_x;
    double scale;
    double nu_src_max_eff;
    double nu_to_eps;
    double soft_scale;
    double sigma_val;
    unsigned int N_soft_eff;
    unsigned int tau_size;
    unsigned int i;
    unsigned int ID_RH;
    unsigned int ID_THETA;
    unsigned int ID_SOFT;
    unsigned int ID_GAMMA;
    size_t rh_soft_base;
    size_t rh_mu_base;
    size_t idx_soft;
    size_t idx_mu;
    const struct internal_abs_store *ia_store;
    int use_fast_sigma;
    struct internal_abs_eval_workspace ws;

    status = -1;
    init_internal_abs_eval_workspace(&ws);

    /* Basic pointer/component checks. */
    if (pt == NULL) {
        return -1;
    }

    comp_id = parse_internal_abs_component(seed_photons_name);
    comp = get_internal_abs_component_ptr(pt, comp_id);
    if ((comp_id == INTABS_COMP_INVALID) || (comp == NULL)) {
        return -1;
    }

    /* All integration dimensions must be strictly positive. */
    if ((N_soft == 0U) || (N_hard == 0U) || (N_R_H == 0U) || (N_theta == 0U)) {
        return -1;
    }

    ia_store = &(pt->core.internal_abs);
    use_fast_sigma = (comp->use_sigma_gamma_gamma_fast != 0) ? 1 : 0;

    /* In peak mode the soft field is collapsed to one representative point. */
    N_soft_eff = (peak != 0) ? 1U : N_soft;
    if (N_soft_eff == 0U) {
        N_soft_eff = 1U;
    }

    /* Save original position; keep a positive fallback for geometric scales. */
    R_H_saved_input = pt->core.R_H;
    R_H_saved = R_H_saved_input;
    if (R_H_saved <= 0.0) {
        R_H_saved = 1.0;
    }

    /* Select characteristic size of the target photon field. */
    if (comp_id == INTABS_COMP_BLR) {
        R_seed = pt->BLR.R_BLR_out;
    } else if (comp_id == INTABS_COMP_DT) {
        R_seed = pt->DT.R_DT;
    } else {
        // better working with R_corona for disk geometry
        R_seed = pt->Corona.R_Corona;
    }
    if (R_seed <= 0.0) {
        R_seed = R_H_saved;
    }

    /*
     * Corona integration is done in distance from the corona center and keeps
     * track of which side of the corona the blob is located on.
     */
    if (comp_id == INTABS_COMP_CORONA) {
        R_H_saved_eval = fabs(R_H_saved - pt->Corona.R_H_Corona);
        if (R_H_saved_input >= pt->Corona.R_H_Corona) {
            corona_side = 1.0;
        } else {
            //between corona and BH
            corona_side = -1.0;
        }
    } else {
        R_H_saved_eval = R_H_saved;
        corona_side = 1.0;
    }
    if (R_H_saved_eval <= 0.0) {
        R_H_saved_eval = 1.0;
    }

    /* Allocate all temporary grids used by the 3D trapezoidal integration. */
    ws.nu_soft = (double *)calloc((size_t)N_R_H * (size_t)N_soft_eff, sizeof(double));
    ws.n_soft = (double *)calloc((size_t)N_R_H * (size_t)N_soft_eff, sizeof(double));
    ws.eps_soft = (double *)calloc((size_t)N_R_H * (size_t)N_soft_eff, sizeof(double));
    ws.dnu_soft = (double *)calloc((size_t)N_R_H * (size_t)N_soft_eff, sizeof(double));
    ws.mu_grid = (double *)calloc((size_t)N_R_H * (size_t)N_theta, sizeof(double));
    ws.one_minus_mu = (double *)calloc((size_t)N_R_H * (size_t)N_theta, sizeof(double));
    ws.dmu_grid = (double *)calloc((size_t)N_R_H * (size_t)N_theta, sizeof(double));
    ws.R_H_grid = (double *)calloc((size_t)N_R_H, sizeof(double));
    ws.d_rh_grid = (double *)calloc((size_t)N_R_H, sizeof(double));
    ws.nu_soft_ref = (double *)calloc((size_t)N_soft_eff, sizeof(double));
    ws.n_soft_ref = (double *)calloc((size_t)N_soft_eff, sizeof(double));
    ws.nu_soft_tmp = (double *)calloc((size_t)N_soft_eff, sizeof(double));
    ws.n_soft_tmp = (double *)calloc((size_t)N_soft_eff, sizeof(double));

    if ((ws.nu_soft == NULL) || (ws.n_soft == NULL) || (ws.eps_soft == NULL) || (ws.dnu_soft == NULL) ||
        (ws.mu_grid == NULL) || (ws.one_minus_mu == NULL) || (ws.dmu_grid == NULL) ||
        (ws.R_H_grid == NULL) || (ws.d_rh_grid == NULL) || (ws.nu_soft_ref == NULL) ||
        (ws.n_soft_ref == NULL) || (ws.nu_soft_tmp == NULL) || (ws.n_soft_tmp == NULL)) {
        return finalize_internal_abs_eval(pt, comp, R_H_saved_input, -1, &ws);
    }

    /* Build one reference seed spectrum close to the source field scale. */
    distance_blob_from_seed_field_geom_center = R_seed/1000;
    if (comp_id == INTABS_COMP_CORONA) {
        /*
         * In corona mode, this quantity is interpreted as distance from the
         * corona center along the jet axis, not as absolute R_H.
         */
        distance_blob_from_seed_field_geom_center = R_seed;
        /*
         * Convert center-relative distance to absolute jet coordinate:
         *   R_H = R_H_corona_center +/- distance_from_center
         * `corona_side` preserves which side of the corona the blob is on.
         */
        R_H_sample = pt->Corona.R_H_Corona + corona_side * distance_blob_from_seed_field_geom_center;
        pt->core.R_H = fmax(R_H_sample, 0.0);
    } else {
        /* BLR/DT path: integration coordinate is already absolute R_H. */
        pt->core.R_H = distance_blob_from_seed_field_geom_center;
    }
    if (sample_seed_field(pt, comp_id, N_soft_eff, peak, ws.nu_soft_ref, ws.n_soft_ref) < 0) {
        return finalize_internal_abs_eval(pt, comp, R_H_saved_input, -1, &ws);
    }

    /* Conversion nu -> dimensionless epsilon = h nu / (m_e c^2). */
    nu_to_eps = HPLANCK / MEC2;

    /*
     * If caller did not provide nu_min, estimate it from the highest soft
     * frequency so the pair-production threshold can be reached.
     */
    if (nu_min > 0.0) {
        nu_min_eff = nu_min;
    } else {
        nu_soft_max = ws.nu_soft_ref[0];
        for (i = 1; i < N_soft_eff; ++i) {
            if (ws.nu_soft_ref[i] > nu_soft_max) {
                nu_soft_max = ws.nu_soft_ref[i];
            }
        }
        if (nu_soft_max > 0.0) {
            nu_min_eff = 1.0e40 / nu_soft_max;
        } else {
            nu_min_eff = 1.0e20;
        }
    }
    if (nu_min_eff <= 0.0) {
        nu_min_eff = 1.0e20;
    }

    /* Ensure hard-photon range is valid; collapse to one point if inverted. */
    nu_src_max_eff = nu_src_max;
    if (nu_src_max_eff <= 0.0) {
        nu_src_max_eff = nu_min_eff;
    }

    if (nu_src_max_eff < nu_min_eff) {
        tau_size = 1U;
    } else {
        tau_size = N_hard;
    }

    /* Allocate/resize component output arrays nu_tau and tau. */
    if (ensure_tau_arrays(comp, tau_size) < 0) {
        return finalize_internal_abs_eval(pt, comp, R_H_saved_input, -1, &ws);
    }

    /* Build the hard-photon grid where tau will be stored (log-spaced). */
    if (tau_size == 1U) {
        comp->nu_tau[0] = nu_min_eff;
    } else {
        for (i = 0; i < tau_size; ++i) {
            comp->nu_tau[i] = pow(10.0,
                                  log10(nu_min_eff) +
                                      (log10(nu_src_max_eff) - log10(nu_min_eff)) * (double)i / (double)(tau_size - 1U));
        }
    }

    /* Build propagation-distance grid (R_H) over up to 3 decades. */
    R_H_ref = R_H_saved_eval;
    if (N_R_H == 1U) {
        ws.R_H_grid[0] = R_H_ref;
    } else {
        for (ID_RH = 0; ID_RH < N_R_H; ++ID_RH) {
            ws.R_H_grid[ID_RH] = pow(10.0, 3.0 * (double)ID_RH / (double)(N_R_H - 1U)) * R_H_ref;
        }
    }

    /*
     * For each R_H:
     * - get soft field (either full resampling or 1/R^2 extrapolation),
     * - convert soft photons to epsilon and dnu bins,
     * - build angular grid limits from source geometry.
     */
    double denom,mu;
    for (ID_RH = 0; ID_RH < N_R_H; ++ID_RH) {
        distance_blob_from_seed_field_geom_center = ws.R_H_grid[ID_RH];
        if (comp_id == INTABS_COMP_CORONA) {
            /*
             * Corona grid is built in distance-from-center. Map each sample to
             * absolute R_H before evaluating geometry-dependent seed fields.
             */
            R_H_sample = pt->Corona.R_H_Corona + corona_side * distance_blob_from_seed_field_geom_center;
            pt->core.R_H = fmax(R_H_sample, 0.0);
        } else {
            pt->core.R_H = distance_blob_from_seed_field_geom_center;
        }
        rh_soft_base = ((size_t)ID_RH) * ((size_t)N_soft_eff);

        if (use_R_H_profile_extrapolation != 0) {

            if (comp_id != INTABS_COMP_CORONA){
                if (distance_blob_from_seed_field_geom_center <= R_seed) {
                    scale = 1.0;
                } else {
                    denom = sqrt(distance_blob_from_seed_field_geom_center * distance_blob_from_seed_field_geom_center + R_seed *R_seed);
                    if (denom > 0.0){
                        mu = distance_blob_from_seed_field_geom_center / denom;
                    }
                    else{
                        mu = 0.0;
                    }
                    scale = (1-mu);
                    }
            }else{
                denom = sqrt(distance_blob_from_seed_field_geom_center * distance_blob_from_seed_field_geom_center + R_seed *R_seed);
                if (denom > 0.0){
                    mu = distance_blob_from_seed_field_geom_center / denom;
                }
                else{
                    mu = 0.0;
                }
                scale = (1-mu)*pi;
            }
            

            for (ID_SOFT = 0; ID_SOFT < N_soft_eff; ++ID_SOFT) {
                idx_soft = rh_soft_base + (size_t)ID_SOFT;
                ws.nu_soft[idx_soft] = ws.nu_soft_ref[ID_SOFT];
                ws.n_soft[idx_soft] = ws.n_soft_ref[ID_SOFT] * scale;
            }
        } else {
            if (sample_seed_field(pt, comp_id, N_soft_eff, peak, ws.nu_soft_tmp, ws.n_soft_tmp) < 0) {
                return finalize_internal_abs_eval(pt, comp, R_H_saved_input, -1, &ws);
            }

            for (ID_SOFT = 0; ID_SOFT < N_soft_eff; ++ID_SOFT) {
                idx_soft = rh_soft_base + (size_t)ID_SOFT;
                ws.nu_soft[idx_soft] = ws.nu_soft_tmp[ID_SOFT];
                ws.n_soft[idx_soft] = ws.n_soft_tmp[ID_SOFT];
            }
        }

        for (ID_SOFT = 0; ID_SOFT < N_soft_eff; ++ID_SOFT) {
            idx_soft = rh_soft_base + (size_t)ID_SOFT;
            ws.eps_soft[idx_soft] = ws.nu_soft[idx_soft] * nu_to_eps;
            if (ID_SOFT == 0U) {
                ws.dnu_soft[idx_soft] = 0.0;
            } else {
                ws.dnu_soft[idx_soft] = ws.nu_soft[idx_soft] - ws.nu_soft[idx_soft - 1U];
            }
        }

        mu_max = 1.0;
        if (distance_blob_from_seed_field_geom_center < R_seed) {
            mu_min = -1.0;
        } else {
            ratio = R_seed / distance_blob_from_seed_field_geom_center;
            if (ratio > 1.0) {
                ratio = 1.0;
            }
            if (ratio < 0.0) {
                ratio = 0.0;
            }
            mu_min = sqrt(1.0 - ratio * ratio);
        }

        rh_mu_base = ((size_t)ID_RH) * ((size_t)N_theta);
        if (N_theta == 1U) {
            ws.mu_grid[rh_mu_base] = mu_min;
            ws.one_minus_mu[rh_mu_base] = fmax(1.0 - mu_min, INTABS_MIN_ONE_MINUS_MU);
            ws.dmu_grid[rh_mu_base] = 0.0;
        } else {
            for (ID_THETA = 0; ID_THETA < N_theta; ++ID_THETA) {
                idx_mu = rh_mu_base + (size_t)ID_THETA;
                ws.mu_grid[idx_mu] = mu_min + (mu_max - mu_min) * (double)ID_THETA / (double)(N_theta - 1U);
                ws.one_minus_mu[idx_mu] = fmax(1.0 - ws.mu_grid[idx_mu], INTABS_MIN_ONE_MINUS_MU);
                if (ID_THETA == 0U) {
                    ws.dmu_grid[idx_mu] = 0.0;
                } else {
                    ws.dmu_grid[idx_mu] = ws.mu_grid[idx_mu] - ws.mu_grid[idx_mu - 1U];
                }
            }
        }
    }

    /* Precompute R_H bin widths for trapezoidal integration along the path. */
    ws.d_rh_grid[0] = 0.0;
    for (ID_RH = 1; ID_RH < N_R_H; ++ID_RH) {
        ws.d_rh_grid[ID_RH] = ws.R_H_grid[ID_RH] - ws.R_H_grid[ID_RH - 1U];
    }

    /*
     * Nested integration order for each gamma energy:
     * 1) integrate over soft frequency,
     * 2) integrate over angle mu,
     * 3) integrate over path length R_H.
     */
    for (ID_GAMMA = 0; ID_GAMMA < tau_size; ++ID_GAMMA) {
        eps_gamma = comp->nu_tau[ID_GAMMA] * nu_to_eps;
        tau_gamma = 0.0;
        prev_rh_integral = 0.0;

        for (ID_RH = 0; ID_RH < N_R_H; ++ID_RH) {
            rh_mu_base = ((size_t)ID_RH) * ((size_t)N_theta);
            rh_soft_base = ((size_t)ID_RH) * ((size_t)N_soft_eff);
            mu_integral = 0.0;
            prev_mu_integral = 0.0;

            for (ID_THETA = 0; ID_THETA < N_theta; ++ID_THETA) {
                idx_mu = rh_mu_base + (size_t)ID_THETA;
                one_minus_mu = ws.one_minus_mu[idx_mu];
                soft_scale = eps_gamma * one_minus_mu * 0.5;

                nu_integral = 0.0;
                prev_integrand = 0.0;
                for (ID_SOFT = 0; ID_SOFT < N_soft_eff; ++ID_SOFT) {
                    idx_soft = rh_soft_base + (size_t)ID_SOFT;
                    s_value = soft_scale * ws.eps_soft[idx_soft];
                    integrand = 0.0;
                    if (s_value >= 1.0) {
                        if (use_fast_sigma != 0) {
                            sigma_val = sigma_gamma_gamma_fast(ia_store, s_value);
                        } else {
                            sigma_val = sigma_gamma_gamma(s_value);
                        }
                        integrand = sigma_val * ws.n_soft[idx_soft] * one_minus_mu;
                    }

                    if (ID_SOFT > 0U) {
                        nu_integral += 0.5 * (prev_integrand + integrand) * ws.dnu_soft[idx_soft];
                    }
                    prev_integrand = integrand;
                }

                if (ID_THETA > 0U) {
                    mu_integral += 0.5 * (prev_mu_integral + nu_integral) * ws.dmu_grid[idx_mu];
                }
                prev_mu_integral = nu_integral;
            }

            if (ID_RH > 0U) {
                tau_gamma += 0.5 * (prev_rh_integral + mu_integral) * ws.d_rh_grid[ID_RH];
            }
            prev_rh_integral = mu_integral;
        }

        comp->tau[ID_GAMMA] = 2.0 * pi * tau_gamma;
    }

    /* Persist computed tau table and the settings used to produce it. */
    comp->is_enabled = 1;
    comp->is_valid = 1;
    comp->use_R_H_profile_extrapolation = use_R_H_profile_extrapolation;
    comp->peak_mode = peak;
    comp->N_soft = N_soft_eff;
    comp->N_hard = N_hard;
    comp->N_R_H = N_R_H;
    comp->N_theta = N_theta;
    comp->nu_min = nu_min_eff;
    comp->nu_src_max = nu_src_max_eff;

    status = (int)tau_size;
    return finalize_internal_abs_eval(pt, comp, R_H_saved_input, status, &ws);
}

/*
 * Recompute tau tables for all components currently marked as enabled,
 * preserving per-component numerical settings whenever available.
 */
void recompute_internal_absorption_tau(struct blob *pt) {
    struct internal_abs_component *comp_blr;
    struct internal_abs_component *comp_dt;
    struct internal_abs_component *comp_corona;
    double nu_src_max;
    int status;

    if (pt == NULL) {
        return;
    }

    nu_src_max = pt->core.nu_stop_grid;
    if (nu_src_max <= 0.0) {
        nu_src_max = 1.0e30;
    }

    comp_blr = &(pt->core.internal_abs.BLR);
    if (comp_blr->is_enabled != 0) {
        status = eval_internal_abs_tau(pt,
                                       "BLR",
                                       comp_blr->nu_min,
                                       sanitize_grid_size(comp_blr->N_soft, 50U),
                                       sanitize_grid_size(comp_blr->N_hard, 50U),
                                       sanitize_grid_size(comp_blr->N_R_H, 50U),
                                       sanitize_grid_size(comp_blr->N_theta, 50U),
                                       comp_blr->use_R_H_profile_extrapolation,
                                       comp_blr->peak_mode,
                                       nu_src_max);
        if (status < 0) {
            comp_blr->is_valid = 0;
        }
    }

    comp_dt = &(pt->core.internal_abs.DT);
    if (comp_dt->is_enabled != 0) {
        status = eval_internal_abs_tau(pt,
                                       "DT",
                                       comp_dt->nu_min,
                                       sanitize_grid_size(comp_dt->N_soft, 50U),
                                       sanitize_grid_size(comp_dt->N_hard, 50U),
                                       sanitize_grid_size(comp_dt->N_R_H, 50U),
                                       sanitize_grid_size(comp_dt->N_theta, 50U),
                                       comp_dt->use_R_H_profile_extrapolation,
                                       comp_dt->peak_mode,
                                       nu_src_max);
        if (status < 0) {
            comp_dt->is_valid = 0;
        }
    }

    comp_corona = &(pt->core.internal_abs.Corona);
    if (comp_corona->is_enabled != 0) {
        status = eval_internal_abs_tau(pt,
                                       "Corona",
                                       comp_corona->nu_min,
                                       sanitize_grid_size(comp_corona->N_soft, 50U),
                                       sanitize_grid_size(comp_corona->N_hard, 50U),
                                       sanitize_grid_size(comp_corona->N_R_H, 50U),
                                       sanitize_grid_size(comp_corona->N_theta, 50U),
                                       comp_corona->use_R_H_profile_extrapolation,
                                       comp_corona->peak_mode,
                                       nu_src_max);
        if (status < 0) {
            comp_corona->is_valid = 0;
        }
    }
}

/*
 * Return total internal opacity at `nu_obs` by summing BLR/DT/Corona
 * interpolated tau contributions. Invalid/non-finite totals are clamped to 0.
 */
double get_internal_abs_tau_at_nu(struct blob *pt, double nu_obs) {
    double tau_tot;

    if (pt == NULL) {
        return 0.0;
    }

    tau_tot = 0.0;
    tau_tot += interp_tau_component(&(pt->core.internal_abs.BLR), nu_obs);
    tau_tot += interp_tau_component(&(pt->core.internal_abs.DT), nu_obs);
    tau_tot += interp_tau_component(&(pt->core.internal_abs.Corona), nu_obs);

    if (!isfinite(tau_tot) || (tau_tot < 0.0)) {
        return 0.0;
    }

    return tau_tot;
}
