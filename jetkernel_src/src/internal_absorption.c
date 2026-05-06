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
#define INTABS_PEAK_SAMPLES 1U

typedef enum {
    INTABS_COMP_INVALID = 0,
    INTABS_COMP_BLR = 1,
    INTABS_COMP_DT = 2,
    INTABS_COMP_CORONA = 3,
    INTABS_COMP_TOTAL = 4
} intabs_comp_t;

/* Per-evaluation guard: build disk seed spectrum once on first sampling call. */
static int intabs_disk_seed_built = 0;

static void reset_intabs_seed_build_guard(void) {
    intabs_disk_seed_built = 0;
}

static void ensure_intabs_disk_seed_built(struct blob *pt_cloned) {
    if ((pt_cloned != NULL) && (intabs_disk_seed_built == 0)) {
        Build_I_nu_Disk(pt_cloned);
        intabs_disk_seed_built = 1;
    }
}

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
    init_internal_abs_component(&(pt->core.internal_abs.Total));
}

/* Release all per-component tau caches stored in the blob. */
void free_internal_abs_store(struct blob *pt) {
    if (pt == NULL) {
        return;
    }
    free_internal_abs_component(&(pt->core.internal_abs.BLR));
    free_internal_abs_component(&(pt->core.internal_abs.DT));
    free_internal_abs_component(&(pt->core.internal_abs.Corona));
    free_internal_abs_component(&(pt->core.internal_abs.Total));
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
static struct internal_abs_component *get_internal_abs_component_ptr(struct blob *pt_cloned, intabs_comp_t comp_id) {
    if (pt_cloned == NULL) {
        return NULL;
    }

    if (comp_id == INTABS_COMP_BLR) {
        return &(pt_cloned->core.internal_abs.BLR);
    }
    if (comp_id == INTABS_COMP_DT) {
        return &(pt_cloned->core.internal_abs.DT);
    }
    if (comp_id == INTABS_COMP_CORONA) {
        return &(pt_cloned->core.internal_abs.Corona);
    }
    if (comp_id == INTABS_COMP_TOTAL) {
        return &(pt_cloned->core.internal_abs.Total);
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
 * Shared solver exit path: restore `pt_cloned->core.R_H`, free temporaries,
 * and invalidate component cache on errors.
 */
static int finalize_internal_abs_eval(struct blob *pt_cloned,
                                      struct internal_abs_component *comp,
                                      double R_H_saved,
                                      int status,
                                      struct internal_abs_eval_workspace *ws) {
    if (pt_cloned != NULL) {
        pt_cloned->core.R_H = R_H_saved;
    }
    free_internal_abs_eval_workspace(ws);
    if ((status < 0) && (comp != NULL)) {
        comp->is_valid = 0;
    }
    return status;
}

static int finalize_internal_abs_total_eval(struct blob *pt_cloned,
                                            struct internal_abs_component *comp_tot,
                                            double R_H_saved_input,
                                            int status,
                                            double *nu_soft_ref[3],
                                            double *n_soft_ref[3],
                                            double *nu_soft_tmp[3],
                                            double *n_soft_tmp[3],
                                            double *nu_soft_common,
                                            double *eps_soft_common,
                                            double *dnu_soft_common,
                                            double *x_grid,
                                            double *dx_grid,
                                            double *mu_grid,
                                            double *dmu_grid,
                                            double *one_minus_mu_grid,
                                            double *sum_soft_mu,
                                            double *comp_soft_interp,
                                            double *comp_mu_min_grid,
                                            double *tmp_resampled) {
    unsigned int c;

    if (pt_cloned != NULL) {
        pt_cloned->core.R_H = R_H_saved_input;
    }

    if ((status < 0) && (comp_tot != NULL)) {
        comp_tot->is_valid = 0;
        comp_tot->is_enabled = 0;
    }

    for (c = 0U; c < 3U; ++c) {
        if (nu_soft_ref[c] != NULL) {
            free(nu_soft_ref[c]);
        }
        if (n_soft_ref[c] != NULL) {
            free(n_soft_ref[c]);
        }
        if (nu_soft_tmp[c] != NULL) {
            free(nu_soft_tmp[c]);
        }
        if (n_soft_tmp[c] != NULL) {
            free(n_soft_tmp[c]);
        }
    }

    if (nu_soft_common != NULL) {
        free(nu_soft_common);
    }
    if (eps_soft_common != NULL) {
        free(eps_soft_common);
    }
    if (dnu_soft_common != NULL) {
        free(dnu_soft_common);
    }
    if (x_grid != NULL) {
        free(x_grid);
    }
    if (dx_grid != NULL) {
        free(dx_grid);
    }
    if (mu_grid != NULL) {
        free(mu_grid);
    }
    if (dmu_grid != NULL) {
        free(dmu_grid);
    }
    if (one_minus_mu_grid != NULL) {
        free(one_minus_mu_grid);
    }
    if (sum_soft_mu != NULL) {
        free(sum_soft_mu);
    }
    if (comp_soft_interp != NULL) {
        free(comp_soft_interp);
    }
    if (comp_mu_min_grid != NULL) {
        free(comp_mu_min_grid);
    }
    if (tmp_resampled != NULL) {
        free(tmp_resampled);
    }

    return status;
}

/*
 * Build the target seed-photon intensity field and expose the sampled DRF
 * arrays plus active frequency bounds for the requested component.
 */
static void build_seed_spectrum(struct blob *pt_cloned,
                                intabs_comp_t comp_id,
                                double **nu_grid,
                                double **n_grid,
                                double *nu_start,
                                double *nu_stop) {

    //pt_cloned->core.theta_n_int=10;
    //pt_cloned->core.l_n_int=10;
    ensure_intabs_disk_seed_built(pt_cloned);
    
    if (comp_id == INTABS_COMP_BLR) {
        
        Build_I_nu_BLR(pt_cloned);
        *nu_grid = pt_cloned->BLR.spec.nu_DRF;
        *n_grid = pt_cloned->BLR.spec.n_nu_DRF;
        *nu_start = pt_cloned->BLR.spec.nu_min_DRF;
        *nu_stop = pt_cloned->BLR.spec.nu_max_DRF;
    } else if (comp_id == INTABS_COMP_DT) {
        Build_I_nu_DT(pt_cloned);
        *nu_grid = pt_cloned->DT.spec.nu_DRF;
        *n_grid = pt_cloned->DT.spec.n_nu_DRF;
        *nu_start = pt_cloned->DT.spec.nu_min_DRF;
        *nu_stop = pt_cloned->DT.spec.nu_max_DRF;
    } else {
        Build_I_nu_Corona(pt_cloned);
        *nu_grid = pt_cloned->Corona.spec.nu_DRF;
        *n_grid = pt_cloned->Corona.spec.n_nu_DRF;
        *nu_start = pt_cloned->Corona.spec.nu_min_DRF;
        *nu_stop = pt_cloned->Corona.spec.nu_max_DRF;
    }
}

/*
 * Sample the current seed field on a compact positive-frequency grid.
 *
 * Modes:
 * - `peak != 0`: return a one-point representation at the seed-field peak
 *   frequency with normalization derived from the integrated spectrum.
 * - `peak == 0`: filter very weak tails, then log-resample to `N_soft` points
 *   using log-log interpolation.
 */
static int sample_seed_field(struct blob *pt_cloned,
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

    if ((pt_cloned == NULL) || (nu_out == NULL) || (n_out == NULL) || (N_soft == 0)) {
        return -1;
    }

    build_seed_spectrum(pt_cloned, comp_id, &nu_grid, &n_grid, &nu_start, &nu_stop);

    size_grid = pt_cloned->core.nu_grid_size;
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
        i_max = 0U;
        for (i = 1U; i < n_sel; ++i) {
            if (y[i] > y[i_max]) {
                i_max = i;
            }
        }

        integral = 0.0;
        for (i = 1U; i < n_sel; ++i) {
            dnu = x[i] - x[i - 1U];
            integral += 0.5 * (y[i - 1U] + y[i]) * dnu;
        }

        nu_out[0] = x[i_max];
        if ((integral > 0.0) && (x[i_max] > 0.0) && (y[i_max] > 0.0)) {
            n_out[0] = integral / x[i_max];
        } else {
            n_out[0] = y[i_max];
        }
        if ((!isfinite(n_out[0])) || (n_out[0] <= 0.0)) {
            n_out[0] = INTABS_MIN_Y;
        }
        for (i = 1U; i < N_soft; ++i) {
            nu_out[i] = nu_out[0];
            n_out[i] = INTABS_MIN_Y;
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

static double get_component_seed_radius(const struct blob *pt_cloned, intabs_comp_t comp_id, double fallback_radius) {
    double radius;

    radius = fallback_radius;
    if ((pt_cloned != NULL) && (comp_id == INTABS_COMP_BLR)) {
        radius =  pt_cloned->BLR.R_BLR_in;
    } else if ((pt_cloned != NULL) && (comp_id == INTABS_COMP_DT)) {
        radius = pt_cloned->DT.R_DT;
    } else if ((pt_cloned != NULL) && (comp_id == INTABS_COMP_CORONA)) {
        radius = pt_cloned->Corona.R_Corona;
    }

    if (radius <= 0.0) {
        radius = fallback_radius;
    }
    if (radius <= 0.0) {
        radius = 1.0;
    }

    return radius;
}

struct internal_abs_geometry_ctx {
    double R_seed;
    double R_H_ref;
    double corona_side;
};

static void resolve_component_geometry_ctx(const struct blob *pt_cloned,
                                           intabs_comp_t comp_id,
                                           double R_H_saved_input,
                                           double R_H_saved_fallback,
                                           struct internal_abs_geometry_ctx *ctx) {
    if (ctx == NULL) {
        return;
    }

    ctx->R_seed = get_component_seed_radius(pt_cloned, comp_id, R_H_saved_fallback);

    if ((pt_cloned != NULL) && (comp_id == INTABS_COMP_CORONA)) {
        ctx->R_H_ref = fabs(R_H_saved_fallback - pt_cloned->Corona.R_H_Corona);
        if (R_H_saved_input >= pt_cloned->Corona.R_H_Corona) {
            ctx->corona_side = 1.0;
        } else {
            ctx->corona_side = -1.0;
        }
    } else {
        ctx->R_H_ref = R_H_saved_fallback;
        ctx->corona_side = 1.0;
    }

    if (ctx->R_H_ref <= 0.0) {
        ctx->R_H_ref = 1.0;
    }
}

static double get_seed_reference_distance(intabs_comp_t comp_id, const struct internal_abs_geometry_ctx *ctx) {
    double distance_ref;

    if (ctx == NULL) {
        return 1.0;
    }

    if (comp_id == INTABS_COMP_CORONA) {
        distance_ref = ctx->R_seed;
    } else {
        distance_ref = 1;
    }

    if (distance_ref <= 0.0) {
        distance_ref = ctx->R_H_ref;
    }
    if (distance_ref <= 0.0) {
        distance_ref = 1.0;
    }

    return distance_ref;
}

static void set_component_sampling_position(struct blob *pt_cloned,
                                            intabs_comp_t comp_id,
                                            double distance_from_center,
                                            double corona_side) {
    double R_H_sample;

    if (pt_cloned == NULL) {
        return;
    }

    if (comp_id == INTABS_COMP_CORONA) {
        R_H_sample = pt_cloned->Corona.R_H_Corona + corona_side * distance_from_center;
        pt_cloned->core.R_H = fmax(R_H_sample, 0.0);
    } else {
        pt_cloned->core.R_H = distance_from_center;
    }
}


static double compute_seed_scale_extrapolated(const struct blob *pt_cloned, intabs_comp_t comp_id, double distance_from_center, double R_seed) {
    double denom;
    double mu;
    double scale;
    double R_BLR_eff,BLR_corr,Delta_BLR;
    if ((distance_from_center <= 0.0) || (R_seed <= 0.0)) {
        return 1.0;
    }

    if (comp_id == INTABS_COMP_DT) {
        if (distance_from_center <= R_seed) {
            return 1.0;
        }

        denom = sqrt(distance_from_center * distance_from_center + R_seed * R_seed);
        if (denom > 0.0) {
            mu = distance_from_center / denom;
        } else {
            mu = 0.0;
        }
        scale = 1.0 - mu;
        return (scale > 0.0) ? scale : INTABS_MIN_Y;
    }
    if (comp_id == INTABS_COMP_BLR) {
        if (distance_from_center <= R_seed) {
            return 1.0;
        }

        denom = sqrt(distance_from_center * distance_from_center + R_seed * R_seed);
        if (denom > 0.0) {
            mu = distance_from_center / denom;
        } else {
            mu = 0.0;
        }
       
        scale = (1.0 - mu);
        return (scale > 0.0) ? scale : INTABS_MIN_Y;
    }
    if (comp_id == INTABS_COMP_CORONA){
        denom = sqrt(distance_from_center * distance_from_center + R_seed * R_seed);
        if (denom > 0.0) {
            mu = distance_from_center / denom;
        } else {
            mu = 0.0;
        }
        scale = (1.0 - mu) * pi;
        return (scale > 0.0) ? scale : INTABS_MIN_Y;
    }
    return 1.0;
}

static double compute_mu_min_for_geometry(double distance_from_center, double R_seed) {
    double ratio;

    if (distance_from_center <= 0.0) {
        return -1.0;
    }

    if (distance_from_center < R_seed) {
        return -1.0;
    }

    ratio = R_seed / distance_from_center;
    if (ratio > 1.0) {
        ratio = 1.0;
    }
    if (ratio < 0.0) {
        ratio = 0.0;
    }

    return sqrt(1.0 - ratio * ratio);
}

static void resample_soft_field_to_common_grid(const double *nu_src,
                                               const double *n_src,
                                               unsigned int src_size,
                                               const double *nu_dst,
                                               unsigned int dst_size,
                                               double *n_dst) {
    unsigned int i;
    unsigned int left;
    unsigned int best_idx;
    double best_diff;
    double diff;
    double log_x;
    double log_x0;
    double log_x1;
    double log_y0;
    double log_y1;
    double t;

    if ((n_dst == NULL) || (nu_dst == NULL) || (dst_size == 0U)) {
        return;
    }

    for (i = 0; i < dst_size; ++i) {
        n_dst[i] = INTABS_MIN_Y;
    }

    if ((nu_src == NULL) || (n_src == NULL) || (src_size == 0U)) {
        return;
    }

    if (src_size == 1U) {
        best_idx = 0U;
        best_diff = fabs(nu_dst[0] - nu_src[0]);
        for (i = 1U; i < dst_size; ++i) {
            diff = fabs(nu_dst[i] - nu_src[0]);
            if (diff < best_diff) {
                best_diff = diff;
                best_idx = i;
            }
        }
        if (n_src[0] > 0.0) {
            n_dst[best_idx] = n_src[0];
        }
        return;
    }

    left = 0U;
    for (i = 0U; i < dst_size; ++i) {
        if ((nu_dst[i] < nu_src[0]) || (nu_dst[i] > nu_src[src_size - 1U])) {
            n_dst[i] = INTABS_MIN_Y;
            continue;
        }

        while (((left + 1U) < src_size) && (nu_dst[i] > nu_src[left + 1U])) {
            left += 1U;
        }

        if (left >= (src_size - 1U)) {
            n_dst[i] = INTABS_MIN_Y;
            continue;
        }

        if ((nu_src[left] <= 0.0) || (nu_src[left + 1U] <= nu_src[left]) || (n_src[left] <= 0.0) || (n_src[left + 1U] <= 0.0)) {
            n_dst[i] = INTABS_MIN_Y;
            continue;
        }

        log_x = log10(nu_dst[i]);
        log_x0 = log10(nu_src[left]);
        log_x1 = log10(nu_src[left + 1U]);
        log_y0 = log10(n_src[left]);
        log_y1 = log10(n_src[left + 1U]);

        if (log_x1 == log_x0) {
            n_dst[i] = n_src[left];
        } else {
            t = (log_x - log_x0) / (log_x1 - log_x0);
            n_dst[i] = pow(10.0, log_y0 + t * (log_y1 - log_y0));
        }

        if ((!isfinite(n_dst[i])) || (n_dst[i] <= 0.0)) {
            n_dst[i] = INTABS_MIN_Y;
        }
    }
}

/*
 * Convert one soft-frequency grid to dimensionless energies and bin widths.
 * For a single-bin grid, keep a finite width to avoid zeroing the integral.
 */
static void build_soft_eps_dnu(const double *nu_soft,
                               unsigned int N_soft,
                               double nu_to_eps,
                               double *eps_soft,
                               double *dnu_soft) {
    unsigned int ID_SOFT;

    if ((nu_soft == NULL) || (eps_soft == NULL) || (dnu_soft == NULL) || (N_soft == 0U)) {
        return;
    }

    for (ID_SOFT = 0U; ID_SOFT < N_soft; ++ID_SOFT) {
        eps_soft[ID_SOFT] = nu_soft[ID_SOFT] * nu_to_eps;
        if ((N_soft == 1U) && (ID_SOFT == 0U)) {
            dnu_soft[ID_SOFT] = fmax(nu_soft[ID_SOFT], INTABS_MIN_Y);
        } else if (ID_SOFT == 0U) {
            dnu_soft[ID_SOFT] = 0.0;
        } else {
            dnu_soft[ID_SOFT] = nu_soft[ID_SOFT] - nu_soft[ID_SOFT - 1U];
        }
    }
}

/*
 * Shared 3D trapezoidal integration:
 * nu_soft -> mu -> R_H for each hard-photon energy.
 *
 * `eps_soft`/`dnu_soft` can be either global (size N_soft) or per-R_H
 * (size N_R_H * N_soft), selected via `soft_grid_per_rh`.
 *
 * Soft densities can be provided either as:
 * - `n_soft_rh`: per-R_H arrays (size N_R_H * N_soft), or
 * - `n_soft_mu`: per-(R_H,mu) arrays (size N_R_H * N_theta * N_soft).
 */
static void integrate_tau_grid(const struct internal_abs_store *ia_store,
                               int use_fast_sigma,
                               const double *nu_tau,
                               unsigned int tau_size,
                               const double *eps_soft,
                               const double *dnu_soft,
                               int soft_grid_per_rh,
                               const double *n_soft_rh,
                               const double *n_soft_mu,
                               const double *one_minus_mu,
                               const double *dmu_grid,
                               const double *d_rh_grid,
                               unsigned int N_soft,
                               unsigned int N_theta,
                               unsigned int N_R_H,
                               double nu_to_eps,
                               double *tau_out) {
    unsigned int ID_GAMMA;
    unsigned int ID_RH;
    unsigned int ID_THETA;
    unsigned int ID_SOFT;
    size_t rh_mu_base;
    size_t rh_soft_base;
    size_t idx_mu;
    size_t idx_soft_energy;
    size_t idx_soft_n;
    double eps_gamma;
    double one_minus_mu_val;
    double soft_scale;
    double s_value;
    double sigma_val;
    double integrand;
    double prev_integrand;
    double nu_integral;
    double prev_mu_integral;
    double mu_integral;
    double prev_rh_integral;
    double tau_gamma;
    double n_soft_val;

    if ((nu_tau == NULL) || (tau_out == NULL) || (eps_soft == NULL) || (dnu_soft == NULL) ||
        (one_minus_mu == NULL) || (dmu_grid == NULL) || (d_rh_grid == NULL) ||
        (N_soft == 0U) || (N_theta == 0U) || (N_R_H == 0U) || (tau_size == 0U)) {
        return;
    }
    if ((n_soft_rh == NULL) && (n_soft_mu == NULL)) {
        return;
    }

    for (ID_GAMMA = 0U; ID_GAMMA < tau_size; ++ID_GAMMA) {
        eps_gamma = nu_tau[ID_GAMMA] * nu_to_eps;
        tau_gamma = 0.0;
        prev_rh_integral = 0.0;

        for (ID_RH = 0U; ID_RH < N_R_H; ++ID_RH) {
            rh_mu_base = ((size_t)ID_RH) * ((size_t)N_theta);
            rh_soft_base = ((size_t)ID_RH) * ((size_t)N_soft);
            mu_integral = 0.0;
            prev_mu_integral = 0.0;

            for (ID_THETA = 0U; ID_THETA < N_theta; ++ID_THETA) {
                idx_mu = rh_mu_base + (size_t)ID_THETA;
                one_minus_mu_val = one_minus_mu[idx_mu];
                soft_scale = eps_gamma * one_minus_mu_val * 0.5;

                nu_integral = 0.0;
                prev_integrand = 0.0;
                for (ID_SOFT = 0U; ID_SOFT < N_soft; ++ID_SOFT) {
                    if (soft_grid_per_rh != 0) {
                        idx_soft_energy = rh_soft_base + (size_t)ID_SOFT;
                    } else {
                        idx_soft_energy = (size_t)ID_SOFT;
                    }

                    s_value = soft_scale * eps_soft[idx_soft_energy];
                    integrand = 0.0;
                    if (s_value >= 1.0) {
                        if (use_fast_sigma != 0) {
                            sigma_val = sigma_gamma_gamma_fast(ia_store, s_value);
                        } else {
                            sigma_val = sigma_gamma_gamma(s_value);
                        }

                        if (n_soft_mu != NULL) {
                            idx_soft_n = (rh_mu_base + (size_t)ID_THETA) * ((size_t)N_soft) + (size_t)ID_SOFT;
                            n_soft_val = n_soft_mu[idx_soft_n];
                        } else {
                            idx_soft_n = rh_soft_base + (size_t)ID_SOFT;
                            n_soft_val = n_soft_rh[idx_soft_n];
                        }
                        integrand = sigma_val * n_soft_val * one_minus_mu_val;
                    }

                    if ((N_soft == 1U) && (ID_SOFT == 0U)) {
                        nu_integral += integrand * dnu_soft[idx_soft_energy];
                    } else if (ID_SOFT > 0U) {
                        nu_integral += 0.5 * (prev_integrand + integrand) * dnu_soft[idx_soft_energy];
                    }
                    prev_integrand = integrand;
                }

                if (ID_THETA > 0U) {
                    mu_integral += 0.5 * (prev_mu_integral + nu_integral) * dmu_grid[idx_mu];
                }
                prev_mu_integral = nu_integral;
            }

            if (ID_RH > 0U) {
                tau_gamma += 0.5 * (prev_rh_integral + mu_integral) * d_rh_grid[ID_RH];
            }
            prev_rh_integral = mu_integral;
        }

        tau_out[ID_GAMMA] = 2.0 * pi * tau_gamma;
    }
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
 * Build one combined internal-absorption tau table by summing all enabled
 * seed fields first, then running a single 3D integration.
 *
 * Integration flow:
 * 1) discover enabled components and gather one reference soft spectrum per
 *    component at a geometry-aware reference position;
 * 2) build the common hard-photon grid (nu_tau) and common soft grid;
 * 3) for each R_H and mu sample, combine component soft fields into one
 *    effective soft density grid;
 * 4) integrate with the shared nu_soft -> mu -> R_H kernel and store tau.
 */
static int eval_internal_abs_tau_total(struct blob *pt_cloned,
                                       double nu_min,
                                       unsigned int N_soft,
                                       unsigned int N_hard,
                                       unsigned int N_R_H,
                                       unsigned int N_theta,
                                       double nu_src_max) {
    int status;
    unsigned int n_enabled;
    unsigned int c;
    unsigned int i;
    unsigned int ID_RH;
    unsigned int ID_THETA;
    unsigned int ID_SOFT;
    intabs_comp_t comp_ids[3];
    struct internal_abs_component *comp_cfg[3];
    struct internal_abs_component *comp_tot;
    double comp_R_seed[3];
    double comp_R_H_ref[3];
    double comp_corona_side[3];
    double comp_mu_min;
    int comp_use_extrapolation[3];
    int comp_peak_mode[3];
    unsigned int comp_N_soft_eff[3];
    double *nu_soft_ref[3];
    double *n_soft_ref[3];
    double *nu_soft_tmp[3];
    double *n_soft_tmp[3];
    double R_H_saved_input;
    double R_H_saved;
    double distance_ref;
    double distance_from_center;
    double nu_min_eff;
    double nu_src_max_eff;
    double nu_soft_min_common;
    double nu_soft_max_common;
    double nu_soft_max_ref;
    double log_nu_min;
    double log_nu_max;
    int use_fast_sigma;
    int use_extrapolation_all;
    double mu_min_global;
    double mu_max;
    double sigma_val;
    double nu_to_eps;
    unsigned int tau_size;
    size_t idx_mu;
    size_t idx_sum_base;
    size_t idx_comp_soft;
    size_t idx_comp_mu;
    const struct internal_abs_store *ia_store;
    double *nu_soft_common;
    double *eps_soft_common;
    double *dnu_soft_common;
    double *x_grid;
    double *dx_grid;
    double *mu_grid;
    double *dmu_grid;
    double *one_minus_mu_grid;
    double *sum_soft_mu;
    double *comp_soft_interp;
    double *comp_mu_min_grid;
    double *tmp_resampled;
    struct internal_abs_geometry_ctx geom_ctx;

    /* Step 0: initialize pointers/counters so every early return is safe. */
    status = -1;
    n_enabled = 0U;
    ia_store = NULL;
    nu_soft_common = NULL;
    eps_soft_common = NULL;
    dnu_soft_common = NULL;
    x_grid = NULL;
    dx_grid = NULL;
    mu_grid = NULL;
    dmu_grid = NULL;
    one_minus_mu_grid = NULL;
    sum_soft_mu = NULL;
    comp_soft_interp = NULL;
    comp_mu_min_grid = NULL;
    tmp_resampled = NULL;

    for (c = 0U; c < 3U; ++c) {
        nu_soft_ref[c] = NULL;
        n_soft_ref[c] = NULL;
        nu_soft_tmp[c] = NULL;
        n_soft_tmp[c] = NULL;
    }

    /* Step 1: validate basic inputs and clear the cached Total slot. */
    if ((pt_cloned == NULL) || (N_soft == 0U) || (N_hard == 0U) || (N_R_H == 0U) || (N_theta == 0U)) {
        return -1;
    }
    reset_intabs_seed_build_guard();

    comp_tot = &(pt_cloned->core.internal_abs.Total);
    comp_tot->is_enabled = 0;
    comp_tot->is_valid = 0;

    /* Step 2: collect the enabled seed components to be combined. */
    if (pt_cloned->core.internal_abs.BLR.is_enabled) {
        comp_ids[n_enabled] = INTABS_COMP_BLR;
        comp_cfg[n_enabled] = &(pt_cloned->core.internal_abs.BLR);
        n_enabled += 1U;
    }
    if (pt_cloned->core.internal_abs.DT.is_enabled) {
        comp_ids[n_enabled] = INTABS_COMP_DT;
        comp_cfg[n_enabled] = &(pt_cloned->core.internal_abs.DT);
        n_enabled += 1U;
    }
    if (pt_cloned->core.internal_abs.Corona.is_enabled) {
        comp_ids[n_enabled] = INTABS_COMP_CORONA;
        comp_cfg[n_enabled] = &(pt_cloned->core.internal_abs.Corona);
        n_enabled += 1U;
    }

    if (n_enabled == 0U) {
        return 0;
    }

    ia_store = &(pt_cloned->core.internal_abs);

    R_H_saved_input = pt_cloned->core.R_H;
    R_H_saved = R_H_saved_input;
    if (R_H_saved <= 0.0) {
        R_H_saved = 1.0;
    }

    /* Aggregated flags: downgraded if any enabled component requests it. */
    use_fast_sigma = 1;
    use_extrapolation_all = 1;
    nu_soft_min_common = 0.0;
    nu_soft_max_common = 0.0;
    nu_soft_max_ref = 0.0;

    /* Step 3: prepare each component reference spectrum and geometry context. */
    for (c = 0U; c < n_enabled; ++c) {
        /* Keep each component's own sampling/extrapolation mode. */
        comp_peak_mode[c] = (comp_cfg[c]->peak_mode != 0) ? 1 : 0;
        comp_use_extrapolation[c] = (comp_cfg[c]->use_R_H_profile_extrapolation != 0) ? 1 : 0;
        comp_N_soft_eff[c] = (comp_peak_mode[c] != 0) ? INTABS_PEAK_SAMPLES : N_soft;
        if (comp_N_soft_eff[c] == 0U) {
            comp_N_soft_eff[c] = 1U;
        }

        if (comp_cfg[c]->use_sigma_gamma_gamma_fast == 0) {
            use_fast_sigma = 0;
        }
        if (comp_use_extrapolation[c] == 0) {
            use_extrapolation_all = 0;
        }

        resolve_component_geometry_ctx(pt_cloned, comp_ids[c], R_H_saved_input, R_H_saved, &geom_ctx);
        comp_R_seed[c] = geom_ctx.R_seed;
        comp_R_H_ref[c] = geom_ctx.R_H_ref;
        comp_corona_side[c] = geom_ctx.corona_side;

        nu_soft_ref[c] = (double *)calloc((size_t)comp_N_soft_eff[c], sizeof(double));
        n_soft_ref[c] = (double *)calloc((size_t)comp_N_soft_eff[c], sizeof(double));
        nu_soft_tmp[c] = (double *)calloc((size_t)comp_N_soft_eff[c], sizeof(double));
        n_soft_tmp[c] = (double *)calloc((size_t)comp_N_soft_eff[c], sizeof(double));
        if ((nu_soft_ref[c] == NULL) || (n_soft_ref[c] == NULL) || (nu_soft_tmp[c] == NULL) || (n_soft_tmp[c] == NULL)) {
            return finalize_internal_abs_total_eval(pt_cloned, comp_tot, R_H_saved_input, -1,
                                                    nu_soft_ref, n_soft_ref, nu_soft_tmp, n_soft_tmp,
                                                    nu_soft_common, eps_soft_common, dnu_soft_common,
                                                    x_grid, dx_grid, mu_grid, dmu_grid, one_minus_mu_grid,
                                                    sum_soft_mu, comp_soft_interp, comp_mu_min_grid, tmp_resampled);
        }

        distance_ref = get_seed_reference_distance(comp_ids[c], &geom_ctx);

        set_component_sampling_position(pt_cloned, comp_ids[c], distance_ref, comp_corona_side[c]);
        /*
         * Reference spectrum at a representative distance:
         * - defines common soft-grid bounds across components,
         * - provides the base shape for optional R_H extrapolation.
         */
        if (sample_seed_field(pt_cloned, comp_ids[c], comp_N_soft_eff[c], comp_peak_mode[c], nu_soft_ref[c], n_soft_ref[c]) < 0) {
            return finalize_internal_abs_total_eval(pt_cloned, comp_tot, R_H_saved_input, -1,
                                                    nu_soft_ref, n_soft_ref, nu_soft_tmp, n_soft_tmp,
                                                    nu_soft_common, eps_soft_common, dnu_soft_common,
                                                    x_grid, dx_grid, mu_grid, dmu_grid, one_minus_mu_grid,
                                                    sum_soft_mu, comp_soft_interp, comp_mu_min_grid, tmp_resampled);
        }

        for (i = 0U; i < comp_N_soft_eff[c]; ++i) {
            if (nu_soft_ref[c][i] > 0.0) {
                if ((nu_soft_min_common <= 0.0) || (nu_soft_ref[c][i] < nu_soft_min_common)) {
                    nu_soft_min_common = nu_soft_ref[c][i];
                }
                if (nu_soft_ref[c][i] > nu_soft_max_common) {
                    nu_soft_max_common = nu_soft_ref[c][i];
                }
                if (nu_soft_ref[c][i] > nu_soft_max_ref) {
                    nu_soft_max_ref = nu_soft_ref[c][i];
                }
            }
        }
    }

    if ((nu_soft_min_common <= 0.0) || (nu_soft_max_common <= 0.0)) {
        return finalize_internal_abs_total_eval(pt_cloned, comp_tot, R_H_saved_input, -1,
                                                nu_soft_ref, n_soft_ref, nu_soft_tmp, n_soft_tmp,
                                                nu_soft_common, eps_soft_common, dnu_soft_common,
                                                x_grid, dx_grid, mu_grid, dmu_grid, one_minus_mu_grid,
                                                sum_soft_mu, comp_soft_interp, comp_mu_min_grid, tmp_resampled);
    }

    /* Step 4: define the hard-photon integration range and output grid size. */
    /*
     * If nu_min is not provided, infer it from the highest reference soft
     * frequency so the gamma-gamma threshold region is reachable.
     */
    if (nu_min > 0.0) {
        nu_min_eff = nu_min;
    } else if (nu_soft_max_ref > 0.0) {
        nu_min_eff = 1.0e40 / nu_soft_max_ref;
    } else {
        nu_min_eff = 1.0e20;
    }
    if (nu_min_eff <= 0.0) {
        nu_min_eff = 1.0e20;
    }

    nu_src_max_eff = nu_src_max;
    if (nu_src_max_eff <= 0.0) {
        nu_src_max_eff = nu_min_eff;
    }

    if (nu_src_max_eff < nu_min_eff) {
        tau_size = 1U;
    } else {
        tau_size = N_hard;
    }

    if (ensure_tau_arrays(comp_tot, tau_size) < 0) {
        return finalize_internal_abs_total_eval(pt_cloned, comp_tot, R_H_saved_input, -1,
                                                nu_soft_ref, n_soft_ref, nu_soft_tmp, n_soft_tmp,
                                                nu_soft_common, eps_soft_common, dnu_soft_common,
                                                x_grid, dx_grid, mu_grid, dmu_grid, one_minus_mu_grid,
                                                sum_soft_mu, comp_soft_interp, comp_mu_min_grid, tmp_resampled);
    }

    /* Hard-photon grid where total tau is stored (log-spaced if size > 1). */
    if (tau_size == 1U) {
        comp_tot->nu_tau[0] = nu_min_eff;
    } else {
        for (i = 0U; i < tau_size; ++i) {
            comp_tot->nu_tau[i] = pow(10.0,
                                      log10(nu_min_eff) +
                                          (log10(nu_src_max_eff) - log10(nu_min_eff)) * (double)i / (double)(tau_size - 1U));
        }
    }

    /* Step 5: allocate shared grids used by the combined integration path. */
    nu_soft_common = (double *)calloc((size_t)N_soft, sizeof(double));
    eps_soft_common = (double *)calloc((size_t)N_soft, sizeof(double));
    dnu_soft_common = (double *)calloc((size_t)N_soft, sizeof(double));
    x_grid = (double *)calloc((size_t)N_R_H, sizeof(double));
    dx_grid = (double *)calloc((size_t)N_R_H, sizeof(double));
    mu_grid = (double *)calloc((size_t)N_R_H * (size_t)N_theta, sizeof(double));
    dmu_grid = (double *)calloc((size_t)N_R_H * (size_t)N_theta, sizeof(double));
    one_minus_mu_grid = (double *)calloc((size_t)N_R_H * (size_t)N_theta, sizeof(double));
    sum_soft_mu = (double *)calloc((size_t)N_R_H * (size_t)N_theta * (size_t)N_soft, sizeof(double));
    comp_soft_interp = (double *)calloc((size_t)n_enabled * (size_t)N_R_H * (size_t)N_soft, sizeof(double));
    comp_mu_min_grid = (double *)calloc((size_t)n_enabled * (size_t)N_R_H, sizeof(double));
    tmp_resampled = (double *)calloc((size_t)N_soft, sizeof(double));

    if ((nu_soft_common == NULL) || (eps_soft_common == NULL) || (dnu_soft_common == NULL) || (x_grid == NULL) || (dx_grid == NULL) ||
        (mu_grid == NULL) || (dmu_grid == NULL) || (one_minus_mu_grid == NULL) || (sum_soft_mu == NULL) ||
        (comp_soft_interp == NULL) || (comp_mu_min_grid == NULL) || (tmp_resampled == NULL)) {
        return finalize_internal_abs_total_eval(pt_cloned, comp_tot, R_H_saved_input, -1,
                                                nu_soft_ref, n_soft_ref, nu_soft_tmp, n_soft_tmp,
                                                nu_soft_common, eps_soft_common, dnu_soft_common,
                                                x_grid, dx_grid, mu_grid, dmu_grid, one_minus_mu_grid,
                                                sum_soft_mu, comp_soft_interp, comp_mu_min_grid, tmp_resampled);
    }

    /* Step 6: build the common soft-frequency grid and path-length grid. */
    if (N_soft == 1U) {
        nu_soft_common[0] = nu_soft_min_common;
    } else if (nu_soft_max_common > nu_soft_min_common) {
        log_nu_min = log10(nu_soft_min_common);
        log_nu_max = log10(nu_soft_max_common);
        for (ID_SOFT = 0U; ID_SOFT < N_soft; ++ID_SOFT) {
            nu_soft_common[ID_SOFT] = pow(10.0,
                                          log_nu_min + (log_nu_max - log_nu_min) * (double)ID_SOFT / (double)(N_soft - 1U));
        }
    } else {
        for (ID_SOFT = 0U; ID_SOFT < N_soft; ++ID_SOFT) {
            nu_soft_common[ID_SOFT] = nu_soft_min_common;
        }
    }

    nu_to_eps = HPLANCK / MEC2;
    build_soft_eps_dnu(nu_soft_common, N_soft, nu_to_eps, eps_soft_common, dnu_soft_common);

    /*
     * Path grid in normalized distance units x in [1, 1e3].
     * Each component maps x -> physical distance via its own comp_R_H_ref[c].
     */
    if (N_R_H == 1U) {
        x_grid[0] = 1.0;
    } else {
        for (ID_RH = 0U; ID_RH < N_R_H; ++ID_RH) {
            x_grid[ID_RH] = pow(10.0, 3.0 * (double)ID_RH / (double)(N_R_H - 1U));
        }
    }
    dx_grid[0] = 0.0;
    for (ID_RH = 1U; ID_RH < N_R_H; ++ID_RH) {
        dx_grid[ID_RH] = x_grid[ID_RH] - x_grid[ID_RH - 1U];
    }

    /*
     * Step 7: for each (R_H, mu), combine component soft fields into a single
     * effective soft density table `sum_soft_mu`.
     */
    for (ID_RH = 0U; ID_RH < N_R_H; ++ID_RH) {
        for (c = 0U; c < n_enabled; ++c) {
            distance_from_center = x_grid[ID_RH] * comp_R_H_ref[c];
            set_component_sampling_position(pt_cloned, comp_ids[c], distance_from_center, comp_corona_side[c]);

            if (comp_use_extrapolation[c] != 0) {
                /* Reuse reference shape, only applying geometric scaling. */
                sigma_val = compute_seed_scale_extrapolated(pt_cloned, comp_ids[c], distance_from_center, comp_R_seed[c]);
                for (ID_SOFT = 0U; ID_SOFT < comp_N_soft_eff[c]; ++ID_SOFT) {
                    nu_soft_tmp[c][ID_SOFT] = nu_soft_ref[c][ID_SOFT];
                    n_soft_tmp[c][ID_SOFT] = n_soft_ref[c][ID_SOFT] * sigma_val;
                }
            } else {
                /* Full local recomputation of the soft field at this R_H. */
                if (sample_seed_field(pt_cloned, comp_ids[c], comp_N_soft_eff[c], comp_peak_mode[c], nu_soft_tmp[c], n_soft_tmp[c]) < 0) {
                    return finalize_internal_abs_total_eval(pt_cloned, comp_tot, R_H_saved_input, -1,
                                                            nu_soft_ref, n_soft_ref, nu_soft_tmp, n_soft_tmp,
                                                            nu_soft_common, eps_soft_common, dnu_soft_common,
                                                            x_grid, dx_grid, mu_grid, dmu_grid, one_minus_mu_grid,
                                                            sum_soft_mu, comp_soft_interp, comp_mu_min_grid, tmp_resampled);
                }
            }

            /*
             * Resample each component on the shared nu grid so components can
             * be added point-by-point before the final integration.
             */
            resample_soft_field_to_common_grid(nu_soft_tmp[c], n_soft_tmp[c], comp_N_soft_eff[c], nu_soft_common, N_soft, tmp_resampled);

            idx_comp_soft = (((size_t)c) * ((size_t)N_R_H) + (size_t)ID_RH) * ((size_t)N_soft);
            for (ID_SOFT = 0U; ID_SOFT < N_soft; ++ID_SOFT) {
                comp_soft_interp[idx_comp_soft + (size_t)ID_SOFT] = tmp_resampled[ID_SOFT];
            }

            idx_comp_mu = ((size_t)c) * ((size_t)N_R_H) + (size_t)ID_RH;
            comp_mu_min_grid[idx_comp_mu] = compute_mu_min_for_geometry(distance_from_center, comp_R_seed[c]);
        }

        /*
         * Global mu-grid lower bound: most permissive acceptance among enabled
         * components at this R_H (smaller mu_min means wider angular support).
         */
        mu_min_global = 1.0;
        for (c = 0U; c < n_enabled; ++c) {
            idx_comp_mu = ((size_t)c) * ((size_t)N_R_H) + (size_t)ID_RH;
            comp_mu_min = comp_mu_min_grid[idx_comp_mu];
            if (comp_mu_min < mu_min_global) {
                mu_min_global = comp_mu_min;
            }
        }
        if (mu_min_global < -1.0) {
            mu_min_global = -1.0;
        }
        mu_max = 1.0;

        for (ID_THETA = 0U; ID_THETA < N_theta; ++ID_THETA) {
            idx_mu = ((size_t)ID_RH) * ((size_t)N_theta) + (size_t)ID_THETA;
            if (N_theta == 1U) {
                mu_grid[idx_mu] = mu_min_global;
                dmu_grid[idx_mu] = 0.0;
            } else {
                mu_grid[idx_mu] = mu_min_global + (mu_max - mu_min_global) * (double)ID_THETA / (double)(N_theta - 1U);
                if (ID_THETA == 0U) {
                    dmu_grid[idx_mu] = 0.0;
                } else {
                    dmu_grid[idx_mu] = mu_grid[idx_mu] - mu_grid[idx_mu - 1U];
                }
            }
            one_minus_mu_grid[idx_mu] = fmax(1.0 - mu_grid[idx_mu], INTABS_MIN_ONE_MINUS_MU);
        }

        for (ID_THETA = 0U; ID_THETA < N_theta; ++ID_THETA) {
            idx_mu = ((size_t)ID_RH) * ((size_t)N_theta) + (size_t)ID_THETA;
            idx_sum_base = (((size_t)ID_RH) * ((size_t)N_theta) + (size_t)ID_THETA) * ((size_t)N_soft);

            for (ID_SOFT = 0U; ID_SOFT < N_soft; ++ID_SOFT) {
                sum_soft_mu[idx_sum_base + (size_t)ID_SOFT] = 0.0;
            }

            /*
             * Build the combined soft density at fixed (R_H, mu):
             * only components visible at this mu contribute, then each
             * contribution is weighted by its geometry-dependent R_H scale.
             */
            for (c = 0U; c < n_enabled; ++c) {
                idx_comp_mu = ((size_t)c) * ((size_t)N_R_H) + (size_t)ID_RH;
                if (mu_grid[idx_mu] < comp_mu_min_grid[idx_comp_mu]) {
                    continue;
                }

                idx_comp_soft = (((size_t)c) * ((size_t)N_R_H) + (size_t)ID_RH) * ((size_t)N_soft);
                for (ID_SOFT = 0U; ID_SOFT < N_soft; ++ID_SOFT) {
                    sum_soft_mu[idx_sum_base + (size_t)ID_SOFT] +=
                        comp_soft_interp[idx_comp_soft + (size_t)ID_SOFT] * comp_R_H_ref[c];
                }
            }
        }
    }

    /* Step 8: run the shared nu_soft -> mu -> R_H integration kernel. */
    integrate_tau_grid(ia_store,
                       use_fast_sigma,
                       comp_tot->nu_tau,
                       tau_size,
                       eps_soft_common,
                       dnu_soft_common,
                       0,
                       NULL,
                       sum_soft_mu,
                       one_minus_mu_grid,
                       dmu_grid,
                       dx_grid,
                       N_soft,
                       N_theta,
                       N_R_H,
                       nu_to_eps,
                       comp_tot->tau);

    /* Step 9: persist outputs/settings in the Total component cache. */
    comp_tot->is_enabled = 1;
    comp_tot->is_valid = 1;
    comp_tot->use_R_H_profile_extrapolation = (use_extrapolation_all != 0) ? 1 : 0;
    comp_tot->use_sigma_gamma_gamma_fast = (use_fast_sigma != 0) ? 1 : 0;
    comp_tot->peak_mode = 0;
    comp_tot->N_soft = N_soft;
    comp_tot->N_hard = N_hard;
    comp_tot->N_R_H = N_R_H;
    comp_tot->N_theta = N_theta;
    comp_tot->nu_min = nu_min_eff;
    comp_tot->nu_src_max = nu_src_max_eff;

    status = (int)tau_size;
    return finalize_internal_abs_total_eval(pt_cloned, comp_tot, R_H_saved_input, status,
                                            nu_soft_ref, n_soft_ref, nu_soft_tmp, n_soft_tmp,
                                            nu_soft_common, eps_soft_common, dnu_soft_common,
                                            x_grid, dx_grid, mu_grid, dmu_grid, one_minus_mu_grid,
                                            sum_soft_mu, comp_soft_interp, comp_mu_min_grid, tmp_resampled);
}

/*
 * Core internal-absorption integration routine for one seed component.
 *
 * Inputs control the quadrature grids:
 * - `N_soft`: soft-photon frequency samples (or one peak sample in peak mode),
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
 * 6) store tau and cache setup into `pt_cloned->core.internal_abs.<component>`.
 *
 * During evaluation `pt_cloned->core.R_H` is temporarily changed for sampling and
 * restored on every return path by `finalize_internal_abs_eval`.
 */
int eval_internal_abs_tau(struct blob *pt_cloned,
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
    double distance_blob_from_seed_field_geom_center;
    double R_H_ref;
    double nu_min_eff;
    double nu_soft_max;
    double mu_min;
    double mu_max;
    double scale;
    double nu_src_max_eff;
    double nu_to_eps;
    unsigned int N_soft_eff;
    unsigned int tau_size;
    unsigned int i;
    unsigned int ID_RH;
    unsigned int ID_THETA;
    unsigned int ID_SOFT;
    size_t rh_soft_base;
    size_t rh_mu_base;
    size_t idx_soft;
    size_t idx_mu;
    const struct internal_abs_store *ia_store;
    int use_fast_sigma;
    struct internal_abs_eval_workspace ws;
    struct internal_abs_geometry_ctx geom_ctx;

    /* Step 0: initialize return state and workspace ownership. */
    status = -1;
    init_internal_abs_eval_workspace(&ws);

    /* Step 1: validate pointers, component name, and integration dimensions. */
    if (pt_cloned == NULL) {
        return -1;
    }

    /*
     * Single-component evaluations invalidate the cached total table,
     * because the total cache is only consistent after a dedicated
     * combined recomputation.
     */
    pt_cloned->core.internal_abs.Total.is_valid = 0;
    pt_cloned->core.internal_abs.Total.is_enabled = 0;

    comp_id = parse_internal_abs_component(seed_photons_name);
    comp = get_internal_abs_component_ptr(pt_cloned, comp_id);
    if ((comp_id == INTABS_COMP_INVALID) || (comp == NULL)) {
        return -1;
    }

    /* All integration dimensions must be strictly positive. */
    if ((N_soft == 0U) || (N_hard == 0U) || (N_R_H == 0U) || (N_theta == 0U)) {
        return -1;
    }

    reset_intabs_seed_build_guard();

    ia_store = &(pt_cloned->core.internal_abs);
    use_fast_sigma = (comp->use_sigma_gamma_gamma_fast != 0) ? 1 : 0;

    /* Step 2: normalize sampling mode and resolve geometry reference scales. */
    /* In peak mode collapse the soft field to one representative peak sample. */
    N_soft_eff = (peak != 0) ? INTABS_PEAK_SAMPLES : N_soft;
    if (N_soft_eff == 0U) {
        N_soft_eff = 1U;
    }

    /* Save original position; keep a positive fallback for geometric scales. */
    R_H_saved_input = pt_cloned->core.R_H;
    R_H_saved = R_H_saved_input;
    if (R_H_saved <= 0.0) {
        R_H_saved = 1.0;
    }

    resolve_component_geometry_ctx(pt_cloned, comp_id, R_H_saved_input, R_H_saved, &geom_ctx);

    /* Step 3: allocate temporary grids used by the 3D trapezoidal integration. */
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
        return finalize_internal_abs_eval(pt_cloned, comp, R_H_saved_input, -1, &ws);
    }

    /* Step 4: build one reference seed spectrum near the source-field scale. */
    distance_blob_from_seed_field_geom_center = get_seed_reference_distance(comp_id, &geom_ctx);
    set_component_sampling_position(pt_cloned, comp_id, distance_blob_from_seed_field_geom_center, geom_ctx.corona_side);
    if (sample_seed_field(pt_cloned, comp_id, N_soft_eff, peak, ws.nu_soft_ref, ws.n_soft_ref) < 0) {
        return finalize_internal_abs_eval(pt_cloned, comp, R_H_saved_input, -1, &ws);
    }

    /* Conversion nu -> dimensionless epsilon = h nu / (m_e c^2). */
    nu_to_eps = HPLANCK / MEC2;

    /*
     * Step 5: set the hard-photon range; infer nu_min from soft photons if
     * caller did not provide it.
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
        return finalize_internal_abs_eval(pt_cloned, comp, R_H_saved_input, -1, &ws);
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

    /* Step 6: build propagation grid in R_H (up to three decades). */
    R_H_ref = geom_ctx.R_H_ref;
    if (N_R_H == 1U) {
        ws.R_H_grid[0] = R_H_ref;
    } else {
        for (ID_RH = 0; ID_RH < N_R_H; ++ID_RH) {
            ws.R_H_grid[ID_RH] = pow(10.0, 3.0 * (double)ID_RH / (double)(N_R_H - 1U)) * R_H_ref;
        }
    }

    /*
     * Step 7: for each R_H sample, build soft-field slices and mu-grid geometry.
     * Soft fields come either from direct sampling or extrapolated scaling.
     */
    for (ID_RH = 0; ID_RH < N_R_H; ++ID_RH) {
        distance_blob_from_seed_field_geom_center = ws.R_H_grid[ID_RH];
        set_component_sampling_position(pt_cloned, comp_id, distance_blob_from_seed_field_geom_center, geom_ctx.corona_side);
        rh_soft_base = ((size_t)ID_RH) * ((size_t)N_soft_eff);

        if (use_R_H_profile_extrapolation != 0) {
            scale = compute_seed_scale_extrapolated(pt_cloned,
                                                    comp_id,
                                                    distance_blob_from_seed_field_geom_center,
                                                    geom_ctx.R_seed);

            for (ID_SOFT = 0; ID_SOFT < N_soft_eff; ++ID_SOFT) {
                idx_soft = rh_soft_base + (size_t)ID_SOFT;
                ws.nu_soft[idx_soft] = ws.nu_soft_ref[ID_SOFT];
                ws.n_soft[idx_soft] = ws.n_soft_ref[ID_SOFT] * scale;
            }
        } else {
            if (sample_seed_field(pt_cloned, comp_id, N_soft_eff, peak, ws.nu_soft_tmp, ws.n_soft_tmp) < 0) {
                return finalize_internal_abs_eval(pt_cloned, comp, R_H_saved_input, -1, &ws);
            }

            for (ID_SOFT = 0; ID_SOFT < N_soft_eff; ++ID_SOFT) {
                idx_soft = rh_soft_base + (size_t)ID_SOFT;
                ws.nu_soft[idx_soft] = ws.nu_soft_tmp[ID_SOFT];
                ws.n_soft[idx_soft] = ws.n_soft_tmp[ID_SOFT];
            }
        }

        build_soft_eps_dnu(ws.nu_soft + rh_soft_base,
                           N_soft_eff,
                           nu_to_eps,
                           ws.eps_soft + rh_soft_base,
                           ws.dnu_soft + rh_soft_base);

        mu_max = 1.0;
        mu_min = compute_mu_min_for_geometry(distance_blob_from_seed_field_geom_center, geom_ctx.R_seed);

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

    /* Step 8: integrate tau with the shared nu_soft -> mu -> R_H kernel. */
    integrate_tau_grid(ia_store,
                       use_fast_sigma,
                       comp->nu_tau,
                       tau_size,
                       ws.eps_soft,
                       ws.dnu_soft,
                       1,
                       ws.n_soft,
                       NULL,
                       ws.one_minus_mu,
                       ws.dmu_grid,
                       ws.d_rh_grid,
                       N_soft_eff,
                       N_theta,
                       N_R_H,
                       nu_to_eps,
                       comp->tau);

    /* Step 9: persist computed tau and the settings used to build it. */
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
    return finalize_internal_abs_eval(pt_cloned, comp, R_H_saved_input, status, &ws);
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
