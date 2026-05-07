#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include "Blazar_SED.h"

/**
 * \file Blazar_SED.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief lettura input file
 * chimata sotto funzioni
 *
 */

//=========================================================================================

struct internal_abs_async_ctx {
    struct blob *pt_worker;
    int worker_status;
};

static const char *particle_type_to_str(particle_type_t particle) {
    switch (particle) {
        case PARTICLE_ELECTRONS:
            return "electrons";
        case PARTICLE_PROTONS:
            return "protons";
        case PARTICLE_SECONDARIES_EL:
            return "secondaries_el";
        case PARTICLE_PRIMARIES_EL:
            return "primaries_el";
        default:
            return "unknown";
    }
}

static int internal_abs_enabled_on_blob(const struct blob *pt) {
    if (pt == NULL) {
        return 0;
    }
    return (pt->core.internal_abs.BLR.is_enabled ||
            pt->core.internal_abs.DT.is_enabled ||
            pt->core.internal_abs.Corona.is_enabled);
}

static int internal_abs_eval_valid_on_blob(const struct blob *pt) {
    if (pt == NULL) {
        return 0;
    }
 
    if (pt->core.internal_abs.BLR.is_enabled && (pt->core.internal_abs.BLR.is_valid == 0)) {
        return 0;
    }
    if (pt->core.internal_abs.DT.is_enabled && (pt->core.internal_abs.DT.is_valid == 0)) {
        return 0;
    }
    if (pt->core.internal_abs.Corona.is_enabled && (pt->core.internal_abs.Corona.is_valid == 0)) {
        return 0;
    }
    return 1;
}

static void copy_internal_abs_component_config(struct internal_abs_component *dst, const struct internal_abs_component *src) {
    if ((dst == NULL) || (src == NULL)) {
        return;
    }

    dst->is_enabled = src->is_enabled;
    dst->is_valid = 0;
    dst->use_R_H_profile_extrapolation = src->use_R_H_profile_extrapolation;
    dst->use_sigma_gamma_gamma_fast = src->use_sigma_gamma_gamma_fast;
    dst->peak_mode = src->peak_mode;
    dst->N_soft = src->N_soft;
    dst->N_hard = src->N_hard;
    dst->N_R_H = src->N_R_H;
    dst->N_theta = src->N_theta;
    dst->tau_size = 0U;
    dst->nu_min = src->nu_min;
    dst->nu_src_max = 0.0;
    dst->nu_tau = NULL;
    dst->tau = NULL;
}

static int alloc_internal_abs_component_arrays(struct internal_abs_component *comp, unsigned int tau_size) {
    if (comp == NULL) {
        return -1;
    }

    if (tau_size == 0U) {
        if (comp->nu_tau != NULL) {
            free(comp->nu_tau);
            comp->nu_tau = NULL;
        }
        if (comp->tau != NULL) {
            free(comp->tau);
            comp->tau = NULL;
        }
        comp->tau_size = 0U;
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
            comp->tau_size = 0U;
            return -1;
        }
    }

    comp->tau_size = tau_size;
    return 0;
}

static int merge_internal_abs_component_result(struct internal_abs_component *dst, const struct internal_abs_component *src) {
    if ((dst == NULL) || (src == NULL)) {
        return -1;
    }

    dst->is_enabled = src->is_enabled;
    dst->use_R_H_profile_extrapolation = src->use_R_H_profile_extrapolation;
    dst->use_sigma_gamma_gamma_fast = src->use_sigma_gamma_gamma_fast;
    dst->peak_mode = src->peak_mode;
    dst->N_soft = src->N_soft;
    dst->N_hard = src->N_hard;
    dst->N_R_H = src->N_R_H;
    dst->N_theta = src->N_theta;
    dst->nu_min = src->nu_min;
    dst->nu_src_max = src->nu_src_max;

    if ((src->is_enabled == 0) || (src->is_valid == 0) || (src->tau_size == 0U) || (src->nu_tau == NULL) || (src->tau == NULL)) {
        dst->is_valid = 0;
        return 0;
    }

    if (alloc_internal_abs_component_arrays(dst, src->tau_size) < 0) {
        dst->is_valid = 0;
        return -1;
    }

    memcpy(dst->nu_tau, src->nu_tau, ((size_t)src->tau_size) * sizeof(double));
    memcpy(dst->tau, src->tau, ((size_t)src->tau_size) * sizeof(double));
    dst->tau_size = src->tau_size;
    dst->is_valid = 1;

    return 0;
}

static int merge_internal_abs_result(struct blob *dst, const struct blob *src) {
    if ((dst == NULL) || (src == NULL)) {
        return -1;
    }

    if (merge_internal_abs_component_result(&(dst->core.internal_abs.BLR), &(src->core.internal_abs.BLR)) < 0) {
        return -1;
    }
    if (merge_internal_abs_component_result(&(dst->core.internal_abs.DT), &(src->core.internal_abs.DT)) < 0) {
        return -1;
    }
    if (merge_internal_abs_component_result(&(dst->core.internal_abs.Corona), &(src->core.internal_abs.Corona)) < 0) {
        return -1;
    }
    return 0;
}

static int get_internal_abs_component_by_name(struct blob *pt, const char *seed_photons_name, struct internal_abs_component **out) {
    if ((pt == NULL) || (seed_photons_name == NULL) || (out == NULL)) {
        return -1;
    }
    if (strcmp(seed_photons_name, "BLR") == 0) {
        *out = &(pt->core.internal_abs.BLR);
        return 0;
    }
    if (strcmp(seed_photons_name, "DT") == 0) {
        *out = &(pt->core.internal_abs.DT);
        return 0;
    }
    if (strcmp(seed_photons_name, "Corona") == 0) {
        *out = &(pt->core.internal_abs.Corona);
        return 0;
    }
    return -1;
}

/*
 * Build a worker blob used by isolated IA evaluations.
 *
 * The worker starts as a shallow copy of src, then internal-absorption
 * storage is re-initialized so the solver can allocate private tau buffers.
 */
static struct blob *make_internal_abs_worker_blob(const struct blob *src) {
    struct blob *worker;

    if (src == NULL) {
        return NULL;
    }

    worker = (struct blob *)malloc(sizeof(struct blob));
    if (worker == NULL) {
        return NULL;
    }

    *worker = *src;
    reset_internal_abs_store(worker);
    copy_internal_abs_component_config(&(worker->core.internal_abs.BLR), &(src->core.internal_abs.BLR));
    copy_internal_abs_component_config(&(worker->core.internal_abs.DT), &(src->core.internal_abs.DT));
    copy_internal_abs_component_config(&(worker->core.internal_abs.Corona), &(src->core.internal_abs.Corona));

    return worker;
}

static void *run_internal_abs_async(void *arg) {
    struct internal_abs_async_ctx *ctx;

    ctx = (struct internal_abs_async_ctx *)arg;
    if ((ctx == NULL) || (ctx->pt_worker == NULL)) {
        return NULL;
    }

    recompute_internal_absorption_tau(ctx->pt_worker);
    ctx->worker_status = internal_abs_eval_valid_on_blob(ctx->pt_worker) ? 0 : -1;
    return NULL;
}

/*
 * Isolated internal-absorption evaluation.
 *
 * This wrapper runs eval_internal_abs_tau() on a worker copy of the blob and
 * merges back only the requested IA component result (tau grid + config flags)
 * into the live blob. The main goal is to avoid side effects on live seed-field
 * caches/buffers that can happen with direct live-blob evaluation.
 */
int eval_internal_abs_tau_isolated(struct blob *pt,
                                   const char *seed_photons_name,
                                   double nu_min,
                                   unsigned int N_soft,
                                   unsigned int N_hard,
                                   unsigned int N_R_H,
                                   unsigned int N_theta,
                                   int use_R_H_profile_extrapolation,
                                   int peak,
                                   double nu_src_max,
                                   double R_H_override) {
    int status;
    struct blob *worker;
    struct internal_abs_component *dst_comp;
    struct internal_abs_component *src_comp;

    if (pt == NULL) {
        return -1;
    }

    dst_comp = NULL;
    if (get_internal_abs_component_by_name(pt, seed_photons_name, &dst_comp) < 0) {
        return -1;
    }

    st_gamma(1.0);
    worker = make_internal_abs_worker_blob(pt);
    if (worker == NULL) {
        dst_comp->is_valid = 0;
        return -1;
    }

    if (isfinite(R_H_override) && (R_H_override > 0.0)) {
        worker->core.R_H = R_H_override;
    }

    status = eval_internal_abs_tau(worker,
                                   seed_photons_name,
                                   nu_min,
                                   N_soft,
                                   N_hard,
                                   N_R_H,
                                   N_theta,
                                   use_R_H_profile_extrapolation,
                                   peak,
                                   nu_src_max);
    if (status >= 0) {
        src_comp = NULL;
        if (get_internal_abs_component_by_name(worker, seed_photons_name, &src_comp) < 0) {
            dst_comp->is_valid = 0;
            status = -1;
        } else if (merge_internal_abs_component_result(dst_comp, src_comp) < 0) {
            dst_comp->is_valid = 0;
            status = -1;
        }
    } else {
        dst_comp->is_valid = 0;
    }

    free_internal_abs_store(worker);
    free(worker);

    return status;
}

void show_blob(struct blob pt ) {
    printf("verbose=%d\n", pt.core.verbose);
    printf("path=%s\n", pt.core.path);
    printf("STEM=%s\n", pt.core.STEM);
    printf("do_Sync=%d\n", pt.core.do_Sync);
    printf("do_SSC=%d\n", pt.core.do_SSC);
    printf("MODE=%s\n", pt.core.MODE);
    printf("PARTICLE=%s\n", particle_type_to_str(pt.core.PARTICLE));
    printf("nu_seed_size=%d\n", pt.core.nu_seed_size);
    printf("nu_IC_size=%d\n", pt.core.nu_IC_size);
    printf("nu_start_Sync=%e\n", pt.Sync.spec.nu_min);
    printf("nu_stop_Sync=%e\n", pt.Sync.spec.nu_max);
    printf("nu_start_SSC=%e\n", pt.SSC.spec.nu_min);
    printf("nu_stop_SSC=%e\n", pt.SSC.spec.nu_max);
    printf("nu_grid_size=%d\n", pt.core.nu_grid_size);
    printf("nu_start_grid=%e\n", pt.core.nu_start_grid);
    printf("nu_stop_grid=%e\n", pt.core.nu_stop_grid);
    printf("B=%e\n", pt.core.B);
    printf("R=%e\n", pt.core.R);
    printf("BulkFactor=%e\n", pt.core.BulkFactor);
    printf("theta=%e\n", pt.core.theta);
    printf("z_cosm=%e\n", pt.core.z_cosm);
    printf("NH_pp=%e\n", pt.PP_gamma.NH_pp);
    printf("N=%e\n", pt.emitters.N);
    printf("Norm_distr=%d\n", pt.emitters.Norm_distr);
    printf("DISTR=%s\n", pt.core.DISTR);
    printf("gmin=%e\n", pt.emitters.gmin);
    printf("gmax=%e\n", pt.emitters.gmax);
    printf("do_EC_Disk=%d\n", pt.core.do_EC_Disk);
    printf("do_EC_BLR=%d\n", pt.core.do_EC_BLR);
    printf("do_EC_DT=%d\n", pt.core.do_EC_DT);
    printf("do_EC_Corona=%d\n", pt.core.do_EC_Corona);
    printf("disk type =%s\n", pt.core.disk_type);
    printf("nu_start_EC_BLR %e\n", pt.BLR.ec.spec.nu_min);
    printf("nu_stop_EC_BLR %e\n", pt.BLR.ec.spec.nu_max);
    printf("Lum Diks %e\n", pt.Disk.L_Disk);
    printf("tau BLR %e\n", pt.BLR.tau_BLR);
    printf("R_inner_Sw %e (Rs)\n", pt.Disk.R_inner_Sw);
    printf("R_ext_Sw %e (Rs)\n", pt.Disk.R_ext_Sw);
    printf("accr eff %e \n", pt.Disk.accr_eff);
    printf("T disk max (T max for MultiBB) %e\n", pt.Disk.T_Disk);
    printf("dist disk BLR (cm))%e\n", pt.BLR.R_BLR_in);
    printf("test array=%e\n", pt.SSC.spec.nuFnu_obs[0]);
    printf("nu_start_EC_DT %e\n", pt.DT.ec.spec.nu_min);
    printf("nu_stop_EC_DT %e\n", pt.DT.ec.spec.nu_max);
    printf("T_DT (T Dusty Torus) %e\n", pt.DT.T_DT);
    printf("dist disk DT (cm))%e\n", pt.DT.R_DT);
    printf("tau DT %e\n", pt.DT.tau_DT);
    printf("nu_start_EC_Corona %e\n", pt.Corona.ec.spec.nu_min);
    printf("nu_stop_EC_Corona %e\n", pt.Corona.ec.spec.nu_max);
    printf("L_Corona %e\n", pt.Corona.L_Corona);
    printf("R_Corona %e\n", pt.Corona.R_Corona);
    printf("R_H_Corona %e\n", pt.Corona.R_H_Corona);
    printf("alpha_Corona %e\n", pt.Corona.alpha_Corona);
    printf("nu_cut_low_Corona %e\n", pt.Corona.nu_cut_low_Corona);
    printf("nu_cut_Corona %e\n", pt.Corona.nu_cut_Corona);
}

void show_temp_ev(  struct temp_ev pt_ev){
    printf("do_Sync_cooling=%d\n", pt_ev.do_Sync_cooling);
    printf("do_Compton_cooling=%d\n", pt_ev.do_Compton_cooling);
    


    printf("L_inj %e (erg/s)\n", pt_ev.L_inj );
    printf("Diff_coeff %e (1/s), t_D=1/Diff_coeff %e (s)\n", pt_ev.Diff_Coeff, pt_ev.t_D0 );
    printf("Acc_coeff %e (1/s),  t_A=1/Acc_coeff %e (s)\n",  pt_ev.Acc_Coeff, pt_ev.t_A0 );
    //printf("T_esc_Coeff %e (R/c)\n", pt_ev.T_esc_Coeff_R_by_c_acc );
    printf("Esc_index %e\n", pt_ev.Esc_Index_acc );
    printf("T_start_Acc %e (s)\n", pt_ev.TStart_Acc );
    printf("T_stop_Acc  %e (s)\n", pt_ev.TStop_Acc );
    printf("T_start_inj %e (s)\n", pt_ev.TStart_Inj );
    printf("T_stop_inj  %e (s)\n", pt_ev.TStop_Inj );
    printf("Num. out file  %d \n", pt_ev.NUM_SET) ;
    printf("T size %d \n", pt_ev.T_SIZE);
}
//=========================================================================================


//=========================================================================================
struct temp_ev MakeTempEv() {
    struct temp_ev ev_root;
    //ev_root.t_unit=; //unit time in light crossing time
    //double t_acc;

    ev_root.do_Sync_cooling = 1;
    ev_root.do_Compton_cooling = 0;
    ev_root.do_Expansion = 0;
    ev_root.do_Adiabatic_cooling = 1;
    ev_root.T_COUNTER=0;



    ev_root.L_inj=1e39;
    //double *T_esc;
    ev_root.t_D0=1.0E4;
    ev_root.t_DA0=ev_root.t_D0*0.5;
    ev_root.t_A0=1.0E3;
    ev_root.Diff_Coeff=1.0/ev_root.t_D0;
    ev_root.Acc_Coeff=1.0/ev_root.t_A0;
    ev_root.Diff_Index=2.0;
    ev_root.Acc_Index=1.0;
    ev_root.Esc_Index_acc=0;
    ev_root.Esc_Index_rad=0.0;
    ev_root.m_B=1.0;
    ev_root.B_rad=1.0;
    ev_root.B_acc=1.0;
    ev_root.B_t=1.0;
    //ev_root.m_R=1.0;
    ev_root.T_esc_Coeff_R_by_c_acc=2.0;
    ev_root.T_esc_Coeff_R_by_c_rad=2.0;
    
    ev_root.TStart_Inj=0.0;
    ev_root.TStop_Inj=0.0;
    ev_root.TStart_Acc=0.0;
    ev_root.TStop_Acc=0.0;
    ev_root.Inj_temp_slope=0.0;
    ev_root.NUM_SET=50;
    ev_root.T_SIZE=1000;
    ev_root.duration=3e4;
    ev_root.E_acc_max=1E200;
    ev_root.Delta_R_acc=1E13;
    //ev_root.R_jet=1E13;
    ev_root.v_exp_by_c=1;
    //ev_root.R_jet_exp=1E13;
    ev_root.t_jet_exp=1E5;
    ev_root.R_jet_t=1E16;
    ev_root.R_H_jet_t=1E17;
    ev_root.R_H_rad_start=1E17;
    ev_root.R_rad_start=1E16;;
    ev_root.gmin_griglia = 1.0e1;
    ev_root.gmax_griglia = 1.0e8;
    ev_root.gamma_grid_size =1E4;
    ev_root.Q_inj_jetset_gamma_grid_size=1E2;

    ev_root.Lambda_max_Turb = 1e15;
    ev_root.Lambda_choer_Turb_factor=0.1;
    ev_root.Gamma_Max_Turb_L_max=Larmor_radius_to_gamma(ev_root.Lambda_max_Turb,0.1, 1.0);
    ev_root.Gamma_Max_Turb_L_coher=Larmor_radius_to_gamma(ev_root.Lambda_max_Turb*ev_root.Lambda_choer_Turb_factor,0.1, 1.0);
    ev_root.LOG_SET = 0;
    ev_root.Q_inj=NULL;
    ev_root.gamma=NULL;
    ev_root.Q_inj_jetset=NULL;
    ev_root.gamma_inj_jetset=NULL;
    //ev_root.N_gamma=NULL;
    
    ev_root.N_rad_gamma=NULL;
    ev_root.N_acc_gamma=NULL;
    ev_root.N_time=NULL;
    ev_root.T_esc_acc=NULL;
    ev_root.T_esc_rad=NULL;
    //ev_root.T_esc_ad_rad=NULL;
    ev_root.T_inj_profile=NULL;
    ev_root.T_acc_profile=NULL;
    return ev_root;
}



struct blob MakeBlob() {

    struct blob spettro_root;

    spettro_root.core.N_THREADS = 0;
    spettro_root.core.spec_array_size=static_spec_arr_size;

    spettro_root.core.WRITE_TO_FILE=0;
    spettro_root.core.BESSEL_TABLE_DONE=0;
    spettro_root.core.verbose = 0;
    sprintf(spettro_root.core.path, "./");
    sprintf(spettro_root.core.STEM, "TEST");

    spettro_root.core.PARTICLE = PARTICLE_ELECTRONS;
    spettro_root.core.do_Sync = 1;
    spettro_root.core.Sync_kernel=1;
    spettro_root.core.do_SSC = 1;
    spettro_root.core.do_IC=1;
    spettro_root.PP_gamma.do_pp_gamma=0;
    spettro_root.Bremss_ep.do_bremss_ep=0;
    //spettro_root.PP_gamma.set_pp_racc_elec = 0;
    //spettro_root.PP_gamma.set_pp_racc_gamma = 0;
    //spettro_root.PP_gamma.set_pp_racc_nu_mu = 0;
    spettro_root.PP_gamma.pp_racc_elec = 1.0;
    spettro_root.PP_gamma.pp_racc_gamma = 1.0;
    spettro_root.PP_gamma.pp_racc_nu_mu = 1.0;
    spettro_root.PP_gamma.E_th_pp_delta_approx=0.1;
    spettro_root.PP_gamma.E_pp_x_delta_approx=0.001;
    
    spettro_root.core.IC_adaptive_e_binning =0;
    spettro_root.core.do_IC_down_scattering =0;
    spettro_root.core.bulk_compton = 0;
    sprintf(spettro_root.core.MODE, "fast");
    //GRID SIZE FOR SEED
    spettro_root.core.nu_seed_size = 200;
    //GRID SIZE FOR IC
    spettro_root.core.nu_IC_size = 100;
    spettro_root.emitters.gamma_grid_size = 1000;
    spettro_root.emitters.gamma_custom_grid_size=1000;
    spettro_root.Sync.spec.nu_min = 1e8;
    spettro_root.Sync.spec.nu_max = 1e20;
    spettro_root.SSC.spec.nu_min = 1e16;
    spettro_root.SSC.spec.nu_max = 1e27;
    //GRID SIZE FOR INTERP
    spettro_root.core.nu_grid_size = 200;
    spettro_root.core.nu_start_grid = 1e8;
    spettro_root.core.nu_stop_grid = 1e27;
    spettro_root.core.emiss_lim=1.0E-120;
    spettro_root.core.B = 0.1;
    spettro_root.Sync.sin_psi = 1.0;
    spettro_root.core.R = 1e15;
    spettro_root.core.R_escape= 1E15;
    spettro_root.core.h_sh = 1;
    spettro_root.core.R_ext_sh =  0;
    spettro_root.core.R_sh =  spettro_root.core.R;
    sprintf(spettro_root.core.GEOMETRY, "spherical");
    sprintf(spettro_root.core.BEAMING_EXPR, "delta");
    spettro_root.core.BulkFactor = 10;
    spettro_root.core.beta_Gamma=eval_beta_gamma(spettro_root.core.BulkFactor);
    spettro_root.core.theta = 3.5;
    spettro_root.core.beam_obj=10.0;
    spettro_root.core.z_cosm = 0.1;
    spettro_root.PP_gamma.NH_pp = 0.1;
    spettro_root.emitters.NH_cold_to_rel_e = 0.1;
    spettro_root.emitters.N = 10;
    spettro_root.emitters.Norm_distr = 1;
    //spettro_root.Norm_distr_L_e_Sync=-1.0;
    spettro_root.emitters.Distr_e_done = 0;
    spettro_root.emitters.do_equilibrium = 0;
    sprintf(spettro_root.core.DISTR, "lp");
    spettro_root.emitters.grid_bounded_to_gamma=1;
    spettro_root.emitters.gmin = 1.0e1;
    spettro_root.emitters.gmax = 1.0e5;
    spettro_root.emitters.gmin_secondaries=spettro_root.emitters.gmin;
    spettro_root.emitters.gmax_secondaries=spettro_root.emitters.gmax*mp_by_me;
    spettro_root.emitters.gmin_griglia = -1.0;
    spettro_root.emitters.gmax_griglia = -1.0;;
    spettro_root.emitters.gamma_cooling_eq=0;
    spettro_root.emitters.T_esc_e_primaries=1;
    spettro_root.emitters.T_esc_e_secondaries=1;

    spettro_root.core.EC_stat=0; 
    spettro_root.core.EC_stat_orig=0;
    //spettro_root.EC_factor=1.0;
    spettro_root.core.do_EC_Disk = 0;
    spettro_root.core.do_EC_BLR = 0;
    spettro_root.core.do_EC_DT = 0;
    spettro_root.core.do_EC_Corona = 0;
    spettro_root.core.do_EC_CMB=0;

    spettro_root.core.do_EC_Star=0;
    spettro_root.core.do_Disk=0;
    spettro_root.core.do_DT=0;
    spettro_root.core.do_Corona=0;

    spettro_root.core.nu_planck_min_factor=1E-4;
    spettro_root.core.nu_planck_max_factor=1E2;
    spettro_root.core.mono_planck_min_factor=0.5;
    spettro_root.core.mono_planck_max_factor=2.0;
    sprintf(spettro_root.core.disk_type, "BB");
    spettro_root.core.R_H=1E17;
    spettro_root.core.R_H_orig = 1E17;
    spettro_root.core.R_H_scale_factor=1.0;
    spettro_root.core.R_ext_emit_factor=1.0;
    //spettro_root.EC_theta_lim=5.0;
    spettro_root.Disk.M_BH = 1E9;

    spettro_root.core.theta_n_int=50;
    spettro_root.core.l_n_int=50;
    spettro_root.Disk.ec.spec.nu_min = 1e13;
    spettro_root.Disk.ec.spec.nu_max = 1e26;
    spettro_root.BLR.ec.spec.nu_min = 1e13;
    spettro_root.BLR.ec.spec.nu_max = 1e26;
    spettro_root.DT.ec.spec.nu_min = 1e13;
    spettro_root.BLR.ec.spec.nu_max = 1e26;
    spettro_root.Corona.ec.spec.nu_min = 1e13;
    spettro_root.Corona.ec.spec.nu_max = 1e30;
    spettro_root.CMB.ec.spec.nu_min = 1e13;
    spettro_root.CMB.ec.spec.nu_max = 1e30;


    spettro_root.Disk.L_Disk = 1e47;
    spettro_root.BLR.tau_BLR = 1e-1;
    spettro_root.Disk.R_inner_Sw = 3.0;
    spettro_root.Disk.R_ext_Sw = 500;
    spettro_root.Disk.T_Disk = 1e5;
    spettro_root.CMB.T_CMB_0=2.725;
    spettro_root.Disk.accr_eff = 0.08;
    spettro_root.BLR.R_BLR_in = 1e18;
    spettro_root.BLR.R_BLR_out=spettro_root.BLR.R_BLR_in*2;
    spettro_root.DT.T_DT = 100;
    spettro_root.DT.R_DT = 5.0e18;
    spettro_root.DT.tau_DT = 1e-1;
    spettro_root.Corona.L_Corona = 1e45;
    spettro_root.Corona.R_Corona = 1e15;
    spettro_root.Corona.R_H_Corona = 0.0;
    spettro_root.Corona.alpha_Corona = 1.0;
    spettro_root.Corona.nu_cut_low_Corona = 0.0;
    spettro_root.Corona.nu_cut_Corona = 1e20;
    spettro_root.Corona.f_Corona_norm = 1.0;
    spettro_root.Star.L_Star = 1e33;
    spettro_root.Star.R_H_Star = 1e14;
    spettro_root.Star.T_Star =6000.;
    spettro_root.Star.theta_Star=90;

    spettro_root.emitters.gam=NULL;
    spettro_root.emitters.Q_inj_e_second=NULL;
    spettro_root.emitters.Q_inj_e_primaries=NULL;
    spettro_root.emitters.Q_inj_e=NULL;


    spettro_root.emitters.Ne=NULL;
    spettro_root.emitters.Ne_custom=NULL;
    //spettro_root.Ne_IC=NULL;
    spettro_root.emitters.Ne_jetset=NULL;
    //spettro_root.Ne_stat=NULL;
    
    spettro_root.emitters.griglia_gamma_Ne_log=NULL;
    spettro_root.emitters.gamma_e_custom=NULL;
    spettro_root.emitters.log_of_griglia_gamma_Ne_log=NULL;
    spettro_root.emitters.griglia_gamma_jetset_Ne_log=NULL;

    spettro_root.emitters.Np=NULL;
    spettro_root.emitters.Np_jetset=NULL;
    spettro_root.emitters.Np_custom=NULL;
     
    spettro_root.emitters.griglia_gamma_Np_log=NULL;
    spettro_root.emitters.griglia_gamma_jetset_Np_log=NULL;
    spettro_root.emitters.gamma_p_custom=NULL;
    spettro_root.emitters.Integrand_over_gamma_grid=NULL;

    reset_internal_abs_store(&spettro_root);
    init_sigma_gamma_gamma_table(&spettro_root);
    
    return spettro_root;
}




//=========================================================================================
void set_seed_freq_start(struct blob *pt_base){
    //pt_base->Sync.spec.nu_min = min(1e6, pt_base->core.nu_start_grid);
    pt_base->Sync.spec.nu_min = 1E6;
    pt_base->Sync.NU_INT_STOP_Sync_SSC=0;
    
    pt_base->Sync.spec.nu_max = 1E20;
    pt_base->SSC.spec.nu_min = 1E14;
    pt_base->SSC.NU_INT_STOP_COMPTON_SSC=0;

    //pt_base->SSC.spec.nu_max = max(1e30, pt_base->core.nu_stop_grid);
    pt_base->SSC.spec.nu_max=1E30;

    pt_base->Disk.ec.spec.nu_min = 1E13;
    //pt_base->Disk.ec.spec.nu_max =  max(1e30, pt_base->core.nu_stop_grid);
    pt_base->Disk.ec.spec.nu_max =1E30;
    pt_base->Disk.spec.NU_INT_MAX=0;
    
    pt_base->BLR.ec.spec.nu_min = 1E13;
    //pt_base->BLR.ec.spec.nu_max =  max(1e30, pt_base->core.nu_stop_grid);
    pt_base->BLR.ec.spec.nu_max = 1E30;
    pt_base->BLR.spec.NU_INT_MAX=0;
    
    pt_base->DT.ec.spec.nu_min = 1E13;
    pt_base->DT.spec.NU_INT_MAX=0;
    pt_base->Corona.ec.spec.nu_min = 1E13;
    pt_base->Corona.ec.spec.nu_max = 1E30;
    pt_base->Corona.spec.NU_INT_MAX=0;
    pt_base->CMB.ec.spec.nu_min = 1E13;
    pt_base->CMB.spec.NU_INT_MAX=0;
    pt_base->Star.spec.NU_INT_MAX=0;
    //pt_base->CMB.ec.spec.nu_max =  max(1e30, pt_base->core.nu_stop_grid);
    pt_base->CMB.ec.spec.nu_max = 1E30;
}


//=========================================================================================
void InitRadiative(struct blob *pt_base,unsigned int update_EC){
    //========================================================
    // Geometry Setup
    //========================================================
    pt_base->core.R_ext_sh =  pt_base->core.R_sh*(1+pt_base->core.h_sh);
    set_R_Sync(pt_base);
    pt_base->core.Vol_region = V_region(pt_base);
    pt_base->core.Surf_region = S_sphere(pt_base);
    SetBeaming(pt_base);
    pt_base->core.beta_Gamma=eval_beta_gamma(pt_base->core.BulkFactor);
    

    
    //========================================================
    // Synchrotron Parameter Initialization
    //========================================================
    pt_base->Sync.nu_B = (q_esu * pt_base->core.B) / (2 * pi * me_g * vluce_cm);
    pt_base->Sync.UB = pow(pt_base->core.B, 2.0) / (8.0 * pi); /*dens. ener. B */
    
    if (pt_base->core.verbose>0) {
        printf("gmin %e   gmax %e \n", pt_base->emitters.gmin, pt_base->emitters.gmax);
        printf("UB=%e \n", pt_base->Sync.UB);
        printf("nu_B_non_rel=%e \n", pt_base->Sync.nu_B);
        printf("beaming factor =%e\n", pt_base->core.beam_obj);
    }
    
    //COSTANTI PER ALFA=FIXED E KERNEL DELTA O KERNEL 2
    pt_base->Sync.C1_Sync_K53 = pow(3, 0.5) * pow(q_esu, 3.0) * pt_base->Sync.sin_psi;
    pt_base->Sync.C1_Sync_K53 *= pt_base->core.B / (MEC2) * one_by_four_pi;
    pt_base->Sync.C2_Sync_K53 = 2.0/(3*pt_base->Sync.nu_B* pt_base->Sync.sin_psi);
    
    pt_base->Sync.C1_Sync_K_AVE= 4*pi*pow(3, 0.5) * pow(q_esu, 2.0)*pt_base->Sync.nu_B/(vluce_cm) * one_by_four_pi;
    pt_base->Sync.C2_Sync_K_AVE=1.0/(3*pt_base->Sync.nu_B);

    pt_base->Sync.C3_Sync_K53 = -1.0 * pow(3, 0.5) * pow(q_esu, 3.0) / (8 * pi * MEC2 * me_g);

    pt_base->Sync.COST_Sync_COOLING = SIGTH * vluce_cm/MEC2;


    //==================================
    //  Bessel Fucntion Setup
    //==================================
    //exit(1);
    if (pt_base->core.BESSEL_TABLE_DONE == 0){
        printf("Bessel Functions\n");
        tabella_Bessel(pt_base);
    }
    
    //========================================================
    // Compton Parameter Initialization
    //========================================================
    pt_base->core.COST_IC_K1 = 3.0 * SIGTH * vluce_cm / 4.0;
    pt_base->core.COST_IC_COOLING = (4.0/3.0) * SIGTH * vluce_cm*HPLANCK/MEC2;


    //========================================================
    // pp parameters initialization
    //========================================================
    //pt_base->PP_gamma.set_pp_racc_elec = 0;
    //pt_base->PP_gamma.set_pp_racc_gamma = 0;
    //pt_base->PP_gamma.set_pp_racc_nu_mu = 0;
    //pt_base->PP_gamma.pp_racc_elec = 1.0;
    //pt_base->PP_gamma.pp_racc_gamma = 1.0;
    //pt_base->PP_gamma.pp_racc_nu_mu = 1.0;

    //========================================================
    // EC  Initialization
    //========================================================
    pt_base->core.R_H_orig=pt_base->core.R_H;
    pt_base->core.EC_stat_orig = pt_base->core.EC_stat;
    if (update_EC>0){
        if (pt_base->core.do_EC_Disk == 1 || pt_base->core.do_EC_BLR == 1 || pt_base->core.do_EC_DT == 1  || pt_base->core.do_EC_Corona == 1 || pt_base->core.do_EC_Star == 1 || pt_base->core.do_EC_CMB == 1 || pt_base->core.do_Disk==1 || pt_base->core.do_DT==1 || pt_base->core.do_Corona==1 || pt_base->core.do_Star==1) 
            {
                spectra_External_Fields(1, pt_base, 1);
        }
    }
    //========================================================

}

void Init(struct blob *pt_base, double luminosity_distance) {
    // if luminosity_distance is negative is evaluated internally
    // otherwise the passed value is used

    //struct spettro *pt_base;
    //double (*pf) (struct spettro *, double);
    
    unsigned int i;
    //char * ENV;
    pt_base->core.SYSPATH=getenv("BLAZARSED");
    set_seed_freq_start(pt_base);

    //pt_base->core.emiss_lim=1.0E-120;

    //sprintf(ENV,'%s',getenv("BLAZARSED"));
    //return;
    //printf("CIAO =%s\n",pt_base->core.SYSPATH);

    //sprintf(pt_base->core.SYSPATH,'%s', ENV);
    //return;
    if (pt_base->core.verbose) {
        printf("SYSPATH =%s\n", pt_base->core.SYSPATH);
        printf("STEM=%s\n", pt_base->core.STEM);
        printf("PATH =%s\n", pt_base->core.path);
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>> Satic Case Initilization <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    }

    //======================================
    // Arrays SetUp
    //======================================
    if (pt_base->core.nu_seed_size>=pt_base->core.spec_array_size){
    	pt_base->core.nu_seed_size=pt_base->core.spec_array_size-1;
    	if (pt_base->core.verbose){
            printf("!!! Warning nu_seed_size  was gt spec_array size \n");
            printf("now set  to spec_array size %d \n ",pt_base->core.spec_array_size-1);
        }
    }
    if (pt_base->core.nu_IC_size>=pt_base->core.spec_array_size){
    	pt_base->core.nu_IC_size=pt_base->core.spec_array_size-1;
    	if (pt_base->core.verbose){
            printf("!!! Warning  nu_IC_size  was gt spec_array size \n ");
            printf("now set  to spec_array size %d\n ",pt_base->core.spec_array_size-1);
        }
    }
    for (i = 0; i < static_spec_arr_size; i++) {
        pt_base->SSC.q_comp[i] = 0.0;
        pt_base->Sync.spec.j_nu[i] = 0.0;
        pt_base->SSC.spec.j_nu[i] = 0.0;
        pt_base->Disk.ec.spec.j_nu[i] = 0.0;
        pt_base->BLR.ec.spec.j_nu[i] = 0.0;
        pt_base->DT.ec.spec.j_nu[i] = 0.0;
        pt_base->Corona.ec.spec.j_nu[i] = 0.0;
        pt_base->Star.ec.spec.j_nu[i] = 0.0;
        pt_base->CMB.ec.spec.j_nu[i] = 0.0;
        pt_base->Sync.alfa_Sync[i] = 0.0;
        pt_base->Sync.spec.I_nu[i] = 0.0;
        pt_base->PP_neutrino.spec_tot.j_nu[i]=0.0;
        pt_base->PP_neutrino.spec_mu.j_nu[i]=0.0;
        pt_base->PP_neutrino.spec_e.j_nu[i]=0.0;
        pt_base->PP_gamma.spec.j_nu[i]=0.0;
        pt_base->Bremss_ep.spec.j_nu[i]=0.0;
    }

    for (i = 0; i < static_spec_arr_size; i++){
        pt_base->Sync.spec.nuFnu_obs[i]=0.0;
        pt_base->SSC.spec.nuFnu_obs[i]=0.0;
        pt_base->Disk.ec.spec.nuFnu_obs[i]=0;
        pt_base->BLR.ec.spec.nuFnu_obs[i]=0;
        pt_base->DT.ec.spec.nuFnu_obs[i]=0;
        pt_base->Corona.ec.spec.nuFnu_obs[i]=0;
        pt_base->Star.ec.spec.nuFnu_obs[i]=0;
        pt_base->CMB.ec.spec.nuFnu_obs[i]=0;
        pt_base->Disk.spec.nuFnu_obs[i]=0;
        pt_base->DT.spec.nuFnu_obs[i]=0;
        pt_base->Corona.spec.nuFnu_obs[i]=0;
        pt_base->Star.spec.nuFnu_obs[i]=0;
        pt_base->PP_gamma.spec.nuFnu_obs[i]=0;
        pt_base->PP_neutrino.spec_tot.nuFnu_obs[i]=0;
        pt_base->PP_neutrino.spec_mu.nuFnu_obs[i]=0;
        pt_base->PP_neutrino.spec_e.nuFnu_obs[i]=0;
        pt_base->Bremss_ep.spec.nuFnu_obs[i]=0;
    }


    //set file number counter
    pt_base->core.OUT_FILE = 1;

    if (luminosity_distance<0){

        pt_base->core.dist = dist_lum_cm(pt_base->core.z_cosm);
    }
    else{
        pt_base->core.dist = luminosity_distance;
    }
    InitRadiative(pt_base,1);
    pt_base->Disk.R_Sw = eval_R_Sw(pt_base->Disk.M_BH);
    pt_base->Disk.R_ext = pt_base->Disk.R_ext_Sw * pt_base->Disk.R_Sw;

   

    if (pt_base->core.verbose) {
        printf("Distanza rigorosa=%e in Mpc \n", pt_base->core.dist/(1.0e6*1.0e2));
        printf("Distanza rigorosa=%e in cm \n", pt_base->core.dist);
    }

    pt_base->emitters.Distr_e_done = 0;
    pt_base->emitters.Distr_e_pp_done = 0;

    if (pt_base->core.verbose) {     
        printf("******************************  Geometry  *********************************\n");
        printf("Volume Geom.=%e\n", pt_base->core.Vol_region);
    }
    
    if (pt_base->core.PARTICLE == PARTICLE_ELECTRONS) {
        if (pt_base->emitters.do_equilibrium == 1) {
            InitNeEquilibrium(pt_base);
        } else {
            InitNe(pt_base);
        }
        pt_base->emitters.N_tot_e_Sferic = pt_base->core.Vol_region * pt_base->emitters.N_e;  
        if (pt_base->core.verbose) {
            FindNe_NpGp(pt_base);
            EvalU_e(pt_base);     
            printf("********************       Leptonic Scenario       ********************\n");
            if (pt_base->emitters.do_equilibrium == 1) {
                printf("equilibrium mode=ON\n");
            }
            printf("type of distr=%d\n", pt_base->emitters.TIPO_DISTR);
            printf("*******  Leptonic Energetic   **********\n");
            printf("N_e=%e Ne/Ne_0=%e\n", pt_base->emitters.N_e, pt_base->emitters.N / pt_base->emitters.N_0e);
            printf("Total number of electrons    =%e\n", pt_base->emitters.N_tot_e_Sferic);
            printf("Gamma_p of N(gamma)*gamma^2 = %e\n", pt_base->emitters.Gamma_p2);
            printf("Gamma_p of N(gamma)*gamma^3 = %e\n", pt_base->emitters.Gamma_p3);
            printf("Peak of  N(gamma)*gamma^2 = %e\n", pt_base->emitters.Np2);
            printf("Peak of  N(gamma)*gamma^3 = %e\n", pt_base->emitters.Np3);
            printf("U_e   blob rest frame =%e erg/cm^3\n", pt_base->emitters.U_e);
            printf("U_e/U_b =%e\n", pt_base->emitters.U_e / pt_base->Sync.UB);
            printf("E_tot (electron)  blob rest frame =%e erg     \n", pt_base->emitters.E_tot_e);
            printf("************************************************************************\n");
        }
    }
    //}else if (strcmp(pt_base->core.PARTICLE, "electrons-equilibrium") == 0){
    //    InitNeEquilibirum(pt_base);
    //    pt_base->emitters.N_tot_e_Sferic = pt_base->core.Vol_region * pt_base->emitters.N_e;
    //    FindNe_NpGp(pt_base);
    //    EvalU_e(pt_base);
    else if (pt_base->core.PARTICLE == PARTICLE_PROTONS) {
        Init_Np_Ne_pp(pt_base);        
        pt_base->emitters.N_tot_p_Sferic = pt_base->core.Vol_region * pt_base->emitters.N_p;             
        //EvalU_p(pt_base);             
        pt_base->emitters.N_tot_e_Sferic = pt_base->core.Vol_region * pt_base->emitters.N_e_pp;
        //EvalU_e(pt_base);
        //FindNe_NpGp(pt_base);
        if (pt_base->core.verbose) {
            EvalU_p(pt_base);
            EvalU_e(pt_base);
            FindNe_NpGp(pt_base);           
            printf("***********************       Hadronic Scenario           ********************\n");
            printf("****** Generate Np and Ne form secondaries **************\n");
          
            printf("******************      Hadronic Energetic    **********\n");
            printf("N_p=%e N_p/N_0p=%e\n", pt_base->emitters.N, pt_base->emitters.N / pt_base->emitters.N_0p);
            printf("Total number of p    =%e\n", pt_base->emitters.N_tot_p_Sferic);
            printf("U_p   blob rest frame =%e erg/cm^3\n", pt_base->emitters.U_p);
            printf("U_p/U_b =%e\n", pt_base->emitters.U_p / pt_base->Sync.UB);
            printf("E_tot (protons)  blob rest frame =%e erg     \n", pt_base->emitters.E_tot_p);
            printf("*******  Scondaries Leptonic Energetic  ****************\n");
            printf("In this case N_0e=1, leptons come from pp, N_e_pp=%e\n", pt_base->emitters.N_e_pp);
            printf("Total number of secondary electrons    =%e\n", pt_base->emitters.N_tot_e_Sferic);
            printf("U_e   blob rest frame =%e erg/cm^3\n", pt_base->emitters.U_e);
            printf("U_e/U_b =%e\n", pt_base->emitters.U_e / pt_base->Sync.UB);
            printf("E_tot (electron)  blob rest frame =%e erg     \n", pt_base->emitters.E_tot_e);
            printf("Gamma_p of N(gamma)*gamma^2 = %e\n", pt_base->emitters.Gamma_p2);
            printf("Gamma_p of N(gamma)*gamma^3 = %e\n", pt_base->emitters.Gamma_p3);
            printf("Peak of  N(gamma)*gamma^2 = %e\n", pt_base->emitters.Np2);
            printf("Peak of  N(gamma)*gamma^3 = %e\n", pt_base->emitters.Np3);
        }
    }

}
 
void Run_SED(struct blob *pt_base){
    double nuFnu_obs_ref_EC;
    int ia_enabled;
    int ia_async_started;
    int ia_thread_join_ok;
    int ia_isolated_eval_ok;
    pthread_t ia_thread;
    struct blob *ia_worker_blob;
    struct internal_abs_async_ctx ia_ctx;

    if (pt_base->core.verbose) {
        printf("STEM=%s\n", pt_base->core.STEM);
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>> RUN      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    }

    ia_enabled = internal_abs_enabled_on_blob(pt_base);
    ia_async_started = 0;
    ia_thread_join_ok = 0;
    ia_isolated_eval_ok = 0;
    ia_worker_blob = NULL;
    ia_ctx.pt_worker = NULL;
    ia_ctx.worker_status = -1;

    if (ia_enabled) {
        /*
         * st_gamma() uses a lazy static initializer in func_math.c.
         * Force one serial call before starting a parallel IA worker.
         */
        st_gamma(1.0);
        ia_worker_blob = make_internal_abs_worker_blob(pt_base);
        if (ia_worker_blob != NULL) {
            ia_ctx.pt_worker = ia_worker_blob;
            if (pthread_create(&ia_thread, NULL, run_internal_abs_async, &ia_ctx) == 0) {
                ia_async_started = 1;
            } else {
                recompute_internal_absorption_tau(ia_worker_blob);
                ia_ctx.worker_status = internal_abs_eval_valid_on_blob(ia_worker_blob) ? 0 : -1;
            }
        }
    }
    //==================================================
    // Evaluate hadronic pp Spectrum
    //==================================================
    if ((pt_base->core.PARTICLE == PARTICLE_PROTONS) && pt_base->PP_gamma.do_pp_gamma) {

        spettro_pp_gamma(1, pt_base);
    }

    if ((pt_base->core.PARTICLE == PARTICLE_PROTONS) && pt_base->PP_neutrino.do_pp_neutrino) {
        spettro_pp_neutrino(1,pt_base);
    }

    //==================================================
    // Evaluate Bremms sp Spectrum
    //==================================================
    if (pt_base->Bremss_ep.do_bremss_ep) {
        spettro_bremss_ep(1, pt_base);
    }


    //==================================================
    // Evaluate Synchrotron Spectrum
    //==================================================
    if (pt_base->core.do_Sync != 0) {
        spettro_sincrotrone(1, pt_base);
    }



    //==================================================
    // Evaluate SSC Spectrum
    //==================================================
    if (pt_base->core.do_SSC && pt_base->core.do_IC) {
       spettro_compton(1, pt_base);
    }


    //==================================================
    // Evaluate EC Spectrum
    //==================================================
	if (pt_base->core.do_IC) {
		if (pt_base->core.do_EC_Disk == 1 || pt_base->core.do_EC_BLR == 1 || pt_base->core.do_EC_DT == 1  || pt_base->core.do_EC_Corona == 1 || pt_base->core.do_EC_Star == 1 || pt_base->core.do_EC_CMB == 1 || pt_base->core.do_Disk==1 || pt_base->core.do_DT==1 || pt_base->core.do_Corona==1 || pt_base->core.do_Star==1) 
        {
                if (pt_base->core.do_EC_Star == 1) {
                    
                    pt_base->core.EC = 4;
                    spettro_EC(1, pt_base);
                }
                if (pt_base->core.do_EC_Disk == 1 || pt_base->core.do_Disk==1) {
                   
                    pt_base->core.EC = 1;
                    if (set_condition_EC_correction(pt_base, pt_base->Disk.R_inner) > 0)
                    {
                        pt_base->core.R_H = pt_base->Disk.R_inner/10;
                        Build_I_nu_Disk(pt_base);
                        spettro_EC(1, pt_base);
                        nuFnu_obs_ref_EC = get_EC_reference(pt_base, pt_base->Disk.ec.spec.nuFnu_obs);
                        pt_base->core.R_H = pt_base->core.R_H_orig;
                        Build_I_nu_Disk(pt_base);
                    }
                    spettro_EC(1, pt_base);
                    if (set_condition_EC_correction(pt_base, pt_base->Disk.R_inner) > 0){
                        update_EC_for_bp(pt_base, nuFnu_obs_ref_EC, pt_base->Disk.R_inner, pt_base->core.nu_IC_size, pt_base->Disk.ec.spec.nuFnu_obs, pt_base->Disk.ec.spec.nu_obs);
                    }
                }
                if (pt_base->core.do_EC_BLR == 1) {
                    pt_base->core.EC = 2;
                    if (set_condition_EC_correction(pt_base, pt_base->BLR.R_BLR_out) > 0)
                    {
                        pt_base->core.R_H = max(1,pt_base->BLR.R_BLR_in/1E10);
                        Build_I_nu_BLR(pt_base);
                        spettro_EC(1, pt_base);
                        nuFnu_obs_ref_EC = get_EC_reference(pt_base, pt_base->BLR.ec.spec.nuFnu_obs);
                        pt_base->core.R_H = pt_base->core.R_H_orig;
                        Build_I_nu_BLR(pt_base);
                    }
                    spettro_EC(1, pt_base);
                    if (set_condition_EC_correction(pt_base, pt_base->BLR.R_BLR_out) > 0){
                        update_EC_for_bp(pt_base, nuFnu_obs_ref_EC, pt_base->BLR.R_BLR_out, pt_base->core.nu_IC_size, pt_base->BLR.ec.spec.nuFnu_obs, pt_base->BLR.ec.spec.nu_obs);
                    }
                }
                if (pt_base->core.do_EC_DT == 1) {
                    pt_base->core.EC = 3;
                    //printf("RUN 1 R_H=%e c=%d , EC_stat=%d\n", pt_base->core.R_H, set_condition_EC_correction(pt_base, pt_base->DT.R_DT), pt_base->core.EC_stat);
                    if (set_condition_EC_correction(pt_base, pt_base->DT.R_DT) > 0)
                    {
                        pt_base->core.R_H = 0;
                        Build_I_nu_DT(pt_base);
                        spettro_EC(1, pt_base);
                        nuFnu_obs_ref_EC = get_EC_reference(pt_base, pt_base->DT.ec.spec.nuFnu_obs);
                        pt_base->core.R_H = pt_base->core.R_H_orig;
                        Build_I_nu_DT(pt_base);
                    }
                    //printf("RUN 2 R_H=%e c=%d , EC_stat=%d\n", pt_base->core.R_H, set_condition_EC_correction(pt_base, pt_base->DT.R_DT), pt_base->core.EC_stat);
                    //printf("RUN 3 R_H=%e c=%d \n", pt_base->core.R_H, set_condition_EC_correction(pt_base, pt_base->DT.R_DT));
                    spettro_EC(1, pt_base);
                    //printf("RUN 4 R_H=%e c=%d \n", pt_base->core.R_H,set_condition_EC_correction(pt_base, pt_base->DT.R_DT) );
                    if (set_condition_EC_correction(pt_base, pt_base->DT.R_DT) > 0)
                    {
                        update_EC_for_bp(pt_base, nuFnu_obs_ref_EC, pt_base->DT.R_DT, pt_base->core.nu_IC_size, pt_base->DT.ec.spec.nuFnu_obs, pt_base->DT.ec.spec.nu_obs);
                    }
                    //printf("RUN 4 R_H=%e c=%d \n", pt_base->core.R_H, set_condition_EC_correction(pt_base, pt_base->DT.R_DT));
                }
                if (pt_base->core.do_EC_CMB == 1) {
                    
                    pt_base->core.EC = 5;
                    spettro_EC(1, pt_base);
                }
                if (pt_base->core.do_EC_Corona == 1) {
                    pt_base->core.EC = 6;
                    if (set_condition_EC_correction(pt_base, pt_base->Corona.R_Corona) > 0)
                    {
                        double R_blob_Corona_ref, R_H_Corona_ref;
                        R_blob_Corona_ref = max(1, pt_base->Corona.R_Corona / 1E10);
                        if (pt_base->core.R_H >= pt_base->Corona.R_H_Corona){
                            R_H_Corona_ref = pt_base->Corona.R_H_Corona + R_blob_Corona_ref;
                        }
                        else{
                            R_H_Corona_ref = pt_base->Corona.R_H_Corona - R_blob_Corona_ref;
                        }
                        pt_base->core.R_H = max(R_H_Corona_ref, 0.0);
                        Build_I_nu_Corona(pt_base);
                        spettro_EC(1, pt_base);
                        nuFnu_obs_ref_EC = get_EC_reference(pt_base, pt_base->Corona.ec.spec.nuFnu_obs);
                        pt_base->core.R_H = pt_base->core.R_H_orig;
                        Build_I_nu_Corona(pt_base);
                    }
                    spettro_EC(1, pt_base);
                    if (set_condition_EC_correction(pt_base, pt_base->Corona.R_Corona) > 0)
                    {
                        update_EC_for_bp(pt_base, nuFnu_obs_ref_EC, pt_base->Corona.R_Corona, pt_base->core.nu_IC_size, pt_base->Corona.ec.spec.nuFnu_obs, pt_base->Corona.ec.spec.nu_obs);
                    }
                }
              
            }
        //printf("=>done\n");
    }
    //==================================================
    //Sum Up all the Spectral Components
    //==================================================
    if (ia_enabled) {
        if (ia_worker_blob != NULL) {
            if (ia_async_started) {
                if (pthread_join(ia_thread, NULL) == 0) {
                    ia_thread_join_ok = 1;
                }
            } else {
                ia_thread_join_ok = 1;
            }

            if (ia_thread_join_ok && (ia_ctx.worker_status == 0)) {
                if (merge_internal_abs_result(pt_base, ia_worker_blob) == 0) {
                    ia_isolated_eval_ok = 1;
                }
            }

            free_internal_abs_store(ia_worker_blob);
            free(ia_worker_blob);
            ia_worker_blob = NULL;
        }

        /*
         * No live-blob fallback: keep isolated-only IA behavior.
         * If isolated evaluation failed, invalidate IA components.
         */
        if (ia_isolated_eval_ok == 0) {
            if (pt_base->core.internal_abs.BLR.is_enabled) {
                pt_base->core.internal_abs.BLR.is_valid = 0;
            }
            if (pt_base->core.internal_abs.DT.is_enabled) {
                pt_base->core.internal_abs.DT.is_valid = 0;
            }
            if (pt_base->core.internal_abs.Corona.is_enabled) {
                pt_base->core.internal_abs.Corona.is_valid = 0;
            }
        }
    }
    common_grid_spectra(1, pt_base);

    //==================================================
    // Energetic
    //==================================================
    //printf("Energetic computation (output to file)\n");


    //EnergeticOutput(pt_base);
    //CoolingRates(pt_base);


}

//=========================================================================================

//==================================================
//Funtions To access Ne and Spectral components form Python
//==================================================
double get_array(double * arr, unsigned int id, unsigned int size){
	if ((id >=0) && (id <= size)){
		return arr[id];
	}
	else{
        printf("exceeded array size in get_spectral_array\n");
        exit(0);
	}
}



double get_spectral_array(double * arr, struct blob * pt, unsigned int id){
	if ((id >=0) && (id <= pt->core.nu_grid_size)){
		return arr[id];
	}
	else{
        printf("exceeded array size in get_spectral_array\n");
        exit(0);
	}
}

void set_spectral_array(double *arr, struct blob *pt, unsigned int id, double val)
{
    if ((id >= 0) && (id <= pt->core.nu_grid_size))
    {
        arr[id]=val;
    }
    else
    {
        printf("exceeded array size in get_spectral_array\n");
        exit(0);
    }
}

double get_elec_array(double * arr, struct blob *pt, unsigned int id){
	if ((id>=0) && (id<=pt->emitters.gamma_grid_size)){
		return arr[id];
	}
	else{
        printf("exceeded array size in get_elec_array\n");
        exit(0);
	}
}

double get_Q_inj_array(double *arr, struct temp_ev *pt_ev, unsigned int id)
{
    if ((id >= 0) && (id <= pt_ev->gamma_grid_size))
    {
        return arr[id];
    }
    else
    {
        printf("exceeded array size in get_Q_inj_array\n");
        exit(0);
    }
}

double get_temp_ev_N_gamma_array(double *arr, struct temp_ev *pt_ev, unsigned int row, unsigned int col)
{
    if (((col >= 0) && (col <= pt_ev->gamma_grid_size)) && ((row >= 0) && (row <= pt_ev->NUM_SET)))
    {
        return arr[row * pt_ev->gamma_grid_size + col];
    }
    else
    {
        printf("exceeded array size in get_temp_ev_N_gamma_array\n");
        exit(0);
    }
}

double get_temp_ev_N_time_array(double *arr, struct temp_ev *pt_ev, unsigned int id)
{
    if ((id >= 0) && (id <= pt_ev->NUM_SET))
    {
            return arr[id];
        }
        else
        {
            printf("exceeded array size in get_temp_ev_N_time_array\n");
            exit(0);
        }
    }

double get_temp_ev_gamma_array(double *arr, struct temp_ev *pt_ev, unsigned int id)
    {
        if ((id >= 0) && (id <= pt_ev->gamma_grid_size))
        {
            return arr[id];
        }
        else
        {
            printf("exceeded array size in get_temp_ev_gamma_array\n");
            exit(0);
        }
    }

void set_elec_array(double * arr,struct blob *pt, double val, unsigned int id){
    if ((id>=0) && (id<=pt->emitters.gamma_grid_size)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_elec_array\n");
            exit(0);
        }
}

void set_q_inj_user_array(double * arr,struct temp_ev *pt, double val, unsigned int id){
    if ((id>=0) && (id<=pt->Q_inj_jetset_gamma_grid_size)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_elec_array\n");
            exit(0);
        }
}

void set_temp_ev_Time_array(double * arr,struct temp_ev *pt, double val, unsigned int id){
    if ((id>=0) && (id<=pt->T_SIZE)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_T_inj_profile\n");
            exit(0);
        }
}


void set_elec_custom_array(double * arr, struct blob *pt,double val, unsigned int id){
    if ((id>=0) && (id<=pt->emitters.gamma_custom_grid_size)){
           arr[id]=val;
        }
        else{
            printf("exceeded array size in set_elec_custom_array\n");
            exit(0);
        }
}

void set_bessel_table(double *arr, struct blob *pt, double val, unsigned int id)
{
    if ((id >= 0) && (id <= static_bess_table_size))
    {
        arr[id] = val;
    }
    else
    {
        printf("exceeded array size in set_bessel_table\n");
        exit(0);
    }
}

double get_temp_ev_array_static(double *arr, unsigned int id){
    if ((id >= 0) && (id <= static_ev_arr_grid_size))
    {
        return arr[id];
    }
    else
    {
        printf("exceeded array size get_temp_ev_array_static\n");
        exit(0);
    }
}
//=========================================================================================
void SetBeaming(struct blob *pt){

	if (strcmp(pt->core.BEAMING_EXPR, "delta") == 0) {
	        pt->core.beam_obj = pt->core.beam_obj;
            pt->core.BulkFactor = pt->core.beam_obj;
    }

	else if (strcmp(pt->core.BEAMING_EXPR, "bulk_theta") == 0) {
		pt->core.beam_obj = get_beaming(pt->core.BulkFactor,pt->core.theta);
	}

	else {
		printf("BEAMING_EXPR variable set to wrong value, posible delta or bulk_theta \n");
		exit(0);
	}

	if (pt->core.verbose) {
	     printf("beaming set to  %e\n",pt->core.beam_obj);
	}
}
