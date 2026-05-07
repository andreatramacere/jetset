import numpy as np
from .base_class import TestBase


class TestInternalAbsorption(TestBase):

    def integration_suite(self, plot=False):
        self.test_internal_absorption_enable_eval_remove(plot=plot)
        self.test_internal_absorption_serialization(plot=plot)

    def _build_internal_abs_jet(self):
        from jetset.jet_model import Jet

        j = Jet(name='compact_int_abs', emitters_distribution='bkn', beaming_expr='bulk_theta')
        j.add_EC_component(EC_components_list=['EC_DT', 'EC_BLR', 'EC_Corona'], disk_type='BB')

        j.parameters.z_cosm.val = 0.03
        j.parameters.L_Disk.val = 2E45
        j.parameters.R_H.val = 1E18
        j.parameters.tau_DT.val = 0.1
        j.parameters.tau_BLR.val = 0.1
        j.parameters.L_Corona.val = 5E44
        j.parameters.R_Corona.val = 5E15
        j.parameters.R_H_Corona.val = 2E17
        j.parameters.alpha_Corona.val = 1.1
        j.parameters.nu_cut_Corona.val = 1E20
        j.parameters.B.val = 0.2

        j.set_gamma_grid_size(120)
        j.set_IC_nu_size(80)
        return j

    @staticmethod
    def _assert_total_cache_disabled(jet):
        total = jet._blob.core.internal_abs.Total
        assert int(total.is_enabled) == 0
        assert int(total.is_valid) == 0

    def test_internal_absorption_enable_eval_remove(self, plot=False):
        j = self._build_internal_abs_jet()
        nu = np.logspace(20, 29, 120)

        y_no_ia = np.asarray(j.eval(nu=nu, get_model=True), dtype=float)
        assert np.all(np.isfinite(y_no_ia))
        r_h_ref = float(j._blob.core.R_H)

        j.enable_internal_absorption('DT', N_soft=12, N_hard=12, N_R_H=10, N_theta=10, use_sigma_gamma_gamma_fast=True)
        j.enable_internal_absorption('BLR', N_soft=12, N_hard=12, N_R_H=10, N_theta=10)
        j.enable_internal_absorption('Corona', N_soft=12, N_hard=12, N_R_H=10, N_theta=10)
        assert 'DT' in j._internal_absorption_comp.keys()
        assert 'BLR' in j._internal_absorption_comp.keys()
        assert 'Corona' in j._internal_absorption_comp.keys()
        assert j._internal_absorption_comp['DT']['pars']['use_sigma_gamma_gamma_fast'] is True
        assert int(j._get_internal_abs_component_on_blob('DT').use_sigma_gamma_gamma_fast) == 1
        assert int(j._get_internal_abs_component_on_blob('BLR').use_sigma_gamma_gamma_fast) == 0
        assert int(j._get_internal_abs_component_on_blob('Corona').use_sigma_gamma_gamma_fast) == 0
        assert int(j._get_internal_abs_component_on_blob('DT').is_valid) == 0
        assert int(j._get_internal_abs_component_on_blob('BLR').is_valid) == 0
        assert int(j._get_internal_abs_component_on_blob('Corona').is_valid) == 0

        tau_dt, nu_dt = j.eval_internal_absorption(comp='DT', peak=False)
        tau_blr, nu_blr = j.eval_internal_absorption(comp='BLR', peak=False)
        tau_corona, nu_corona = j.eval_internal_absorption(comp='Corona', peak=False)
        np.testing.assert_allclose(float(j._blob.core.R_H), r_h_ref, rtol=1E-12, atol=0.0)
        self._assert_total_cache_disabled(j)

        for comp_name, tau_arr, nu_arr in (
            ('DT', tau_dt, nu_dt),
            ('BLR', tau_blr, nu_blr),
            ('Corona', tau_corona, nu_corona),
        ):
            tau_arr = np.asarray(tau_arr, dtype=float)
            nu_arr = np.asarray(nu_arr, dtype=float)
            assert tau_arr.ndim == 1
            assert nu_arr.ndim == 1
            assert tau_arr.size > 0
            assert tau_arr.size == nu_arr.size
            assert np.all(np.isfinite(tau_arr))
            assert np.all(np.isfinite(nu_arr))
            assert np.all(tau_arr >= 0.0)
            c_comp = j._get_internal_abs_component_on_blob(comp_name)
            assert int(c_comp.is_enabled) == 1
            assert int(c_comp.is_valid) == 1
            assert int(c_comp.tau_size) > 0
            assert int(c_comp.peak_mode) == 0

        raw_tau_blr_peak, _ = j._internal_absorption_comp['BLR']['obj'].eval_tau_photons(
            nu_src=np.logspace(20, 29, 80),
            R_H=1.7 * j.parameters.R_H.val,
            peak=True,
            use_R_H_profile_extrapolation=j._internal_absorption_comp['BLR']['pars']['use_R_H_profile_extrapolation'],
        )
        np.testing.assert_allclose(float(j._blob.core.R_H), r_h_ref, rtol=1E-12, atol=0.0)
        self._assert_total_cache_disabled(j)
        raw_tau_blr_peak = np.asarray(raw_tau_blr_peak, dtype=float)
        assert raw_tau_blr_peak.size > 0
        assert np.all(np.isfinite(raw_tau_blr_peak))
        assert np.any(raw_tau_blr_peak > 0.0)
        c_blr_peak = j._get_internal_abs_component_on_blob('BLR')
        assert int(c_blr_peak.peak_mode) == 1
        assert int(c_blr_peak.N_soft) == 1
        assert int(c_blr_peak.is_valid) == 1

        y_ia = np.asarray(j.eval(nu=nu, get_model=True), dtype=float)
        assert np.all(np.isfinite(y_ia))
        np.testing.assert_allclose(float(j._blob.core.R_H), r_h_ref, rtol=1E-12, atol=0.0)
        self._assert_total_cache_disabled(j)
        assert int(j._get_internal_abs_component_on_blob('DT').is_valid) == 1
        assert int(j._get_internal_abs_component_on_blob('BLR').is_valid) == 1
        assert int(j._get_internal_abs_component_on_blob('Corona').is_valid) == 1

        m = (y_no_ia > 0) & (nu >= 1E22)
        assert np.any(m)
        ratio = y_ia[m] / y_no_ia[m]
        assert np.all(ratio <= 1.0 + 1E-4)
        assert np.any(ratio < 1.0 - 1E-3)

        j.remove_internal_absorption('DT')
        tau_none, nu_none = j.eval_internal_absorption(comp='DT')
        assert tau_none is None
        assert nu_none is None
        assert 'DT' not in j._internal_absorption_comp.keys()
        assert 'BLR' in j._internal_absorption_comp.keys()
        assert 'Corona' in j._internal_absorption_comp.keys()
        assert int(j._get_internal_abs_component_on_blob('DT').is_enabled) == 0
        assert int(j._get_internal_abs_component_on_blob('DT').is_valid) == 0
        self._assert_total_cache_disabled(j)

        j.remove_internal_absorption('BLR')
        j.remove_internal_absorption('Corona')
        assert len(j._internal_absorption_comp.keys()) == 0
        assert int(j._get_internal_abs_component_on_blob('BLR').is_enabled) == 0
        assert int(j._get_internal_abs_component_on_blob('BLR').is_valid) == 0
        assert int(j._get_internal_abs_component_on_blob('Corona').is_enabled) == 0
        assert int(j._get_internal_abs_component_on_blob('Corona').is_valid) == 0
        self._assert_total_cache_disabled(j)
 
    def test_internal_absorption_serialization(self, plot=False):
        from jetset.jet_model import Jet
        from jetset.internal_absorption import BlazarSED

        j = self._build_internal_abs_jet()
        assert hasattr(BlazarSED, 'eval_internal_abs_tau_isolated')

        dt_cfg = dict(
            comp='DT',
            nu_min=1E21,
            N_soft=10,
            N_hard=11,
            N_R_H=9,
            N_theta=8,
            use_R_H_profile_extrapolation=True,
            use_sigma_gamma_gamma_fast=True
        )
        blr_cfg = dict(
            comp='BLR',
            nu_min=1E20,
            N_soft=9,
            N_hard=10,
            N_R_H=8,
            N_theta=7,
            use_R_H_profile_extrapolation=False,
            use_sigma_gamma_gamma_fast=False
        )
        corona_cfg = dict(
            comp='Corona',
            nu_min=1E20,
            N_soft=8,
            N_hard=9,
            N_R_H=7,
            N_theta=6,
            use_R_H_profile_extrapolation=True,
            use_sigma_gamma_gamma_fast=True
        )

        j.enable_internal_absorption(**dt_cfg)
        j.enable_internal_absorption(**blr_cfg)
        j.enable_internal_absorption(**corona_cfg)
        j.save_model('test_internal_absorption.pkl')

        new_j = Jet.load_model('test_internal_absorption.pkl')
        assert set(new_j._internal_absorption_comp.keys()) == {'DT', 'BLR', 'Corona'}
        r_h_ref = float(new_j._blob.core.R_H)

        for cfg in (dt_cfg, blr_cfg, corona_cfg):
            comp = cfg['comp']
            p = new_j._internal_absorption_comp[comp]['pars']
            assert p['comp'] == comp
            assert p['N_soft'] == cfg['N_soft']
            assert p['N_hard'] == cfg['N_hard']
            assert p['N_R_H'] == cfg['N_R_H']
            assert p['N_theta'] == cfg['N_theta']
            assert p['use_R_H_profile_extrapolation'] == cfg['use_R_H_profile_extrapolation']
            assert p['use_sigma_gamma_gamma_fast'] == cfg['use_sigma_gamma_gamma_fast']
            np.testing.assert_allclose(p['nu_min'], cfg['nu_min'], rtol=1E-12)
            assert int(new_j._get_internal_abs_component_on_blob(comp).use_sigma_gamma_gamma_fast) == int(cfg['use_sigma_gamma_gamma_fast'])

            tau, nu_tau = new_j.eval_internal_absorption(comp=comp)
            tau = np.asarray(tau, dtype=float)
            nu_tau = np.asarray(nu_tau, dtype=float)
            assert tau.size > 0
            assert tau.size == nu_tau.size
            assert np.all(np.isfinite(tau))
            assert np.all(np.isfinite(nu_tau))
            assert np.all(tau >= 0.0)
            c_comp = new_j._get_internal_abs_component_on_blob(comp)
            assert int(c_comp.is_enabled) == 1
            assert int(c_comp.is_valid) == 1
            assert int(c_comp.peak_mode) == 0
            assert int(c_comp.N_soft) == cfg['N_soft']
            np.testing.assert_allclose(float(new_j._blob.core.R_H), r_h_ref, rtol=1E-12, atol=0.0)

            custom_nu = np.logspace(21, 28, 37)
            tau_custom, nu_custom = new_j.eval_internal_absorption(comp=comp, nu=custom_nu)
            tau_custom = np.asarray(tau_custom, dtype=float)
            nu_custom = np.asarray(nu_custom, dtype=float)
            assert tau_custom.size == custom_nu.size
            np.testing.assert_allclose(nu_custom, custom_nu, rtol=0, atol=0)
            assert np.all(np.isfinite(tau_custom))
            assert np.all(tau_custom >= 0.0)
            np.testing.assert_allclose(float(new_j._blob.core.R_H), r_h_ref, rtol=1E-12, atol=0.0)

        c_dt = new_j._get_internal_abs_component_on_blob('DT')
        c_blr = new_j._get_internal_abs_component_on_blob('BLR')
        c_corona = new_j._get_internal_abs_component_on_blob('Corona')
        assert int(c_dt.is_enabled) == 1
        assert int(c_blr.is_enabled) == 1
        assert int(c_corona.is_enabled) == 1
        self._assert_total_cache_disabled(new_j)

        y = np.asarray(new_j.eval(nu=np.logspace(20, 29, 80), get_model=True), dtype=float)
        assert np.all(np.isfinite(y))
        np.testing.assert_allclose(float(new_j._blob.core.R_H), r_h_ref, rtol=1E-12, atol=0.0)
        assert int(new_j._get_internal_abs_component_on_blob('DT').is_valid) == 1
        assert int(new_j._get_internal_abs_component_on_blob('BLR').is_valid) == 1
        assert int(new_j._get_internal_abs_component_on_blob('Corona').is_valid) == 1
        self._assert_total_cache_disabled(new_j)
