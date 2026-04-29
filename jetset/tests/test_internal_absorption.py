import numpy as np
from .base_class import TestBase


class TestInternalAbsorption(TestBase):

    def integration_suite(self, plot=False):
        self.test_internal_absorption_enable_eval_remove(plot=plot)
        self.test_internal_absorption_serialization(plot=plot)

    def _build_internal_abs_jet(self):
        from jetset.jet_model import Jet

        j = Jet(name='compact_int_abs', emitters_distribution='bkn', beaming_expr='bulk_theta')
        j.add_EC_component(EC_components_list=['EC_DT', 'EC_BLR'], disk_type='BB')

        j.parameters.z_cosm.val = 0.03
        j.parameters.L_Disk.val = 2E45
        j.parameters.R_H.val = 1E18
        j.parameters.tau_DT.val = 0.1
        j.parameters.tau_BLR.val = 0.1
        j.parameters.B.val = 0.2

        j.set_gamma_grid_size(120)
        j.set_IC_nu_size(80)
        return j

    def test_internal_absorption_enable_eval_remove(self, plot=False):
        j = self._build_internal_abs_jet()
        nu = np.logspace(20, 29, 120)

        y_no_ia = np.asarray(j.eval(nu=nu, get_model=True), dtype=float)
        assert np.all(np.isfinite(y_no_ia))

        j.enable_internal_absorption('DT', N_soft=12, N_hard=12, N_R_H=10, N_theta=10)
        j.enable_internal_absorption('BLR', N_soft=12, N_hard=12, N_R_H=10, N_theta=10)
        assert 'DT' in j._internal_absorption_comp.keys()
        assert 'BLR' in j._internal_absorption_comp.keys()

        tau_dt, nu_dt = j.eval_internal_absorption(comp='DT', peak=False)
        tau_blr, nu_blr = j.eval_internal_absorption(comp='BLR', peak=False)

        for tau_arr, nu_arr in ((tau_dt, nu_dt), (tau_blr, nu_blr)):
            tau_arr = np.asarray(tau_arr, dtype=float)
            nu_arr = np.asarray(nu_arr, dtype=float)
            assert tau_arr.ndim == 1
            assert nu_arr.ndim == 1
            assert tau_arr.size > 0
            assert tau_arr.size == nu_arr.size
            assert np.all(np.isfinite(tau_arr))
            assert np.all(np.isfinite(nu_arr))
            assert np.all(tau_arr >= 0.0)

        y_ia = np.asarray(j.eval(nu=nu, get_model=True), dtype=float)
        assert np.all(np.isfinite(y_ia))

        m = y_no_ia > 0
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

    def test_internal_absorption_serialization(self, plot=False):
        from jetset.jet_model import Jet

        j = self._build_internal_abs_jet()
        dt_cfg = dict(
            comp='DT',
            nu_min=1E21,
            N_soft=10,
            N_hard=11,
            N_R_H=9,
            N_theta=8,
            use_R_H_profile_extrapolation=True
        )
        blr_cfg = dict(
            comp='BLR',
            nu_min=1E20,
            N_soft=9,
            N_hard=10,
            N_R_H=8,
            N_theta=7,
            use_R_H_profile_extrapolation=False
        )

        j.enable_internal_absorption(**dt_cfg)
        j.enable_internal_absorption(**blr_cfg)
        j.save_model('test_internal_absorption.pkl')

        new_j = Jet.load_model('test_internal_absorption.pkl')
        assert set(new_j._internal_absorption_comp.keys()) == {'DT', 'BLR'}

        for cfg in (dt_cfg, blr_cfg):
            comp = cfg['comp']
            p = new_j._internal_absorption_comp[comp]['pars']
            assert p['comp'] == comp
            assert p['N_soft'] == cfg['N_soft']
            assert p['N_hard'] == cfg['N_hard']
            assert p['N_R_H'] == cfg['N_R_H']
            assert p['N_theta'] == cfg['N_theta']
            assert p['use_R_H_profile_extrapolation'] == cfg['use_R_H_profile_extrapolation']
            np.testing.assert_allclose(p['nu_min'], cfg['nu_min'], rtol=1E-12)

            tau, nu_tau = new_j.eval_internal_absorption(comp=comp)
            tau = np.asarray(tau, dtype=float)
            nu_tau = np.asarray(nu_tau, dtype=float)
            assert tau.size > 0
            assert tau.size == nu_tau.size
            assert np.all(np.isfinite(tau))
            assert np.all(np.isfinite(nu_tau))
            assert np.all(tau >= 0.0)

            custom_nu = np.logspace(21, 28, 37)
            tau_custom, nu_custom = new_j.eval_internal_absorption(comp=comp, nu=custom_nu)
            tau_custom = np.asarray(tau_custom, dtype=float)
            nu_custom = np.asarray(nu_custom, dtype=float)
            assert tau_custom.size == custom_nu.size
            np.testing.assert_allclose(nu_custom, custom_nu, rtol=0, atol=0)
            assert np.all(np.isfinite(tau_custom))
            assert np.all(tau_custom >= 0.0)

        c_dt = new_j._get_internal_abs_component_on_blob('DT')
        c_blr = new_j._get_internal_abs_component_on_blob('BLR')
        assert int(c_dt.is_enabled) == 1
        assert int(c_blr.is_enabled) == 1

        y = np.asarray(new_j.eval(nu=np.logspace(20, 29, 80), get_model=True), dtype=float)
        assert np.all(np.isfinite(y))
