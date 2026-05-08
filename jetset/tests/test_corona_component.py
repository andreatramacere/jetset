import numpy as np

from .base_class import TestBase


class TestCoronaComponent(TestBase):

    def integration_suite(self, plot=False):
        self.test_corona_component(plot=plot)
        self.test_corona_low_energy_cutoff_parameter()

    def test_corona_component(self, plot=False):
        from jetset.jet_model import Jet

        j = Jet(name='test_corona_component', emitters_distribution='lp')
        j.add_EC_component(['EC_Corona'])

        j.parameters.L_Corona.val = 5E44
        j.parameters.R_Corona.val = 5E14
        j.parameters.alpha_Corona.val = 1.1
        j.parameters.nu_cut_Corona.val = 1E20

        nu = np.logspace(12, 30, 120)
        j.parameters.R_H.val = 1E16
        j.parameters.R_H_Corona.val = 5E15
        y_far_above = np.asarray(j.eval(nu=nu, get_model=True), dtype=float)
        assert np.all(np.isfinite(y_far_above))
        assert np.any(y_far_above > j._blob.core.emiss_lim)

        j.parameters.R_H.val = 1E15
        j.parameters.R_H_Corona.val = 6E15
        y_far_below = np.asarray(j.eval(nu=nu, get_model=True), dtype=float)
        assert np.all(np.isfinite(y_far_below))
        assert np.any(y_far_below > j._blob.core.emiss_lim)

        assert j.get_spectral_component_by_name('Corona', verbose=False) is not None
        assert j.get_spectral_component_by_name('EC_Corona', verbose=False) is not None

        n_max = int(j._blob.Corona.spec.NU_INT_MAX) + 1
        nu_drf = np.asarray(j._blob.Corona.spec.nu_DRF, dtype=float)[:n_max]
        l_nu_drf = np.asarray(j._blob.Corona.spec.L_nu_DRF, dtype=float)[:n_max]
        f_nu = l_nu_drf / j.parameters.L_Corona.val
        area = np.trapz(f_nu, nu_drf)
        np.testing.assert_allclose(area, 1.0, rtol=5E-2, atol=1E-3)

    def test_corona_low_energy_cutoff_parameter(self):
        from jetset.jet_model import Jet

        j = Jet(name='test_corona_low_energy_cutoff', emitters_distribution='lp')
        j.add_EC_component(['EC_Corona'])

        j.parameters.L_Corona.val = 1E45
        j.parameters.R_Corona.val = 5E15
        j.parameters.alpha_Corona.val = 1.1
        j.parameters.nu_cut_Corona.val = 1E20
        j.parameters.R_H.val = 1E16
        j.parameters.R_H_Corona.val = 5E15

        j.parameters.nu_cut_low_Corona.val = 0.0
        j.eval(nu=np.logspace(12, 30, 120), get_model=True)
        n_max = int(j._blob.Corona.spec.NU_INT_MAX) + 1
        nu_no_cut = np.asarray(j._blob.Corona.spec.nu_DRF, dtype=float)[:n_max]
        l_no_cut = np.asarray(j._blob.Corona.spec.L_nu_DRF, dtype=float)[:n_max]

        j.parameters.nu_cut_low_Corona.val = 1E18
        j.eval(nu=np.logspace(12, 30, 120), get_model=True)
        n_max = int(j._blob.Corona.spec.NU_INT_MAX) + 1
        nu_with_cut = np.asarray(j._blob.Corona.spec.nu_DRF, dtype=float)[:n_max]
        l_with_cut = np.asarray(j._blob.Corona.spec.L_nu_DRF, dtype=float)[:n_max]

        low_nu = 1E16
        high_nu = 1E19
        l_low_no_cut = np.interp(low_nu, nu_no_cut, l_no_cut)
        l_high_no_cut = np.interp(high_nu, nu_no_cut, l_no_cut)
        l_low_with_cut = np.interp(low_nu, nu_with_cut, l_with_cut)
        l_high_with_cut = np.interp(high_nu, nu_with_cut, l_with_cut)

        ratio_no_cut = l_low_no_cut / l_high_no_cut
        ratio_with_cut = l_low_with_cut / l_high_with_cut
        assert ratio_with_cut < ratio_no_cut * 1E-4
