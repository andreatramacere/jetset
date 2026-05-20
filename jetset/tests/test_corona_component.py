import numpy as np

from .base_class import TestBase


class TestCoronaComponent(TestBase):

    @staticmethod
    def _get_corona_src_sed_arrays(jet):
        corona = jet.get_spectral_component_by_name('Corona', verbose=False)
        assert corona is not None

        # Read source-frame luminosity arrays through the public spectral-component API.
        corona.fill_SED(skip_zeros=True)
        nu_src = np.asarray(corona.SED.nu_src.value, dtype=float)
        nu_lnu_src = np.asarray(corona.SED.nuLnu_src.value, dtype=float)

        m = np.isfinite(nu_src) & np.isfinite(nu_lnu_src) & (nu_src > 0.0) & (nu_lnu_src > 0.0)
        nu_src = nu_src[m]
        nu_lnu_src = nu_lnu_src[m]
        assert nu_src.size > 1
        return nu_src, nu_lnu_src

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

        nu_src, nu_lnu_src = self._get_corona_src_sed_arrays(j)
        l_nu_src = nu_lnu_src / nu_src
        f_nu = l_nu_src / j.parameters.L_Corona.val
        area = np.trapezoid(f_nu, nu_src)
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
        nu_no_cut, nu_lnu_no_cut = self._get_corona_src_sed_arrays(j)
        l_no_cut = nu_lnu_no_cut / nu_no_cut

        j.parameters.nu_cut_low_Corona.val = 1E18
        j.eval(nu=np.logspace(12, 30, 120), get_model=True)
        nu_with_cut, nu_lnu_with_cut = self._get_corona_src_sed_arrays(j)
        l_with_cut = nu_lnu_with_cut / nu_with_cut

        low_nu = 1E16
        high_nu = 1E19
        l_low_no_cut = np.interp(low_nu, nu_no_cut, l_no_cut)
        l_high_no_cut = np.interp(high_nu, nu_no_cut, l_no_cut)
        l_low_with_cut = np.interp(low_nu, nu_with_cut, l_with_cut)
        l_high_with_cut = np.interp(high_nu, nu_with_cut, l_with_cut)

        ratio_no_cut = l_low_no_cut / l_high_no_cut
        ratio_with_cut = l_low_with_cut / l_high_with_cut
        assert ratio_with_cut < ratio_no_cut * 1E-4
