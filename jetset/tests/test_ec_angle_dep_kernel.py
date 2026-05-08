import numpy as np
import pytest

from .base_class import TestBase


class TestECAngleDepKernel(TestBase):

    def integration_suite(self, plot=False):
        self.test_ec_dt_angle_dep_kernel_blob_frame(plot=plot)

    @staticmethod
    def _build_ec_dt_jet(kernel_mode):
        from jetset.jet_model import Jet

        j = Jet(
            name='ec_dt_angle_dep_kernel_test',
            emitters_distribution='lp',
            beaming_expr='bulk_theta',
            verbose=False,
        )
        j.add_EC_component(EC_components_list=['EC_DT'], disk_type='BB')
        j.set_external_field_transf('blob')

        j.parameters.BulkFactor.val = 10.0
        j.parameters.theta.val = 1.0
        j.parameters.R_H.val = 2.0 * j.parameters.R_DT.val

        j.set_gamma_grid_size(121)
        j.set_seed_nu_size(100)
        j.set_IC_nu_size(100)

        if not hasattr(j._blob.core, 'EC_kernel'):
            pytest.skip('EC kernel mode not exposed in this build')

        j._blob.core.EC_kernel = int(kernel_mode)
        j._blob.core.EC_angle_n_phi = 32
        return j

    def test_ec_dt_angle_dep_kernel_blob_frame(self, plot=False):
        nu_eval = np.logspace(18, 29, 160)

        jet_iso = self._build_ec_dt_jet(kernel_mode=0)
        jet_ang = self._build_ec_dt_jet(kernel_mode=1)

        jet_iso.eval(nu=nu_eval, get_model=True)
        jet_ang.eval(nu=nu_eval, get_model=True)

        comp_iso = jet_iso.get_spectral_component_by_name('EC_DT', verbose=False)
        comp_ang = jet_ang.get_spectral_component_by_name('EC_DT', verbose=False)
        assert comp_iso is not None
        assert comp_ang is not None

        nu_iso = np.asarray(comp_iso.SED.nu, dtype=float)
        nu_ang = np.asarray(comp_ang.SED.nu, dtype=float)
        y_iso = np.asarray(comp_iso.SED.nuFnu, dtype=float)
        y_ang = np.asarray(comp_ang.SED.nuFnu, dtype=float)

        assert np.all(np.isfinite(nu_iso))
        assert np.all(np.isfinite(nu_ang))
        assert np.all(np.isfinite(y_iso))
        assert np.all(np.isfinite(y_ang))

        if nu_iso.shape != nu_ang.shape or np.any(np.abs((nu_iso - nu_ang) / nu_iso) > 1e-10):
            y_ang = np.interp(nu_iso, nu_ang, y_ang, left=np.nan, right=np.nan)

        assert np.all(np.isfinite(y_ang))
        assert np.any(y_iso > 0.0)

        peak_iso = np.nanmax(y_iso)
        assert peak_iso > 0.0

        mask = y_iso > (peak_iso * 1e-2)
        assert np.any(mask)

        rel_diff = np.abs(y_ang[mask] - y_iso[mask]) / y_iso[mask]
        assert np.nanmax(rel_diff) < 0.10
