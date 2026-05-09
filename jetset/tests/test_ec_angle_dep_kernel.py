import numpy as np
import pytest

from .base_class import TestBase


class TestECAngleDepKernel(TestBase):

    def integration_suite(self, plot=False):
        self.test_ec_dt_blob_vs_disk_notebook_latest(plot=plot)

    @staticmethod
    def _build_ec_dt_jet():
        from jetset.jet_model import Jet

        j = Jet(
            name='ec_dt_angle_dep_kernel_test',
            electron_distribution='bkn',
            electron_distribution_log_values=False,
            beaming_expr='bulk_theta',
            verbose=False,
        )
        j.add_EC_component(['EC_DT'], disk_type='BB')
        j.set_par('N', val=1.0)
        j.set_par('p', val=2.5)
        j.set_par('p_1', val=3.0)
        j.set_par('gamma_break', val=1.0e3)
        j.set_par('gmin', val=1.0)
        j.set_par('gmax', val=1.0e4)
        j.set_par('R', val=1.0e15)
        j.set_par('B', val=0.1)
        j.set_par('BulkFactor', val=10.0)
        j.set_par('theta', val=20.0)
        j.set_par('z_cosm', val=0.1)
        j.set_par('L_Disk', val=1.0e45)
        j.set_par('tau_DT', val=0.1)
        j.set_par('R_DT', val=1.0e18)
        j.set_par('T_DT', val=600.0)
        j.parameters.R_H.val = j.parameters.R_DT.val / 1000.0
        j.electron_distribution.update()

        j.set_gamma_grid_size(200)
        j.set_nu_grid(1.0e14, 1.0e29, 100)
        j._blob.core.theta_size_seed_fields = 200
        if hasattr(j.spectral_components, 'Sync'):
            j.spectral_components.Sync.state = 'off'
        if hasattr(j.spectral_components, 'SSC'):
            j.spectral_components.SSC.state = 'off'
        if hasattr(j.spectral_components, 'EC_DT'):
            j.spectral_components.EC_DT.state = 'on'

        if not hasattr(j._blob.core, 'EC_kernel'):
            pytest.skip('EC kernel mode not exposed in this build')
        if not hasattr(j._blob.core, 'theta_size_seed_fields'):
            pytest.skip('theta_size_seed_fields not exposed in Python interface')

        j._blob.core.theta_size_seed_fields = 20
        j._blob.core.EC_angle_n_phi = 10
        return j

    @staticmethod
    def _eval_blob_to_disk_ratio(jet, disk_kernel, blob_kernel):
        jet._blob.core.EC_kernel = int(disk_kernel)
        jet.set_external_field_transf('disk')
        jet.eval()
        comp_disk = jet.get_spectral_component_by_name('EC_DT', verbose=False)
        assert comp_disk is not None
        y_disk = np.asarray(comp_disk.SED.nuFnu, dtype=float)
        assert np.all(np.isfinite(y_disk))

        mask = (y_disk > np.nanmax(y_disk) * 1e-2) & (y_disk > 0.0)
        assert np.any(mask)

        jet._blob.core.EC_kernel = int(blob_kernel)
        jet.set_external_field_transf('blob')
        jet.eval()
        comp_blob = jet.get_spectral_component_by_name('EC_DT', verbose=False)
        assert comp_blob is not None
        y_blob = np.asarray(comp_blob.SED.nuFnu, dtype=float)
        assert np.all(np.isfinite(y_blob))

        ratio = y_blob[mask] / y_disk[mask]
        assert np.all(np.isfinite(ratio))
        return ratio

    @staticmethod
    def _assert_ratio_near_unity(ratio, label):
        rmin = float(np.nanmin(ratio))
        rmax = float(np.nanmax(ratio))
        assert rmin > 0.9, f'{label}: min ratio={rmin}'
        assert rmax < 1.1, f'{label}: max ratio={rmax}'

    def test_ec_dt_blob_vs_disk_notebook_latest(self, plot=False):
        jet = self._build_ec_dt_jet()

        
        ratio_k1_vs_k1 = self._eval_blob_to_disk_ratio(jet, disk_kernel=1, blob_kernel=1)
        #self._assert_ratio_near_unity(ratio_k1_vs_k1, 'k1 disk vs k1 blob')

       
        ratio_k0_vs_k1 = self._eval_blob_to_disk_ratio(jet, disk_kernel=0, blob_kernel=1)
        #self._assert_ratio_near_unity(ratio_k0_vs_k1, 'k0 disk vs k1 blob')
        print(ratio_k1_vs_k1)
        print(ratio_k0_vs_k1)
