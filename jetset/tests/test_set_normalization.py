import numpy as np

from jetset.jet_model import Jet
from jetset.jetkernel import jetkernel as BlazarSED
from jetset.jet_emitters_factory import InjEmittersFactory


class TestSetNormalization:

    @staticmethod
    def _make_jet(sync_only=False):
        j = Jet(verbose=False)
        j.parameters.z_cosm.val = 0.1

        if sync_only:
            ssc = j.get_spectral_component_by_name('SSC', verbose=False)
            if ssc is not None:
                ssc.state = 'off'

        return j

    @staticmethod
    def _integration_bounds(j):
        gmin = max(1.1 * j.parameters.gmin.val, 1.0)
        gmax = 0.5 * j.parameters.gmax.val
        if gmax <= gmin:
            gmax = j.parameters.gmax.val
        return gmin, gmax

    def test_set_N_from_U_emitters(self):
        j = self._make_jet()
        gmin, gmax = self._integration_bounds(j)

        target_u = 1e-4
        j.set_N_from_U_emitters(target_u, gmin=gmin, gmax=gmax)

        u_out = j.emitters_distribution.eval_U(gmin=gmin, gmax=gmax)
        np.testing.assert_allclose(u_out, target_u, rtol=1e-2)

    def test_set_N_from_U_vol_emitters(self):
        j = self._make_jet()
        gmin, gmax = self._integration_bounds(j)

        j.set_blob()
        target_u = 3e-5
        target_u_vol = target_u * j._blob.core.Vol_region

        j.set_N_from_U_vol_emitters(target_u_vol, gmin=gmin, gmax=gmax)

        u_out = j.emitters_distribution.eval_U(gmin=gmin, gmax=gmax)
        np.testing.assert_allclose(u_out * j._blob.core.Vol_region, target_u_vol, rtol=1e-2)

    def test_set_N_from_L_sync(self):
        j = self._make_jet()
        target_l_sync = 5e43

        j.set_N_from_L_sync(target_l_sync)

        j.set_blob()
        l_sync_out = BlazarSED.Power_Sync_Electron(j._blob) * j.get_beaming() ** 4
        np.testing.assert_allclose(l_sync_out, target_l_sync, rtol=1e-2)

    def test_set_N_from_F_sync(self):
        j = self._make_jet()
        target_f_sync = 3e-12

        j.set_N_from_F_sync(target_f_sync)

        j.set_blob()
        l_sync_out = BlazarSED.Power_Sync_Electron(j._blob) * j.get_beaming() ** 4
        f_sync_out = l_sync_out / (4.0 * np.pi * j.get_DL_cm() ** 2)
        np.testing.assert_allclose(f_sync_out, target_f_sync, rtol=1e-2)

    def test_set_N_from_nuLnu(self):
        j = self._make_jet(sync_only=True)
        nu_src = 1e14
        target_nu_lnu = 1e43

        j.set_N_from_nuLnu(target_nu_lnu, nu_src)

        nu_obs = nu_src / (1.0 + j.parameters.z_cosm.val)
        model = np.atleast_1d(j.eval(nu=np.array([nu_obs]), fill_SED=False, get_model=True))
        nu_lnu_out = model[0] * 4.0 * np.pi * j.get_DL_cm() ** 2
        np.testing.assert_allclose(nu_lnu_out, target_nu_lnu, rtol=2e-2)

    def test_set_N_from_nuFnu(self):
        j = self._make_jet(sync_only=True)
        nu_obs = 1e14
        target_nu_fnu = 1e-14

        j.set_N_from_nuFnu(target_nu_fnu, nu_obs)

        model = np.atleast_1d(j.eval(nu=np.array([nu_obs]), fill_SED=False, get_model=True))
        np.testing.assert_allclose(model[0], target_nu_fnu, rtol=2e-2)


class TestSetNormalizationLeptonicEquilibrium:

    @staticmethod
    def _make_eq_jet(sync_only=False):
        q_inj = InjEmittersFactory().create_inj_emitters(
            'pl',
            emitters_type='electrons',
            normalize=False
        )
        q_inj.parameters.p.val = 2.2
        q_inj.parameters.gmin.val = 10.0
        q_inj.parameters.gmax.val = 1e6
        q_inj.parameters.Q.val = 1e-5

        j = Jet(emitters_distribution=q_inj, emitters_type='electrons', verbose=False)
        j.parameters.z_cosm.val = 0.1
        j.parameters.T_esc_e_primaries.val = 2.0
        j.eval()

        if sync_only:
            ssc = j.get_spectral_component_by_name('SSC', verbose=False)
            if ssc is not None:
                ssc.state = 'off'

        return j

    @staticmethod
    def _integration_bounds(j):
        gmin = max(1.1 * j.parameters.gmin.val, 1.0)
        gmax = 0.5 * j.parameters.gmax.val
        if gmax <= gmin:
            gmax = j.parameters.gmax.val
        return gmin, gmax

    def test_set_N_from_U_emitters_eq(self):
        j = self._make_eq_jet()
        gmin, gmax = self._integration_bounds(j)

        target_u = 1e-4
        j.set_N_from_U_emitters(target_u, gmin=gmin, gmax=gmax)
        j.eval()

        u_out = j.emitters_distribution.eval_U(gmin=gmin, gmax=gmax)
        np.testing.assert_allclose(u_out, target_u, rtol=5e-2)

    def test_set_N_from_U_vol_emitters_eq(self):
        j = self._make_eq_jet()
        gmin, gmax = self._integration_bounds(j)

        j.set_blob()
        target_u = 5e-5
        target_u_vol = target_u * j._blob.core.Vol_region

        j.set_N_from_U_vol_emitters(target_u_vol, gmin=gmin, gmax=gmax)
        j.eval()

        u_out = j.emitters_distribution.eval_U(gmin=gmin, gmax=gmax)
        np.testing.assert_allclose(u_out * j._blob.core.Vol_region, target_u_vol, rtol=5e-2)

    def test_set_N_from_L_sync_eq(self):
        j = self._make_eq_jet()
        target_l_sync = 5e43

        j.set_N_from_L_sync(target_l_sync)

        j.set_blob()
        l_sync_out = BlazarSED.Power_Sync_Electron(j._blob) * j.get_beaming() ** 4
        np.testing.assert_allclose(l_sync_out, target_l_sync, rtol=2e-2)

    def test_set_N_from_F_sync_eq(self):
        j = self._make_eq_jet()
        target_f_sync = 3e-12

        j.set_N_from_F_sync(target_f_sync)

        j.set_blob()
        l_sync_out = BlazarSED.Power_Sync_Electron(j._blob) * j.get_beaming() ** 4
        f_sync_out = l_sync_out / (4.0 * np.pi * j.get_DL_cm() ** 2)
        np.testing.assert_allclose(f_sync_out, target_f_sync, rtol=2e-2)

    def test_set_N_from_nuLnu_eq(self):
        j = self._make_eq_jet(sync_only=True)
        nu_src = 1e14
        target_nu_lnu = 1e43

        j.set_N_from_nuLnu(target_nu_lnu, nu_src)

        nu_obs = nu_src / (1.0 + j.parameters.z_cosm.val)
        model = np.atleast_1d(j.eval(nu=np.array([nu_obs]), fill_SED=False, get_model=True))
        nu_lnu_out = model[0] * 4.0 * np.pi * j.get_DL_cm() ** 2
        np.testing.assert_allclose(nu_lnu_out, target_nu_lnu, rtol=5e-2)

    def test_set_N_from_nuFnu_eq(self):
        j = self._make_eq_jet(sync_only=True)
        nu_obs = 1e14
        target_nu_fnu = 1e-14

        j.set_N_from_nuFnu(target_nu_fnu, nu_obs)

        model = np.atleast_1d(j.eval(nu=np.array([nu_obs]), fill_SED=False, get_model=True))
        np.testing.assert_allclose(model[0], target_nu_fnu, rtol=5e-2)
