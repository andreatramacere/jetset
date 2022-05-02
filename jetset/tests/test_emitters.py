import  pytest
from .base_class import TestBase

class TestEmitters(TestBase):

    def integration_suite(self,plot=False):
        self._all(plot=plot)

    def test_custom_emitters(self,plot=True):
        from jetset.jet_model import Jet

        from jetset.jet_emitters import EmittersDistribution
        import numpy as np
        def distr_func_bkn(gamma_break, gamma, s1, s2):
            return np.power(gamma, -s1) * (1. + (gamma / gamma_break)) ** (-(s2 - s1))

        n_e = EmittersDistribution('custom_bkn',spectral_type='bkn')
        n_e.add_par('gamma_break', par_type='turn-over-energy', val=1E3, vmin=1., vmax=None, unit='lorentz-factor')
        n_e.add_par('s1', par_type='LE_spectral_slope', val=2.5, vmin=-10., vmax=10, unit='')
        n_e.add_par('s2', par_type='LE_spectral_slope', val=3.2, vmin=-10., vmax=10, unit='')
        n_e.set_distr_func(distr_func_bkn)
        n_e.parameters.show_pars()
        n_e.parameters.s1.val = 2.0
        n_e.parameters.s2.val = 3.5
        if plot is True:
            n_e.plot()

        my_jet= Jet(emitters_distribution=n_e)
        my_jet.Norm_distr = True
        my_jet.parameters.N.val = 5E4
        my_jet.eval()
        np.testing.assert_allclose(my_jet.emitters_distribution.eval_N(), my_jet.parameters.N.val, rtol=1E-5)
        print(my_jet.emitters_distribution.eval_N(), my_jet.parameters.N.val)
        print(n_e.eval_N(), my_jet.parameters.N.val)
        assert (my_jet.emitters_distribution.emitters_type=='electrons')
        my_jet.save_model('test_jet_custom_emitters.pkl')
        my_jet = Jet.load_model('test_jet_custom_emitters.pkl')
        my_jet.eval()


    def test_custom_emitters_log_values(self,plot=True):
        from jetset.jet_model import Jet

        from jetset.jet_emitters import EmittersDistribution
        import numpy as np
        def distr_func_bkn(gamma_break, gamma, s1, s2):
            return np.power(gamma, -s1) * (1. + (gamma / gamma_break)) ** (-(s2 - s1))

        n_e = EmittersDistribution('custom_bkn',spectral_type='bkn',log_values=True)
        n_e.add_par('gamma_break', par_type='turn-over-energy', val=3, vmin=np.log10(2), vmax=np.log10(1E6), unit='lorentz-factor',log=True)
        n_e.add_par('s1', par_type='LE_spectral_slope', val=2.5, vmin=-10., vmax=10, unit='')
        n_e.add_par('s2', par_type='LE_spectral_slope', val=3.2, vmin=-10., vmax=10, unit='')
        n_e.set_distr_func(distr_func_bkn)
        n_e.parameters.show_pars()
        n_e.parameters.s1.val = 2.0
        n_e.parameters.s2.val = 3.5
        if plot is True:
            n_e.plot()

        my_jet= Jet(emitters_distribution=n_e)
        my_jet.Norm_distr = True
        my_jet.parameters.N.val = 5E4
        my_jet.eval()
        np.testing.assert_allclose(my_jet.emitters_distribution.eval_N(), my_jet.parameters.N.val, rtol=1E-5)
        print(my_jet.emitters_distribution.eval_N(), my_jet.parameters.N.val)
        print(n_e.eval_N(), my_jet.parameters.N.val)
        assert (my_jet.emitters_distribution.emitters_type=='electrons')
        my_jet.save_model('test_jet_custom_emitters.pkl')
        my_jet = Jet.load_model('test_jet_custom_emitters.pkl')
        my_jet.eval()

    def test_emitters_log_values(self,plot=True):
        from jetset.jet_model import Jet

        from jetset.jet_emitters_factory import  EmittersFactory
        ef=EmittersFactory()
        for distr_name in EmittersFactory.available_distributions_list():
            n_e=ef.create_emitters(name=distr_name,log_values=True)
            my_jet = Jet(emitters_distribution=n_e)
            my_jet.eval()
            my_jet.save_model('test_jet_custom_emitters_log.pkl')
            my_jet = Jet.load_model('test_jet_custom_emitters_log.pkl')
            my_jet.eval()

    def test_custom_emitters_array(self,plot=True):
        from jetset.jet_model import Jet
        from jetset.jet_emitters import EmittersArrayDistribution
        import numpy as np

        # gamma array
        gamma = np.logspace(1, 8, 500)

        # gamma array this is n(\gamma) in 1/cm^3/gamma
        n_gamma = gamma ** -2 * 1E-5 * np.exp(-gamma / 1E5)

        N1 = np.trapz(n_gamma, gamma)

        n_distr = EmittersArrayDistribution(name='array_distr', emitters_type='protons', gamma_array=gamma, n_gamma_array=n_gamma,normalize=False)

        N2 = np.trapz(n_distr._array_n_gamma, n_distr._array_gamma)

        j = Jet(emitters_distribution=n_distr, verbose=False)

        j.parameters.z_cosm.val = z = 0.001
        j.parameters.beam_obj.val = 1
        j.parameters.N.val = 1
        j.parameters.NH_pp.val = 1
        j.parameters.B.val = 0.01
        j.parameters.R.val = 1E18
        j.set_IC_nu_size(100)
        j.gamma_grid_size = 200

        N3 = np.trapz(j.emitters_distribution.n_gamma_p, j.emitters_distribution.gamma_p)

        np.testing.assert_allclose(N1, N2, rtol=1E-5)
        np.testing.assert_allclose(N1, N3, rtol=1E-2)
        np.testing.assert_allclose(N1, j.emitters_distribution.eval_N(), rtol=1E-2)
        assert (j.emitters_distribution.emitters_type == 'protons')

        j.eval()
        j.save_model('test_jet_custom_emitters_array.pkl')
        j = Jet.load_model('test_jet_custom_emitters_array.pkl')
        j.eval()


    def test_emitters_U(self,plot=True):
        from jetset.jet_model import Jet
        from jetset.jetkernel import jetkernel
        import numpy as np

        j = Jet(emitters_distribution='plc', verbose=False, emitters_type='protons')
        j.parameters.z_cosm.val = z = 0.001
        j.parameters.beam_obj.val = 1
        j.parameters.gamma_cut.val = 1000 / (jetkernel.MPC2_TeV)
        j.parameters.NH_pp.val = 1
        j.parameters.N.val = 1
        j.parameters.p.val = 2.0
        j.parameters.B.val = 1.0
        j.parameters.R.val = 1E18
        j.parameters.gmin.val = 1
        j.parameters.gmax.val = 1E8
        j.set_emiss_lim(1E-60)
        j.set_IC_nu_size(100)
        j.gamma_grid_size = 200
        gmin=1.0/jetkernel.MPC2_TeV
        j.set_N_from_U_emitters(1.0, gmin=gmin)
        m = j.emitters_distribution.gamma_p>gmin
        N1 = jetkernel.MPC2*np.trapz(j.emitters_distribution.n_gamma_p[m]*j.emitters_distribution.gamma_p[m],j.emitters_distribution.gamma_p[m])
        N2 = j.emitters_distribution.eval_U(gmin=gmin)
        np.testing.assert_allclose(N1, N2, rtol=1E-5)