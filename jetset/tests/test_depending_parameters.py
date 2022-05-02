import pytest
from .base_class import TestBase

class TestDependingParameters(TestBase):

    def integration_suite(self,plot=False):
        self._all(plot=plot)

    def test_dep_par(self,plot=False):
        from jetset.jet_emitters import EmittersDistribution
        import numpy as np

        def distr_func_bkn(gamma_break, gamma, s1, s2):
            return np.power(gamma, -s1) * (1. + (gamma / gamma_break)) ** (-(s2 - s1))

        n_e_bkn = EmittersDistribution('bkn', spectral_type='bkn')
        n_e_bkn.add_par('gamma_break', par_type='turn-over-energy', val=1E3, vmin=1., vmax=None, unit='lorentz-factor')
        n_e_bkn.add_par('s1', par_type='LE_spectral_slope', val=2.5, vmin=-10., vmax=10, unit='')
        n_e_bkn.add_par('s2', par_type='LE_spectral_slope', val=3.2, vmin=-10., vmax=10, unit='')
        n_e_bkn.set_distr_func(distr_func_bkn)
        n_e_bkn.parameters.show_pars()
        n_e_bkn.parameters.s1.val = 2.0
        n_e_bkn.parameters.s2.val = 3.5
        n_e_bkn.update()
        n_e_bkn.parameters.show_pars()

        from jetset.jet_model import Jet

        j = Jet(emitters_distribution=n_e_bkn)

        # def par_func(s1):
        #    return s1+1
        j.make_dependent_par(par='s2', depends_on=['s1'], par_expr='s1+1')
        print('here')
        j.parameters.s1.val = 3
        print('done')
        np.testing.assert_allclose(j.parameters.s2.val, j.parameters.s1.val + 1)
        j.save_model('jet.pkl')
        new_jet = Jet.load_model('jet.pkl')
        print('here')
        new_jet.parameters.s1.val = 2
        print('done')

        np.testing.assert_allclose(new_jet.parameters.s2.val, new_jet.parameters.s1.val + 1)
        j.eval()
        new_jet.show_model()


    def test_dep_par_log_values(self,plot=False):
        from jetset.jet_emitters import EmittersDistribution
        import numpy as np

        def distr_func_bkn(gamma_break, gamma, s1, s2):
            return np.power(gamma, -s1) * (1. + (gamma / gamma_break)) ** (-(s2 - s1))

        n_e_bkn = EmittersDistribution('bkn', spectral_type='bkn')
        n_e_bkn.add_par('gamma_break', par_type='turn-over-energy', val=3, vmin=0.1, vmax=None, unit='lorentz-factor',log=True)
        n_e_bkn.add_par('s1', par_type='LE_spectral_slope', val=2.5, vmin=-10., vmax=10, unit='')
        n_e_bkn.add_par('s2', par_type='LE_spectral_slope', val=3.2, vmin=-10., vmax=10, unit='')
        n_e_bkn.set_distr_func(distr_func_bkn)
        n_e_bkn.parameters.show_pars()
        n_e_bkn.parameters.s1.val = 2.0
        n_e_bkn.parameters.s2.val = 3.5
        n_e_bkn.update()
        n_e_bkn.parameters.show_pars()

        from jetset.jet_model import Jet

        j = Jet(emitters_distribution=n_e_bkn)

        # def par_func(s1):
        #    return s1+1

        j.make_dependent_par(par='gamma_break', depends_on=['gmax'], par_expr='gmax / 10')
        j.parameters.gmax.val = 1E5
        print('par',j.parameters.gamma_break.val, j.parameters.gmax.val / 10)
        np.testing.assert_allclose(j.parameters.gamma_break.val, np.log10(j.parameters.gmax.val/ 10))
        j.save_model('jet.pkl')
        new_jet = Jet.load_model('jet.pkl')
        new_jet.parameters.gmax.val = 1E4
        np.testing.assert_allclose(new_jet.parameters.gamma_break.val, np.log10(new_jet.parameters.gmax.val/ 10))
        j.eval()
        new_jet.show_model()

        j = Jet(emitters_distribution='lppl',emitters_distribution_log_values=True)
        j.show_pars()
        j.parameters.gmax.val = 5
        j.make_dependent_par(par='gamma0_log_parab', depends_on=['gmax'], par_expr='gmax / 10')
        j.show_pars()  
        
        print('par',j.parameters.gamma0_log_parab.val, j.parameters.gmax.val - 1)
        np.testing.assert_allclose(j.parameters.gamma0_log_parab.val, j.parameters.gmax.val - 1)
        print(j.parameters.gamma0_log_parab.val, j.parameters.gmax.val -1)
        j.save_model('jet.pkl')
        new_jet = Jet.load_model('jet.pkl')
        new_jet.parameters.gmax.val = 4
        np.testing.assert_allclose(new_jet.parameters.gamma0_log_parab.val, new_jet.parameters.gmax.val -1)
        j.eval()
        new_jet.show_model()


    def test_dep_par_jet(self,plot=False):
        import numpy as np

        from jetset.jet_model import Jet
        j = Jet(emitters_distribution='plc')
        j.add_user_par(name='B0',units='G',val=1E-5,val_min=0,val_max=None)
        j.make_dependent_par(par='B',depends_on=['B0'],par_expr='B0*5')
        np.testing.assert_allclose(j.parameters.B.val, j.parameters.B0.val*5)

        j.save_model('test.pkl')
        new_j=Jet.load_model('test.pkl')
        new_j.parameters.B0.val=1
        np.testing.assert_allclose(new_j.parameters.B.val, new_j.parameters.B0.val*5)


    def test_dep_par_composite_model(self,plot=False):
        import numpy as np
        from jetset.model_manager import FitModel
        from jetset.jet_model import Jet

        jet = Jet(emitters_distribution='plc')
        fit_model = FitModel(jet=jet, name='SSC-best-fit-lsb', template=None)
        fit_model.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
        fit_model.jet_leptonic.parameters.R_H.val = 5E17
        fit_model.jet_leptonic.parameters.R_H.frozen = False
        fit_model.jet_leptonic.parameters.R_H.fit_range = [1E15, 1E19]
        fit_model.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
        fit_model.jet_leptonic.add_user_par(name='B0', units='G', val=1E3, val_min=0, val_max=None)
        fit_model.jet_leptonic.add_user_par(name='R0', units='cm', val=5E13, val_min=0, val_max=None)
        fit_model.jet_leptonic.add_user_par(name='m_B', val=1, val_min=1, val_max=2)

        fit_model.jet_leptonic.parameters.R0.frozen = True
        fit_model.jet_leptonic.parameters.B0.frozen = True

        def par_func(R0,B0,R_H,m_B): 
            return B0*np.power((R0/R_H),m_B)

        fit_model.jet_leptonic.make_dependent_par(par='B', depends_on=['B0', 'R0', 'R_H','m_B'], par_expr=par_func)

        B0=fit_model.jet_leptonic.parameters.B0.val
        R0 = fit_model.jet_leptonic.parameters.R0.val
        R_H = fit_model.jet_leptonic.parameters.R_H.val
        m_B= fit_model.jet_leptonic.parameters.m_B.val

        np.testing.assert_allclose(fit_model.jet_leptonic.parameters.B.val, par_func(B0,R0,R_H,m_B))

        fit_model.save_model('test_composite.pkl')
        fit_model.parameters._build_par_table()
        t=fit_model.parameters._par_table
        t.sort(t.colnames)
        
        for new_fit_model in [FitModel.load_model('test_composite.pkl'),fit_model.clone()]:
            
            new_fit_model.parameters._build_par_table()
            _t=new_fit_model.parameters._par_table
            _t.sort(_t.colnames)

            assert(all([_t[i]['frozen']==t[i]['frozen'] for i in range(len(_t))]))
            assert(all([_t[i]['val']==t[i]['val'] for i in range(len(_t))]))

            new_fit_model.jet_leptonic.parameters.B0.val=1E4

            B0 = new_fit_model.jet_leptonic.parameters.B0.val
            R0 = new_fit_model.jet_leptonic.parameters.R0.val
            R_H = new_fit_model.jet_leptonic.parameters.R_H.val
            
            np.testing.assert_allclose(new_fit_model.jet_leptonic.parameters.B.val, par_func(B0,R0,R_H,m_B))

            