import pytest
import  numpy as np
from .base_class import TestBase


class TestJets(TestBase):

    def integration_suite(self,plot=False):
        self.test_build_bessel(plot=plot)
        self.test_jet(plot=plot)
        self.test_set_N_from_nuFnu(plot=plot)
        self.test_EC(plot=plot)
        self.test_hadronic_jet(plot=plot)

    def test_hadronic_jet(self,plot=False):
        print('--------> test_hadronic_jet',plot)
        from jetset.jet_model import Jet
        j = Jet(proton_distribution='plc')
        j.parameters.gmin.val = 2
        j.parameters.gmax.val = 1E8
        j.parameters.NH_pp.val = 1E10
        j.parameters.N.val = 1E1
        j.parameters.B.val = 80

        j.parameters.p.val = 2.5
        j.eval()
        j.show_model()

        j.eval()
        #print('--------> j.energetic_report')
        #j.energetic_report(verbose=False)
        #assert('U_p_cold' not in j.energetic_dict.keys())
        #np.testing.assert_allclose(j.energetic_dict['U_p'],j.emitters_distribution.eval_U(),rtol=1E-2)

        sum1=j.spectral_components.Sum.SED.nuFnu
        if plot is True:
            j.plot_model()
        j.save_model('test_jet_hadronic.pkl')
        j_new = Jet.load_model('test_jet_hadronic.pkl')
        j_new.eval()
        sum2 = j_new.spectral_components.Sum.SED.nuFnu
        np.testing.assert_allclose(sum2,sum1, rtol=1E-5)
        #print('-------->  j_new.energetic_report')
        #j_new.energetic_report(verbose=False)
        #assert('U_p_cold' not in j.energetic_dict.keys())   
        #np.testing.assert_allclose(j.energetic_dict['U_p'],j.emitters_distribution.eval_U(),rtol=1E-2)
        

    def test_build_bessel(self,plot=False):
        print('--------> test_build_bessel',plot)
        from jetset.jet_model import Jet
        Jet().eval()

    def test_jet(self,plot=False):
        print('--------> test_jet',plot)
        from jetset.jet_model import Jet
        j=Jet()
        j.eval()
       


        sum1 = j.spectral_components.Sum.SED.nuFnu
        if plot is True:
            j.plot_model()
            j.emitters_distribution.plot()
            j.emitters_distribution.plot2p()
            j.emitters_distribution.plot3p()
            j.emitters_distribution.plot3p(energy_unit='eV')
            j.emitters_distribution.plot3p(energy_unit='erg')
        j.save_model('test_jet.pkl')
        j_new=Jet.load_model('test_jet.pkl')
        j_new.eval()
        sum2 = j_new.spectral_components.Sum.SED.nuFnu
        np.testing.assert_allclose(sum2, sum1, rtol=1E-5)
        
        j=j_new
        j.energetic_report()
        assert('U_p_target' not in j.energetic_dict.keys())
        np.testing.assert_allclose(j.energetic_dict['U_e'],j.emitters_distribution.eval_U(),rtol=1E-2)

    def test_set_N_from_nuFnu(self,plot=False):
        print('--------> test_set_N_from_nuFnu',plot)
        from jetset.jet_model import Jet
        from jetset.jetkernel import jetkernel
        import numpy as np
        nu = 1E15
        nuFnu = 1E-15
        j = Jet()
        j.set_N_from_nuFnu(nuFnu_obs=nuFnu,nu_obs=nu)
        y = j.eval(nu=[nu], get_model=True)
        np.testing.assert_allclose(y, nuFnu, rtol=1E-2)

       
    def test_EC(self,plot=False):
        print('--------> test_EC', plot)