import pytest
import numpy as np
from .base_class import TestBase
from jetset.jet_model import Jet

class TestJets(TestBase):

    def integration_suite(self,plot=False):
        self.test_jet(plot=plot)
        self.test_build_bessel(plot=plot)
        
        self.test_set_N_from_nuFnu(plot=plot)
        
    def test_build_bessel(self,plot=False):
        print('--------> test_build_bessel',plot)
        
        Jet().eval()

    def test_jet(self,plot=False):
        print('--------> test_jet',plot)
       
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

        j.energetic_report()
        assert('U_p_target' not in j.energetic_dict.keys())
        np.testing.assert_allclose(j.energetic_dict['U_e'],j.emitters_distribution.eval_U(),rtol=1E-2)
        j.save_model('test_jet.pkl')
        
        j_new=Jet.load_model('test_jet.pkl')
        j_new.eval()
        sum2 = j_new.spectral_components.Sum.SED.nuFnu
        np.testing.assert_allclose(sum2, sum1, rtol=1E-5)
        j_new.energetic_report()
        assert('U_p_target' not in j_new.energetic_dict.keys())
        np.testing.assert_allclose(j_new.energetic_dict['U_e'],j_new.emitters_distribution.eval_U(),rtol=1E-2)

    def test_set_N_from_nuFnu(self,plot=False):
        print('--------> test_set_N_from_nuFnu',plot)
      
        nu = 1E15
        nuFnu = 1E-15
        j = Jet()
        j.set_N_from_nuFnu(nuFnu_obs=nuFnu,nu_obs=nu)
        y = j.eval(nu=[nu], get_model=True)
        np.testing.assert_allclose(y, nuFnu, rtol=1E-2)

    def test_synch_pol(self):
        p=2.0
        jet=Jet(name='test',electron_distribution='plc')
        jet.spectral_components.SSC.state='off'
        jet.parameters.p.val=p
        nu_range=np.logspace(10,22,50)
        p_nu,nuF_nu=jet.eval_synch_pol(nu_range)
        th_pol_power_law=(p+1)/(p+7/3)
        np.testing.assert_allclose(p_nu[:10],th_pol_power_law,rtol=0.1,atol=0.1)

class TestJetHadronic(TestBase):

    def integration_suite(self,plot=False):
        self.test_hadronic_jet(plot=plot)

    def test_hadronic_jet(self,plot=False):
        hadronic_func(plot)

def hadronic_func(plot):
    print('--------> test_hadronic_jet',plot)
        
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
    j.energetic_report(verbose=False)
    assert('U_p_cold' not in j.energetic_dict.keys())
    np.testing.assert_allclose(j.energetic_dict['U_p'],j.emitters_distribution.eval_U(),rtol=1E-2)

    sum1=j.spectral_components.Sum.SED.nuFnu
    if plot is True:
        j.plot_model()
    
    j.save_model('test_jet_hadronic.pkl')
    
    
    j_new = Jet.load_model('test_jet_hadronic.pkl')
    j_new.eval()
    sum2 = j_new.spectral_components.Sum.SED.nuFnu
    np.testing.assert_allclose(sum2,sum1, rtol=1E-5)
    j_new.energetic_report(verbose=False)
    assert('U_p_cold' not in j_new.energetic_dict.keys())   
    np.testing.assert_allclose(j_new.energetic_dict['U_p'],j_new.emitters_distribution.eval_U(),rtol=1E-2)