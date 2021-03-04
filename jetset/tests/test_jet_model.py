import pytest
import  numpy as np

def test_hadronic_jet(plot=False):
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
    sum1=j.spectral_components.Sum.SED.nuFnu
    if plot is True:
        j.plot_model()
    j.save_model('test_jet_hadronic.pkl')
    j_new = Jet.load_model('test_jet_hadronic.pkl')
    j_new.eval()
    sum2 = j_new.spectral_components.Sum.SED.nuFnu
    np.testing.assert_allclose(sum2,sum1, rtol=1E-5)

def test_build_bessel():
    from jetset.jet_model import Jet
    Jet().eval()




@pytest.mark.users
def test_jet(plot=False):
    print('--------> test_jet',plot)
    from jetset.jet_model import Jet
    j=Jet()
    j.eval()
    j.energetic_report()
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

def test_set_N_from_nuFnu(plot=False):
    from jetset.jet_model import Jet
    from jetset.jetkernel import jetkernel
    import numpy as np
    nu = 1E15
    nuFnu = 1E-15
    j = Jet()
    j.set_N_from_nuFnu(nuFnu_obs=nuFnu,nu_obs=nu)
    y = j.eval(nu=[nu], get_model=True)
    np.testing.assert_allclose(y, nuFnu, rtol=1E-2)