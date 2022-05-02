import pytest
import numpy as np
def test(plot=False):
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
    print('--------> j.energetic_report')
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
    print('-------->  j_new.energetic_report')
    j_new.energetic_report(verbose=False)
    assert('U_p_cold' not in j.energetic_dict.keys())   
    np.testing.assert_allclose(j.energetic_dict['U_p'],j.emitters_distribution.eval_U(),rtol=1E-2)