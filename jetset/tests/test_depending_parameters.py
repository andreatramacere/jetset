import pytest


def test_dep_par_distr(plot=False):
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

    def par_func(s1=n_e_bkn.parameters.s1):
        return s1.val + 1

    n_e_bkn.parameters.s2.make_dependent(par_func)

    np.testing.assert_allclose(n_e_bkn.parameters.s2.val , n_e_bkn.parameters.s1.val+1)

    from jetset.jet_model import Jet

    j = Jet(emitters_distribution=n_e_bkn)

    def par_func(s1=j.parameters.s1):
        return s1.val + 1

    j.parameters.s2.make_dependent(par_func)
    j.parameters.s1.val=3
    np.testing.assert_allclose(j.parameters.s2.val, j.parameters.s1.val + 1)
    j.save_model('jet.pkl')
    j_new=Jet.load_model('jet.pkl')
    j_new.parameters.s1.val = 2
    np.testing.assert_allclose(j_new.parameters.s2.val, j_new.parameters.s1.val + 1)
    j_new.eval()
    j_new.show_model()