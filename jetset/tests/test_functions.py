import pytest

from .test_mcmc import *
from .test_jet_model import *
from .test_data import *
from .test_phenom_constr import *
from .test_model_fit import *
from .test_emitters import *

def test_full(plot=False):
    from jetset.plot_sedfit import plt
    plt.ioff()
    test_jet(plot)
    test_hadronic_model(plot)
    sed_data=test_data_loader()
    print('done')
    my_shape=test_spectral_indices(sed_data)
    print('done')
    my_shape=test_sed_shaper(my_shape)
    print('done')
    prefit_jet, my_shape = test_model_constr(my_shape)
    print('done')
    jet_lsb, model_minimizer_lsb, fit_model_lsb = test_model_fit_lsb(sed_data,my_shape)
    jet_minuit,model_minimizer_minuit = test_model_fit_minuit(sed_data,my_shape)

def test_short(plot=False):
    test_jet(plot)
    sed_data = test_data_loader(plot)
    print('done')
    my_shape = test_spectral_indices(sed_data,plot)
    print('done')
    prefit_jet, my_shape = test_sed_shaper(my_shape,plot)
    print('done')
    prefit_jet = test_model_constr(my_shape,plot)
    print('done')
    jet_lsb, model_minimizer_lsb,fit_model_lsb = test_model_fit_lsb(sed_data, my_shape,plot)
    print('TEST PASSED: OK')


@pytest.mark.users
def test_users():
    test_short(plot=False)