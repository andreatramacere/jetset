import pytest

from .test_mcmc import *
from .test_jet_model import *
from .test_data import *
from .test_phenom_constr import *
from .test_model_fit import *
from .test_emitters import *
from .test_ebl import *
from .test_mcmc import test_emcee
from .test_depending_parameters import *

def test_full(plot=False):
    from jetset.plot_sedfit import plt
    plt.ioff()
    test_jet(plot=plot)
    test_hadronic_jet(plot=plot)
    fit_model, model_minimizer,data = test_model_fit (minimizer='lsb', sed_number=1)
    fit_model, model_minimizer,data = test_model_fit(minimizer='minuit', sed_number=1)
    fit_model, model_minimizer,data = test_model_fit(minimizer='lsb', sed_number=2)
    test_ebl(plot=plot)
    test_ebl_jet(plot=plot)
    test_ebl_jet_fit(plot=plot)
    test_emcee(plot=plot)
    test_dep_par(plot=plot)
    test_dep_par_jet(plot=plot)
    test_dep_par_composite_model(plot=plot)
    test_hadronic_jet(plot=plot)

def test_short(plot=False):
    test_jet(plot)
    jet_lsb, model_minimizer_lsb,sed_data = test_model_fit (minimizer='lsb', sed_number=1)
    print('TEST PASSED: OK')




@pytest.mark.integration
def test_custom(plot=False):
    test_model_fit(plot=plot,sed_number=2,minimizer='lsb')
    test_model_fit_dep_pars(plot=plot,sed_number=2,minimizer='lsb')
    test_emcee(plot=plot)
    test_custom_emitters(plot=plot)
    test_custom_emitters_array(plot=plot)
    test_ebl_jet(plot=plot)
    test_ebl(plot=plot)
    test_ebl_jet_fit(plot=plot)
    test_dep_par(plot=plot)
    test_dep_par_jet(plot=plot)
    test_dep_par_composite_model(plot=plot)
    test_hadronic_jet(plot=plot)