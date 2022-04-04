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
from .test_composite_model import *

@pytest.mark.integration
class TestFull:

    def test_emitters_log_values(self, plot=False):
        test_emitters_log_values(plot=plot)

    def test_custom_emitters(self, plot=False):
        test_custom_emitters(plot=plot)

    def test_custom_emitters_array(self, plot=False):
        test_custom_emitters_array(plot=plot)

    def test_custom_emitters_log_values(self, plot=False):
        test_custom_emitters_log_values(plot=plot)

    def test_emitters_U(self,plot=False):
        test_emitters_U(plot=plot)

    def test_set_N_from_nuFnu(self, plot=False):
        test_set_N_from_nuFnu(plot=plot)

    def test_model_fit_lsb(self, plot=False):
        fit_model, model_minimizer, data = test_model_fit(plot=plot, sed_number=1, minimizer='lsb')


    def test_model_fit_lsb_1(self, plot=False):
        fit_model, model_minimizer, data = test_model_fit(plot=plot, sed_number=2, minimizer='lsb')


    def test_model_fit_minuit(self, plot=False):
        fit_model, model_minimizer, data = test_model_fit(plot=plot, sed_number=2, minimizer='minuit')


    def test_model_fit_dep_pars(self, plot=False):
        test_model_fit_dep_pars(plot=plot, sed_number=2, minimizer='lsb')


    def test_emcee(self, plot=False):
        test_emcee(plot=plot)


    def test_ebl_jet(self, plot=False):
        test_ebl_jet(plot=plot)


    def test_ebl(self, plot=False):
        test_ebl(plot=plot)


    def test_ebl_jet_fit(self, plot=False):
        test_ebl_jet_fit(plot=plot)


    def test_dep_par(self, plot=False):
        test_dep_par(plot=plot)


    def test_dep_par_jet(self, plot=False):
        test_dep_par_jet(plot=plot)

    def test_dep_par_log_values(self,plot=False):
        test_dep_par_log_values(plot=plot)

    def test_dep_par_composite_model(self, plot=False):
        test_dep_par_composite_model(plot=plot)


    def test_hadronic_jet(self, plot=False):
        test_hadronic_jet(plot=plot)

    def test_composite_model_pars(self, plot=False):
        test_composite_model_pars(plot=plot)