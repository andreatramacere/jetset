import pytest

from .test_mcmc import TestEmcee
from .test_jet_model import TestJets
#from .test_data import *
from .test_phenom_constr import TestPhenomenologyConstr
from .test_model_fit import TestModelFit
from .test_emitters import TestEmitters
from .test_ebl import TestEBL
from .test_mcmc import TestEmcee
from .test_depending_parameters import TestDependingParameters
from .test_composite_model import *

@pytest.fixture
def plot():
   input = False
   return input

def test(plot=plot):
    
    t=TestJets()
    t.integration_suite()

    t=TestEmitters()
    t.integration_suite()

    t=TestDependingParameters()
    t.integration_suite()
    
    t=TestCompositeModel()
    t.integration_suite()

    t=TestEBL()
    t.integration_suite()

    t=TestPhenomenologyConstr()
    phenom_dict=t.integration_suite(sed_number=1)

    t=TestModelFit()
    fit_dict=t.integration_suite(phenom_dict=phenom_dict,use_ebl=False,use_dep_pars=False)
    
    t=TestEmcee()
    t.integration_suite(fit_dict,plot=plot)


