import pytest
from .base_class import TestBase

from .test_mcmc import TestEmcee
from .test_jet_model import TestJets
from .test_phenom_constr import TestPhenomenologyConstr
from .test_model_fit import TestModelFit
from .test_emitters import TestEmitters
from .test_ebl import TestEBL
from .test_mcmc import TestEmcee
from .test_depending_parameters import TestDependingParameters
from .test_composite_model import TestCompositeModel
from .test_temp_ev import TestTempEv

@pytest.fixture
def plot():
   input = False
   return input

class TestIntegration(TestBase):
  
   def test_jet(self,plot=plot):
      
      t=TestJets()
      t.integration_suite()

   def test_emitters(self,plot=plot):
      t=TestEmitters()
      t.integration_suite()

   def test_dep_pars(self,plot=plot):
      t=TestDependingParameters()
      t.integration_suite()
   
   def test_composit_model(self,plot=plot):
      t=TestCompositeModel()
      t.integration_suite()

   def test_composit_ebl(self,plot=plot):
      t=TestEBL()
      t.integration_suite()

      
   def test_model_fit(self,phenom_dict=None,plot=plot):
      t=TestPhenomenologyConstr()
      phenom_dict=t.integration_suite(sed_number=1)
      t=TestModelFit()
      fit_dict=t.integration_suite(sed_number=None,phenom_dict=phenom_dict,use_ebl=False,use_dep_pars=False,skip_minuit=True)
      return fit_dict  

   def test_emcee(self,fit_dict=None,plot=plot):
      t=TestEmcee()
      if fit_dict is None:
         sed_number=1
      else:
         sed_number=None
      t.integration_suite(fit_dict=fit_dict,sed_number=sed_number,plot=plot)

   def test_temp_ev(self,plot=plot):
      t=TestTempEv() 
      t.integration_suite(plot=plot)  


