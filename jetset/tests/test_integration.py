import pytest
import os
from .base_class import TestBase

from .test_mcmc import TestEmcee
from .test_jet_model import TestJets,hadronic_func
from .test_model_fit import TestModelFit
from .test_emitters import TestEmitters
from .test_ebl import TestEBL
from .test_mcmc import TestEmcee
from .test_depending_parameters import TestDependingParameters
from .test_composite_model import TestCompositeModel
from .test_temp_ev import TestTempEv
from .test_galactic import TestGalactic

@pytest.fixture
def plot():
   input = False
   return input

class TestIntegration(TestBase):

   def test_jet(self,plot=plot):
      t=TestJets()
      t.integration_suite()
   
   @pytest.mark.skipif(os.getenv('WF_ENV')=='CONDA', reason="not running with conda") 
   def test_jet_hadronic(self,plot=plot):
       hadronic_func(plot)
   #   #t=TestJetHadronic()
   #   #t.test_hadronic_jet(plot=plot)
   
   @pytest.mark.skipif(os.getenv('WF_ENV')=='CONDA', reason="not running with conda") 
   def test_galactic(self,plot=plot):
      t=TestGalactic()
      t.integration_suite(plot=plot)

   def test_emitters(self,plot=plot):
      t=TestEmitters()
      t.integration_suite(plot=plot)

   def test_dep_pars(self,plot=plot):
      t=TestDependingParameters()
      t.integration_suite(plot=plot)
   
   def test_composit_model(self,plot=plot):
      t=TestCompositeModel()
      t.integration_suite(plot=plot)

   def test_composit_ebl(self,plot=plot):
      t=TestEBL()
      t.integration_suite(plot=plot)

      
   def test_model_fit(self,phenom_dict=None,plot=plot):
      from .test_phenom_constr import prepare_asset
      phenom_dict=prepare_asset(sed_number=1)
      t=TestModelFit()
      t.integration_suite(sed_number=None,phenom_dict=phenom_dict,use_ebl=False,use_dep_pars=False,skip_minuit=True,plot=plot)
   

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


