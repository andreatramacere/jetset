import pytest
from .base_class import TestBase

from .test_emcee import TestEmcee
from .test_ultranest import TestUltranest
from .test_jet_model import TestJets,hadronic_func
from .test_model_fit import TestModelFit
from .test_emitters import TestEmitters
from .test_ebl import TestEBL
from .test_depending_parameters import TestDependingParameters
from .test_composite_model import TestCompositeModel
from .test_temp_ev import TestTempEv
from .test_galactic import TestGalactic
from .test_set_normalization import TestSetNormalization, TestSetNormalizationLeptonicEquilibrium
from .test_internal_absorption import TestInternalAbsorption
from .test_corona_component import TestCoronaComponent
from .test_data import TestData
from .test_phenom_constr import TestPhenomenologyConstr

@pytest.fixture
def plot():
   input = False
   return input

class TestIntegration(TestBase):

   def test_jet(self,plot=plot):
      t=TestJets()
      t.integration_suite()

   def test_jet_synch_pol(self):
      t=TestJets()
      t.test_synch_pol()
   
   #@pytest.mark.skipif(os.getenv('WF_ENV')=='CONDA', reason="not running with conda") 
   def test_jet_hadronic(self,plot=plot):
       hadronic_func(plot)

   #@pytest.mark.skipif(os.getenv('WF_ENV')=='CONDA', reason="not running with conda") 
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

   def test_composit_ebl_fit(self,plot=plot):
      t=TestEBL()
      t.test_ebl_jet_fit(plot=plot,sed_number=2,minimizer='lsb')

      
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

   def test_ultranest(self,fit_dict=None,plot=plot):
      t=TestUltranest()
      if fit_dict is None:
         sed_number=1
      else:
         sed_number=None
      t.integration_suite(fit_dict=fit_dict,sed_number=sed_number,plot=plot)

   def test_temp_ev(self,plot=plot):
      t=TestTempEv() 
      t.integration_suite(plot=plot)

   def test_set_normalization(self):
      t=TestSetNormalization()
      t.test_set_N_from_U_emitters()
      t.test_set_N_from_U_vol_emitters()
      t.test_set_N_from_L_sync()
      t.test_set_N_from_F_sync()
      t.test_set_N_from_nuLnu()
      t.test_set_N_from_nuFnu()

   def test_set_normalization_leptonic_equilibrium(self):
      t=TestSetNormalizationLeptonicEquilibrium()
      t.test_set_N_from_U_emitters_eq()
      t.test_set_N_from_U_vol_emitters_eq()
      t.test_set_N_from_L_sync_eq()
      t.test_set_N_from_F_sync_eq()
      t.test_set_N_from_nuLnu_eq()
      t.test_set_N_from_nuFnu_eq()

   def test_internal_absorption(self,plot=plot):
      t=TestInternalAbsorption()
      t.integration_suite(plot=plot)

   def test_corona_component(self,plot=plot):
      t=TestCoronaComponent()
      t.integration_suite(plot=plot)

   def test_data(self,plot=plot):
      t=TestData()
      t.integration_suite(plot=plot)

   def test_phenom_constr(self,plot=plot):
      t=TestPhenomenologyConstr()
      t.integration_suite(sed_number=1,plot=plot)
