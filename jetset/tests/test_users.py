import sys
import pytest
from .base_class import TestBase


from .test_jet_model import TestJets,hadronic_func
from .test_emitters import TestEmitters
from .test_ebl import TestEBL
from .test_depending_parameters import TestDependingParameters
from .test_composite_model import TestCompositeModel

@pytest.fixture
def plot():
   input = False
   return input


class TestUser(TestBase):
   
   def test_emitters(self,plot=plot):
      t=TestEmitters()
      t.integration_suite()

   def test_jet(self,plot=plot):
      t=TestJets()
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

   def test_jet_hadronic(self,plot=plot):
      hadronic_func(plot)
      #t=TestJetHadronic()
      #t.test_hadronic_jet(plot=plot)