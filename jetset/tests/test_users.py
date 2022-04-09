import sys
import pytest
from .test_composite_model import TestCompositeModel
from .test_jet_model import TestJets
from .test_depending_parameters import TestDependingParameters

@pytest.fixture
def plot():
   input = False
   return input


def test(plot):
   for TestClass in [TestJets,TestDependingParameters,TestCompositeModel]:
      t=TestClass()
      t.test_all(plot=plot)
   