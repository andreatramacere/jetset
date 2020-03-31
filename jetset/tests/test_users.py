import sys
import pytest
from .test_functions import *

@pytest.fixture
def plot():
   input = False
   return input


def test_foo(plot):
   test_jet(plot)