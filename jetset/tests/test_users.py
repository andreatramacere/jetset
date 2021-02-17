import sys
import pytest
from .test_functions import test_short

@pytest.fixture
def plot():
   input = False
   return input


@pytest.mark.users
def test_users():
    test_short(plot=False)