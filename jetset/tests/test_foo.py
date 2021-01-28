import pytest

def test_my_foo():
    from jetset.jet_model import Jet
    j = Jet()
    j.eval()
    j.energetic_report()
