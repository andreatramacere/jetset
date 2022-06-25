import pytest
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from .base_class import TestBase

class TestGalactic(TestBase):

    def integration_suite(self,plot=False):
        self.test_galactic_leptonic_unbeamed(plot=plot)
        self.test_galactic_hadronic_unbeamed(plot=plot)
        self.test_galactic_hadronic_beamed(plot=plot)

    def test_galactic_leptonic_unbeamed(self,plot=False):
        from jetset.jet_model import GalacticUnbeamed

        pwn=GalacticUnbeamed(emitters_distribution='bkn',verbose=False,emitters_type='electrons',distance=2*u.pc,name='pwn')
        pwn.add_EC_component('EC_CMB')
        pwn.parameters.R.val=(2*u.pc).to('cm').value
        pwn.parameters.B.val=1E-4
        pwn.parameters.gamma_break.val=1E8
        pwn.parameters.gmax.val=8E9
        pwn.parameters.gmin.val=5E5
        pwn.parameters.p.val=3.2
        pwn.parameters.p_1.val=3.8

        pwn.parameters.N.val=1E-8
        pwn.show_model()

        pwn.eval()
        if plot is True:
            p=pwn.plot_model(frame='src')
        
        pwn.energetic_report()
        assert('NH_cold_to_rel_e'  in pwn.energetic_dict.keys())
        _par_array=pwn._build_energetic_dict()
        np.testing.assert_allclose(pwn.energetic_dict['U_Synch'],pwn.energetic_dict['U_Synch_DRF'],rtol=1E-3)
        pwn.save_model('pwn.pkl')
        GalacticUnbeamed.load_model('pwn.pkl')

    def test_galactic_hadronic_unbeamed(self,plot=False):
        from jetset.jet_model import GalacticUnbeamed

        gal_hadronic=GalacticUnbeamed(emitters_distribution='plc',verbose=False,emitters_type='protons',distance=2*u.pc,name='gal_hadronic_unbeamed')
        gal_hadronic.parameters.R.val=1E18
        gal_hadronic.parameters.N.val=1000
        gal_hadronic.parameters.B.val=1E-3
        gal_hadronic.show_model()

        gal_hadronic.eval()
        if plot is  True:
            p=gal_hadronic.plot_model(frame='src')
            p.setlim(y_min=1E33)
        
        gal_hadronic.energetic_report()
        gal_hadronic.save_model('gal_hadronic.pkl')
        GalacticUnbeamed.load_model('gal_hadronic.pkl')

    def test_galactic_hadronic_beamed(self,plot=False):
        from jetset.jet_model import GalacticBeamed

        gal_hadronic=GalacticBeamed(emitters_distribution='plc',verbose=False,emitters_type='protons',distance=2*u.pc,name='gal_hadronic_beamed')
        gal_hadronic.parameters.R.val=1E18
        gal_hadronic.parameters.N.val=1000
        gal_hadronic.parameters.B.val=1E-3

        gal_hadronic.show_model()

        gal_hadronic.eval()
        
        if plot is True:
            p=gal_hadronic.plot_model(frame='src')
            p.setlim(y_min=1E36)

        gal_hadronic.energetic_report()
        gal_hadronic.save_model('gal_hadronic_beamed.pkl')
        GalacticBeamed.load_model('gal_hadronic_beamed.pkl')
