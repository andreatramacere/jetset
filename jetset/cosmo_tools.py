from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)


__author__ = "Andrea Tramacere"




__all__=['Cosmo']


class  Cosmo (object):

    def __init__(self,astropy_cosmo=None,DL_cm=None):

        _c = None

        if DL_cm is not None and astropy_cosmo is not None:
            raise  RuntimeError('Either you provide an astropy comsology objet, or luminosity distance in cm, or nothing')

        elif astropy_cosmo is None and DL_cm is None:
            from astropy.cosmology import Planck13 as cosmo

            _c = cosmo

        elif astropy_cosmo is not None and DL_cm is None:

            _c = astropy_cosmo

        elif astropy_cosmo is None and DL_cm is not None:

            _c=None
        else:
            raise RuntimeError('Either you provide an astropy comsology objet, or luminosity distance in cm, or nothing')


        self._c = _c
        self._DL_cm = DL_cm


    def get_DL_cm(self,z):
        if self._c is not None:
            _d= self.cosmo.luminosity_distance( z ).to('cm').value
        else:
            _d = self._DL_cm

        return _d
    
    




