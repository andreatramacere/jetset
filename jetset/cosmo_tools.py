from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)


__author__ = "Andrea Tramacere"

from astropy.units import  Unit as u


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

    def __repr__(self):
        if self._c is not None:
            s= '%s' % self._c
        elif self._DL_cm is not None:
            x=self._DL_cm*u('cm')

            s= 'cosmology is not defined, the luminosity distance has been set to %s'%x
        else:
            s='you have to define either an astropy comoslogy model or set a luminosity distance in cm'

        return s

    def get_DL_cm(self,z=None):
        if self._c is not None:
            _d= self._c.luminosity_distance( z ).to('cm').value
        else:
            _d = self._DL_cm

        return _d
    
    




