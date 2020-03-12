from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)




__author__ = "Andrea Tramacere"

from .model_parameters import ModelParameterArray, ModelParameter
from .utils import safe_run



__all__=['JetParameter','JetModelDictionaryPar']


class JetModelDictionaryPar(object):

    def __init__(self,
                 ptype=None,
                 vmin=None,
                 vmax=None,
                 punit=None,
                 froz=False,
                 log=False,
                 val=None,
                 is_in_blob=True,
                 allowed_values=None ):

        self.ptype =ptype
        self.val=val
        self.vmin =vmin
        self.vmax =vmax
        self.punit =punit
        self.froz =froz
        self.log =log
        self.is_in_blob=is_in_blob
        self.allowed_values =allowed_values


class JetParameter(ModelParameter):
    """
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to  handles SSC/EC parameters,
    overriding the :meth:`.ModelParameter.set` in order to propagate the
    parameter value to the BlazarSED object instance
    """
    def __init__(self,jet, is_in_blob=True, **keywords):

        self._jet=jet

        self.allowed_par_types = ['emitters_density',
                                  'target_density',
                                  'Disk',
                                  'BLR',
                                  'DT',
                                  'region_size',
                                  'region_position',
                                  'electron_energy',
                                  'LE_spectral_slope',
                                  'HE_spectral_slope',
                                  'high-energy-cut-off',
                                  'low-energy-cut-off',
                                  'spectral_curvature',
                                  'turn-over-energy',
                                  'magnetic_field',
                                  'beaming',
                                  'jet-viewing-angle',
                                  'jet-bulk-factor',
                                  'redshift']

        self._is_in_blob=is_in_blob
        super(JetParameter,self).__init__(  **keywords)

        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val_to_blob(self.name, val)

    @safe_run
    def set(self,**keywords):
        """
        overrides the  :meth:`.ModelParameter.set` method in order to propagate the
        parameter value to the BlazarSED object instance
        """
        super(JetParameter,self).set(**keywords )

        if 'val' in keywords.keys():
            self.assign_val_to_blob(self.name, keywords['val'])



    def assign_val_to_blob(self, name, val):
        """
        sets the :class:`.JetParameter` value in the BlazarSED object
        """

        if self._is_in_blob is True:
            if hasattr(self._jet._blob,name):
                b=getattr(self._jet._blob,name)



                if type(b)==int:
                    if self._val.islog is True:
                        val=10**val
                    val=int(val)

                elif type(b)==float:
                    if self._val._islog is True:
                        val=10**val
                    val=float(val)

                elif type(b)==str:
                    val=val

                setattr(self._jet._blob,name,val)
