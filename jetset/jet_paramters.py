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
                 allowed_values=None ):

        self.ptype =ptype
        self.vmin =vmin
        self.vmax =vmax
        self.punit =punit
        self.froz =froz
        self.log =log
        self.allowed_values =allowed_values


class JetParameter(ModelParameter):
    """

    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to  handles SSC/EC parameters,
    overriding the :meth:`.ModelParameter.set` in order to propagate the
    parameter value to the BlazarSED object instance


    """
    def __init__(self,jet,**keywords):

        self._jet=jet

        self.par_type_list = ['proton_density',
                              'electron_density',
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

        if 'par_type' in keywords.keys() and keywords['par_type'] not in self.par_type_list:
            msg= "blob_parameter type %s not allowed"%keywords['par_type'] + "\n please choose among %s"%self.par_type_list
            raise ValueError("%s"%msg)



        super(JetParameter,self).__init__(  **keywords)


        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val(self.name,val)



    #OVERRIDES Base Method

    @safe_run
    def set(self,**keywords):
        """
        overrides the  :meth:`.ModelParameter.set` method in order to propagate the
        parameter value to the BlazarSED object instance
        """
        super(JetParameter,self).set(**keywords )


        if 'val' in keywords.keys():

            self.assign_val(self.name,keywords['val'])




    def assign_val(self,name,val):
        """
        sets the :class:`.JetParameter` value in the BlazarSED object
        """
        #TODO improve this with allowed values: Done!
        #if name=='disk_type':

        #    if val not in allowed_disk_type:

        #        raise RuntimeError('Disk type %s, not in allowed '%val,allowed_disk_type)
        #print('-> assign val',name,val)
        b=getattr(self._jet._blob,name)


        #print('1 name',name,val,self._val.islog)
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
        #print('2 name',name, val, self._val.islog)
        setattr(self._jet._blob,name,val)