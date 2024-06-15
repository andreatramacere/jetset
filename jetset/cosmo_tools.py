__author__ = "Andrea Tramacere"

from astropy.units import  Unit as u
from astropy.units.quantity import Quantity
from astropy import cosmology 
from astropy.table import Table
from astropy.cosmology import Cosmology

import warnings

__all__=['Cosmo']


class  Cosmo (object):


    def __init__(self,astropy_cosmo=None,DL_cm=None,verbose=False):
        
        _c = None
        self._c_name=None
        if DL_cm is not None and astropy_cosmo is not None:
            raise  RuntimeError('Either you provide an astropy cosmology objet, or luminosity distance in cm, or nothing')

        elif  astropy_cosmo is None and DL_cm is None:
            _c = self._get_custom_cosmology()
            if verbose is True:
                print("cosmology set to custom",str(_c))
            

        elif astropy_cosmo is not None and DL_cm is None:
            _c = astropy_cosmo
            
        elif astropy_cosmo is None and DL_cm is not None:

            _c=None
            if verbose is True:
                print("using cosmo without z and only DL, should be used only for galactic objects!!")
                print("z will be fixed to zero")
        else:
            raise RuntimeError('Either you provide an astropy cosmology objet, or luminosity distance in cm, or nothing')


        self._c = _c
        if DL_cm is not None:
            if isinstance(DL_cm,Quantity):
                if u(DL_cm)=='cm':
                    pass
                else:
                    DL_cm=DL_cm.to('cm')
            else:
                DL_cm=DL_cm*u('cm')
        self._DL_cm = DL_cm

    def __repr__(self):
        if self._c is not None:
            s= '%s' % self._c
        elif self._DL_cm is not None:
            x=self._DL_cm

            s= 'cosmology is not defined, the luminosity distance has been set to %s'%x
        else:
            s='you have to define either an astropy cosmology model or set a luminosity distance in cm'

        return s

    def get_DL_cm(self,z=None):
        if self._c is not None:
            #THIS IS FIXING THE ERROR WITH PICKLED COSMO
            #TODO: open issue on astropy!
            try:
                _d= self._c.luminosity_distance( z ).to('cm').value
            except:
                _d = self._c.luminosity_distance(z)
                _d = _d.value*u(str(_d.unit))
                _d = _d.to('cm').value
        else:
            _d = self._DL_cm.value

        return _d
    
    def _serialize_model(self):
        _model = {}
        if self._c is not None:
            _model['_astropy_cosmo']=Table(self._c.to_format('astropy.table'))
        else:
            _model['_astropy_cosmo']=None
            _model['_DL_cm'] = self._DL_cm.value
        return   _model 
          
    def __getstate__(self):
        return  self._serialize_model()

    def __setstate__(self,state):
        astropy_cosmo,DL_cm=self._decode_model(state)
        self.__init__(astropy_cosmo=astropy_cosmo,DL_cm=DL_cm)
    
    @classmethod
    def from_model(cls,model):
        astropy_cosmo,DL_cm  = cls._decode_model(model)
        return cls(astropy_cosmo,DL_cm)

    @staticmethod
    def _get_custom_cosmology():
        try:
            return cosmology.Planck13
        except:
            return cosmology.FlatLambdaCDM(H0=67.8, Om0=0.307)


    @staticmethod
    def _decode_model(model):
        DL_cm=None
        astropy_cosmo=None
        try:
            if '_astropy_cosmo' in model.keys():
                if model['_astropy_cosmo'] is not None:
                    try:
                        astropy_cosmo=Cosmology.from_format(Table(model['_astropy_cosmo']),format='astropy.table')
                    except Exception as e:
                        warnings.warn('unable to get astropy.cosmology from loaded instance, setting to Planck13')
                        astropy_cosmo = cosmology.Planck13
            if '_DL_cm' in model.keys():
                DL_cm=model['_DL_cm']
        except Exception as e:
            warnings.warn('failed to decode saved astropy model, reason: %s'%str(e))       
        return astropy_cosmo,DL_cm  

