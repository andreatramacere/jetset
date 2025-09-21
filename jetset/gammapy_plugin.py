__author__ = "Andrea Tramacere"

import os

try:
    from gammapy.modeling.models import (
    SpectralModel,
)
    from gammapy.modeling.parameter import Parameter,Parameters
    from gammapy.estimators import FluxPoints
    from gammapy.datasets import FluxPointsDataset
    from gammapy.modeling import Fit

    gammapy_installed = True

except:
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
    
    if on_rtd is True:
        SpectralModel=object
        pass
    else:
        raise  ImportError('to use gammapy plugin you need to install gammapy: https://docs.gammapy.org/0.19/getting-started/install.html')

import astropy.units as u
import  numpy as np
from .model_parameters import SettingDependentParError
from .model_manager import FitModel

__all__=['GammapyJetsetModel','GammapyJetsetModelFactory']

def string_to_int(s: str) -> int:
    """Convert a string into a unique integer using base-256 encoding."""
    result = 0
    for ch in s:
        result = result * 256 + ord(ch)
    return result

def int_to_string(n: int) -> str:
    """Recover the original string from the integer representation."""
    chars = []
    while n > 0:
        chars.append(chr(n % 256))
        n //= 256
    return ''.join(reversed(chars))

class GammapyJetsetModel(SpectralModel):
    

    def __init__(self,jetset_model,clone=True):
       
        if clone is True:
            _jetset_model = jetset_model.clone()
        else:
            _jetset_model= jetset_model
        
        self._jetset_model=_jetset_model
    
         
        parameters = []
        self._parameter_values_string = {}
        self._parameters_hash_dict={}
        for ID,p in enumerate(self._jetset_model.parameters.par_array):
            
            if isinstance(self._jetset_model,FitModel):
                gp_name=f'{p.name}_{p.model.name}'
                
            else:
                gp_name=f'{p.name}'

            self._parameters_hash_dict[gp_name]=p
            if type(p.val) is str:
                self._parameter_values_string[gp_name]=p.val
                parameter = Parameter(gp_name, string_to_int(p.val), frozen=p.frozen)
            else:
                parameter = Parameter(gp_name, p.val, frozen=p.frozen)
            
            if _jetset_model.parameters.par_array[ID].units is not None:
                try:
                    parameter.unit = p.units
                except:
                    parameter.unit = ''
            else:
                  parameter.unit = ''


            if p.val_min is not None:
                parameter.min=p.val_min
            
            if p.fit_range_min is not None:
                 parameter.min=p.fit_range_min

            if p.val_max is not None:
                parameter.max=p.val_max
            
            if p.fit_range_max is not None:
                 parameter.max=p.fit_range_max

            parameter._is_jetset_dependent=p._is_dependent
            parameters.append(parameter)
            
        self.default_parameters = Parameters(parameters)
        self.tag=_jetset_model.name
        super(GammapyJetsetModel, self).__init__()

    def evaluate(self,energy=None,**kwargs):
        print(kwargs)
        if energy is None:
            el1=np.log10( self._jetset_model.nu_min)
            el2=np.log10( self._jetset_model.nu_max)
            energy=(np.logspace(el1,el2,self._jetset_model.nu_size)*u.Hz).to('eV',equivalencies=u.spectral())
        
        nu = energy.to("Hz", equivalencies=u.spectral())

        for p in self.parameters:
            if p.name not in kwargs.keys()  and not p._is_jetset_dependent:
                if p.name in self._parameter_values_string:
                    #self._jetset_model.set_par(p.name ,val=self._parameter_values_string[p.name])
                    self._parameters_hash_dict[p.name].set(val=self._parameter_values_string[p.name])
                else:
                    #self._jetset_model.set_par(p.name ,val=p.value)
                    self._parameters_hash_dict[p.name].set(val=p.value)

        for k,v in kwargs.items():
            try:
                if k in self._parameter_values_string:
                    #self._jetset_model.set_par(k ,val=self._parameter_values_string[k])
                    self._parameters_hash_dict[k].set(val=self._parameter_values_string[k])
                else:
                    #self._jetset_model.set_par(k ,val=v.value)
                    self._parameters_hash_dict[k].set(val=v.value)
            except SettingDependentParError:
                pass
            except Exception as e:
                raise(e)
          

        self._jetset_model.eval(nu=nu.value)
        if isinstance(self._jetset_model,FitModel):
            nuFnu = self._jetset_model.SED.nuFnu
            raw_flux = u.Quantity(np.asarray(nuFnu), unit=u.erg / (u.cm**2 * u.s))
            _spec = (raw_flux.to(u.eV / (u.cm**2 * u.s)) / (energy.to(u.eV)**2)).to(1 / (u.cm**2 * u.s * u.eV))
        else:
            _spec= self._jetset_model.spectral_components.Sum.SED.nuFnu.to("eV/(cm2 s)")/(energy.to('eV')**2)
        return _spec.to("1 / (cm2 eV s)")
    
    @property
    def jetset_model(self):
        return self._jetset_model
    

def GammapyJetsetModelFactory(jetset_model,clone=True):
    return GammapyJetsetModel(jetset_model,clone=clone)
