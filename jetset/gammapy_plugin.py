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

__all__=['GammapyJetsetModel','GammapyJetsetModelFactory']

class GammapyJetsetModel(SpectralModel):
    

    def __init__(self,jetset_model,clone=True):
       
        if clone is True:
            _jetset_model = jetset_model.clone()
        else:
            _jetset_model= jetset_model
        
        self._jetset_model=_jetset_model
    
        self._jetset_model.add_user_par(name='fake_norm',units='',val=1,val_min=0,val_max=None)
        self._jetset_model.parameters.fake_norm.frozen=True
        parameters = []
        
        for ID,p in enumerate(self._jetset_model.parameters.par_array):
            #print(p.name)
            if p.name=='fake_norm':
                is_norm=True
            else:
                is_norm=False

            parameter = Parameter(p.name, p.val, is_norm=is_norm,frozen=p.frozen)
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

            parameters.append(parameter)
        self.default_parameters = Parameters(parameters)
        self.tag=_jetset_model.name
        super(GammapyJetsetModel, self).__init__()

    def evaluate(self,energy=None,**kwargs):

        if energy is None:
            el1=np.log10( self._jetset_model.nu_min)
            el2=np.log10( self._jetset_model.nu_max)
            energy=(np.logspace(el1,el2,self._jetset_model.nu_size)*u.Hz).to('eV',equivalencies=u.spectral())
        
        nu = energy.to("Hz", equivalencies=u.spectral())

        for p in self.parameters:
            if p.name not in kwargs.keys():
                self._jetset_model.set_par(p.name ,val=p.value)

        for k,v in kwargs.items():
            self._jetset_model.set_par(k,val=v.value)

        self._jetset_model.eval(nu=nu.value)
        _spec= self._jetset_model.spectral_components.Sum.SED.nuFnu.to('eV cm-2 s-1')/(energy.to('eV')**2)
        return _spec.to("1 / (cm2 eV s)")
    
    @property
    def jetset_model(self):
        return self._jetset_model
    

def GammapyJetsetModelFactory(jetset_model,clone=True):
    return GammapyJetsetModel(jetset_model,clone=clone)