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
        pass
    else:

        raise  ImportError('to use gammapy plugin you need to install gammapy: https://docs.gammapy.org/0.19/getting-started/install.html')

import astropy.units as u
import  numpy as np

def GammapyJetsetModelFactory(jetmodel,clone=True):
    
    if clone is True:
        _jetset_model = jetmodel.clone()
    else:
        _jetset_model= jetmodel
    

    class GammapyJetsetModel(SpectralModel):
        """    """

        for ID,p in enumerate(_jetset_model.parameters.par_array):
            #print(p.name)
            exec("%s=Parameter(p.name,p.val, frozen=p.frozen)"%p.name)
            if _jetset_model.parameters.par_array[ID].units is not None:
                try:
                    exec("%s.unit = '%s'"%(p.name,_jetset_model.parameters.par_array[ID].units))
                except:
                    exec("%s.unit = '' "%(p.name))
            else:
                exec("%s.unit = '' "%(p.name))

            if _jetset_model.parameters.par_array[ID].val_min is not None:
                exec("%s.min = %s"%(p.name,_jetset_model.parameters.par_array[ID].val_min))

            if _jetset_model.parameters.par_array[ID].val_max is not None:
                exec("%s.max = %s"%(p.name,_jetset_model.parameters.par_array[ID].val_max)) 
       
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
    
    gammapy_obj=GammapyJetsetModel()
    setattr(gammapy_obj,'_jetset_model',_jetset_model)
    setattr(gammapy_obj,'tag',_jetset_model.name)
    return gammapy_obj