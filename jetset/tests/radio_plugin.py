from jetset.analytical_model import AnalyticalParameter
from jetset.model_parameters import ModelParameterArray
from jetset.base_model import Model
from jetset.spectral_shapes import SED



class RadioSpectrum(Model):
    def __init__(self,nu_size=100,cosmo=None,**keywords):
        """
        """
       
        super(RadioSpectrum,self).__init__(
              nu_size=nu_size,
              cosmo=cosmo,
              nu_min=1E6,
              nu_max=1E12,
              **keywords)

        self.name='radio_spectrum'
        self.model_type='radio_spectrum'



        self.SED = SED(name=self.model_type)

        self.parameters = ModelParameterArray()

        self.parameters.add_par(AnalyticalParameter(self,name='alpha_radio',par_type='spectral-slope',val=0.0,val_min=-10.,val_max=10.,units=''))
        self.parameters.add_par(AnalyticalParameter(self,name='nu_ssa',par_type='turn-over freq',val=1E10,val_min=1E6,val_max=1E12,units='Hz'))
        self.parameters.add_par(AnalyticalParameter(self,name='nu_cut',par_type='',val=1E12,val_min=1E6,val_max=1E12,units='Hz'))
        self.parameters.add_par(AnalyticalParameter(self,name='nuFnu_p',par_type='flux-const',val=1E-13,val_min=1E-30,val_max=1E-5,units='erg cm2 s-1'))


    def _func(self,nu,nu_ssa,s2,nu_cut):
        return np.power(nu,2.5)*(1.+(nu/nu_ssa))**(-(s2+2.5))*np.exp(-(nu/nu_cut))*nu

    def lin_func(self,nu):
        s2=self.parameters.get_par_by_name('alpha_radio').val
        nu_ssa=self.parameters.get_par_by_name('nu_ssa').val
        nuFnu_p=self.parameters.get_par_by_name('nuFnu_p').val
        nu_cut=self.parameters.get_par_by_name('nu_cut').val
        res= self._func(nu,nu_ssa,s2,nu_cut)
        _Norm=nuFnu_p/self._func(nu_cut,nu_ssa,s2,nu_cut)
        return res*_Norm
