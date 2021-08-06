__author__ = "Andrea Tramacere"


from .model_parameters import ModelParameterArray, ModelParameter

import numpy as np
from astropy import constants as const
from .spectral_shapes import SED
from .base_model import  Model


__all__=['AnalyticalParameter','Disk']

class AnalyticalParameter(ModelParameter):
    """
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to analytical parameters.
    """
    def __init__(self,polymodel,**keywords):
        
        self.polymodel=polymodel

        self.allowed_par_types=['Temperature','peak freq','peak flux']

        super(AnalyticalParameter,self).__init__(  **keywords)

        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val(self.name,val)
            
        
        
    def set(self,**keywords):
        super(AnalyticalParameter,self).set(**keywords )
        
        """
        sets a parameter value checking for physical boundaries 
        """
        
        if 'val' in keywords.keys():    
            self.assign_val(self.name,keywords['val']) 
        
    
    def assign_val(self,name,val):
        """
        assigns the paramter valut to the :class:`Template` object
        """
        setattr(self.polymodel,name,val)






class Disk(Model):
    def __init__(self,cosmo,z,nu_size=100,name='disk',**keywords):
        """
        """
        
        super(Disk,self).__init__(  **keywords)
        
        self.z=z

        self.name=name

        self.model_type='analytical'
            
        self.SED = SED(name=self.model_type)
    
        self.parameters = ModelParameterArray()
      
        self.T_Disk=10000
        self.parameters.add_par(AnalyticalParameter(self,name='T_Disk',par_type='Temperature',val=10000,val_min=0.,val_max=None,units=''))
        
       
        self.nuFnu_p_D=1E-10
        self.parameters.add_par(AnalyticalParameter(self,name='nuFnu_p',par_type='peak flux',val=1E-10,val_min=0.,val_max=1E-3,units='erg cm^-2 s^-1',log=True))
        
        self.HPLANCK=const.h.cgs.value
        self.SIGSB=const.sigma_sb.cgs.value
        self.vluce_cm=const.c.cgs.value
        self.K_boltz=const.k_B.cgs.value


        self.cosmo=cosmo

        self.DL=self.cosmo.get_DL_cm(self.z)
        
        self.L_D=self.get_L_D()
        self.nu_p_Disk=self.get_nu_p()
        
        

    def set_disk_pars(self,fit_model,model_name):
        
        self.L_Disk=self.get_L_D()
          
        self.L_Disk_err=self.L_Disk*(fit_model.parameters.get(model_name,'nuFnu_p','best_fit_err')/
                                     fit_model.parameters.get(model_name,'nuFnu_p','best_fit_val'))

                        
        self.T_Disk=fit_model.parameters.get(model_name,'T_Disk','best_fit_val')
        self.T_Disk_err=fit_model.parameters.get(model_name,'T_Disk','best_fit_err')

        self.nu_p_Disk=self.get_nu_p()
        self.nu_p_Disk_err=self.nu_p_Disk*(self.T_Disk_err/self.T_Disk)

    
    def get_nu_p(self):
        return self.T_Disk*(1.39*5.879e10)

    
    def get_L_D(self):
        return 4*np.pi*self.DL*self.DL*self.nuFnu_p
        
        
    #def log_func(self,log_nu):
    #    return np.log10(self.lin_func(np.power(10,log_nu)))

    
    def lin_func(self,nu):
        nu=nu*(1+self.z)
        a = 2 * self.HPLANCK * nu*nu*nu*nu / (self.vluce_cm* self.vluce_cm)
        a *= 1.0 / (np.exp((self.HPLANCK * nu) / (self.K_boltz * self.T_Disk)) - 1)
        return a/(self.SIGSB*self.T_Disk*self.T_Disk*self.T_Disk*self.T_Disk/np.pi)*self.nuFnu_p
    
    
    def log_func(self,log_nu):
        
        return np.log10(self.lin_func(np.power(10,log_nu)))