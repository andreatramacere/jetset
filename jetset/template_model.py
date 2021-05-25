__author__ = "Andrea Tramacere"



from .data_loader import log_to_lin, lin_to_log

from scipy.interpolate import interp1d


import numpy as np

import os

from .spectral_shapes import SED
from .plot_sedfit import PlotSED,PlotSpecComp

from .model_parameters import ModelParameter,ModelParameterArray
from .base_model import  Model

__all__=['SpectralTemplateLogLog', 'TemplateParameter', 'BigBlueBumpTemplateLogLog', 'HostGalaxyTemplateLogLog']

class TemplateParameter(ModelParameter):
    """
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to Template parameters.
    """
    def __init__(self,template,**keywords):
        
        self.template=template

        self.allowed_par_types=['nu-scale','nuFnu-scale']


        super(TemplateParameter,self).__init__(  **keywords)

        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val(self.name,val)
            
        
        

    def set(self,**keywords):
        super(TemplateParameter,self).set(**keywords )
        
        """
        sets a parameter value checking for physical boundaries 
        """

        if 'val' in keywords.keys():
            self.assign_val(self.name,keywords['val']) 
        
    
    def assign_val(self,name,val):

        setattr(self.template,name,val)
        


class SpectralTemplateLogLog(Model):
    """
    Class to handle spectral templates
    """
    def __init__(self,template_type,cosmo,z=None,nu_size=100,name='TemplateModel'):
        """
        """
        super(SpectralTemplateLogLog, self).__init__(name=name)


        
        if z is not None  :
            self.z=z
            self.z_scale=np.log10(1+z)
        else:
            self.z=0.0
            self.z_scale=0.0

        self._scale = 'log-log'

        self.allowed_templates= self.get_allowed_template_name()
        
        self.nu_size=nu_size
        
        self.name = template_type
        
        self.model_type='template'
            
        self.SED = SED(name=self.name)
    
        self.parameters = ModelParameterArray()
      
        self.cosmo=cosmo

        self.DL=self.cosmo.get_DL_cm(self.z)
        self.flux_plot_lim=1E-30

    @staticmethod
    def get_allowed_template_name():
        return ['BBB','host_galaxy','user_defined']

    @classmethod
    def template_factory(cls,template_type,cosmo,z=None,nu_size=100,name='TemplateModel'):
        if template_type == 'BBB':
            return BigBlueBumpTemplateLogLog(template_type, cosmo, z=z, nu_size=nu_size, name=name)


        elif template_type == 'host_galaxy':

            return HostGalaxyTemplateLogLog(template_type, cosmo, z=z, nu_size=nu_size, name=name)

        else:
            raise ValueError("Wrong template type=%s, allowed=" % (template_type, cls.get_allowed_template_name()))

    def plot_model(self,plot_obj=None,clean=False,label=None,sed_data=None,color=None, density=False,frame='obs'):
        plot_obj = self._set_up_plot(plot_obj, sed_data, frame, density)

        if clean is True:
            plot_obj.clean_model_lines()


        if label is None:
            label=self.name

        plot_obj.add_model_plot(self.SED, line_style='-', label=label, flim=self.flux_plot_lim,color=color, density=density, frame=frame)

        return plot_obj


    def get_redshift(self):
        return self.z

    def set_Lum(self,nuFnu_p):
        self.L_D=self.get_L_D(nuFnu_p)

    
    def get_Lum(self,nuFnu_p):
        return 4*np.pi*self.DL*self.DL*nuFnu_p
        


    def load_teplate_from_file(self, file_path):
        """
        Loads a template from a file
        """
        xy = np.loadtxt(file_path, usecols=(0, 1), unpack=False, dtype=float)
        xy = xy[np.argsort(xy[:, 0])]
        self.nu_template = xy[:, 0]
        self.nuFnu_template = xy[:, 1]
        
   
    def log_func(self,nu_log):
        
        x_shift=getattr(self,self.x_scale)
        y_shift=getattr(self,self.y_scale)
        if np.shape(nu_log)==():
                nu=np.array([nu_log])
        
        
        x_log=nu_log-x_shift

        model=np.ones(x_log.size)+np.log10(self.flux_plot_lim)
            
        msk = x_log >self.nu_template.min() 
        msk*= x_log <self.nu_template.max()
            
        
        model[msk] = self.interp_func(x_log[msk])+y_shift
        
        return model
    
    
    def eval(self,fill_SED=True,nu=None,get_model=False,loglog=False):
        """
        Evaluates the Template for the current parameters values
        """    

        if nu is None:

            nu = np.copy(self.nu_template)


        if loglog is False:
            log_nu=np.log10(nu)

            lin_nu=nu
        else:
            log_nu=nu

            lin_nu=np.power(10.,log_nu)
        
        
        log_model= self.log_func(log_nu)
        
        model=np.power(10.,log_model)

        if fill_SED is True:
            self._fill(lin_nu = lin_nu, lin_model = model)

        if get_model is  True:
            if loglog is False:
            
                return model
            else:
                
                return log_model

        else:
            return None
            
  
class BigBlueBumpTemplateLogLog(SpectralTemplateLogLog):


    def __init__(self,template_type,cosmo,z=None,nu_size=100,name='TemplateModel'):

        super(BigBlueBumpTemplateLogLog, self).__init__(template_type,
                                                     cosmo,
                                                     z=z,
                                                     nu_size=nu_size,
                                                     name=name)

        self.Templates_dir = os.path.dirname(__file__) + '/Spectral_Templates_Repo'



        file_path = '%s/QSO_BBB_Template.dat' % self.Templates_dir

        self.load_teplate_from_file(file_path)
        self.nu_p_template = 10 ** self.nu_template[self.nuFnu_template.argmax()]
        # set nu_template according to redshift
        self.nu_template = self.nu_template - self.z_scale

        # add parameter
        self.nuFnu_p_BBB = 1.0
        self.nu_scale = 0.0
        self.parameters.add_par(
            TemplateParameter(self, name='nuFnu_p_BBB', par_type='nuFnu-scale', val=0.0, val_min=-20.0, val_max=20.0,
                              units='erg cm^-2 s^-1'))
        self.parameters.add_par(
            TemplateParameter(self, name='nu_scale', par_type='nu-scale', val=0.0, val_min=-2.0, val_max=2.0,
                              units='Hz'))
        self.y_scale = 'nuFnu_p_BBB'
        self.x_scale = 'nu_scale'
        self.nuLnu_p_BBB = self.get_Lum(self.nuFnu_p_BBB)
        self.T_BBB = self.get_T_BBB()
        self.interp_func = interp1d(self.nu_template, self.nuFnu_template)

    def get_T_BBB(self):
        return self.nu_p_template / (1.39 * 5.879e10)

    def set_BBB_pars(self, fit_model, model_name):
        self.DL = self.cosmo.get_DL_cm(self.z)

        nuFnu_p, nuFnu_p_err = log_to_lin(
            log_val=fit_model.parameters.get(model_name, 'nuFnu_p_BBB', 'best_fit_val'),
            log_err=fit_model.parameters.get(model_name, 'nuFnu_p_BBB', 'best_fit_err'))
        err_rel = nuFnu_p_err / nuFnu_p
        self.nuLnu_p_BBB = self.get_Lum(nuFnu_p)
        self.nuLnu_p_BBB_err = self.nuLnu_p_BBB * err_rel

        self.nu_p, self.nu_p_err = log_to_lin(
            log_val=fit_model.parameters.get(model_name, 'nu_scale', 'best_fit_val'),
            log_err=fit_model.parameters.get(model_name, 'nu_scale', 'best_fit_err'))
        err_rel = self.nu_p_err / self.nu_p
        self.nu_p = self.nu_p_template * self.nu_p
        self.nu_p_err = self.nu_p * err_rel

        self.T_Disk = self.get_T_BBB()
        self.T_Disk_err = self.T_Disk * err_rel


class HostGalaxyTemplateLogLog(SpectralTemplateLogLog):

    def __init__(self, template_type, cosmo, z=None, nu_size=100, name='TemplateModel'):

        super(HostGalaxyTemplateLogLog, self).__init__(template_type,
                                                     cosmo,
                                                     z=z,
                                                     nu_size=nu_size,
                                                     name=name)



        self.Templates_dir = os.path.dirname(__file__) + '/Spectral_Templates_Repo'

        file_path = '%s/HostgalaxyTemplate.dat' % self.Templates_dir

        self.load_teplate_from_file(file_path)
        self.nu_p_template = 10 ** self.nu_template[self.nuFnu_template.argmax()]
        # set nu_template according to redshift
        self.nu_template = self.nu_template - self.z_scale

        # add parameter
        self.nuFnu_p_host = 1.0
        self.nu_scale = 0.0

        self.parameters.add_par(
            TemplateParameter(self, name='nuFnu_p_host', par_type='nuFnu-scale', val=0.0, val_min=-20.0, val_max=20.0,
                              units='erg cm^-2 s^-1'))
        self.parameters.add_par(
            TemplateParameter(self, name='nu_scale', par_type='nu-scale', val=0.0, val_min=-2.0, val_max=2.0,
                              units='Hz'))

        self.y_scale = 'nuFnu_p_host'
        self.x_scale = 'nu_scale'
        self.nuLnu_p_host = self.get_Lum(self.nuFnu_p_host)

        self.interp_func = interp1d(self.nu_template, self.nuFnu_template)



    def set_host_pars(self, fit_model, model_name):
        self.DL = self.cosmo.get_DL_cm(self.z)

        nuFnu_p, nuFnu_p_err = log_to_lin(
            log_val=fit_model.parameters.get(model_name, 'nuFnu_p_host', 'best_fit_val'),
            log_err=fit_model.parameters.get(model_name, 'nuFnu_p_host', 'best_fit_err'))
        err_rel = nuFnu_p_err / nuFnu_p

        self.nuLnu_p_host = self.get_Lum(nuFnu_p)
        self.nuLnu_p_host_err = self.nuLnu_p_host * err_rel

        err_rel = self.nuLnu_p_host_err / self.nuLnu_p_host
        self.nuLnu_p_host = 4 * np.pi * self.DL * self.DL * self.nuFnu_p_host
        self.nuLnu_p_host_err = self.nuLnu_p_host * err_rel

        self.nu_p, self.nu_p_err = log_to_lin(
            log_val=fit_model.parameters.get(model_name, 'nu_scale', 'best_fit_val'),
            log_err=fit_model.parameters.get(model_name, 'nu_scale', 'best_fit_err'))
        err_rel = self.nu_p_err / self.nu_p
        self.nu_p = self.nu_p * self.nu_p
        self.nu_p_err = self.nu_p * err_rel

    def load_teplate_from_file(self, file_path):
        """
        Loads a template from a file
        """
        xy = np.loadtxt(file_path, usecols=(0, 1), unpack=False, dtype=float)
        xy = xy[np.argsort(xy[:, 0])]
        self.nu_template = xy[:, 0]
        self.nuFnu_template = xy[:, 1]