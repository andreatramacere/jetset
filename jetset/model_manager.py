from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"


import numpy as np

from  . import minimizer

#import sed_models_dic as Model_dic

from .jet_model  import Jet

#from template_model import Template

#from output import WorkPlace

from .model_parameters import ModelParameterArray

from .spectral_shapes import  SED
   
from .base_model import  Model

from .plot_sedfit import  PlotSED

from .utils import  clean_var_name

from .jet_model import Jet

import copy

import  pickle

__all__=['FitModel']

class ModelComoponentContainer(object):

    def __init__(self):
        pass

class FitModel(Model):
    
    
    """
 
    This class creates an interface to handle  a fit Model. The relevant class members are 
    
        - :class:`.SED` object storing the SED data points
        - `model_pars`, a :class:`.par_array_builder` object returned by the :class:`.SSC_param` constructor 
        - :class:`.jet_builder` object 
    
    Parameters
    
    :param jet: (:class:`.par_array_builder` object)
     
    :param elec_distr_type: (str)  name of the electron distribution model
        
    :param template: (:class:`.spectral_shapes.template` object)
    
    
    Members
    
    :ivar SED: (:class:`.spectral_shapes.SED`) object to store the SED data points
    
    
    :ivar model_pars: (:class:`.par_array_builder` object)
     

  
    """
    
    def __init__(self,elec_distr=None,jet=None,name='no-name',out_dir=None,flag=None,template=None,loglog_poly=None,analytical=None,nu_size=100,cosmo=None,  **keywords):
 
        """
        Constructor
        
        args: 
           
        """
        
       
        super(FitModel,self).__init__(  **keywords)
       
        if  jet is not None and elec_distr is not None:
            #!! warning or error?
            
            print ("you can't provide both elec_distr and jet, only one")
            
            raise RuntimeError
        
          
       


        self.sed_data=None
        
        self.nu_min_fit=1E6
        
        self.nu_max_fit=1E30
        
        self.fitname=None
        
        self.name=name

        self.SED=SED(name=self.name)

        self.nu_min=1E6
        self.nu_max=1E30
        self.nu_size=nu_size

        self.flux_plot_lim=1E-30

        self.components_list=[]

        self.components=ModelComoponentContainer()

        self.parameters=ModelParameterArray()

        if elec_distr is not None:
            jet=Jet(cosmo=cosmo,name=flag, electron_distribution=elec_distr, jet_workplace=None)
            
            self.add_component(jet)
        
        if jet is not None:
            self.add_component(jet)

        if jet is not None:
            self.cosmo = jet.cosmo
        else:
            self.cosmo = cosmo


        if template is not None:
            self.add_component(template)
           
        
        if loglog_poly is not None:
            self.add_component(loglog_poly)
        
        
        if analytical is not None:
            self.add_component(analytical)




    def plot_model(self,plot_obj=None,clean=False,sed_data=None):
        if plot_obj is None:
            plot_obj=PlotSED(sed_data=sed_data)


        if clean==True:
            plot_obj.clean_model_lines()

        line_style='--'


        for mc in self.components_list:
            comp_label = mc.name
            #print ('comp_label',comp_label,mc.SED)
            try:
                plot_obj.add_model_plot(mc.SED, line_style=line_style,label=comp_label,flim=self.flux_plot_lim)
            except Exception as e:
                #print('name',e)
                pass

            if hasattr(mc,'spectral_components_list'):
                for c in mc.spectral_components_list:

                    comp_label = c.name
                    if comp_label!='Sum':
                        #print('comp comp_label', comp_label,c.SED)
                        plot_obj.add_model_plot(c.SED, line_style=line_style, label=comp_label, flim=self.flux_plot_lim)

        line_style = '-'
        #print('comp_label', self.name)
        plot_obj.add_model_plot(self.SED, line_style=line_style, label=self.name, flim=self.flux_plot_lim,fit_range=np.log10([self.nu_min_fit,self.nu_max_fit]))

        plot_obj.add_residual_plot(data=sed_data, model=self,fit_range=np.log10([self.nu_min_fit,self.nu_max_fit]))
        return plot_obj




    def set_nu_grid(self,nu_min=None,nu_max=None,nu_size=None):
        if nu_size is not None:
            self.nu_size=nu_size
        
        if nu_min is not None:
            self.nu_min=nu_min
        
        if nu_max is not None:
            self.nu_max=nu_max
        
        for model_comp in self.components_list:
        
            if nu_size is not None:
                model_comp.nu_size=nu_size
            
            if nu_min is not None:
                model_comp.nu_min=nu_min
            
            if nu_max is not None:
                model_comp.nu_max=nu_max

    def add_component(self, moldel_comp):
        
        self.components_list.append(moldel_comp)
        
        for par in moldel_comp.parameters.par_array:
            
            #if moldel_comp.model_type=='jet':
                
            #    fit_range=set_param_rage(par.name,par.val)
                
            #    par.set(fit_range=fit_range)
                
            self.parameters.add_par(par)

        #print('-->', self.components, moldel_comp.name, moldel_comp)
        setattr(self.components,  clean_var_name(moldel_comp.name), moldel_comp)
    
    def show_pars(self):
        """
        shows all the parameters of the fit_array
        
        """
        
        self.parameters.show_pars()

    def show_model(self):
        for c in self.components_list:
            c.show_model()

    def freeze(self,par_name):
        self.set(par_name,'frozen')


    def free(self,par_name):
        self.set(par_name,'free')
    
    def set(self,par_name,*args,**kw):
        """
        sets the value of a specific parameter of the parameter array
        
        :param par_name: 
        
        .. note::
            see  the documentation of :class:`.par_array_builder` for arguments and keywords
        """
        #print "Model in model manager",args,kw
        
        self.parameters.set(par_name,*args,**kw)

    def set_par(self,par_name,val):
        """
        shortcut to :class:`ModelParametersArray.set` method
        set a parameter value

        :param par_name: (srt), name of the parameter
        :param val: parameter value

        """

        self.parameters.set(par_name, val=val)
  
    def get(self,par_name,*args):
        """
        returns the value of a specific keyword for a specific parameter of the parameter array
        
        :param par_name: 
        
        .. note::
            see  the documentation of :class:`.par_array_builder` for arguments and keywords
        """
        #print "Model in model manager",args,kw
        
        return self.parameters.get(par_name,*args)
        
        
    def get_val(self,par_name):
        """
        returns  the value of a specific parameter of the parameter array
        
        :param par_name: 
        
        .. note::
            see  the documentation of :class:`.par_array_builder` for arguments and keywords
        """
        #print "Model in model manager",args,kw
        
        return self.parameters.get(par_name,'val')
  
  
    def fit(self,sed_data,nu_min,nu_max,fitname=None):
        """
        shortcut to call :func:`.minimizer.fit_SED` 
        
        :param SEDdata: SEDdata object
        :param nu_min: minimun frequency for the fit range interval
        :param nu_max: maximum frequency for the fit range interval
        """
        self.sed_data=sed_data
        
        self.nu_min_fit=nu_min
        
        self.nu_max_fit=nu_max
        
        self.fitname=fitname
        
        return minimizer.fit_SED(self,self.sed_data, self.nu_min_fit, self.nu_max_fit,self.fitname)
  


    
    
    
    def eval(self,nu=None,fill_SED=True,get_model=False,loglog=False,label=None,phys_output=False):
        """
        evaluates the SED for the current parameters and fills the :class:`.SED` member
        """
        
        
        if nu is None:
            #print ("--->", self.nu_min,self.nu_max,self.nu_size)
            
            x1=np.log10(self.nu_min)

            x2=np.log10(self.nu_max)
            
            lin_nu=np.logspace(x1,x2,self.nu_size)

            log_nu=np.log10(lin_nu)
            
            model=np.zeros(lin_nu.size)
                 
          
        else:
            
            if np.shape(nu)==():
 
                nu=np.array([nu])
            
            if loglog==True:
                lin_nu=np.power(10.,nu)
                log_nu=nu
            else:
                log_nu=np.log10(nu)
                lin_nu=nu
        
            
            
        #print lin_nu 
        model=np.zeros(lin_nu.size)
        
        
        
        for model_comp in self.components_list:
            
            #print "model",model_comp.name
            
            if loglog==False:
                model+= model_comp.eval(nu=lin_nu,fill_SED=fill_SED,get_model=True,loglog=loglog)
            else:
                model+= np.power(10.,model_comp.eval(nu=log_nu,fill_SED=fill_SED,get_model=True,loglog=loglog))
         
        
    
        if fill_SED==True:
 
            self.SED.fill(nu=lin_nu,nuFnu=model)
            #TODO ADD cosmo properly to all components
            #self.SED.fill_nuLnu(z=self.jet_obj.get_par_by_type('redshift').val, dl=self.jet_obj.get_DL_cm())
            
        if get_model==True:
            
            if loglog==True:
                model=np.log10(model)
                    
            
            return model        
        
        else:
            
            return None



    @classmethod
    def _build_serializable(cls):
        return cls()

    def save_model(self,file_name):
        c=self._build_serializable()

        c._serialized_model_list=[]
        for _m in self.components_list:
            if isinstance(_m,Jet):
                c._serialized_model_list.append([_m._serialize_model(),_m.model_type])
            else:
                c._serialized_model_list.append([_m,_m.model_type])


        _keep_member_list = ['model_type',
                             'name',
                             # TODO 'parameters' removed  form _keep_member_list because we need to improve parmeter serilization in order to handle _blob member
                             #'parameters'
                             'SED',
                             '_scale',
                             'nu_size',
                             'nu_min',
                             'nu_max',
                             'sed_data',
                             'nu_min_fit',
                             'nu_max_fit',
                             'fitname',
                             'flux_plot_lim',
                             'cosmo']

        for km in _keep_member_list:
            setattr(c,km,getattr(self,km))


        pickle.dump(c, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

        #self.components = _components
        #self.components_list = _components_list

    @classmethod
    def load_model(cls, file_name):
        #c=cls()
        c = pickle.load(open(file_name, "rb"))
        for _m in c._serialized_model_list:

            print(_m)
            if _m[1]=='jet':
                j = Jet()
                j._decode_model( _m[0])

                c.add_component(j)
            else:
                c.add_component(_m[0])
        c.eval()
        return c