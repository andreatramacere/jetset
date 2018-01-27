"""
Module: model_template
========================

Overview
--------
   
This module  provides an interface handle the template model




Classes relations
---------------------------------------
.. figure::  classes_template_model.png
   :align:   center  


Classes and Inheritance Structure
----------------------------------------------
.. inheritance-diagram:: BlazarSEDFit.template_model


  


Module API
-----------

Summary
---------
.. autosummary::

   Template
   TemplateParameter
    
Module API

"""

'''
Created on 2013 1 27

@author: orion
'''
from cosmo_tools import Cosmo

from data_loader import log_to_lin, lin_to_log

from scipy.interpolate import interp1d

from scipy import loadtxt,log10,array,power,linspace,shape,zeros,ones

import numpy as np

import os

from spectral_shapes import SED

from model_parameters import ModelParameter,ModelParameterArray


__all__=['Template','TemplateParameter']

class TemplateParameter(ModelParameter):
    """
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to Template parameters.
    """
    def __init__(self,template,**keywords):
        
        self.template=template

        self.par_type_list=['nu-scale','nuFnu-scale']
        
        if 'par_type' in keywords.keys() and keywords['par_type'] not in self.par_type_list:
            print "blob_parameter type %s not allowed"%self.par_type
            print "please choose among ",self.par_type_list
            return



        super(TemplateParameter,self).__init__(  **keywords)
        
       
        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val(self.name,val)
            
        
        
    #OVERRIDES Base Method
    def set(self,**keywords):
        super(TemplateParameter,self).set(**keywords )
        
        """
        sets a parameter value checking for physical boundaries 
        """
        
        if 'val' in keywords.keys():
            self.assign_val(self.name,keywords['val']) 
        
    
    def assign_val(self,name,val):
        """
        assigns the paramter valut to the :class:`Template` object
        """
        setattr(self.template,name,val)
        


class Template(object):
    """
    Class to handle spectral templates
    """
    def __init__(self,template_type,z=None,nu_size=100):
        """
        """
        
        self.Templates_dir=os.path.dirname(__file__)+'/Spectral_Templates_Repo'
        
        if z is not None  :
            self.z=z
            self.z_scale=log10(1+z)
        else:
            self.z=0.0
            self.z_scale=0.0
            
        self.allowed_templates=['Disk','host-galaxy','user-defined']
        
        self.nu_size=nu_size
        
        self.name = template_type
        
        self.model_type='template'
            
        self.SED = SED(name=self.name)
    
        self.parameters = ModelParameterArray()
      
        self.cosmo_eval=Cosmo(units='cm')

        self.DL=self.cosmo_eval.DL(self.z)
        
        

        
        
        #BBB TEMPLATE
        if template_type=='BBB':
            self.__scale='log-log'
            
            file_path='%s/QSO_BBB_Template.dat'%self.Templates_dir
            
            self.load_teplate_from_file(file_path)
            self.nu_p_template= 10**self.nu_template[self.nuFnu_template.argmax()]
            #set nu_template according to redshift
            self.nu_template =  self.nu_template-self.z_scale  
            
            #add parameter
            self.nuFnu_p_BBB=1.0
            self.nu_scale=0.0
            self.parameters.add_par(TemplateParameter(self,name='nuFnu_p_BBB',par_type='nuFnu-scale',val=0.0,val_min=-20.0,val_max=20.0,units='erg cm^-2 s^-1'))
            self.parameters.add_par(TemplateParameter(self,name='nu_scale',par_type='nu-scale',val=0.0,val_min=-2.0,val_max=2.0,units='Hz'))
            self.y_scale='nuFnu_p_BBB'          
            self.x_scale='nu_scale'          
            self.nuLnu_p_BBB=self.get_Lum(self.nuFnu_p_BBB)        
            self.T_BBB=self.get_T_BBB() 
            self.interp_func=interp1d(self.nu_template,self.nuFnu_template)
        
       
        elif template_type=='host-galaxy':
            
            self.__scale='log-log'
            
            file_path='%s/HostgalaxyTemplate.dat'%self.Templates_dir
            
            self.load_teplate_from_file(file_path)
            self.nu_p_template= 10**self.nu_template[self.nuFnu_template.argmax()]
            #set nu_template according to redshift
            self.nu_template =  self.nu_template-self.z_scale  
            
            
            #add parameter
            self.nuFnu_p_host=1.0
            self.nu_scale=0.0
            
            self.parameters.add_par(TemplateParameter(self,name='nuFnu_p_host',par_type='nuFnu-scale',val=0.0,val_min=-20.0,val_max=20.0,units='erg cm^-2 s^-1'))
            self.parameters.add_par(TemplateParameter(self,name='nu_scale',par_type='nu-scale',val=0.0,val_min=-2.0,val_max=2.0,units='Hz'))
            
            self.y_scale='nuFnu_p_host'  
            self.x_scale='nu_scale'          
            self.nuLnu_p_host=self.get_Lum(self.nuFnu_p_host)        

            self.interp_func=interp1d(self.nu_template,self.nuFnu_template)
            
        else:
            print"Wrong template type=%s, allowed="%(template_type,self.allowed_templates)
    
    
    
    def PlotModel(self,Plot,clean=False,autoscale=False):
        if Plot is not None  :
            if clean==True:
                Plot.clean_model_lines()
           
            Plot.add_model_plot(self.SED,autoscale=autoscale,line_style='-')
        else:
            print "the plot window is not defined"
     

    def get_T_BBB(self):
        return self.nu_p_template/(1.39*5.879e10)
    
    
    def set_BBB_pars(self,fit_model):
        
        self.DL=self.cosmo_eval.DL(self.z)
        
        
        nuFnu_p,nuFnu_p_err=log_to_lin(log_val=fit_model.parameters.get('nuFnu_p_BBB','best_fit_val'),log_err=fit_model.parameters.get('nuFnu_p_BBB','best_fit_err'))
        err_rel=nuFnu_p_err/nuFnu_p
        self.nuLnu_p_BBB=self.get_Lum(nuFnu_p)   
        self.nuLnu_p_BBB_err=self.nuLnu_p_BBB*err_rel
       
        self.nu_p,self.nu_p_err=log_to_lin(log_val=fit_model.parameters.get('nu_scale','best_fit_val'),log_err=fit_model.parameters.get('nu_scale','best_fit_err'))
        err_rel=self.nu_p_err/ self.nu_p
        self.nu_p=self.nu_p_template*self.nu_p
        self.nu_p_err=self.nu_p*err_rel
    
        self.T_Disk=self.get_T_BBB()
        self.T_Disk_err= self.T_Disk*err_rel
        
        
    def set_host_pars(self,fit_model):
        
        self.DL=self.cosmo_eval.DL(self.z)
        
        
        nuFnu_p,nuFnu_p_err=log_to_lin(log_val=fit_model.parameters.get('nuFnu_p_host','best_fit_val'),log_err=fit_model.parameters.get('nuFnu_p_host','best_fit_err'))
        err_rel=nuFnu_p_err/nuFnu_p
        
        self.nuLnu_p_host=self.get_Lum(nuFnu_p)                
        self.nuLnu_p_host_err=self.nuLnu_p_host*err_rel
        
        err_rel=self.nuLnu_p_host_err/ self.nuLnu_p_host
        self.nuLnu_p_host=4*np.pi* self.DL* self.DL*self.nuFnu_p_host
        self.nuLnu_p_host_err=self.nuLnu_p_host*err_rel
      
        self.nu_p,self.nu_p_err=log_to_lin(log_val=fit_model.parameters.get('nu_scale','best_fit_val'),log_err=fit_model.parameters.get('nu_scale','best_fit_err'))
        err_rel=self.nu_p_err/ self.nu_p
        self.nu_p=self.nu_p*self.nu_p
        self.nu_p_err=self.nu_p*err_rel
    
    
    def set_Lum(self,nuFnu_p):
        self.L_D=self.get_L_D(nuFnu_p)

    
    def get_Lum(self,nuFnu_p):
        return 4*np.pi*self.DL*self.DL*nuFnu_p
        

    def load_teplate_from_file(self,file_path):
        """
        Loads a template from a file
        """
        x,y=loadtxt(file_path, usecols = (0,1), unpack=True, dtype=float)
            
        my_set=zip(x,y)
        my_set.sort(key=lambda tup: tup[0])
    
        self.nu_template,self.nuFnu_template= array(zip(*my_set))

        
   
    def log_func(self,nu_log):
        
        x_shift=getattr(self,self.x_scale)
        y_shift=getattr(self,self.y_scale)
        if shape(nu_log)==():
                nu=array([nu_log])
        
        
        x_log=nu_log-x_shift

        model=ones(x_log.size)*-20.0 
            
        msk = x_log >self.nu_template.min() 
        msk*= x_log <self.nu_template.max()
            
        
        model[msk] = self.interp_func(x_log[msk])+y_shift
        
        return model
    
    
    def eval(self,fill_SED=True,nu=None,get_model=False,loglog=False):
        """
        Evaluates the Template for the current parameters values
        """    
        if loglog==False:
            log_nu=log10(nu)
            
            lin_nu=nu
        else:
            log_nu=nu

            lin_nu=power(10.,log_nu)
        
        
        log_model= self.log_func(log_nu)
        
        model=power(10.,log_model)

        if fill_SED==True:
            
            self.SED.fill(nu=lin_nu, nuFnu=model)
            #print nu.size,nu
        
        if get_model==True:
            if loglog==False:
            
                return model
            else:
                
                return log_model

        else:
            return None
            
  

#-----------------------------------------------
