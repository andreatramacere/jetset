from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"


from .model_parameters import ModelParameterArray, ModelParameter
from .spectral_shapes import SED
import numpy as np

__all__=['Model']
class Model(object):
    
    
    def __init__(self,name='no-name',nu_size=100):
        
        self.model_type=None
        
        self.name=name

        self.SED = SED(name=self.name)
    
        self.parameters = ModelParameterArray()    
        
        self._scale='lin-lin'
        
        self.nu_size=nu_size
        
        self.nu_min=None
        
        self.nu_max=None
        
        
   
    def eval(self,fill_SED=True,nu=None,get_model=False,loglog=False,plot=None,label=None):
        
        if nu is None:
            print("--->", self.nu_min,self.nu_max,self.nu_size)
            
            x1=np.log10(self.nu_min)

            x2=np.log10(self.nu_max)
            
            lin_nu=np.logspace(x1,x2,self.nu_size)
            
            model=np.zeros(lin_nu.size)
            #print "x1",x1     
            
            
        else:
            

            if np.shape(nu)==():
 
                nu=np.array([nu])
            
            
            
           
            
        if loglog==False:
            log_nu=np.log10(nu)
            
            lin_nu=nu
            
            model=self.lin_func(lin_nu)
        
        else:
            log_nu=nu

            lin_nu=np.power(10.,log_nu)
    
            log_model= self.log_func(log_nu)
            
            model=np.power(10.,log_model)
        
        
        if fill_SED==True:
            #print "base model", self.model_type,self.nu_min,self.nu_max
            #print lin_nu[0]
            self.SED.fill(nu=lin_nu,nuFnu=model)
        
        
        #print model[0]
        
        if plot is not None:
                if label is None:
                    label= self.name
                    
                self.PlotModel(plot, clean=True, label=self.name)
        
        if get_model==True:
            if loglog==False:
            
                return model
            else:
                
                return log_model
        else:
            return None
        
    def PlotModel(self,Plot,clean=False,autoscale=False,label=None):
        #print "PlotModel base model"
        if Plot is not None:
            if clean==True:
                Plot.clean_model_lines()
           
            if label is None:
                    label= self.name
           
            Plot.add_model_plot(self.SED,autoscale=autoscale,line_style='-',label=label)
        
            #print self.SED.nu
            
            
        else:
            print ("the plot window is not defined")
            
    def set_nu_grid(self,nu_min=None,nu_max=None,nu_size=None):
        if nu_size is not None:
            self.nu_size=nu_size
        
        if nu_min is not None:
            self.nu_min=nu_min
        
        if nu_max is not None:
            self.nu_max=nu_max
    
    def lin_func(self,nu):
        pass   
    
    def log_func(self,log_nu):
        pass



    def get_residuals(self, data, log_log=False,filter_UL=True):

        model = self.eval(nu=data['nu_data'], fill_SED=False, get_model=True, loglog=False)

        # print "loglog",loglog
        if filter_UL ==True:
            msk=data['UL']==False
        else:
            msk=np.ones(data.size,dt=True)

        residuals = (data['nuFnu_data'] - model) /  data['dnuFnu_data']

        nu_residuals=data['nu_data']


        if log_log == False:
            return nu_residuals[msk], residuals[msk]
        else:
            return  np.log10(nu_residuals[msk]),  residuals[msk]