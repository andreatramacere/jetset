#!/usr/bin/env python
"""
===================================================================
Moudule: spectral_shapes
===================================================================

This module contains all the classes necessary to handle the spectral
shapes, including spectral templates




Classes and Inheritance Structure
-------------------------------------------------------------------

.. inheritance-diagram:: BlazarSEDFit.spectral_shapes
   


Classes relations
----------------------------------------------

.. figure::  classes_spectral_shapes.png    
   :align:   center     


  
.. autosummary::
   SED
   poly_shape
   
    
Module API
-------------------------------------------------------------------

"""
from __future__ import absolute_import

from numpy import  array,zeros,log10

class SED(object):
    """
    Class handling the SED 
    """
    def __init__(self,name=None,nu=None,nuFnu=None,nu_residuals=None,residuals=None):
        self.name=name
        
        self.nu=nu
        self.nuFnu=nuFnu
        
        self.nu_residuals=nu_residuals
        self.residuals=residuals
        
    def get_model_points(self,log_log=False):
    
 
        if log_log==False:
            
            return self.nu,self.nuFnu
        
        else:
        
            msk=self.nuFnu>0
            
            y=log10(self.nuFnu[msk])

            x=log10(self.nu[msk])
        
        return x,y
            
        
    def get_residuals(self,log_log=False): 
        if log_log==False:
            return self.nu_residuals, self.residuals
        else:
            return log10(self.nu_residuals), self.residuals

    def fill(self,nu=None,nuFnu=None,nu_residuals=None,residuals=None):
        
        if nu is not None:
            self.nu=nu
        
        if nuFnu is not None:
            self.nuFnu=nuFnu
         
        if residuals is not None:
            
            self.residuals=residuals
         
        if nu_residuals is not None:
            
            self.nu_residuals=nu_residuals
            
            
            
class poly_shape(object):
    """
    Class for log-log polynomial shapes
    """
    def __init__(self,name=None,nu=None,nuFnu=None):
        self.name=name
        self.nu=nu
        self.nuFnu=nuFnu
    
    def get_model_points(self):
        return self.nu,self.nuFnu
    
        

