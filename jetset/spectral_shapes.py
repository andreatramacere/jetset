from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

from numpy import  array,zeros,log10
from astropy import units


__all__=['SED']





class SED(object):
    """
    Class handling the SED 
    """
    def __init__(self,name=None,nu=None,nuFnu=None,nu_residuals=None,residuals=None):
        self.name=name

        self._nu_units=units.Hz
        self._nuFnu_units=units.erg/units.cm**2/units.s

        #if nu is not None:
        self.nu=nu

        #if nuFnu  is not None:
        self.nuFnu=nuFnu
        
        self.nu_residuals=nu_residuals
        self.residuals=residuals


    @property
    def nu(self):
        return self._nu

    @nu.setter
    def nu(self,nu):
        if nu is None:
            self._nu=nu
        else:
            self._nu=nu*self._nu_units


    @property
    def nuFnu(self):
        return self._nuFnu


    @nuFnu.setter
    def nuFnu(self,nuFnu):
        if nuFnu is None:
            self._nuFnu=nuFnu
        else:
            self._nuFnu = nuFnu * self._nuFnu_units



    def get_model_points(self,log_log=False):
    
 
        if log_log==False:
            
            return self.nu.value,self.nuFnu.value
        
        else:
        
            msk=self.nuFnu.value>0
            
            y=log10(self.nuFnu[msk].value)

            x=log10(self.nu[msk].value)
        
        return x,y

    def get_residuals(self, log_log=False):

        residuals = self.residuals
        nu_residuals = self.nu_residuals

        if log_log == False:
            return nu_residuals, residuals
        else:
            return nu_residuals, log10(residuals)

    def fill(self,nu=None,nuFnu=None,nu_residuals=None,residuals=None):
        
        #if nu is not None:
        self.nu=nu
        
        #if nuFnu is not None:
        self.nuFnu=nuFnu
         
        if residuals is not None:
            
            self.residuals=residuals
         
        if nu_residuals is not None:
            
            self.nu_residuals=nu_residuals
            
            
#
# class poly_shape(object):
#     """
#     Class for log-log polynomial shapes
#     """
#     def __init__(self,name=None,nu=None,nuFnu=None):
#         self.name=name
#         self.nu=nu
#         self.nuFnu=nuFnu
#
#     def get_model_points(self):
#         return self.nu,self.nuFnu
#
#

