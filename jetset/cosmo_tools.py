"""
Module: cosmo_tools
===================================================================

This module contains the class Cosmo, implementing cosmological
calculations following Hogg 1999 (http://arxiv.org/abs/astroph/9905116)



..



Classes and Inheritance Structure
-------------------------------------------------------------------

.. inheritance-diagram:: BlazarSEDFit.cosmo_tools
    


  
.. autosummary::
    Cosmo
   
Classes relations
----------------------------------------------

.. figure::  classes_cosmo_tools.png    
   :align:   center     


    
Module API
-------------------------------------------------------------------

"""

from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"


import numpy as np
from scipy import integrate
#import jet_wrapper
from .jetkernel import jetkernel as BlazarSED

from .output import section_separator


__all__=['Cosmo']


class Cosmo(object):
    """
    
    This class sets-up a cosmological framework to evalute comosmological distances and volume.
    Calculations are done following Hogg 1999 

    :param units: (str), units for distances "cm", or "Gpc", or "Mpc"
    :param Omega_m: (float), if not provided taken from BlazarSED code value
    :param Omega_l: (float), if not provided taken from BlazarSED code value
    :param Omega_k: (float), if not provided taken from BlazarSED code value
        
    """
    def __init__(self, units=None, Omega_m=None, Omega_l=None, Omega_k=None, H_0=None,steps_DC=None,steps_dVC=None):
        if Omega_m is None:
            self._Omega_m = BlazarSED.Omega_matter
        else:
            self._Omega_m = Omega_m
        
        if  Omega_l is None:
            self._Omega_l = BlazarSED.Omega_lambda
        else:
            self._Omega_l = Omega_l
        
        if Omega_k is None:
            self._Omega_k = BlazarSED.Omega_k
        else:
            self._Omega_k = Omega_k
        
        if H_0 is None:
            self._H_0 = BlazarSED.H_0
        else:
            self._H_0=H_0
        
        self._units=units

        if steps_DC is None:
            self._steps_DC=10000
        else:
            self._steps_DC=steps_DC
        
        if steps_dVC is None:
            self._steps_dVC=10000
        else:
            self._steps_dVC=steps_dVC

        if self._units is None or self._units=='cm':
            self._D_H= BlazarSED.vluce_km/self._H_0* BlazarSED.parsec*1E6*1E2
        
        elif self.units=='Mpc':
            print ('Mpc')
            self._D_H= BlazarSED.vluce_km/self._H_0
        
        elif self.units=='Gpc':
            print ("Gpc")
            self._D_H= BlazarSED.vluce_km/self._H_0/1E3
    
    
    
    def show(self):
        """
        shows the set-up of the comological framework

        """
        
        print (section_separator)
        print ("Cosmological set-up")
        print ("units= ",self._units)
        print ("H0 = ",self._H0)
        print ("Omega matter =",self._Omega_m)
        print ("Omega lambda = ",self._Omega_l)
        print ("Omega k = ",self._Omega_k)
        print (section_separator)



    def EZ(self,z):
        z2=(1+z)*(1+z)
        z3=z2*(1+z)
        return np.sqrt(self._Omega_m*z3 + self._Omega_k*z2 + self._Omega_l)
        


    def DC(self,z):
        """
        evaluates the comoving (line-of-sight) distance
        
        :param z: (float) redhsift
        
        :returns: comoving (line-of-sight) distance
        """
        z_grid=np.linspace(0,z,self._steps_DC)
        f=1.0/self.EZ(z_grid)
        integral = integrate.simps(f,z_grid)
        return self._D_H*integral



    def DM(self,z):
        """
        evaluates the comoving (transverse) distance
        
        :param z: (float) redhsift
        
        :returns: comoving (transverse) distance
        """
        return self.DC(z)
    
    
    
    def DA(self,z):
        """
        evaluates the angular distance
        
        :param z: (float) redhsift
        
        :returns: angular distance
        """
        return self.DM(z)/(1+z)
    
    def DL(self,z):
        """
        evaluates the luminosity distance

        :param z: (float) redhsift
    
        :returns: luminosity distance
        
        """
        return self.DM(z)*(1+z)
    
    def four_pi_dVC(self,z):
        a=self.DA(z)
        return np.pi*4.0*self._D_H*a*a*(1+z)*(1+z)/self.EZ(z)

    
    
    
    def VC(self,z2,z1=None):
        """
        evaluates the comoving volume  in a spehrical shell between z1 and  z2
        
        :param z2: (float) redhsift
        :param z1: (float) redhsift, optional, if not provided
                    the returned volume will be the sphere with a z=z2                    

        :returns: comoving volume
        """
        if z1 is not None:
            r2=self.DC(z2)
            r1=self.DC(z1)
            
            a1=r1*r1*r1
            a2=r2*r2*r2
            
            return 4./3*np.pi*(a2-a1)
        else:
            a=self.DC(z2)
            return 4./3*np.pi*a*a*a