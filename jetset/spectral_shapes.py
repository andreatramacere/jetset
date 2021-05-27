
__author__ = "Andrea Tramacere"

import numpy as np
from astropy import units
from .frame_converter import  convert_nuFnu_to_nuLnu_src, convert_nu_to_src, convert_nu_src_to_nu_blob, convert_nuLnu_src_to_nuLnu_blob
from .utils import *

__all__=['SED']





class SED(object):
    """
    Class handling the SED 
    """
    def __init__(self,
                 name=None,
                 nu=None,
                 nuFnu=None,
                 nu_residuals=None,
                 residuals=None,
                 nu_src_residuals=None,
                 nuLnu_src_residuals=None,
                 dl=None,
                 z=None,
                 log_log=False,
                 beaming=None):

        if beaming is None:
            beaming =1

        self.beaming=beaming

        self.name=name

        self._nu_units=units.Hz
        self._nuFnu_units=units.erg/units.cm**2/units.s

        self._nu_src_units=units.Hz
        self._nuLnu_src_units = units.erg/units.s
        self._nuLnu_blob_units  = units.erg/units.s


        self.nu=(nu)

        self.nuFnu=(nuFnu)

        self._loglog=log_log

        if z is not None and dl is not None:
            # calling setter do not change
            self.nu_src= (z)

            # calling setter do not change
            self.nuLnu_src = (z,dl)

        self.nu_residuals=nu_residuals
        self.residuals=residuals

        self.nu_src_residuals = nu_src_residuals
        self.nuLnu_src_residuals = nuLnu_src_residuals

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

    @property
    def nu_src(self):
        return self._nu_src

    @nu_src.setter
    def nu_src(self, z):
        #print('->self._nu',self._nu)
        if self._nu is None:
            self._nu_src = self._nu
        else:
            if self._loglog is True:
                self._nu_src =np.log10( convert_nu_to_src(10**(self._nu.value),z,in_frame='obs') )* self._nu_units
            else:
                self._nu_src = convert_nu_to_src(self._nu.value,z,in_frame='obs') * self._nu_units

        #print('->self._nu_src', self._nu_src)

    @property
    def nuLnu_src(self):
        return self._nuLnu

    @nuLnu_src.setter
    def nuLnu_src(self, t):
        z,dl=t
        #print('2')
        #print('->',t,z,dl,self._loglog)
        if self._nuFnu is None:
            self._nuLnu = None
        else:
            if self._loglog is True:
                self._nuLnu =np.log10( convert_nuFnu_to_nuLnu_src(10**(self._nuFnu.value),z,'obs',dl)) * self._nuLnu_src_units
            else:
                self._nuLnu = convert_nuFnu_to_nuLnu_src(self._nuFnu.value, z, 'obs', dl) * self._nuLnu_src_units

    @property
    def nuLnu_blob(self):
        return convert_nuLnu_src_to_nuLnu_blob(self._nuLnu, beaming=self.beaming, in_frame='src')* self._nuLnu_src_units

    @property
    def nu_blob(self):
        return convert_nu_src_to_nu_blob(self._nu_src, beaming=self.beaming, in_frame='src')* self._nu_units

    def get_model_points(self,log_log=False,frame='obs'):

        check_frame(frame)
        if frame == 'obs':
            x, y = self.nu.value, self.nuFnu.value
        elif frame == 'src':
            x, y = self.nu_src.value, self.nuLnu_src.value
        elif frame == 'blob':
            x,y = self.nu_blob.value, self.nuLnu_blob.value
        else:
            unexpected_behaviour()

        if log_log is True:
            x=np.copy(x)
            y=np.copy(y)
            msk = y>0
            x = np.log10(x)
            y[~msk] = -1E10
            y[msk]= np.log10(y[msk])

        return x, y

    def get_residuals(self, log_log=False):

        residuals = self.residuals
        nu_residuals = self.nu_residuals

        if log_log == False:
            return nu_residuals, residuals
        else:
            return nu_residuals, np.log10(residuals)

    def fill(self,nu=None,nuFnu=None,nu_residuals=None,residuals=None,log_log=False):

        self._loglog=log_log
        #if nu is not None:
        self.nu=(nu)
        
        #if nuFnu is not None:
        self.nuFnu=(nuFnu)
         
        if residuals is not None:
            
            self.residuals=residuals
         
        if nu_residuals is not None:
            
            self.nu_residuals=nu_residuals

    def fill_nuLnu(self, nu_src_residuals=None, nuLnu_src_residuals=None, z=None, dl=None):


        if z is not None and dl is not None:

            #calling setter do not change
            self.nu_src = (z)

            # calling setter do not change
            self.nuLnu_src = (z,dl)

        if nuLnu_src_residuals is not None:
            self.nuLnu_src_residuals = nuLnu_src_residuals

        if nu_src_residuals is not None:
            self.nu_src_residuals = nu_src_residuals
            
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

