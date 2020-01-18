from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)




__author__ = "Andrea Tramacere"





import os
import json
import  sys
import pickle
import six

import numpy as np
from numpy import log10,array,zeros,power,shape

from scipy.interpolate import interp1d
from astropy.table import Table
import warnings
import inspect



on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd == True:
    try:
        from jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from jetkernel import jetkernel as BlazarSED

from . import spectral_shapes

from .model_parameters import ModelParameterArray, ModelParameter
from .base_model import  Model

from .output import makedir,WorkPlace,clean_dir

from .jetkernel_models_dic import nuFnu_obs_dic,gamma_dic,available_N_distr,N_distr_descr,n_seed_dic,allowed_disk_type

from  .plot_sedfit import PlotSED,PlotPdistr,PlotSpecComp,PlotSeedPhotons,plt

from .cosmo_tools import Cosmo

from .utils import *


__all__=['Jet','JetParameter','JetSpecComponent','ElectronDistribution','build_emitting_region_dic',
         'build_ExtFields_dic']



class NoTraceBackWithLineNumber(Exception):
    def __init__(self, msg):
        try:
            ln = sys.exc_info()[-1].tb_lineno
        except AttributeError:
            ln = inspect.currentframe().f_back.f_lineno
        self.args = "{0.__name__} (line {1}): {2}".format(type(self), ln, msg),
        sys.exit(self)

def new_version_warning():
    m = '\n\n' + '*'*80 + '\n'
    m+= 'Something wrong has happened. Please, look at the exception message. Note also that\n'
    m+= 'starting from version 1.1.0, the R parameter as default is linear \nand not logarithmic, please update your scripts\n'
    m+= 'Also the format of the jet_model has changed, now it is a binary file.\n'
    m+= '*' * 80 + '\n'
    warnings.warn(m)


def old_model_warning():
    m = '\n\n' + '*'*80 + '\n'
    m+= 'you are loading a model supported for version<1.1.0, starting from version 1.1.0 \n'
    m+= 'the saved model has changed,  plase update to the new model the new format, \n'
    m += 'by saving it with this version\n'
    m+= '*' * 80 + '\n'
    warnings.warn(m)

class JetkerneltException(NoTraceBackWithLineNumber):

    def __init__(self, message='Jeset  exception', debug_message=''):
        super(JetkerneltException, self).__init__(message)
        self.message=message
        self.debug_message=debug_message


def str_hook(pairs):
    new_pairs = []
    for key, value in pairs:
        if isinstance(value, unicode):
            value = value.encode('utf-8')
        if isinstance(key, unicode):
            key = key.encode('utf-8')
        new_pairs.append((key, value))
    return dict(new_pairs)



def safe_run(func):

    def func_wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
           message =  'the jetkernel failed\n'
           message += '\n exception message: '
           message += '%s'%e
           new_version_warning()

           raise JetkerneltException(message=message)


    return func_wrapper


class JetModelDictionaryPar(object):

    def __init__(self,
                 ptype=None,
                 vmin=None,
                 vmax=None,
                 punit=None,
                 froz=False,
                 log=False,
                 allowed_values=None ):

        self.ptype =ptype
        self.vmin =vmin
        self.vmax =vmax
        self.punit =punit
        self.froz =froz
        self.log =log
        self.allowed_values =allowed_values



class JetParameter(ModelParameter):
    """
    
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to  handles SSC/EC parameters, 
    overriding the :meth:`.ModelParameter.set` in order to propagate the
    parameter value to the BlazarSED object instance
           
 
    """
    def __init__(self,blob,**keywords):
    
        self._blob=blob
        
        self.par_type_list=['electron_density',
                            'Disk',
                            'BLR',
                            'DT',
                            'region_size',
                            'region_position',
                            'electron_energy',
                            'LE_spectral_slope',
                            'HE_spectral_slope',
                            'high-energy-cut-off',
                            'low-energy-cut-off',
                            'spectral_curvature',
                            'turn-over-energy',
                            'magnetic_field',
                            'beaming',
                            'jet-viewing-angle',
                            'jet-bulk-factor',
                            'redshift']
        
        if 'par_type' in keywords.keys() and keywords['par_type'] not in self.par_type_list:
            msg= "blob_parameter type %s not allowed"%keywords['par_type'] + "\n please choose among %s"%self.par_type_list
            raise ValueError("%s"%msg)
            
        
 
        super(JetParameter,self).__init__(  **keywords)
        
       
        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val(self.name,val)
            
        
        
    #OVERRIDES Base Method

    @safe_run
    def set(self,**keywords):
        """        
        overrides the  :meth:`.ModelParameter.set` method in order to propagate the
        parameter value to the BlazarSED object instance
        """
        super(JetParameter,self).set(**keywords )
        
        
        if 'val' in keywords.keys():
            
            self.assign_val(self.name,keywords['val']) 




    def assign_val(self,name,val):
        """
        sets the :class:`.JetParameter` value in the BlazarSED object
        """
        #TODO improve this with allowed values: Done!
        #if name=='disk_type':

        #    if val not in allowed_disk_type:

        #        raise RuntimeError('Disk type %s, not in allowed '%val,allowed_disk_type)
        #print('-> assign val',name,val)
        b=getattr(self._blob,name)
        
        
        #print('1 name',name,val,self._val.islog)
        if type(b)==int:
            if self._val.islog is True:
                val=10**val
            val=int(val)

        elif type(b)==float:
            if self._val._islog is True:
                val=10**val
            val=float(val)

        elif type(b)==str:
            val=val
        #print('2 name',name, val, self._val.islog)
        setattr(self._blob,name,val)
        
        
    


def build_emitting_region_dic(beaming_expr='delta'):
    """
    
    Builds a dictionary to init the :class:`.JetParameter` 
    objects   for the emitting region:
        
        - **R**, the radius of the emitting region in cm
        
        - **B**, the magnetic field in G
        
        - **beaming**, the beaming factor
        
        - **z**, the redshift

   
    """
    
    model_dic={}
    model_dic['R']=JetModelDictionaryPar(ptype='region_size',vmin=1E3,vmax=1E30,punit='cm',froz=False,log=False)
    #    ['region_size',0,30,'cm',False,True]
    model_dic['R_H'] = JetModelDictionaryPar(ptype='region_position',vmin=0,vmax=None,punit='cm',froz=True)
    #['region_position', 0, None, 'cm']
    model_dic['B']=JetModelDictionaryPar(ptype='magnetic_field',vmin=0,vmax=None,punit='G')
    #['magnetic_field',0,None,'G']
    
    if beaming_expr=='bulk_theta':
        model_dic['theta']=JetModelDictionaryPar(ptype='jet-viewing-angle',vmin=0,vmax=None,punit='deg')
        #['jet-viewing-angle',0.0,None,'deg']
        model_dic['BulkFactor']=JetModelDictionaryPar(ptype='jet-bulk-factor',vmin=1.0,vmax=None,punit='Lorentz-factor')
        #['jet-bulk-factor',1.0,None,'Lorentz-factor']
    elif beaming_expr=='delta' or beaming_expr=='':
        model_dic['beam_obj']=JetModelDictionaryPar(ptype='beaming',vmin=1E-4,vmax=None,punit='Lorentz-factor')
        #['beaming', 1, None, '']
    else:
        raise  RuntimeError('''wrong beaming_expr="%s" value, allowed 'delta' or 'bulk_theta' '''%(beaming_expr))
    
    model_dic['z_cosm']=JetModelDictionaryPar(ptype='redshift',vmin=0,vmax=None,punit='')
    # ['redshift',0,None,'']
    
    
    return model_dic
    

def BLR_constraints(L_Disk):
    r1_min = 1E17 * (L_Disk / 1E45) ** 0.5
    r2_min = r1_min * 1.01
    r2_max = r1_min * 2
    print('-> BLR constraints',L_Disk,r1_min, r2_min, r2_max)
    return r1_min, r2_min, r2_max



def DT_constraints(L_Disk):
    return 1E18 * (L_Disk / 1E45) ** 0.5

def build_ExtFields_dic(EC_model_list,allowed_EC_components_list,):
    """

        """
    

        
    
    model_dic={}

    for EC_model in EC_model_list:

        #print('EC_model',EC_model)
        #if EC_model not in allowed_EC_components_list:
        #   raise RuntimeError("EC model %s not allowed"%EC_model,"please choose among ", allowed_EC_components_list)
     
     
        if 'Disk' in EC_model:

            model_dic['L_Disk']=JetModelDictionaryPar(ptype='Disk',vmin=0,vmax=None,punit='erg/s')
            #['Disk',0,None,'erg/s']
            model_dic['R_inner_Sw']=JetModelDictionaryPar(ptype='Disk',vmin=0,vmax=None,punit='Sw. radii')
            # ['Disk',0,None,'Sw. radii']
            model_dic['R_ext_Sw']=JetModelDictionaryPar(ptype='Disk',vmin=0,vmax=None,punit='Sw. radii')
            # ['Disk',0,None,'Sw. radii']
            model_dic['T_Disk']=JetModelDictionaryPar(ptype='Disk',vmin=0,vmax=None,punit='K')
            # ['Disk',0,None,'K']
            model_dic['accr_eff']=JetModelDictionaryPar(ptype='Disk',vmin=0,vmax=None,punit='')
            #['Disk',0,None,'']
            model_dic['disk_type']=JetModelDictionaryPar(ptype='Disk',vmin=None,vmax=None,punit='',froz=True,allowed_values=allowed_disk_type)
            #['Disk',None,None,'',True,False,allowed_disk_type]
            model_dic['M_BH'] = JetModelDictionaryPar(ptype='Disk',vmin=0,vmax=None,punit='M_sun')
            #['Disk', 0, None, 'M_sun']

        if 'BLR' in EC_model:
            #r1_BLR_min, r2_BLR_min, r2_BLR_max = BLR_constraints(L_Disk)
            model_dic['tau_BLR']=JetModelDictionaryPar(ptype='BLR',vmin=0,vmax=1.0,punit='')
            #['BLR',0.0,1.0,'']
            model_dic['R_BLR_in']=JetModelDictionaryPar(ptype='BLR',vmin=0,vmax=None,punit='cm',froz=True)
            #['BLR',0,None,'cm',True]
            model_dic['R_BLR_out']=JetModelDictionaryPar(ptype='BLR',vmin=0,vmax=None,punit='cm',froz=True)
            #['BLR',0,None,'cm',True]
    
    
            
        if 'DT' in EC_model:
            model_dic['T_DT']=JetModelDictionaryPar(ptype='DT',vmin=0,vmax=None,punit='K')
            #['DT',0.0,None,'K']
            model_dic['R_DT']=JetModelDictionaryPar(ptype='DT',vmin=0,vmax=None,punit='cm')
            #['DT',0,None,'cm',True]
            model_dic['tau_DT']=JetModelDictionaryPar(ptype='DT',vmin=0,vmax=1.0,punit='')
            #['DT',0.0,1.0,'']

    return model_dic








class ArrayElectronDistribution(object):

    def __init__(self, e_array, n_array, gamma_grid_size=None):
        self.e_array = e_array
        self.n_array = n_array
        _size = e_array.size

        if n_array.size != _size:
            raise RuntimeError('e_array and n_array must have same size')
        self.size = _size
        if gamma_grid_size is None:
            self.gamma_grid_size = _size * 2 + 1


class ElectronDistribution(object):

    def __init__(self,name,jet,gamma_grid_size=None,log_values=False):

        self.elec_models_list = available_N_distr

        if name == 'from_array':
            pass
        elif name not in self.elec_models_list:
            raise Warning("electron distribution model %s not allowed" % name)
            print("please choose among: ", self.elec_models_list)
            return None
        else:
            pass

        self._name=name
        self._log_values=log_values


        self._jet=jet
        set_str_attr(jet._blob,'DISTR',name)

        if gamma_grid_size is not None:
            self.set_grid_size(gamma_grid_size)

        else:
            self._set_blob()
            self._fill()


    @staticmethod
    def available_distributions():
        for k in N_distr_descr.keys():
            print('%s: %s'%(k,N_distr_descr[k]))


    @classmethod
    def from_array(cls,jet,custom_Ne):
        name='from_array'


        BlazarSED.build_Ne_custom(jet._blob,custom_Ne.size)
        Ne_custom_ptr = getattr(jet._blob, 'Ne_custom')
        gamma_custom_ptr = getattr(jet._blob,'gamma_e_custom')
        for ID in range(custom_Ne.size):
            BlazarSED.set_elec_custom_array(gamma_custom_ptr,jet._blob,custom_Ne.e_array[ID], ID)
            BlazarSED.set_elec_custom_array(Ne_custom_ptr, jet._blob, custom_Ne.n_array[ID],ID)

        setattr(jet._blob,'gmin',custom_Ne.e_array[0] )
        setattr(jet._blob,'gmax',custom_Ne.e_array[-1] )

        return cls(name,jet,custom_Ne.gamma_grid_size)


    def update(self):
        self._set_blob()
        self._fill()

    def set_grid_size(self,gamma_grid_size):
        setattr(self._jet._blob,'gamma_grid_size' ,gamma_grid_size)
        self._set_blob()
        self._fill()

    def _set_blob(self):
        #BlazarSED.MakeNe(self._jet._blob)
        BlazarSED.InitNe(self._jet._blob)
        self._N_name, self._gamma_name = gamma_dic['electron_distr']
        self.Ne_ptr = getattr(self._jet._blob, self._N_name)
        self.gamma_ptr = getattr(self._jet._blob, self._gamma_name)


    def _fill(self):
        size = self._jet._blob.gamma_grid_size
        self.gamma = zeros(size)
        self.n_gamma = zeros(size)

        for ID in range(size):
            self.gamma[ID] = BlazarSED.get_elec_array(self.gamma_ptr, self._jet._blob, ID)
            self.n_gamma[ID] = BlazarSED.get_elec_array(self.Ne_ptr, self._jet._blob, ID)


    def plot(self, p=None, y_min=None,y_max=None,x_min=None,x_max=None):

        self.update()
        if p is None:
            p=PlotPdistr()

        p.plot_distr(self.gamma,self.n_gamma,y_min=y_min,y_max=y_max,x_min=x_min,x_max=x_max)

        return p


    def plot3p(self, p=None,y_min=None,y_max=None,x_min=None,x_max=None):
        self.update()

        if p is None:
            p = PlotPdistr()

        p.plot_distr3p(self.gamma,self.n_gamma,y_min=y_min,y_max=y_max,x_min=x_min,x_max=x_max)

        return p

    def set_bounds(self,a,b,log_val=False):
        if log_val == False:
            return [a,b]

        else:
            return np.log10([a,b])

    def _build_electron_distribution_dic(self,electron_distribution_name):
        """
        Builds the dictionary to init the :class:`.JetParameter`
        objects   for the electron  distribution:

        The following :class:`.JetParameter`: objects the *do not* depend
        on the type of electron distribution

                - N, particle density in cm^-3

                - gmin, the minimum value of the electron Lorentz factor

                - gmax, the maximum value of the electron Lorentz factor

        The following :class:`.JetParameter`: objects *depend* on the type of electron distribution:

            - **power law**, electron_distribution='pl'

               - p

            - **broken power-law**, electron_distribution= **'bkn'**

                - p
                - p_1
                - gamma_break

            - **log-parabola**, electron_distribution= **'lp'**

                - r
                - s
                - gamma0_log_parab (fixed)

            - **log-parabola** with a low energy power-law tail, electron_distribution= **'lppl'**

                - r
                - s
                - gamma0_log_parab

            - **log-parabola** defined by peak energy, electron_distribution= **'lpep'**

                - r
                - s
                - gammap_log_parab,

            - **power-law cut-off**, lectron_distribution= **'plc'**

                - p
                - gamma_cut

        """

        model_dic = {}
        model_dic['N'] = JetModelDictionaryPar(ptype='electron_density',vmin=0,vmax=None,punit='cm^-3')
        #['electron_density', 0, None, 'cm^-3']

        a_h,b_h=self.set_bounds(1,1E15,log_val=self._log_values)
        a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
        a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)
        model_dic['gmin'] = JetModelDictionaryPar(ptype='low-energy-cut-off',vmin=a_l,vmax=b_l,punit='lorentz-factor',log=self._log_values)
        #['low-energy-cut-off', a_l, b_l, 'Lorentz-factor',False,self._log_values]
        model_dic['gmax'] = JetModelDictionaryPar(ptype='high-energy-cut-off',vmin=a_h,vmax=b_h,punit='lorentz-factor',log=self._log_values)
        #['high-energy-cut-off', a_h, b_h, 'Lorentz-factor',False,self._log_values]

        if electron_distribution_name == 'pl':
            model_dic['p'] = JetModelDictionaryPar(ptype='HE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['HE_spectral_slope', -10, 10, '']

        if electron_distribution_name == 'bkn':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['p_1'] = JetModelDictionaryPar(ptype='HE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['HE_spectral_slope', -10, 10, '']
            model_dic['gamma_break'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if electron_distribution_name == 'lp':
            model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
            # ['LE_spectral_slope', -10, 10, '']
            model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature', vmin=-15., vmax=15., punit='')
            # ['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t,punit='lorentz-factor', log=self._log_values,froz=True)

        if electron_distribution_name == 'lppl':
            model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            #['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if electron_distribution_name == 'lpep':
            model_dic['r'] =  JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            #['spectral_curvature', -15, 15, '']
            model_dic['gammap_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if electron_distribution_name == 'plc':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['gamma_cut'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if electron_distribution_name == 'spitkov':
            model_dic['spit_index'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['spit_temp'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['spit_gamma_th'] =JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']

        if electron_distribution_name == 'lppl_pile_up':
            model_dic['s'] =JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['r'] =JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            # ['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['gamma_inj'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']

            model_dic['gamma_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
           # ['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['ratio_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=0.0,vmax=None,punit='')
            #['turn-over-energy', 0, None, '']

            model_dic['alpha_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=0.0,vmax=10.0,punit='')
            #['turn-over-energy', 0.0, 10, '']


        if electron_distribution_name == 'bkn_pile_up':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['p_1'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            # ['HE_spectral_slope', -10, 10, '']
            model_dic['gamma_break'] =  JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']

            model_dic['gamma_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['gamma_pile_up_cut'] =JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']

            model_dic['alpha_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=0.0,vmax=10.0,punit='')
            #['turn-over-energy', 0.0, 10, '']

        if electron_distribution_name == 'from_array':
            model_dic.pop('gmin')
            model_dic.pop('gmax')

        return model_dic




class JetSeedPhotons(object):
    """

    """
    def __init__(self,name,blob_object,var_name=None):
        self.name = name

        self._blob_object = blob_object
        self._n_name, self._nu_name = n_seed_dic[self.name]

        self.n_ptr = getattr(blob_object, self._n_name)

        self.nu_ptr = getattr(blob_object, self._nu_name)
        #self.SED = spectral_shapes.SED(name=self.name)
        if var_name is not None:
            self._var_name=var_name

        self.fill(emiss_lim=self._blob_object.emiss_lim)

    def fill(self,log_log=False,emiss_lim=0):
        self.nu,self.n=self.get_spectral_points(log_log=log_log,emiss_lim=emiss_lim)

    def get_spectral_points(self,log_log=False,emiss_lim=0):

        #try:

        size=self._blob_object.nu_seed_size
        x=zeros(size)
        y=zeros(size)

        for i in range(size):
            x[i]=BlazarSED.get_spectral_array(self.nu_ptr,self._blob_object,i)
            y[i]=BlazarSED.get_spectral_array(self.n_ptr,self._blob_object,i)

            #print("->%e %e"%(x[i],y[i]))

        msk_nan=np.isnan(x)
        msk_nan+=np.isnan(y)
        #print('emiss lim',self.get_emiss_lim())
        x[msk_nan]=0.
        y[msk_nan]=emiss_lim

        msk=y<emiss_lim


        y[msk]=emiss_lim



        if log_log==True:
            msk = y <= 0.
            y[msk] = emiss_lim

            #x=x[msk]
            #    y=y[msk]

            x=log10(x)
            y=log10(y)



        return x,y

        #except:
        #    raise RuntimeError ('model evaluation failed in get_spectral_points')


    def plot(self, y_min=None,y_max=None):
        self.fill(emiss_lim=self._blob_object.emiss_lim)
        p=PlotSeedPhotons()
        p.plot(nu=self.nu,nuFnu=self.n,y_min=y_min,y_max=y_max)

        return p


class JetSpecComponent(object):
    """
    """

    def __repr__(self):
        return str(self.show())

    #def __str__(self):
    #    return str(self.show())


    def __init__(self,jet_obj,name,blob_object,var_name=None,state_dict=None,state=None):
        
        self.name=name
        self.jet_obj=jet_obj

        self._blob_object=blob_object
        self._nuFnu_name, self._nu_name=nuFnu_obs_dic[self.name]
        
        self.nuFnu_ptr=getattr(blob_object,self._nuFnu_name)
        
        self.nu_ptr=getattr(blob_object,self._nu_name)
    
        self.SED=spectral_shapes.SED(name=self.name)
        self.seed_field=None
        #print('->',name,n_seed_dic.keys())
        if name in n_seed_dic.keys():
            self.seed_field=JetSeedPhotons(name,blob_object)


        if var_name is not None:
            self._var_name=var_name

            if state_dict is None:
                self._state_dict = dict()
                self._state_dict['on'] = 1
                self._state_dict['off'] = 0
            else:
                self._state_dict=state_dict
            self.state='on'
        else:
            self._state_dict = {}
            self._var_name=None
            self._state='on'

        if state is not None and self._state_dict != {}:
            self.state=state

    def get_emiss_lim(self,seed=False):
        return self._blob_object.emiss_lim


    def fill_SED(self,log_log=False):

        x,y=self.get_SED_points( log_log=log_log)

        self.SED.fill(nu=x,nuFnu=y,log_log=log_log)
        self.SED.fill_nuLnu(z=self.jet_obj.get_par_by_type('redshift').val,dl=self.jet_obj.get_DL_cm())

        if self.seed_field is not None:
            self.seed_field.fill(log_log=log_log)



    def get_SED_points(self, log_log=False):

        size = self._blob_object.nu_grid_size
        x = zeros(size)
        y = zeros(size)

        for i in range(size):
            x[i] = BlazarSED.get_spectral_array(self.nu_ptr, self._blob_object, i)
            y[i] = BlazarSED.get_spectral_array(self.nuFnu_ptr, self._blob_object, i)


        msk_nan = np.isnan(x)
        msk_nan += np.isnan(y)

        x[msk_nan] = 0.
        y[msk_nan] = self.get_emiss_lim()

        msk = y < self.get_emiss_lim()

        y[msk] = self.get_emiss_lim()

        if log_log == True:
            msk = y <= 0.
            y[msk] = self.get_emiss_lim()

            x = log10(x)
            y = log10(y)

        return x, y




    def update(self):
        size = self._blob.nu_grid_size
        x = zeros(size)
        y = zeros(size)

        for i in range(size):
            x[i] = BlazarSED.get_spectral_array(self.nu_ptr, self._blob, i)
            y[i] = BlazarSED.get_spectral_array(self.nuFnu_ptr, self._blob, i)


    def show(self):
        print('name     :',self.name)
        print('var name :',self._var_name)
        print('state    :',self._state)
        if self._state_dict is not None:
            print('allowed states :',[k for k in self._state_dict.keys()])


    @property
    def state(self,):
        return self._state

    @state.setter
    def state(self, val):
        if self._state_dict!={}:
            if val not in self._state_dict.keys():
                raise RuntimeError('val', val, 'not in allowed', self._state_dict.keys())
            self._state = val
            if self._var_name is not None:
                setattr(self._blob_object, self._var_name, self._state_dict[val])
        else:
            raise Warning('the state of the spectral component',self.name,' can not be changed')

    def get_var_state(self,):
        if self._var_name is not None:
            return  getattr(self._blob_object,self._var_name)
        else:
            return None



    def plot(self, y_min=None,y_max=None):
        p=PlotSpecComp()
        p.plot(nu=self.SED.nu.value,nuFnu=self.SED.nuFnu.value,y_min=y_min,y_max=y_max)

        return p



class TempEvol(object):

    def __init__(self,out_dir='temp_ev',flag='tests',clean_work_dir=True):

        self.build_TempEv()

        self.set_path(out_dir, clean_work_dir=clean_work_dir)

        self.set_flag(flag)

    def build_TempEv(self,duration=1E5,
                     TStart_Acc=0.,
                     TStop_Acc=0.,
                     TStart_Inj=0.,
                     TStop_Inj=0.,
                     T_esc_Coeff=1E60,
                     Esc_Index=0.,
                     Acc_Index=1.,
                     Diff_Index=2.,
                     T_SIZE=5000,
                     NUM_SET=50,
                     Lambda_max_Turb=1E30):




        self._temp_ev = BlazarSED.MakeTempEv()

        self._temp_ev.duration = duration
        self._temp_ev.TStart_Acc = TStart_Acc
        self._temp_ev.TStop_Acc = TStop_Acc
        self._temp_ev.TStart_Inj = TStart_Inj
        self._temp_ev.TStop_Inj = TStop_Inj
        self._temp_ev.T_esc_Coeff = T_esc_Coeff
        self._temp_ev.Esc_Index = Esc_Index
        self._temp_ev.Acc_Index = Acc_Index
        self._temp_ev.Diff_Index = Diff_Index
        self._temp_ev.T_SIZE = T_SIZE
        self._temp_ev.NUM_SET = NUM_SET
        self._temp_ev.Lambda_max_Turb=Lambda_max_Turb



    def set_path(self, path, clean_work_dir=True):
        if path.endswith('/'):
            pass
        else:
            path += '/'

        set_str_attr(self._temp_ev, 'path', path)
        # set_str(self._blob.path,path)
        makedir(path, clean_work_dir=clean_work_dir)



    def set_flag(self,flag):
        self._temp_ev.STEM=flag


    def run(self,jet):
        BlazarSED.Run_temp_evolution(jet._blob, self._temp_ev)


class SpecCompList(object):

    def __init__(self,sc_list):
        self._sc_list=sc_list
        self._table=None


    def show(self):
        for sc in self._sc_list:
            sc.show()

    def __repr__(self):
        return str(self.show())

    def build_table(self, restframe='obs'):

        _names = ['nu']
        _cols=[]

        check_frame(restframe)
        if restframe=='obs':
           _cols.append(self._sc_list[0].SED.nu)
        elif restframe=='src':
            _cols.append(self._sc_list[0].SED.nu_src)
        else:
            unexpetced_behaviour()

        for ID,sc in enumerate(self._sc_list):
            _names.append(sc.name)
            if restframe == 'obs':
                _cols.append(sc.SED.nuFnu)
            else:
                _cols.append(sc.SED.nuLnu_src)



        _meta=dict(src_name=sc.jet_obj.name)
        _meta['redshift']=sc.jet_obj.get_par_by_type('redshift').val
        _meta['restframe']= restframe
        self._table = Table(_cols, names=_names,meta=_meta)

    @property
    def table(self):
        return self._table



class Jet(Model):
    """

    This class allows to build a ``Jet`` model providing the interface to the
    C code, giving  full access to the physical parameters and
    providing the methods to run the code.
    A :class:`Jet` object  will store the
    the physical parameters in  the ::py:attr:`Jet.parameters`  that is :class:`.ModelParameterArray` class,
    i.e. a collection of :class:`JetParameter` objects.
    All the physical parameters are  also accessible as attributes of
    the  ::py:attr:`Jet.parameters`

    **Examples**




    """
    def __repr__(self):
        return str(self.show_model())

    #def __str__(self):
    #    return str(self.show_model())

    def __init__(self,
                 cosmo=None,
                 name='tests',
                 electron_distribution='pl',
                 electron_distribution_log_values=False,
                 beaming_expr='delta',
                 jet_workplace=None,
                 verbose=None,
                 clean_work_dir=True,
                 **keywords):

        """

        Parameters
        ----------
        cosmo
        name
        electron_distribution
        electron_distribution_log_values
        beaming_expr
        jet_workplace
        verbose
        clean_work_dir
        keywords
        """
        super(Jet,self).__init__(  **keywords)

        #print('cosmo',cosmo)
        if cosmo is not None:
            self.cosmo=cosmo
        else:
            self.cosmo= Cosmo()
        #print('cosmo', self.cosmo)
        self.name = name

        self.model_type='jet'

        self._scale='lin-lin'

        self._blob = self.build_blob(verbose=verbose)

        if jet_workplace is None:
            jet_workplace=WorkPlace()
            out_dir= jet_workplace.out_dir + '/' + self.name + '_jet_prod/'
        else:
            out_dir=jet_workplace.out_dir+'/'+self.name+'_jet_prod/'

        self.set_path(out_dir,clean_work_dir=clean_work_dir)

        self.set_flag(self.name)


        self._allowed_EC_components_list=['EC_BLR',
                                          'DT',
                                          'EC_DT',
                                          'Star',
                                          'EC_Start',
                                          'CMB',
                                          'EC_CMB',
                                          'Disk',
                                          'EC_Disk',
                                          'All']

        self.EC_components_list =[]
        self.spectral_components_list=[]


        self.spectral_components=SpecCompList(self.spectral_components_list)
        self.add_basic_components()


        self.SED=self.get_spectral_component_by_name('Sum').SED

        self.parameters = ModelParameterArray()

        self._emitting_region_dic = None
        self._electron_distribution_dic= None
        self._external_photon_fields_dic= None




        self.set_electron_distribution(electron_distribution,electron_distribution_log_values)

        self._blob.adaptive_e_binning=0

        self.set_emitting_region(beaming_expr)


        self.flux_plot_lim=1E-30
        self.set_emiss_lim(1E-120)

        self._IC_states = {}
        self._IC_states['on'] = 1
        self._IC_states['off'] = 0

        self._external_field_transf = {}
        self._external_field_transf['blob'] = 0
        self._external_field_transf['disk'] = 1


    @staticmethod
    def available_electron_distributions():
        ElectronDistribution.available_distributions()

    def build_blob(self,verbose=None):

        blob = BlazarSED.MakeBlob()
        # temp_ev=BlazarSED.MakeTempEv()
        #print('blob',blob.griglia_gamma_Ne_log)
        if verbose is None:
            blob.verbose = 0
        else:
            blob.verbose = verbose

        set_str_attr(blob, 'path', './')
        # set_str(blob.path,'./')

        set_str_attr(blob, 'MODE', 'custom')
        # blob.MODE='custom'
        blob.gamma_grid_size = 1000

        blob.nu_IC_size = 50
        blob.nu_seed_size = 100

        blob.do_Sync = 2

        blob.do_SSC = 1

        blob.R = 5.0e15

        blob.B = 0.1

        blob.z_cosm = 0.1

        blob.BulkFactor = 10
        blob.theta = 0.1

        blob.N = 100

        blob.L_Disk=1E45

        blob.L_DT=1E45

        blob.gmin = 2

        blob.gmax = 1e6

        blob.nu_start_Sync = 1e6
        blob.nu_stop_Sync = 1e20

        blob.nu_start_SSC = 1e14
        blob.nu_stop_SSC = 1e30

        blob.nu_start_grid = 1e6
        blob.nu_stop_grid = 1e30

        return blob

    def _serialize_model(self):
        _model = {}
        _model['electron_distribution'] = self._electron_distribution_name
        _model['electron_distribution_log_values'] = self._electron_distribution_log_values
        _model['beaming_expr'] = self._beaming_expr
        _model['spectral_components_name'] = self.get_spectral_component_names_list()
        _model['spectral_components_state'] = [c.state for c in self.spectral_components_list]
        _model['EC_components_name'] = self.EC_components_list
        _model['basic_components_name'] = self.basic_components_list
        _model['cosmo'] = self.cosmo

        _model['pars'] = {}
        _model['pars']=self.parameters._serialize_pars()

        _model['internal_pars'] = {}
        _model['internal_pars']['nu_size'] = self.nu_size
        _model['internal_pars']['nu_seed_size'] = self.nu_seed_size
        _model['internal_pars']['gamma_grid_size'] = self.gamma_grid_size
        _model['internal_pars']['IC_nu_size']=self.IC_nu_size
        _model['internal_pars']['nu_min']=self.nu_min
        _model['internal_pars']['nu_max']=self.nu_max
        _model['internal_pars']['Norm_distr'] = self.Norm_distr
        return _model

    def save_model(self,file_name):
        pickle.dump(self._serialize_model(), open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    def _decode_model(self,_model):

        self.cosmo = _model['cosmo']
        self.model_type = 'jet'

        self.set_blob()
        self.parameters = ModelParameterArray()

        self.set_electron_distribution(str(_model['electron_distribution']))

        for c in self.basic_components_list:
            if c not in _model['basic_components_name']:
                self.del_spectral_component(c)

        self.add_EC_component(_model['EC_components_name'])

        for ID, c in enumerate(_model['spectral_components_name']):
            comp = getattr(self.spectral_components, c)
            if comp._state_dict != {}:
                comp.state = _model['spectral_components_state'][ID]

        self.SED = self.get_spectral_component_by_name('Sum').SED

        self.set_emitting_region(str(_model['beaming_expr']))
        self.set_electron_distribution(str(_model['electron_distribution']))

        _par_dict = _model['pars']
        self.parameters._decode_pars(_par_dict)

        _par_dict = _model['internal_pars']
        for k in _par_dict.keys():
            #print ('set', k,_par_dict[k])
            setattr(self,k,_par_dict[str(k)])
            #self.set_par(par_name=str(k), val=_par_dict[str(k)])

    @classmethod
    def load_old_model(cls, file_name):

        old_model_warning()
        jet = cls()
        with open(file_name, 'r') as infile:
            _model = json.load(infile)

        # print ('_model',_model)

        jet.model_type = 'jet'

        jet.init_BlazarSED()
        jet.parameters = ModelParameterArray()

        jet.set_electron_distribution(str(_model['electron_distribution']))

        for c in jet.basic_components_list:
            if c not in _model['basic_components_name']:
                jet.del_spectral_component(c)

        jet.add_EC_component(_model['EC_components_name'])

        for ID, c in enumerate(_model['spectral_components_name']):
            comp = getattr(jet.spectral_components, c)
            if comp._state_dict != {}:
                comp.state = _model['spectral_components_state'][ID]

        jet.SED = jet.get_spectral_component_by_name('Sum').SED

        jet.set_emitting_region(str(_model['beaming_expr']))
        jet.set_electron_distribution(str(_model['electron_distribution']))
        _par_dict = _model['pars']
        jet.show_pars()
        for k in _par_dict.keys():
            # print ('set', k,_par_dict[k])
            jet.set_par(par_name=str(k), val=_par_dict[str(k)])

        jet.eval()
        return jet

    @classmethod
    @safe_run
    def load_model(cls,file_name):



        _model=pickle.load( open(file_name, "rb" ) )
        jet = cls()
        jet._decode_model(_model)
        jet.show_pars()
        jet.eval()
        return jet

    def set_electron_distribution(self,name=None,log_values=False):

        self._electron_distribution_log_values=log_values

        if self._electron_distribution_dic is not None:
            self.del_par_from_dic(self._electron_distribution_dic)

        if isinstance(name, ArrayElectronDistribution):
            self._electron_distribution_name='from_array'
            #print('ciccio')
            self.electron_distribution=ElectronDistribution.from_array(self,name)

        elif isinstance(name, ElectronDistribution):
            self.electron_distribution = name

        else:



            self.electron_distribution = ElectronDistribution(name, self,log_values=log_values)

        if self.electron_distribution is not None:
            self._electron_distribution_name = name
            self._electron_distribution_dic = self.electron_distribution._build_electron_distribution_dic(
                self._electron_distribution_name)

            self.add_par_from_dic(self._electron_distribution_dic)

        else:
            raise RuntimeError('name for electron distribution was not valid')




    def show_spectral_components(self):

        print ("Spectral components for Jet model:%s"%(self.name))

        for comp in self.spectral_components_list:

            print ("comp: %s "%(comp.name))

        print()




    def get_spectral_component_by_name(self,name,verbose=True):
        for i in range(len(self.spectral_components_list)):
            if self.spectral_components_list[i].name==name:
                return self.spectral_components_list[i]
        else:
            if verbose==True:
                print ("no spectral components with name %s found"%name)
                print ("names in array are:")
                self.show_spectral_components()

            return None



    def list_spectral_components(self):
        for i in range(len(self.spectral_components_list)):
            print (self.spectral_components_list[i].name)

    def get_spectral_component_names_list(self):
        _l=[]
        for i in range(len(self.spectral_components_list)):
            _l.append(self.spectral_components_list[i].name)
        return _l



    def del_spectral_component(self,name):
        print('deleting spectral component', name)
        if name in self.EC_components_list:
            self.del_EC_component(name)
        else:
            self._del_spectral_component(name)
        if name in self.basic_components_list:
            self.basic_components_list.remove(name)

    def _del_spectral_component(self, name, verbose=True):



        comp=self.get_spectral_component_by_name(name,verbose=verbose)

        if comp._state_dict!={}:

            comp.state='off'

        if comp is not None:

            self.spectral_components_list.remove(comp)

            self.spectral_components.__delattr__(name)



    def _add_spectral_component(self, name, var_name=None, state_dict=None,state=None):

        self.spectral_components_list.append(JetSpecComponent(self,name, self._blob, var_name=var_name, state_dict=state_dict,state=state))
        setattr(self.spectral_components,name,self.spectral_components_list[-1])

    def _update_spectral_components(self):
        _l=[]
        for ID,s in enumerate(self.spectral_components_list):
            self.spectral_components_list[ID]=JetSpecComponent(self,s.name, self._blob, var_name=s._var_name, state_dict=s._state_dict,state=s.state)
            setattr(self.spectral_components, self.spectral_components_list[ID].name, self.spectral_components_list[ID])




    def set_emitting_region(self,beaming_expr):

        if  self._emitting_region_dic is not None:
            self.del_par_from_dic(self._emitting_region_dic)

        self._beaming_expr=beaming_expr

        set_str_attr(self._blob,'BEAMING_EXPR',beaming_expr)

        #setattr(self._blob,'BEAMING_EXPR',beaming_expr)

        self._emitting_region_dic=build_emitting_region_dic(beaming_expr=beaming_expr)

        self.add_par_from_dic(self._emitting_region_dic)


    def set_N_from_Ue(self,U_e):
        N = 1.0
        setattr(self._blob, 'N', N)
        gamma_grid_size = self._blob.gamma_grid_size
        self.electron_distribution.set_grid_size(100)
        self.set_blob()
        BlazarSED.EvalU_e(self._blob)
        ratio = self._blob.U_e/ U_e
        self.electron_distribution.set_grid_size(gamma_grid_size)
        self.set_par('N', val=N / ratio)

    def set_N_from_Le(self,L_e):
        gamma_grid_size = self._blob.gamma_grid_size
        self.electron_distribution.set_grid_size(100)
        self.set_blob()
        U_e=L_e/    self._blob.Vol_sphere
        self.electron_distribution.set_grid_size(gamma_grid_size)
        self.set_N_from_Ue(U_e)


    def set_N_from_L_sync(self,L_sync):

        N = 1.0
        setattr(self._blob, 'N', N)
        gamma_grid_size = self._blob.gamma_grid_size
        self.electron_distribution.set_grid_size(100)
        self.set_blob()
        delta = self._blob.beam_obj
        ratio = (BlazarSED.Power_Sync_Electron(self._blob)* delta ** 4)/L_sync
        self.electron_distribution.set_grid_size(gamma_grid_size)

        # print 'N',N/ratio
        self.set_par('N', val=N / ratio)


    def set_N_from_F_sync(self, F_sync):
        self.set_blob()
        DL = self._blob.dist
        L = F_sync * DL * DL * 4.0 * np.pi
        self.set_N_from_L_sync(L)

    def set_N_from_nuLnu(self,nuLnu_src, nu_src):
        """
        sets the normalization of N to match the rest frame luminosity L_0, at a given frequency nu_0
        """
        N = 1.0
        setattr(self._blob, 'N', N)
        gamma_grid_size = self._blob.gamma_grid_size
        self.electron_distribution.set_grid_size(100)



        self.set_blob()

        delta = self._blob.beam_obj
        nu_blob = nu_src / delta

        L_out = BlazarSED.Lum_Sync_at_nu(self._blob, nu_blob) * delta ** 4


        ratio = (L_out / nuLnu_src)
        self.electron_distribution.set_grid_size(gamma_grid_size)

        #print 'N',N/ratio
        self.set_par('N', val=N/ratio)


    def set_N_from_nuFnu(self, nuFnu_obs, nu_obs):
        """
        sets the normalization of N to match the observed flux nu0F_nu0 at a given frequency nu_0
        """

        self.set_blob()
        DL = self._blob.dist

        L = nuFnu_obs * DL * DL * 4.0 * np.pi

        nu_rest = nu_obs * (1 + self._blob.z_cosm)

        self.set_N_from_nuLnu( L, nu_rest)



    def set_B_eq(self, nuFnu_obs, nu_obs, B_min=1E-9,B_max=1.0,N_pts=20,plot=False):
        """
        returns equipartiont B
        """
        # print B_min


        b_grid = np.logspace(np.log10(B_min), np.log10(B_max), N_pts)
        print ('B grid min ',B_min)
        print ('B grid max ',B_max)
        print ('grid points',N_pts)
        U_e = np.zeros(N_pts)
        U_B = np.zeros(N_pts)
        N = np.zeros(N_pts)
        self.set_par('B', b_grid[0])
        self.set_blob()

        for ID, b in enumerate(b_grid):
            self.set_par('B', b)
            self.set_par('N', 1.0)
            # print 'B_eq',ID
            self.set_N_from_nuFnu(nuFnu_obs, nu_obs)
            N[ID]=self.get_par_by_name('N').val
            self.set_blob()
            #
            U_e[ID] = self._blob.U_e
            U_B[ID] = self._blob.UB
            # delta=Jet.get_beaming()
            # print "check L_in=%4.4e L_out=%4.4e"%(L_0,(L_0/delta**4)/BlazarSED.Power_Sync_Electron(Jet._Jet__blob))

        ID_min = np.argmin(U_B + U_e)

        if plot==True:
            #import  pylab as plt
            fig, ax = plt.subplots()
            ax.plot(b_grid,U_e)
            ax.plot(b_grid,U_B)
            ax.plot(b_grid,U_B+U_e)
            ax.scatter(b_grid[ID_min],U_e[ID_min]+U_B[ID_min])

            ax.semilogy()
            ax.semilogx()
            plt.show()


        self.set_par('B', val=b_grid[ID_min])
        self.set_par('N', val=N[ID_min])
        print ('setting B to ',b_grid[ID_min])
        print ('setting N to ',N[ID_min])
        return b_grid[ID_min],b_grid,U_B,U_e



    def add_basic_components(self):
        self.basic_components_list=['Sum','Sync','SSC']

        self._add_spectral_component('Sum')
        self._add_spectral_component('Sync', var_name='do_Sync', state_dict=dict((('on', 1), ('off', 0), ('self-abs', 2))),state='self-abs')

        self._add_spectral_component('SSC', var_name='do_SSC', state_dict=dict((('on', 1), ('off', 0))))

    def add_sync_component(self,state='self-abs'):
        self._add_spectral_component('Sync', var_name='do_Sync',
                                     state_dict=dict((('on', 1), ('off', 0), ('self-abs', 2))), state=state)

    def add_SSC_component(self,state='on'):
        self._add_spectral_component('SSC', var_name='do_SSC', state_dict=dict((('on', 1), ('off', 0))),state=state)

    def del_EC_component(self,EC_components_list):
        if isinstance(EC_components_list, six.string_types):
            EC_components_list = [EC_components_list]

        #print(EC_components_list)

        if 'All' in EC_components_list:
            EC_components_list=self._allowed_EC_components_list[::]


        for EC_component in EC_components_list:


            if EC_component not in self._allowed_EC_components_list:
                raise RuntimeError("EC_component %s not allowed" % EC_component, "please choose among ",
                                   self._allowed_EC_components_list)

            if EC_component=='Disk':
                if self.get_spectral_component_by_name('Disk', verbose=False) is not None:
                    self._del_spectral_component('Disk', verbose=False)
                    self.EC_components_list.remove('Disk')

                if self.get_spectral_component_by_name('EC_Disk',verbose=False) is not None:
                    self._del_spectral_component('EC_Disk')
                    self._blob.do_EC_Disk = 0
                    self.EC_components_list.remove('EC_Disk')

                if self.get_spectral_component_by_name('EC_BLR', verbose=False) is not None:
                    self._blob.do_EC_BLR=0
                    self._del_spectral_component('EC_BLR', verbose=False)
                    self.EC_components_list.remove('EC_BLR')


            if EC_component=='EC_Disk':
                if self.get_spectral_component_by_name('EC_Disk', verbose=False) is not None:
                    self._blob.do_EC_Disk=0
                    self._del_spectral_component('EC_Disk', verbose=False)
                    self.EC_components_list.remove('EC_Disk')


            if EC_component=='EC_BLR':
                if self.get_spectral_component_by_name('EC_BLR', verbose=False) is not None:
                    self._blob.do_EC_BLR=0
                    self._del_spectral_component('EC_BLR', verbose=False)
                    self.EC_components_list.remove('EC_BLR')

            if EC_component=='DT':
                if self.get_spectral_component_by_name('DT', verbose=False) is not None:
                    self._del_spectral_component('DT', verbose=False)
                    self.EC_components_list.remove('DT')
                if self.get_spectral_component_by_name('EC_DT', verbose=False) is not None:
                    self._blob.do_EC_DT = 0
                    self._del_spectral_component('EC_DT', verbose=False)
                    self.EC_components_list.remove('EC_DT')

            if EC_component=='EC_DT':
                if self.get_spectral_component_by_name('EC_DT', verbose=False) is not None:
                    self._blob.do_EC_DT=0
                    self._del_spectral_component('EC_DT', verbose=False)
                    self.EC_components_list.remove('EC_DT')

            if EC_component=='EC_CMB':
                if self.get_spectral_component_by_name('EC_CMB', verbose=False) is not None:
                    self._blob.do_EC_CMB=0
                    self._del_spectral_component('EC_CMB', verbose=False)
                    self.EC_components_list.remove('EC_CMB')

            #if EC_component=='EC_CMB_stat':
            #    if self.get_spectral_component_by_name('EC_CMB_stat', verbose=False) is not None:
            #        self._blob.do_EC_CMB_stat=0
            #        self._del_spectral_component('EC_CMB_stat', verbose=False)
            #        self.EC_components_list.remove('EC_CMB_stat')

        self.del_par_from_dic(build_ExtFields_dic(EC_components_list,self._allowed_EC_components_list))






    def add_EC_component(self,EC_components_list):

        #self._blob.R_BLR_in, self._blob.R_BLR_out, _=BLR_constraints(self._blob.L_Disk)
        #self._blob.R_BLR_DT=DT_constraints(self._blob.L_Disk)

        if isinstance(EC_components_list, six.string_types):
            EC_components_list=[EC_components_list]

        if 'All' in EC_components_list:
            EC_components_list=self._allowed_EC_components_list[::]

        for EC_component in EC_components_list:
            if EC_component not in self._allowed_EC_components_list:
                raise RuntimeError("EC_component %s not allowed" % EC_component, "please choose among ",
                                   self._allowed_EC_components_list)

            if EC_component == 'Disk':
                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_Disk':
                #self._blob.do_EC_Disk=1
                if self.get_spectral_component_by_name('EC_Disk',verbose=False) is None:
                    self._add_spectral_component('EC_Disk', var_name='do_EC_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_Disk')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_BLR':
                #self._blob.do_EC_BLR=1
                if self.get_spectral_component_by_name('EC_BLR',verbose=False) is None:
                    self._add_spectral_component('EC_BLR', var_name='do_EC_BLR', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_BLR')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    # TODO add state
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')


            if EC_component == 'DT':
                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self._add_spectral_component('DT',var_name='do_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('DT')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    # TODO add state
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')


            if EC_component=='EC_DT':
                #self._blob.do_EC_DT=1
                if self.get_spectral_component_by_name('EC_DT',verbose=False) is None:
                    self._add_spectral_component('EC_DT', var_name='do_EC_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_DT')

                #TODO add state
                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self._add_spectral_component('DT',var_name='do_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('DT')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    # TODO add state
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_CMB':
                #self._blob.do_EC_CMB=1
                if self.get_spectral_component_by_name('EC_CMB',verbose=False) is None:
                    self._add_spectral_component('EC_CMB', var_name='do_EC_CMB', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_CMB')

            #if EC_component == 'EC_CMB_stat':
            #    #self._blob.do_EC_CMB_stat = 1
            #    if self.get_spectral_component_by_name('EC_CMB_stat', verbose=False) is  None:
            #        self._add_spectral_component('EC_CMB_stat', var_name='do_EC_CMB_stat', state_dict=dict((('on', 1), ('off', 0))))
            #        self.EC_components_list.append('EC_CMB_stat')



        self.add_par_from_dic(build_ExtFields_dic(self.EC_components_list,self._allowed_EC_components_list))




    def del_par_from_dic(self,model_dic):
        """
        """
        for key in model_dic.keys():

            par=self.parameters.get_par_by_name(key)

            if par is not None:
                self.parameters.del_par(par)




    def add_par_from_dic(self,model_dic):
        """
        add the :class:`.JetParameter` object to the :class:`.ModelParameterArray`
        usign the dictionaries built by the :func:`build_emitting_region_dic`
        and :func:`build_electron_distribution_dic`
        """
        for key in model_dic.keys():

            pname=key
            #log=False
            #allowed_values=None
            p_test=self.parameters.get_par_by_name(pname)

            #TODO add keys to dict and remove length checking

            #print('-> add_par_from_dic', key)
            if p_test is None:

                #print('-> par name', pname, 'log',model_dic[key].log)
                pval=getattr(self._blob,pname)

                #if len(model_dic[key])>5:
                log=model_dic[key].log

                #if len(model_dic[key])>6:
                allowed_values=model_dic[key].allowed_values

                #IMPORTANT
                #This has to stay here, because the parameter, even if log, is stored as linear
                if log is True:
                    #print('-> is log',log)
                    pval = np.log10(pval)
                    #print('-> log_pval', pval)

                ptype=model_dic[key].ptype
                vmin=model_dic[key].vmin
                vmax=model_dic[key].vmax
                punit=model_dic[key].punit

                #froz=False

                #if len(model_dic[key])>4:
                froz=model_dic[key].froz


                self.parameters.add_par(JetParameter(self._blob,name=pname,
                                                     par_type=ptype,
                                                     val=pval,
                                                     val_min=vmin,
                                                     val_max=vmax,
                                                     units=punit,
                                                     frozen=froz,
                                                     log=log,
                                                     allowed_values=allowed_values))



    #PROTECTED MEMBER ACCESS
    #!! controlla i paramteri ....
    def get_electron_distribution_name(self):
        return self.electron_distribution._name



    def get_DL_cm(self,eval=False):

        if eval is True:
            self.set_blob()

        return self._blob.dist


    @safe_run
    def get_beaming(self,theta=None,bulk_factor=None,beaming=None):

        if theta is None:
            theta=self._blob.theta

        if bulk_factor is None:
            bulk_factor=self._blob.BulkFactor

        if beaming is None:
            beaming=self._blob.beam_obj

        BlazarSED.SetBeaming(self._blob)

        return self._blob.beam_obj

    def set_flag(self,flag):
        self._blob.STEM=flag

    def get_flag(self):
        return self._blob.STEM

    def get_path(self):
        return self._blob.path

    def set_path(self,path,clean_work_dir=True):
        if path.endswith('/'):
            pass
        else:
            path+='/'

        set_str_attr(self._blob,'path',path)
        #set_str(self._blob.path,path)
        makedir(path,clean_work_dir=clean_work_dir)

    #def set_SSC_mode(self,val):
    #    self._blob.do_SSC=val

    #def get_SSC_mode(self):
    #    return self._blob.do_SSC





    def set_IC_mode(self,val):

        if val not in self._IC_states.keys():
            raise RuntimeError('val',val,'not in allowed values',self._IC_states.keys())

        self._blob.do_IC=self._IC_states[val]

    def get_IC_mode(self):
        return dict(map(reversed, self._IC_states.items()))[self._blob.do_IC]



    def set_external_field_transf(self,val):
        if val not in self._external_field_transf.keys():
            raise RuntimeError('val',val,'not in allowed values',self._external_field_transf.keys())

        self._blob.EC_stat=self._external_field_transf[val]

    def get_external_field_transf(self):
        return dict(map(reversed, self._external_field_transf.items()))[self._blob.EC_stat]

    def set_emiss_lim(self,val):
        self._blob.emiss_lim=val

    def get_emiss_lim(self):
        return self._blob.emiss_lim


    #def set_disk_type(self,disk_type='MultiBB'):

    #    if disk_type in allowed_disk_type:

    #        self._blob.disk_type = disk_type
    #    else:
    #        raise RuntimeError('Disk type %s, not in allowed '%disk_type,allowed_disk_type)

        #BlazarSED.build_photons(self._blob)

    @property
    def IC_nu_size(self):
        return self._blob.nu_IC_size

    @IC_nu_size.setter
    def IC_nu_size(self, val):
        self.set_IC_nu_size(val)

    def get_IC_nu_size(self):
        return self._blob.nu_IC_size

    def set_IC_nu_size(self, val):
        if val > 1000:
            raise RuntimeError('value can not exceed 1000')
        self._blob.nu_IC_size = val

    @property
    def nu_seed_size(self):
        return self._blob.nu_seed_size

    def get_seed_nu_size(self):
        return self._blob.nu_seed_size

    @nu_seed_size.setter
    def nu_seed_size(self,val):
        self.set_seed_nu_size(val)

    def set_seed_nu_size(self,val):
        if val>1000:
            raise RuntimeError('value can not exceed 1000')
        self._blob.nu_seed_size=val





    def set_gamma_grid_size(self,val):
        self.electron_distribution.set_grid_size(gamma_grid_size=val)

    @property
    def gamma_grid_size(self):
        return self._blob.gamma_grid_size

    @gamma_grid_size.setter
    def gamma_grid_size(self,val):
        self.electron_distribution.set_grid_size(gamma_grid_size=val)

    @property
    def nu_min(self):
        return self._get_nu_min_grid()

    @nu_min.setter
    def nu_min(self, val):
        if hasattr(self, '_blob'):
            self._set_nu_min_grid(val)

    def _get_nu_min_grid(self):
        return  self._blob.nu_start_grid


    def _set_nu_min_grid(self, val):
        self._blob.nu_start_grid=val
    @property
    def Norm_distr(self):
        return self._blob.Norm_distr

    @Norm_distr.setter
    def Norm_distr(self, val):
        if val == 1:
            self._blob.Norm_distr = val
        elif val == 0:
            self._blob.Norm_distr = val
        else:
            raise RuntimeError('value', val, 'not allowed, allowed 0 or 1')

    def switch_Norm_distr_ON(self):
        self._blob.Norm_distr=1

    def switch_Norm_distr_OFF(self):
        self._blob.Norm_distr=0

    @property
    def nu_max(self):
        return self._get_nu_max_grid()

    @nu_max.setter
    def nu_max(self, val):
        if hasattr(self, '_blob'):
            self._set_nu_max_grid(val)

    def _set_nu_max_grid(self, val):
        self._blob.nu_stop_grid=val

    def _get_nu_max_grid(self):
        return  self._blob.nu_stop_grid

    def set_nu_grid_size(self,size):
        self.nu_size=size

    @property
    def nu_size(self):
        return self._get_nu_grid_size()

    @nu_size.setter
    def nu_size(self, val):
        if hasattr(self, '_blob'):
            self._set_nu_grid_size(val)

    def _set_nu_grid_size(self, val):
        self._blob.nu_grid_size=val
        #BlazarSED.build_photons(self._blob)

    def _get_nu_grid_size(self):
        return  self._blob.nu_grid_size



    def set_verbosity(self,val):
        self._blob.verbose=val

    def get_verbosity(self):
        return  self._blob.verbose

    def debug_synch(self):
        print ("nu stop synch", self._blob.nu_stop_Sync)
        print ("nu stop synch ssc", self._blob.nu_stop_Sync_ssc)
        print ("ID MAX SYNCH", self._blob.NU_INT_STOP_Sync_SSC)

    def debug_SSC(self):
        print ("nu start SSC", self._blob.nu_start_SSC)
        print ("nu stop SSC", self._blob.nu_stop_SSC)
        print ("ID MAX SSC", self._blob.NU_INT_STOP_COMPTON_SSC)

    def set_par(self,par_name,val):
        """
        shortcut to :class:`ModelParametersArray.set` method
        set a parameter value

        :param par_name: (srt), name of the parameter
        :param val: parameter value

        """

        self.parameters.set(par_name, val=val)



    def get_par_by_type(self,par_type):
        """

        get parameter by type

        """
        for param in self.parameters.par_array:
            if param.par_type==par_type:
                return param

        return None

    def get_par_by_name(self,par_name):
        """

        get parameter by type

        """
        for param in self.parameters.par_array:
            if param.name==par_name:
                return param

        return None


    def show_pars(self,sort_key='par type'):
        #self.parameters.par_array.sort(key=lambda x: x.name, reverse=False)
        self.parameters.show_pars(sort_key=sort_key)


    def show_electron_distribution(self):
        print(
            "-------------------------------------------------------------------------------------------------------------------")
        print('electron distribution:')
        print(" type: %s  " % (self._electron_distribution_name))
        print(" electron energy grid size: ", self.gamma_grid_size)
        print(" gmin grid : %e" % self._blob.gmin_griglia)
        print(" gmax grid : %e" % self._blob.gmax_griglia)
        print(" normalization ", self.Norm_distr>0)
        print(" log-values ", self._electron_distribution_log_values)
        print('')
        self.parameters.par_array.sort(key=lambda x: x.name, reverse=False)

        self.parameters.show_pars(names_list=self._electron_distribution_dic.keys())

    def show_model(self):
        """
        shortcut to :class:`ModelParametersArray.show_pars` method
        shows all the paramters in the model

        """
        print("")
        print("-------------------------------------------------------------------------------------------------------------------")

        print ("jet model description")
        print(
            "-------------------------------------------------------------------------------------------------------------------")
        print("name: %s  " % (self.name))
        print('')
        print('electron distribution:')
        print(" type: %s  " % (self._electron_distribution_name))
        print (" electron energy grid size: ",self.gamma_grid_size)
        print (" gmin grid : %e"%self._blob.gmin_griglia)
        print (" gmax grid : %e"%self._blob.gmax_griglia)
        print(" normalization ", self.Norm_distr>0)
        print(" log-values ", self._electron_distribution_log_values)
        print('')
        print('radiative fields:')
        print (" seed photons grid size: ", self.nu_seed_size)
        print (" IC emission grid size: ", self.get_IC_nu_size())
        print (' source emissivity lower bound :  %e' % self._blob.emiss_lim)
        print(' spectral components:')
        for _s in self.spectral_components_list:
            print("   name:%s,"%_s.name, 'state:', _s.state)
        print('external fields transformation method:', self.get_external_field_transf())
        #print('disk type:', self._blob.disk_type)
        print ('')
        print ('SED info:')
        print (' nu grid size :%d' % self._get_nu_grid_size())
        print (' nu mix (Hz): %e' % self._get_nu_min_grid())
        print (' nu max (Hz): %e' % self._get_nu_max_grid())
        print('')
        print('flux plot lower bound   :  %e' % self.flux_plot_lim)
        print('')

        self.show_pars()

        print("-------------------------------------------------------------------------------------------------------------------")

    def plot_model(self,plot_obj=None,clean=False,label=None,comp=None,sed_data=None,color=None,auto_label=True,line_style='-',frame='obs'):
        if plot_obj is None:
            plot_obj=PlotSED(sed_data=sed_data, frame= frame)


        if clean==True:
            plot_obj.clean_model_lines()

        #line_style='-'

        if comp is not None:
            c = self.get_spectral_component_by_name(comp)

            if label is not None:
                comp_label = label
            elif label is None and auto_label is True:
                comp_label = c.name
            else:
                comp_label=None
            if c.state!='off':
                plot_obj.add_model_plot(c.SED, line_style=line_style, label=comp_label,flim=self.flux_plot_lim,color=color,auto_label=auto_label)

        else:
            for c in self.spectral_components_list:
                comp_label = c.name
                if auto_label is not True:
                    comp_label=label
                #print('comp label',comp_label)
                if c.state != 'off' and c.name!='Sum':
                    plot_obj.add_model_plot(c.SED, line_style=line_style, label=comp_label,flim=self.flux_plot_lim,auto_label=auto_label,color=color)


            c=self.get_spectral_component_by_name('Sum')
            if label is not None:
                comp_label = label
            else:
                comp_label='Sum'

            plot_obj.add_model_plot(c.SED, line_style='--', label=comp_label, flim=self.flux_plot_lim,color=color)

        return plot_obj

    @safe_run
    def set_blob(self):
        BlazarSED.Init(self._blob,self.cosmo.get_DL_cm(self.parameters.z_cosm.val))

    @safe_run
    def set_external_fields(self):
        BlazarSED.Init(self._blob,self.cosmo.get_DL_cm(self.parameters.z_cosm.val))
        BlazarSED.spectra_External_Fields(1,self._blob)

    @safe_run
    def eval(self,init=True,fill_SED=True,nu=None,get_model=False,loglog=False,plot=None,label=None,phys_output=False):
        """
        Runs the BlazarSED  code for the current `JetModel` instance.

        :param init: (boolean), "defualt=True" initializes the BlazarSED code
            for the current `Jet` instance parameters values.

        """

        if self.electron_distribution is None:
            raise  RuntimeError('electron distribution not defined')


        if init==True:

            BlazarSED.Init(self._blob,self.cosmo.get_DL_cm(self.parameters.z_cosm.val) )
            #self.set_electron_distribution()

            #TODO investigate if this is necessary!!!
            self._update_spectral_components()

        BlazarSED.Run_SED(self._blob)


        if phys_output==True:
            BlazarSED.EnergeticOutput(self._blob)

        nu_sed_sum,nuFnu_sed_sum= self.spectral_components.Sum.get_SED_points()
        #self.get_SED_points()
        #print('nu_sed_sum,nuFnu_sed_sum',nu_sed_sum,nuFnu_sed_sum)


        if fill_SED==True:
            #TODO check if this is not usefule!!!
            #self.SED.fill(nu=nu_sed_sum,nuFnu=nuFnu_sed_sum)

            for i in range(len(self.spectral_components_list)):

                #print ('fill name',self.spectral_components_list[i].name)
                #nu_sed,nuFnu_sed= self.get_SED_points(name=self.spectral_components_list[i].name)
                #self.spectral_components_list[i].SED.fill(nu=nu_sed, nuFnu=nuFnu_sed)

                self.spectral_components_list[i].fill_SED()

        if nu is None:



            if loglog==True:

                model=log10(nuFnu_sed_sum)
            else:
                model=nuFnu_sed_sum

        else :

            if shape(nu)==():

                nu=array([nu])

            if loglog==False:
                nu_log=log10(nu)
            else:
                nu_log=nu


            nu_sed_log=np.log10(nu_sed_sum)
            nuFnu_sed_log=np.log10(nuFnu_sed_sum)

            msk= nu_log > nu_sed_log.min()


            msk*=  nu_log < nu_sed_log.max()



            if loglog==False:
                model=zeros(nu_log.size)

                nuFnu_log=np.log10(nuFnu_sed_sum)

                f_interp =interp1d(nu_sed_log,nuFnu_sed_log,bounds_error=True)


                model[msk] = power(10.0,f_interp(nu_log[msk]))

            else:
                model=zeros(nu_log.size) + np.log10(self.flux_plot_lim)

                f_interp =interp1d(nu_sed_log,nuFnu_sed_log,bounds_error=True)

                model[msk] = f_interp(nu_log[msk])


        if plot is not None:
            if label is None:
                label= self.name

            self.plot_model(plot, clean=True, label=label)


        if get_model==True:
            return model
        else:
            return None

    @safe_run
    def get_SED_points(self,log_log=False,name='Sum'):

        try:
            spec_comp=self.get_spectral_component_by_name(name)

            #print ('spec_comp',spec_comp.name,spec_comp._nu_name)
            nuFnu_ptr=spec_comp.nuFnu_ptr
            nu_ptr=spec_comp.nu_ptr

            size=self._blob.nu_grid_size
            x=zeros(size)
            y=zeros(size)

            for i in range(size):
                x[i]=BlazarSED.get_spectral_array(nu_ptr,self._blob,i)
                y[i]=BlazarSED.get_spectral_array(nuFnu_ptr,self._blob,i)

                #print("%s %e %e"%(name,x[i],y[i]))

            msk_nan=np.isnan(x)
            msk_nan+=np.isnan(y)
            #print('emiss lim',self.get_emiss_lim())
            x[msk_nan]=0.
            y[msk_nan]=self.get_emiss_lim()

            msk=y<self.get_emiss_lim()


            y[msk]=self.get_emiss_lim()



            if log_log==True:
                msk = y <= 0.
                y[msk] = self.get_emiss_lim()

                #x=x[msk]
                #    y=y[msk]

                x=log10(x)
                y=log10(y)



            return x,y

        except:
            raise RuntimeError ('model evaluation failed in get_SED_points')

    @safe_run
    def energetic_report(self,write_file=False,getstring=True,wd=None,name=None,verbose=True):
        self.energetic_dict={}

        _energetic = BlazarSED.EnergeticOutput(self._blob,0)
        _par_array=ModelParameterArray()

        _name = [i for i in _energetic.__class__.__dict__.keys() if i[:1] != '_']
        units=''
        for _n in _name:
            if _n[0]=='L':
                par_type='Lum. blob rest. frme.'
                units='erg/s'
            elif _n[0] == 'U' and 'DRF' not in _n:
                par_type = 'Energy dens. blob rest. frame'
                units = 'erg/cm^3'
            elif _n[0] == 'U' and 'DRF'  in _n:
                par_type = 'Energy dens. disk rest. frame'
                units = 'erg/cm^3'
            elif _n[0] == 'j':
                par_type = 'jet Lum.'
                units = 'erg/s'
            else:
                warnings.warn('enegetic name %s not understood'%_n)

            self.energetic_dict[_n]=getattr(_energetic, _n)

            _par_array.add_par(ModelParameter(name=_n, val=getattr(_energetic, _n), units=units,par_type=par_type))

            self.energetic_report_table = _par_array.par_table
            self.energetic_report_table.remove_columns(['log','frozen','phys. bound. min','phys. bound. max'])
            self.energetic_report_table.rename_column('par type','type')

        self._energetic_report = self.energetic_report_table.pformat_all()

        if  verbose is True:
            print("-----------------------------------------------------------------------------------------")
            print("jet eneregetic report:")
            self.energetic_report_table.pprint_all()
            print("-----------------------------------------------------------------------------------------")





        if write_file==True:

            if wd is None:
                wd = self.wd

            if name is None:
                name = 'energetic_report_%s' % self.name + '.txt'

            outname = '%s/%s' % (wd, name)

            outfile = open(outname, 'w')

            for text in self._energetic_report:
                print(text, file=outfile)

            outfile.close()



    def get_SED_peak(self,peak_name=None,freq_range=None,log_log=False):

        if peak_name is not None and freq_range is not None:
            print ("either you provide peak_name or freq_range")
            raise ValueError
        elif   peak_name is not None:
            try:
                if log_log==False:
                    return getattr(self._blob,peak_name)
                else:
                     return np.log10(getattr(self._blob,peak_name) )
            except:
                print ("peak name %s, not found, check name"%peak_name)
                raise ValueError

        else:
            x,y=self.spectral_components.Sum.get_SED_points(log_log=log_log)
            msk1=x>freq_range[0]
            msk2=x<freq_range[1]
            #print x,freq_range,x[msk1*msk2]
            y_m= y[msk1*msk2].max()
            x_id= np.argmax(y[msk1*msk2])
            return x[msk1*msk2][x_id],y_m

    def get_component_peak(self,comp_name=None,log_log=False):
        comp = self.get_spectral_component_by_name(comp_name)

        ID = np.argmax(comp.SED.nuFnu.value)
        x_p, y_p = comp.SED.nu[ID].value, comp.SED.nuFnu[ID].value

        if log_log is True:
            x_p, y_p=np.log10([x_p,y_p])

        return x_p, y_p



def set_str_attr(obj,name,val):
    #print('set obj', obj,'name',name ,'to', val)
    try:

        try:
            setattr(obj, name,val)
        except:
            setattr(obj, name, val.encode('ascii'))
    except Exception as e:
        raise RuntimeError('error setting attr',name,'execption:',e)







        
