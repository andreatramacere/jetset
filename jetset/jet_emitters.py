from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)




__author__ = "Andrea Tramacere"



import os
import numpy as np
from numpy import log10,array,zeros,power,shape

from .jetkernel_models_dic import gamma_dic_e ,gamma_dic_p, available_N_distr, N_distr_descr, available_emitters_type
from .plot_sedfit import PlotPdistr
from .jet_paramters import *
from .utils import set_str_attr

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd == True:
    try:
        from jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from jetkernel import jetkernel as BlazarSED


__all__=[ 'EmittersDistribution', 'ArrayDistribution']


class ArrayDistribution(object):

    def __init__(self, e_array, n_array, gamma_grid_size=None):
        self.e_array = e_array
        self.n_array = n_array
        _size = e_array.size

        if n_array.size != _size:
            raise RuntimeError('e_array and n_array must have same size')
        self.size = _size
        if gamma_grid_size is None:
            self.gamma_grid_size = _size * 2 + 1


class EmittersDistribution(object):

    def __init__(self,name,jet,gamma_grid_size=None,log_values=False,emitters_type='electrons'):
        #print('==> __init__ EmittersDistribution 1')
        self.available_emitters_type = available_emitters_type
        self.available_N_distr = available_N_distr

        self._check_emitters_type(emitters_type)
        self.emitters_type=emitters_type

        self._check_N_distr_name(name)
        self.available_N_distr = available_N_distr

        self._name=name
        self._log_values=log_values


        self._jet=jet
        set_str_attr(jet._blob, 'DISTR',name)
        set_str_attr(jet._blob, 'PARTICLE', emitters_type)
        #print('==> __init__ EmittersDistribution 2')
        if gamma_grid_size is not None:
            self.set_grid_size(gamma_grid_size)
            print('==>,EmittersDistribution 2')
        else:
            self._set_blob()
            #print('==>,EmittersDistribution 3')
            self._fill()
            #print('==>,EmittersDistribution 4')

    def _check_emitters_type(self,emitters_type):
        if emitters_type not in available_emitters_type:
            raise RuntimeError ("type of emitters :%s, not allowed" % emitters_type, "please choose among: ",
                          self.available_emitters_type)
        else:
            pass

    def _check_N_distr_name(self,name):
        if name == 'from_array':
            pass
        elif name not in self.available_N_distr:
            raise RuntimeError ("electron distribution model %s not allowed" % name, "please choose among: ", self.available_N_distr)
        else:
            pass

    @staticmethod
    def available_distributions():
        for k in N_distr_descr.keys():
            print('%s: %s'%(k,N_distr_descr[k]))


    @classmethod
    def from_array(cls,jet,custom_Ne,emitters_type='electrons'):
        #TODO update for protons

        if emitters_type !='electrons':
            raise RuntimeError('not implemented yet for',emitters_type, 'only for electrons')

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


        if self.emitters_type=='electrons':
            BlazarSED.InitNe(self._jet._blob)
            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self.Ne_ptr = getattr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = getattr(self._jet._blob, self._gammae_name)

            #print("==> _set_blob 1", self._Ne_name, self._gammae_name, self.Ne_ptr)

        elif self.emitters_type=='protons':
            BlazarSED.Init_Np_Ne_pp(self._jet._blob)
            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self._Np_name, self._gammap_name = gamma_dic_p['proton_distr']
            self.Np_ptr = getattr(self._jet._blob, self._Np_name)
            self.p_gamma_ptr = getattr(self._jet._blob, self._gammap_name)
            self.Ne_ptr = getattr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = getattr(self._jet._blob, self._gammae_name)
            #print("==> _set_blob 1", self._Np_name, self._gammap_name, self.Np_ptr)
        else:
            raise RuntimeError('')


    def _fill(self):
        #print ("==> __fill 1")
        size = self._jet._blob.gamma_grid_size
        self.gamma_e = zeros(size)
        self.n_gamma_e = zeros(size)
        self.gamma_p = zeros(size)
        self.n_gamma_p = zeros(size)

        #print("==> __fill 2")
        for ID in range(size):
            self.gamma_e[ID] = BlazarSED.get_elec_array(self.e_gamma_ptr, self._jet._blob, ID)
            self.n_gamma_e[ID] = BlazarSED.get_elec_array(self.Ne_ptr, self._jet._blob, ID)


        if self.emitters_type=='protons':
            #print("==> __fill 3", self.p_gamma_ptr, self.Np_ptr)
            for ID in range(size):
                self.gamma_p[ID] = BlazarSED.get_elec_array(self.p_gamma_ptr, self._jet._blob, ID)
                #print("==> __fill 4")
                self.n_gamma_p[ID] = BlazarSED.get_elec_array(self.Np_ptr, self._jet._blob, ID)

        #print("==> __fill 5")


    def plot(self, p=None, y_min=None,y_max=None,x_min=None,x_max=None,energy_scale='gamma'):

        self.update()
        if p is None:
            p=PlotPdistr()

        p.plot_distr(self.gamma_e,
                     self.n_gamma_e,
                     y_min=y_min,
                     y_max=y_max,
                     x_min=x_min,
                     x_max=x_max,
                     particle='electrons',
                     energy_scale=energy_scale)

        if  self.emitters_type=='protons':
            p.plot_distr(self.gamma_p,
                         self.n_gamma_p,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='protons',
                         energy_scale=energy_scale)

        return p

    def plot2p(self, p=None, y_min=None, y_max=None, x_min=None, x_max=None,energy_scale='gamma'):
        self.update()

        if p is None:
            p = PlotPdistr()

        p.plot_distr2p(self.gamma_e,
                       self.n_gamma_e,
                       y_min=y_min,
                       y_max=y_max,
                       x_min=x_min,
                       x_max=x_max,
                       particle='electrons',
                       energy_scale=energy_scale)

        if self.emitters_type == 'protons':
            p.plot_distr2p(self.gamma_p,
                           self.n_gamma_p,
                           y_min=y_min,
                           y_max=y_max,
                           x_min=x_min,
                           x_max=x_max,
                           particle='protons',
                           energy_scale=energy_scale)

        return p


    def plot3p(self, p=None,y_min=None,y_max=None,x_min=None,x_max=None,energy_scale='gamma'):
        self.update()

        if p is None:
            p = PlotPdistr()

        p.plot_distr3p(self.gamma_e,
                       self.n_gamma_e,
                       y_min=y_min,
                       y_max=y_max,
                       x_min=x_min,
                       x_max=x_max,
                       particle='electrons',
                       energy_scale=energy_scale)

        if  self.emitters_type=='protons':
            p.plot_distr3p(self.gamma_p,
                           self.n_gamma_p,
                           y_min=y_min,
                           y_max=y_max,
                           x_min=x_min,
                           x_max=x_max,
                           particle='protons',
                           energy_scale=energy_scale)

        return p

    def set_bounds(self,a,b,log_val=False):
        if log_val == False:
            return [a,b]

        else:
            return np.log10([a,b])

    def _build_emitters_distribution_dic(self, distribution_name, emitters_type='electrons'):
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
        if emitters_type =='electrons':
            model_dic['N'] = JetModelDictionaryPar(ptype='electron_density',vmin=0,vmax=None,punit='cm^-3')
            #['electron_density', 0, None, 'cm^-3']

        elif emitters_type =='protons':
            model_dic['N'] = JetModelDictionaryPar(ptype='proton_density', vmin=0, vmax=None, punit='cm^-3')
            #model_dic['N_e_pp'] = JetModelDictionaryPar(ptype='electron_density', vmin=0, vmax=None, punit='cm^-3',froz='False')
            model_dic['NH_pp'] = JetModelDictionaryPar(ptype='proton_density', vmin=0, vmax=None, punit='cm^-3',froz ='False')

        a_h,b_h=self.set_bounds(1,1E15,log_val=self._log_values)
        a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
        a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)
        model_dic['gmin'] = JetModelDictionaryPar(ptype='low-energy-cut-off',vmin=a_l,vmax=b_l,punit='lorentz-factor',log=self._log_values)
        #['low-energy-cut-off', a_l, b_l, 'Lorentz-factor',False,self._log_values]
        model_dic['gmax'] = JetModelDictionaryPar(ptype='high-energy-cut-off',vmin=a_h,vmax=b_h,punit='lorentz-factor',log=self._log_values)
        #['high-energy-cut-off', a_h, b_h, 'Lorentz-factor',False,self._log_values]

        if distribution_name == 'pl':
            model_dic['p'] = JetModelDictionaryPar(ptype='HE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['HE_spectral_slope', -10, 10, '']

        if distribution_name == 'bkn':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['p_1'] = JetModelDictionaryPar(ptype='HE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['HE_spectral_slope', -10, 10, '']
            model_dic['gamma_break'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if distribution_name == 'lp':
            model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
            # ['LE_spectral_slope', -10, 10, '']
            model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature', vmin=-15., vmax=15., punit='')
            # ['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t,punit='lorentz-factor', log=self._log_values,froz=True)

        if distribution_name == 'lppl':
            model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            #['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if distribution_name == 'lpep':
            model_dic['r'] =  JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            #['spectral_curvature', -15, 15, '']
            model_dic['gammap_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if distribution_name == 'plc':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['gamma_cut'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]

        if distribution_name == 'spitkov':
            model_dic['spit_index'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            #['LE_spectral_slope', -10, 10, '']
            model_dic['spit_temp'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['spit_gamma_th'] =JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            #['turn-over-energy', 1, None, 'Lorentz-factor']

        if distribution_name == 'lppl_pile_up':
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


        if distribution_name == 'bkn_pile_up':
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

        if distribution_name == 'from_array':
            model_dic.pop('gmin')
            model_dic.pop('gmax')

        return model_dic



#
# class ArrayElectronDistribution(object):
#
#     def __init__(self, e_array, n_array, gamma_grid_size=None):
#         self.e_array = e_array
#         self.n_array = n_array
#         _size = e_array.size
#
#         if n_array.size != _size:
#             raise RuntimeError('e_array and n_array must have same size')
#         self.size = _size
#         if gamma_grid_size is None:
#             self.gamma_grid_size = _size * 2 + 1
#
#
# class ElectronDistribution(object):
#
#     def __init__(self,name,jet,gamma_grid_size=None,log_values=False):
#
#         self.elec_models_list = available_N_distr
#
#         if name == 'from_array':
#             pass
#         elif name not in self.elec_models_list:
#             raise Warning("electron distribution model %s not allowed" % name)
#             print("please choose among: ", self.elec_models_list)
#             return None
#         else:
#             pass
#
#         self._name=name
#         self._log_values=log_values
#
#
#         self._jet=jet
#         set_str_attr(jet._blob,'DISTR',name)
#
#         if gamma_grid_size is not None:
#             self.set_grid_size(gamma_grid_size)
#
#         else:
#             self._set_blob()
#             self._fill()
#
#
#     @staticmethod
#     def available_distributions():
#         for k in N_distr_descr.keys():
#             print('%s: %s'%(k,N_distr_descr[k]))
#
#
#     @classmethod
#     def from_array(cls,jet,custom_Ne):
#         name='from_array'
#
#
#         BlazarSED.build_Ne_custom(jet._blob,custom_Ne.size)
#         Ne_custom_ptr = getattr(jet._blob, 'Ne_custom')
#         gamma_custom_ptr = getattr(jet._blob,'gamma_e_custom')
#         for ID in range(custom_Ne.size):
#             BlazarSED.set_elec_custom_array(gamma_custom_ptr,jet._blob,custom_Ne.e_array[ID], ID)
#             BlazarSED.set_elec_custom_array(Ne_custom_ptr, jet._blob, custom_Ne.n_array[ID],ID)
#
#         setattr(jet._blob,'gmin',custom_Ne.e_array[0] )
#         setattr(jet._blob,'gmax',custom_Ne.e_array[-1] )
#
#         return cls(name,jet,custom_Ne.gamma_grid_size)
#
#
#     def update(self):
#         self._set_blob()
#         self._fill()
#
#     def set_grid_size(self,gamma_grid_size):
#         setattr(self._jet._blob,'gamma_grid_size' ,gamma_grid_size)
#         self._set_blob()
#         self._fill()
#
#     def _set_blob(self):
#         #BlazarSED.MakeNe(self._jet._blob)
#         BlazarSED.InitNe(self._jet._blob)
#         self._N_name, self._gamma_name = gamma_dic['electron_distr']
#         self.Ne_ptr = getattr(self._jet._blob, self._N_name)
#         self.gamma_ptr = getattr(self._jet._blob, self._gamma_name)
#
#
#     def _fill(self):
#         size = self._jet._blob.gamma_grid_size
#         self.gamma = zeros(size)
#         self.n_gamma = zeros(size)
#
#         for ID in range(size):
#             self.gamma[ID] = BlazarSED.get_elec_array(self.gamma_ptr, self._jet._blob, ID)
#             self.n_gamma[ID] = BlazarSED.get_elec_array(self.Ne_ptr, self._jet._blob, ID)
#
#
#     def plot(self, p=None, y_min=None,y_max=None,x_min=None,x_max=None):
#
#         self.update()
#         if p is None:
#             p=PlotPdistr()
#
#         p.plot_distr(self.gamma,self.n_gamma,y_min=y_min,y_max=y_max,x_min=x_min,x_max=x_max)
#
#         return p
#
#
#     def plot3p(self, p=None,y_min=None,y_max=None,x_min=None,x_max=None):
#         self.update()
#
#         if p is None:
#             p = PlotPdistr()
#
#         p.plot_distr3p(self.gamma,self.n_gamma,y_min=y_min,y_max=y_max,x_min=x_min,x_max=x_max)
#
#         return p
#
#     def set_bounds(self,a,b,log_val=False):
#         if log_val == False:
#             return [a,b]
#
#         else:
#             return np.log10([a,b])
#
#     def _build_electron_distribution_dic(self,electron_distribution_name):
#         """
#         Builds the dictionary to init the :class:`.JetParameter`
#         objects   for the electron  distribution:
#
#         The following :class:`.JetParameter`: objects the *do not* depend
#         on the type of electron distribution
#
#                 - N, particle density in cm^-3
#
#                 - gmin, the minimum value of the electron Lorentz factor
#
#                 - gmax, the maximum value of the electron Lorentz factor
#
#         The following :class:`.JetParameter`: objects *depend* on the type of electron distribution:
#
#             - **power law**, electron_distribution='pl'
#
#                - p
#
#             - **broken power-law**, electron_distribution= **'bkn'**
#
#                 - p
#                 - p_1
#                 - gamma_break
#
#             - **log-parabola**, electron_distribution= **'lp'**
#
#                 - r
#                 - s
#                 - gamma0_log_parab (fixed)
#
#             - **log-parabola** with a low energy power-law tail, electron_distribution= **'lppl'**
#
#                 - r
#                 - s
#                 - gamma0_log_parab
#
#             - **log-parabola** defined by peak energy, electron_distribution= **'lpep'**
#
#                 - r
#                 - s
#                 - gammap_log_parab,
#
#             - **power-law cut-off**, lectron_distribution= **'plc'**
#
#                 - p
#                 - gamma_cut
#
#         """
#
#         model_dic = {}
#         model_dic['N'] = JetModelDictionaryPar(ptype='electron_density', vmin=0, vmax=None, punit='cm^-3')
#         #['electron_density', 0, None, 'cm^-3']
#
#         a_h,b_h=self.set_bounds(1,1E15,log_val=self._log_values)
#         a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
#         a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)
#         model_dic['gmin'] = JetModelDictionaryPar(ptype='low-energy-cut-off', vmin=a_l, vmax=b_l, punit='lorentz-factor', log=self._log_values)
#         #['low-energy-cut-off', a_l, b_l, 'Lorentz-factor',False,self._log_values]
#         model_dic['gmax'] = JetModelDictionaryPar(ptype='high-energy-cut-off', vmin=a_h, vmax=b_h, punit='lorentz-factor', log=self._log_values)
#         #['high-energy-cut-off', a_h, b_h, 'Lorentz-factor',False,self._log_values]
#
#         if electron_distribution_name == 'pl':
#             model_dic['p'] = JetModelDictionaryPar(ptype='HE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['HE_spectral_slope', -10, 10, '']
#
#         if electron_distribution_name == 'bkn':
#             model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['LE_spectral_slope', -10, 10, '']
#             model_dic['p_1'] = JetModelDictionaryPar(ptype='HE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['HE_spectral_slope', -10, 10, '']
#             model_dic['gamma_break'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]
#
#         if electron_distribution_name == 'lp':
#             model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             # ['LE_spectral_slope', -10, 10, '']
#             model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature', vmin=-15., vmax=15., punit='')
#             # ['spectral_curvature', -15, 15, '']
#             model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values, froz=True)
#
#         if electron_distribution_name == 'lppl':
#             model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['LE_spectral_slope', -10, 10, '']
#             model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature', vmin=-15., vmax=15., punit='')
#             #['spectral_curvature', -15, 15, '']
#             model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]
#
#         if electron_distribution_name == 'lpep':
#             model_dic['r'] =  JetModelDictionaryPar(ptype='spectral_curvature', vmin=-15., vmax=15., punit='')
#             #['spectral_curvature', -15, 15, '']
#             model_dic['gammap_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]
#
#         if electron_distribution_name == 'plc':
#             model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['LE_spectral_slope', -10, 10, '']
#             model_dic['gamma_cut'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', a_t, b_t, 'Lorentz-factor',False,self._log_values]
#
#         if electron_distribution_name == 'spitkov':
#             model_dic['spit_index'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['LE_spectral_slope', -10, 10, '']
#             model_dic['spit_temp'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', 1, None, 'Lorentz-factor']
#             model_dic['spit_gamma_th'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', 1, None, 'Lorentz-factor']
#
#         if electron_distribution_name == 'lppl_pile_up':
#             model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['LE_spectral_slope', -10, 10, '']
#             model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature', vmin=-15., vmax=15., punit='')
#             # ['spectral_curvature', -15, 15, '']
#             model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', 1, None, 'Lorentz-factor']
#             model_dic['gamma_inj'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', 1, None, 'Lorentz-factor']
#
#             model_dic['gamma_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#            # ['turn-over-energy', 1, None, 'Lorentz-factor']
#             model_dic['ratio_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=0.0, vmax=None, punit='')
#             #['turn-over-energy', 0, None, '']
#
#             model_dic['alpha_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=0.0, vmax=10.0, punit='')
#             #['turn-over-energy', 0.0, 10, '']
#
#
#         if electron_distribution_name == 'bkn_pile_up':
#             model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             #['LE_spectral_slope', -10, 10, '']
#             model_dic['p_1'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
#             # ['HE_spectral_slope', -10, 10, '']
#             model_dic['gamma_break'] =  JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', 1, None, 'Lorentz-factor']
#
#             model_dic['gamma_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', 1, None, 'Lorentz-factor']
#             model_dic['gamma_pile_up_cut'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t, punit='lorentz-factor', log=self._log_values)
#             #['turn-over-energy', 1, None, 'Lorentz-factor']
#
#             model_dic['alpha_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=0.0, vmax=10.0, punit='')
#             #['turn-over-energy', 0.0, 10, '']
#
#         if electron_distribution_name == 'from_array':
#             model_dic.pop('gmin')
#             model_dic.pop('gmax')
#
#         return model_dic
